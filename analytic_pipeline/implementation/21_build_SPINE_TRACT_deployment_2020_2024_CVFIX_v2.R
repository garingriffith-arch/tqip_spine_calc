# ==============================================================================
# SPINE-TRACT Stage 21: Locked 2020-2024 deployment models
# Executed CV-parser revision v2; XGBoost as_pandas=FALSE.
# Model development source: repaired-DX Stage 3A-R and Stage 3C-R.
#
# DOES NOT ALTER manuscript temporal models or re-estimate published metrics.
# Five R XGBoost classifiers + three Python XGBoost joint quantile boosters.
# 2024-selected structural configuration is LOCKED per endpoint.
# 3-fold internal CV over the entire eligible 2020-2024 cohort selects only the
# number of boosting rounds for the all-years deployment model.
# CV is for deployment fitting, NOT a manuscript validation performance estimate.
#
# OUTPUT: <project>/outputs/stage6_deployment_2020_2024/artifacts/
#          models/*.json, encoders/*_encoder.rds, qc/*.csv, artifact_manifest.csv
# No patient-level data or predictions are written.
#
# Required R packages: data.table, Matrix, xgboost, reticulate, digest
# Python (duration stage): xgboost >= 2.0 in SPINE_TRACT_PY_ENV (default
#                        tbi-tract-xgb-gpu).
#
# Environment variables (optional; examples in the accompanying README):
#   SPINE_TRACT_PROJECT_DIR         workspace root (defaults to known local path)
#   SPINE_TRACT_DEPLOY_STAGE        classification | duration | all [all]
#   SPINE_TRACT_DEPLOY_DEVICE       cuda | cpu [cuda]
#   SPINE_TRACT_PY_ENV              conda environment [tbi-tract-xgb-gpu]
#   SPINE_TRACT_DEPLOY_MAX_ROUNDS   max CV boosting rounds [6000]
#   SPINE_TRACT_DEPLOY_NFOLD        folds, >=3 [3]
#   SPINE_TRACT_DEPLOY_THREADS      R/Python threads [8]
#   SPINE_TRACT_DEPLOY_FORCE        1 to overwrite completed endpoints [0]
#   SPINE_TRACT_DEPLOY_SEED         random seed [20261019]
# ==============================================================================

options(stringsAsFactors = FALSE)
required <- c("data.table", "Matrix", "xgboost", "reticulate", "digest")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x
get_int <- function(k, default) {
  z <- suppressWarnings(as.integer(Sys.getenv(k, unset = as.character(default))))
  if (length(z) != 1L || is.na(z)) stop("Invalid integer environment variable: ", k)
  z
}

PROJECT_DIR <- Sys.getenv("SPINE_TRACT_PROJECT_DIR", unset =
  "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization")
STAGE <- tolower(Sys.getenv("SPINE_TRACT_DEPLOY_STAGE", unset = "all"))
DEVICE <- tolower(Sys.getenv("SPINE_TRACT_DEPLOY_DEVICE", unset = "cuda"))
PY_ENV <- Sys.getenv("SPINE_TRACT_PY_ENV", unset = "tbi-tract-xgb-gpu")
MAX_ROUNDS <- get_int("SPINE_TRACT_DEPLOY_MAX_ROUNDS", 6000L)
NFOLD <- get_int("SPINE_TRACT_DEPLOY_NFOLD", 3L)
THREADS <- get_int("SPINE_TRACT_DEPLOY_THREADS", 8L)
SEED <- get_int("SPINE_TRACT_DEPLOY_SEED", 20261019L)
FORCE <- identical(Sys.getenv("SPINE_TRACT_DEPLOY_FORCE", unset = "0"), "1")
if (!STAGE %in% c("all", "classification", "duration")) stop("Invalid SPINE_TRACT_DEPLOY_STAGE")
if (!DEVICE %in% c("cpu", "cuda")) stop("Device must be cpu or cuda")
if (NFOLD < 3L || THREADS < 1L || MAX_ROUNDS < 500L) stop("Invalid CV parameters")
if (!dir.exists(PROJECT_DIR)) stop("Project directory does not exist: ", PROJECT_DIR)

DATA_FILE <- file.path(PROJECT_DIR, "data", "derived", "stage2_predictor_adjudication",
                       "spine_tract_frozen_primary_2020_2024.rds")
SPEC_FILE <- file.path(PROJECT_DIR, "data", "derived", "stage2_predictor_adjudication",
                       "spine_tract_predictor_sets.rds")
CLASS_SOURCE <- file.path(PROJECT_DIR, "outputs", "stage3_final_classification_REPAIRED_DX",
                          "rolling_origin_gpu", "14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv")
DURATION_SOURCE <- file.path(PROJECT_DIR, "outputs", "stage3_duration_quantiles_REPAIRED_DX",
                             "rolling_origin_gpu", "04_SELECTED_DURATION_CONFIGS.csv")
OUT <- file.path(PROJECT_DIR, "outputs", "stage6_deployment_2020_2024", "artifacts")
MODELS_DIR <- file.path(OUT, "models")
ENC_DIR <- file.path(OUT, "encoders")
QC_DIR <- file.path(OUT, "qc")
for (f in c(MODELS_DIR, ENC_DIR, QC_DIR)) dir.create(f, recursive = TRUE, showWarnings = FALSE)

needed <- c(DATA_FILE, SPEC_FILE, if (STAGE != "duration") CLASS_SOURCE,
            if (STAGE != "classification") DURATION_SOURCE)
if (any(!file.exists(needed))) stop("Missing prerequisites:\n", paste(needed[!file.exists(needed)], collapse = "\n"))

cat("\nSPINE-TRACT FINAL ALL-YEARS DEPLOYMENT\n")
cat("  Workspace: ", PROJECT_DIR, "\n  Stage: ", STAGE,
    "\n  Device: ", DEVICE, "\n  Internal CV: ", NFOLD,
    " folds, max ", MAX_ROUNDS, " rounds\n", sep = "")
cat("  Output: ", OUT, "\n\n", sep = "")

DT <- as.data.table(readRDS(DATA_FILE))
SPEC <- readRDS(SPEC_FILE)
stopifnot(nrow(DT) == 424017L,
          all(c("inc_key", "admission_year") %in% names(DT)))
DT[, admission_year := as.integer(admission_year)]
expected_by_year <- c(`2020`=75237L, `2021`=83076L, `2022`=85301L,
                      `2023`=89754L, `2024`=90649L)
observed_by_year <- table(factor(DT$admission_year, levels = 2020:2024))
if (!identical(as.integer(observed_by_year), unname(expected_by_year))) {
  stop("Year-specific repaired primary cohort counts do not match the locked manuscript")
}
pg <- SPEC$predictor_groups
backbone <- SPEC$predictor_sets$clinical_backbone
CORE <- unique(c(setdiff(backbone, pg$spine_granular), pg$spine_coarse,
                 pg$chest_detail, pg$pupil))
CONTEXT <- unique(c(CORE, pg$social_context))
if (length(CORE) != 54L || length(CONTEXT) != 57L ||
    !"dx_tl_soft_tissue_any" %in% CORE ||
    !setequal(setdiff(CONTEXT, CORE), pg$social_context)) {
  stop("Final repaired 54/57 predictor lock is not reconstructible")
}
if (length(setdiff(CONTEXT, names(DT)))) stop("Missing predictors: ",
  paste(setdiff(CONTEXT, names(DT)), collapse = ", "))

# These definitions intentionally match the repaired Stage 3A and Stage 3C source.
DT[, discharge3_final := factor(fcase(
  as.character(discharge4) == "Home/home health", "Home/home health",
  as.character(discharge4) %in% c("Inpatient rehabilitation", "Other post-acute facility"),
    "Post-acute facility",
  as.character(discharge4) == "Death/hospice", "Death/hospice",
  default = NA_character_),
  levels = c("Home/home health", "Post-acute facility", "Death/hospice"))]
DT[, hlos3_final := factor(as.character(hlos3), levels = c("<=7 days", "8-19 days", ">=20 days"))]
DT[, icu3_final := factor(as.character(icu3), levels = c("No ICU", "ICU 1-7 days", "ICU >=8 days"))]
DT[, vent3_final := factor(as.character(vent3),
  levels = c("No ventilation", "Ventilation 1-7 days", "Ventilation >=8 days"))]

CLASS_SPECS <- list(
  discharge3_final = list(type="multiclass", policy="final_context",
                         predictors=CONTEXT, levels=levels(DT$discharge3_final)),
  hlos3_final = list(type="multiclass", policy="final_context",
                     predictors=CONTEXT, levels=levels(DT$hlos3_final)),
  icu3_final = list(type="multiclass", policy="final_core",
                    predictors=CORE, levels=levels(DT$icu3_final)),
  vent3_final = list(type="multiclass", policy="final_core",
                     predictors=CORE, levels=levels(DT$vent3_final)),
  operative_spine_intervention = list(type="binary", policy="final_core",
                                     predictors=CORE, levels=c("No operation", "Recorded operation"))
)
DURATION_SPECS <- list(
  hospital_los = list(y_var="hospital_los_days", lower_bound=0),
  icu_los_conditional = list(y_var="icu_los_days", lower_bound=1),
  ventilator_days_conditional = list(y_var="ventilator_days_clean", lower_bound=1)
)

# Frozen preprocessing; no reference levels are dropped in the one-hot encoder.
normalize_category <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | z == "" | z %in% c("NA", "NaN", "<NA>")] <- "__UNKNOWN__"
  z
}
fit_encoder <- function(d, predictors) {
  enc <- setNames(vector("list", length(predictors)), predictors)
  for (v in predictors) {
    x <- d[[v]]
    if (is.numeric(x) || is.integer(x)) {
      enc[[v]] <- list(type="numeric")
    } else {
      z <- normalize_category(x)
      enc[[v]] <- list(type="categorical",
                       levels=unique(c(sort(unique(z)), "__UNKNOWN__", "__OTHER__")))
    }
  }
  enc
}
encode_features <- function(d, predictors, encoder, dense = FALSE) {
  n <- nrow(d)
  blocks <- vector("list", length(predictors))
  for (i in seq_along(predictors)) {
    v <- predictors[[i]]
    info <- encoder[[v]]
    if (is.null(info)) stop("Encoder missing predictor ", v)
    if (identical(info$type, "numeric")) {
      if (dense) {
        m <- matrix(as.numeric(d[[v]]), nrow=n, ncol=1L)
      } else {
        m <- Matrix::Matrix(as.numeric(d[[v]]), ncol=1L, sparse=TRUE)
      }
      colnames(m) <- v
    } else {
      z <- normalize_category(d[[v]])
      z[!z %in% info$levels] <- "__OTHER__"
      f <- factor(z, levels=info$levels)
      if (anyNA(f)) stop("Failed to encode categorical predictor: ", v)
      if (dense) {
        m <- matrix(0, nrow=n, ncol=length(info$levels))
        m[cbind(seq_len(n), as.integer(f))] <- 1
      } else {
        m <- Matrix::sparseMatrix(i=seq_len(n), j=as.integer(f), x=1,
                                  dims=c(n, length(info$levels)))
      }
      colnames(m) <- paste0(v, "__", make.names(info$levels, unique=TRUE))
    }
    blocks[[i]] <- m
  }
  X <- do.call(cbind, blocks)
  if (dense) storage.mode(X) <- "double"
  X
}

# Preserve both Stage 3 per-predictor encoder and TBI-style metadata fields.
make_metadata <- function(id, spec, encoder, feature_names, rounds, config, n_rows,
                          type, source_file) {
  cats <- names(encoder)[vapply(encoder, function(z) identical(z$type, "categorical"), logical(1))]
  nums <- names(encoder)[vapply(encoder, function(z) identical(z$type, "numeric"), logical(1))]
  list(
    endpoint_id=id,
    endpoint_label=gsub("_", " ", id),
    type=type,
    predictors=spec$predictors,
    encoder=list(numeric_vars=nums, categorical_vars=cats,
                 categorical_levels=setNames(lapply(encoder[cats], `[[`, "levels"), cats),
                 per_predictor=encoder),
    feature_names=feature_names,
    predictor_policy=spec$policy,
    outcome_levels=spec$levels %||% NULL,
    outcome_var=spec$y_var %||% id,
    quantile_alphas=if (identical(type, "joint_quantile")) c(.10,.50,.90) else NULL,
    transformation=spec$transformation %||% NULL,
    lower_bound=spec$lower_bound %||% NULL,
    selected_config=as.list(config),
    selected_rounds=as.integer(rounds),
    training_n=as.integer(n_rows),
    training_years=2020:2024,
    method="Locked 2024 structural config; 3-fold all-years internal CV for boosting rounds; fit all eligible years",
    seed=SEED,
    runtime=if (identical(type, "joint_quantile")) "Python XGBoost >=2 (reticulate)" else "R XGBoost",
    selected_configuration_source=source_file,
    built_at=as.character(Sys.time())
  )
}

sha256 <- function(path) digest::digest(object=path, algo="sha256", file=TRUE, serialize=FALSE)
source_hash <- list(
  frozen_dataset=unname(tools::md5sum(DATA_FILE)),
  predictor_spec=unname(tools::md5sum(SPEC_FILE)),
  classifier_selected=if (file.exists(CLASS_SOURCE)) unname(tools::md5sum(CLASS_SOURCE)) else NA_character_,
  duration_selected=if (file.exists(DURATION_SOURCE)) unname(tools::md5sum(DURATION_SOURCE)) else NA_character_
)
write.csv(as.data.frame(source_hash), file.path(QC_DIR, "source_file_checksums.csv"), row.names=FALSE)

get_config <- function(path, id, label_col) {
  tab <- data.table::fread(path)
  must <- c(label_col, "evaluation_year", "max_depth", "min_child_weight",
            "subsample", "colsample_bytree", "lambda", "best_iteration")
  if (length(setdiff(must, names(tab)))) stop("Missing selected-config columns in ", path, ": ",
                                             paste(setdiff(must, names(tab)), collapse=", "))
  row <- tab[get(label_col) == id & evaluation_year == 2024L]
  if (nrow(row) != 1L) stop("Expected exactly one 2024 selected configuration for ", id)
  numeric_cols <- setdiff(must, c(label_col, "evaluation_year"))
  if (!all(is.finite(as.numeric(unlist(row[, ..numeric_cols]))))) {
    stop("Non-finite selected configuration for ", id)
  }
  row
}

qc_file <- function(id) file.path(QC_DIR, paste0(id, "_serialization_qc.csv"))
model_file <- function(id) file.path(MODELS_DIR, paste0(id, ".json"))
meta_file <- function(id) file.path(ENC_DIR, paste0(id, "_encoder.rds"))

is_completed <- function(id) {
  paths <- c(model_file(id), meta_file(id), qc_file(id))
  if (FORCE) return(FALSE)
  if (all(!file.exists(paths))) return(FALSE)
  if (!all(file.exists(paths))) stop("Partial artifact files found for ", id,
    ". Remove all three corresponding files to rerun, or set SPINE_TRACT_DEPLOY_FORCE=1")
  qc <- data.table::fread(qc_file(id))
  if (nrow(qc) != 1L || !all(c("model_sha256", "meta_sha256", "max_abs_logit_or_quantile_difference") %in% names(qc))) {
    stop("Malformed existing QC for ", id)
  }
  if (qc$model_sha256 != sha256(model_file(id)) ||
      qc$meta_sha256 != sha256(meta_file(id)) ||
      !is.finite(qc$max_abs_logit_or_quantile_difference) ||
      qc$max_abs_logit_or_quantile_difference > 1e-5) {
    stop("Artifact checksum or parity QC mismatch for ", id, "; do not reuse it")
  }
  TRUE
}

record_qc <- function(id, kind, n, rounds, selected, delta, model_path, metadata_path) {
  if (!is.finite(delta) || delta > 1e-5) stop("Serialization parity failure: ", id, " (", delta, ")")
  obj <- data.table::data.table(
    endpoint_id=id, type=kind, training_n=n, training_years="2020-2024",
    selected_2024_config_id=as.character(selected$config_id %||% NA_character_),
    selected_rounds=as.integer(rounds), internal_cv_folds=NFOLD,
    max_abs_logit_or_quantile_difference=delta,
    model_sha256=sha256(model_path), meta_sha256=sha256(metadata_path),
    model_size_bytes=file.info(model_path)$size,
    model_format="XGBoost native JSON", built_at=as.character(Sys.time()),
    source_file_md5=source_hash$frozen_dataset)
  data.table::fwrite(obj, qc_file(id))
}

select_rounds_from_R_cv <- function(dmat, params, id) {
  set.seed(SEED)
  cv <- xgboost::xgb.cv(params=params, data=dmat, nrounds=MAX_ROUNDS,
                       nfold=NFOLD, stratified=TRUE, early_stopping_rounds=75L,
                       prediction=FALSE, verbose=0, showsd=TRUE)
  log <- as.data.frame(cv$evaluation_log)
  test_cols <- grep("^test_.*_mean$", names(log), value=TRUE)
  if (length(test_cols) != 1L) stop("Unexpected classification CV metric columns for ", id,
                                    ": ", paste(names(log), collapse=", "))
  rounds <- which.min(log[[test_cols]])
  if (rounds > MAX_ROUNDS - 15L) stop("Classification CV hit the boosting-round cap for ", id,
    ". Increase SPINE_TRACT_DEPLOY_MAX_ROUNDS; do not publish capped model")
  data.table::fwrite(as.data.table(log), file.path(QC_DIR, paste0(id, "_cv_history.csv")))
  as.integer(rounds)
}

fit_classification <- function(id) {
  cat("\n===== CLASSIFICATION: ", id, " =====\n", sep="")
  if (is_completed(id)) { cat("Verified files already present; skipping\n"); return(invisible(NULL)) }
  sp <- CLASS_SPECS[[id]]
  if (is.null(sp)) stop("Unknown classification endpoint")
  selected <- get_config(CLASS_SOURCE, id, "endpoint")
  if (sp$type == "multiclass") {
    y <- as.integer(DT[[id]]) - 1L
    ok <- !is.na(y)
    if (!all(sort(unique(y[ok])) == 0:2)) stop("Invalid 3-class outcome for ", id)
  } else {
    y <- suppressWarnings(as.integer(as.character(DT[[id]])))
    ok <- !is.na(y) & y %in% 0:1
    if (!identical(sort(unique(y[ok])), 0:1)) stop("Invalid binary outcome for ", id)
  }
  if (sum(ok) < 1000L) stop("Unexpectedly small deployment cohort for ", id)
  predictors <- sp$predictors
  enc <- fit_encoder(DT[ok], predictors)
  X <- encode_features(DT[ok], predictors, enc, dense=FALSE)
  feature_names <- colnames(X)
  if (anyDuplicated(feature_names)) stop("Duplicate encoded feature names for ", id)
  dmat <- xgboost::xgb.DMatrix(data=X, label=y[ok], missing=NA_real_)
  params <- list(
    objective=if (sp$type == "multiclass") "multi:softprob" else "binary:logistic",
    eval_metric=if (sp$type == "multiclass") "mlogloss" else "logloss",
    eta=0.05, max_depth=as.integer(selected$max_depth),
    min_child_weight=as.numeric(selected$min_child_weight),
    subsample=as.numeric(selected$subsample),
    colsample_bytree=as.numeric(selected$colsample_bytree),
    lambda=as.numeric(selected$lambda), nthread=THREADS, verbosity=0,
    tree_method=if (DEVICE == "cuda") "gpu_hist" else "hist")
  if (sp$type == "multiclass") params$num_class <- length(sp$levels)
  rounds <- select_rounds_from_R_cv(dmat, params, id)
  cat("All-years internal CV selected boosting rounds: ", rounds, "\n")
  set.seed(SEED + 10000L)
  fit <- xgboost::xgb.train(params=params, data=dmat, nrounds=rounds, verbose=0)
  set.seed(SEED + 21L)
  test_idx <- sample.int(nrow(X), min(128L, nrow(X)))
  dtest <- xgboost::xgb.DMatrix(X[test_idx,,drop=FALSE], missing=NA_real_)
  before <- as.numeric(predict(fit, dtest))
  path <- model_file(id)
  xgboost::xgb.save(fit, fname=path)
  reloaded <- xgboost::xgb.load(path)
  after <- as.numeric(predict(reloaded, dtest))
  if (length(before) != length(after)) stop("Serialized classification output shape changed: ", id)
  delta <- max(abs(before - after))
  md <- make_metadata(id, sp, enc, feature_names, rounds, selected,
                      sum(ok), sp$type, CLASS_SOURCE)
  md$source_model_engine_version <- as.character(utils::packageVersion("xgboost"))
  saveRDS(md, meta_file(id), version=3, compress=TRUE)
  record_qc(id, sp$type, sum(ok), rounds, selected, delta, path, meta_file(id))
  cat("JSON/RDS parity OK, maximum difference: ", format(delta, scientific=TRUE), "\n")
  rm(X,dmat,dtest,fit,reloaded,md);gc(verbose=FALSE)
}

quantile_matrix <- function(p, n) {
  if (inherits(p, "python.builtin.object")) p <- reticulate::py_to_r(p)
  if (is.matrix(p) || is.data.frame(p)) {
    m <- as.matrix(p)
  } else {
    m <- matrix(as.numeric(p), nrow=n, ncol=3L, byrow=TRUE)
  }
  if (!identical(as.integer(dim(m)), c(as.integer(n),3L))) stop("Invalid joint quantile shape: ",
    paste(dim(m),collapse="x"))
  storage.mode(m) <- "double"
  m
}

# Python XGBoost CV uses dictionary output when as_pandas = FALSE.
# Inspect the test mean metric rather than allowing reticulate to flatten it.
select_rounds_from_python_cv <- function(xgb, dmat, params, id) {
  cat("Running three-fold joint-quantile CV for ", id, "\n")
  cv <- xgb$cv(
    params = params,
    dtrain = dmat,
    num_boost_round = as.integer(MAX_ROUNDS),
    nfold = as.integer(NFOLD),
    stratified = FALSE,
    early_stopping_rounds = as.integer(100L),
    seed = as.integer(SEED),
    shuffle = TRUE,
    verbose_eval = FALSE,
    as_pandas = FALSE
  )
  raw <- tryCatch(reticulate::py_to_r(cv), error = function(e) cv)
  if (is.data.frame(raw)) {
    history <- data.table::as.data.table(raw)
  } else if (is.list(raw) && !is.null(names(raw)) &&
             any(grepl("^(train|test)[-_.]", names(raw)))) {
    numeric_columns <- lapply(raw, function(z) {
      as.numeric(unlist(z, recursive = TRUE, use.names = FALSE))
    })
    if (!length(numeric_columns) || length(unique(lengths(numeric_columns))) != 1L) {
      stop("Inconsistent joint-quantile CV series lengths for ", id)
    }
    history <- data.table::as.data.table(numeric_columns, check.names = FALSE)
  } else if (is.matrix(raw) && is.numeric(raw) && ncol(raw) == 4L) {
    history <- data.table::as.data.table(raw)
    data.table::setnames(history, c(
      "train-quantile-mean", "train-quantile-std",
      "test-quantile-mean", "test-quantile-std"
    ))
  } else {
    stop("Unrecognized Python XGBoost CV return structure for ", id,
         ". R class: ", paste(class(raw), collapse = "/"),
         "; names: ", paste(head(names(raw), 8L), collapse = ", "))
  }
  normalized <- gsub("[._]", "-", names(history))
  test_cols <- which(grepl("^test-.*-mean$", normalized))
  if (length(test_cols) != 1L) {
    stop("Expected exactly one test mean CV metric for ", id,
         "; observed columns: ", paste(names(history), collapse = ", "))
  }
  metric_name <- names(history)[test_cols]
  scores <- suppressWarnings(as.numeric(history[[test_cols]]))
  if (!length(scores) || any(!is.finite(scores))) {
    stop("Non-finite or empty test CV metric history for ", id)
  }
  rounds <- as.integer(which.min(scores))
  if (rounds < 1L || rounds > nrow(history)) {
    stop("Invalid optimal boosting rounds for ", id)
  }
  if (rounds > MAX_ROUNDS - 15L) {
    stop("Joint-quantile CV hit the round cap for ", id,
         ". Increase SPINE_TRACT_DEPLOY_MAX_ROUNDS; do not publish capped model")
  }
  history[, boosting_round := seq_len(.N)]
  data.table::setcolorder(history, c("boosting_round",
                                    setdiff(names(history), "boosting_round")))
  data.table::fwrite(history, file.path(QC_DIR, paste0(id, "_cv_history.csv")))
  cat("Selected ", rounds, " rounds (minimum test CV metric ",
      format(scores[[rounds]], digits = 8), "; column ", metric_name,
      ")\n", sep = "")
  rounds
}

fit_duration <- function(id, pyxgb) {
  cat("\n===== JOINT QUANTILE DURATION: ", id, " =====\n", sep="")
  if (is_completed(id)) { cat("Verified files already present; skipping\n"); return(invisible(NULL)) }
  sp <- DURATION_SPECS[[id]]
  if (is.null(sp)) stop("Unknown duration endpoint")
  selected <- get_config(DURATION_SOURCE, id, "outcome")
  y <- suppressWarnings(as.numeric(as.character(DT[[sp$y_var]])))
  ok <- is.finite(y) & if (sp$lower_bound == 0) y >= 0 else y > 0
  if (sum(ok) < 1000L) stop("Unexpectedly small eligible duration cohort: ", id)
  sp$predictors <- CORE
  sp$policy <- "final_core"
  sp$levels <- NULL
  sp$transformation <- "log1p"
  enc <- fit_encoder(DT[ok], CORE)
  X <- encode_features(DT[ok], CORE, enc, dense=TRUE)
  feature_names <- colnames(X)
  if (anyDuplicated(feature_names)) stop("Duplicate encoded features: ", id)
  ylog <- log1p(y[ok])
  dmat <- pyxgb$DMatrix(X, label=ylog, missing=NA_real_)
  params <- reticulate::dict(
    objective="reg:quantileerror", quantile_alpha=c(0.10,0.50,0.90),
    eta=0.05, max_depth=as.integer(selected$max_depth),
    min_child_weight=as.numeric(selected$min_child_weight),
    subsample=as.numeric(selected$subsample),
    colsample_bytree=as.numeric(selected$colsample_bytree),
    lambda=as.numeric(selected$lambda), tree_method="hist", device=DEVICE,
    nthread=as.integer(THREADS), seed=as.integer(SEED + 20000L))
  rounds <- select_rounds_from_python_cv(pyxgb,dmat,params,id)
  cat("All-years internal CV selected boosting rounds: ",rounds,"\n")
  fit <- pyxgb$train(params=params, dtrain=dmat,
                      num_boost_round=as.integer(rounds), verbose_eval=FALSE)
  if (DEVICE == "cuda") {
    fitted_config <- as.character(fit$save_config())
    if (!grepl('"device"[[:space:]]*:[[:space:]]*"cuda', fitted_config)) {
      stop("CUDA was requested but Python XGBoost did not report CUDA for ", id)
    }
  }
  set.seed(SEED + 22L)
  test_idx <- sample.int(nrow(X), min(128L,nrow(X)))
  dtest <- pyxgb$DMatrix(X[test_idx,,drop=FALSE], missing=NA_real_)
  before <- quantile_matrix(fit$predict(dtest),length(test_idx))
  path <- model_file(id)
  fit$save_model(path)
  loaded <- pyxgb$Booster(model_file=path)
  after <- quantile_matrix(loaded$predict(dtest),length(test_idx))
  delta <- max(abs(before-after))
  md <- make_metadata(id,sp,enc,feature_names,rounds,selected,
                      sum(ok),"joint_quantile",DURATION_SOURCE)
  md$source_model_engine_version <- as.character(pyxgb$`__version__`)
  saveRDS(md,meta_file(id),version=3,compress=TRUE)
  record_qc(id,"joint_quantile",sum(ok),rounds,selected,delta,path,meta_file(id))
  cat("Python JSON reload parity OK, maximum difference: ",
      format(delta,scientific=TRUE),"\n")
  rm(X,dmat,dtest,fit,loaded,md);gc(verbose=FALSE)
}

if (STAGE %in% c("all","classification")) {
  for (id in names(CLASS_SPECS)) fit_classification(id)
}
if (STAGE %in% c("all","duration")) {
  # Only initialize Python if at least one duration artifact needs creating.
  todo <- names(DURATION_SPECS)[!vapply(names(DURATION_SPECS), is_completed, logical(1))]
  if (length(todo)) {
    reticulate::use_condaenv(PY_ENV, required=TRUE)
    pyxgb <- reticulate::import("xgboost",convert=TRUE)
    major <- as.integer(strsplit(as.character(pyxgb$`__version__`),"\\.")[[1L]][1L])
    if (!is.finite(major) || major < 2L) stop("Python xgboost >=2 is required for joint quantile regression")
    for (id in todo) fit_duration(id,pyxgb)
  }
}

qc_paths <- list.files(QC_DIR,pattern="_serialization_qc\\.csv$",full.names=TRUE)
if (!length(qc_paths)) stop("No verified model outputs found")
all_qc <- data.table::rbindlist(lapply(qc_paths,data.table::fread),fill=TRUE)
if (anyDuplicated(all_qc$endpoint_id)) stop("Duplicate QC endpoint IDs")
expected_ids <- c(names(CLASS_SPECS),names(DURATION_SPECS))
if (!all(all_qc$endpoint_id %in% expected_ids)) stop("Unknown artifact QC endpoint")
for (id in all_qc$endpoint_id) if (!is_completed(id) && !FORCE) stop("QC check failed: ", id)
all_qc <- all_qc[match(intersect(expected_ids,all_qc$endpoint_id),all_qc$endpoint_id)]
data.table::fwrite(all_qc,file.path(OUT,"artifact_manifest.csv"))
if (STAGE == "all" && !setequal(all_qc$endpoint_id,expected_ids)) {
  stop("All-stage run incomplete: expected 8 exported models, found ",nrow(all_qc))
}

notes <- c(
  "SPINE-TRACT deployment artifacts, full eligible 2020-2024 cohort",
  "Five R XGBoost multiclass/binary models and three Python XGBoost joint quantile models.",
  "Structural configurations: locked from original 2024 forward temporal evaluation.",
  paste0("Boosting rounds: ",NFOLD,"-fold internal all-years cross-validation; not independent validation."),
  "Duration models jointly predict log1p Q10/Q50/Q90, postprocessed using expm1,",
  "feasible lower bound (0 for HLOS; 1 for conditional ICU/vent), then row-wise sorting.",
  "Native numeric missing routing, explicit __UNKNOWN__ categorical, __OTHER__ unseen category.",
  "JSON model parity checked after saving/reloading; QC is not external validation.",
  "RDS encoders are prediction specifications and contain no patient-level data.",
  "DO NOT COPY THE LEGACY ROOT data/model_bundle.rds INTO THE NEW APPLICATION.",
  "The Python joint quantile boosters REQUIRE Python XGBoost >=2 inference runtime.",
  "Individual-level ACS PUF data and prediction caches MUST NOT be committed to GitHub.",
  "Review ACS source-data and derived-model artifact sharing terms before any public upload.",
  "See artifact_manifest.csv and qc/*_cv_history.csv and *_serialization_qc.csv."
)
writeLines(notes,file.path(OUT,"DEPLOYMENT_README.txt"))
cat("\nFinished. Verified models: ",nrow(all_qc)," / 8\n",sep="")
cat("Staged at: ",OUT,"\n",sep="")
print(all_qc[,.(endpoint_id,training_n,selected_rounds,max_abs_logit_or_quantile_difference)])