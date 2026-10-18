# ==============================================================================
# 04r_run_SPINE_TRACT_stage2C_REPAIRED_DX_gpu_v1.R
# SPINE-TRACT Stage 2C-R — repaired-DX cumulative finalist adjudication v1
#
# Forward temporal predictor-family adjudication screen.
#
# PURPOSE
#   Compare prespecified predictor-family variants using the same temporal
#   architecture and standardized XGBoost training settings.
#
# IMPORTANT
#   - This is NOT final hyperparameter tuning.
#   - This is NOT the final manuscript performance analysis.
#   - The goal is to determine whether predictor families add stable,
#     clinically meaningful value across forward temporal evaluations.
#   - Borderline families should be adjudicated clinically and can undergo
#     targeted retuning later; no universal delta-AUROC cutoff is imposed.
#
# TEMPORAL FOLDS
#   Eval 2022: train 2020, tune 2021, refit 2020-2021, evaluate 2022
#   Eval 2023: train 2020-2021, tune 2022, refit 2020-2022, evaluate 2023
#   Eval 2024: train 2020-2022, tune 2023, refit 2020-2023, evaluate 2024
#
# HARDWARE
#   Requires a CUDA-enabled XGBoost R binary and HARD-FAILS if the trained
#   booster is not on cuda:0. CPU/data preparation uses ~80% logical threads.
#
# RUN AFTER
#   02_build_spine_tract_candidate_dataset.R
# ==============================================================================

rm(list = ls())
gc()

required_pkgs <- c("data.table", "xgboost", "Matrix")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ", paste(missing_pkgs, collapse = ", "),
    ". Install before running Stage 2B.",
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(xgboost)
  library(Matrix)
})

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

DATA_DIR <- file.path(PROJECT_DIR, "data", "derived", "stage2_predictor_adjudication")
OUT_DIR <- file.path(PROJECT_DIR, "outputs", "stage2_predictor_adjudication", "stage2C_REPAIRED_DX_v1")
PRED_DIR <- file.path(OUT_DIR, "predictions")
MODEL_DIR <- file.path(OUT_DIR, "models")

for (d in c(OUT_DIR, PRED_DIR, MODEL_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

DATA_FILE <- file.path(DATA_DIR, "spine_tract_frozen_primary_2020_2024.rds")
SPEC_FILE <- file.path(DATA_DIR, "spine_tract_predictor_sets.rds")

if (!file.exists(DATA_FILE)) stop("Frozen Stage 2A dataset not found: ", DATA_FILE, call. = FALSE)
if (!file.exists(SPEC_FILE)) stop("Predictor specification not found: ", SPEC_FILE, call. = FALSE)

TARGET_CPU_FRACTION <- 0.80
LOGICAL_CORES <- parallel::detectCores(logical = TRUE)
if (is.na(LOGICAL_CORES) || LOGICAL_CORES < 1L) LOGICAL_CORES <- 1L
CPU_THREADS <- max(1L, min(LOGICAL_CORES, as.integer(round(LOGICAL_CORES * TARGET_CPU_FRACTION))))
data.table::setDTthreads(CPU_THREADS)

MODEL_SEED <- 20260929L
MAX_ROUNDS <- 3000L
EARLY_STOP <- 100L
PRINT_EVERY <- 100L

# Standardized family-screen parameters, modeled after the final TBI-TRACT
# architecture. Final model tuning occurs only after predictor lock.
ETA <- 0.05
SUBSAMPLE <- 0.85
COLSAMPLE <- 0.85
MIN_CHILD_WEIGHT <- 10

dt <- as.data.table(readRDS(DATA_FILE))
spec <- readRDS(SPEC_FILE)
predictor_sets <- spec$predictor_sets
predictor_groups <- spec$predictor_groups

clinical_backbone <- predictor_sets$clinical_backbone

coarse_full <- unique(c(
  setdiff(
    clinical_backbone,
    predictor_groups$spine_granular
  ),
  predictor_groups$spine_coarse
))

predictor_sets <- list(
  coarse_plus_chest_pupil = unique(c(
    coarse_full,
    predictor_groups$chest_detail,
    predictor_groups$pupil
  )),

  coarse_plus_chest_pupil_no_soft = setdiff(
    unique(c(
      coarse_full,
      predictor_groups$chest_detail,
      predictor_groups$pupil
    )),
    "dx_tl_soft_tissue_any"
  ),

  granular_plus_chest_pupil = unique(c(
    clinical_backbone,
    predictor_groups$chest_detail,
    predictor_groups$pupil
  )),

  coarse_plus_chest_pupil_social = unique(c(
    coarse_full,
    predictor_groups$chest_detail,
    predictor_groups$pupil,
    predictor_groups$social_context
  ))
)

predictor_manifest <- rbindlist(
  lapply(
    names(predictor_sets),
    function(nm) {
      data.table(
        variant = nm,
        n_predictors = length(predictor_sets[[nm]]),
        predictor = predictor_sets[[nm]]
      )
    }
  )
)

fwrite(
  predictor_manifest,
  file.path(OUT_DIR, "00b_STAGE2C_PREDICTOR_MANIFEST.csv")
)

if (nrow(dt) != 424017L) {
  stop("Frozen primary dataset N is not 424,017. Observed: ", nrow(dt), call. = FALSE)
}

repair_expected <- data.table(
  predictor = c(
    "dx_thoracic_vertebral_fracture",
    "dx_lumbar_vertebral_fracture",
    "dx_tl_spinal_cord_injury_any",
    "dx_tl_root_cauda_any",
    "dx_tl_soft_tissue_any",
    "dx_cervical_spine_injury",
    "dx_tbi_any"
  ),
  expected_positive_n = c(
    218002L,
    259812L,
    21661L,
    2215L,
    19988L,
    65376L,
    92755L
  )
)

repair_preflight <- repair_expected[
  ,
  .(
    predictor,
    expected_positive_n,
    observed_positive_n = vapply(
      predictor,
      function(v) sum(dt[[v]] == 1L, na.rm = TRUE),
      integer(1)
    ),
    missing_n = vapply(
      predictor,
      function(v) sum(is.na(dt[[v]])),
      integer(1)
    )
  )
]

repair_preflight[
  ,
  concordant :=
    observed_positive_n == expected_positive_n &
      missing_n == 0L
]

fwrite(
  repair_preflight,
  file.path(OUT_DIR, "00_REPAIRED_DX_PREFLIGHT.csv")
)

if (!all(repair_preflight$concordant)) {
  print(repair_preflight)
  stop(
    "HARD FAIL: repaired diagnosis phenotype preflight failed.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. CUDA hard-fail detection
# ------------------------------------------------------------------------------

config_to_text <- function(model) {
  cfg <- xgb.config(model)
  if (is.character(cfg)) {
    paste(cfg, collapse = "")
  } else if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::toJSON(cfg, auto_unbox = TRUE, null = "null")
  } else {
    paste(capture.output(str(cfg)), collapse = " ")
  }
}

assert_cuda_booster <- function(model, context = "XGBoost booster") {
  txt <- tolower(config_to_text(model))
  ver <- packageVersion("xgboost")

  if (ver < package_version("2.0.0")) {
    gpu_ok <- grepl("grow_gpu_hist|gpu_hist", txt, perl = TRUE)
  } else {
    gpu_ok <- grepl(
      '"device"\\s*:\\s*"cuda(:0)?"',
      txt,
      perl = TRUE
    ) || grepl("grow_gpu_hist", txt, perl = TRUE)
  }

  if (!gpu_ok) {
    stop(
      "\nHARD FAIL: ", context, " is not using the GPU.\n",
      "Detected xgboost version: ", as.character(ver), "\n",
      "For 1.x we require gpu_hist; for 2.x+ we require CUDA device configuration.\n",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

gpu_params_for_version <- function(cpu_threads) {
  ver <- packageVersion("xgboost")

  if (ver < package_version("2.0.0")) {
    list(
      mode = paste0(
        "R CUDA legacy: xgboost ", as.character(ver),
        " | tree_method=gpu_hist | gpu_id=0"
      ),
      extra = list(
        tree_method = "gpu_hist",
        predictor = "gpu_predictor",
        gpu_id = 0L
      )
    )
  } else {
    list(
      mode = paste0(
        "R CUDA modern: xgboost ", as.character(ver),
        " | tree_method=hist | device=cuda:0"
      ),
      extra = list(
        tree_method = "hist",
        device = "cuda:0"
      )
    )
  }
}

detect_xgb_compute <- function(cpu_threads) {
  set.seed(1L)
  x <- matrix(rnorm(20000 * 32), nrow = 20000, ncol = 32)
  y <- rbinom(20000, 1, 0.4)
  d <- xgb.DMatrix(x, label = y, missing = NA_real_)

  compute <- gpu_params_for_version(cpu_threads)

  params <- c(
    list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = 5,
      eta = 0.1,
      nthread = cpu_threads,
      verbosity = 1
    ),
    compute$extra
  )

  probe <- tryCatch(
    xgb.train(
      params = params,
      data = d,
      nrounds = 20,
      verbose = 0
    ),
    error = function(e) {
      stop(
        "\nCUDA XGBoost R smoke test failed:\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )

  assert_cuda_booster(probe, "CUDA smoke-test booster")
  cfg_txt <- config_to_text(probe)

  rm(probe, d, x, y)
  invisible(gc())

  compute$config <- cfg_txt
  compute
}

compute <- detect_xgb_compute(CPU_THREADS)

gpu_query <- tryCatch(
  system2(
    "nvidia-smi",
    c("--query-gpu=name,memory.total,driver_version,utilization.gpu",
      "--format=csv,noheader"),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) paste("nvidia-smi unavailable:", conditionMessage(e))
)

hardware <- data.table(
  logical_cores = LOGICAL_CORES,
  cpu_threads = CPU_THREADS,
  xgboost_version = as.character(packageVersion("xgboost")),
  compute_mode = compute$mode,
  gpu_query = paste(gpu_query, collapse = " | ")
)
fwrite(hardware, file.path(OUT_DIR, "00_compute_audit.csv"))

cat("\n============================================================\n")
cat("SPINE-TRACT STAGE 2B: R-NATIVE CUDA PREDICTOR ADJUDICATION\n")
cat("============================================================\n")
cat("N:            ", format(nrow(dt), big.mark = ","), "\n", sep = "")
cat("CPU threads:  ", CPU_THREADS, "/", LOGICAL_CORES, "\n", sep = "")
cat("XGBoost:      ", as.character(packageVersion("xgboost")), "\n", sep = "")
cat("Compute:      ", compute$mode, "\n", sep = "")
cat("GPU:          ", paste(gpu_query, collapse = " | "), "\n", sep = "")
cat("Max rounds:   ", MAX_ROUNDS, "\n", sep = "")
cat("Early stop:   ", EARLY_STOP, "\n", sep = "")

# ------------------------------------------------------------------------------
# 2. Metrics
# ------------------------------------------------------------------------------

clamp_prob <- function(p, eps = 1e-7) pmin(pmax(p, eps), 1 - eps)

fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  r <- rank(p, ties.method = "average")
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fast_auprc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])
  npos <- sum(y == 1L)
  if (npos == 0L) return(NA_real_)
  ord <- order(p, decreasing = TRUE)
  ys <- y[ord]
  precision <- cumsum(ys == 1L) / seq_along(ys)
  mean(precision[ys == 1L])
}

binary_brier <- function(y, p) mean((y - p)^2)
binary_logloss <- function(y, p) {
  pp <- clamp_prob(p)
  -mean(y * log(pp) + (1 - y) * log(1 - pp))
}

calibration_is <- function(y, p) {
  pp <- clamp_prob(p)
  lp <- qlogis(pp)

  fit_slope <- tryCatch(
    glm(y ~ lp, family = binomial()),
    error = function(e) NULL
  )
  slope <- if (is.null(fit_slope)) NA_real_ else unname(coef(fit_slope)[2])

  fit_int <- tryCatch(
    glm(y ~ offset(lp), family = binomial()),
    error = function(e) NULL
  )
  intercept <- if (is.null(fit_int)) NA_real_ else unname(coef(fit_int)[1])

  c(intercept = intercept, slope = slope)
}

multiclass_metrics <- function(y0, prob, levels, endpoint, variant, fold_year, best_iter) {
  k <- length(levels)
  if (ncol(prob) != k) stop("Prediction matrix class mismatch.", call. = FALSE)

  valid <- !is.na(y0)
  y <- y0[valid]
  p <- prob[valid, , drop = FALSE]

  true_p <- p[cbind(seq_along(y), y + 1L)]
  mlogloss <- -mean(log(pmax(true_p, 1e-7)))

  onehot <- matrix(0, nrow = length(y), ncol = k)
  onehot[cbind(seq_along(y), y + 1L)] <- 1
  brier <- mean(rowSums((onehot - p)^2))

  class_rows <- rbindlist(lapply(seq_len(k), function(j) {
    yy <- as.integer(y == (j - 1L))
    pp <- p[, j]
    cal <- calibration_is(yy, pp)

    data.table(
      endpoint = endpoint,
      variant = variant,
      evaluation_year = fold_year,
      metric_level = "class",
      class = levels[j],
      n = length(yy),
      prevalence = mean(yy),
      auroc = fast_auc(yy, pp),
      auprc = fast_auprc(yy, pp),
      pr_lift = if (mean(yy) > 0) fast_auprc(yy, pp) / mean(yy) else NA_real_,
      brier = binary_brier(yy, pp),
      logloss = binary_logloss(yy, pp),
      calibration_intercept = unname(cal["intercept"]),
      calibration_slope = unname(cal["slope"]),
      multiclass_logloss = mlogloss,
      multiclass_brier = brier,
      best_iteration = best_iter
    )
  }))

  macro <- data.table(
    endpoint = endpoint,
    variant = variant,
    evaluation_year = fold_year,
    metric_level = "macro",
    class = "MACRO",
    n = max(class_rows$n),
    prevalence = NA_real_,
    auroc = mean(class_rows$auroc, na.rm = TRUE),
    auprc = mean(class_rows$auprc, na.rm = TRUE),
    pr_lift = mean(class_rows$pr_lift, na.rm = TRUE),
    brier = mean(class_rows$brier, na.rm = TRUE),
    logloss = mean(class_rows$logloss, na.rm = TRUE),
    calibration_intercept =
      median(class_rows$calibration_intercept, na.rm = TRUE),
    calibration_slope =
      median(class_rows$calibration_slope, na.rm = TRUE),
    multiclass_logloss =
      unique(class_rows$multiclass_logloss)[1],
    multiclass_brier =
      unique(class_rows$multiclass_brier)[1],
    best_iteration = best_iter
  )

  rbind(class_rows, macro, fill = TRUE)
}

binary_metrics <- function(y, p, endpoint, variant, fold_year, best_iter) {
  cal <- calibration_is(y, p)
  prev <- mean(y == 1L)

  data.table(
    endpoint = endpoint,
    variant = variant,
    evaluation_year = fold_year,
    metric_level = "binary",
    class = "Event",
    n = length(y),
    prevalence = prev,
    auroc = fast_auc(y, p),
    auprc = fast_auprc(y, p),
    pr_lift = if (prev > 0) fast_auprc(y, p) / prev else NA_real_,
    brier = binary_brier(y, p),
    logloss = binary_logloss(y, p),
    calibration_intercept = unname(cal["intercept"]),
    calibration_slope = unname(cal["slope"]),
    multiclass_logloss = NA_real_,
    multiclass_brier = NA_real_,
    best_iteration = best_iter
  )
}

# ------------------------------------------------------------------------------
# 3. Fold-specific encoder
# ------------------------------------------------------------------------------

is_numeric_predictor <- function(x) is.numeric(x) || is.integer(x)

fit_encoder <- function(train_dt, predictors) {
  out <- list()
  for (v in predictors) {
    x <- train_dt[[v]]
    if (is_numeric_predictor(x)) {
      out[[v]] <- list(type = "numeric")
    } else {
      z <- as.character(x)
      z[is.na(z) | z == ""] <- "__UNKNOWN__"
      lev <- sort(unique(z))
      lev <- unique(c(lev, "__UNKNOWN__", "__OTHER__"))
      out[[v]] <- list(type = "categorical", levels = lev)
    }
  }
  out
}

encode_with_spec <- function(d, predictors, enc) {
  blocks <- vector("list", length(predictors))
  block_names <- character(length(predictors))

  for (i in seq_along(predictors)) {
    v <- predictors[i]
    info <- enc[[v]]
    x <- d[[v]]

    if (info$type == "numeric") {
      m <- Matrix(as.numeric(x), ncol = 1, sparse = TRUE)
      colnames(m) <- v
      blocks[[i]] <- m
      block_names[i] <- v
    } else {
      z <- as.character(x)
      z[is.na(z) | z == ""] <- "__UNKNOWN__"
      z[!z %in% info$levels] <- "__OTHER__"
      f <- factor(z, levels = info$levels)

      ii <- seq_along(f)
      jj <- as.integer(f)
      m <- sparseMatrix(
        i = ii,
        j = jj,
        x = 1,
        dims = c(length(f), length(info$levels))
      )
      colnames(m) <- paste0(v, "__", make.names(info$levels, unique = TRUE))
      blocks[[i]] <- m
      block_names[i] <- v
    }
  }

  X <- do.call(cbind, blocks)
  X
}

# ------------------------------------------------------------------------------
# 4. Model configuration
# ------------------------------------------------------------------------------

base_multiclass_params <- function(k) {
  c(
    list(
      objective = "multi:softprob",
      eval_metric = "mlogloss",
      num_class = k,
      max_depth = 4L,
      eta = ETA,
      subsample = SUBSAMPLE,
      colsample_bytree = COLSAMPLE,
      min_child_weight = MIN_CHILD_WEIGHT,
      lambda = 1,
      alpha = 0,
      nthread = CPU_THREADS
    ),
    compute$extra
  )
}

base_binary_params <- function() {
  c(
    list(
      objective = "binary:logistic",
      eval_metric = "logloss",
      max_depth = 3L,
      eta = ETA,
      subsample = SUBSAMPLE,
      colsample_bytree = COLSAMPLE,
      min_child_weight = MIN_CHILD_WEIGHT,
      lambda = 1,
      alpha = 0,
      nthread = CPU_THREADS
    ),
    compute$extra
  )
}

extract_best_iteration <- function(model, fallback) {
  b <- tryCatch(model$best_iteration, error = function(e) NULL)
  if (!is.null(b) && length(b) == 1L && is.finite(b)) return(as.integer(b))
  a <- tryCatch(xgb.attr(model, "best_iteration"), error = function(e) NULL)
  if (!is.null(a) && length(a) == 1L && !is.na(a)) return(as.integer(a) + 1L)
  as.integer(fallback)
}

folds <- list(
  `2022` = list(train = 2020L, tune = 2021L, dev = 2020:2021, eval = 2022L),
  `2023` = list(train = 2020:2021, tune = 2022L, dev = 2020:2022, eval = 2023L),
  `2024` = list(train = 2020:2022, tune = 2023L, dev = 2020:2023, eval = 2024L)
)

endpoint_specs <- list(
  discharge4 = list(
    type = "multiclass",
    levels = c(
      "Home/home health",
      "Inpatient rehabilitation",
      "Other post-acute facility",
      "Death/hospice"
    )
  ),
  hlos3 = list(
    type = "multiclass",
    levels = c("<=7 days", "8-19 days", ">=20 days")
  ),
  icu3 = list(
    type = "multiclass",
    levels = c("No ICU", "ICU 1-7 days", "ICU >=8 days")
  ),
  vent3 = list(
    type = "multiclass",
    levels = c(
      "No ventilation",
      "Ventilation 1-7 days",
      "Ventilation >=8 days"
    )
  ),
  operative_spine_intervention = list(
    type = "binary",
    levels = c("0", "1")
  )
)

# Endpoint-specific allowable variants.
# Social context is only tested for disposition/HLOS.
# Airway state is deliberately NOT tested for ventilation and is restricted
# to distal/non-airway endpoints because it can reveal support already underway.
variant_endpoints <- list(
  coarse_plus_chest_pupil = names(endpoint_specs),
  coarse_plus_chest_pupil_no_soft = names(endpoint_specs),
  granular_plus_chest_pupil = names(endpoint_specs),
  coarse_plus_chest_pupil_social = c("discharge4", "hlos3")
)

# ------------------------------------------------------------------------------
# 5. Fit one temporal fold
# ------------------------------------------------------------------------------

fit_one <- function(endpoint, variant, fold) {
  endpoint_spec <- endpoint_specs[[endpoint]]
  predictors <- predictor_sets[[variant]]

  train_idx <- which(dt$admission_year %in% fold$train)
  tune_idx <- which(dt$admission_year %in% fold$tune)
  dev_idx <- which(dt$admission_year %in% fold$dev)
  eval_idx <- which(dt$admission_year %in% fold$eval)

  # Outcome-specific complete outcome only; predictor missingness is retained.
  if (endpoint_spec$type == "multiclass") {
    y_factor <- factor(as.character(dt[[endpoint]]), levels = endpoint_spec$levels)
    y_all <- as.integer(y_factor) - 1L
    train_idx <- train_idx[!is.na(y_all[train_idx])]
    tune_idx <- tune_idx[!is.na(y_all[tune_idx])]
    dev_idx <- dev_idx[!is.na(y_all[dev_idx])]
    eval_idx <- eval_idx[!is.na(y_all[eval_idx])]
  } else {
    y_all <- as.integer(dt[[endpoint]])
    train_idx <- train_idx[!is.na(y_all[train_idx]) & y_all[train_idx] %in% c(0L, 1L)]
    tune_idx <- tune_idx[!is.na(y_all[tune_idx]) & y_all[tune_idx] %in% c(0L, 1L)]
    dev_idx <- dev_idx[!is.na(y_all[dev_idx]) & y_all[dev_idx] %in% c(0L, 1L)]
    eval_idx <- eval_idx[!is.na(y_all[eval_idx]) & y_all[eval_idx] %in% c(0L, 1L)]
  }

  if (length(train_idx) == 0L || length(tune_idx) == 0L ||
      length(dev_idx) == 0L || length(eval_idx) == 0L) {
    stop("Empty temporal partition for ", endpoint, " / ", variant, call. = FALSE)
  }

  # Fit encoder only from training-period data.
  enc <- fit_encoder(dt[train_idx], predictors)

  X_train <- encode_with_spec(dt[train_idx], predictors, enc)
  X_tune <- encode_with_spec(dt[tune_idx], predictors, enc)
  X_dev <- encode_with_spec(dt[dev_idx], predictors, enc)
  X_eval <- encode_with_spec(dt[eval_idx], predictors, enc)

  if (!identical(colnames(X_train), colnames(X_tune)) ||
      !identical(colnames(X_train), colnames(X_dev)) ||
      !identical(colnames(X_train), colnames(X_eval))) {
    stop("Encoder produced inconsistent feature columns.", call. = FALSE)
  }

  dtrain <- xgb.DMatrix(X_train, label = y_all[train_idx], missing = NA_real_)
  dtune <- xgb.DMatrix(X_tune, label = y_all[tune_idx], missing = NA_real_)

  params <- if (endpoint_spec$type == "multiclass") {
    base_multiclass_params(length(endpoint_spec$levels))
  } else {
    base_binary_params()
  }

  set.seed(MODEL_SEED)

  # R xgboost API compatibility:
  #   xgboost 1.x uses watchlist=
  #   newer R builds accept evals=
  # The historical Windows GPU build installed for this workstation is 1.7.2,
  # so using evals= would be swallowed into params and early stopping would fail.
  if (packageVersion("xgboost") < package_version("2.0.0")) {
    es <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = MAX_ROUNDS,
      watchlist = list(train = dtrain, tune = dtune),
      early_stopping_rounds = EARLY_STOP,
      maximize = FALSE,
      verbose = 1,
      print_every_n = PRINT_EVERY
    )
  } else {
    es <- xgb.train(
      params = params,
      data = dtrain,
      nrounds = MAX_ROUNDS,
      evals = list(train = dtrain, tune = dtune),
      early_stopping_rounds = EARLY_STOP,
      maximize = FALSE,
      verbose = 1,
      print_every_n = PRINT_EVERY
    )
  }

  assert_cuda_booster(
    es,
    paste0(endpoint, " / ", variant, " / eval ", fold$eval, " early-stopping fit")
  )

  best_iter <- extract_best_iteration(es, MAX_ROUNDS)

  hit_ceiling <- best_iter >= (MAX_ROUNDS - 5L)
  if (hit_ceiling) {
    warning(
      endpoint, " / ", variant, " / eval ", fold$eval,
      " selected the boosting ceiling (", best_iter, "/", MAX_ROUNDS, ")."
    )
  }

  rm(dtrain, dtune, X_train, X_tune, es)
  gc()

  ddev <- xgb.DMatrix(X_dev, label = y_all[dev_idx], missing = NA_real_)
  deval <- xgb.DMatrix(X_eval, label = y_all[eval_idx], missing = NA_real_)

  set.seed(MODEL_SEED)
  fit <- xgb.train(
    params = params,
    data = ddev,
    nrounds = best_iter,
    verbose = 0
  )

  assert_cuda_booster(
    fit,
    paste0(endpoint, " / ", variant, " / eval ", fold$eval, " refit")
  )

  raw <- predict(fit, deval)

  if (endpoint_spec$type == "multiclass") {
    if (is.matrix(raw) || is.data.frame(raw)) {
      prob <- as.matrix(raw)
    } else {
      prob <- matrix(
        as.numeric(raw),
        ncol = length(endpoint_spec$levels),
        byrow = TRUE
      )
    }

    met <- multiclass_metrics(
      y_all[eval_idx],
      prob,
      endpoint_spec$levels,
      endpoint,
      variant,
      fold$eval,
      best_iter
    )

    pred_dt <- data.table(
      admission_year = dt$admission_year[eval_idx],
      inc_key = dt$inc_key[eval_idx],
      endpoint = endpoint,
      variant = variant,
      y = y_all[eval_idx]
    )
    for (j in seq_along(endpoint_spec$levels)) {
      pred_dt[, paste0("p_", make.names(endpoint_spec$levels[j])) := prob[, j]]
    }
  } else {
    p <- as.numeric(raw)
    met <- binary_metrics(
      y_all[eval_idx],
      p,
      endpoint,
      variant,
      fold$eval,
      best_iter
    )
    pred_dt <- data.table(
      admission_year = dt$admission_year[eval_idx],
      inc_key = dt$inc_key[eval_idx],
      endpoint = endpoint,
      variant = variant,
      y = y_all[eval_idx],
      pred = p
    )
  }

  model_slug <- paste(endpoint, variant, fold$eval, sep = "__")

  fwrite(
    pred_dt,
    file.path(PRED_DIR, paste0(model_slug, "_predictions.csv"))
  )

  saveRDS(
    list(
      model = fit,
      endpoint = endpoint,
      variant = variant,
      evaluation_year = fold$eval,
      predictors = predictors,
      feature_names = colnames(X_dev),
      encoder = enc,
      best_iteration = best_iter,
      compute_mode = compute$mode
    ),
    file.path(MODEL_DIR, paste0(model_slug, "_model.rds")),
    compress = FALSE
  )

  qc <- data.table(
    endpoint = endpoint,
    variant = variant,
    evaluation_year = fold$eval,
    train_years = paste(fold$train, collapse = ","),
    tune_years = paste(fold$tune, collapse = ","),
    dev_refit_years = paste(fold$dev, collapse = ","),
    n_train = length(train_idx),
    n_tune = length(tune_idx),
    n_dev_refit = length(dev_idx),
    n_eval = length(eval_idx),
    n_predictors = length(predictors),
    n_encoded_features = ncol(X_dev),
    best_iteration = best_iter,
    hit_round_ceiling = hit_ceiling
  )

  rm(ddev, deval, X_dev, X_eval, fit, raw)
  gc()

  list(metrics = met, qc = qc)
}

# ------------------------------------------------------------------------------
# 6. Run adjudication screen
# ------------------------------------------------------------------------------

run_plan <- rbindlist(lapply(names(variant_endpoints), function(v) {
  data.table(
    variant = v,
    endpoint = variant_endpoints[[v]]
  )
}))

# Reference first so delta calculations always have a baseline.
run_plan[, variant_order := match(variant, names(variant_endpoints))]
setorder(run_plan, variant_order, endpoint)
fwrite(run_plan[, !"variant_order"], file.path(OUT_DIR, "01_run_plan.csv"))

all_metrics <- list()
all_qc <- list()
runtime_rows <- list()
counter <- 0L
total_runs <- nrow(run_plan) * length(folds)

for (i in seq_len(nrow(run_plan))) {
  variant <- run_plan$variant[i]
  endpoint <- run_plan$endpoint[i]

  for (fy in names(folds)) {
    counter <- counter + 1L
    cat("\n============================================================\n")
    cat("[", counter, "/", total_runs, "] ", endpoint, " | ", variant,
        " | eval ", fy, "\n", sep = "")
    cat("============================================================\n")

    t0 <- Sys.time()
    res <- fit_one(endpoint, variant, folds[[fy]])
    elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

    all_metrics[[length(all_metrics) + 1L]] <- res$metrics
    all_qc[[length(all_qc) + 1L]] <- res$qc
    runtime_rows[[length(runtime_rows) + 1L]] <- data.table(
      endpoint = endpoint,
      variant = variant,
      evaluation_year = as.integer(fy),
      runtime_minutes = elapsed
    )

    # Incremental checkpoints so a long GPU run is recoverable.
    fwrite(
      rbindlist(all_metrics, fill = TRUE),
      file.path(OUT_DIR, "02_metrics_checkpoint.csv")
    )
    fwrite(
      rbindlist(all_qc, fill = TRUE),
      file.path(OUT_DIR, "03_training_QC_checkpoint.csv")
    )
    fwrite(
      rbindlist(runtime_rows, fill = TRUE),
      file.path(OUT_DIR, "04_runtime_checkpoint.csv")
    )
  }
}

metrics <- rbindlist(all_metrics, fill = TRUE)
training_qc <- rbindlist(all_qc, fill = TRUE)
runtime <- rbindlist(runtime_rows, fill = TRUE)

fwrite(metrics, file.path(OUT_DIR, "05_predictor_adjudication_metrics.csv"))
fwrite(training_qc, file.path(OUT_DIR, "06_predictor_adjudication_training_QC.csv"))
fwrite(runtime, file.path(OUT_DIR, "07_predictor_adjudication_runtime.csv"))

# ------------------------------------------------------------------------------
# 7. Delta summaries against pragmatic coarse finalist
# ------------------------------------------------------------------------------

headline <- metrics[
  metric_level %in% c("macro", "binary")
]

headline_key_qc <- headline[
  ,
  .N,
  by = .(
    endpoint,
    variant,
    evaluation_year
  )
]

if (
  nrow(headline_key_qc) != 51L ||
    any(headline_key_qc$N != 1L)
) {
  fwrite(
    headline_key_qc,
    file.path(OUT_DIR, "08a_HEADLINE_KEY_QC_FAILED.csv")
  )
  stop(
    "Expected exactly 51 unique headline rows.",
    call. = FALSE
  )
}

reference <- headline[
  variant == "coarse_plus_chest_pupil",
  .(
    endpoint,
    evaluation_year,
    ref_auroc = auroc,
    ref_auprc = auprc,
    ref_pr_lift = pr_lift,
    ref_brier = brier,
    ref_logloss = fifelse(
      !is.na(multiclass_logloss),
      multiclass_logloss,
      logloss
    ),
    ref_cal_intercept = calibration_intercept,
    ref_cal_slope = calibration_slope
  )
]

reference_qc <- reference[
  ,
  .N,
  by = .(
    endpoint,
    evaluation_year
  )
]

if (
  nrow(reference_qc) != 15L ||
    any(reference_qc$N != 1L)
) {
  stop(
    "Pragmatic coarse reference is not unique by endpoint/year.",
    call. = FALSE
  )
}

compare <- merge(
  headline,
  reference,
  by = c(
    "endpoint",
    "evaluation_year"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(compare) != nrow(headline)) {
  stop(
    "Delta merge changed row count.",
    call. = FALSE
  )
}

compare[
  ,
  observed_logloss := fifelse(
    !is.na(multiclass_logloss),
    multiclass_logloss,
    logloss
  )
]

compare[, `:=`(
  delta_auroc = auroc - ref_auroc,
  delta_auprc = auprc - ref_auprc,
  delta_pr_lift = pr_lift - ref_pr_lift,
  delta_brier = brier - ref_brier,
  delta_logloss = observed_logloss - ref_logloss,
  delta_abs_calibration_intercept =
    abs(calibration_intercept) -
      abs(ref_cal_intercept),
  delta_abs_calibration_slope_error =
    abs(calibration_slope - 1) -
      abs(ref_cal_slope - 1)
)]

fwrite(
  compare,
  file.path(
    OUT_DIR,
    "08_foldwise_deltas_vs_coarse_plus_chest_pupil.csv"
  )
)

summary <- compare[
  variant != "coarse_plus_chest_pupil",
  .(
    n_temporal_evaluations = .N,
    median_delta_auroc =
      median(delta_auroc, na.rm = TRUE),
    worst_delta_auroc =
      min(delta_auroc, na.rm = TRUE),
    median_delta_auprc =
      median(delta_auprc, na.rm = TRUE),
    worst_delta_auprc =
      min(delta_auprc, na.rm = TRUE),
    median_delta_brier =
      median(delta_brier, na.rm = TRUE),
    median_delta_logloss =
      median(delta_logloss, na.rm = TRUE),
    median_delta_abs_calibration_intercept =
      median(
        delta_abs_calibration_intercept,
        na.rm = TRUE
      ),
    median_delta_abs_calibration_slope_error =
      median(
        delta_abs_calibration_slope_error,
        na.rm = TRUE
      )
  ),
  by = .(
    variant,
    endpoint
  )
]

fwrite(
  summary,
  file.path(
    OUT_DIR,
    "09_stage2C_finalist_delta_summary.csv"
  )
)

# Decision-oriented table: granular-vs-coarse, soft-tissue ablation, context.
decision_table <- dcast(
  headline[
    ,
    .(
      endpoint,
      evaluation_year,
      variant,
      auroc,
      auprc,
      brier
    )
  ],
  endpoint + evaluation_year ~ variant,
  value.var = c(
    "auroc",
    "auprc",
    "brier"
  )
)

fwrite(
  decision_table,
  file.path(
    OUT_DIR,
    "10_stage2C_decision_table_wide.csv"
  )
)

# Compact predictor counts.
predictor_counts <- unique(
  predictor_manifest[
    ,
    .(
      variant,
      n_predictors
    )
  ]
)

fwrite(
  predictor_counts,
  file.path(
    OUT_DIR,
    "11_stage2C_predictor_counts.csv"
  )
)

# Class-specific comparison to pragmatic coarse finalist.
class_ref <- metrics[
  variant == "coarse_plus_chest_pupil" &
    metric_level == "class",
  .(
    endpoint,
    evaluation_year,
    class,
    ref_class_auroc = auroc,
    ref_class_auprc = auprc,
    ref_class_brier = brier
  )
]

class_compare <- merge(
  metrics[
    metric_level == "class"
  ],
  class_ref,
  by = c(
    "endpoint",
    "evaluation_year",
    "class"
  ),
  all.x = TRUE,
  sort = FALSE
)

class_compare[, `:=`(
  delta_class_auroc =
    auroc - ref_class_auroc,
  delta_class_auprc =
    auprc - ref_class_auprc,
  delta_class_brier =
    brier - ref_class_brier
)]

fwrite(
  class_compare,
  file.path(
    OUT_DIR,
    "12_stage2C_class_specific_deltas.csv"
  )
)

# ------------------------------------------------------------------------------
# 8. Final README
# ------------------------------------------------------------------------------

readme <- c(
  "SPINE-TRACT Stage 2C-R repaired-DX cumulative finalist adjudication complete.",
  "",
  paste0("Compute mode: ", compute$mode),
  paste0("CPU threads: ", CPU_THREADS, "/", LOGICAL_CORES),
  "",
  "This stage contains 51 GPU fits:",
  "- coarse + chest + pupil: 5 endpoints x 3 folds = 15",
  "- same minus TL soft tissue: 5 x 3 = 15",
  "- granular + chest + pupil: 5 x 3 = 15",
  "- coarse + chest + pupil + social context: discharge/HLOS x 3 = 6",
  "",
  "Decision questions:",
  "1. Does granular spine representation materially outperform coarse after cumulative augmentation?",
  "2. Is TL soft-tissue information still necessary after chest+pupil augmentation?",
  "3. Does social-context augmentation remain useful for discharge/HLOS?",
  "",
  "Review:",
  "09_stage2C_finalist_delta_summary.csv",
  "10_stage2C_decision_table_wide.csv",
  "11_stage2C_predictor_counts.csv",
  "12_stage2C_class_specific_deltas.csv"
)

writeLines(
  readme,
  file.path(
    OUT_DIR,
    "README_STAGE2C_REPAIRED_DX.txt"
  )
)

cat("\n============================================================\n")
cat("STAGE 2C-R REPAIRED-DX COMPLETE\n")
cat("============================================================\n")
cat("Expected completed fits: 51\n")
cat("Headline rows: ", nrow(headline), "\n", sep = "")
cat("Review:\n  ",
    file.path(OUT_DIR, "09_stage2C_finalist_delta_summary.csv"),
    "\n  ",
    file.path(OUT_DIR, "10_stage2C_decision_table_wide.csv"),
    "\n",
    sep = "")
