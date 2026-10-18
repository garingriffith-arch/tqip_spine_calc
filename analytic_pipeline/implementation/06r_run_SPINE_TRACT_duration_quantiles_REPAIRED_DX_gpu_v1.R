# ==============================================================================
# 06r_run_SPINE_TRACT_duration_quantiles_REPAIRED_DX_gpu_v1.R
# SPINE-TRACT Stage 3C-R — repaired-DX Q10/Q50/Q90 duration models
#
# PURPOSE
#   Refit the three final continuous-duration models after repair of the
#   diagnosis phenotypes and lock of the corrected 54-predictor final_core.
#
# OUTCOMES
#   hospital_los                 : unconditional, observed HLOS
#   icu_los_conditional          : among ICU users (ICU LOS > 0)
#   ventilator_days_conditional  : among ventilated patients (vent days > 0)
#
# MODEL
#   XGBoost quantile regression, Q10/Q50/Q90 jointly
#   objective = reg:quantileerror
#   outcome transformation = log1p(days)
#   eta = 0.05
#   exact same prespecified 24-config structural grid as repaired Stage 3A
#
# TEMPORAL DESIGN
#   train 2020      | tune 2021 | refit 2020-21 | eval 2022
#   train 2020-21   | tune 2022 | refit 2020-22 | eval 2023
#   train 2020-22   | tune 2023 | refit 2020-23 | eval 2024
#
# DISPLAY / EVALUATION
#   Predictions are back-transformed with expm1(), truncated to the feasible
#   lower support, and row-wise rearranged to guarantee q10 <= q50 <= q90.
#
# IMPORTANT
#   Requires a Python XGBoost build supporting reg:quantileerror (XGBoost >=2).
#   The script uses the existing conda environment "tbi-tract-xgb-gpu" by
#   default; override with SPINE_TRACT_PY_ENV if needed.
# ==============================================================================

rm(list = ls())
gc()

required <- c("data.table", "reticulate")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0L) {
  stop(
    "Missing package(s): ",
    paste(missing, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(reticulate)
})

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

DATA_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication",
  "spine_tract_frozen_primary_2020_2024.rds"
)

SPEC_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication",
  "spine_tract_predictor_sets.rds"
)

CLASS_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification_REPAIRED_DX",
  "rolling_origin_gpu"
)

GRID_FILE <- file.path(
  CLASS_DIR,
  "02_PRESPECIFIED_SEARCH_GRID.csv"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_duration_quantiles_REPAIRED_DX",
  "rolling_origin_gpu"
)

GRID_DIR <- file.path(OUT_DIR, "grid_search")
PRED_DIR <- file.path(OUT_DIR, "predictions")
MODEL_DIR <- file.path(OUT_DIR, "models")

for (d in c(OUT_DIR, GRID_DIR, PRED_DIR, MODEL_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

needed <- c(DATA_FILE, SPEC_FILE, GRID_FILE)
if (!all(file.exists(needed))) {
  stop(
    "Missing prerequisite file(s):\n",
    paste(needed[!file.exists(needed)], collapse = "\n"),
    call. = FALSE
  )
}

MODEL_SEED <- 20261002L
ETA <- 0.05
MAX_ROUNDS <- 3000L
EARLY_STOP <- 100L

logical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(logical_cores) || logical_cores < 2L) logical_cores <- 2L

CPU_THREADS <- max(
  1L,
  min(
    logical_cores - 1L,
    as.integer(round(0.80 * logical_cores))
  )
)

setDTthreads(CPU_THREADS)

# ------------------------------------------------------------------------------
# 1. Python XGBoost CUDA preflight
# ------------------------------------------------------------------------------

PY_ENV <- Sys.getenv(
  "SPINE_TRACT_PY_ENV",
  unset = "tbi-tract-xgb-gpu"
)

reticulate::use_condaenv(
  PY_ENV,
  required = TRUE
)

py_cfg <- reticulate::py_config()

xgb <- reticulate::import(
  "xgboost",
  convert = TRUE
)

xgb_version <- as.character(xgb$`__version__`)

major_version <- suppressWarnings(
  as.integer(
    strsplit(xgb_version, "\\.")[[1L]][1L]
  )
)

if (is.na(major_version) || major_version < 2L) {
  stop(
    "Python XGBoost >=2 is required for reg:quantileerror. Observed version: ",
    xgb_version,
    call. = FALSE
  )
}

probe_x <- matrix(
  rnorm(3000),
  nrow = 300L,
  ncol = 10L
)
probe_y <- abs(rnorm(300L))

probe_d <- xgb$DMatrix(
  probe_x,
  label = probe_y
)

probe_fit <- xgb$train(
  params = reticulate::dict(
    objective = "reg:quantileerror",
    quantile_alpha = c(0.10, 0.50, 0.90),
    tree_method = "hist",
    device = "cuda",
    eta = 0.1,
    max_depth = 2L,
    nthread = as.integer(CPU_THREADS),
    seed = MODEL_SEED
  ),
  dtrain = probe_d,
  num_boost_round = 3L,
  verbose_eval = FALSE
)

probe_cfg <- as.character(probe_fit$save_config())

if (!grepl(
  '"device"[[:space:]]*:[[:space:]]*"cuda',
  probe_cfg,
  ignore.case = TRUE
)) {
  stop(
    "CUDA was not verified in the Python XGBoost quantile probe.",
    call. = FALSE
  )
}

fwrite(
  data.table(
    python = as.character(py_cfg$python)[1L],
    python_xgboost_version = xgb_version,
    device = "cuda",
    logical_cpu_threads = logical_cores,
    allocated_host_threads = CPU_THREADS
  ),
  file.path(OUT_DIR, "00_COMPUTE_ENVIRONMENT.csv")
)

rm(probe_x, probe_y, probe_d, probe_fit)
gc()

# ------------------------------------------------------------------------------
# 2. Data / repaired predictor lock
# ------------------------------------------------------------------------------

dt <- as.data.table(readRDS(DATA_FILE))
spec <- readRDS(SPEC_FILE)
search_grid <- fread(GRID_FILE)

dt[, inc_key := as.character(inc_key)]
dt[, admission_year := as.integer(admission_year)]

if (nrow(dt) != 424017L) {
  stop(
    "Expected repaired primary N=424017; observed ",
    nrow(dt),
    ".",
    call. = FALSE
  )
}

required_grid_cols <- c(
  "config_id",
  "max_depth",
  "min_child_weight",
  "subsample",
  "colsample_bytree",
  "lambda"
)

if (
  nrow(search_grid) != 24L ||
    !all(required_grid_cols %in% names(search_grid))
) {
  stop(
    "Repaired Stage 3A prespecified grid is missing or malformed.",
    call. = FALSE
  )
}

search_grid <- search_grid[, ..required_grid_cols]

pg <- spec$predictor_groups
clinical_backbone <- spec$predictor_sets$clinical_backbone

final_core <- unique(c(
  setdiff(
    clinical_backbone,
    pg$spine_granular
  ),
  pg$spine_coarse,
  pg$chest_detail,
  pg$pupil
))

if (length(final_core) != 54L) {
  stop(
    "Expected corrected final_core length 54; observed ",
    length(final_core),
    ".",
    call. = FALSE
  )
}

if (!"dx_tl_soft_tissue_any" %in% final_core) {
  stop(
    "Corrected final_core is missing dx_tl_soft_tissue_any.",
    call. = FALSE
  )
}

missing_predictors <- setdiff(final_core, names(dt))
if (length(missing_predictors) > 0L) {
  stop(
    "Frozen dataset lacks final_core predictor(s): ",
    paste(missing_predictors, collapse = ", "),
    call. = FALSE
  )
}

fwrite(
  data.table(
    predictor_policy = "final_core",
    n_predictors = length(final_core),
    predictor = final_core
  ),
  file.path(OUT_DIR, "01_FINAL_DURATION_PREDICTOR_MANIFEST.csv")
)

# ------------------------------------------------------------------------------
# 3. Encoder
# ------------------------------------------------------------------------------

is_num <- function(x) {
  is.numeric(x) || is.integer(x)
}

clean_cat <- function(x) {
  z <- trimws(as.character(x))
  z[
    is.na(z) |
      z == "" |
      z %in% c("NA", "NaN", "<NA>")
  ] <- "__UNKNOWN__"
  z
}

fit_encoder <- function(d, predictors) {
  enc <- vector("list", length(predictors))
  names(enc) <- predictors

  for (v in predictors) {
    x <- d[[v]]

    if (is_num(x)) {
      enc[[v]] <- list(type = "numeric")
    } else {
      z <- clean_cat(x)
      lev <- unique(c(
        sort(unique(z)),
        "__UNKNOWN__",
        "__OTHER__"
      ))

      enc[[v]] <- list(
        type = "categorical",
        levels = lev
      )
    }
  }

  enc
}

encode_dense <- function(d, predictors, enc) {
  n <- nrow(d)
  blocks <- vector("list", length(predictors))

  for (i in seq_along(predictors)) {
    v <- predictors[i]
    info <- enc[[v]]
    x <- d[[v]]

    if (info$type == "numeric") {
      m <- matrix(
        as.numeric(x),
        nrow = n,
        ncol = 1L
      )
      colnames(m) <- v
      blocks[[i]] <- m
    } else {
      z <- clean_cat(x)
      z[!z %in% info$levels] <- "__OTHER__"

      f <- factor(
        z,
        levels = info$levels
      )

      m <- matrix(
        0,
        nrow = n,
        ncol = length(info$levels)
      )

      m[
        cbind(
          seq_len(n),
          as.integer(f)
        )
      ] <- 1

      colnames(m) <- paste0(
        v,
        "__",
        make.names(info$levels, unique = TRUE)
      )

      blocks[[i]] <- m
    }
  }

  X <- do.call(cbind, blocks)
  storage.mode(X) <- "double"
  X
}

# ------------------------------------------------------------------------------
# 4. Quantile helpers / metrics
# ------------------------------------------------------------------------------

extract_best_iteration <- function(model, fallback) {
  out <- tryCatch(
    as.numeric(
      reticulate::py_to_r(
        model$best_iteration
      )
    ),
    error = function(e) NA_real_
  )

  if (
    length(out) >= 1L &&
      is.finite(out[1L])
  ) {
    return(as.integer(out[1L] + 1L))
  }

  as.integer(fallback)
}

extract_best_score <- function(model) {
  out <- tryCatch(
    as.numeric(
      reticulate::py_to_r(
        model$best_score
      )
    ),
    error = function(e) NA_real_
  )

  if (
    length(out) >= 1L &&
      is.finite(out[1L])
  ) {
    return(out[1L])
  }

  NA_real_
}

as_quantile_matrix <- function(pred, n, k = 3L) {
  p <- tryCatch(
    reticulate::py_to_r(pred),
    error = function(e) pred
  )

  if (
    is.matrix(p) ||
      is.data.frame(p)
  ) {
    p <- as.matrix(p)
  } else {
    p <- matrix(
      as.numeric(p),
      nrow = n,
      ncol = k,
      byrow = TRUE
    )
  }

  if (
    nrow(p) != n ||
      ncol(p) != k
  ) {
    stop(
      "Unexpected quantile prediction shape: ",
      paste(dim(p), collapse = " x "),
      call. = FALSE
    )
  }

  storage.mode(p) <- "double"
  p
}

rearrange_three <- function(p) {
  lo <- pmin(p[, 1L], p[, 2L], p[, 3L])
  hi <- pmax(p[, 1L], p[, 2L], p[, 3L])
  mid <- rowSums(p) - lo - hi

  cbind(
    q10 = lo,
    q50 = mid,
    q90 = hi
  )
}

pinball_loss <- function(y, q, alpha) {
  e <- y - q

  mean(
    ifelse(
      e >= 0,
      alpha * e,
      (alpha - 1) * e
    ),
    na.rm = TRUE
  )
}

duration_metrics <- function(
    y,
    q10,
    q50,
    q90,
    benchmark_median
) {
  mae <- mean(
    abs(y - q50),
    na.rm = TRUE
  )

  benchmark_mae <- mean(
    abs(y - benchmark_median),
    na.rm = TRUE
  )

  data.table(
    N = length(y),
    MAE_days = mae,
    constant_median_MAE_days = benchmark_mae,
    improvement_vs_constant_median_pct =
      if (
        is.finite(benchmark_mae) &&
          benchmark_mae > 0
      ) {
        100 * (
          benchmark_mae -
            mae
        ) /
          benchmark_mae
      } else {
        NA_real_
      },
    median_bias_days =
      median(
        q50 - y,
        na.rm = TRUE
      ),
    pinball_q10 =
      pinball_loss(
        y,
        q10,
        0.10
      ),
    pinball_q50 =
      pinball_loss(
        y,
        q50,
        0.50
      ),
    pinball_q90 =
      pinball_loss(
        y,
        q90,
        0.90
      ),
    mean_pinball =
      mean(c(
        pinball_loss(
          y,
          q10,
          0.10
        ),
        pinball_loss(
          y,
          q50,
          0.50
        ),
        pinball_loss(
          y,
          q90,
          0.90
        )
      )),
    empirical_q10 =
      mean(
        y <= q10,
        na.rm = TRUE
      ),
    empirical_q50 =
      mean(
        y <= q50,
        na.rm = TRUE
      ),
    empirical_q90 =
      mean(
        y <= q90,
        na.rm = TRUE
      ),
    q10_q90_coverage =
      mean(
        y >= q10 &
          y <= q90,
        na.rm = TRUE
      ),
    median_q10_q90_width =
      median(
        q90 - q10,
        na.rm = TRUE
      )
  )
}

# ------------------------------------------------------------------------------
# 5. Duration specs and rolling folds
# ------------------------------------------------------------------------------

duration_specs <- list(
  hospital_los = list(
    y_var = "hospital_los_days",
    lower_bound = 0,
    eligible = function(y) {
      is.finite(y) &
        y >= 0
    }
  ),
  icu_los_conditional = list(
    y_var = "icu_los_days",
    lower_bound = 1,
    eligible = function(y) {
      is.finite(y) &
        y > 0
    }
  ),
  ventilator_days_conditional = list(
    y_var = "ventilator_days_clean",
    lower_bound = 1,
    eligible = function(y) {
      is.finite(y) &
        y > 0
    }
  )
)

for (nm in names(duration_specs)) {
  if (!duration_specs[[nm]]$y_var %in% names(dt)) {
    stop(
      "Frozen dataset lacks duration outcome: ",
      duration_specs[[nm]]$y_var,
      call. = FALSE
    )
  }
}

folds <- list(
  `2022` = list(
    train = 2020L,
    tune = 2021L,
    dev = 2020:2021,
    eval = 2022L
  ),
  `2023` = list(
    train = 2020:2021,
    tune = 2022L,
    dev = 2020:2022,
    eval = 2023L
  ),
  `2024` = list(
    train = 2020:2022,
    tune = 2023L,
    dev = 2020:2023,
    eval = 2024L
  )
)

# ------------------------------------------------------------------------------
# 6. One duration endpoint x temporal fold
# ------------------------------------------------------------------------------

fit_one_duration <- function(id, fold) {
  sp <- duration_specs[[id]]

  y_raw <- suppressWarnings(
    as.numeric(
      as.character(
        dt[[sp$y_var]]
      )
    )
  )

  eligible <- sp$eligible(y_raw)

  train_idx <- which(
    dt$admission_year %in% fold$train &
      eligible
  )
  tune_idx <- which(
    dt$admission_year %in% fold$tune &
      eligible
  )
  dev_idx <- which(
    dt$admission_year %in% fold$dev &
      eligible
  )
  eval_idx <- which(
    dt$admission_year %in% fold$eval &
      eligible
  )

  if (
    min(
      length(train_idx),
      length(tune_idx),
      length(dev_idx),
      length(eval_idx)
    ) <
      100L
  ) {
    stop(
      "Unexpectedly small temporal duration partition for ",
      id,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  y_train <- log1p(y_raw[train_idx])
  y_tune <- log1p(y_raw[tune_idx])

  enc_train <- fit_encoder(
    dt[train_idx],
    final_core
  )

  X_train <- encode_dense(
    dt[train_idx],
    final_core,
    enc_train
  )
  X_tune <- encode_dense(
    dt[tune_idx],
    final_core,
    enc_train
  )

  if (!identical(colnames(X_train), colnames(X_tune))) {
    stop(
      "Duration tuning encoder mismatch for ",
      id,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  dtrain <- xgb$DMatrix(
    X_train,
    label = y_train,
    missing = NA_real_
  )
  dtune <- xgb$DMatrix(
    X_tune,
    label = y_tune,
    missing = NA_real_
  )

  grid_rows <- vector(
    "list",
    nrow(search_grid)
  )

  for (i in seq_len(nrow(search_grid))) {
    g <- search_grid[i]

    params <- reticulate::dict(
      objective = "reg:quantileerror",
      quantile_alpha = c(
        0.10,
        0.50,
        0.90
      ),
      eta = ETA,
      max_depth =
        as.integer(
          g$max_depth
        ),
      min_child_weight =
        as.numeric(
          g$min_child_weight
        ),
      subsample =
        as.numeric(
          g$subsample
        ),
      colsample_bytree =
        as.numeric(
          g$colsample_bytree
        ),
      lambda =
        as.numeric(
          g$lambda
        ),
      tree_method = "hist",
      device = "cuda",
      nthread =
        as.integer(
          CPU_THREADS
        ),
      seed =
        as.integer(
          MODEL_SEED +
            fold$eval +
            i
        )
    )

    cat(
      "\n",
      id,
      " | eval ",
      fold$eval,
      " | ",
      g$config_id,
      " (",
      i,
      "/",
      nrow(search_grid),
      ")\n",
      sep = ""
    )

    fit <- xgb$train(
      params = params,
      dtrain = dtrain,
      num_boost_round =
        as.integer(
          MAX_ROUNDS
        ),
      evals = list(
        reticulate::tuple(
          dtrain,
          "train"
        ),
        reticulate::tuple(
          dtune,
          "tune"
        )
      ),
      early_stopping_rounds =
        as.integer(
          EARLY_STOP
        ),
      maximize = FALSE,
      verbose_eval = FALSE
    )

    grid_rows[[i]] <- cbind(
      data.table(
        outcome = id,
        evaluation_year = fold$eval,
        best_iteration =
          extract_best_iteration(
            fit,
            MAX_ROUNDS
          ),
        tune_quantile_score =
          extract_best_score(
            fit
          )
      ),
      g
    )

    rm(fit)
    gc()
  }

  grid_table <- rbindlist(
    grid_rows,
    fill = TRUE
  )

  if (
    all(
      !is.finite(
        grid_table$tune_quantile_score
      )
    )
  ) {
    stop(
      "All duration tuning scores are non-finite for ",
      id,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  setorder(
    grid_table,
    tune_quantile_score,
    best_iteration,
    max_depth,
    -lambda
  )

  selected <- copy(
    grid_table[1L]
  )

  fwrite(
    grid_table,
    file.path(
      GRID_DIR,
      paste0(
        id,
        "__eval",
        fold$eval,
        "__grid.csv"
      )
    )
  )

  rm(
    X_train,
    X_tune,
    dtrain,
    dtune,
    enc_train,
    grid_rows
  )
  gc()

  # Refit on all pre-evaluation development years.
  enc_dev <- fit_encoder(
    dt[dev_idx],
    final_core
  )

  X_dev <- encode_dense(
    dt[dev_idx],
    final_core,
    enc_dev
  )
  X_eval <- encode_dense(
    dt[eval_idx],
    final_core,
    enc_dev
  )

  if (!identical(colnames(X_dev), colnames(X_eval))) {
    stop(
      "Duration refit encoder mismatch for ",
      id,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  y_dev_log <- log1p(
    y_raw[dev_idx]
  )

  ddev <- xgb$DMatrix(
    X_dev,
    label = y_dev_log,
    missing = NA_real_
  )

  deval <- xgb$DMatrix(
    X_eval,
    missing = NA_real_
  )

  final_params <- reticulate::dict(
    objective = "reg:quantileerror",
    quantile_alpha = c(
      0.10,
      0.50,
      0.90
    ),
    eta = ETA,
    max_depth =
      as.integer(
        selected$max_depth
      ),
    min_child_weight =
      as.numeric(
        selected$min_child_weight
      ),
    subsample =
      as.numeric(
        selected$subsample
      ),
    colsample_bytree =
      as.numeric(
        selected$colsample_bytree
      ),
    lambda =
      as.numeric(
        selected$lambda
      ),
    tree_method = "hist",
    device = "cuda",
    nthread =
      as.integer(
        CPU_THREADS
      ),
    seed =
      as.integer(
        MODEL_SEED +
          10000L +
          fold$eval
      )
  )

  final_fit <- xgb$train(
    params = final_params,
    dtrain = ddev,
    num_boost_round =
      as.integer(
        selected$best_iteration
      ),
    verbose_eval = FALSE
  )

  raw_pred <- as_quantile_matrix(
    final_fit$predict(deval),
    n = length(eval_idx),
    k = 3L
  )

  raw_crossing_rate <- mean(
    raw_pred[, 1L] >
      raw_pred[, 2L] |
      raw_pred[, 2L] >
        raw_pred[, 3L]
  )

  pred_log <- rearrange_three(
    raw_pred
  )

  pred_days <- expm1(
    pred_log
  )

  pred_days <- pmax(
    pred_days,
    sp$lower_bound
  )

  # Back-transform can introduce tiny numerical disorder after truncation.
  pred_days <- rearrange_three(
    pred_days
  )

  y_eval <- y_raw[eval_idx]
  benchmark_median <- median(
    y_raw[dev_idx],
    na.rm = TRUE
  )

  met <- duration_metrics(
    y = y_eval,
    q10 = pred_days[, "q10"],
    q50 = pred_days[, "q50"],
    q90 = pred_days[, "q90"],
    benchmark_median =
      benchmark_median
  )

  met[, `:=`(
    outcome = id,
    evaluation_year = fold$eval,
    predictor_policy = "final_core",
    n_predictors = length(final_core),
    selected_config =
      selected$config_id,
    selected_rounds =
      as.integer(
        selected$best_iteration
      ),
    tune_quantile_score =
      as.numeric(
        selected$tune_quantile_score
      ),
    raw_quantile_crossing_rate =
      raw_crossing_rate,
    benchmark_development_median_days =
      benchmark_median
  )]

  pred_dt <- data.table(
    admission_year =
      dt$admission_year[eval_idx],
    inc_key =
      dt$inc_key[eval_idx],
    observed_days =
      y_eval,
    q10_days =
      pred_days[, "q10"],
    q50_median_days =
      pred_days[, "q50"],
    q90_days =
      pred_days[, "q90"]
  )

  fwrite(
    pred_dt,
    file.path(
      PRED_DIR,
      paste0(
        id,
        "__eval",
        fold$eval,
        "__Q10_Q50_Q90_predictions.csv"
      )
    )
  )

  saveRDS(
    list(
      model = final_fit,
      encoder = enc_dev,
      outcome = id,
      evaluation_year = fold$eval,
      predictors = final_core,
      selected_config = selected,
      lower_bound = sp$lower_bound,
      transformation = "log1p"
    ),
    file.path(
      MODEL_DIR,
      paste0(
        id,
        "__eval",
        fold$eval,
        "__model.rds"
      )
    ),
    compress = FALSE
  )

  qc <- data.table(
    outcome = id,
    evaluation_year = fold$eval,
    n_train = length(train_idx),
    n_tune = length(tune_idx),
    n_dev_refit = length(dev_idx),
    n_eval = length(eval_idx),
    selected_config =
      selected$config_id,
    selected_rounds =
      as.integer(
        selected$best_iteration
      ),
    hit_round_ceiling =
      as.integer(
        selected$best_iteration
      ) >=
        (MAX_ROUNDS - 5L),
    raw_quantile_crossing_rate =
      raw_crossing_rate,
    n_predictors =
      length(final_core),
    n_encoded_features =
      ncol(X_dev)
  )

  rm(
    X_dev,
    X_eval,
    ddev,
    deval,
    final_fit,
    raw_pred,
    pred_log,
    pred_days,
    pred_dt,
    enc_dev
  )
  gc()

  list(
    metrics = met,
    selected = selected,
    qc = qc
  )
}

# ------------------------------------------------------------------------------
# 7. Run 3 outcomes x 3 folds = 9 final temporal evaluations
# ------------------------------------------------------------------------------

run_plan <- CJ(
  outcome = names(duration_specs),
  evaluation_year =
    as.integer(
      names(folds)
    ),
  unique = TRUE
)

run_plan[
  ,
  outcome_order :=
    match(
      outcome,
      names(duration_specs)
    )
]
setorder(
  run_plan,
  outcome_order,
  evaluation_year
)
run_plan[, outcome_order := NULL]

fwrite(
  run_plan,
  file.path(OUT_DIR, "02_RUN_PLAN.csv")
)

metric_rows <- list()
selected_rows <- list()
qc_rows <- list()

for (i in seq_len(nrow(run_plan))) {
  id <- run_plan$outcome[i]
  year <- as.character(
    run_plan$evaluation_year[i]
  )

  cat(
    "\n============================================================\n",
    "[",
    i,
    "/",
    nrow(run_plan),
    "] DURATION: ",
    id,
    " | eval ",
    year,
    "\n============================================================\n",
    sep = ""
  )

  res <- fit_one_duration(
    id,
    folds[[year]]
  )

  metric_rows[[length(metric_rows) + 1L]] <-
    res$metrics

  s <- copy(
    res$selected
  )
  s[, `:=`(
    outcome = id,
    evaluation_year =
      as.integer(year)
  )]
  selected_rows[[length(selected_rows) + 1L]] <- s

  qc_rows[[length(qc_rows) + 1L]] <-
    res$qc

  fwrite(
    rbindlist(metric_rows, fill = TRUE),
    file.path(
      OUT_DIR,
      "05_FINAL_DURATION_TEMPORAL_METRICS_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(selected_rows, fill = TRUE),
    file.path(
      OUT_DIR,
      "04_SELECTED_DURATION_CONFIGS_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(qc_rows, fill = TRUE),
    file.path(
      OUT_DIR,
      "06_DURATION_TRAINING_QC_CHECKPOINT.csv"
    )
  )
}

metrics <- rbindlist(
  metric_rows,
  fill = TRUE
)

selected <- rbindlist(
  selected_rows,
  fill = TRUE
)

qc <- rbindlist(
  qc_rows,
  fill = TRUE
)

fwrite(
  selected,
  file.path(
    OUT_DIR,
    "04_SELECTED_DURATION_CONFIGS.csv"
  )
)

fwrite(
  metrics,
  file.path(
    OUT_DIR,
    "05_FINAL_DURATION_TEMPORAL_METRICS.csv"
  )
)

fwrite(
  qc,
  file.path(
    OUT_DIR,
    "06_DURATION_TRAINING_QC.csv"
  )
)

temporal_summary <- metrics[
  ,
  .(
    min_MAE_days =
      min(
        MAE_days,
        na.rm = TRUE
      ),
    median_MAE_days =
      median(
        MAE_days,
        na.rm = TRUE
      ),
    max_MAE_days =
      max(
        MAE_days,
        na.rm = TRUE
      ),
    median_improvement_vs_constant_pct =
      median(
        improvement_vs_constant_median_pct,
        na.rm = TRUE
      ),
    median_bias_days =
      median(
        median_bias_days,
        na.rm = TRUE
      ),
    median_q10_empirical =
      median(
        empirical_q10,
        na.rm = TRUE
      ),
    median_q50_empirical =
      median(
        empirical_q50,
        na.rm = TRUE
      ),
    median_q90_empirical =
      median(
        empirical_q90,
        na.rm = TRUE
      ),
    median_q10_q90_coverage =
      median(
        q10_q90_coverage,
        na.rm = TRUE
      ),
    median_q10_q90_width =
      median(
        median_q10_q90_width,
        na.rm = TRUE
      )
  ),
  by = outcome
]

fwrite(
  temporal_summary,
  file.path(
    OUT_DIR,
    "07_DURATION_TEMPORAL_SUMMARY.csv"
  )
)

if (
  any(
    qc$hit_round_ceiling,
    na.rm = TRUE
  )
) {
  warning(
    "At least one duration model hit the 3000-round ceiling. Review 06_DURATION_TRAINING_QC.csv before final lock."
  )
}

writeLines(
  c(
    "SPINE-TRACT REPAIRED-DX DURATION QUANTILES COMPLETE",
    "",
    "Predictors: corrected 54-variable final_core.",
    "Outcomes: HLOS unconditional; ICU LOS conditional on ICU use; ventilator days conditional on ventilation.",
    "Quantiles: Q10/Q50/Q90 on log1p(days), back-transformed and row-wise rearranged.",
    "Temporal folds: evaluation 2022, 2023, 2024.",
    "Each endpoint/fold was fully retuned over the same 24 structural configurations used by repaired Stage 3A.",
    "",
    "Review first:",
    "  05_FINAL_DURATION_TEMPORAL_METRICS.csv",
    "  07_DURATION_TEMPORAL_SUMMARY.csv",
    "  06_DURATION_TRAINING_QC.csv",
    "  04_SELECTED_DURATION_CONFIGS.csv"
  ),
  file.path(
    OUT_DIR,
    "README_RESULTS.txt"
  )
)

cat(
  "\n============================================================\n"
)
cat(
  "REPAIRED-DX DURATION QUANTILES COMPLETE\n"
)
cat(
  "============================================================\n"
)
print(
  metrics[
    order(
      outcome,
      evaluation_year
    )
  ]
)
cat(
  "\nOutput:\n  ",
  OUT_DIR,
  "\n",
  sep = ""
)
