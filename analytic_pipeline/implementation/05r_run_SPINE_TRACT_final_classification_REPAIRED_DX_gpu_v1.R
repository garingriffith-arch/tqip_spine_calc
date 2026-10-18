# ==============================================================================
# 05r_run_SPINE_TRACT_final_classification_REPAIRED_DX_gpu_v1.R
# SPINE-TRACT Stage 3A-R — final repaired-DX rolling-origin classification
#
# PURPOSE
#   Refit the final manuscript-facing SPINE-TRACT classification models after
#   repair of the diagnosis-derived phenotypes and Stage 2B/2C re-adjudication.
#
# LOCKED PREDICTOR POLICY
#   final_core (54 predictors):
#     repaired clinical backbone
#     + coarse spine representation INCLUDING dx_tl_soft_tissue_any
#     + detailed chest injury
#     + pupillary response
#
#   final_context (57 predictors):
#     final_core + race + ethnicity + payer
#
#   Endpoints:
#     discharge3_final               -> final_context
#     hlos3_final                    -> final_context
#     icu3_final                     -> final_core
#     vent3_final                    -> final_core
#     operative_spine_intervention   -> final_core
#
# TEMPORAL DESIGN
#   Fold 1: train 2020      | tune 2021 | refit 2020-21 | eval 2022
#   Fold 2: train 2020-21   | tune 2022 | refit 2020-22 | eval 2023
#   Fold 3: train 2020-22   | tune 2023 | refit 2020-23 | eval 2024
#
# HYPERPARAMETER SEARCH
#   Exact 24-configuration prespecified structural strategy used in TBI-TRACT:
#     max_depth            2,3,4,5
#     min_child_weight     5,10,20
#     subsample            0.75,0.90,1.00
#     colsample_bytree     0.75,0.90,1.00
#     lambda               1,5
#     + historical heuristic row depth=3, child=10, subsample=.85,
#       colsample=.85, lambda=1
#   23 full-grid configurations are sampled with fixed seed before fitting.
#   eta fixed at 0.05. Boosting rounds selected by early stopping in the
#   fold-specific tuning year.
#
# IMPORTANT
#   - 2024 is development-informed forward temporal evaluation, not external
#     validation.
#   - No evaluation-year recalibration.
#   - Numeric predictor missingness remains native to XGBoost.
#   - Categorical missingness uses explicit __UNKNOWN__; unseen levels use
#     __OTHER__.
#   - This writes to a NEW repaired-DX output directory. Old Stage 3 outputs are
#     preserved for audit only.
# ==============================================================================

rm(list = ls())
gc()

required <- c("data.table", "Matrix", "xgboost")
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
  library(Matrix)
  library(xgboost)
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

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification_REPAIRED_DX",
  "rolling_origin_gpu"
)
GRID_DIR <- file.path(OUT_DIR, "grid_search")
PRED_DIR <- file.path(OUT_DIR, "predictions")
MODEL_DIR <- file.path(OUT_DIR, "models")

for (d in c(OUT_DIR, GRID_DIR, PRED_DIR, MODEL_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

if (!file.exists(DATA_FILE) || !file.exists(SPEC_FILE)) {
  stop(
    "Missing repaired frozen dataset or predictor specification.",
    call. = FALSE
  )
}

MODEL_SEED <- 20261002L
GRID_SEED <- 20260907L
ETA <- 0.05
MAX_ROUNDS <- 3000L
EARLY_STOP <- 75L
PRINT_EVERY <- 100L

logical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(logical_cores) || logical_cores < 1L) logical_cores <- 1L
CPU_THREADS <- max(
  1L,
  min(
    logical_cores,
    as.integer(round(logical_cores * 0.80))
  )
)
data.table::setDTthreads(CPU_THREADS)

# ------------------------------------------------------------------------------
# 1. Load repaired data and reconstruct locked 54/57 predictor sets
# ------------------------------------------------------------------------------

dt <- as.data.table(readRDS(DATA_FILE))
spec <- readRDS(SPEC_FILE)

dt[, inc_key := as.character(inc_key)]
dt[, admission_year := as.integer(admission_year)]

if (nrow(dt) != 424017L) {
  stop(
    "Expected repaired frozen primary N=424017; observed ",
    nrow(dt),
    ".",
    call. = FALSE
  )
}

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

final_context <- unique(c(
  final_core,
  pg$social_context
))

if (length(final_core) != 54L) {
  stop(
    "Expected corrected final_core length 54; observed ",
    length(final_core),
    ".",
    call. = FALSE
  )
}

if (length(final_context) != 57L) {
  stop(
    "Expected corrected final_context length 57; observed ",
    length(final_context),
    ".",
    call. = FALSE
  )
}

if (!"dx_tl_soft_tissue_any" %in% final_core) {
  stop(
    "HARD FAIL: corrected final_core does not include dx_tl_soft_tissue_any.",
    call. = FALSE
  )
}

missing_predictors <- setdiff(
  unique(c(final_core, final_context)),
  names(dt)
)
if (length(missing_predictors) > 0L) {
  stop(
    "Frozen dataset lacks locked predictor(s): ",
    paste(missing_predictors, collapse = ", "),
    call. = FALSE
  )
}

# Diagnosis repair hard preflight.
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
      function(v) {
        sum(dt[[v]] == 1L, na.rm = TRUE)
      },
      integer(1)
    ),
    missing_n = vapply(
      predictor,
      function(v) {
        sum(is.na(dt[[v]]))
      },
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
    "HARD FAIL: repaired diagnosis preflight failed.",
    call. = FALSE
  )
}

predictor_manifest <- rbindlist(
  list(
    data.table(
      predictor_policy = "final_core",
      n_predictors = length(final_core),
      predictor = final_core
    ),
    data.table(
      predictor_policy = "final_context",
      n_predictors = length(final_context),
      predictor = final_context
    )
  )
)

fwrite(
  predictor_manifest,
  file.path(OUT_DIR, "01_FINAL_PREDICTOR_MANIFEST.csv")
)

# ------------------------------------------------------------------------------
# 2. Construct final endpoint variables
# ------------------------------------------------------------------------------

if (!"discharge4" %in% names(dt)) {
  stop("Frozen dataset lacks discharge4.", call. = FALSE)
}

dt[, discharge3_final := factor(
  fcase(
    as.character(discharge4) == "Home/home health",
      "Home/home health",
    as.character(discharge4) %in% c(
      "Inpatient rehabilitation",
      "Other post-acute facility"
    ),
      "Post-acute facility",
    as.character(discharge4) == "Death/hospice",
      "Death/hospice",
    default = NA_character_
  ),
  levels = c(
    "Home/home health",
    "Post-acute facility",
    "Death/hospice"
  )
)]

dt[, hlos3_final := factor(
  as.character(hlos3),
  levels = c(
    "<=7 days",
    "8-19 days",
    ">=20 days"
  )
)]

dt[, icu3_final := factor(
  as.character(icu3),
  levels = c(
    "No ICU",
    "ICU 1-7 days",
    "ICU >=8 days"
  )
)]

dt[, vent3_final := factor(
  as.character(vent3),
  levels = c(
    "No ventilation",
    "Ventilation 1-7 days",
    "Ventilation >=8 days"
  )
)]

endpoint_specs <- list(
  discharge3_final = list(
    type = "multiclass",
    levels = c(
      "Home/home health",
      "Post-acute facility",
      "Death/hospice"
    ),
    predictors = final_context,
    policy = "final_context"
  ),
  hlos3_final = list(
    type = "multiclass",
    levels = c(
      "<=7 days",
      "8-19 days",
      ">=20 days"
    ),
    predictors = final_context,
    policy = "final_context"
  ),
  icu3_final = list(
    type = "multiclass",
    levels = c(
      "No ICU",
      "ICU 1-7 days",
      "ICU >=8 days"
    ),
    predictors = final_core,
    policy = "final_core"
  ),
  vent3_final = list(
    type = "multiclass",
    levels = c(
      "No ventilation",
      "Ventilation 1-7 days",
      "Ventilation >=8 days"
    ),
    predictors = final_core,
    policy = "final_core"
  ),
  operative_spine_intervention = list(
    type = "binary",
    predictors = final_core,
    policy = "final_core"
  )
)

# ------------------------------------------------------------------------------
# 3. Prespecified 24-configuration structural grid
# ------------------------------------------------------------------------------

full_grid <- CJ(
  max_depth = c(2L, 3L, 4L, 5L),
  min_child_weight = c(5, 10, 20),
  subsample = c(0.75, 0.90, 1.00),
  colsample_bytree = c(0.75, 0.90, 1.00),
  lambda = c(1, 5),
  unique = TRUE
)

default_row <- data.table(
  max_depth = 3L,
  min_child_weight = 10,
  subsample = 0.85,
  colsample_bytree = 0.85,
  lambda = 1
)

set.seed(GRID_SEED)
sampled_grid <- full_grid[
  sample(.N, size = 23L, replace = FALSE)
]

search_grid <- unique(
  rbind(
    default_row,
    sampled_grid,
    fill = TRUE
  )
)

# Defensive completion if the default happened to be sampled.
if (nrow(search_grid) < 24L) {
  remaining <- fsetdiff(
    full_grid,
    search_grid[
      subsample != 0.85 |
        colsample_bytree != 0.85
    ]
  )
  need <- 24L - nrow(search_grid)
  search_grid <- rbind(
    search_grid,
    remaining[seq_len(need)],
    fill = TRUE
  )
}

search_grid[, config_id := sprintf("cfg_%02d", seq_len(.N))]

setcolorder(
  search_grid,
  c(
    "config_id",
    "max_depth",
    "min_child_weight",
    "subsample",
    "colsample_bytree",
    "lambda"
  )
)

if (nrow(search_grid) != 24L) {
  stop(
    "Expected exactly 24 prespecified grid configurations; observed ",
    nrow(search_grid),
    ".",
    call. = FALSE
  )
}

fwrite(
  search_grid,
  file.path(OUT_DIR, "02_PRESPECIFIED_SEARCH_GRID.csv")
)

# ------------------------------------------------------------------------------
# 4. GPU detection / compatibility
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

gpu_params_for_version <- function() {
  ver <- packageVersion("xgboost")

  if (ver < package_version("2.0.0")) {
    list(
      mode = paste0(
        "R CUDA legacy xgboost ",
        as.character(ver),
        " | gpu_hist | gpu_id=0"
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
        "R CUDA modern xgboost ",
        as.character(ver),
        " | hist | cuda:0"
      ),
      extra = list(
        tree_method = "hist",
        device = "cuda:0"
      )
    )
  }
}

assert_cuda_booster <- function(model, context) {
  txt <- tolower(config_to_text(model))
  ver <- packageVersion("xgboost")

  ok <- if (ver < package_version("2.0.0")) {
    grepl("grow_gpu_hist|gpu_hist", txt, perl = TRUE)
  } else {
    grepl(
      '"device"\\s*:\\s*"cuda(:0)?"',
      txt,
      perl = TRUE
    ) ||
      grepl("grow_gpu_hist", txt, perl = TRUE)
  }

  if (!ok) {
    stop(
      "HARD FAIL: GPU not verified for ",
      context,
      ".",
      call. = FALSE
    )
  }

  invisible(TRUE)
}

compute <- gpu_params_for_version()

# Small hard-fail GPU probe.
set.seed(1L)
probe_x <- matrix(
  rnorm(20000),
  nrow = 1000L,
  ncol = 20L
)
probe_y <- rbinom(1000L, 1L, 0.3)
probe_d <- xgb.DMatrix(probe_x, label = probe_y)

probe_params <- c(
  list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth = 2L,
    eta = 0.2,
    nthread = CPU_THREADS,
    verbosity = 0
  ),
  compute$extra
)

probe_fit <- xgb.train(
  params = probe_params,
  data = probe_d,
  nrounds = 3L,
  verbose = 0
)

assert_cuda_booster(probe_fit, "startup GPU probe")

fwrite(
  data.table(
    xgboost_version = as.character(packageVersion("xgboost")),
    compute_mode = compute$mode,
    logical_cores = logical_cores,
    allocated_cpu_threads = CPU_THREADS
  ),
  file.path(OUT_DIR, "03_COMPUTE_ENVIRONMENT.csv")
)

rm(probe_x, probe_y, probe_d, probe_fit)
gc()

# ------------------------------------------------------------------------------
# 5. Temporal encoder
# ------------------------------------------------------------------------------

is_numeric_predictor <- function(x) {
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

    if (is_numeric_predictor(x)) {
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

encode_with <- function(d, predictors, enc) {
  blocks <- vector("list", length(predictors))

  for (i in seq_along(predictors)) {
    v <- predictors[i]
    info <- enc[[v]]
    x <- d[[v]]

    if (info$type == "numeric") {
      m <- Matrix(
        as.numeric(x),
        ncol = 1L,
        sparse = TRUE
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

      m <- sparseMatrix(
        i = seq_along(f),
        j = as.integer(f),
        x = 1,
        dims = c(
          length(f),
          length(info$levels)
        )
      )

      colnames(m) <- paste0(
        v,
        "__",
        make.names(info$levels, unique = TRUE)
      )

      blocks[[i]] <- m
    }
  }

  do.call(cbind, blocks)
}

# ------------------------------------------------------------------------------
# 6. Metrics
# ------------------------------------------------------------------------------

clamp_prob <- function(p, eps = 1e-7) {
  pmin(
    pmax(as.numeric(p), eps),
    1 - eps
  )
}

fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p) & is.finite(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  n1 <- as.double(sum(y == 1L))
  n0 <- as.double(sum(y == 0L))

  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  r <- rank(p, ties.method = "average")

  (
    sum(r[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}

fast_auprc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p) & is.finite(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  npos <- sum(y == 1L)
  if (npos == 0L) {
    return(NA_real_)
  }

  ord <- order(p, decreasing = TRUE)
  ys <- y[ord]

  precision <- cumsum(ys == 1L) / seq_along(ys)

  mean(
    precision[ys == 1L]
  )
}

binary_calibration <- function(y, p) {
  pp <- clamp_prob(p)
  lp <- qlogis(pp)

  if (length(unique(y)) < 2L) {
    return(
      c(
        intercept = NA_real_,
        slope = NA_real_
      )
    )
  }

  fit_i <- tryCatch(
    suppressWarnings(
      glm(
        y ~ offset(lp),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  fit_s <- tryCatch(
    suppressWarnings(
      glm(
        y ~ lp,
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  c(
    intercept =
      if (is.null(fit_i)) NA_real_ else unname(coef(fit_i)[1L]),
    slope =
      if (is.null(fit_s)) NA_real_ else unname(coef(fit_s)[2L])
  )
}

binary_metrics <- function(
    y,
    p,
    endpoint,
    eval_year,
    best_iter
) {
  pp <- clamp_prob(p)
  prev <- mean(y == 1L)
  cal <- binary_calibration(y, p)
  pr <- fast_auprc(y, p)

  data.table(
    endpoint = endpoint,
    evaluation_year = eval_year,
    metric_level = "binary",
    class = "Event",
    n = length(y),
    prevalence = prev,
    auroc = fast_auc(y, p),
    auprc = pr,
    pr_lift =
      if (is.finite(prev) && prev > 0) pr / prev else NA_real_,
    brier = mean((y - p)^2),
    logloss = -mean(
      y * log(pp) +
        (1 - y) * log(1 - pp)
    ),
    calibration_intercept = unname(cal["intercept"]),
    calibration_slope = unname(cal["slope"]),
    multiclass_logloss = NA_real_,
    multiclass_brier = NA_real_,
    best_iteration = best_iter
  )
}

multiclass_metrics <- function(
    y0,
    prob,
    levels,
    endpoint,
    eval_year,
    best_iter
) {
  k <- length(levels)
  y1 <- y0 + 1L

  class_rows <- rbindlist(
    lapply(
      seq_len(k),
      function(j) {
        yj <- as.integer(y1 == j)
        pj <- prob[, j]
        pp <- clamp_prob(pj)
        prev <- mean(yj == 1L)
        cal <- binary_calibration(yj, pj)
        pr <- fast_auprc(yj, pj)

        data.table(
          endpoint = endpoint,
          evaluation_year = eval_year,
          metric_level = "class",
          class = levels[j],
          n = length(yj),
          prevalence = prev,
          auroc = fast_auc(yj, pj),
          auprc = pr,
          pr_lift =
            if (is.finite(prev) && prev > 0) pr / prev else NA_real_,
          brier = mean((yj - pj)^2),
          logloss = -mean(
            yj * log(pp) +
              (1 - yj) * log(1 - pp)
          ),
          calibration_intercept = unname(cal["intercept"]),
          calibration_slope = unname(cal["slope"]),
          multiclass_logloss = NA_real_,
          multiclass_brier = NA_real_,
          best_iteration = best_iter
        )
      }
    ),
    fill = TRUE
  )

  ptrue <- prob[
    cbind(
      seq_along(y1),
      y1
    )
  ]
  ptrue <- pmax(ptrue, 1e-15)

  y_onehot <- matrix(
    0,
    nrow = length(y1),
    ncol = k
  )
  y_onehot[
    cbind(
      seq_along(y1),
      y1
    )
  ] <- 1

  mc_logloss <- -mean(log(ptrue))
  mc_brier <- mean(
    rowSums(
      (prob - y_onehot)^2
    )
  )

  class_rows[
    ,
    `:=`(
      multiclass_logloss = mc_logloss,
      multiclass_brier = mc_brier
    )
  ]

  macro <- data.table(
    endpoint = endpoint,
    evaluation_year = eval_year,
    metric_level = "macro",
    class = "MACRO",
    n = length(y1),
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
    multiclass_logloss = mc_logloss,
    multiclass_brier = mc_brier,
    best_iteration = best_iter
  )

  rbindlist(
    list(
      class_rows,
      macro
    ),
    fill = TRUE
  )
}

# ------------------------------------------------------------------------------
# 7. XGBoost helpers
# ------------------------------------------------------------------------------

make_params <- function(spec, g) {
  base <- list(
    eval_metric =
      if (spec$type == "multiclass") "mlogloss" else "logloss",
    eta = ETA,
    max_depth = as.integer(g$max_depth),
    min_child_weight = as.numeric(g$min_child_weight),
    subsample = as.numeric(g$subsample),
    colsample_bytree = as.numeric(g$colsample_bytree),
    lambda = as.numeric(g$lambda),
    nthread = CPU_THREADS,
    verbosity = 0
  )

  if (spec$type == "multiclass") {
    base$objective <- "multi:softprob"
    base$num_class <- length(spec$levels)
  } else {
    base$objective <- "binary:logistic"
  }

  c(base, compute$extra)
}

extract_best_iteration <- function(model, fallback) {
  b <- tryCatch(model$best_iteration, error = function(e) NULL)

  if (
    !is.null(b) &&
      length(b) == 1L &&
      is.finite(b)
  ) {
    return(as.integer(b))
  }

  a <- tryCatch(
    xgb.attr(model, "best_iteration"),
    error = function(e) NULL
  )

  if (
    !is.null(a) &&
      length(a) == 1L &&
      !is.na(a)
  ) {
    return(as.integer(a) + 1L)
  }

  as.integer(fallback)
}

extract_best_score <- function(model) {
  # First try native attribute.
  a <- tryCatch(
    xgb.attr(model, "best_score"),
    error = function(e) NULL
  )

  if (
    !is.null(a) &&
      length(a) >= 1L &&
      !is.na(a[1L])
  ) {
    z <- suppressWarnings(as.numeric(a[1L]))
    if (is.finite(z)) return(z)
  }

  # Then evaluation log.
  log <- tryCatch(
    model$evaluation_log,
    error = function(e) NULL
  )

  if (
    !is.null(log) &&
      nrow(log) > 0L
  ) {
    tune_cols <- grep(
      "^tune_",
      names(log),
      value = TRUE
    )
    if (length(tune_cols) >= 1L) {
      return(
        min(
          as.numeric(log[[tune_cols[1L]]]),
          na.rm = TRUE
        )
      )
    }
  }

  NA_real_
}

# ------------------------------------------------------------------------------
# 8. Rolling-origin folds
# ------------------------------------------------------------------------------

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
# 9. Fit one endpoint x fold: full grid tune -> refit -> evaluation
# ------------------------------------------------------------------------------

fit_one <- function(endpoint, fold) {
  spec_ep <- endpoint_specs[[endpoint]]
  predictors <- spec_ep$predictors

  if (spec_ep$type == "multiclass") {
    yf <- factor(
      as.character(dt[[endpoint]]),
      levels = spec_ep$levels
    )
    y_all <- as.integer(yf) - 1L
  } else {
    y_all <- as.integer(dt[[endpoint]])
  }

  train_idx <- which(dt$admission_year %in% fold$train)
  tune_idx <- which(dt$admission_year %in% fold$tune)
  dev_idx <- which(dt$admission_year %in% fold$dev)
  eval_idx <- which(dt$admission_year %in% fold$eval)

  if (spec_ep$type == "multiclass") {
    train_idx <- train_idx[!is.na(y_all[train_idx])]
    tune_idx <- tune_idx[!is.na(y_all[tune_idx])]
    dev_idx <- dev_idx[!is.na(y_all[dev_idx])]
    eval_idx <- eval_idx[!is.na(y_all[eval_idx])]
  } else {
    train_idx <- train_idx[
      !is.na(y_all[train_idx]) &
        y_all[train_idx] %in% c(0L, 1L)
    ]
    tune_idx <- tune_idx[
      !is.na(y_all[tune_idx]) &
        y_all[tune_idx] %in% c(0L, 1L)
    ]
    dev_idx <- dev_idx[
      !is.na(y_all[dev_idx]) &
        y_all[dev_idx] %in% c(0L, 1L)
    ]
    eval_idx <- eval_idx[
      !is.na(y_all[eval_idx]) &
        y_all[eval_idx] %in% c(0L, 1L)
    ]
  }

  if (
    length(train_idx) == 0L ||
      length(tune_idx) == 0L ||
      length(dev_idx) == 0L ||
      length(eval_idx) == 0L
  ) {
    stop(
      "Empty temporal partition for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  # Training-period encoder only for tuning.
  enc_train <- fit_encoder(
    dt[train_idx],
    predictors
  )

  X_train <- encode_with(
    dt[train_idx],
    predictors,
    enc_train
  )
  X_tune <- encode_with(
    dt[tune_idx],
    predictors,
    enc_train
  )

  if (!identical(colnames(X_train), colnames(X_tune))) {
    stop(
      "Tuning encoder column mismatch for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  dtrain <- xgb.DMatrix(
    X_train,
    label = y_all[train_idx],
    missing = NA_real_
  )
  dtune <- xgb.DMatrix(
    X_tune,
    label = y_all[tune_idx],
    missing = NA_real_
  )

  grid_rows <- vector("list", nrow(search_grid))

  for (i in seq_len(nrow(search_grid))) {
    g <- search_grid[i]
    params <- make_params(spec_ep, g)

    cat(
      "\n",
      endpoint,
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

    set.seed(MODEL_SEED + fold$eval + i)

    if (packageVersion("xgboost") < package_version("2.0.0")) {
      fit_es <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = MAX_ROUNDS,
        watchlist = list(
          train = dtrain,
          tune = dtune
        ),
        early_stopping_rounds = EARLY_STOP,
        maximize = FALSE,
        verbose = 0
      )
    } else {
      fit_es <- xgb.train(
        params = params,
        data = dtrain,
        nrounds = MAX_ROUNDS,
        evals = list(
          train = dtrain,
          tune = dtune
        ),
        early_stopping_rounds = EARLY_STOP,
        maximize = FALSE,
        verbose = 0
      )
    }

    if (i == 1L) {
      assert_cuda_booster(
        fit_es,
        paste0(
          endpoint,
          " eval ",
          fold$eval,
          " tuning"
        )
      )
    }

    grid_rows[[i]] <- cbind(
      data.table(
        endpoint = endpoint,
        predictor_policy = spec_ep$policy,
        evaluation_year = fold$eval,
        best_iteration =
          extract_best_iteration(
            fit_es,
            MAX_ROUNDS
          ),
        tune_score =
          extract_best_score(
            fit_es
          )
      ),
      g
    )

    rm(fit_es)
    gc()
  }

  grid_table <- rbindlist(
    grid_rows,
    fill = TRUE
  )

  setorder(
    grid_table,
    tune_score,
    best_iteration,
    max_depth,
    -lambda
  )

  selected <- copy(grid_table[1L])

  if (!is.finite(selected$tune_score)) {
    stop(
      "No finite tuning score for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  fwrite(
    grid_table,
    file.path(
      GRID_DIR,
      paste0(
        endpoint,
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

  # Refit encoder on all pre-evaluation development years.
  enc_dev <- fit_encoder(
    dt[dev_idx],
    predictors
  )

  X_dev <- encode_with(
    dt[dev_idx],
    predictors,
    enc_dev
  )
  X_eval <- encode_with(
    dt[eval_idx],
    predictors,
    enc_dev
  )

  if (!identical(colnames(X_dev), colnames(X_eval))) {
    stop(
      "Refit encoder column mismatch for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  ddev <- xgb.DMatrix(
    X_dev,
    label = y_all[dev_idx],
    missing = NA_real_
  )
  deval <- xgb.DMatrix(
    X_eval,
    missing = NA_real_
  )

  final_params <- make_params(
    spec_ep,
    selected
  )

  set.seed(MODEL_SEED + 10000L + fold$eval)

  fit <- xgb.train(
    params = final_params,
    data = ddev,
    nrounds =
      as.integer(
        selected$best_iteration
      ),
    verbose = 0
  )

  assert_cuda_booster(
    fit,
    paste0(
      endpoint,
      " eval ",
      fold$eval,
      " final refit"
    )
  )

  raw <- predict(
    fit,
    deval
  )

  if (spec_ep$type == "multiclass") {
    k <- length(spec_ep$levels)

    prob <- if (
      is.matrix(raw) ||
        is.data.frame(raw)
    ) {
      as.matrix(raw)
    } else {
      matrix(
        as.numeric(raw),
        ncol = k,
        byrow = TRUE
      )
    }

    if (
      nrow(prob) != length(eval_idx) ||
        ncol(prob) != k
    ) {
      stop(
        "Unexpected multiclass prediction shape for ",
        endpoint,
        " eval ",
        fold$eval,
        ".",
        call. = FALSE
      )
    }

    met <- multiclass_metrics(
      y0 = y_all[eval_idx],
      prob = prob,
      levels = spec_ep$levels,
      endpoint = endpoint,
      eval_year = fold$eval,
      best_iter =
        as.integer(
          selected$best_iteration
        )
    )

    pred_dt <- data.table(
      admission_year =
        dt$admission_year[eval_idx],
      inc_key =
        dt$inc_key[eval_idx],
      observed_class =
        spec_ep$levels[
          y_all[eval_idx] + 1L
        ]
    )

    for (j in seq_along(spec_ep$levels)) {
      pred_dt[
        ,
        paste0(
          "p_",
          make.names(
            spec_ep$levels[j]
          )
        ) :=
          prob[, j]
      ]
    }
  } else {
    p <- as.numeric(raw)

    met <- binary_metrics(
      y = y_all[eval_idx],
      p = p,
      endpoint = endpoint,
      eval_year = fold$eval,
      best_iter =
        as.integer(
          selected$best_iteration
        )
    )

    pred_dt <- data.table(
      admission_year =
        dt$admission_year[eval_idx],
      inc_key =
        dt$inc_key[eval_idx],
      y =
        y_all[eval_idx],
      pred = p
    )
  }

  fwrite(
    pred_dt,
    file.path(
      PRED_DIR,
      paste0(
        endpoint,
        "__eval",
        fold$eval,
        "__predictions.csv"
      )
    )
  )

  saveRDS(
    list(
      model = fit,
      endpoint = endpoint,
      evaluation_year = fold$eval,
      predictor_policy = spec_ep$policy,
      predictors = predictors,
      encoder = enc_dev,
      feature_names = colnames(X_dev),
      selected_config = selected,
      compute_mode = compute$mode
    ),
    file.path(
      MODEL_DIR,
      paste0(
        endpoint,
        "__eval",
        fold$eval,
        "__model.rds"
      )
    ),
    compress = FALSE
  )

  qc <- data.table(
    endpoint = endpoint,
    predictor_policy = spec_ep$policy,
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
    selected_config = selected$config_id,
    selected_rounds =
      as.integer(
        selected$best_iteration
      ),
    tune_score =
      as.numeric(
        selected$tune_score
      ),
    hit_round_ceiling =
      as.integer(
        selected$best_iteration
      ) >=
        (MAX_ROUNDS - 5L)
  )

  rm(
    X_dev,
    X_eval,
    ddev,
    deval,
    fit,
    raw,
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
# 10. Run all 5 endpoints x 3 folds
# ------------------------------------------------------------------------------

run_plan <- CJ(
  endpoint = names(endpoint_specs),
  evaluation_year = as.integer(names(folds)),
  unique = TRUE
)

run_plan[
  ,
  endpoint_order :=
    match(
      endpoint,
      names(endpoint_specs)
    )
]
setorder(
  run_plan,
  endpoint_order,
  evaluation_year
)
run_plan[, endpoint_order := NULL]

fwrite(
  run_plan,
  file.path(OUT_DIR, "04_RUN_PLAN.csv")
)

metrics_rows <- list()
selected_rows <- list()
qc_rows <- list()
runtime_rows <- list()

for (i in seq_len(nrow(run_plan))) {
  endpoint <- run_plan$endpoint[i]
  year <- as.character(
    run_plan$evaluation_year[i]
  )

  cat(
    "\n============================================================\n",
    "[",
    i,
    "/",
    nrow(run_plan),
    "] FINAL CLASSIFICATION: ",
    endpoint,
    " | eval ",
    year,
    "\n============================================================\n",
    sep = ""
  )

  t0 <- Sys.time()

  res <- fit_one(
    endpoint,
    folds[[year]]
  )

  metrics_rows[[length(metrics_rows) + 1L]] <-
    res$metrics

  s <- copy(res$selected)
  s[
    ,
    `:=`(
      endpoint = endpoint,
      evaluation_year = as.integer(year)
    )
  ]
  selected_rows[[length(selected_rows) + 1L]] <- s

  qc_rows[[length(qc_rows) + 1L]] <- res$qc

  runtime_rows[[length(runtime_rows) + 1L]] <- data.table(
    endpoint = endpoint,
    evaluation_year = as.integer(year),
    runtime_minutes =
      as.numeric(
        difftime(
          Sys.time(),
          t0,
          units = "mins"
        )
      )
  )

  fwrite(
    rbindlist(metrics_rows, fill = TRUE),
    file.path(
      OUT_DIR,
      "05_FINAL_CLASSIFICATION_METRICS_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(selected_rows, fill = TRUE),
    file.path(
      OUT_DIR,
      "14_FINAL_SELECTED_CONFIGS_BY_FOLD_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(qc_rows, fill = TRUE),
    file.path(
      OUT_DIR,
      "15_FINAL_TRAINING_QC_CHECKPOINT.csv"
    )
  )
}

metrics <- rbindlist(
  metrics_rows,
  fill = TRUE
)
selected <- rbindlist(
  selected_rows,
  fill = TRUE
)
training_qc <- rbindlist(
  qc_rows,
  fill = TRUE
)
runtime <- rbindlist(
  runtime_rows,
  fill = TRUE
)

# Selected config must be unique for all 15 endpoint/fold combinations.
selected_key_qc <- selected[
  ,
  .N,
  by = .(
    endpoint,
    evaluation_year
  )
]

if (
  nrow(selected_key_qc) != 15L ||
    any(selected_key_qc$N != 1L)
) {
  fwrite(
    selected_key_qc,
    file.path(
      OUT_DIR,
      "13_SELECTED_CONFIG_KEY_QC_FAILED.csv"
    )
  )
  stop(
    "Selected config uniqueness QC failed.",
    call. = FALSE
  )
}

# Headline rows: 4 multiclass macro x 3 folds + surgery binary x 3 = 15.
headline <- metrics[
  metric_level %in%
    c(
      "macro",
      "binary"
    )
]

headline_qc <- headline[
  ,
  .N,
  by = .(
    endpoint,
    evaluation_year
  )
]

if (
  nrow(headline_qc) != 15L ||
    any(headline_qc$N != 1L)
) {
  fwrite(
    headline_qc,
    file.path(
      OUT_DIR,
      "13_HEADLINE_KEY_QC_FAILED.csv"
    )
  )
  stop(
    "Headline metric uniqueness QC failed.",
    call. = FALSE
  )
}

fwrite(
  metrics,
  file.path(
    OUT_DIR,
    "05_FINAL_CLASSIFICATION_METRICS.csv"
  )
)

fwrite(
  headline,
  file.path(
    OUT_DIR,
    "06_FINAL_CLASSIFICATION_HEADLINE_METRICS.csv"
  )
)

fwrite(
  training_qc,
  file.path(
    OUT_DIR,
    "07_FINAL_CLASSIFICATION_TRAINING_QC.csv"
  )
)

fwrite(
  runtime,
  file.path(
    OUT_DIR,
    "08_FINAL_CLASSIFICATION_RUNTIME.csv"
  )
)

# Preserve downstream-compatible filename.
fwrite(
  selected,
  file.path(
    OUT_DIR,
    "14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv"
  )
)

# Compact 2024 class-specific table.
metrics_2024 <- metrics[
  evaluation_year == 2024L
]

fwrite(
  metrics_2024,
  file.path(
    OUT_DIR,
    "09_2024_FINAL_CLASSIFICATION_METRICS.csv"
  )
)

# Temporal range summary.
temporal_summary <- headline[
  ,
  .(
    min_auroc =
      min(auroc, na.rm = TRUE),
    median_auroc =
      median(auroc, na.rm = TRUE),
    max_auroc =
      max(auroc, na.rm = TRUE),
    min_auprc =
      min(auprc, na.rm = TRUE),
    median_auprc =
      median(auprc, na.rm = TRUE),
    max_auprc =
      max(auprc, na.rm = TRUE),
    min_calibration_slope =
      min(calibration_slope, na.rm = TRUE),
    median_calibration_slope =
      median(calibration_slope, na.rm = TRUE),
    max_calibration_slope =
      max(calibration_slope, na.rm = TRUE)
  ),
  by = endpoint
]

fwrite(
  temporal_summary,
  file.path(
    OUT_DIR,
    "10_FINAL_CLASSIFICATION_TEMPORAL_SUMMARY.csv"
  )
)

# Predictor lock note.
lock_note <- data.table(
  item = c(
    "Primary population",
    "final_core",
    "final_context",
    "Disposition predictor policy",
    "HLOS trajectory predictor policy",
    "ICU trajectory predictor policy",
    "Ventilation trajectory predictor policy",
    "Operation predictor policy",
    "TL soft tissue",
    "Social context",
    "2024 interpretation"
  ),
  specification = c(
    "Direct presentations only; repaired frozen primary N=424,017",
    "54 predictors: coarse spine including TL soft tissue + detailed chest + pupil",
    "57 predictors: final_core + race + ethnicity + payer",
    "final_context",
    "final_context",
    "final_core",
    "final_core",
    "final_core",
    "Retained after repaired-DX Stage 2C",
    "Disposition/HLOS only; contextual rather than biological",
    "Development-informed forward temporal evaluation"
  )
)

fwrite(
  lock_note,
  file.path(
    OUT_DIR,
    "11_FINAL_CLASSIFICATION_LOCK_NOTE.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT FINAL REPAIRED-DX CLASSIFICATION COMPLETE",
    "",
    "Predictor policy:",
    "  final_core = 54",
    "  final_context = 57",
    "",
    "Temporal evaluation:",
    "  2022, 2023, 2024 rolling-origin folds",
    "",
    "Each endpoint/fold was fully tuned across the prespecified 24-config structural grid.",
    "No evaluation-year recalibration was performed.",
    "",
    "Review first:",
    "  06_FINAL_CLASSIFICATION_HEADLINE_METRICS.csv",
    "  09_2024_FINAL_CLASSIFICATION_METRICS.csv",
    "  10_FINAL_CLASSIFICATION_TEMPORAL_SUMMARY.csv",
    "  14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv",
    "  07_FINAL_CLASSIFICATION_TRAINING_QC.csv",
    "",
    "After this passes, proceed to repaired-DX duration quantiles and ridge comparator."
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
  "FINAL REPAIRED-DX CLASSIFICATION COMPLETE\n"
)
cat(
  "============================================================\n"
)
cat(
  "final_core predictors: ",
  length(final_core),
  "\n",
  sep = ""
)
cat(
  "final_context predictors: ",
  length(final_context),
  "\n\n",
  sep = ""
)
print(
  headline[
    order(
      endpoint,
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
