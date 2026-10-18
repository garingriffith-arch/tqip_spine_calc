# ==============================================================================
# 09s_extend_SPINE_TRACT_ridge_vent2024_lambda_only_v1.R
# SPINE-TRACT repaired-DX Ridge comparator — TARGETED boundary continuation
#
# PURPOSE
#   Continue ONLY the one unresolved Ridge tuning path:
#     endpoint = vent3_final
#     evaluation year = 2024
#
# WHY
#   v1.2 successfully forced all 250 requested glmnet lambdas to be evaluated.
#   Fourteen of 15 endpoint/fold fits had an interior optimum. Only vent3_final
#   / eval2024 selected the smallest available lambda.
#
#   The observed v1.2 floor was exactly lambda_max * 1e-6 because glmnet's
#   internal eps control limits automatically generated paths to that ratio.
#   This script avoids that limitation by supplying an EXPLICIT decreasing
#   lambda continuation from the v1.2 floor to lambda_max * 1e-10.
#
# DESIGN
#   Training for lambda path: 2020-2022
#   Tuning year:              2023
#   Final refit:              2020-2023
#   Evaluation:               2024
#
# OUTPUT
#   Replaces ONLY vent3_final/eval2024 in copies of the v1.2 Ridge metrics/QC,
#   then builds the final XGBoost-vs-Ridge comparator.
#
#   No other Ridge fit is rerun.
# ==============================================================================

rm(list = ls())
gc()

required <- c(
  "data.table",
  "Matrix",
  "glmnet"
)

missing <- required[
  !vapply(
    required,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
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
  library(glmnet)
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

XGB_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification_REPAIRED_DX",
  "rolling_origin_gpu"
)

XGB_METRICS_FILE <- file.path(
  XGB_DIR,
  "05_FINAL_CLASSIFICATION_METRICS.csv"
)

RIDGE_V12_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_ridge_comparator_REPAIRED_DX_v1_2"
)

RIDGE_V12_METRICS <- file.path(
  RIDGE_V12_DIR,
  "01_RIDGE_TEMPORAL_METRICS.csv"
)

RIDGE_V12_QC <- file.path(
  RIDGE_V12_DIR,
  "02_RIDGE_TRAINING_QC.csv"
)

UNRESOLVED_FILE <- file.path(
  RIDGE_V12_DIR,
  "02b_UNRESOLVED_LAMBDA_BOUNDARY.csv"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_ridge_comparator_REPAIRED_DX_FINAL"
)

PRED_DIR <- file.path(
  OUT_DIR,
  "predictions"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  PRED_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

needed <- c(
  DATA_FILE,
  SPEC_FILE,
  XGB_METRICS_FILE,
  RIDGE_V12_METRICS,
  RIDGE_V12_QC,
  UNRESOLVED_FILE
)

if (!all(file.exists(needed))) {
  stop(
    "Missing prerequisite file(s):\n",
    paste(
      needed[!file.exists(needed)],
      collapse = "\n"
    ),
    call. = FALSE
  )
}

SEED <- 20261002L
TARGET_ENDPOINT <- "vent3_final"
TARGET_YEAR <- 2024L
EXTENDED_FLOOR_RATIO <- 1e-10
N_EXTENSION_LAMBDAS <- 120L

logical_cores <- parallel::detectCores(logical = TRUE)

if (
  is.na(logical_cores) ||
    logical_cores < 2L
) {
  logical_cores <- 2L
}

N_THREADS <- max(
  1L,
  min(
    logical_cores - 1L,
    as.integer(
      round(
        logical_cores *
          0.80
      )
    )
  )
)

setDTthreads(
  N_THREADS
)

# Prevent glmnet from terminating an explicitly supplied path early.
glmnet::glmnet.control(
  fdev = 0,
  devmax = 1
)

# ------------------------------------------------------------------------------
# 1. Verify that v1.2 has EXACTLY the expected single unresolved fold
# ------------------------------------------------------------------------------

old_metrics <- fread(
  RIDGE_V12_METRICS
)

old_qc <- fread(
  RIDGE_V12_QC
)

unresolved <- fread(
  UNRESOLVED_FILE
)

if (
  nrow(unresolved) != 1L ||
    unresolved$endpoint[1L] !=
      TARGET_ENDPOINT ||
    unresolved$evaluation_year[1L] !=
      TARGET_YEAR
) {
  print(unresolved)

  stop(
    "Expected exactly one unresolved boundary: vent3_final / eval2024.",
    call. = FALSE
  )
}

if (
  nrow(old_qc) != 15L ||
    !all(
      old_qc$full_requested_lambda_path_completed ==
        TRUE
    )
) {
  stop(
    "v1.2 QC is not the expected completed 15-fold run.",
    call. = FALSE
  )
}

target_qc_old <- old_qc[
  endpoint ==
    TARGET_ENDPOINT &
    evaluation_year ==
      TARGET_YEAR
]

if (
  nrow(target_qc_old) != 1L ||
    target_qc_old$selected_lambda_at_path_boundary !=
      TRUE
) {
  stop(
    "Target v1.2 QC row is not the expected boundary case.",
    call. = FALSE
  )
}

lambda_max <- as.numeric(
  target_qc_old$lambda_path_max
)

old_floor <- as.numeric(
  target_qc_old$lambda_path_min
)

if (
  !is.finite(lambda_max) ||
    !is.finite(old_floor) ||
    lambda_max <= 0 ||
    old_floor <= 0
) {
  stop(
    "Invalid v1.2 lambda bounds.",
    call. = FALSE
  )
}

old_floor_ratio <- old_floor /
  lambda_max

cat(
  "\nConfirmed only unresolved fold:\n",
  "  endpoint: ",
  TARGET_ENDPOINT,
  "\n  evaluation: ",
  TARGET_YEAR,
  "\n  lambda_max: ",
  signif(lambda_max, 8),
  "\n  old floor: ",
  signif(old_floor, 8),
  "\n  old floor ratio: ",
  format(
    old_floor_ratio,
    scientific = TRUE
  ),
  "\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 2. Load data and reconstruct corrected final_core
# ------------------------------------------------------------------------------

dt <- as.data.table(
  readRDS(
    DATA_FILE
  )
)

spec <- readRDS(
  SPEC_FILE
)

xgb_metrics <- fread(
  XGB_METRICS_FILE
)

dt[
  ,
  inc_key :=
    as.character(
      inc_key
    )
]

dt[
  ,
  admission_year :=
    as.integer(
      admission_year
    )
]

if (nrow(dt) != 424017L) {
  stop(
    "Expected repaired primary N=424017.",
    call. = FALSE
  )
}

pg <- spec$predictor_groups
clinical_backbone <-
  spec$predictor_sets$clinical_backbone

final_core <- unique(c(
  setdiff(
    clinical_backbone,
    pg$spine_granular
  ),
  pg$spine_coarse,
  pg$chest_detail,
  pg$pupil
))

if (
  length(final_core) != 54L ||
    !"dx_tl_soft_tissue_any" %in%
      final_core
) {
  stop(
    "Corrected 54-predictor final_core reconstruction failed.",
    call. = FALSE
  )
}

dt[
  ,
  vent3_final :=
    factor(
      as.character(
        vent3
      ),
      levels = c(
        "No ventilation",
        "Ventilation 1-7 days",
        "Ventilation >=8 days"
      )
    )
]

levels_target <- levels(
  dt$vent3_final
)

# ------------------------------------------------------------------------------
# 3. Same v1.2 Ridge encoder
# ------------------------------------------------------------------------------

is_num <- function(x) {
  is.numeric(x) ||
    is.integer(x)
}

clean_cat <- function(x) {
  z <- trimws(
    as.character(x)
  )

  z[
    is.na(z) |
      z == "" |
      z %in%
        c(
          "NA",
          "NaN",
          "<NA>"
        )
  ] <- "__UNKNOWN__"

  z
}

fit_ridge_encoder <- function(
    d,
    predictors
) {
  enc <- vector(
    "list",
    length(predictors)
  )

  names(enc) <- predictors

  for (v in predictors) {
    x <- d[[v]]

    if (is_num(x)) {
      med <- median(
        as.numeric(x),
        na.rm = TRUE
      )

      if (!is.finite(med)) {
        med <- 0
      }

      enc[[v]] <- list(
        type = "numeric",
        median = med
      )
    } else {
      z <- clean_cat(x)

      lev <- unique(c(
        sort(
          unique(z)
        ),
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

encode_ridge <- function(
    d,
    predictors,
    enc
) {
  n <- nrow(d)
  blocks <- list()
  block_i <- 0L

  for (v in predictors) {
    info <- enc[[v]]
    x <- d[[v]]

    if (info$type == "numeric") {
      z <- as.numeric(x)

      miss <- as.integer(
        is.na(z) |
          !is.finite(z)
      )

      z[
        is.na(z) |
          !is.finite(z)
      ] <- info$median

      block_i <- block_i + 1L

      m <- Matrix(
        z,
        nrow = n,
        ncol = 1L,
        sparse = TRUE
      )

      colnames(m) <- v
      blocks[[block_i]] <- m

      block_i <- block_i + 1L

      mm <- Matrix(
        miss,
        nrow = n,
        ncol = 1L,
        sparse = TRUE
      )

      colnames(mm) <- paste0(
        v,
        "__MISSING"
      )

      blocks[[block_i]] <- mm
    } else {
      z <- clean_cat(x)

      z[
        !z %in%
          info$levels
      ] <- "__OTHER__"

      f <- factor(
        z,
        levels = info$levels
      )

      m <- sparseMatrix(
        i = seq_len(n),
        j = as.integer(f),
        x = 1,
        dims = c(
          n,
          length(info$levels)
        )
      )

      colnames(m) <- paste0(
        v,
        "__",
        make.names(
          info$levels,
          unique = TRUE
        )
      )

      block_i <- block_i + 1L
      blocks[[block_i]] <- m
    }
  }

  do.call(
    cbind,
    blocks
  )
}

# ------------------------------------------------------------------------------
# 4. Metric helpers
# ------------------------------------------------------------------------------

clamp_prob <- function(
    p,
    eps = 1e-7
) {
  pmin(
    pmax(
      as.numeric(p),
      eps
    ),
    1 - eps
  )
}

fast_auc <- function(y, p) {
  ok <- !is.na(y) &
    !is.na(p) &
    is.finite(p)

  y <- as.integer(
    y[ok]
  )

  p <- as.numeric(
    p[ok]
  )

  n1 <- as.double(
    sum(
      y == 1L
    )
  )

  n0 <- as.double(
    sum(
      y == 0L
    )
  )

  if (
    n1 == 0 ||
      n0 == 0
  ) {
    return(NA_real_)
  }

  r <- rank(
    p,
    ties.method = "average"
  )

  (
    sum(
      r[y == 1L]
    ) -
      n1 *
        (
          n1 + 1
        ) /
        2
  ) /
    (
      n1 *
        n0
    )
}

fast_auprc <- function(y, p) {
  ok <- !is.na(y) &
    !is.na(p) &
    is.finite(p)

  y <- as.integer(
    y[ok]
  )

  p <- as.numeric(
    p[ok]
  )

  npos <- sum(
    y == 1L
  )

  if (npos == 0L) {
    return(NA_real_)
  }

  ord <- order(
    p,
    decreasing = TRUE
  )

  ys <- y[ord]

  precision <- cumsum(
    ys == 1L
  ) /
    seq_along(ys)

  mean(
    precision[
      ys == 1L
    ]
  )
}

binary_calibration <- function(y, p) {
  pp <- clamp_prob(p)
  lp <- qlogis(pp)

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
      if (is.null(fit_i)) {
        NA_real_
      } else {
        unname(
          coef(fit_i)[1L]
        )
      },
    slope =
      if (is.null(fit_s)) {
        NA_real_
      } else {
        unname(
          coef(fit_s)[2L]
        )
      }
  )
}

multiclass_metric_rows <- function(
    y_factor,
    prob,
    levels,
    endpoint,
    eval_year
) {
  y_index <- as.integer(
    factor(
      y_factor,
      levels = levels
    )
  )

  class_rows <- rbindlist(
    lapply(
      seq_along(levels),
      function(j) {
        y <- as.integer(
          y_index == j
        )

        p <- prob[, j]
        pp <- clamp_prob(p)
        prev <- mean(
          y == 1L
        )
        pr <- fast_auprc(
          y,
          p
        )
        cal <- binary_calibration(
          y,
          p
        )

        data.table(
          endpoint = endpoint,
          evaluation_year = eval_year,
          metric_level = "class",
          class = levels[j],
          n = length(y),
          prevalence = prev,
          auroc =
            fast_auc(
              y,
              p
            ),
          auprc = pr,
          pr_lift =
            if (
              is.finite(prev) &&
                prev > 0
            ) {
              pr /
                prev
            } else {
              NA_real_
            },
          brier =
            mean(
              (
                y -
                  p
              )^2
            ),
          logloss =
            -mean(
              y *
                log(pp) +
                (
                  1 -
                    y
                ) *
                  log(
                    1 -
                      pp
                  )
            ),
          calibration_intercept =
            unname(
              cal["intercept"]
            ),
          calibration_slope =
            unname(
              cal["slope"]
            ),
          multiclass_logloss =
            NA_real_,
          multiclass_brier =
            NA_real_
        )
      }
    )
  )

  ptrue <- pmax(
    prob[
      cbind(
        seq_along(y_index),
        y_index
      )
    ],
    1e-15
  )

  mc_logloss <- -mean(
    log(ptrue)
  )

  onehot <- matrix(
    0,
    nrow = length(y_index),
    ncol = length(levels)
  )

  onehot[
    cbind(
      seq_along(y_index),
      y_index
    )
  ] <- 1

  mc_brier <- mean(
    rowSums(
      (
        prob -
          onehot
      )^2
    )
  )

  class_rows[
    ,
    `:=`(
      multiclass_logloss =
        mc_logloss,
      multiclass_brier =
        mc_brier
    )
  ]

  macro <- data.table(
    endpoint = endpoint,
    evaluation_year = eval_year,
    metric_level = "macro",
    class = "MACRO",
    n = length(y_index),
    prevalence = NA_real_,
    auroc =
      mean(
        class_rows$auroc,
        na.rm = TRUE
      ),
    auprc =
      mean(
        class_rows$auprc,
        na.rm = TRUE
      ),
    pr_lift =
      mean(
        class_rows$pr_lift,
        na.rm = TRUE
      ),
    brier =
      mean(
        class_rows$brier,
        na.rm = TRUE
      ),
    logloss =
      mean(
        class_rows$logloss,
        na.rm = TRUE
      ),
    calibration_intercept =
      median(
        class_rows$calibration_intercept,
        na.rm = TRUE
      ),
    calibration_slope =
      median(
        class_rows$calibration_slope,
        na.rm = TRUE
      ),
    multiclass_logloss =
      mc_logloss,
    multiclass_brier =
      mc_brier
  )

  rbindlist(
    list(
      class_rows,
      macro
    )
  )
}

# ------------------------------------------------------------------------------
# 5. Targeted lambda continuation: vent3_final / eval2024 only
# ------------------------------------------------------------------------------

y_all <- dt$vent3_final

train_idx <- which(
  dt$admission_year %in%
    2020:2022 &
    !is.na(y_all)
)

tune_idx <- which(
  dt$admission_year ==
    2023L &
    !is.na(y_all)
)

dev_idx <- which(
  dt$admission_year %in%
    2020:2023 &
    !is.na(y_all)
)

eval_idx <- which(
  dt$admission_year ==
    2024L &
    !is.na(y_all)
)

enc_train <- fit_ridge_encoder(
  dt[train_idx],
  final_core
)

X_train <- encode_ridge(
  dt[train_idx],
  final_core,
  enc_train
)

X_tune <- encode_ridge(
  dt[tune_idx],
  final_core,
  enc_train
)

# Explicit continuation starts at the old floor and extends four additional
# orders of magnitude to lambda_max * 1e-10.
extended_floor <-
  lambda_max *
    EXTENDED_FLOOR_RATIO

lambda_extension <- exp(
  seq(
    log(old_floor),
    log(extended_floor),
    length.out =
      N_EXTENSION_LAMBDAS
  )
)

# Guarantee strict decreasing order and no accidental duplicate.
lambda_extension <- sort(
  unique(
    as.numeric(
      lambda_extension
    )
  ),
  decreasing = TRUE
)

if (
  min(lambda_extension) >
    lambda_max *
      EXTENDED_FLOOR_RATIO *
      1.01
) {
  stop(
    "Explicit continuation did not reach requested near-zero floor.",
    call. = FALSE
  )
}

set.seed(
  SEED +
    TARGET_YEAR
)

fit_extension <- glmnet(
  x = X_train,
  y =
    factor(
      y_all[train_idx],
      levels = levels_target
    ),
  family = "multinomial",
  alpha = 0,
  lambda =
    lambda_extension,
  standardize = TRUE,
  type.multinomial = "grouped"
)

actual_lambda <- as.numeric(
  fit_extension$lambda
)

if (
  length(actual_lambda) !=
    length(lambda_extension)
) {
  stop(
    "glmnet did not return every explicitly supplied continuation lambda.",
    call. = FALSE
  )
}

pred_tune <- predict(
  fit_extension,
  newx = X_tune,
  type = "response"
)

y_tune_index <- as.integer(
  factor(
    y_all[tune_idx],
    levels = levels_target
  )
)

losses <- vapply(
  seq_along(actual_lambda),
  function(k) {
    p <- pred_tune[
      ,
      ,
      k,
      drop = FALSE
    ]

    p <- matrix(
      p,
      nrow = length(tune_idx),
      ncol = length(levels_target)
    )

    ptrue <- pmax(
      p[
        cbind(
          seq_along(
            y_tune_index
          ),
          y_tune_index
        )
      ],
      1e-15
    )

    -mean(
      log(ptrue)
    )
  },
  numeric(1)
)

extension_table <- data.table(
  path_index =
    seq_along(
      actual_lambda
    ),
  lambda =
    actual_lambda,
  lambda_ratio_to_original_max =
    actual_lambda /
      lambda_max,
  tuning_logloss =
    losses
)

fwrite(
  extension_table,
  file.path(
    OUT_DIR,
    "01_VENT2024_LAMBDA_CONTINUATION.csv"
  )
)

best_j <- which.min(
  losses
)

selected_lambda <-
  actual_lambda[best_j]

selected_ratio <-
  selected_lambda /
    lambda_max

selected_at_extension_floor <-
  best_j ==
    length(actual_lambda)

cat(
  "\nTargeted continuation result:\n",
  "  old tuning logloss: ",
  signif(
    target_qc_old$tuning_logloss,
    10
  ),
  "\n  new tuning logloss: ",
  signif(
    losses[best_j],
    10
  ),
  "\n  selected lambda: ",
  signif(
    selected_lambda,
    10
  ),
  "\n  selected ratio to original lambda_max: ",
  format(
    selected_ratio,
    scientific = TRUE
  ),
  "\n  selected at extension floor: ",
  selected_at_extension_floor,
  "\n",
  sep = ""
)

rm(
  fit_extension,
  pred_tune,
  X_train,
  X_tune,
  enc_train
)
gc()

# ------------------------------------------------------------------------------
# 6. Refit target on 2020-2023 and evaluate 2024
# ------------------------------------------------------------------------------

enc_dev <- fit_ridge_encoder(
  dt[dev_idx],
  final_core
)

X_dev <- encode_ridge(
  dt[dev_idx],
  final_core,
  enc_dev
)

X_eval <- encode_ridge(
  dt[eval_idx],
  final_core,
  enc_dev
)

final_fit <- glmnet(
  x = X_dev,
  y =
    factor(
      y_all[dev_idx],
      levels = levels_target
    ),
  family = "multinomial",
  alpha = 0,
  lambda =
    selected_lambda,
  standardize = TRUE,
  type.multinomial = "grouped"
)

pred <- predict(
  final_fit,
  newx = X_eval,
  type = "response",
  s = selected_lambda
)

prob <- matrix(
  pred,
  nrow = length(eval_idx),
  ncol = length(levels_target)
)

colnames(prob) <- levels_target

target_metrics <- multiclass_metric_rows(
  y_factor =
    y_all[eval_idx],
  prob = prob,
  levels = levels_target,
  endpoint = TARGET_ENDPOINT,
  eval_year = TARGET_YEAR
)

pred_dt <- data.table(
  admission_year =
    dt$admission_year[eval_idx],
  inc_key =
    dt$inc_key[eval_idx],
  observed_class =
    as.character(
      y_all[eval_idx]
    )
)

for (j in seq_along(levels_target)) {
  pred_dt[
    ,
    paste0(
      "p_",
      make.names(
        levels_target[j]
      )
    ) :=
      prob[, j]
  ]
}

fwrite(
  pred_dt,
  file.path(
    PRED_DIR,
    "vent3_final__eval2024__ridge_predictions.csv"
  )
)

# Replace ONLY target metric rows in v1.2.
final_ridge_metrics <- rbindlist(
  list(
    old_metrics[
      !(
        endpoint ==
          TARGET_ENDPOINT &
          evaluation_year ==
            TARGET_YEAR
      )
    ],
    target_metrics
  ),
  fill = TRUE,
  use.names = TRUE
)

setorder(
  final_ridge_metrics,
  endpoint,
  evaluation_year,
  metric_level,
  class
)

target_qc_new <- copy(
  target_qc_old
)

target_qc_new[
  ,
  `:=`(
    selected_lambda =
      selected_lambda,
    tuning_logloss =
      losses[best_j],
    lambda_path_index =
      best_j,
    lambda_path_n =
      length(actual_lambda),
    lambda_path_max =
      lambda_max,
    lambda_path_min =
      min(actual_lambda),
    selected_lambda_ratio_to_max =
      selected_ratio,
    selected_lambda_at_path_boundary =
      selected_at_extension_floor,
    selected_at_near_zero_floor =
      selected_ratio <=
        1.1 *
          EXTENDED_FLOOR_RATIO,
    full_requested_lambda_path_completed =
      TRUE,
    n_encoded_features =
      ncol(X_dev)
  )
]

target_qc_new[
  ,
  continuation_from_v12 :=
    TRUE
]

final_qc <- copy(
  old_qc[
    !(
      endpoint ==
        TARGET_ENDPOINT &
        evaluation_year ==
          TARGET_YEAR
    )
  ]
)

final_qc[
  ,
  continuation_from_v12 :=
    FALSE
]

final_qc <- rbindlist(
  list(
    final_qc,
    target_qc_new
  ),
  fill = TRUE,
  use.names = TRUE
)

setorder(
  final_qc,
  endpoint,
  evaluation_year
)

fwrite(
  final_ridge_metrics,
  file.path(
    OUT_DIR,
    "02_FINAL_RIDGE_TEMPORAL_METRICS.csv"
  )
)

fwrite(
  final_qc,
  file.path(
    OUT_DIR,
    "03_FINAL_RIDGE_TRAINING_QC.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Final fair XGBoost-vs-Ridge comparison
# ------------------------------------------------------------------------------

join_keys <- c(
  "endpoint",
  "evaluation_year",
  "metric_level",
  "class"
)

xgb_keep <- intersect(
  c(
    join_keys,
    "n",
    "prevalence",
    "auroc",
    "auprc",
    "pr_lift",
    "brier",
    "logloss",
    "calibration_intercept",
    "calibration_slope",
    "multiclass_logloss",
    "multiclass_brier"
  ),
  names(xgb_metrics)
)

xgb_for_compare <-
  xgb_metrics[
    ,
    ..xgb_keep
  ]

setnames(
  xgb_for_compare,
  setdiff(
    names(xgb_for_compare),
    join_keys
  ),
  paste0(
    "xgb_",
    setdiff(
      names(xgb_for_compare),
      join_keys
    )
  )
)

ridge_for_compare <- copy(
  final_ridge_metrics
)

setnames(
  ridge_for_compare,
  setdiff(
    names(ridge_for_compare),
    join_keys
  ),
  paste0(
    "ridge_",
    setdiff(
      names(ridge_for_compare),
      join_keys
    )
  )
)

comparison <- merge(
  xgb_for_compare,
  ridge_for_compare,
  by = join_keys,
  all = TRUE,
  sort = FALSE
)

comparison[, `:=`(
  delta_auroc_xgb_minus_ridge =
    xgb_auroc -
      ridge_auroc,
  delta_auprc_xgb_minus_ridge =
    xgb_auprc -
      ridge_auprc,
  delta_brier_xgb_minus_ridge =
    xgb_brier -
      ridge_brier,
  delta_logloss_xgb_minus_ridge =
    xgb_logloss -
      ridge_logloss
)]

fwrite(
  comparison,
  file.path(
    OUT_DIR,
    "04_XGBOOST_VS_RIDGE_TEMPORAL_COMPARISON.csv"
  )
)

cohort_qc <- comparison[
  !is.na(xgb_n) &
    !is.na(ridge_n),
  .(
    rows_compared = .N,
    n_mismatch_rows =
      sum(
        xgb_n !=
          ridge_n
      ),
    prevalence_mismatch_rows =
      sum(
        abs(
          xgb_prevalence -
            ridge_prevalence
        ) >
          1e-12,
        na.rm = TRUE
      )
  )
]

fwrite(
  cohort_qc,
  file.path(
    OUT_DIR,
    "05_RIDGE_XGBOOST_COHORT_CONCORDANCE_QC.csv"
  )
)

if (
  cohort_qc$n_mismatch_rows >
    0L ||
    cohort_qc$prevalence_mismatch_rows >
      0L
) {
  stop(
    "Final Ridge/XGBoost cohort concordance failed.",
    call. = FALSE
  )
}

headline <- comparison[
  metric_level %in%
    c(
      "macro",
      "binary"
    )
]

fwrite(
  headline,
  file.path(
    OUT_DIR,
    "06_FINAL_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv"
  )
)

# Compact boundary-resolution record.
resolution <- data.table(
  endpoint =
    TARGET_ENDPOINT,
  evaluation_year =
    TARGET_YEAR,
  original_lambda_max =
    lambda_max,
  v12_lambda_floor =
    old_floor,
  v12_floor_ratio =
    old_floor_ratio,
  extension_floor =
    min(actual_lambda),
  extension_floor_ratio =
    min(actual_lambda) /
      lambda_max,
  selected_lambda =
    selected_lambda,
  selected_lambda_ratio =
    selected_ratio,
  selected_at_extension_floor =
    selected_at_extension_floor,
  v12_tuning_logloss =
    target_qc_old$tuning_logloss,
  final_tuning_logloss =
    losses[best_j],
  delta_tuning_logloss =
    losses[best_j] -
      target_qc_old$tuning_logloss
)

fwrite(
  resolution,
  file.path(
    OUT_DIR,
    "07_VENT2024_BOUNDARY_RESOLUTION.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT FINAL REPAIRED-DX RIDGE COMPARATOR COMPLETE",
    "",
    "Fourteen of 15 endpoint/fold Ridge models use the fully traversed v1.2 paths.",
    "The only unresolved v1.2 fold (vent3_final/eval2024) was continued explicitly below the glmnet automatic path floor to lambda_max*1e-10.",
    "No other model was rerun.",
    "",
    "A floor selection after this explicit continuation is interpreted as an effectively near-unpenalized Ridge optimum, not an accidentally truncated lambda path.",
    "",
    "Review:",
    "  06_FINAL_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv",
    "  05_RIDGE_XGBOOST_COHORT_CONCORDANCE_QC.csv",
    "  03_FINAL_RIDGE_TRAINING_QC.csv",
    "  07_VENT2024_BOUNDARY_RESOLUTION.csv"
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
  "TARGETED RIDGE BOUNDARY CONTINUATION COMPLETE\n"
)
cat(
  "============================================================\n\n"
)
print(
  resolution
)
cat(
  "\n2024 headline comparator:\n"
)
print(
  headline[
    evaluation_year ==
      2024L
  ]
)
cat(
  "\nOutput:\n  ",
  OUT_DIR,
  "\n",
  sep = ""
)
