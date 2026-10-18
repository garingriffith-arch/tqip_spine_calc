# ==============================================================================
# 09r_run_SPINE_TRACT_fair_ridge_REPAIRED_DX_v1_2.R
# SPINE-TRACT Stage 3B-R v1.2 — fair rolling-origin Ridge comparator
#
# PURPOSE
#   Fit a conventional penalized-regression comparator using the SAME corrected
#   endpoint-specific predictor policy and SAME temporal folds as the final
#   repaired-DX XGBoost models.
#
# PREDICTOR POLICY
#   final_context (57): discharge and HLOS trajectory
#   final_core    (54): ICU, ventilation, recorded TL operation
#
# RIDGE PREPROCESSING
#   - numeric: train-period median imputation + explicit missingness indicator
#   - categorical: explicit __UNKNOWN__; unseen levels -> __OTHER__; one-hot
#   - alpha = 0
#   - lambda selected on the pre-evaluation temporal tuning year
#   - refit on all pre-evaluation development years at selected lambda
#
# MODEL FAMILY
#   multinomial Ridge for 3-class trajectories
#   binomial Ridge for recorded TL operation
#
# IMPORTANT
#   This is a comparator, not a second model-selection exercise.
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

# glmnet ordinarily truncates a lambda path early when fractional deviance
# improvement becomes very small. That caused v1.1 to return only 172-215 of
# the requested 250 lambdas and left most temporal folds at the path minimum.
# Disable that INTERNAL path truncation for this comparator so the prespecified
# lambda range is actually evaluated.
glmnet::glmnet.control(
  fdev = 0,
  devmax = 1
)

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

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_ridge_comparator_REPAIRED_DX_v1_2"
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
  XGB_METRICS_FILE
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

setDTthreads(N_THREADS)

# ------------------------------------------------------------------------------
# 1. Data / locked predictor sets
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
    "Expected repaired primary N=424017; observed ",
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

if (
  length(final_core) != 54L ||
    length(final_context) != 57L
) {
  stop(
    "Corrected predictor lock mismatch: final_core=",
    length(final_core),
    ", final_context=",
    length(final_context),
    ".",
    call. = FALSE
  )
}

# Final 3-class disposition.
dt[
  ,
  discharge3_final :=
    factor(
      fcase(
        as.character(discharge4) ==
          "Home/home health",
        "Home/home health",
        as.character(discharge4) %in%
          c(
            "Inpatient rehabilitation",
            "Other post-acute facility"
          ),
        "Post-acute facility",
        as.character(discharge4) ==
          "Death/hospice",
        "Death/hospice",
        default = NA_character_
      ),
      levels = c(
        "Home/home health",
        "Post-acute facility",
        "Death/hospice"
      )
    )
]

dt[
  ,
  hlos3_final :=
    factor(
      as.character(hlos3),
      levels = c(
        "<=7 days",
        "8-19 days",
        ">=20 days"
      )
    )
]

dt[
  ,
  icu3_final :=
    factor(
      as.character(icu3),
      levels = c(
        "No ICU",
        "ICU 1-7 days",
        "ICU >=8 days"
      )
    )
]

dt[
  ,
  vent3_final :=
    factor(
      as.character(vent3),
      levels = c(
        "No ventilation",
        "Ventilation 1-7 days",
        "Ventilation >=8 days"
      )
    )
]

endpoint_specs <- list(
  discharge3_final = list(
    type = "multinomial",
    levels = levels(
      dt$discharge3_final
    ),
    predictors = final_context,
    policy = "final_context"
  ),
  hlos3_final = list(
    type = "multinomial",
    levels = levels(
      dt$hlos3_final
    ),
    predictors = final_context,
    policy = "final_context"
  ),
  icu3_final = list(
    type = "multinomial",
    levels = levels(
      dt$icu3_final
    ),
    predictors = final_core,
    policy = "final_core"
  ),
  vent3_final = list(
    type = "multinomial",
    levels = levels(
      dt$vent3_final
    ),
    predictors = final_core,
    policy = "final_core"
  ),
  operative_spine_intervention = list(
    type = "binomial",
    predictors = final_core,
    policy = "final_core"
  )
)

# ------------------------------------------------------------------------------
# 2. Ridge encoder
# ------------------------------------------------------------------------------

is_num <- function(x) {
  is.numeric(x) || is.integer(x)
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
      z[!z %in% info$levels] <-
        "__OTHER__"

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
# 3. Metrics
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

  if (
    length(
      unique(y)
    ) <
      2L
  ) {
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

binary_metric_row <- function(
    y,
    p,
    endpoint,
    eval_year
) {
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
    metric_level = "binary",
    class = "Event",
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
        pr / prev
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
    multiclass_logloss = NA_real_,
    multiclass_brier = NA_real_
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
              pr / prev
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
          multiclass_logloss = NA_real_,
          multiclass_brier = NA_real_
        )
      }
    )
  )

  ptrue <- pmax(
    prob[
      cbind(
        seq_along(
          y_index
        ),
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
# 4. Temporal folds
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
# 5. Fit one endpoint / fold
# ------------------------------------------------------------------------------

fit_one <- function(
    endpoint,
    fold
) {
  ep <- endpoint_specs[[endpoint]]
  predictors <- ep$predictors

  if (ep$type == "multinomial") {
    y_all <- factor(
      as.character(
        dt[[endpoint]]
      ),
      levels = ep$levels
    )
  } else {
    y_all <- as.integer(
      dt[[endpoint]]
    )
  }

  train_idx <- which(
    dt$admission_year %in%
      fold$train
  )
  tune_idx <- which(
    dt$admission_year %in%
      fold$tune
  )
  dev_idx <- which(
    dt$admission_year %in%
      fold$dev
  )
  eval_idx <- which(
    dt$admission_year %in%
      fold$eval
  )

  if (ep$type == "multinomial") {
    train_idx <- train_idx[
      !is.na(
        y_all[train_idx]
      )
    ]
    tune_idx <- tune_idx[
      !is.na(
        y_all[tune_idx]
      )
    ]
    dev_idx <- dev_idx[
      !is.na(
        y_all[dev_idx]
      )
    ]
    eval_idx <- eval_idx[
      !is.na(
        y_all[eval_idx]
      )
    ]
  } else {
    train_idx <- train_idx[
      !is.na(
        y_all[train_idx]
      ) &
        y_all[train_idx] %in%
          c(
            0L,
            1L
          )
    ]
    tune_idx <- tune_idx[
      !is.na(
        y_all[tune_idx]
      ) &
        y_all[tune_idx] %in%
          c(
            0L,
            1L
          )
    ]
    dev_idx <- dev_idx[
      !is.na(
        y_all[dev_idx]
      ) &
        y_all[dev_idx] %in%
          c(
            0L,
            1L
          )
    ]
    eval_idx <- eval_idx[
      !is.na(
        y_all[eval_idx]
      ) &
        y_all[eval_idx] %in%
          c(
            0L,
            1L
          )
    ]
  }

  # Temporal tuning encoder fit ONLY on training years.
  enc_train <- fit_ridge_encoder(
    dt[train_idx],
    predictors
  )

  X_train <- encode_ridge(
    dt[train_idx],
    predictors,
    enc_train
  )

  X_tune <- encode_ridge(
    dt[tune_idx],
    predictors,
    enc_train
  )

  if (!identical(
    colnames(X_train),
    colnames(X_tune)
  )) {
    stop(
      "Ridge tuning encoder mismatch for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  set.seed(
    SEED +
      fold$eval
  )

  if (ep$type == "multinomial") {
    fit_path <- glmnet(
      x = X_train,
      y =
        factor(
          y_all[train_idx],
          levels = ep$levels
        ),
      family = "multinomial",
      alpha = 0,
      standardize = TRUE,
      type.multinomial = "grouped",
      nlambda = 250L,
      lambda.min.ratio = 1e-8
    )

    lambda_path <- fit_path$lambda

    pred_tune <- predict(
      fit_path,
      newx = X_tune,
      type = "response"
    )

    y_tune_index <- as.integer(
      factor(
        y_all[tune_idx],
        levels = ep$levels
      )
    )

    losses <- vapply(
      seq_along(lambda_path),
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
          ncol = length(ep$levels)
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
  } else {
    fit_path <- glmnet(
      x = X_train,
      y =
        y_all[train_idx],
      family = "binomial",
      alpha = 0,
      standardize = TRUE,
      nlambda = 250L,
      lambda.min.ratio = 1e-8
    )

    lambda_path <- fit_path$lambda

    pred_tune <- predict(
      fit_path,
      newx = X_tune,
      type = "response"
    )

    yy <- y_all[tune_idx]

    losses <- vapply(
      seq_along(lambda_path),
      function(k) {
        p <- clamp_prob(
          pred_tune[, k]
        )

        -mean(
          yy *
            log(p) +
            (
              1 -
                yy
            ) *
              log(
                1 -
                  p
              )
        )
      },
      numeric(1)
    )
  }

  best_j <- which.min(
    losses
  )

  selected_lambda <-
    lambda_path[best_j]

  lambda_boundary <- best_j %in%
    c(
      1L,
      length(lambda_path)
    )

  lambda_path_max <- max(
    lambda_path,
    na.rm = TRUE
  )

  lambda_path_min <- min(
    lambda_path,
    na.rm = TRUE
  )

  selected_lambda_ratio_to_max <-
    selected_lambda /
      lambda_path_max

  selected_at_near_zero_floor <-
    best_j ==
      length(lambda_path) &&
      selected_lambda_ratio_to_max <=
        1.1e-8

  rm(
    fit_path,
    pred_tune,
    X_train,
    X_tune,
    enc_train
  )
  gc()

  # Refit encoder on all pre-evaluation years.
  enc_dev <- fit_ridge_encoder(
    dt[dev_idx],
    predictors
  )

  X_dev <- encode_ridge(
    dt[dev_idx],
    predictors,
    enc_dev
  )

  X_eval <- encode_ridge(
    dt[eval_idx],
    predictors,
    enc_dev
  )

  if (!identical(
    colnames(X_dev),
    colnames(X_eval)
  )) {
    stop(
      "Ridge refit encoder mismatch for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  if (ep$type == "multinomial") {
    final_fit <- glmnet(
      x = X_dev,
      y =
        factor(
          y_all[dev_idx],
          levels = ep$levels
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
      ncol = length(ep$levels)
    )

    colnames(prob) <- ep$levels

    metrics <- multiclass_metric_rows(
      y_factor =
        y_all[eval_idx],
      prob = prob,
      levels = ep$levels,
      endpoint = endpoint,
      eval_year = fold$eval
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

    for (j in seq_along(ep$levels)) {
      pred_dt[
        ,
        paste0(
          "p_",
          make.names(
            ep$levels[j]
          )
        ) :=
          prob[, j]
      ]
    }
  } else {
    final_fit <- glmnet(
      x = X_dev,
      y =
        y_all[dev_idx],
      family = "binomial",
      alpha = 0,
      lambda =
        selected_lambda,
      standardize = TRUE
    )

    p <- as.numeric(
      predict(
        final_fit,
        newx = X_eval,
        type = "response",
        s = selected_lambda
      )
    )

    metrics <- binary_metric_row(
      y =
        y_all[eval_idx],
      p = p,
      endpoint = endpoint,
      eval_year = fold$eval
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
        "__ridge_predictions.csv"
      )
    )
  )

  qc <- data.table(
    endpoint = endpoint,
    predictor_policy = ep$policy,
    evaluation_year = fold$eval,
    n_train = length(train_idx),
    n_tune = length(tune_idx),
    n_dev_refit = length(dev_idx),
    n_eval = length(eval_idx),
    selected_lambda =
      selected_lambda,
    tuning_logloss =
      losses[best_j],
    lambda_path_index =
      best_j,
    lambda_path_n =
      length(lambda_path),
    lambda_path_max =
      lambda_path_max,
    lambda_path_min =
      lambda_path_min,
    selected_lambda_ratio_to_max =
      selected_lambda_ratio_to_max,
    selected_lambda_at_path_boundary =
      lambda_boundary,
    selected_at_near_zero_floor =
      selected_at_near_zero_floor,
    full_requested_lambda_path_completed =
      length(lambda_path) ==
        250L,
    n_predictors =
      length(predictors),
    n_encoded_features =
      ncol(X_dev)
  )

  rm(
    X_dev,
    X_eval,
    final_fit,
    pred_dt,
    enc_dev
  )
  gc()

  list(
    metrics = metrics,
    qc = qc
  )
}

# ------------------------------------------------------------------------------
# 6. Run comparator
# ------------------------------------------------------------------------------

run_plan <- CJ(
  endpoint =
    names(
      endpoint_specs
    ),
  evaluation_year =
    as.integer(
      names(folds)
    ),
  unique = TRUE
)

run_plan[
  ,
  endpoint_order :=
    match(
      endpoint,
      names(
        endpoint_specs
      )
    )
]

setorder(
  run_plan,
  endpoint_order,
  evaluation_year
)

run_plan[
  ,
  endpoint_order := NULL
]

fwrite(
  run_plan,
  file.path(
    OUT_DIR,
    "00_RUN_PLAN.csv"
  )
)

metric_rows <- list()
qc_rows <- list()

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
    "] RIDGE: ",
    endpoint,
    " | eval ",
    year,
    "\n============================================================\n",
    sep = ""
  )

  res <- fit_one(
    endpoint,
    folds[[year]]
  )

  metric_rows[[length(metric_rows) + 1L]] <- res$metrics

  qc_rows[[length(qc_rows) + 1L]] <- res$qc

  fwrite(
    rbindlist(
      metric_rows,
      fill = TRUE
    ),
    file.path(
      OUT_DIR,
      "01_RIDGE_TEMPORAL_METRICS_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(
      qc_rows,
      fill = TRUE
    ),
    file.path(
      OUT_DIR,
      "02_RIDGE_TRAINING_QC_CHECKPOINT.csv"
    )
  )
}

ridge_metrics <- rbindlist(
  metric_rows,
  fill = TRUE
)

ridge_qc <- rbindlist(
  qc_rows,
  fill = TRUE
)

fwrite(
  ridge_metrics,
  file.path(
    OUT_DIR,
    "01_RIDGE_TEMPORAL_METRICS.csv"
  )
)

fwrite(
  ridge_qc,
  file.path(
    OUT_DIR,
    "02_RIDGE_TRAINING_QC.csv"
  )
)

if (
  any(
    ridge_qc$full_requested_lambda_path_completed !=
      TRUE
  )
) {
  stop(
    "HARD FAIL: glmnet did not traverse all 250 requested lambdas despite path-truncation controls. Review 02_RIDGE_TRAINING_QC.csv.",
    call. = FALSE
  )
}

# A minimum-boundary selection is acceptable ONLY when the minimum is the
# prespecified near-zero floor (<=1.1e-8 of lambda_max). In that case the tuned
# comparator is effectively the least-regularized ridge model supported by the
# prespecified search range rather than an accidentally truncated path.
bad_boundary <- ridge_qc[
  selected_lambda_at_path_boundary ==
    TRUE &
    selected_at_near_zero_floor !=
      TRUE
]

if (nrow(bad_boundary) > 0L) {
  fwrite(
    bad_boundary,
    file.path(
      OUT_DIR,
      "02b_UNRESOLVED_LAMBDA_BOUNDARY.csv"
    )
  )

  stop(
    "HARD FAIL: at least one Ridge fit remains boundary-limited above the prespecified near-zero floor. Review 02b_UNRESOLVED_LAMBDA_BOUNDARY.csv.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 7. Fair XGBoost-vs-Ridge comparison
# ------------------------------------------------------------------------------

join_keys <- c(
  "endpoint",
  "evaluation_year",
  "metric_level",
  "class"
)

xgb_for_compare <- copy(
  xgb_metrics
)

keep_xgb <- intersect(
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
  names(xgb_for_compare)
)

xgb_for_compare <-
  xgb_for_compare[
    ,
    ..keep_xgb
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
  ridge_metrics
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
    "03_XGBOOST_VS_RIDGE_TEMPORAL_COMPARISON.csv"
  )
)

cohort_qc <- comparison[
  !is.na(xgb_n) &
    !is.na(ridge_n),
  .(
    rows_compared = .N,
    n_mismatch_rows =
      sum(
        xgb_n != ridge_n
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
    "04_RIDGE_XGBOOST_COHORT_CONCORDANCE_QC.csv"
  )
)

if (
  cohort_qc$n_mismatch_rows > 0L ||
    cohort_qc$prevalence_mismatch_rows > 0L
) {
  stop(
    "Ridge/XGBoost evaluation cohort concordance failed. Review file 04.",
    call. = FALSE
  )
}

headline_compare <- comparison[
  metric_level %in%
    c(
      "macro",
      "binary"
    )
]

fwrite(
  headline_compare,
  file.path(
    OUT_DIR,
    "05_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT REPAIRED-DX FAIR RIDGE COMPARATOR COMPLETE",
    "",
    "Ridge used the same corrected endpoint-specific predictor policy and rolling temporal folds as XGBoost.",
    "Numeric predictors: temporal-train median imputation + explicit missing indicator.",
    "Categorical predictors: explicit Unknown/Other + one-hot encoding.",
    "alpha=0; lambda selected on the fold-specific tuning year.",
    "v1.2 forces glmnet to traverse all 250 requested lambdas (fdev=0, devmax=1) with lambda.min.ratio=1e-8; a minimum-boundary selection is accepted only at that prespecified near-zero floor.",
    "",
    "Review first:",
    "  05_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv",
    "  03_XGBOOST_VS_RIDGE_TEMPORAL_COMPARISON.csv",
    "  04_RIDGE_XGBOOST_COHORT_CONCORDANCE_QC.csv",
    "  02_RIDGE_TRAINING_QC.csv"
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
  "REPAIRED-DX FAIR RIDGE COMPARATOR COMPLETE\n"
)
cat(
  "============================================================\n"
)
print(
  headline_compare[
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
