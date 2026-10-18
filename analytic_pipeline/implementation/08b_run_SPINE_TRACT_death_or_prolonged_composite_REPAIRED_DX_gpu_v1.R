# ==============================================================================
# 08b_run_SPINE_TRACT_death_or_prolonged_composite_REPAIRED_DX_gpu_v1.R
# SPINE-TRACT targeted early-death composite sensitivity
#
# PURPOSE
#   Fit the one remaining prespecified heavy robustness analysis after the
#   repaired-DX primary model freeze:
#
#     1) Death/hospice OR hospital LOS >=20 days
#     2) Death/hospice OR ICU LOS >=8 days
#     3) Death/hospice OR ventilation >=8 days
#
# These are SENSITIVITY outcomes only. They do not replace the locked primary
# HLOS/ICU/ventilation trajectory models.
#
# PREDICTOR POLICY
#   death_or_hlos20 -> corrected final_context (57 predictors)
#   death_or_icu8   -> corrected final_core    (54 predictors)
#   death_or_vent8  -> corrected final_core    (54 predictors)
#
# TEMPORAL DESIGN
#   Fold 1: train 2020      | tune 2021 | refit 2020-21 | eval 2022
#   Fold 2: train 2020-21   | tune 2022 | refit 2020-22 | eval 2023
#   Fold 3: train 2020-22   | tune 2023 | refit 2020-23 | eval 2024
#
# HYPERPARAMETER SEARCH
#   Reuses the exact 24-configuration prespecified Stage 3A structural grid.
#   eta = 0.05; early stopping selects fold-specific boosting rounds.
#
# IMPORTANT
#   - No primary model is changed.
#   - 2024 remains development-informed forward temporal evaluation.
#   - Death/hospice is treated as a positive composite event irrespective of
#     observed resource duration.
#   - Among non-deaths, the resource component is positive only when the locked
#     categorical endpoint is observed in the prolonged class.
# ==============================================================================

rm(list = ls())
gc()

required <- c(
  "data.table",
  "Matrix",
  "xgboost"
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
  "stage4_death_or_prolonged_composite_REPAIRED_DX",
  "rolling_origin_gpu"
)

GRID_DIR <- file.path(
  OUT_DIR,
  "grid_search"
)

MODEL_DIR <- file.path(
  OUT_DIR,
  "models"
)

PRED_DIR <- file.path(
  OUT_DIR,
  "predictions"
)

for (d in c(
  OUT_DIR,
  GRID_DIR,
  MODEL_DIR,
  PRED_DIR
)) {
  dir.create(
    d,
    recursive = TRUE,
    showWarnings = FALSE
  )
}

needed <- c(
  DATA_FILE,
  SPEC_FILE,
  GRID_FILE
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

MODEL_SEED <- 20261002L
ETA <- 0.05
MAX_ROUNDS <- 3000L
EARLY_STOP <- 75L

logical_cores <- parallel::detectCores(
  logical = TRUE
)

if (
  is.na(logical_cores) ||
    logical_cores < 1L
) {
  logical_cores <- 1L
}

CPU_THREADS <- max(
  1L,
  min(
    logical_cores,
    as.integer(
      round(
        logical_cores *
          0.80
      )
    )
  )
)

data.table::setDTthreads(
  CPU_THREADS
)

# ------------------------------------------------------------------------------
# 1. Load repaired primary dataset and predictor lock
# ------------------------------------------------------------------------------

dt <- as.data.table(
  readRDS(
    DATA_FILE
  )
)

spec <- readRDS(
  SPEC_FILE
)

grid <- fread(
  GRID_FILE
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

required_grid_cols <- c(
  "config_id",
  "max_depth",
  "min_child_weight",
  "subsample",
  "colsample_bytree",
  "lambda"
)

if (
  nrow(grid) != 24L ||
    !all(
      required_grid_cols %in%
        names(grid)
    )
) {
  stop(
    "Prespecified Stage 3A search grid is malformed.",
    call. = FALSE
  )
}

grid <- grid[
  ,
  ..required_grid_cols
]

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

final_context <- unique(c(
  final_core,
  pg$social_context
))

if (
  length(final_core) != 54L ||
    length(final_context) != 57L ||
    !"dx_tl_soft_tissue_any" %in%
      final_core
) {
  stop(
    "Corrected 54/57 predictor lock could not be reconstructed.",
    call. = FALSE
  )
}

# Repaired-DX preflight.
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

preflight <- repair_expected[
  ,
  .(
    predictor,
    expected_positive_n,
    observed_positive_n =
      vapply(
        predictor,
        function(v) {
          sum(
            dt[[v]] == 1L,
            na.rm = TRUE
          )
        },
        integer(1)
      ),
    missing_n =
      vapply(
        predictor,
        function(v) {
          sum(
            is.na(
              dt[[v]]
            )
          )
        },
        integer(1)
      )
  )
]

preflight[
  ,
  concordant :=
    observed_positive_n ==
      expected_positive_n &
      missing_n ==
        0L
]

fwrite(
  preflight,
  file.path(
    OUT_DIR,
    "00_REPAIRED_DX_PREFLIGHT.csv"
  )
)

if (!all(preflight$concordant)) {
  stop(
    "Repaired diagnosis preflight failed.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Define composite sensitivity outcomes
# ------------------------------------------------------------------------------

death_hospice <- as.integer(
  as.character(
    dt$discharge4
  ) ==
    "Death/hospice"
)

hlos_obs <- as.character(
  dt$hlos3
)

icu_obs <- as.character(
  dt$icu3
)

vent_obs <- as.character(
  dt$vent3
)

dt[
  ,
  death_or_hlos20 :=
    fcase(
      death_hospice ==
        1L,
      1L,

      death_hospice ==
        0L &
        !is.na(hlos_obs) &
        hlos_obs ==
          ">=20 days",
      1L,

      death_hospice ==
        0L &
        !is.na(hlos_obs),
      0L,

      default =
        NA_integer_
    )
]

dt[
  ,
  death_or_icu8 :=
    fcase(
      death_hospice ==
        1L,
      1L,

      death_hospice ==
        0L &
        !is.na(icu_obs) &
        icu_obs ==
          "ICU >=8 days",
      1L,

      death_hospice ==
        0L &
        !is.na(icu_obs),
      0L,

      default =
        NA_integer_
    )
]

dt[
  ,
  death_or_vent8 :=
    fcase(
      death_hospice ==
        1L,
      1L,

      death_hospice ==
        0L &
        !is.na(vent_obs) &
        vent_obs ==
          "Ventilation >=8 days",
      1L,

      death_hospice ==
        0L &
        !is.na(vent_obs),
      0L,

      default =
        NA_integer_
    )
]

endpoint_specs <- list(
  death_or_hlos20 = list(
    label =
      "Death/hospice or hospital LOS >=20 days",
    predictors =
      final_context,
    policy =
      "final_context"
  ),

  death_or_icu8 = list(
    label =
      "Death/hospice or ICU LOS >=8 days",
    predictors =
      final_core,
    policy =
      "final_core"
  ),

  death_or_vent8 = list(
    label =
      "Death/hospice or ventilation >=8 days",
    predictors =
      final_core,
    policy =
      "final_core"
  )
)

endpoint_counts <- rbindlist(
  lapply(
    names(endpoint_specs),
    function(ep) {
      dt[
        ,
        .(
          N =
            sum(
              !is.na(
                get(ep)
              )
            ),
          events =
            sum(
              get(ep) ==
                1L,
              na.rm = TRUE
            ),
          prevalence =
            mean(
              get(ep) ==
                1L,
              na.rm = TRUE
            )
        ),
        by =
          admission_year
      ][
        ,
        endpoint :=
          ep
      ]
    }
  ),
  fill = TRUE
)

setcolorder(
  endpoint_counts,
  c(
    "endpoint",
    "admission_year",
    "N",
    "events",
    "prevalence"
  )
)

fwrite(
  endpoint_counts,
  file.path(
    OUT_DIR,
    "01_COMPOSITE_ENDPOINT_COUNTS_BY_YEAR.csv"
  )
)

# ------------------------------------------------------------------------------
# 3. GPU / encoding / metric helpers
# ------------------------------------------------------------------------------

config_to_text <- function(model) {
  cfg <- xgb.config(model)

  if (is.character(cfg)) {
    paste(
      cfg,
      collapse = ""
    )
  } else if (
    requireNamespace(
      "jsonlite",
      quietly = TRUE
    )
  ) {
    jsonlite::toJSON(
      cfg,
      auto_unbox = TRUE,
      null = "null"
    )
  } else {
    paste(
      capture.output(
        str(cfg)
      ),
      collapse = " "
    )
  }
}

gpu_params_for_version <- function() {
  ver <- packageVersion(
    "xgboost"
  )

  if (
    ver <
      package_version(
        "2.0.0"
      )
  ) {
    list(
      mode = paste0(
        "R CUDA legacy xgboost ",
        as.character(ver),
        " | gpu_hist | gpu_id=0"
      ),
      extra = list(
        tree_method =
          "gpu_hist",
        predictor =
          "gpu_predictor",
        gpu_id =
          0L
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
        tree_method =
          "hist",
        device =
          "cuda:0"
      )
    )
  }
}

assert_cuda_booster <- function(
    model,
    context
) {
  txt <- tolower(
    config_to_text(
      model
    )
  )

  ver <- packageVersion(
    "xgboost"
  )

  ok <- if (
    ver <
      package_version(
        "2.0.0"
      )
  ) {
    grepl(
      "grow_gpu_hist|gpu_hist",
      txt,
      perl = TRUE
    )
  } else {
    grepl(
      '"device"\\s*:\\s*"cuda(:0)?"',
      txt,
      perl = TRUE
    ) ||
      grepl(
        "grow_gpu_hist",
        txt,
        perl = TRUE
      )
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

is_numeric_predictor <- function(x) {
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

fit_encoder <- function(
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

    if (
      is_numeric_predictor(x)
    ) {
      enc[[v]] <- list(
        type = "numeric"
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

encode_with <- function(
    d,
    predictors,
    enc
) {
  blocks <- vector(
    "list",
    length(predictors)
  )

  for (i in seq_along(predictors)) {
    v <- predictors[i]
    info <- enc[[v]]
    x <- d[[v]]

    if (
      info$type ==
        "numeric"
    ) {
      m <- Matrix(
        as.numeric(x),
        ncol = 1L,
        sparse = TRUE
      )

      colnames(m) <- v
      blocks[[i]] <- m
    } else {
      z <- clean_cat(x)

      z[
        !z %in%
          info$levels
      ] <- "__OTHER__"

      f <- factor(
        z,
        levels =
          info$levels
      )

      m <- sparseMatrix(
        i =
          seq_along(f),
        j =
          as.integer(f),
        x =
          1,
        dims =
          c(
            length(f),
            length(
              info$levels
            )
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

      blocks[[i]] <- m
    }
  }

  do.call(
    cbind,
    blocks
  )
}

clamp_prob <- function(
    p,
    eps = 1e-7
) {
  pmin(
    pmax(
      as.numeric(p),
      eps
    ),
    1 -
      eps
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
      y ==
        1L
    )
  )

  n0 <- as.double(
    sum(
      y ==
        0L
    )
  )

  if (
    n1 == 0 ||
      n0 == 0
  ) {
    return(
      NA_real_
    )
  }

  r <- rank(
    p,
    ties.method =
      "average"
  )

  (
    sum(
      r[
        y ==
          1L
      ]
    ) -
      n1 *
        (
          n1 +
            1
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
    y ==
      1L
  )

  if (
    npos ==
      0L
  ) {
    return(
      NA_real_
    )
  }

  ord <- order(
    p,
    decreasing =
      TRUE
  )

  yy <- y[ord]

  precision <- cumsum(
    yy ==
      1L
  ) /
    seq_along(yy)

  mean(
    precision[
      yy ==
        1L
    ]
  )
}

calibration <- function(y, p) {
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
        intercept =
          NA_real_,
        slope =
          NA_real_
      )
    )
  }

  fi <- tryCatch(
    suppressWarnings(
      glm(
        y ~ offset(lp),
        family =
          binomial()
      )
    ),
    error =
      function(e) NULL
  )

  fs <- tryCatch(
    suppressWarnings(
      glm(
        y ~ lp,
        family =
          binomial()
      )
    ),
    error =
      function(e) NULL
  )

  c(
    intercept =
      if (
        is.null(fi)
      ) {
        NA_real_
      } else {
        unname(
          coef(fi)[1L]
        )
      },
    slope =
      if (
        is.null(fs)
      ) {
        NA_real_
      } else {
        unname(
          coef(fs)[2L]
        )
      }
  )
}

binary_metrics <- function(
    y,
    p
) {
  pp <- clamp_prob(p)
  prev <- mean(
    y ==
      1L
  )

  cal <- calibration(
    y,
    p
  )

  pr <- fast_auprc(
    y,
    p
  )

  data.table(
    N = length(y),
    events =
      sum(
        y ==
          1L
      ),
    prevalence =
      prev,
    AUROC =
      fast_auc(
        y,
        p
      ),
    AUPRC =
      pr,
    PR_lift =
      if (
        is.finite(prev) &&
          prev >
            0
      ) {
        pr /
          prev
      } else {
        NA_real_
      },
    Brier =
      mean(
        (
          y -
            p
        )^2
      ),
    log_loss =
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
      )
  )
}

extract_best_iteration <- function(
    model,
    fallback
) {
  b <- tryCatch(
    model$best_iteration,
    error =
      function(e) NULL
  )

  if (
    !is.null(b) &&
      length(b) ==
        1L &&
      is.finite(b)
  ) {
    return(
      as.integer(b)
    )
  }

  a <- tryCatch(
    xgb.attr(
      model,
      "best_iteration"
    ),
    error =
      function(e) NULL
  )

  if (
    !is.null(a) &&
      length(a) ==
        1L &&
      !is.na(a)
  ) {
    return(
      as.integer(a) +
        1L
    )
  }

  as.integer(
    fallback
  )
}

extract_best_score <- function(
    model
) {
  a <- tryCatch(
    xgb.attr(
      model,
      "best_score"
    ),
    error =
      function(e) NULL
  )

  if (
    !is.null(a) &&
      length(a) >=
        1L &&
      !is.na(
        a[1L]
      )
  ) {
    z <- suppressWarnings(
      as.numeric(
        a[1L]
      )
    )

    if (
      is.finite(z)
    ) {
      return(z)
    }
  }

  log <- tryCatch(
    model$evaluation_log,
    error =
      function(e) NULL
  )

  if (
    !is.null(log) &&
      nrow(log) >
        0L
  ) {
    tune_cols <- grep(
      "^tune_",
      names(log),
      value =
        TRUE
    )

    if (
      length(tune_cols) >=
        1L
    ) {
      return(
        min(
          as.numeric(
            log[[
              tune_cols[1L]
            ]]
          ),
          na.rm =
            TRUE
        )
      )
    }
  }

  NA_real_
}

make_params <- function(g) {
  c(
    list(
      objective =
        "binary:logistic",
      eval_metric =
        "logloss",
      eta =
        ETA,
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
      nthread =
        CPU_THREADS,
      verbosity =
        0
    ),
    compute$extra
  )
}

# GPU hard probe.
set.seed(1L)

probe_x <- matrix(
  rnorm(20000),
  nrow =
    1000L,
  ncol =
    20L
)

probe_y <- rbinom(
  1000L,
  1L,
  0.3
)

probe_d <- xgb.DMatrix(
  probe_x,
  label =
    probe_y
)

probe_fit <- xgb.train(
  params = c(
    list(
      objective =
        "binary:logistic",
      eval_metric =
        "logloss",
      max_depth =
        2L,
      eta =
        0.2,
      nthread =
        CPU_THREADS,
      verbosity =
        0
    ),
    compute$extra
  ),
  data =
    probe_d,
  nrounds =
    3L,
  verbose =
    0
)

assert_cuda_booster(
  probe_fit,
  "startup GPU probe"
)

fwrite(
  data.table(
    xgboost_version =
      as.character(
        packageVersion(
          "xgboost"
        )
      ),
    compute_mode =
      compute$mode,
    logical_cores =
      logical_cores,
    allocated_cpu_threads =
      CPU_THREADS
  ),
  file.path(
    OUT_DIR,
    "02_COMPUTE_ENVIRONMENT.csv"
  )
)

rm(
  probe_x,
  probe_y,
  probe_d,
  probe_fit
)

gc()

# ------------------------------------------------------------------------------
# 4. Rolling-origin folds
# ------------------------------------------------------------------------------

folds <- list(
  `2022` = list(
    train =
      2020L,
    tune =
      2021L,
    dev =
      2020:2021,
    eval =
      2022L
  ),

  `2023` = list(
    train =
      2020:2021,
    tune =
      2022L,
    dev =
      2020:2022,
    eval =
      2023L
  ),

  `2024` = list(
    train =
      2020:2022,
    tune =
      2023L,
    dev =
      2020:2023,
    eval =
      2024L
  )
)

# ------------------------------------------------------------------------------
# 5. Fit one composite endpoint x fold
# ------------------------------------------------------------------------------

fit_one <- function(
    endpoint,
    fold
) {
  ep <- endpoint_specs[[endpoint]]
  predictors <- ep$predictors

  y_all <- as.integer(
    dt[[endpoint]]
  )

  train_idx <- which(
    dt$admission_year %in%
      fold$train &
      !is.na(y_all)
  )

  tune_idx <- which(
    dt$admission_year %in%
      fold$tune &
      !is.na(y_all)
  )

  dev_idx <- which(
    dt$admission_year %in%
      fold$dev &
      !is.na(y_all)
  )

  eval_idx <- which(
    dt$admission_year ==
      fold$eval &
      !is.na(y_all)
  )

  if (
    min(
      length(train_idx),
      length(tune_idx),
      length(dev_idx),
      length(eval_idx)
    ) <=
      0L
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

  if (
    !identical(
      colnames(X_train),
      colnames(X_tune)
    )
  ) {
    stop(
      "Tuning encoder mismatch for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  dtrain <- xgb.DMatrix(
    X_train,
    label =
      y_all[train_idx],
    missing =
      NA_real_
  )

  dtune <- xgb.DMatrix(
    X_tune,
    label =
      y_all[tune_idx],
    missing =
      NA_real_
  )

  grid_rows <- vector(
    "list",
    nrow(grid)
  )

  for (
    i in seq_len(
      nrow(grid)
    )
  ) {
    g <- grid[i]

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
      nrow(grid),
      ")\n",
      sep = ""
    )

    set.seed(
      MODEL_SEED +
        fold$eval +
        i
    )

    params <- make_params(g)

    if (
      packageVersion(
        "xgboost"
      ) <
        package_version(
          "2.0.0"
        )
    ) {
      fit_es <- xgb.train(
        params =
          params,
        data =
          dtrain,
        nrounds =
          MAX_ROUNDS,
        watchlist = list(
          train =
            dtrain,
          tune =
            dtune
        ),
        early_stopping_rounds =
          EARLY_STOP,
        maximize =
          FALSE,
        verbose =
          0
      )
    } else {
      fit_es <- xgb.train(
        params =
          params,
        data =
          dtrain,
        nrounds =
          MAX_ROUNDS,
        evals = list(
          train =
            dtrain,
          tune =
            dtune
        ),
        early_stopping_rounds =
          EARLY_STOP,
        maximize =
          FALSE,
        verbose =
          0
      )
    }

    if (
      i ==
        1L
    ) {
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
        endpoint =
          endpoint,
        endpoint_label =
          ep$label,
        predictor_policy =
          ep$policy,
        evaluation_year =
          fold$eval,
        best_iteration =
          extract_best_iteration(
            fit_es,
            MAX_ROUNDS
          ),
        tune_logloss =
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
    fill =
      TRUE
  )

  if (
    all(
      !is.finite(
        grid_table$tune_logloss
      )
    )
  ) {
    stop(
      "All tuning scores non-finite for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  setorder(
    grid_table,
    tune_logloss,
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

  # Refit on all pre-evaluation development years.
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

  if (
    !identical(
      colnames(X_dev),
      colnames(X_eval)
    )
  ) {
    stop(
      "Refit encoder mismatch for ",
      endpoint,
      " eval ",
      fold$eval,
      ".",
      call. = FALSE
    )
  }

  ddev <- xgb.DMatrix(
    X_dev,
    label =
      y_all[dev_idx],
    missing =
      NA_real_
  )

  deval <- xgb.DMatrix(
    X_eval,
    missing =
      NA_real_
  )

  set.seed(
    MODEL_SEED +
      10000L +
      fold$eval
  )

  fit <- xgb.train(
    params =
      make_params(
        selected
      ),
    data =
      ddev,
    nrounds =
      as.integer(
        selected$best_iteration
      ),
    verbose =
      0
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

  p <- as.numeric(
    predict(
      fit,
      deval
    )
  )

  y_eval <- y_all[eval_idx]

  met <- binary_metrics(
    y_eval,
    p
  )

  met[, `:=`(
    endpoint =
      endpoint,
    endpoint_label =
      ep$label,
    predictor_policy =
      ep$policy,
    evaluation_year =
      fold$eval,
    selected_config =
      selected$config_id,
    selected_rounds =
      as.integer(
        selected$best_iteration
      ),
    tune_logloss =
      as.numeric(
        selected$tune_logloss
      )
  )]

  pred_dt <- data.table(
    admission_year =
      dt$admission_year[eval_idx],
    inc_key =
      dt$inc_key[eval_idx],
    y =
      y_eval,
    pred =
      p
  )

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
      model =
        fit,
      endpoint =
        endpoint,
      endpoint_label =
        ep$label,
      evaluation_year =
        fold$eval,
      predictor_policy =
        ep$policy,
      predictors =
        predictors,
      encoder =
        enc_dev,
      feature_names =
        colnames(X_dev),
      selected_config =
        selected,
      sensitivity_only =
        TRUE
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
    compress =
      FALSE
  )

  qc <- data.table(
    endpoint =
      endpoint,
    endpoint_label =
      ep$label,
    predictor_policy =
      ep$policy,
    evaluation_year =
      fold$eval,
    train_years =
      paste(
        fold$train,
        collapse = ","
      ),
    tune_years =
      paste(
        fold$tune,
        collapse = ","
      ),
    dev_refit_years =
      paste(
        fold$dev,
        collapse = ","
      ),
    n_train =
      length(train_idx),
    n_tune =
      length(tune_idx),
    n_dev_refit =
      length(dev_idx),
    n_eval =
      length(eval_idx),
    n_predictors =
      length(predictors),
    n_encoded_features =
      ncol(X_dev),
    selected_config =
      selected$config_id,
    selected_rounds =
      as.integer(
        selected$best_iteration
      ),
    tune_logloss =
      as.numeric(
        selected$tune_logloss
      ),
    hit_round_ceiling =
      as.integer(
        selected$best_iteration
      ) >=
        (
          MAX_ROUNDS -
            5L
        )
  )

  rm(
    X_dev,
    X_eval,
    ddev,
    deval,
    fit,
    p,
    pred_dt,
    enc_dev
  )

  gc()

  list(
    metrics =
      met,
    selected =
      selected,
    qc =
      qc
  )
}

# ------------------------------------------------------------------------------
# 6. Run 3 composites x 3 folds = 9 final sensitivity evaluations
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
  unique =
    TRUE
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
  endpoint_order :=
    NULL
]

fwrite(
  run_plan,
  file.path(
    OUT_DIR,
    "03_RUN_PLAN.csv"
  )
)

metric_rows <- list()
selected_rows <- list()
qc_rows <- list()
runtime_rows <- list()

for (
  i in seq_len(
    nrow(
      run_plan
    )
  )
) {
  endpoint <-
    run_plan$endpoint[i]

  year <- as.character(
    run_plan$evaluation_year[i]
  )

  cat(
    "\n============================================================\n",
    "[",
    i,
    "/",
    nrow(run_plan),
    "] COMPOSITE SENSITIVITY: ",
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

  metric_rows[[length(metric_rows) + 1L]] <- res$metrics

  s <- copy(
    res$selected
  )

  s[
    ,
    `:=`(
      endpoint =
        endpoint,
      evaluation_year =
        as.integer(year)
    )
  ]

  selected_rows[[length(selected_rows) + 1L]] <- s

  qc_rows[[length(qc_rows) + 1L]] <- res$qc

  runtime_rows[[length(runtime_rows) + 1L]] <- data.table(
    endpoint =
      endpoint,
    evaluation_year =
      as.integer(year),
    runtime_minutes =
      as.numeric(
        difftime(
          Sys.time(),
          t0,
          units =
            "mins"
        )
      )
  )

  fwrite(
    rbindlist(
      metric_rows,
      fill =
        TRUE
    ),
    file.path(
      OUT_DIR,
      "04_COMPOSITE_TEMPORAL_METRICS_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(
      selected_rows,
      fill =
        TRUE
    ),
    file.path(
      OUT_DIR,
      "05_SELECTED_CONFIGS_CHECKPOINT.csv"
    )
  )

  fwrite(
    rbindlist(
      qc_rows,
      fill =
        TRUE
    ),
    file.path(
      OUT_DIR,
      "06_TRAINING_QC_CHECKPOINT.csv"
    )
  )
}

metrics <- rbindlist(
  metric_rows,
  fill =
    TRUE
)

selected <- rbindlist(
  selected_rows,
  fill =
    TRUE
)

qc <- rbindlist(
  qc_rows,
  fill =
    TRUE
)

runtime <- rbindlist(
  runtime_rows,
  fill =
    TRUE
)

if (
  nrow(metrics) !=
    9L ||
    nrow(qc) !=
      9L
) {
  stop(
    "Expected 9 completed composite sensitivity evaluations.",
    call. = FALSE
  )
}

if (
  any(
    qc$hit_round_ceiling,
    na.rm =
      TRUE
  )
) {
  warning(
    "At least one composite model hit the 3000-round ceiling. Review 06_TRAINING_QC.csv."
  )
}

fwrite(
  metrics,
  file.path(
    OUT_DIR,
    "04_COMPOSITE_TEMPORAL_METRICS.csv"
  )
)

fwrite(
  selected,
  file.path(
    OUT_DIR,
    "05_SELECTED_CONFIGS.csv"
  )
)

fwrite(
  qc,
  file.path(
    OUT_DIR,
    "06_TRAINING_QC.csv"
  )
)

fwrite(
  runtime,
  file.path(
    OUT_DIR,
    "07_RUNTIME.csv"
  )
)

# 2024 manuscript-facing compact table.
metrics_2024 <- metrics[
  evaluation_year ==
    2024L
]

fwrite(
  metrics_2024,
  file.path(
    OUT_DIR,
    "08_2024_COMPOSITE_SENSITIVITY.csv"
  )
)

temporal_summary <- metrics[
  ,
  .(
    prevalence_min =
      min(
        prevalence,
        na.rm = TRUE
      ),
    prevalence_max =
      max(
        prevalence,
        na.rm = TRUE
      ),
    AUROC_min =
      min(
        AUROC,
        na.rm = TRUE
      ),
    AUROC_median =
      median(
        AUROC,
        na.rm = TRUE
      ),
    AUROC_max =
      max(
        AUROC,
        na.rm = TRUE
      ),
    AUPRC_min =
      min(
        AUPRC,
        na.rm = TRUE
      ),
    AUPRC_median =
      median(
        AUPRC,
        na.rm = TRUE
      ),
    AUPRC_max =
      max(
        AUPRC,
        na.rm = TRUE
      ),
    max_abs_calibration_intercept =
      max(
        abs(
          calibration_intercept
        ),
        na.rm = TRUE
      ),
    max_abs_calibration_slope_error =
      max(
        abs(
          calibration_slope -
            1
        ),
        na.rm = TRUE
      )
  ),
  by = .(
    endpoint,
    endpoint_label,
    predictor_policy
  )
]

fwrite(
  temporal_summary,
  file.path(
    OUT_DIR,
    "09_COMPOSITE_TEMPORAL_SUMMARY.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT DEATH-OR-PROLONGED COMPOSITE SENSITIVITY COMPLETE",
    "",
    "Sensitivity outcomes only; primary model definitions are unchanged.",
    "",
    "Composite outcomes:",
    "  Death/hospice OR HLOS >=20 days",
    "  Death/hospice OR ICU LOS >=8 days",
    "  Death/hospice OR ventilation >=8 days",
    "",
    "Predictor policy:",
    "  HLOS composite -> final_context (57)",
    "  ICU/vent composites -> final_core (54)",
    "",
    "Temporal design:",
    "  rolling evaluation 2022, 2023, 2024",
    "  exact 24-config Stage 3A structural grid",
    "",
    "Review first:",
    "  08_2024_COMPOSITE_SENSITIVITY.csv",
    "  09_COMPOSITE_TEMPORAL_SUMMARY.csv",
    "  06_TRAINING_QC.csv",
    "  05_SELECTED_CONFIGS.csv"
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
  "DEATH-OR-PROLONGED COMPOSITE SENSITIVITY COMPLETE\n"
)
cat(
  "============================================================\n\n"
)

print(
  metrics[
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
