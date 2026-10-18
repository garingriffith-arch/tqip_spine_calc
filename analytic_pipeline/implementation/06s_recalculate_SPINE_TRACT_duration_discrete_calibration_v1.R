# ==============================================================================
# 06s_recalculate_SPINE_TRACT_duration_discrete_calibration_v1.R
# SPINE-TRACT Stage 3C-R — NO-REFIT duration calibration postprocessing
#
# PURPOSE
#   Recalculate interval coverage and quantile-calibration diagnostics for the
#   already-fitted repaired-DX duration models using a numerical equality
#   tolerance appropriate for integer-day outcomes.
#
# WHY
#   ICU Q10 predictions cluster at values such as 2.0000000595 days. Raw
#   floating-point comparison treats an observed 2-day stay as BELOW that
#   estimate, artificially depressing central-range coverage. A tolerance of
#   1e-6 day (~0.086 seconds) treats these values as numerically equivalent.
#
#   For discrete outcomes, exact empirical P(Y <= Q_alpha) need not equal alpha
#   because probability mass can occur at the predicted integer quantile.
#   Therefore this script additionally reports:
#
#     P(Y < Q_alpha - tol)
#     P(Y <= Q_alpha + tol)
#
#   A nominal alpha can legitimately fall between those two quantities.
#
# THIS SCRIPT DOES NOT REFIT OR RETUNE ANY MODEL.
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

DURATION_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_duration_quantiles_REPAIRED_DX",
  "rolling_origin_gpu"
)

PRED_DIR <- file.path(
  DURATION_DIR,
  "predictions"
)

METRICS_FILE <- file.path(
  DURATION_DIR,
  "05_FINAL_DURATION_TEMPORAL_METRICS.csv"
)

if (!file.exists(METRICS_FILE) || !dir.exists(PRED_DIR)) {
  stop(
    "Missing repaired-DX duration metrics or prediction directory.",
    call. = FALSE
  )
}

TOL <- 1e-6

metrics <- fread(
  METRICS_FILE
)

pred_files <- list.files(
  PRED_DIR,
  pattern = "__eval[0-9]{4}__Q10_Q50_Q90_predictions\\.csv$",
  full.names = TRUE
)

if (length(pred_files) != 9L) {
  stop(
    "Expected 9 repaired-DX duration prediction files; found ",
    length(pred_files),
    ".",
    call. = FALSE
  )
}

parse_file <- function(path) {
  b <- basename(path)

  outcome <- sub(
    "__eval[0-9]{4}__Q10_Q50_Q90_predictions\\.csv$",
    "",
    b
  )

  year <- as.integer(
    sub(
      ".*__eval([0-9]{4})__Q10_Q50_Q90_predictions\\.csv$",
      "\\1",
      b
    )
  )

  list(
    outcome = outcome,
    evaluation_year = year
  )
}

diag_rows <- list()

for (f in pred_files) {
  id <- parse_file(f)
  d <- fread(f)

  needed <- c(
    "observed_days",
    "q10_days",
    "q50_median_days",
    "q90_days"
  )

  if (!all(needed %in% names(d))) {
    stop(
      "Prediction file missing expected columns: ",
      f,
      call. = FALSE
    )
  }

  y <- as.numeric(d$observed_days)
  q10 <- as.numeric(d$q10_days)
  q50 <- as.numeric(d$q50_median_days)
  q90 <- as.numeric(d$q90_days)

  ok <- is.finite(y) &
    is.finite(q10) &
    is.finite(q50) &
    is.finite(q90)

  y <- y[ok]
  q10 <- q10[ok]
  q50 <- q50[ok]
  q90 <- q90[ok]

  row <- data.table(
    outcome = id$outcome,
    evaluation_year = id$evaluation_year,
    N = length(y),

    tolerance_days = TOL,

    p_y_lt_q10 =
      mean(y < q10 - TOL),
    p_y_le_q10 =
      mean(y <= q10 + TOL),

    p_y_lt_q50 =
      mean(y < q50 - TOL),
    p_y_le_q50 =
      mean(y <= q50 + TOL),

    p_y_lt_q90 =
      mean(y < q90 - TOL),
    p_y_le_q90 =
      mean(y <= q90 + TOL),

    q10_nominal_bracketed =
      mean(y < q10 - TOL) <= 0.10 &&
        mean(y <= q10 + TOL) >= 0.10,

    q50_nominal_bracketed =
      mean(y < q50 - TOL) <= 0.50 &&
        mean(y <= q50 + TOL) >= 0.50,

    q90_nominal_bracketed =
      mean(y < q90 - TOL) <= 0.90 &&
        mean(y <= q90 + TOL) >= 0.90,

    q10_q90_coverage_tolerance =
      mean(
        y + TOL >= q10 &
          y - TOL <= q90
      ),

    lower_tail_miss_tolerance =
      mean(
        y + TOL < q10
      ),

    upper_tail_miss_tolerance =
      mean(
        y - TOL > q90
      ),

    q10_exact_integer_proximity_n =
      sum(
        abs(
          q10 -
            round(q10)
        ) <= TOL
      ),

    q50_exact_integer_proximity_n =
      sum(
        abs(
          q50 -
            round(q50)
        ) <= TOL
      ),

    q90_exact_integer_proximity_n =
      sum(
        abs(
          q90 -
            round(q90)
        ) <= TOL
      )
  )

  diag_rows[[length(diag_rows) + 1L]] <- row
}

diag <- rbindlist(
  diag_rows,
  fill = TRUE
)

setorder(
  diag,
  outcome,
  evaluation_year
)

fwrite(
  diag,
  file.path(
    DURATION_DIR,
    "08_DURATION_DISCRETE_QUANTILE_CALIBRATION.csv"
  )
)

corrected <- merge(
  metrics,
  diag[
    ,
    .(
      outcome,
      evaluation_year,
      p_y_lt_q10,
      p_y_le_q10,
      p_y_lt_q50,
      p_y_le_q50,
      p_y_lt_q90,
      p_y_le_q90,
      q10_nominal_bracketed,
      q50_nominal_bracketed,
      q90_nominal_bracketed,
      q10_q90_coverage_tolerance,
      lower_tail_miss_tolerance,
      upper_tail_miss_tolerance
    )
  ],
  by = c(
    "outcome",
    "evaluation_year"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(corrected) != nrow(metrics)) {
  stop(
    "Corrected metric merge changed row count.",
    call. = FALSE
  )
}

corrected[
  ,
  q10_q90_coverage_raw_float :=
    q10_q90_coverage
]

corrected[
  ,
  q10_q90_coverage :=
    q10_q90_coverage_tolerance
]

fwrite(
  corrected,
  file.path(
    DURATION_DIR,
    "05b_FINAL_DURATION_TEMPORAL_METRICS_DISCRETE_TOLERANCE.csv"
  )
)

summary <- corrected[
  ,
  .(
    min_MAE_days =
      min(MAE_days, na.rm = TRUE),
    median_MAE_days =
      median(MAE_days, na.rm = TRUE),
    max_MAE_days =
      max(MAE_days, na.rm = TRUE),

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

    median_p_y_lt_q10 =
      median(
        p_y_lt_q10,
        na.rm = TRUE
      ),

    median_p_y_le_q10 =
      median(
        p_y_le_q10,
        na.rm = TRUE
      ),

    median_p_y_lt_q50 =
      median(
        p_y_lt_q50,
        na.rm = TRUE
      ),

    median_p_y_le_q50 =
      median(
        p_y_le_q50,
        na.rm = TRUE
      ),

    median_p_y_lt_q90 =
      median(
        p_y_lt_q90,
        na.rm = TRUE
      ),

    median_p_y_le_q90 =
      median(
        p_y_le_q90,
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
  summary,
  file.path(
    DURATION_DIR,
    "07b_DURATION_TEMPORAL_SUMMARY_DISCRETE_TOLERANCE.csv"
  )
)

cat(
  "\n============================================================\n"
)
cat(
  "NO-REFIT DISCRETE-DURATION CALIBRATION RECALCULATION COMPLETE\n"
)
cat(
  "============================================================\n"
)
cat(
  "Numerical tolerance: ",
  format(TOL, scientific = TRUE),
  " day\n\n",
  sep = ""
)
print(
  corrected[
    ,
    .(
      outcome,
      evaluation_year,
      MAE_days,
      q10_q90_coverage_raw_float,
      q10_q90_coverage,
      p_y_lt_q10,
      p_y_le_q10,
      p_y_lt_q90,
      p_y_le_q90
    )
  ][
    order(
      outcome,
      evaluation_year
    )
  ]
)
