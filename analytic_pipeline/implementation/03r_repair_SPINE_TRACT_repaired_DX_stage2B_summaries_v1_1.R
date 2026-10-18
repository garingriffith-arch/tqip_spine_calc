# ==============================================================================
# 03r_repair_SPINE_TRACT_repaired_DX_stage2B_summaries_v1_1.R
# SPINE-TRACT Stage 2B-R — POSTPROCESSING REPAIR ONLY
#
# PURPOSE
#   Repair the summary-generation bug that occurs AFTER all 135 repaired-DX GPU
#   fits have completed.
#
# ROOT CAUSE
#   multiclass_metrics() created the "macro" row with:
#
#     class_rows[, .(
#       endpoint = endpoint,
#       variant = variant,
#       evaluation_year = fold_year,
#       ...
#     )]
#
#   Inside data.table j, endpoint and variant referred to full columns of
#   class_rows rather than scalar values. Therefore the supposed macro summary
#   was recycled to K rows (4 for discharge, 3 for HLOS/ICU/vent) instead of one.
#
#   The model fits, class-specific metrics, prediction caches, and model objects
#   are valid. Only the redundant macro rows and downstream delta-summary merge
#   are affected.
#
# THIS SCRIPT:
#   - DOES NOT FIT ANY MODEL;
#   - reads the already-written repaired-DX Stage 2B metrics;
#   - reconstructs exactly one macro row per multiclass endpoint/variant/fold
#     from the correct class rows;
#   - preserves binary surgery rows unchanged;
#   - backs up the buggy metrics CSV;
#   - rewrites the canonical metrics and downstream files 08-11;
#   - verifies uniqueness before every merge.
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

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage2_predictor_adjudication",
  "gpu_ablation_REPAIRED_DX_v1"
)

METRICS_FILE <- file.path(
  OUT_DIR,
  "05_predictor_adjudication_metrics.csv"
)

if (!file.exists(METRICS_FILE)) {
  stop(
    "Stage 2B-R metrics file not found:\n  ",
    METRICS_FILE,
    "\nThe GPU fits must finish through creation of file 05 before this repair can run.",
    call. = FALSE
  )
}

metrics_raw <- fread(
  METRICS_FILE
)

required_cols <- c(
  "endpoint",
  "variant",
  "evaluation_year",
  "metric_level",
  "class",
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
  "multiclass_brier",
  "best_iteration"
)

missing_cols <- setdiff(
  required_cols,
  names(metrics_raw)
)

if (length(missing_cols) > 0L) {
  stop(
    "Metrics file is missing expected column(s): ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. Confirm this is the known macro-row replication bug
# ------------------------------------------------------------------------------

macro_bug_qc <- metrics_raw[
  metric_level == "macro",
  .(
    macro_rows = .N,
    unique_class_values =
      uniqueN(class),
    unique_best_iterations =
      uniqueN(best_iteration)
  ),
  by = .(
    endpoint,
    variant,
    evaluation_year
  )
]

fwrite(
  macro_bug_qc,
  file.path(
    OUT_DIR,
    "05a_BUGGY_MACRO_ROW_QC.csv"
  )
)

# The original screen contains 135 fits:
#   45 endpoint x variant combinations across 3 temporal folds.
expected_fit_keys <- 135L

fit_key_count <- uniqueN(
  metrics_raw[
    ,
    .(
      endpoint,
      variant,
      evaluation_year
    )
  ]
)

if (fit_key_count != expected_fit_keys) {
  stop(
    "Expected 135 completed endpoint/variant/fold fit keys; observed ",
    fit_key_count,
    ". This is not just the known macro-summary bug; inspect the run before proceeding.",
    call. = FALSE
  )
}

# Every fit must have either class rows (multiclass) or one binary row.
fit_structure <- metrics_raw[
  ,
  .(
    class_rows =
      sum(metric_level == "class"),
    binary_rows =
      sum(metric_level == "binary"),
    macro_rows =
      sum(metric_level == "macro")
  ),
  by = .(
    endpoint,
    variant,
    evaluation_year
  )
]

bad_structure <- fit_structure[
  !(
    (
      class_rows > 0L &
        binary_rows == 0L
    ) |
      (
        class_rows == 0L &
          binary_rows == 1L
      )
  )
]

if (nrow(bad_structure) > 0L) {
  fwrite(
    bad_structure,
    file.path(
      OUT_DIR,
      "05b_UNEXPECTED_METRIC_STRUCTURE.csv"
    )
  )

  stop(
    "Unexpected metric-row structure detected. See 05b_UNEXPECTED_METRIC_STRUCTURE.csv.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 2. Reconstruct exactly one macro row from correct class-specific rows
# ------------------------------------------------------------------------------

class_rows <- metrics_raw[
  metric_level == "class"
]

binary_rows <- metrics_raw[
  metric_level == "binary"
]

macro_fixed <- class_rows[
  ,
  .(
    metric_level = "macro",
    class = "MACRO",
    n = max(n),
    prevalence = NA_real_,
    auroc =
      mean(
        auroc,
        na.rm = TRUE
      ),
    auprc =
      mean(
        auprc,
        na.rm = TRUE
      ),
    pr_lift =
      mean(
        pr_lift,
        na.rm = TRUE
      ),
    brier =
      mean(
        brier,
        na.rm = TRUE
      ),
    logloss =
      mean(
        logloss,
        na.rm = TRUE
      ),
    calibration_intercept =
      median(
        calibration_intercept,
        na.rm = TRUE
      ),
    calibration_slope =
      median(
        calibration_slope,
        na.rm = TRUE
      ),
    multiclass_logloss =
      unique(
        multiclass_logloss[
          !is.na(
            multiclass_logloss
          )
        ]
      )[1L],
    multiclass_brier =
      unique(
        multiclass_brier[
          !is.na(
            multiclass_brier
          )
        ]
      )[1L],
    best_iteration =
      unique(
        best_iteration
      )[1L]
  ),
  by = .(
    endpoint,
    variant,
    evaluation_year
  )
]

metrics_fixed <- rbindlist(
  list(
    class_rows,
    macro_fixed,
    binary_rows
  ),
  fill = TRUE,
  use.names = TRUE
)

setorder(
  metrics_fixed,
  variant,
  endpoint,
  evaluation_year,
  metric_level,
  class
)

# One headline row per completed fit.
headline <- metrics_fixed[
  metric_level %in%
    c(
      "macro",
      "binary"
    )
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
  nrow(headline_key_qc) !=
    expected_fit_keys ||
    any(
      headline_key_qc$N !=
        1L
    )
) {
  fwrite(
    headline_key_qc,
    file.path(
      OUT_DIR,
      "05c_HEADLINE_KEY_QC_FAILED.csv"
    )
  )

  stop(
    "Headline uniqueness check failed. No canonical files were overwritten.",
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 3. Backup the buggy canonical metrics, then replace with corrected metrics
# ------------------------------------------------------------------------------

backup_file <- file.path(
  OUT_DIR,
  "05_predictor_adjudication_metrics_BUGGY_MACRO_BACKUP.csv"
)

if (!file.exists(backup_file)) {
  ok <- file.copy(
    METRICS_FILE,
    backup_file,
    overwrite = FALSE
  )

  if (!isTRUE(ok)) {
    stop(
      "Could not create backup of buggy metrics file. Aborting before overwrite.",
      call. = FALSE
    )
  }
}

fwrite(
  metrics_fixed,
  METRICS_FILE
)

fwrite(
  headline,
  file.path(
    OUT_DIR,
    "05d_CORRECTED_HEADLINE_METRICS.csv"
  )
)

# ------------------------------------------------------------------------------
# 4. Rebuild foldwise deltas against clinical backbone
# ------------------------------------------------------------------------------

reference <- headline[
  variant == "clinical_backbone",
  .(
    endpoint,
    evaluation_year,
    ref_auroc = auroc,
    ref_auprc = auprc,
    ref_pr_lift = pr_lift,
    ref_brier = brier,
    ref_logloss =
      fifelse(
        !is.na(
          multiclass_logloss
        ),
        multiclass_logloss,
        logloss
      ),
    ref_cal_intercept =
      calibration_intercept,
    ref_cal_slope =
      calibration_slope
  )
]

reference_key_qc <- reference[
  ,
  .N,
  by = .(
    endpoint,
    evaluation_year
  )
]

if (
  nrow(reference_key_qc) != 15L ||
    any(
      reference_key_qc$N !=
        1L
    )
) {
  fwrite(
    reference_key_qc,
    file.path(
      OUT_DIR,
      "05e_REFERENCE_KEY_QC_FAILED.csv"
    )
  )

  stop(
    "Clinical-backbone reference is not unique by endpoint/year.",
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
    "Delta merge changed row count: headline=",
    nrow(headline),
    ", compare=",
    nrow(compare),
    ".",
    call. = FALSE
  )
}

compare[
  ,
  observed_logloss :=
    fifelse(
      !is.na(
        multiclass_logloss
      ),
      multiclass_logloss,
      logloss
    )
]

compare[, `:=`(
  delta_auroc =
    auroc -
      ref_auroc,
  delta_auprc =
    auprc -
      ref_auprc,
  delta_pr_lift =
    pr_lift -
      ref_pr_lift,
  delta_brier =
    brier -
      ref_brier,
  delta_logloss =
    observed_logloss -
      ref_logloss,
  delta_abs_calibration_intercept =
    abs(
      calibration_intercept
    ) -
      abs(
        ref_cal_intercept
      ),
  delta_abs_calibration_slope_error =
    abs(
      calibration_slope -
        1
    ) -
      abs(
        ref_cal_slope -
          1
      )
)]

fwrite(
  compare,
  file.path(
    OUT_DIR,
    "08_foldwise_deltas_vs_clinical_backbone.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. Predictor-family summaries
# ------------------------------------------------------------------------------

safe_median <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  median(
    x,
    na.rm = TRUE
  )
}

safe_min <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }

  min(
    x,
    na.rm = TRUE
  )
}

summary <- compare[
  variant != "clinical_backbone",
  .(
    n_temporal_evaluations = .N,
    median_delta_auroc =
      safe_median(
        delta_auroc
      ),
    worst_delta_auroc =
      safe_min(
        delta_auroc
      ),
    median_delta_auprc =
      safe_median(
        delta_auprc
      ),
    worst_delta_auprc =
      safe_min(
        delta_auprc
      ),
    median_delta_brier =
      safe_median(
        delta_brier
      ),
    median_delta_logloss =
      safe_median(
        delta_logloss
      ),
    median_delta_abs_calibration_intercept =
      safe_median(
        delta_abs_calibration_intercept
      ),
    median_delta_abs_calibration_slope_error =
      safe_median(
        delta_abs_calibration_slope_error
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
    "09_predictor_family_delta_summary.csv"
  )
)

overall_summary <- summary[
  ,
  .(
    endpoints_tested = .N,
    temporal_evaluations =
      sum(
        n_temporal_evaluations
      ),
    median_of_endpoint_delta_auroc =
      safe_median(
        median_delta_auroc
      ),
    worst_endpoint_delta_auroc =
      safe_min(
        worst_delta_auroc
      ),
    median_of_endpoint_delta_auprc =
      safe_median(
        median_delta_auprc
      ),
    median_of_endpoint_delta_logloss =
      safe_median(
        median_delta_logloss
      ),
    median_of_endpoint_delta_brier =
      safe_median(
        median_delta_brier
      )
  ),
  by = variant
]

fwrite(
  overall_summary,
  file.path(
    OUT_DIR,
    "10_predictor_family_overall_summary.csv"
  )
)

# ------------------------------------------------------------------------------
# 6. Class-specific deltas
# ------------------------------------------------------------------------------

class_ref <- metrics_fixed[
  variant ==
    "clinical_backbone" &
    metric_level ==
      "class",
  .(
    endpoint,
    evaluation_year,
    class,
    ref_class_auroc =
      auroc,
    ref_class_auprc =
      auprc,
    ref_class_brier =
      brier,
    ref_class_cal_intercept =
      calibration_intercept,
    ref_class_cal_slope =
      calibration_slope
  )
]

class_ref_qc <- class_ref[
  ,
  .N,
  by = .(
    endpoint,
    evaluation_year,
    class
  )
]

if (any(class_ref_qc$N != 1L)) {
  stop(
    "Class-specific clinical-backbone reference is not unique.",
    call. = FALSE
  )
}

class_compare <- merge(
  metrics_fixed[
    metric_level ==
      "class"
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
    auroc -
      ref_class_auroc,
  delta_class_auprc =
    auprc -
      ref_class_auprc,
  delta_class_brier =
    brier -
      ref_class_brier
)]

fwrite(
  class_compare,
  file.path(
    OUT_DIR,
    "11_class_specific_deltas.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Final postprocessing QC
# ------------------------------------------------------------------------------

qc <- data.table(
  check = c(
    "Completed fit keys",
    "Corrected headline rows",
    "Clinical-backbone reference rows",
    "Foldwise delta rows",
    "Corrected macro rows",
    "Binary rows",
    "Class rows"
  ),
  observed = c(
    fit_key_count,
    nrow(headline),
    nrow(reference),
    nrow(compare),
    nrow(macro_fixed),
    nrow(binary_rows),
    nrow(class_rows)
  ),
  expected = c(
    135L,
    135L,
    15L,
    135L,
    NA_integer_,
    NA_integer_,
    NA_integer_
  )
)

qc[
  ,
  passed :=
    is.na(expected) |
      observed ==
        expected
]

fwrite(
  qc,
  file.path(
    OUT_DIR,
    "12_POSTPROCESSING_REPAIR_QC.csv"
  )
)

if (!all(qc$passed)) {
  print(qc)
  stop(
    "Postprocessing QC failed.",
    call. = FALSE
  )
}

writeLines(
  c(
    "SPINE-TRACT Stage 2B-R repaired-DX screen: postprocessing repair complete.",
    "",
    "All 135 GPU fits had already completed before the original script failed.",
    "No model was refit by this repair script.",
    "",
    "Root cause:",
    "The multiclass macro calculation recycled endpoint/variant columns and produced K duplicate macro rows per fit.",
    "Class-specific rows, binary rows, model objects, and prediction caches were unaffected.",
    "",
    "Canonical 05_predictor_adjudication_metrics.csv has been replaced with corrected metrics.",
    "The original buggy file is preserved as 05_predictor_adjudication_metrics_BUGGY_MACRO_BACKUP.csv.",
    "",
    "Review next:",
    "09_predictor_family_delta_summary.csv",
    "10_predictor_family_overall_summary.csv",
    "11_class_specific_deltas.csv",
    "12_POSTPROCESSING_REPAIR_QC.csv"
  ),
  file.path(
    OUT_DIR,
    "README_POSTPROCESSING_REPAIR.txt"
  )
)

cat(
  "\n============================================================\n"
)
cat(
  "STAGE 2B-R POSTPROCESSING REPAIR COMPLETE\n"
)
cat(
  "============================================================\n"
)
cat(
  "No GPU models were refit.\n"
)
cat(
  "Corrected headline rows: ",
  nrow(headline),
  " / expected 135\n",
  sep = ""
)
cat(
  "Clinical-backbone reference rows: ",
  nrow(reference),
  " / expected 15\n",
  sep = ""
)
cat(
  "Foldwise delta rows: ",
  nrow(compare),
  " / expected 135\n\n",
  sep = ""
)
cat(
  "Return the entire folder:\n  ",
  OUT_DIR,
  "\n",
  sep = ""
)
