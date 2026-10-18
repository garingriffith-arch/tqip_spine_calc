# ==============================================================================
# 10s_repair_SPINE_TRACT_pooled_OOF_AUROC_REPAIRED_DX_v1.R
# SPINE-TRACT manuscript-output POSTPROCESSING REPAIR ONLY
#
# PURPOSE
#   Repair the pooled 2022-2024 OOF AUROC values in
#   02_CLASSIFICATION_POOLED_OOF_POINT_ESTIMATES.csv.
#
# ROOT CAUSE
#   The original pooled AUROC helper used integer event/non-event counts.
#   With ~265,000 pooled observations, n1*n0 exceeds the 32-bit integer range.
#   That overflows the denominator and produces NaN AUROC values.
#
#   Year-specific 2022/2023/2024 AUROCs and the 2024 bootstrap CIs are NOT
#   affected because those evaluation-year sample sizes remain below the
#   overflow threshold for the reported outputs.
#
# THIS SCRIPT
#   - DOES NOT refit any model;
#   - DOES NOT rerun bootstrap;
#   - reads the already-locked prediction caches;
#   - recalculates pooled 2022-2024 metrics using DOUBLE counts;
#   - backs up the buggy pooled table;
#   - rewrites the canonical pooled table;
#   - writes a QC file.
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

CLASS_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification_REPAIRED_DX",
  "rolling_origin_gpu"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage5_manuscript_outputs_REPAIRED_DX",
  "frozen_predictions_v1"
)

POOL_FILE <- file.path(
  OUT_DIR,
  "02_CLASSIFICATION_POOLED_OOF_POINT_ESTIMATES.csv"
)

if (!file.exists(POOL_FILE)) {
  stop(
    "Canonical pooled manuscript table not found:\n  ",
    POOL_FILE,
    call. = FALSE
  )
}

clamp_prob <- function(p, eps = 1e-7) {
  pmin(
    pmax(
      as.numeric(p),
      eps
    ),
    1 - eps
  )
}

fast_auc_double <- function(y, p) {
  ok <- !is.na(y) &
    !is.na(p) &
    is.finite(p)

  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  n1 <- as.double(
    sum(y == 1L)
  )

  n0 <- as.double(
    sum(y == 0L)
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

  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

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

  yy <- y[ord]

  precision <- cumsum(
    yy == 1L
  ) /
    seq_along(yy)

  mean(
    precision[
      yy == 1L
    ]
  )
}

binary_point_metrics <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  pp <- clamp_prob(p)
  prev <- mean(
    y == 1L
  )

  fit_slope <- tryCatch(
    suppressWarnings(
      glm(
        y ~ qlogis(pp),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  fit_int <- tryCatch(
    suppressWarnings(
      glm(
        y ~ offset(qlogis(pp)),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  data.table(
    n = length(y),
    prevalence = prev,
    auroc =
      fast_auc_double(
        y,
        p
      ),
    auprc =
      fast_auprc(
        y,
        p
      ),
    pr_lift =
      if (prev > 0) {
        fast_auprc(
          y,
          p
        ) /
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
      if (is.null(fit_int)) {
        NA_real_
      } else {
        unname(
          coef(fit_int)[1L]
        )
      },
    calibration_slope =
      if (is.null(fit_slope)) {
        NA_real_
      } else {
        unname(
          coef(fit_slope)[2L]
        )
      }
  )
}

load_binary_output <- function(year, id) {
  pred_dir <- file.path(
    CLASS_DIR,
    "predictions"
  )

  if (
    id %in%
      c(
        "home",
        "postacute",
        "death"
      )
  ) {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "discharge3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    if (id == "home") {
      return(
        data.table(
          admission_year =
            d$admission_year,
          inc_key =
            d$inc_key,
          output_id =
            "disposition_home",
          output_label =
            "Home/home health",
          y =
            as.integer(
              d$observed_class ==
                "Home/home health"
            ),
          p =
            d$p_Home.home.health
        )
      )
    }

    if (id == "postacute") {
      return(
        data.table(
          admission_year =
            d$admission_year,
          inc_key =
            d$inc_key,
          output_id =
            "disposition_postacute",
          output_label =
            "Post-acute facility",
          y =
            as.integer(
              d$observed_class ==
                "Post-acute facility"
            ),
          p =
            d$p_Post.acute.facility
        )
      )
    }

    return(
      data.table(
        admission_year =
          d$admission_year,
        inc_key =
          d$inc_key,
        output_id =
          "disposition_death_hospice",
        output_label =
          "Death/hospice",
        y =
          as.integer(
            d$observed_class ==
              "Death/hospice"
          ),
        p =
          d$p_Death.hospice
      )
    )
  }

  if (id == "hlos20") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "hlos3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    return(
      data.table(
        admission_year =
          d$admission_year,
        inc_key =
          d$inc_key,
        output_id =
          "hlos_ge20",
        output_label =
          "Hospital LOS >=20 days",
        y =
          as.integer(
            d$observed_class ==
              ">=20 days"
          ),
        p =
          d$p_X..20.days
      )
    )
  }

  if (
    id %in%
      c(
        "icu_any",
        "icu8"
      )
  ) {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "icu3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    if (id == "icu_any") {
      return(
        data.table(
          admission_year =
            d$admission_year,
          inc_key =
            d$inc_key,
          output_id =
            "icu_any",
          output_label =
            "Any ICU use",
          y =
            as.integer(
              d$observed_class !=
                "No ICU"
            ),
          p =
            d$p_ICU.1.7.days +
              d$p_ICU...8.days
        )
      )
    }

    return(
      data.table(
        admission_year =
          d$admission_year,
        inc_key =
          d$inc_key,
        output_id =
          "icu_ge8",
        output_label =
          "ICU >=8 days",
        y =
          as.integer(
            d$observed_class ==
              "ICU >=8 days"
          ),
        p =
          d$p_ICU...8.days
      )
    )
  }

  if (
    id %in%
      c(
        "vent_any",
        "vent8"
      )
  ) {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "vent3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    if (id == "vent_any") {
      return(
        data.table(
          admission_year =
            d$admission_year,
          inc_key =
            d$inc_key,
          output_id =
            "vent_any",
          output_label =
            "Any invasive mechanical ventilation",
          y =
            as.integer(
              d$observed_class !=
                "No ventilation"
            ),
          p =
            d$p_Ventilation.1.7.days +
              d$p_Ventilation...8.days
        )
      )
    }

    return(
      data.table(
        admission_year =
          d$admission_year,
        inc_key =
          d$inc_key,
        output_id =
          "vent_ge8",
        output_label =
          "Ventilation >=8 days",
        y =
          as.integer(
            d$observed_class ==
              "Ventilation >=8 days"
          ),
        p =
          d$p_Ventilation...8.days
      )
    )
  }

  if (id == "surgery") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "operative_spine_intervention__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    return(
      data.table(
        admission_year =
          d$admission_year,
        inc_key =
          d$inc_key,
        output_id =
          "recorded_tl_operation",
        output_label =
          "Recorded thoracolumbar operative treatment",
        y =
          as.integer(
            d$y
          ),
        p =
          as.numeric(
            d$pred
          )
      )
    )
  }

  stop(
    "Unknown output id: ",
    id,
    call. = FALSE
  )
}

output_ids <- c(
  "home",
  "postacute",
  "death",
  "hlos20",
  "icu_any",
  "icu8",
  "vent_any",
  "vent8",
  "surgery"
)

rows <- list()

for (i in seq_along(output_ids)) {
  id <- output_ids[i]

  per_year <- lapply(
    2022:2024,
    function(year) {
      load_binary_output(
        year,
        id
      )
    }
  )

  pooled <- rbindlist(
    per_year,
    use.names = TRUE
  )

  m <- binary_point_metrics(
    pooled$y,
    pooled$p
  )

  m[, `:=`(
    output_id =
      pooled$output_id[1L],
    output_label =
      pooled$output_label[1L],
    evaluation_year =
      "Pooled 2022-2024"
  )]

  rows[[i]] <- m
}

corrected <- rbindlist(
  rows,
  fill = TRUE
)

setcolorder(
  corrected,
  c(
    "n",
    "prevalence",
    "auroc",
    "auprc",
    "pr_lift",
    "brier",
    "logloss",
    "calibration_intercept",
    "calibration_slope",
    "output_id",
    "output_label",
    "evaluation_year"
  )
)

if (
  nrow(corrected) !=
    9L ||
    any(
      !is.finite(
        corrected$auroc
      )
    )
) {
  print(corrected)

  stop(
    "Pooled AUROC repair QC failed.",
    call. = FALSE
  )
}

backup <- file.path(
  OUT_DIR,
  "02_CLASSIFICATION_POOLED_OOF_POINT_ESTIMATES_BUGGY_AUROC_BACKUP.csv"
)

if (!file.exists(backup)) {
  ok <- file.copy(
    POOL_FILE,
    backup,
    overwrite = FALSE
  )

  if (!isTRUE(ok)) {
    stop(
      "Could not back up original pooled table.",
      call. = FALSE
    )
  }
}

fwrite(
  corrected,
  POOL_FILE
)

fwrite(
  corrected[
    ,
    .(
      output_id,
      output_label,
      n,
      prevalence,
      auroc
    )
  ],
  file.path(
    OUT_DIR,
    "02a_POOLED_AUROC_REPAIR_QC.csv"
  )
)

cat(
  "\n============================================================\n"
)
cat(
  "POOLED OOF AUROC REPAIR COMPLETE — NO MODEL REFIT\n"
)
cat(
  "============================================================\n\n"
)
print(
  corrected[
    ,
    .(
      output_label,
      n,
      prevalence,
      auroc,
      auprc
    )
  ]
)
