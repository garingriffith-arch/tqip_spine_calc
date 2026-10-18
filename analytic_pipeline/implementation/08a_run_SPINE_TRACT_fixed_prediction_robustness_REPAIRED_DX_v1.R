# ==============================================================================
# 08a_run_SPINE_TRACT_fixed_prediction_robustness_REPAIRED_DX_v1.R
# SPINE-TRACT repaired-DX robustness suite using LOCKED prediction caches
#
# PURPOSE
#   Rebuild manuscript-facing robustness analyses after the diagnosis repair
#   WITHOUT refitting the final models.
#
# INCLUDED
#   1. 2022-2024 overall fixed-prediction performance.
#   2. 2024 + temporal subgroup robustness:
#        - GCS severity
#        - age
#        - sex
#        - clinical-input missingness
#        - thoracic vs lumbar vs both-region injury
#        - fracture-only vs neural injury
#        - SCI vs no SCI
#        - cord vs root/cauda neural pattern
#        - concomitant cervical spine injury
#        - concomitant TBI
#        - race / ethnicity / payer (descriptive subgroup calibration)
#   3. Ventilation performance by baseline airway state.
#   4. Survivor-only and early-death exclusion sensitivities.
#   5. Continuous-duration subgroup / survivor / early-death sensitivity.
#   6. 2024 broad-phenotype classification sensitivity by applying the LOCKED
#      2024 narrow-cohort models to the repaired broad-union dataset.
#   7. Temporal-era performance deltas (2022->2024; 2023->2024).
#
# NOT INCLUDED / INTENTIONALLY NOT REFIT HERE
#   - transfer-in model sensitivity: primary use case is locked to direct
#     presentations; transfer-in rows are outside the frozen modeling dataset.
#   - death-or-prolonged composite model retuning: separate heavy sensitivity,
#     not necessary to regenerate the locked primary manuscript estimates.
#   - facility heterogeneity: only if usable facility identifiers are available.
#
# IMPORTANT
#   - NO primary model is tuned or refit.
#   - Broad-cohort analysis is application of the locked 2024 narrow-trained
#     classification model, not a separately trained broad-cohort model.
#   - Subgroup performance is descriptive and should not be labeled "fairness"
#     from a single metric alone.
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

DATA_DIR <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication"
)

PRIMARY_FILE <- file.path(
  DATA_DIR,
  "spine_tract_frozen_primary_2020_2024.rds"
)

UNION_FILE <- file.path(
  DATA_DIR,
  "spine_tract_frozen_candidate_union_2020_2024.rds"
)

SPEC_FILE <- file.path(
  DATA_DIR,
  "spine_tract_predictor_sets.rds"
)

CLASS_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification_REPAIRED_DX",
  "rolling_origin_gpu"
)

DURATION_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_duration_quantiles_REPAIRED_DX",
  "rolling_origin_gpu"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage4_robustness_REPAIRED_DX",
  "fixed_predictions_v1"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

needed <- c(
  PRIMARY_FILE,
  UNION_FILE,
  SPEC_FILE,
  file.path(
    CLASS_DIR,
    "05_FINAL_CLASSIFICATION_METRICS.csv"
  ),
  file.path(
    DURATION_DIR,
    "05_FINAL_DURATION_TEMPORAL_METRICS.csv"
  )
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

MIN_N <- 500L
MIN_EVENTS <- 25L
DURATION_TOL <- 1e-6

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

num <- function(x) {
  suppressWarnings(
    as.numeric(
      as.character(x)
    )
  )
}

norm_chr <- function(x) {
  z <- trimws(
    as.character(x)
  )

  z[
    is.na(z) |
      z == "" |
      z %in% c(
        "NA",
        "NaN",
        "<NA>"
      )
  ] <- "__UNKNOWN__"

  z
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
    1 - eps
  )
}

fast_auc <- function(y, p) {
  ok <- !is.na(y) &
    !is.na(p) &
    is.finite(p)

  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)

  if (
    n1 == 0L ||
      n0 == 0L
  ) {
    return(NA_real_)
  }

  r <- rank(
    p,
    ties.method = "average"
  )

  (
    sum(r[y == 1L]) -
      n1 * (n1 + 1) / 2
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

  npos <- sum(y == 1L)

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

calibration <- function(y, p) {
  ok <- !is.na(y) &
    !is.na(p) &
    is.finite(p)

  y <- as.integer(y[ok])
  p <- clamp_prob(p[ok])

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

  lp <- qlogis(p)

  fi <- tryCatch(
    suppressWarnings(
      glm(
        y ~ offset(lp),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  fs <- tryCatch(
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
      if (is.null(fi)) {
        NA_real_
      } else {
        unname(
          coef(fi)[1L]
        )
      },
    slope =
      if (is.null(fs)) {
        NA_real_
      } else {
        unname(
          coef(fs)[2L]
        )
      }
  )
}

binary_metrics <- function(y, p) {
  ok <- !is.na(y) &
    !is.na(p) &
    is.finite(p)

  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  if (length(y) == 0L) {
    return(data.table())
  }

  cal <- calibration(
    y,
    p
  )

  pp <- clamp_prob(p)
  prev <- mean(
    y == 1L
  )

  data.table(
    N = length(y),
    events = sum(
      y == 1L
    ),
    prevalence = prev,
    AUROC =
      fast_auc(
        y,
        p
      ),
    AUPRC =
      fast_auprc(
        y,
        p
      ),
    PR_lift =
      if (
        is.finite(prev) &&
          prev > 0
      ) {
        fast_auprc(
          y,
          p
        ) /
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

pinball <- function(
    y,
    q,
    alpha
) {
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
    tol = DURATION_TOL
) {
  ok <- is.finite(y) &
    is.finite(q10) &
    is.finite(q50) &
    is.finite(q90)

  y <- y[ok]
  q10 <- q10[ok]
  q50 <- q50[ok]
  q90 <- q90[ok]

  if (length(y) == 0L) {
    return(data.table())
  }

  data.table(
    N = length(y),
    MAE_days =
      mean(
        abs(
          y -
            q50
        )
      ),
    median_bias_days =
      median(
        q50 -
          y
      ),
    mean_pinball =
      mean(c(
        pinball(
          y,
          q10,
          0.10
        ),
        pinball(
          y,
          q50,
          0.50
        ),
        pinball(
          y,
          q90,
          0.90
        )
      )),
    p_y_lt_q10 =
      mean(
        y <
          q10 -
            tol
      ),
    p_y_le_q10 =
      mean(
        y <=
          q10 +
            tol
      ),
    p_y_lt_q50 =
      mean(
        y <
          q50 -
            tol
      ),
    p_y_le_q50 =
      mean(
        y <=
          q50 +
            tol
      ),
    p_y_lt_q90 =
      mean(
        y <
          q90 -
            tol
      ),
    p_y_le_q90 =
      mean(
        y <=
          q90 +
            tol
      ),
    q10_q90_coverage =
      mean(
        y +
          tol >=
          q10 &
          y -
            tol <=
          q90
      ),
    median_q10_q90_width =
      median(
        q90 -
          q10
      )
  )
}

safe_bind <- function(x) {
  if (length(x) == 0L) {
    return(data.table())
  }

  rbindlist(
    x,
    fill = TRUE,
    use.names = TRUE
  )
}

# ------------------------------------------------------------------------------
# 2. Data and corrected final predictor policy
# ------------------------------------------------------------------------------

dt <- as.data.table(
  readRDS(
    PRIMARY_FILE
  )
)

union <- as.data.table(
  readRDS(
    UNION_FILE
  )
)

spec <- readRDS(
  SPEC_FILE
)

dt[, inc_key := as.character(inc_key)]
union[, inc_key := as.character(inc_key)]

dt[, admission_year := as.integer(admission_year)]
union[, admission_year := as.integer(admission_year)]

if (nrow(dt) != 424017L) {
  stop(
    "Expected repaired primary N=424017; observed ",
    nrow(dt),
    ".",
    call. = FALSE
  )
}

if (nrow(union) != 440380L) {
  stop(
    "Expected repaired broad-union N=440380; observed ",
    nrow(union),
    ".",
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
            dt[[v]] ==
              1L,
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
# 3. Sensitivity strata
# ------------------------------------------------------------------------------

ge <- num(
  dt$gcs_eye_clean
)
gv <- num(
  dt$gcs_verbal_clean
)
gm <- num(
  dt$gcs_motor_clean
)

gt <- rep(
  NA_real_,
  nrow(dt)
)

gok <- !is.na(ge) &
  !is.na(gv) &
  !is.na(gm)

gt[gok] <-
  ge[gok] +
  gv[gok] +
  gm[gok]

dt[
  ,
  gcs_severity_sens :=
    fcase(
      !is.na(gt) &
        gt >= 13 &
        gt <= 15,
      "Mild GCS 13-15",

      !is.na(gt) &
        gt >= 9 &
        gt <= 12,
      "Moderate GCS 9-12",

      !is.na(gt) &
        gt >= 3 &
        gt <= 8,
      "Severe GCS 3-8",

      default =
        "GCS unknown/unclassifiable"
    )
]

agev <- num(
  dt$age_years
)

dt[
  ,
  age_group_sens :=
    fcase(
      !is.na(agev) &
        agev >= 18 &
        agev <= 39,
      "18-39",

      !is.na(agev) &
        agev >= 40 &
        agev <= 64,
      "40-64",

      !is.na(agev) &
        agev >= 65 &
        agev <= 79,
      "65-79",

      !is.na(agev) &
        agev >= 80 &
        agev <= 89,
      "80-89",

      default =
        "Age unknown"
    )
]

# Clinical-input missingness uses final_core only; social/context variables are
# intentionally excluded from the missingness burden definition.
miss_count <- integer(
  nrow(dt)
)

for (v in final_core) {
  x <- dt[[v]]

  if (
    is.numeric(x) ||
      is.integer(x)
  ) {
    miss_count <-
      miss_count +
      as.integer(
        is.na(x) |
          !is.finite(
            as.numeric(x)
          )
      )
  } else {
    z <- norm_chr(x)

    miss_count <-
      miss_count +
      as.integer(
        z %in%
          c(
            "__UNKNOWN__",
            "__OTHER__"
          )
      )
  }
}

dt[
  ,
  clinical_input_missing_count :=
    miss_count
]

dt[
  ,
  missingness_sens :=
    fcase(
      clinical_input_missing_count ==
        0L,
      "0 missing/unknown",

      clinical_input_missing_count ==
        1L,
      "1 missing/unknown",

      clinical_input_missing_count ==
        2L,
      "2 missing/unknown",

      clinical_input_missing_count >=
        3L,
      ">=3 missing/unknown",

      default =
        "Unknown"
    )
]

thoracic_region <-
  dt$dx_thoracic_vertebral_fracture ==
    1L |
  dt$dx_thoracic_spinal_cord_injury ==
    1L |
  dt$dx_thoracic_nerve_root_injury ==
    1L

lumbar_region <-
  dt$dx_lumbar_vertebral_fracture ==
    1L |
  dt$dx_lumbar_spinal_cord_injury ==
    1L |
  dt$dx_lumbar_sacral_nerve_root_injury ==
    1L |
  dt$dx_cauda_equina_injury ==
    1L

dt[
  ,
  spine_region_sens :=
    fcase(
      thoracic_region &
        lumbar_region,
      "Thoracic + lumbar",

      thoracic_region &
        !lumbar_region,
      "Thoracic only",

      !thoracic_region &
        lumbar_region,
      "Lumbar only",

      default =
        "Region unclassifiable"
    )
]

any_fracture <-
  dt$dx_thoracic_vertebral_fracture ==
    1L |
  dt$dx_lumbar_vertebral_fracture ==
    1L

any_neural <-
  dt$dx_tl_spinal_cord_injury_any ==
    1L |
  dt$dx_tl_root_cauda_any ==
    1L

dt[
  ,
  fracture_neural_sens :=
    fcase(
      any_fracture &
        !any_neural,
      "Fracture only",

      any_neural,
      "Neural injury present",

      default =
        "Other qualifying phenotype"
    )
]

dt[
  ,
  sci_sens :=
    fifelse(
      dx_tl_spinal_cord_injury_any ==
        1L,
      "SCI present",
      "No TL SCI"
    )
]

dt[
  ,
  neural_pattern_sens :=
    fcase(
      dx_tl_spinal_cord_injury_any ==
        1L &
        dx_tl_root_cauda_any ==
          1L,
      "Cord + root/cauda",

      dx_tl_spinal_cord_injury_any ==
        1L,
      "Cord injury",

      dx_tl_root_cauda_any ==
        1L,
      "Root/cauda injury",

      default =
        "No neural injury"
    )
]

dt[
  ,
  cervical_sens :=
    fifelse(
      dx_cervical_spine_injury ==
        1L,
      "Concomitant cervical injury",
      "No cervical injury"
    )
]

dt[
  ,
  tbi_sens :=
    fifelse(
      dx_tbi_any ==
        1L,
      "Concomitant TBI",
      "No TBI"
    )
]

intub <- num(
  dt$gcsq_intubated_recovered
)

resp <- norm_chr(
  dt$respiratory_assistance_clean
)

resp_assisted <- grepl(
  "assist",
  tolower(resp)
)

resp_known <- !resp %in%
  c(
    "__UNKNOWN__",
    "__OTHER__"
  )

dt[
  ,
  airway_state_sens :=
    fcase(
      intub ==
        1 |
        resp_assisted,
      "Baseline airway-positive",

      intub ==
        0 &
        resp_known &
        !resp_assisted,
      "Baseline airway-negative",

      default =
        "Baseline airway uncertain"
    )
]

dt[
  ,
  death_hospice :=
    as.integer(
      as.character(discharge4) ==
        "Death/hospice"
    )
]

dt[
  ,
  early_death_le1 :=
    as.integer(
      death_hospice ==
        1L &
        !is.na(hospital_los_days) &
        hospital_los_days <=
          1
    )
]

dt[
  ,
  early_death_le2 :=
    as.integer(
      death_hospice ==
        1L &
        !is.na(hospital_los_days) &
        hospital_los_days <=
          2
    )
]

stratum_counts <- rbindlist(
  list(
    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          gcs_severity_sens
      )
    ][
      ,
      subgroup_domain :=
        "GCS severity"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          age_group_sens
      )
    ][
      ,
      subgroup_domain :=
        "Age"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          norm_chr(
            sex_clean
          )
      )
    ][
      ,
      subgroup_domain :=
        "Sex"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          missingness_sens
      )
    ][
      ,
      subgroup_domain :=
        "Clinical input missingness"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          spine_region_sens
      )
    ][
      ,
      subgroup_domain :=
        "TL region"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          fracture_neural_sens
      )
    ][
      ,
      subgroup_domain :=
        "Fracture vs neural"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          sci_sens
      )
    ][
      ,
      subgroup_domain :=
        "SCI"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          neural_pattern_sens
      )
    ][
      ,
      subgroup_domain :=
        "Neural pattern"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          cervical_sens
      )
    ][
      ,
      subgroup_domain :=
        "Concomitant cervical injury"
    ],

    dt[
      ,
      .(N = .N),
      by = .(
        admission_year,
        subgroup_level =
          tbi_sens
      )
    ][
      ,
      subgroup_domain :=
        "Concomitant TBI"
    ]
  ),
  fill = TRUE
)

fwrite(
  stratum_counts,
  file.path(
    OUT_DIR,
    "01_SENSITIVITY_STRATUM_COUNTS.csv"
  )
)

# ------------------------------------------------------------------------------
# 4. Canonical classification prediction loader
# ------------------------------------------------------------------------------

load_class_prediction <- function(
    year,
    endpoint
) {
  f <- file.path(
    CLASS_DIR,
    "predictions",
    paste0(
      endpoint,
      "__eval",
      year,
      "__predictions.csv"
    )
  )

  if (!file.exists(f)) {
    stop(
      "Missing classification prediction cache: ",
      f,
      call. = FALSE
    )
  }

  d <- fread(f)
  d[, inc_key := as.character(inc_key)]
  d
}

class_output_specs <- data.table(
  output_id = c(
    "disposition_home",
    "disposition_postacute",
    "disposition_death_hospice",

    "hlos_le7",
    "hlos_8_19",
    "hlos_ge20",

    "icu_none",
    "icu_1_7",
    "icu_ge8",
    "icu_any",

    "vent_none",
    "vent_1_7",
    "vent_ge8",
    "vent_any",

    "recorded_tl_operation"
  ),
  output_label = c(
    "Home/home health",
    "Post-acute facility",
    "Death/hospice",

    "Hospital LOS <=7 days",
    "Hospital LOS 8-19 days",
    "Hospital LOS >=20 days",

    "No ICU",
    "ICU 1-7 days",
    "ICU >=8 days",
    "Any ICU use",

    "No ventilation",
    "Ventilation 1-7 days",
    "Ventilation >=8 days",
    "Any invasive mechanical ventilation",

    "Recorded thoracolumbar operative treatment"
  ),
  endpoint = c(
    rep(
      "discharge3_final",
      3
    ),
    rep(
      "hlos3_final",
      3
    ),
    rep(
      "icu3_final",
      4
    ),
    rep(
      "vent3_final",
      4
    ),
    "operative_spine_intervention"
  ),
  target = c(
    "Home/home health",
    "Post-acute facility",
    "Death/hospice",

    "<=7 days",
    "8-19 days",
    ">=20 days",

    "No ICU",
    "ICU 1-7 days",
    "ICU >=8 days",
    "ANY",

    "No ventilation",
    "Ventilation 1-7 days",
    "Ventilation >=8 days",
    "ANY",

    "EVENT"
  ),
  headline = c(
    TRUE,
    TRUE,
    TRUE,

    FALSE,
    FALSE,
    TRUE,

    FALSE,
    FALSE,
    TRUE,
    TRUE,

    FALSE,
    FALSE,
    TRUE,
    TRUE,

    TRUE
  )
)

extract_binary_output <- function(
    pred,
    spec_row
) {
  ep <- spec_row$endpoint
  target <- spec_row$target

  if (
    ep ==
      "operative_spine_intervention"
  ) {
    return(
      data.table(
        admission_year =
          as.integer(
            pred$admission_year
          ),
        inc_key =
          as.character(
            pred$inc_key
          ),
        y =
          as.integer(
            pred$y
          ),
        p =
          as.numeric(
            pred$pred
          )
      )
    )
  }

  obs <- as.character(
    pred$observed_class
  )

  if (
    ep ==
      "discharge3_final"
  ) {
    map <- c(
      "Home/home health" =
        "p_Home.home.health",
      "Post-acute facility" =
        "p_Post.acute.facility",
      "Death/hospice" =
        "p_Death.hospice"
    )

    p <- as.numeric(
      pred[[unname(
            map[target]
          )]]
    )

    y <- as.integer(
      obs ==
        target
    )
  } else if (
    ep ==
      "hlos3_final"
  ) {
    map <- c(
      "<=7 days" =
        "p_X..7.days",
      "8-19 days" =
        "p_X8.19.days",
      ">=20 days" =
        "p_X..20.days"
    )

    p <- as.numeric(
      pred[[unname(
            map[target]
          )]]
    )

    y <- as.integer(
      obs ==
        target
    )
  } else if (
    ep ==
      "icu3_final"
  ) {
    if (target == "ANY") {
      p <-
        as.numeric(
          pred$p_ICU.1.7.days
        ) +
        as.numeric(
          pred$p_ICU...8.days
        )

      y <- as.integer(
        obs !=
          "No ICU"
      )
    } else {
      map <- c(
        "No ICU" =
          "p_No.ICU",
        "ICU 1-7 days" =
          "p_ICU.1.7.days",
        "ICU >=8 days" =
          "p_ICU...8.days"
      )

      p <- as.numeric(
        pred[[unname(
              map[target]
            )]]
      )

      y <- as.integer(
        obs ==
          target
      )
    }
  } else if (
    ep ==
      "vent3_final"
  ) {
    if (target == "ANY") {
      p <-
        as.numeric(
          pred$p_Ventilation.1.7.days
        ) +
        as.numeric(
          pred$p_Ventilation...8.days
        )

      y <- as.integer(
        obs !=
          "No ventilation"
      )
    } else {
      map <- c(
        "No ventilation" =
          "p_No.ventilation",
        "Ventilation 1-7 days" =
          "p_Ventilation.1.7.days",
        "Ventilation >=8 days" =
          "p_Ventilation...8.days"
      )

      p <- as.numeric(
        pred[[unname(
              map[target]
            )]]
      )

      y <- as.integer(
        obs ==
          target
      )
    }
  } else {
    stop(
      "Unknown endpoint in output specification: ",
      ep,
      call. = FALSE
    )
  }

  data.table(
    admission_year =
      as.integer(
        pred$admission_year
      ),
    inc_key =
      as.character(
        pred$inc_key
      ),
    y =
      y,
    p =
      p
  )
}

# ------------------------------------------------------------------------------
# 5. Classification overall + subgroup robustness
# ------------------------------------------------------------------------------

overall_rows <- list()
subgroup_rows <- list()
early_rows <- list()
airway_rows <- list()

subgroup_defs <- list(
  "GCS severity" =
    "gcs_severity_sens",
  "Age" =
    "age_group_sens",
  "Sex" =
    "sex_clean",
  "Clinical input missingness" =
    "missingness_sens",
  "TL region" =
    "spine_region_sens",
  "Fracture vs neural" =
    "fracture_neural_sens",
  "SCI" =
    "sci_sens",
  "Neural pattern" =
    "neural_pattern_sens",
  "Concomitant cervical injury" =
    "cervical_sens",
  "Concomitant TBI" =
    "tbi_sens",
  "Race" =
    "race_clean",
  "Ethnicity" =
    "ethnicity_clean",
  "Payer" =
    "payer_clean"
)

for (year in 2022:2024) {
  year_data <- dt[
    admission_year ==
      year
  ]

  for (
    i in seq_len(
      nrow(
        class_output_specs
      )
    )
  ) {
    sp <- class_output_specs[i]

    pred <- load_class_prediction(
      year,
      sp$endpoint
    )

    d <- extract_binary_output(
      pred,
      sp
    )

    d <- merge(
      d,
      year_data,
      by = c(
        "admission_year",
        "inc_key"
      ),
      all.x = TRUE,
      sort = FALSE
    )

    if (
      nrow(d) !=
        nrow(pred)
    ) {
      stop(
        "Primary prediction-to-dataset merge mismatch for ",
        sp$output_id,
        " / ",
        year,
        ".",
        call. = FALSE
      )
    }

    m <- binary_metrics(
      d$y,
      d$p
    )

    m[, `:=`(
      evaluation_year =
        year,
      output_id =
        sp$output_id,
      output_label =
        sp$output_label,
      headline =
        sp$headline
    )]

    overall_rows[[length(overall_rows) + 1L]] <- m

    for (dom in names(subgroup_defs)) {
      v <- subgroup_defs[[dom]]
      g <- norm_chr(
        d[[v]]
      )

      for (lv in unique(g)) {
        idx <- which(
          g ==
            lv
        )

        if (
          length(idx) <
            MIN_N
        ) {
          next
        }

        yy <- d$y[idx]
        pp <- d$p[idx]

        if (
          sum(
            yy ==
              1L,
            na.rm = TRUE
          ) <
            MIN_EVENTS ||
            sum(
              yy ==
                0L,
              na.rm = TRUE
            ) <
              MIN_EVENTS
        ) {
          next
        }

        sm <- binary_metrics(
          yy,
          pp
        )

        sm[, `:=`(
          evaluation_year =
            year,
          output_id =
            sp$output_id,
          output_label =
            sp$output_label,
          headline =
            sp$headline,
          subgroup_domain =
            dom,
          subgroup_level =
            lv
        )]

        subgroup_rows[[length(subgroup_rows) + 1L]] <- sm
      }
    }

    # Resource / trajectory sensitivity to death truncation.
    if (
      sp$output_id %in%
        c(
          "hlos_ge20",
          "icu_ge8",
          "icu_any",
          "vent_ge8",
          "vent_any"
        )
    ) {
      sensitivity_filters <- list(
        "Overall" =
          rep(
            TRUE,
            nrow(d)
          ),
        "Exclude death/hospice" =
          d$death_hospice !=
            1L,
        "Exclude death/hospice <=1 day" =
          d$early_death_le1 !=
            1L,
        "Exclude death/hospice <=2 days" =
          d$early_death_le2 !=
            1L
      )

      for (
        nm in names(
          sensitivity_filters
        )
      ) {
        idx <- which(
          sensitivity_filters[[nm]]
        )

        yy <- d$y[idx]
        pp <- d$p[idx]

        if (
          length(idx) <
            MIN_N ||
            sum(
              yy ==
                1L
            ) <
              MIN_EVENTS ||
            sum(
              yy ==
                0L
            ) <
              MIN_EVENTS
        ) {
          next
        }

        em <- binary_metrics(
          yy,
          pp
        )

        em[, `:=`(
          evaluation_year =
            year,
          output_id =
            sp$output_id,
          output_label =
            sp$output_label,
          sensitivity_population =
            nm
        )]

        early_rows[[length(early_rows) + 1L]] <- em
      }
    }

    if (
      sp$output_id %in%
        c(
          "vent_none",
          "vent_1_7",
          "vent_ge8",
          "vent_any"
        )
    ) {
      g <- d$airway_state_sens

      for (lv in unique(g)) {
        idx <- which(
          g ==
            lv
        )

        yy <- d$y[idx]
        pp <- d$p[idx]

        if (
          length(idx) <
            MIN_N ||
            sum(
              yy ==
                1L
            ) <
              MIN_EVENTS ||
            sum(
              yy ==
                0L
            ) <
              MIN_EVENTS
        ) {
          next
        }

        am <- binary_metrics(
          yy,
          pp
        )

        am[, `:=`(
          evaluation_year =
            year,
          output_id =
            sp$output_id,
          output_label =
            sp$output_label,
          airway_state =
            lv
        )]

        airway_rows[[length(airway_rows) + 1L]] <- am
      }
    }
  }
}

overall <- safe_bind(
  overall_rows
)

subgroups <- safe_bind(
  subgroup_rows
)

early_sens <- safe_bind(
  early_rows
)

airway_sens <- safe_bind(
  airway_rows
)

fwrite(
  overall,
  file.path(
    OUT_DIR,
    "02_CLASSIFICATION_OVERALL_TEMPORAL.csv"
  )
)

fwrite(
  subgroups,
  file.path(
    OUT_DIR,
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv"
  )
)

fwrite(
  early_sens,
  file.path(
    OUT_DIR,
    "04_EARLY_DEATH_RESOURCE_SENSITIVITY.csv"
  )
)

fwrite(
  airway_sens,
  file.path(
    OUT_DIR,
    "05_VENTILATION_AIRWAY_STATE_SENSITIVITY.csv"
  )
)

# Compact 2024 subgroup range summary.
subgroup_summary_2024 <- subgroups[
  evaluation_year ==
    2024L &
    headline ==
      TRUE,
  .(
    n_strata =
      .N,
    AUROC_min =
      min(
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
    subgroup_domain,
    output_id,
    output_label
  )
]

fwrite(
  subgroup_summary_2024,
  file.path(
    OUT_DIR,
    "06_2024_SUBGROUP_RANGE_SUMMARY.csv"
  )
)

# ------------------------------------------------------------------------------
# 6. Continuous-duration robustness from locked prediction caches
# ------------------------------------------------------------------------------

duration_ids <- c(
  "hospital_los",
  "icu_los_conditional",
  "ventilator_days_conditional"
)

duration_overall_rows <- list()
duration_subgroup_rows <- list()
duration_early_rows <- list()
duration_airway_rows <- list()

duration_subgroup_defs <- list(
  "GCS severity" =
    "gcs_severity_sens",
  "Age" =
    "age_group_sens",
  "Sex" =
    "sex_clean",
  "Clinical input missingness" =
    "missingness_sens",
  "TL region" =
    "spine_region_sens",
  "Fracture vs neural" =
    "fracture_neural_sens",
  "SCI" =
    "sci_sens",
  "Neural pattern" =
    "neural_pattern_sens",
  "Concomitant cervical injury" =
    "cervical_sens",
  "Concomitant TBI" =
    "tbi_sens"
)

for (year in 2022:2024) {
  year_data <- dt[
    admission_year ==
      year
  ]

  for (id in duration_ids) {
    f <- file.path(
      DURATION_DIR,
      "predictions",
      paste0(
        id,
        "__eval",
        year,
        "__Q10_Q50_Q90_predictions.csv"
      )
    )

    if (!file.exists(f)) {
      stop(
        "Missing duration prediction cache: ",
        f,
        call. = FALSE
      )
    }

    d <- fread(f)
    d[, inc_key := as.character(inc_key)]

    d <- merge(
      d,
      year_data,
      by = c(
        "admission_year",
        "inc_key"
      ),
      all.x = TRUE,
      sort = FALSE
    )

    dm <- duration_metrics(
      y =
        as.numeric(
          d$observed_days
        ),
      q10 =
        as.numeric(
          d$q10_days
        ),
      q50 =
        as.numeric(
          d$q50_median_days
        ),
      q90 =
        as.numeric(
          d$q90_days
        )
    )

    dm[, `:=`(
      evaluation_year =
        year,
      outcome =
        id
    )]

    duration_overall_rows[[length(duration_overall_rows) + 1L]] <- dm

    for (dom in names(duration_subgroup_defs)) {
      v <- duration_subgroup_defs[[dom]]
      g <- norm_chr(
        d[[v]]
      )

      for (lv in unique(g)) {
        idx <- which(
          g ==
            lv
        )

        if (
          length(idx) <
            MIN_N
        ) {
          next
        }

        sm <- duration_metrics(
          y =
            as.numeric(
              d$observed_days[idx]
            ),
          q10 =
            as.numeric(
              d$q10_days[idx]
            ),
          q50 =
            as.numeric(
              d$q50_median_days[idx]
            ),
          q90 =
            as.numeric(
              d$q90_days[idx]
            )
        )

        sm[, `:=`(
          evaluation_year =
            year,
          outcome =
            id,
          subgroup_domain =
            dom,
          subgroup_level =
            lv
        )]

        duration_subgroup_rows[[length(duration_subgroup_rows) + 1L]] <- sm
      }
    }

    sensitivity_filters <- list(
      "Overall" =
        rep(
          TRUE,
          nrow(d)
        ),
      "Exclude death/hospice" =
        d$death_hospice !=
          1L,
      "Exclude death/hospice <=1 day" =
        d$early_death_le1 !=
          1L,
      "Exclude death/hospice <=2 days" =
        d$early_death_le2 !=
          1L
    )

    for (
      nm in names(
        sensitivity_filters
      )
    ) {
      idx <- which(
        sensitivity_filters[[nm]]
      )

      if (
        length(idx) <
          MIN_N
      ) {
        next
      }

      em <- duration_metrics(
        y =
          as.numeric(
            d$observed_days[idx]
          ),
        q10 =
          as.numeric(
            d$q10_days[idx]
          ),
        q50 =
          as.numeric(
            d$q50_median_days[idx]
          ),
        q90 =
          as.numeric(
            d$q90_days[idx]
          )
      )

      em[, `:=`(
        evaluation_year =
          year,
        outcome =
          id,
        sensitivity_population =
          nm
      )]

      duration_early_rows[[length(duration_early_rows) + 1L]] <- em
    }

    if (
      id ==
        "ventilator_days_conditional"
    ) {
      g <- d$airway_state_sens

      for (lv in unique(g)) {
        idx <- which(
          g ==
            lv
        )

        if (
          length(idx) <
            MIN_N
        ) {
          next
        }

        am <- duration_metrics(
          y =
            as.numeric(
              d$observed_days[idx]
            ),
          q10 =
            as.numeric(
              d$q10_days[idx]
            ),
          q50 =
            as.numeric(
              d$q50_median_days[idx]
            ),
          q90 =
            as.numeric(
              d$q90_days[idx]
            )
        )

        am[, `:=`(
          evaluation_year =
            year,
          outcome =
            id,
          airway_state =
            lv
        )]

        duration_airway_rows[[length(duration_airway_rows) + 1L]] <- am
      }
    }
  }
}

duration_overall <- safe_bind(
  duration_overall_rows
)

duration_subgroups <- safe_bind(
  duration_subgroup_rows
)

duration_early <- safe_bind(
  duration_early_rows
)

duration_airway <- safe_bind(
  duration_airway_rows
)

fwrite(
  duration_overall,
  file.path(
    OUT_DIR,
    "07_DURATION_OVERALL_TEMPORAL.csv"
  )
)

fwrite(
  duration_subgroups,
  file.path(
    OUT_DIR,
    "08_DURATION_SUBGROUP_ROBUSTNESS.csv"
  )
)

fwrite(
  duration_early,
  file.path(
    OUT_DIR,
    "09_DURATION_EARLY_DEATH_SENSITIVITY.csv"
  )
)

fwrite(
  duration_airway,
  file.path(
    OUT_DIR,
    "10_DURATION_AIRWAY_STATE_SENSITIVITY.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Locked 2024 narrow-trained models applied to repaired broad cohort
# ------------------------------------------------------------------------------

fit_encoder_apply <- function(
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
      z <- norm_chr(x)

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

union[
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

        default =
          NA_character_
      ),
      levels = c(
        "Home/home health",
        "Post-acute facility",
        "Death/hospice"
      )
    )
]

union[
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

union[
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

union[
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

broad_2024 <- union[
  admission_year ==
    2024L
]

broad_rows <- list()

broad_endpoints <- c(
  "discharge3_final",
  "hlos3_final",
  "icu3_final",
  "vent3_final",
  "operative_spine_intervention"
)

for (ep in broad_endpoints) {
  model_file <- file.path(
    CLASS_DIR,
    "models",
    paste0(
      ep,
      "__eval2024__model.rds"
    )
  )

  if (!file.exists(model_file)) {
    stop(
      "Missing locked 2024 model object: ",
      model_file,
      call. = FALSE
    )
  }

  obj <- readRDS(
    model_file
  )

  X <- fit_encoder_apply(
    broad_2024,
    obj$predictors,
    obj$encoder
  )

  if (
    !identical(
      colnames(X),
      obj$feature_names
    )
  ) {
    stop(
      "Broad sensitivity feature mismatch for ",
      ep,
      ".",
      call. = FALSE
    )
  }

  dmat <- xgb.DMatrix(
    X,
    missing = NA_real_
  )

  raw <- predict(
    obj$model,
    dmat
  )

  if (
    ep ==
      "operative_spine_intervention"
  ) {
    y <- as.integer(
      broad_2024$operative_spine_intervention
    )

    p <- as.numeric(raw)

    for (
      population in c(
        "Broad union",
        "Broad-only sensitivity"
      )
    ) {
      idx <- if (
        population ==
          "Broad union"
      ) {
        seq_len(
          nrow(
            broad_2024
          )
        )
      } else {
        which(
          broad_2024$primary_complete_trajectory !=
            TRUE
        )
      }

      if (
        length(idx) >=
          MIN_N &&
          sum(
            y[idx] ==
              1L,
            na.rm = TRUE
          ) >=
            MIN_EVENTS &&
          sum(
            y[idx] ==
              0L,
            na.rm = TRUE
          ) >=
            MIN_EVENTS
      ) {
        m <- binary_metrics(
          y[idx],
          p[idx]
        )

        m[, `:=`(
          evaluation_year =
            2024L,
          population =
            population,
          output_id =
            "recorded_tl_operation",
          output_label =
            "Recorded thoracolumbar operative treatment"
        )]

        broad_rows[[length(broad_rows) + 1L]] <- m
      }
    }
  } else {
    levels_ep <- levels(
      broad_2024[[ep]]
    )

    prob <- if (
      is.matrix(raw) ||
        is.data.frame(raw)
    ) {
      as.matrix(raw)
    } else {
      matrix(
        as.numeric(raw),
        ncol =
          length(levels_ep),
        byrow = TRUE
      )
    }

    if (
      nrow(prob) !=
        nrow(
          broad_2024
        ) ||
        ncol(prob) !=
          length(levels_ep)
    ) {
      stop(
        "Unexpected broad prediction shape for ",
        ep,
        ".",
        call. = FALSE
      )
    }

    for (
      population in c(
        "Broad union",
        "Broad-only sensitivity"
      )
    ) {
      idx <- if (
        population ==
          "Broad union"
      ) {
        seq_len(
          nrow(
            broad_2024
          )
        )
      } else {
        which(
          broad_2024$primary_complete_trajectory !=
            TRUE
        )
      }

      for (
        j in seq_along(
          levels_ep
        )
      ) {
        y <- as.integer(
          as.character(
            broad_2024[[ep]][idx]
          ) ==
            levels_ep[j]
        )

        p <- prob[
          idx,
          j
        ]

        if (
          length(idx) <
            MIN_N ||
            sum(
              y ==
                1L,
              na.rm = TRUE
            ) <
              MIN_EVENTS ||
            sum(
              y ==
                0L,
              na.rm = TRUE
            ) <
              MIN_EVENTS
        ) {
          next
        }

        m <- binary_metrics(
          y,
          p
        )

        m[, `:=`(
          evaluation_year =
            2024L,
          population =
            population,
          output_id =
            paste0(
              ep,
              "__",
              make.names(
                levels_ep[j]
              )
            ),
          output_label =
            levels_ep[j]
        )]

        broad_rows[[length(broad_rows) + 1L]] <- m
      }
    }
  }

  rm(
    obj,
    X,
    dmat,
    raw
  )
  gc()
}

broad_metrics <- safe_bind(
  broad_rows
)

fwrite(
  broad_metrics,
  file.path(
    OUT_DIR,
    "11_BROAD_PHENOTYPE_CLASSIFICATION_2024.csv"
  )
)

# ------------------------------------------------------------------------------
# 8. Temporal-era summary (fixed performance, no alternate retraining)
# ------------------------------------------------------------------------------

headline_overall <- overall[
  headline ==
    TRUE
]

era_wide <- dcast(
  headline_overall[
    ,
    .(
      output_id,
      output_label,
      evaluation_year,
      AUROC,
      AUPRC,
      Brier,
      calibration_intercept,
      calibration_slope
    )
  ],
  output_id +
    output_label ~
    evaluation_year,
  value.var = c(
    "AUROC",
    "AUPRC",
    "Brier",
    "calibration_intercept",
    "calibration_slope"
  )
)

for (
  metric in c(
    "AUROC",
    "AUPRC",
    "Brier",
    "calibration_intercept",
    "calibration_slope"
  )
) {
  c22 <- paste0(
    metric,
    "_2022"
  )
  c23 <- paste0(
    metric,
    "_2023"
  )
  c24 <- paste0(
    metric,
    "_2024"
  )

  if (
    all(
      c(
        c22,
        c23,
        c24
      ) %in%
        names(
          era_wide
        )
    )
  ) {
    era_wide[
      ,
      paste0(
        "delta_2024_minus_2022_",
        metric
      ) :=
        get(c24) -
          get(c22)
    ]

    era_wide[
      ,
      paste0(
        "delta_2024_minus_2023_",
        metric
      ) :=
        get(c24) -
          get(c23)
    ]
  }
}

fwrite(
  era_wide,
  file.path(
    OUT_DIR,
    "12_TEMPORAL_ERA_PERFORMANCE_COMPARISON.csv"
  )
)

# ------------------------------------------------------------------------------
# 9. Completion manifest
# ------------------------------------------------------------------------------

manifest <- data.table(
  sensitivity_domain = c(
    "Repaired phenotype preflight",
    "Age / sex / GCS subgroup robustness",
    "Clinical-input missingness",
    "Spine phenotype subgroups",
    "Cervical / TBI concomitant injury",
    "Race / ethnicity / payer subgroup calibration",
    "Ventilation baseline airway state",
    "Survivor / early-death exclusion",
    "Continuous-duration subgroup robustness",
    "Broad spine phenotype",
    "Temporal-era performance",
    "Transfer-in sensitivity",
    "Death-or-prolonged composite retune",
    "Facility heterogeneity"
  ),
  status = c(
    "COMPLETE",
    "COMPLETE",
    "COMPLETE",
    "COMPLETE",
    "COMPLETE",
    "COMPLETE - DESCRIPTIVE",
    "COMPLETE",
    "COMPLETE",
    "COMPLETE",
    "COMPLETE - LOCKED NARROW MODEL APPLIED TO BROAD 2024",
    "COMPLETE - FIXED TEMPORAL PERFORMANCE; NO ALTERNATE RETRAINING",
    "NOT RUN - OUTSIDE LOCKED DIRECT-PRESENTATION USE CASE",
    "NOT RUN IN NO-REFIT SUITE",
    "NOT RUN - ONLY IF USABLE FACILITY IDS AVAILABLE"
  ),
  output = c(
    "00_REPAIRED_DX_PREFLIGHT.csv",
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "05_VENTILATION_AIRWAY_STATE_SENSITIVITY.csv",
    "04_EARLY_DEATH_RESOURCE_SENSITIVITY.csv",
    "08_DURATION_SUBGROUP_ROBUSTNESS.csv",
    "11_BROAD_PHENOTYPE_CLASSIFICATION_2024.csv",
    "12_TEMPORAL_ERA_PERFORMANCE_COMPARISON.csv",
    NA_character_,
    NA_character_,
    NA_character_
  )
)

fwrite(
  manifest,
  file.path(
    OUT_DIR,
    "13_ROBUSTNESS_COMPLETION_MANIFEST.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT REPAIRED-DX FIXED-PREDICTION ROBUSTNESS COMPLETE",
    "",
    "NO primary model was tuned or refit.",
    "",
    "Review first:",
    "  06_2024_SUBGROUP_RANGE_SUMMARY.csv",
    "  04_EARLY_DEATH_RESOURCE_SENSITIVITY.csv",
    "  05_VENTILATION_AIRWAY_STATE_SENSITIVITY.csv",
    "  08_DURATION_SUBGROUP_ROBUSTNESS.csv",
    "  09_DURATION_EARLY_DEATH_SENSITIVITY.csv",
    "  11_BROAD_PHENOTYPE_CLASSIFICATION_2024.csv",
    "  12_TEMPORAL_ERA_PERFORMANCE_COMPARISON.csv",
    "  13_ROBUSTNESS_COMPLETION_MANIFEST.csv",
    "",
    "Broad phenotype sensitivity applies the LOCKED 2024 narrow-trained model to broad-cohort patients.",
    "Transfer-in sensitivity was not run because the locked use case is direct presentation.",
    "Death-or-prolonged composite model retuning is intentionally separated from this no-refit suite."
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
  "REPAIRED-DX FIXED-PREDICTION ROBUSTNESS COMPLETE\n"
)
cat(
  "============================================================\n"
)
cat(
  "Output:\n  ",
  OUT_DIR,
  "\n",
  sep = ""
)
