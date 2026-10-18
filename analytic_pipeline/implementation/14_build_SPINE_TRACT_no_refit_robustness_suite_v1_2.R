# ==============================================================================
# 14_build_SPINE_TRACT_no_refit_robustness_suite_v1_2.R
# SPINE-TRACT Stage 4A v1.2 — no-refit subgroup / robustness evaluation
#
# PURPOSE
#   Evaluate the already-locked SPINE-TRACT models across clinically relevant
#   subgroups using frozen rolling-origin prediction caches.
#
# NO MODEL IS REFIT.
#
# CLASSIFICATION OUTPUTS
#   - Home/home health
#   - Post-acute facility
#   - Death/hospice
#   - Hospital LOS >=20 days
#   - Any ICU
#   - ICU >=8 days
#   - Any invasive mechanical ventilation
#   - Ventilation >=8 days
#   - Final locked 59-predictor recorded TL operation model
#
# DURATION OUTPUTS
#   - Hospital LOS Q50 / Q10-Q90
#   - ICU LOS conditional on ICU use
#   - Ventilator duration conditional on ventilation
#
# SUBGROUPS
#   - age band
#   - sex
#   - GCS severity from component score
#   - thoracolumbar fracture distribution
#   - neural injury yes/no
#   - neural injury type
#   - concomitant cervical injury
#   - concomitant TBI
#   - baseline-input missingness quartile
#
# EVALUATION WINDOWS
#   - pooled rolling-origin OOF 2022-2024
#   - 2024 development-informed forward temporal evaluation
#
# EARLY-DEATH SENSITIVITY
#   If hospital LOS and disposition fields are available, additionally evaluates
#   all binary outputs after excluding death/hospice patients with HLOS <=1 day
#   and <=2 days. This is an evaluation sensitivity only; models are not refit.
# ==============================================================================

rm(list = ls())
gc()

required <- c("data.table")
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

suppressPackageStartupMessages(
  library(data.table)
)

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

CURRENT_DATA_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication",
  "spine_tract_frozen_primary_2020_2024.rds"
)

CURRENT_SPEC_FILE <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication",
  "spine_tract_predictor_sets.rds"
)

CLASS_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification",
  "rolling_origin_gpu",
  "predictions"
)

SURGERY_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_surgery_composite_temporal_lock",
  "predictions"
)

DURATION_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_duration_quantiles",
  "rolling_origin_R",
  "predictions"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage4_no_refit_robustness"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
)

needed <- c(
  CURRENT_DATA_FILE,
  CURRENT_SPEC_FILE,
  CLASS_DIR,
  SURGERY_DIR,
  DURATION_DIR
)

bad <- needed[
  !file.exists(needed) &
    !dir.exists(needed)
]

if (length(bad) > 0L) {
  stop(
    "Missing required artifact(s):\n",
    paste(bad, collapse = "\n"),
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 1. Current frozen patient-level subgroup data
# ------------------------------------------------------------------------------

dt <- as.data.table(
  readRDS(CURRENT_DATA_FILE)
)
spec <- readRDS(
  CURRENT_SPEC_FILE
)

dt[, inc_key := as.character(inc_key)]
dt[, admission_year := as.integer(admission_year)]

if (nrow(dt) != 424017L) {
  stop(
    "Expected frozen primary N=424017; observed ",
    nrow(dt),
    ".",
    call. = FALSE
  )
}

pg <- spec$predictor_groups
clinical_backbone <- spec$predictor_sets$clinical_backbone

coarse_spine_no_soft <- setdiff(
  pg$spine_coarse,
  "dx_tl_soft_tissue_any"
)

final_core <- unique(c(
  setdiff(
    clinical_backbone,
    pg$spine_granular
  ),
  coarse_spine_no_soft,
  pg$chest_detail,
  pg$pupil
))

if (length(final_core) != 53L) {
  stop(
    "Expected 53 final_core predictors; observed ",
    length(final_core),
    ".",
    call. = FALSE
  )
}

# Robust 0/1 flag interpretation.
# Some frozen diagnosis predictors are factor/character representations rather
# than numeric 0/1. Comparing those directly with 1L can yield NA.
is_yes_flag <- function(x) {
  if (is.logical(x)) {
    return(
      !is.na(x) &
        x
    )
  }

  if (
    is.numeric(x) ||
      is.integer(x)
  ) {
    return(
      !is.na(x) &
        as.numeric(x) == 1
    )
  }

  z <- toupper(
    trimws(
      as.character(x)
    )
  )

  !is.na(z) &
    z %in% c(
      "1",
      "TRUE",
      "T",
      "YES",
      "Y"
    )
}

# Age.
dt[, subgroup_age_band := cut(
  age_years,
  breaks = c(
    18,
    45,
    65,
    80,
    Inf
  ),
  right = FALSE,
  labels = c(
    "18-44",
    "45-64",
    "65-79",
    "80-89"
  )
)]

# Sex.
dt[, subgroup_sex := as.character(
  sex_clean
)]
dt[
  is.na(subgroup_sex) |
    trimws(subgroup_sex) == "",
  subgroup_sex := "Unknown"
]

# GCS severity from components only.
dt[, subgroup_gcs_total_components := fifelse(
  !is.na(gcs_eye_clean) &
    !is.na(gcs_verbal_clean) &
    !is.na(gcs_motor_clean),
  as.numeric(gcs_eye_clean) +
    as.numeric(gcs_verbal_clean) +
    as.numeric(gcs_motor_clean),
  NA_real_
)]

dt[, subgroup_gcs_severity := fifelse(
  is.na(
    subgroup_gcs_total_components
  ),
  "Unknown",
  fifelse(
    subgroup_gcs_total_components <= 8,
    "Severe (3-8)",
    fifelse(
      subgroup_gcs_total_components <= 12,
      "Moderate (9-12)",
      "Mild (13-15)"
    )
  )
)]

# Fracture distribution.
thoracic_fx <- is_yes_flag(
  dt$dx_thoracic_vertebral_fracture
)
lumbar_fx <- is_yes_flag(
  dt$dx_lumbar_vertebral_fracture
)

dt[, subgroup_tl_fracture_distribution := fifelse(
  thoracic_fx &
    lumbar_fx,
  "Thoracic + lumbar fracture",
  fifelse(
    thoracic_fx,
    "Thoracic fracture only",
    fifelse(
      lumbar_fx,
      "Lumbar fracture only",
      "No TL vertebral fracture"
    )
  )
)]

# Neural injury.
cord_injury <- is_yes_flag(
  dt$dx_tl_spinal_cord_injury_any
)
root_cauda_injury <- is_yes_flag(
  dt$dx_tl_root_cauda_any
)

dt[, subgroup_neural_injury := fifelse(
  cord_injury |
    root_cauda_injury,
  "Neural injury",
  "No neural injury"
)]

dt[, subgroup_neural_type := fifelse(
  cord_injury &
    root_cauda_injury,
  "Cord + root/cauda",
  fifelse(
    cord_injury,
    "Cord",
    fifelse(
      root_cauda_injury,
      "Root/cauda",
      "No neural injury"
    )
  )
)]

# Concomitant injuries.
cervical_injury <- is_yes_flag(
  dt$dx_cervical_spine_injury
)
tbi_injury <- is_yes_flag(
  dt$dx_tbi_any
)

dt[, subgroup_cervical_injury := fifelse(
  cervical_injury,
  "Concomitant cervical injury",
  "No cervical injury"
)]

dt[, subgroup_tbi := fifelse(
  tbi_injury,
  "Concomitant TBI",
  "No TBI"
)]

# Missingness burden across the locked clinical core.
is_missing_input <- function(x) {
  if (
    is.numeric(x) ||
      is.integer(x)
  ) {
    return(
      is.na(x)
    )
  }

  z <- toupper(
    trimws(
      as.character(x)
    )
  )

  is.na(z) |
    z == "" |
    z %in% c(
      "UNKNOWN",
      "__UNKNOWN__",
      "UNK",
      "MISSING",
      "NOT KNOWN"
    )
}

missing_matrix <- vapply(
  final_core,
  function(v) {
    as.integer(
      is_missing_input(
        dt[[v]]
      )
    )
  },
  integer(nrow(dt))
)

if (
  is.null(dim(missing_matrix))
) {
  missing_matrix <- matrix(
    missing_matrix,
    ncol = 1L
  )
}

dt[, subgroup_input_missing_count := rowSums(
  missing_matrix
)]

# Missingness strata.
# The frozen cohort has heavy ties in missing-count values, so nominal quartiles
# can collapse (the prior script produced Q1/Q3/Q4 with no Q2). Use transparent
# count-based strata instead.
missing_counts <- sort(
  unique(
    dt$subgroup_input_missing_count
  )
)

if (length(missing_counts) >= 3L) {
  low_cut <- missing_counts[1L]
  mid_cut <- missing_counts[2L]

  dt[, subgroup_missingness_stratum := fifelse(
    subgroup_input_missing_count <= low_cut,
    paste0(
      "Low (",
      low_cut,
      " missing)"
    ),
    fifelse(
      subgroup_input_missing_count <= mid_cut,
      paste0(
        "Intermediate (",
        mid_cut,
        " missing)"
      ),
      paste0(
        "High (>=",
        mid_cut + 1L,
        " missing)"
      )
    )
  )]
} else {
  dt[, subgroup_missingness_stratum := as.character(
    subgroup_input_missing_count
  )]
}

missing_summary <- dt[
  ,
  .(
    N = .N,
    min_missing =
      min(
        subgroup_input_missing_count
      ),
    median_missing =
      median(
        subgroup_input_missing_count
      ),
    max_missing =
      max(
        subgroup_input_missing_count
      )
  ),
  by = subgroup_missingness_stratum
]

fwrite(
  missing_summary,
  file.path(
    OUT_DIR,
    "00_MISSINGNESS_STRATUM_DEFINITION.csv"
  )
)

# Early-death proxy using the frozen discharge4 endpoint retained in Stage 2A.
dt[, early_death_1d := FALSE]
dt[, early_death_2d := FALSE]

if (
  all(
    c(
      "hospital_los_days",
      "discharge4"
    ) %in%
      names(dt)
  )
) {
  death_hospice <- as.character(
    dt$discharge4
  ) == "Death/hospice"

  dt[
    death_hospice &
      !is.na(
        hospital_los_days
      ) &
      hospital_los_days <= 1,
    early_death_1d := TRUE
  ]

  dt[
    death_hospice &
      !is.na(
        hospital_los_days
      ) &
      hospital_los_days <= 2,
    early_death_2d := TRUE
  ]
} else {
  warning(
    "Early-death sensitivity could not be constructed because hospital_los_days or discharge4 is absent."
  )
}

subgroup_vars <- c(
  subgroup_age_band =
    "Age band",
  subgroup_sex =
    "Sex",
  subgroup_gcs_severity =
    "GCS severity",
  subgroup_tl_fracture_distribution =
    "TL fracture distribution",
  subgroup_neural_injury =
    "Neural injury",
  subgroup_neural_type =
    "Neural injury type",
  subgroup_cervical_injury =
    "Concomitant cervical injury",
  subgroup_tbi =
    "Concomitant TBI",
  subgroup_missingness_stratum =
    "Baseline-input missingness stratum"
)

# Subgroup derivation QC. These groups should not be all-unknown.
subgroup_derivation_qc <- rbindlist(
  lapply(
    names(subgroup_vars),
    function(v) {
      x <- as.character(
        dt[[v]]
      )

      data.table(
        subgroup_variable = v,
        subgroup_label =
          unname(
            subgroup_vars[v]
          ),
        N = length(x),
        missing_or_unknown_n =
          sum(
            is.na(x) |
              trimws(x) == ""
          ),
        distinct_nonmissing_levels =
          uniqueN(
            x[
              !is.na(x) &
                trimws(x) != ""
            ]
          )
      )
    }
  ),
  fill = TRUE
)

fwrite(
  subgroup_derivation_qc,
  file.path(
    OUT_DIR,
    "00c_SUBGROUP_DERIVATION_QC.csv"
  )
)

if (
  any(
    subgroup_derivation_qc$
      distinct_nonmissing_levels ==
      0L
  )
) {
  stop(
    "At least one subgroup derivation produced zero nonmissing levels. See 00c_SUBGROUP_DERIVATION_QC.csv.",
    call. = FALSE
  )
}

early_death_qc <- data.table(
  total_primary_n = nrow(dt),
  early_death_1d_n =
    sum(
      dt$early_death_1d
    ),
  early_death_2d_n =
    sum(
      dt$early_death_2d
    )
)

fwrite(
  early_death_qc,
  file.path(
    OUT_DIR,
    "00d_EARLY_DEATH_DEFINITION_QC.csv"
  )
)

# ------------------------------------------------------------------------------
# 2. Binary prediction loaders
# ------------------------------------------------------------------------------

load_binary_output <- function(
    year,
    id
) {
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
        CLASS_DIR,
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
            as.character(d$inc_key),
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
            as.character(d$inc_key),
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
          as.character(d$inc_key),
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
        CLASS_DIR,
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
          as.character(d$inc_key),
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
        CLASS_DIR,
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
            as.character(d$inc_key),
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
          as.character(d$inc_key),
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
        CLASS_DIR,
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
            as.character(d$inc_key),
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
          as.character(d$inc_key),
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
        SURGERY_DIR,
        paste0(
          "G_zero_new_input_composites__eval",
          year,
          ".csv"
        )
      )
    )

    return(
      data.table(
        admission_year =
          d$admission_year,
        inc_key =
          as.character(d$inc_key),
        output_id =
          "recorded_tl_operation",
        output_label =
          "Recorded thoracolumbar operative treatment",
        y =
          as.integer(d$y),
        p =
          as.numeric(d$pred)
      )
    )
  }

  stop(
    "Unknown binary output id: ",
    id,
    call. = FALSE
  )
}

binary_ids <- c(
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
  y <- as.integer(y)
  p <- as.numeric(p)

  n1 <- as.numeric(
    sum(y == 1L)
  )
  n0 <- as.numeric(
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
    sum(r[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) /
    (
      n1 *
        n0
    )
}

fast_auprc <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)

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

binary_metrics <- function(d) {
  # Prediction caches can contain endpoint-specific unknown rows or a sensitivity
  # filter can leave non-finite values in a subgroup. Metrics must be computed
  # only on rows with known observed outcome AND finite prediction.
  y_raw <- suppressWarnings(
    as.integer(d$y)
  )
  p_raw <- suppressWarnings(
    as.numeric(d$p)
  )

  ok <- !is.na(y_raw) &
    y_raw %in% c(0L, 1L) &
    is.finite(p_raw)

  y <- y_raw[ok]
  p <- p_raw[ok]

  n <- length(y)

  if (n == 0L) {
    return(
      data.table(
        N = 0L,
        events = 0L,
        nonevents = 0L,
        prevalence = NA_real_,
        auroc = NA_real_,
        auprc = NA_real_,
        pr_lift = NA_real_,
        brier = NA_real_,
        logloss = NA_real_,
        calibration_slope = NA_real_,
        reportable = FALSE
      )
    )
  }

  pp <- clamp_prob(p)

  events <- sum(
    y == 1L
  )
  nonevents <- sum(
    y == 0L
  )

  prev <- events / n

  # AUROC requires both classes; AUPRC requires at least one event.
  auc <- if (
    events > 0L &&
      nonevents > 0L
  ) {
    fast_auc(
      y,
      p
    )
  } else {
    NA_real_
  }

  pr <- if (
    events > 0L
  ) {
    fast_auprc(
      y,
      p
    )
  } else {
    NA_real_
  }

  lp <- qlogis(pp)

  slope_fit <- if (
    events > 0L &&
      nonevents > 0L
  ) {
    tryCatch(
      suppressWarnings(
        glm(
          y ~ lp,
          family = binomial()
        )
      ),
      error = function(e) NULL
    )
  } else {
    NULL
  }

  data.table(
    N = n,
    events = events,
    nonevents = nonevents,
    prevalence = prev,
    auroc = auc,
    auprc = pr,
    pr_lift =
      if (
        is.finite(prev) &&
          prev > 0 &&
          is.finite(pr)
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
    calibration_slope =
      if (
        is.null(
          slope_fit
        )
      ) {
        NA_real_
      } else {
        unname(
          coef(
            slope_fit
          )[2L]
        )
      },
    reportable =
      n >= 500L &&
        events >= 50L &&
        nonevents >= 50L
  )
}

duration_metrics <- function(d) {
  y_raw <- suppressWarnings(
    as.numeric(d$observed_days)
  )
  q10_raw <- suppressWarnings(
    as.numeric(d$q10_days)
  )
  q50_raw <- suppressWarnings(
    as.numeric(d$q50_median_days)
  )
  q90_raw <- suppressWarnings(
    as.numeric(d$q90_days)
  )

  ok <- is.finite(y_raw) &
    is.finite(q10_raw) &
    is.finite(q50_raw) &
    is.finite(q90_raw)

  y <- y_raw[ok]
  q10 <- q10_raw[ok]
  q50 <- q50_raw[ok]
  q90 <- q90_raw[ok]

  n <- length(y)

  if (n == 0L) {
    return(
      data.table(
        N = 0L,
        median_observed = NA_real_,
        mae = NA_real_,
        median_bias = NA_real_,
        empirical_q10 = NA_real_,
        empirical_q50 = NA_real_,
        empirical_q90 = NA_real_,
        q10_q90_coverage = NA_real_,
        median_range_width = NA_real_,
        reportable = FALSE
      )
    )
  }

  data.table(
    N = n,
    median_observed =
      median(y),
    mae =
      mean(
        abs(
          y -
            q50
        )
      ),
    median_bias =
      median(
        q50 -
          y
      ),
    empirical_q10 =
      mean(
        y <=
          q10
      ),
    empirical_q50 =
      mean(
        y <=
          q50
      ),
    empirical_q90 =
      mean(
        y <=
          q90
      ),
    q10_q90_coverage =
      mean(
        y >=
          q10 &
          y <=
            q90
      ),
    median_range_width =
      median(
        q90 -
          q10
      ),
    reportable =
      n >= 200L
  )
}

# ------------------------------------------------------------------------------
# 4. Build pooled OOF binary predictions and attach subgroup data
# ------------------------------------------------------------------------------

patient_cols <- unique(c(
  "inc_key",
  "admission_year",
  names(subgroup_vars),
  "subgroup_input_missing_count",
  "early_death_1d",
  "early_death_2d"
))

patient_subgroups <- dt[
  ,
  ..patient_cols
]

binary_long <- rbindlist(
  lapply(
    binary_ids,
    function(id) {
      rbindlist(
        lapply(
          2022:2024,
          function(year) {
            load_binary_output(
              year,
              id
            )
          }
        ),
        fill = TRUE
      )
    }
  ),
  fill = TRUE
)

binary_long <- merge(
  binary_long,
  patient_subgroups,
  by = c(
    "inc_key",
    "admission_year"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (
  binary_long[
    is.na(
      subgroup_age_band
    ),
    .N
  ] > 0L
) {
  warning(
    "Some binary prediction rows did not match subgroup data."
  )
}


binary_cache_qc <- binary_long[
  ,
  .(
    rows = .N,
    unknown_y = sum(
      is.na(y) |
        !y %in% c(0L, 1L)
    ),
    nonfinite_prediction = sum(
      !is.finite(
        as.numeric(p)
      )
    )
  ),
  by = .(
    output_id,
    output_label,
    admission_year
  )
]

fwrite(
  binary_cache_qc,
  file.path(
    OUT_DIR,
    "00b_BINARY_PREDICTION_CACHE_QC.csv"
  )
)

# Overall reference metrics by output/window.
overall_rows <- list()

for (
  window in c(
    "Pooled 2022-2024",
    "2024"
  )
) {
  dwin <- if (
    window ==
      "2024"
  ) {
    binary_long[
      admission_year ==
        2024L
    ]
  } else {
    binary_long
  }

  out <- dwin[
    ,
    binary_metrics(.SD),
    by = .(
      output_id,
      output_label
    )
  ]

  out[
    ,
    evaluation_window :=
      window
  ]

  overall_rows[[length(overall_rows) + 1L]] <- out
}

overall_binary <- rbindlist(
  overall_rows,
  fill = TRUE
)

fwrite(
  overall_binary,
  file.path(
    OUT_DIR,
    "01_BINARY_OVERALL_REFERENCE_METRICS.csv"
  )
)

# Subgroup metrics.
subgroup_rows <- list()

for (
  window in c(
    "Pooled 2022-2024",
    "2024"
  )
) {
  dwin <- if (
    window ==
      "2024"
  ) {
    binary_long[
      admission_year ==
        2024L
    ]
  } else {
    binary_long
  }

  for (
    v in names(
      subgroup_vars
    )
  ) {
    label <-
      unname(
        subgroup_vars[v]
      )

    z <- dwin[
      !is.na(
        get(v)
      ),
      binary_metrics(.SD),
      by = .(
        output_id,
        output_label,
        subgroup_level =
          as.character(
            get(v)
          )
      )
    ]

    z[
      ,
      `:=`(
        evaluation_window =
          window,
        subgroup_variable =
          v,
        subgroup_label =
          label
      )
    ]

    subgroup_rows[[length(subgroup_rows) + 1L]] <- z
  }
}

binary_subgroups <- rbindlist(
  subgroup_rows,
  fill = TRUE
)

# Add change from overall reference.
binary_subgroups <- merge(
  binary_subgroups,
  overall_binary[
    ,
    .(
      output_id,
      evaluation_window,
      overall_auroc =
        auroc,
      overall_auprc =
        auprc,
      overall_brier =
        brier,
      overall_calibration_slope =
        calibration_slope
    )
  ],
  by = c(
    "output_id",
    "evaluation_window"
  ),
  all.x = TRUE,
  sort = FALSE
)

binary_subgroups[
  ,
  `:=`(
    delta_auroc_vs_overall =
      auroc -
        overall_auroc,
    delta_auprc_vs_overall =
      auprc -
        overall_auprc,
    delta_brier_vs_overall =
      brier -
        overall_brier
  )
]

fwrite(
  binary_subgroups,
  file.path(
    OUT_DIR,
    "02_BINARY_SUBGROUP_PERFORMANCE.csv"
  )
)

# A compact flag table for materially lower subgroup performance.
binary_flags <- binary_subgroups[
  reportable == TRUE &
    (
      delta_auroc_vs_overall <=
        -0.05 |
        abs(
          calibration_slope -
            1
        ) >=
        0.20
    )
][
  order(
    evaluation_window,
    output_id,
    delta_auroc_vs_overall
  )
]

fwrite(
  binary_flags,
  file.path(
    OUT_DIR,
    "03_BINARY_SUBGROUP_REVIEW_FLAGS.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. Early-death evaluation sensitivity — no refit
# ------------------------------------------------------------------------------

early_rows <- list()

for (
  window in c(
    "Pooled 2022-2024",
    "2024"
  )
) {
  dwin <- if (
    window ==
      "2024"
  ) {
    binary_long[
      admission_year ==
        2024L
    ]
  } else {
    binary_long
  }

  for (
    sens in c(
      "All patients",
      "Exclude death/hospice <=1 day",
      "Exclude death/hospice <=2 days"
    )
  ) {
    ds <- if (
      sens ==
        "Exclude death/hospice <=1 day"
    ) {
      dwin[
        early_death_1d ==
          FALSE
      ]
    } else if (
      sens ==
        "Exclude death/hospice <=2 days"
    ) {
      dwin[
        early_death_2d ==
          FALSE
      ]
    } else {
      dwin
    }

    m <- ds[
      ,
      binary_metrics(.SD),
      by = .(
        output_id,
        output_label
      )
    ]

    m[
      ,
      `:=`(
        evaluation_window =
          window,
        sensitivity =
          sens
      )
    ]

    early_rows[[length(early_rows) + 1L]] <- m
  }
}

early_death_sensitivity <- rbindlist(
  early_rows,
  fill = TRUE
)

fwrite(
  early_death_sensitivity,
  file.path(
    OUT_DIR,
    "04_EARLY_DEATH_EVALUATION_SENSITIVITY.csv"
  )
)

# ------------------------------------------------------------------------------
# 6. Duration subgroup evaluation
# ------------------------------------------------------------------------------

duration_ids <- c(
  "hospital_los",
  "icu_los_conditional",
  "ventilator_days_conditional"
)

duration_long <- rbindlist(
  lapply(
    duration_ids,
    function(id) {
      rbindlist(
        lapply(
          2022:2024,
          function(year) {
            f <- file.path(
              DURATION_DIR,
              paste0(
                id,
                "__eval",
                year,
                "__Q10_Q50_Q90_predictions.csv"
              )
            )

            d <- fread(f)

            d[
              ,
              `:=`(
                inc_key =
                  as.character(
                    inc_key
                  ),
                outcome_id =
                  id
              )
            ]

            d
          }
        ),
        fill = TRUE
      )
    }
  ),
  fill = TRUE
)

duration_long <- merge(
  duration_long,
  patient_subgroups,
  by = c(
    "inc_key",
    "admission_year"
  ),
  all.x = TRUE,
  sort = FALSE
)

duration_overall_rows <- list()

for (
  window in c(
    "Pooled 2022-2024",
    "2024"
  )
) {
  dwin <- if (
    window ==
      "2024"
  ) {
    duration_long[
      admission_year ==
        2024L
    ]
  } else {
    duration_long
  }

  m <- dwin[
    ,
    duration_metrics(.SD),
    by = outcome_id
  ]

  m[
    ,
    evaluation_window :=
      window
  ]

  duration_overall_rows[[length(duration_overall_rows) + 1L]] <- m
}

duration_overall <- rbindlist(
  duration_overall_rows,
  fill = TRUE
)

fwrite(
  duration_overall,
  file.path(
    OUT_DIR,
    "05_DURATION_OVERALL_REFERENCE_METRICS.csv"
  )
)

duration_subgroup_rows <- list()

for (
  window in c(
    "Pooled 2022-2024",
    "2024"
  )
) {
  dwin <- if (
    window ==
      "2024"
  ) {
    duration_long[
      admission_year ==
        2024L
    ]
  } else {
    duration_long
  }

  for (
    v in names(
      subgroup_vars
    )
  ) {
    label <-
      unname(
        subgroup_vars[v]
      )

    z <- dwin[
      !is.na(
        get(v)
      ),
      duration_metrics(.SD),
      by = .(
        outcome_id,
        subgroup_level =
          as.character(
            get(v)
          )
      )
    ]

    z[
      ,
      `:=`(
        evaluation_window =
          window,
        subgroup_variable =
          v,
        subgroup_label =
          label
      )
    ]

    duration_subgroup_rows[[length(duration_subgroup_rows) + 1L]] <- z
  }
}

duration_subgroups <- rbindlist(
  duration_subgroup_rows,
  fill = TRUE
)

duration_subgroups <- merge(
  duration_subgroups,
  duration_overall[
    ,
    .(
      outcome_id,
      evaluation_window,
      overall_mae =
        mae,
      overall_coverage =
        q10_q90_coverage,
      overall_median_bias =
        median_bias
    )
  ],
  by = c(
    "outcome_id",
    "evaluation_window"
  ),
  all.x = TRUE,
  sort = FALSE
)

duration_subgroups[
  ,
  `:=`(
    delta_mae_vs_overall =
      mae -
        overall_mae,
    delta_coverage_vs_overall =
      q10_q90_coverage -
        overall_coverage
  )
]

fwrite(
  duration_subgroups,
  file.path(
    OUT_DIR,
    "06_DURATION_SUBGROUP_PERFORMANCE.csv"
  )
)

duration_flags <- duration_subgroups[
  reportable == TRUE &
    (
      delta_mae_vs_overall >=
        2 |
        abs(
          delta_coverage_vs_overall
        ) >=
        0.10
    )
][
  order(
    evaluation_window,
    outcome_id,
    -delta_mae_vs_overall
  )
]

fwrite(
  duration_flags,
  file.path(
    OUT_DIR,
    "07_DURATION_SUBGROUP_REVIEW_FLAGS.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Subgroup prevalence / composition table
# ------------------------------------------------------------------------------

composition_rows <- list()

for (
  v in names(
    subgroup_vars
  )
) {
  z <- dt[
    ,
    .N,
    by = .(
      admission_year,
      subgroup_level =
        as.character(
          get(v)
        )
    )
  ]

  z[
    ,
    `:=`(
      subgroup_variable =
        v,
      subgroup_label =
        unname(
          subgroup_vars[v]
        )
    )
  ]

  composition_rows[[length(composition_rows) + 1L]] <- z
}

composition <- rbindlist(
  composition_rows,
  fill = TRUE
)

fwrite(
  composition,
  file.path(
    OUT_DIR,
    "08_SUBGROUP_COMPOSITION_BY_YEAR.csv"
  )
)

manifest <- data.table(
  item = c(
    "Model refitting",
    "Evaluation windows",
    "Surgery model",
    "Classification subgroup reportability",
    "Duration subgroup reportability",
    "Early-death analysis"
  ),
  specification = c(
    "None",
    "Pooled rolling-origin OOF 2022-2024 and 2024",
    "Final locked 59-predictor G_zero_new_input_composites model",
    "N>=500, >=50 events, >=50 nonevents",
    "N>=200",
    "Evaluation-only sensitivity excluding death/hospice with HLOS <=1 or <=2 days when available"
  )
)

fwrite(
  manifest,
  file.path(
    OUT_DIR,
    "09_STAGE4A_MANIFEST.csv"
  )
)

cat("\n============================================================\n")
cat("STAGE 4A v1.2 NO-REFIT ROBUSTNESS COMPLETE\n")
cat("============================================================\n")
cat(
  "Output directory:\n  ",
  OUT_DIR,
  "\n\n",
  sep = ""
)
cat(
  "Binary subgroup review flags: ",
  nrow(binary_flags),
  "\n",
  sep = ""
)
cat(
  "Duration subgroup review flags: ",
  nrow(duration_flags),
  "\n",
  sep = ""
)

writeLines(
  c(
    "SPINE-TRACT STAGE 4A NO-REFIT ROBUSTNESS",
    "",
    "No model was refit.",
    "",
    "Review first:",
    "  03_BINARY_SUBGROUP_REVIEW_FLAGS.csv",
    "  02_BINARY_SUBGROUP_PERFORMANCE.csv",
    "  04_EARLY_DEATH_EVALUATION_SENSITIVITY.csv",
    "  07_DURATION_SUBGROUP_REVIEW_FLAGS.csv",
    "  06_DURATION_SUBGROUP_PERFORMANCE.csv",
    "",
    "This stage evaluates subgroup transportability/stability of frozen predictions.",
    "It does not replace prespecified cohort/predictor sensitivity refits such as broad phenotype, airway-state augmentation, or social-context augmentation."
  ),
  file.path(
    OUT_DIR,
    "README_RESULTS.txt"
  )
)
