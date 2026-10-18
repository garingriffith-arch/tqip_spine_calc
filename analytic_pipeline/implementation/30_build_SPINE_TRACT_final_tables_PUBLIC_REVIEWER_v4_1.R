# ==============================================================================
# 30_build_SPINE_TRACT_final_tables_PUBLIC_REVIEWER_v4_1.R
#
# SPINE-TRACT PUBLIC / REVIEWER-FACING TABLE PIPELINE
#
# MAIN
#   Table 1  Cohort characteristics + temporal case-mix comparison
#   Table 2  Final 2024 probability / trajectory performance
#   Table 3  Final 2024 continuous-duration performance
#
# SUPPLEMENT
#   eTable 1  Cohort construction, observation status, endpoint definitions
#   eTable 2  Final predictor dictionary + missingness
#   eTable 3  Predictor adjudication and final architecture decisions
#   eTable 4  Complete 2022/2023/2024 classification performance
#   eTable 5  Temporal duration performance + threshold justification
#   eTable 6  Paper-facing model specification and reproducibility
#   eTable 7  XGBoost-vs-ridge and predictor-family sensitivity
#   eTable 8  Subgroup / airway / early-death / broad-phenotype robustness
#   eTable 9  Recorded thoracolumbar operative-treatment endpoint audit
#
# Reviewer-facing principles:
#   - locked repaired-DX artifacts only
#   - no large-sample P-value column in Table 1
#   - calibration / proper scoring reported beside discrimination
#   - surgery = recorded treatment receipt, not operative indication
#   - transfer-in exclusion is transparent in cohort construction
#   - machine-readable CSV source data are emitted with the DOCX files
#
# This script DOES NOT fit, tune, or recalibrate any model.
# ==============================================================================

rm(list = ls())
gc()

required_packages <- c("data.table", "officer", "flextable", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(officer)
  library(flextable)
  library(scales)
})

# ------------------------------------------------------------------------------
# 1. Paths and artifact discovery
# ------------------------------------------------------------------------------

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

if (!dir.exists(PROJECT_DIR)) {
  stop("Project directory does not exist: ", PROJECT_DIR, call. = FALSE)
}

OUT_DIR <- file.path(PROJECT_DIR, "outputs", "SPINE_TRACT_PUBLICATION_TABLES")
SOURCE_DIR <- file.path(OUT_DIR, "source_data")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(SOURCE_DIR, recursive = TRUE, showWarnings = FALSE)

search_roots <- c(
  file.path(PROJECT_DIR, "outputs"),
  file.path(PROJECT_DIR, "data", "derived"),
  file.path(PROJECT_DIR, "R"),
  file.path(PROJECT_DIR, "reference")
)

all_files <- unique(c(
  unlist(
    lapply(
      search_roots[dir.exists(search_roots)],
      function(root) {
        list.files(
          root,
          recursive = TRUE,
          full.names = TRUE,
          include.dirs = FALSE
        )
      }
    ),
    use.names = FALSE
  ),
  list.files(
    PROJECT_DIR,
    recursive = FALSE,
    full.names = TRUE,
    include.dirs = FALSE
  )
))

find_file <- function(filename, prefer = NULL, required = TRUE) {
  hits <- all_files[basename(all_files) == filename]
  if (!is.null(prefer) && length(hits) > 0L) {
    pref <- hits[grepl(prefer, hits, fixed = TRUE)]
    if (length(pref) > 0L) hits <- pref
  }
  if (length(hits) == 0L) {
    if (isTRUE(required)) {
      stop(
        "Could not find required artifact: ", filename,
        if (!is.null(prefer)) paste0("\nPreferred path contained: ", prefer) else "",
        call. = FALSE
      )
    }
    return(NA_character_)
  }
  hits[which.min(nchar(hits))]
}

files <- list(
  primary = find_file(
    "spine_tract_frozen_primary_2020_2024.rds",
    "stage2_predictor_adjudication"
  ),
  predictor_spec = find_file(
    "spine_tract_predictor_sets.rds",
    "stage2_predictor_adjudication"
  ),
  cohort_flow = find_file("06_cohort_flow_by_year.csv", "gate1_audit_v2_1"),
  operative_categories = find_file("20_operative_category_counts_by_year.csv", "gate1_audit_v2_1"),
  operative_rates = find_file("21_operative_event_rate_by_year.csv", "gate1_audit_v2_1"),
  operative_timing = find_file("22_procedure_timing_completeness.csv", "gate1_audit_v2_1"),
  operative_map = find_file("spine_operative_icd10pcs_map.csv", NULL, required = FALSE),
  matched_procedures = find_file(
    "spine_gate1_v2_matched_operative_procedures_2020_2024.rds",
    "gate1_audit_v2_1",
    required = FALSE
  ),

  classification_temporal = find_file(
    "01_CLASSIFICATION_TEMPORAL_POINT_ESTIMATES.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),
  classification_2024_ci = find_file(
    "04_CLASSIFICATION_2024_BOOTSTRAP_CI.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),
  duration_temporal = find_file(
    "07_DURATION_TEMPORAL_POINT_ESTIMATES.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),
  duration_2024_ci = find_file(
    "09_DURATION_2024_BOOTSTRAP_CI.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),

  ridge = find_file(
    "06_FINAL_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv",
    "stage3_ridge_comparator_REPAIRED_DX_FINAL",
    required = FALSE
  ),
  stage2b = find_file(
    "09_predictor_family_delta_summary.csv",
    "gpu_ablation_REPAIRED_DX_v1",
    required = FALSE
  ),
  stage2c = find_file(
    "09_stage2C_finalist_delta_summary.csv",
    "stage2C_REPAIRED_DX_v1",
    required = FALSE
  ),
  class_configs = find_file(
    "14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv",
    "stage3_final_classification_REPAIRED_DX",
    required = FALSE
  ),
  duration_configs = find_file(
    "04_SELECTED_DURATION_CONFIGS.csv",
    "stage3_duration_quantiles_REPAIRED_DX",
    required = FALSE
  ),

  subgroup = find_file(
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "stage4_robustness_REPAIRED_DX",
    required = FALSE
  ),
  early_death = find_file(
    "04_EARLY_DEATH_RESOURCE_SENSITIVITY.csv",
    "stage4_robustness_REPAIRED_DX",
    required = FALSE
  ),
  airway = find_file(
    "05_VENTILATION_AIRWAY_STATE_SENSITIVITY.csv",
    "stage4_robustness_REPAIRED_DX",
    required = FALSE
  ),
  broad = find_file(
    "11_BROAD_PHENOTYPE_CLASSIFICATION_2024.csv",
    "stage4_robustness_REPAIRED_DX",
    required = FALSE
  ),
  composite_2024 = find_file(
    "08_2024_COMPOSITE_SENSITIVITY.csv",
    "stage4_death_or_prolonged_composite_REPAIRED_DX",
    required = FALSE
  )
)

artifact_manifest <- data.table(
  artifact = names(files),
  path = unlist(files, use.names = FALSE)
)
artifact_manifest[, exists := file.exists(path)]
fwrite(artifact_manifest, file.path(OUT_DIR, "TABLE_SOURCE_ARTIFACT_MANIFEST.csv"))

# ------------------------------------------------------------------------------
# 2. Locked data and final predictor policy
# ------------------------------------------------------------------------------

dt <- as.data.table(readRDS(files$primary))
spec <- readRDS(files$predictor_spec)
dt[, admission_year := as.integer(admission_year)]

if (nrow(dt) != 424017L) {
  stop("Expected repaired primary cohort N=424,017; observed ", nrow(dt), ".", call. = FALSE)
}

pg <- spec$predictor_groups
clinical_backbone <- spec$predictor_sets$clinical_backbone

final_core <- unique(c(
  setdiff(clinical_backbone, pg$spine_granular),
  pg$spine_coarse,
  pg$chest_detail,
  pg$pupil
))
final_context <- unique(c(final_core, pg$social_context))

if (
  length(final_core) != 54L ||
  length(final_context) != 57L ||
  !"dx_tl_soft_tissue_any" %in% final_core
) {
  stop("Corrected 54/57 predictor lock could not be reconstructed.", call. = FALSE)
}

classification_temporal <- fread(files$classification_temporal)
classification_2024_ci <- fread(files$classification_2024_ci)
duration_temporal <- fread(files$duration_temporal)
duration_2024_ci <- fread(files$duration_2024_ci)
cohort_flow <- fread(files$cohort_flow)
operative_categories <- fread(files$operative_categories)
operative_rates <- fread(files$operative_rates)
operative_timing <- fread(files$operative_timing)

# ------------------------------------------------------------------------------
# 3. Publication formatting helpers
# ------------------------------------------------------------------------------

NAVY <- "#17365D"
VERY_LIGHT_GRAY <- "#F2F2F2"
MID_GRAY <- "#7F7F7F"

fmt_int <- function(x) scales::comma(round(x), accuracy = 1)
fmt_pct <- function(x) ifelse(is.na(x), "", scales::percent(x, accuracy = 0.1))
fmt3 <- function(x) ifelse(is.na(x), "", sprintf("%.3f", x))
fmt_signed3 <- function(x) ifelse(is.na(x), "", sprintf("%+.3f", x))
fmt_ci3 <- function(est, lo, hi) {
  ifelse(is.na(est), "", sprintf("%.3f (%.3f-%.3f)", est, lo, hi))
}
fmt_days_ci2 <- function(est, lo, hi) {
  ifelse(is.na(est), "", sprintf("%.2f (%.2f-%.2f)", est, lo, hi))
}

safe_char <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | z == ""] <- "Unknown"
  z
}

first_existing <- function(candidates) {
  hit <- candidates[candidates %in% names(dt)]
  if (length(hit) == 0L) return(NA_character_)
  hit[1L]
}

binary_smd <- function(p_dev, p_2024) {
  denom <- sqrt((p_dev * (1 - p_dev) + p_2024 * (1 - p_2024)) / 2)
  if (!is.finite(denom) || denom == 0) return(0)
  abs((p_2024 - p_dev) / denom)
}

continuous_smd <- function(x_dev, x_2024) {
  x_dev <- as.numeric(x_dev)
  x_2024 <- as.numeric(x_2024)
  sd_pool <- sqrt((var(x_dev, na.rm = TRUE) + var(x_2024, na.rm = TRUE)) / 2)
  if (!is.finite(sd_pool) || sd_pool == 0) return(0)
  abs((mean(x_2024, na.rm = TRUE) - mean(x_dev, na.rm = TRUE)) / sd_pool)
}

pretty_name <- function(x) {
  known <- c(
    "age_years" = "Age",
    "sex_clean" = "Sex",
    "sex" = "Sex",
    "gcs_eye_clean" = "GCS eye",
    "gcs_verbal_clean" = "GCS verbal",
    "gcs_motor_clean" = "GCS motor",
    "gcsq_intubated_recovered" = "GCS qualifier: intubated",
    "gcsq_sedated_paralyzed_recovered" = "GCS qualifier: sedated/paralyzed",
    "gcsq_eye_obstruction_recovered" = "GCS qualifier: eye obstruction",
    "gcsq_unknown_recovered" = "GCS qualifier: unknown/unassessable",
    "sbp_clean" = "Systolic blood pressure",
    "hr_clean" = "Pulse",
    "pulse_clean" = "Pulse",
    "rr_clean" = "Respiratory rate",
    "spo2_clean" = "Oxygen saturation",
    "temperature_c_clean" = "Temperature",
    "temperature_c_recovered" = "Temperature",
    "supplemental_oxygen_clean" = "Supplemental oxygen",
    "supplemental_oxygen_recovered" = "Supplemental oxygen",
    "mechanism_clean" = "Mechanism",
    "primary_ecode" = "Mechanism",
    "pupil_clean" = "Pupillary response",
    "dx_thoracic_vertebral_fracture" = "Thoracic vertebral fracture",
    "dx_lumbar_vertebral_fracture" = "Lumbar vertebral fracture",
    "dx_tl_spinal_cord_injury_any" = "Thoracolumbar spinal cord injury",
    "dx_tl_root_cauda_any" = "Thoracolumbar root/cauda injury",
    "dx_tl_soft_tissue_any" = "Thoracolumbar ligament/dislocation/sprain",
    "dx_cervical_spine_injury" = "Concomitant cervical spine injury",
    "dx_tbi_any" = "Concomitant traumatic brain injury",
    "race_clean" = "Race",
    "ethnicity_clean" = "Ethnicity",
    "payer_clean" = "Payer",
    "pmhx_anticoagulant_therapy" = "Anticoagulant therapy",
    "pmhx_bleeding_disorder" = "Bleeding disorder",
    "pmhx_diabetes" = "Diabetes",
    "pmhx_hypertension" = "Hypertension",
    "pmhx_current_smoker" = "Current smoking",
    "pmhx_functional_dependence" = "Functional dependence",
    "pmhx_copd" = "COPD",
    "pmhx_chf" = "Congestive heart failure"
  )

  out <- unname(
    known[x]
  )

  miss <- is.na(out)

  if (any(miss)) {
    z <- x[miss]
    z <- gsub("^pmhx_", "", z)
    z <- gsub("^dx_", "", z)
    z <- gsub("_any$", "", z)
    z <- gsub("_clean$", "", z)
    z <- gsub("_recovered$", "", z)
    z <- gsub("_", " ", z, fixed = TRUE)

    # Public-facing predictor labels use sentence case rather than carrying
    # analysis-variable title casing into the supplement.
    z <- tolower(z)
    z <- paste0(
      toupper(
        substr(
          z,
          1L,
          1L
        )
      ),
      substr(
        z,
        2L,
        nchar(z)
      )
    )

    z <- gsub("\\bCopd\\b", "COPD", z)
    z <- gsub("\\bChf\\b", "Congestive heart failure", z)
    z <- gsub("\\bTbi\\b", "TBI", z)
    z <- gsub("\\bTl\\b", "Thoracolumbar", z)
    z <- gsub("\\bGcs\\b", "GCS", z)
    z <- gsub("\\bIcu\\b", "ICU", z)
    z <- gsub("\\bLos\\b", "LOS", z)

    out[miss] <- z
  }

  out
}

pretty_symbols <- function(x) {
  z <- as.character(x)
  z <- gsub(">=", "≥", z, fixed = TRUE)
  z <- gsub("<=", "≤", z, fixed = TRUE)
  z <- gsub("Q10-Q90", "Q10–Q90", z, fixed = TRUE)
  z <- gsub("2022-2024", "2022–2024", z, fixed = TRUE)
  z <- gsub("2020-2023", "2020–2023", z, fixed = TRUE)
  z
}

humanize_model_endpoint <- function(x) {
  map <- c(
    "discharge3_final" = "Disposition trajectory",
    "discharge4" = "Disposition",
    "hlos3_final" = "Hospital LOS trajectory",
    "hlos3" = "Hospital LOS",
    "icu3_final" = "ICU trajectory",
    "icu3" = "ICU",
    "vent3_final" = "Ventilation trajectory",
    "vent3" = "Ventilation",
    "operative_spine_intervention" = "Recorded thoracolumbar operation",
    "hospital_los" = "Hospital LOS quantiles",
    "icu_los_conditional" = "ICU LOS quantiles among ICU users",
    "ventilator_days_conditional" = "Ventilator-duration quantiles among ventilated patients"
  )

  out <- unname(
    map[
      as.character(x)
    ]
  )

  out[
    is.na(out)
  ] <- pretty_symbols(
    gsub(
      "_",
      " ",
      as.character(x)[is.na(out)],
      fixed = TRUE
    )
  )

  out
}

classification_domain <- function(output_id) {
  fcase(
    grepl("^disposition", output_id),
    "Disposition",

    grepl("^hlos", output_id),
    "Hospital LOS",

    grepl("^icu", output_id),
    "ICU",

    grepl("^vent", output_id),
    "Ventilation",

    output_id ==
      "recorded_tl_operation",
    "Spine treatment",

    default =
      "Other"
  )
}

humanize_predictor_policy <- function(x) {
  z <- as.character(x)

  fcase(
    z %in%
      c(
        "final_core",
        "clinical_core",
        "Pragmatic clinical"
      ),
    "Clinical core (54 predictors)",

    z %in%
      c(
        "final_context",
        "Pragmatic clinical + race + ethnicity + payer"
      ),
    "Contextual set (57 predictors)",

    default =
      pretty_symbols(z)
  )
}

col_or_na <- function(
    d,
    candidates,
    default = NA
) {
  hit <- candidates[
    candidates %in%
      names(d)
  ]

  if (length(hit) == 0L) {
    return(
      rep(
        default,
        nrow(d)
      )
    )
  }

  d[[hit[1L]]]
}


humanize_subgroup_level <- function(
    domain,
    level
) {
  z <- as.character(level)

  z[
    z %in%
      c(
        "__UNKNOWN__",
        "Unknown",
        "UNKNOWN"
      )
  ] <- "Unknown / not recorded"

  if (domain == "Sex") {
    z[z == "SEX_1"] <- "Male"
    z[z == "SEX_2"] <- "Female"
    z[z == "SEX_3"] <- "Non-binary"
  }

  z
}

predictor_units <- function(v) {
  known <- c(
    "age_years" = "years",
    "sbp_clean" = "mm Hg",
    "hr_clean" = "beats/min",
    "pulse_clean" = "beats/min",
    "rr_clean" = "breaths/min",
    "spo2_clean" = "%",
    "temperature_c_clean" = "degrees C",
    "temperature_c_recovered" = "degrees C",
    "gcs_eye_clean" = "points",
    "gcs_verbal_clean" = "points",
    "gcs_motor_clean" = "points"
  )

  out <- unname(known[v])

  if (is.na(out)) {
    if (predictor_type(v) == "Binary") return("0/1")
    if (predictor_type(v) == "Categorical") return("categorical")
    return("registry-defined units")
  }

  out
}

implementation_levels <- function(v) {
  x <- dt[[v]]
  typ <- predictor_type(v)

  if (typ == "Binary") {
    return("0 = absent/no; 1 = present/yes")
  }

  if (typ == "Categorical") {
    z <- sort(unique(safe_char(x)))

    if (length(z) > 25L) {
      return(
        paste0(
          length(z),
          " observed categories; see source coding"
        )
      )
    }

    return(
      paste(
        z,
        collapse = " | "
      )
    )
  }

  xnum <- suppressWarnings(as.numeric(x))
  xr <- range(xnum[is.finite(xnum)], na.rm = TRUE)

  if (length(xr) == 2L && all(is.finite(xr))) {
    return(
      paste0(
        "Observed range: ",
        signif(xr[1L], 6),
        " to ",
        signif(xr[2L], 6)
      )
    )
  }

  "Numeric"
}

predictor_domain <- function(v) {
  if (grepl("^age|^sex", v)) return("Demographics")
  if (grepl("^gcs|pupil", v)) return("Neurologic examination")
  if (grepl("sbp|hr_|pulse|rr_|spo2|temperature|oxygen", v)) return("Admission physiology")
  if (grepl("mechanism|ecode", v)) return("Injury mechanism")
  if (grepl("^dx_thoracic_vertebral|^dx_lumbar_vertebral|^dx_tl_|root|cauda", v)) {
    return("Thoracolumbar phenotype")
  }
  if (grepl("^dx_", v)) return("Associated injury")
  if (grepl("^pmhx_", v)) return("Pre-existing conditions")
  if (v %in% pg$social_context) return("Social / health-system context")
  "Other"
}

predictor_source <- function(v) {
  dom <- predictor_domain(v)
  if (dom == "Thoracolumbar phenotype") return("ICD-10-CM diagnosis-derived phenotype")
  if (dom == "Associated injury") return("ICD-10-CM concomitant injury phenotype")
  if (dom == "Pre-existing conditions") return("ACS TQP pre-existing-condition field")
  if (dom == "Neurologic examination") return("Admission neurologic examination / qualifier field")
  if (dom == "Admission physiology") return("Admission vital-sign / oxygen-support field")
  if (dom == "Injury mechanism") return("Registry injury-mechanism field")
  if (dom == "Demographics") return("Registry demographic field")
  if (dom == "Social / health-system context") return("Registry demographic / payer field")
  "Derived registry field"
}

predictor_type <- function(v) {
  x <- dt[[v]]
  if (is.numeric(x) || is.integer(x)) {
    vals <- unique(x[!is.na(x)])
    if (length(vals) <= 2L && all(vals %in% c(0, 1))) return("Binary")
    return("Numeric")
  }
  "Categorical"
}

missing_count <- function(x) {
  if (is.numeric(x) || is.integer(x)) {
    return(sum(is.na(x) | !is.finite(as.numeric(x))))
  }
  z <- trimws(as.character(x))
  sum(
    is.na(z) |
    z == "" |
    toupper(z) %in% c("UNKNOWN", "__UNKNOWN__", "UNK")
  )
}

sectionize_display <- function(
    dat,
    group_col
) {
  if (
    is.null(group_col) ||
      !group_col %in%
        names(dat)
  ) {
    return(
      list(
        data =
          copy(dat),
        section_rows =
          integer()
      )
    )
  }

  groups <- unique(
    as.character(
      dat[[group_col]]
    )
  )

  visible_cols <- setdiff(
    names(dat),
    group_col
  )

  pieces <- vector(
    "list",
    length(groups)
  )

  for (i in seq_along(groups)) {
    g <- groups[i]

    body <- copy(
      dat[
        as.character(
          get(group_col)
        ) ==
          g,
        ..visible_cols
      ]
    )

    for (v in names(body)) {
      body[
        ,
        (v) :=
          as.character(
            get(v)
          )
      ]
    }

    section <- as.data.table(
      as.list(
        setNames(
          rep(
            "",
            length(
              visible_cols
            )
          ),
          visible_cols
        )
      )
    )

    section[
      1L,
      (visible_cols[1L]) :=
        g
    ]

    section[
      ,
      `.section_row` :=
        TRUE
    ]

    body[
      ,
      `.section_row` :=
        FALSE
    ]

    pieces[[i]] <- rbindlist(
      list(
        section,
        body
      ),
      use.names =
        TRUE,
      fill =
        TRUE
    )
  }

  out <- rbindlist(
    pieces,
    use.names =
      TRUE,
    fill =
      TRUE
  )

  rows <- which(
    out$.section_row
  )

  out[
    ,
    `.section_row` :=
      NULL
  ]

  list(
    data =
      out,
    section_rows =
      rows
  )
}

make_pub_ft <- function(
    dat,
    header_labels = NULL,
    note = NULL,
    section_col = NULL,
    font_size = 10,
    header_size = 10
) {
  sec <- sectionize_display(
    dat,
    section_col
  )

  display_dat <- sec$data
  section_rows <- sec$section_rows

  ft <- flextable(
    display_dat
  )

  if (!is.null(header_labels)) {
    ft <- set_header_labels(
      ft,
      values =
        header_labels
    )
  }

  ft <- theme_booktabs(
    ft
  )

  ft <- bg(
    ft,
    part =
      "header",
    bg =
      NAVY
  )

  ft <- color(
    ft,
    part =
      "header",
    color =
      "white"
  )

  ft <- bold(
    ft,
    part =
      "header"
  )

  ft <- font(
    ft,
    fontname =
      "Arial",
    part =
      "all"
  )

  ft <- fontsize(
    ft,
    size =
      font_size,
    part =
      "body"
  )

  ft <- fontsize(
    ft,
    size =
      header_size,
    part =
      "header"
  )

  ft <- align(
    ft,
    align =
      "center",
    part =
      "all"
  )

  if (ncol(display_dat) >= 1L) {
    ft <- align(
      ft,
      j =
        1,
      align =
        "left",
      part =
        "body"
    )
  }

  ft <- valign(
    ft,
    valign =
      "center",
    part =
      "all"
  )

  ft <- padding(
    ft,
    padding =
      3,
    part =
      "all"
  )

  if (length(section_rows) > 0L) {
    for (ii in section_rows) {
      ft <- merge_at(
        ft,
        i =
          ii,
        j =
          seq_len(
            ncol(
              display_dat
            )
          ),
        part =
          "body"
      )

      ft <- bg(
        ft,
        i =
          ii,
        bg =
          VERY_LIGHT_GRAY,
        part =
          "body"
      )

      ft <- color(
        ft,
        i =
          ii,
        color =
          NAVY,
        part =
          "body"
      )

      ft <- bold(
        ft,
        i =
          ii,
        part =
          "body"
      )

      ft <- align(
        ft,
        i =
          ii,
        align =
          "left",
        part =
          "body"
      )

      ft <- padding(
        ft,
        i =
          ii,
        padding.top =
          4,
        padding.bottom =
          4,
        part =
          "body"
      )
    }
  }

  ft <- set_table_properties(
    ft,
    layout =
      "autofit",
    width =
      1
  )

  if (!is.null(note)) {
    ft <- add_footer_lines(
      ft,
      values =
        note
    )

    ft <- fontsize(
      ft,
      size =
        9,
      part =
        "footer"
    )

    ft <- italic(
      ft,
      part =
        "footer"
    )

    ft <- color(
      ft,
      color =
        MID_GRAY,
      part =
        "footer"
    )

    ft <- align(
      ft,
      align =
        "left",
      part =
        "footer"
    )
  }

  ft
}

table_title_fpar <- function(
    prefix,
    title
) {
  fpar(
    ftext(
      prefix,
      prop =
        fp_text(
          font.family =
            "Times New Roman",
          font.size =
            12,
          bold =
            TRUE
        )
    ),
    ftext(
      paste0(
        " ",
        title
      ),
      prop =
        fp_text(
          font.family =
            "Times New Roman",
          font.size =
            12,
          bold =
            FALSE
        )
    ),
    fp_p =
      fp_par(
        padding.bottom =
          5,
        keep_with_next =
          TRUE
      )
  )
}

add_table_to_doc <- function(
    doc,
    prefix,
    title,
    ft,
    page_break_before = FALSE
) {
  if (isTRUE(page_break_before)) {
    doc <- body_add_break(
      doc
    )
  }

  doc <- body_add_fpar(
    doc,
    value =
      table_title_fpar(
        prefix,
        title
      )
  )

  doc <- body_add_flextable(
    doc,
    value =
      ft,
    align =
      "center"
  )

  doc
}

landscape_section <- prop_section(
  page_size = page_size(orient = "landscape", width = 11.7, height = 8.3),
  page_margins = page_mar(top = 0.45, bottom = 0.45, left = 0.45, right = 0.45),
  type = "continuous"
)

# ------------------------------------------------------------------------------
# 4. MAIN TABLE 1 — Cohort characteristics + temporal case-mix comparison
# ------------------------------------------------------------------------------

dev <- dt[
  admission_year %in%
    2020:2023
]

eval24 <- dt[
  admission_year ==
    2024L
]

allcohort <- dt

table1_rows <- list()

add_continuous_row <- function(
    domain,
    label,
    var
) {
  x_all <- as.numeric(
    allcohort[[var]]
  )

  x_dev <- as.numeric(
    dev[[var]]
  )

  x_24 <- as.numeric(
    eval24[[var]]
  )

  table1_rows[[length(table1_rows) + 1L]] <<- data.table(
    Domain =
      domain,
    Characteristic =
      label,
    Overall =
      sprintf(
        "%.1f (%.1f)",
        mean(
          x_all,
          na.rm = TRUE
        ),
        sd(
          x_all,
          na.rm = TRUE
        )
      ),
    `Development 2020–2023` =
      sprintf(
        "%.1f (%.1f)",
        mean(
          x_dev,
          na.rm = TRUE
        ),
        sd(
          x_dev,
          na.rm = TRUE
        )
      ),
    `2024 temporal cohort` =
      sprintf(
        "%.1f (%.1f)",
        mean(
          x_24,
          na.rm = TRUE
        ),
        sd(
          x_24,
          na.rm = TRUE
        )
      ),
    `Absolute SMD` =
      sprintf(
        "%.3f",
        continuous_smd(
          x_dev,
          x_24
        )
      )
  )
}

add_binary_row <- function(
    domain,
    label,
    indicator_all,
    indicator_dev,
    indicator_24
) {
  p_all <- mean(
    indicator_all,
    na.rm = TRUE
  )

  p_dev <- mean(
    indicator_dev,
    na.rm = TRUE
  )

  p_24 <- mean(
    indicator_24,
    na.rm = TRUE
  )

  table1_rows[[length(table1_rows) + 1L]] <<- data.table(
    Domain =
      domain,
    Characteristic =
      label,
    Overall =
      fmt_pct(
        p_all
      ),
    `Development 2020–2023` =
      fmt_pct(
        p_dev
      ),
    `2024 temporal cohort` =
      fmt_pct(
        p_24
      ),
    `Absolute SMD` =
      sprintf(
        "%.3f",
        binary_smd(
          p_dev,
          p_24
        )
      )
  )
}

age_var <- first_existing(
  c(
    "age_years",
    "age"
  )
)

if (!is.na(age_var)) {
  add_continuous_row(
    "Demographics",
    "Age, mean (SD), y",
    age_var
  )
}

# Report one sex category rather than duplicating a nearly complementary pair.
sex_var <- first_existing(
  c(
    "sex_clean",
    "sex"
  )
)

if (!is.na(sex_var)) {
  sx_all <- tolower(
    safe_char(
      allcohort[[sex_var]]
    )
  )

  sx_dev <- tolower(
    safe_char(
      dev[[sex_var]]
    )
  )

  sx_24 <- tolower(
    safe_char(
      eval24[[sex_var]]
    )
  )

  add_binary_row(
    "Demographics",
    "Female",
    grepl(
      "female|sex_2|^2$",
      sx_all
    ),
    grepl(
      "female|sex_2|^2$",
      sx_dev
    ),
    grepl(
      "female|sex_2|^2$",
      sx_24
    )
  )
}

gcs_vars <- c(
  "gcs_eye_clean",
  "gcs_verbal_clean",
  "gcs_motor_clean"
)

if (
  all(
    gcs_vars %in%
      names(dt)
  )
) {
  gcs_total_all <-
    as.numeric(
      allcohort$gcs_eye_clean
    ) +
    as.numeric(
      allcohort$gcs_verbal_clean
    ) +
    as.numeric(
      allcohort$gcs_motor_clean
    )

  gcs_total_dev <-
    as.numeric(
      dev$gcs_eye_clean
    ) +
    as.numeric(
      dev$gcs_verbal_clean
    ) +
    as.numeric(
      dev$gcs_motor_clean
    )

  gcs_total_24 <-
    as.numeric(
      eval24$gcs_eye_clean
    ) +
    as.numeric(
      eval24$gcs_verbal_clean
    ) +
    as.numeric(
      eval24$gcs_motor_clean
    )

  add_binary_row(
    "Neurologic status",
    "GCS 13–15",
    gcs_total_all >=
      13 &
      gcs_total_all <=
        15,
    gcs_total_dev >=
      13 &
      gcs_total_dev <=
        15,
    gcs_total_24 >=
      13 &
      gcs_total_24 <=
        15
  )

  add_binary_row(
    "Neurologic status",
    "GCS 9–12",
    gcs_total_all >=
      9 &
      gcs_total_all <=
        12,
    gcs_total_dev >=
      9 &
      gcs_total_dev <=
        12,
    gcs_total_24 >=
      9 &
      gcs_total_24 <=
        12
  )

  add_binary_row(
    "Neurologic status",
    "GCS 3–8",
    gcs_total_all >=
      3 &
      gcs_total_all <=
        8,
    gcs_total_dev >=
      3 &
      gcs_total_dev <=
        8,
    gcs_total_24 >=
      3 &
      gcs_total_24 <=
        8
  )
}

# Selected, explicitly labelled mechanisms. "Other" is not used as a top-N
# catch-all because that would make the section difficult to interpret.
mechanism_var <- first_existing(
  c(
    "mechanism_clean",
    "primary_ecode"
  )
)

if (!is.na(mechanism_var)) {
  ma <- safe_char(
    allcohort[[mechanism_var]]
  )

  md <- safe_char(
    dev[[mechanism_var]]
  )

  m24 <- safe_char(
    eval24[[mechanism_var]]
  )

  candidate_mechanisms <- c(
    "Fall",
    "Transport",
    "Assault",
    "Firearm/explosive"
  )

  available <- candidate_mechanisms[
    candidate_mechanisms %in%
      unique(ma)
  ]

  if (length(available) > 0L) {
    mechanism_counts <- sort(
      table(
        factor(
          ma[
            ma %in%
              available
          ],
          levels =
            available
        )
      ),
      decreasing =
        TRUE
    )

    keep_mechanisms <- names(
      mechanism_counts
    )[
      seq_len(
        min(
          3L,
          length(
            mechanism_counts
          )
        )
      )
    ]

    for (lv in keep_mechanisms) {
      add_binary_row(
        "Selected injury mechanisms",
        lv,
        ma ==
          lv,
        md ==
          lv,
        m24 ==
          lv
      )
    }
  }
}

spine_rows <- c(
  "dx_thoracic_vertebral_fracture",
  "dx_lumbar_vertebral_fracture",
  "dx_tl_spinal_cord_injury_any",
  "dx_tl_root_cauda_any",
  "dx_cervical_spine_injury",
  "dx_tbi_any"
)

for (
  v in spine_rows[
    spine_rows %in%
      names(dt)
  ]
) {
  add_binary_row(
    "Spine and concomitant injury phenotype",
    pretty_name(v),
    allcohort[[v]] ==
      1L,
    dev[[v]] ==
      1L,
    eval24[[v]] ==
      1L
  )
}

if (
  all(
    c(
      "dx_thoracic_vertebral_fracture",
      "dx_lumbar_vertebral_fracture"
    ) %in%
      names(dt)
  )
) {
  # Insert a clinically useful combined-region descriptor.
  add_binary_row(
    "Spine and concomitant injury phenotype",
    "Thoracic + lumbar vertebral fractures",
    allcohort$dx_thoracic_vertebral_fracture ==
      1L &
      allcohort$dx_lumbar_vertebral_fracture ==
        1L,
    dev$dx_thoracic_vertebral_fracture ==
      1L &
      dev$dx_lumbar_vertebral_fracture ==
        1L,
    eval24$dx_thoracic_vertebral_fracture ==
      1L &
      eval24$dx_lumbar_vertebral_fracture ==
        1L
  )
}

# Keep only high-yield baseline conditions in the main table. The complete
# predictor dictionary remains in the supplement.
pmhx_candidates <- intersect(
  c(
    "pmhx_anticoagulant_therapy",
    "pmhx_functional_dependence",
    "pmhx_diabetes"
  ),
  names(dt)
)

for (v in pmhx_candidates) {
  add_binary_row(
    "Selected pre-existing conditions",
    pretty_name(v),
    allcohort[[v]] ==
      1L,
    dev[[v]] ==
      1L,
    eval24[[v]] ==
      1L
  )
}

table1 <- rbindlist(
  table1_rows,
  fill =
    TRUE
)

table1 <- rbindlist(
  list(
    data.table(
      Domain =
        "Cohort",
      Characteristic =
        "N",
      Overall =
        fmt_int(
          nrow(
            allcohort
          )
        ),
      `Development 2020–2023` =
        fmt_int(
          nrow(
            dev
          )
        ),
      `2024 temporal cohort` =
        fmt_int(
          nrow(
            eval24
          )
        ),
      `Absolute SMD` =
        "—"
    ),
    table1
  ),
  use.names =
    TRUE,
  fill =
    TRUE
)

fwrite(
  table1,
  file.path(
    SOURCE_DIR,
    "Table_1_Cohort_Characteristics.csv"
  )
)

ft1 <- make_pub_ft(
  table1,
  section_col =
    "Domain",
  note = paste(
    "Values are percentages unless otherwise indicated.",
    "Development contains 2020–2023 admissions; 2024 is the latest forward temporal evaluation cohort.",
    "Absolute standardized mean differences compare development with 2024; no large-sample P values are reported.",
    "Values with |SMD| ≥0.10 are bolded."
  )
)

# Bold the SMDs that cross the conventional display threshold.
table1_display <- sectionize_display(
  table1,
  "Domain"
)$data

smd_rows <- which(
  suppressWarnings(
    as.numeric(
      table1_display$`Absolute SMD`
    )
  ) >=
    0.10
)

if (length(smd_rows) > 0L) {
  ft1 <- bold(
    ft1,
    i =
      smd_rows,
    j =
      "Absolute SMD",
    part =
      "body"
  )
}

ft1 <- padding(
  ft1,
  padding =
    2,
  part =
    "all"
)

# ------------------------------------------------------------------------------
# 5. MAIN TABLE 2 — Probability / trajectory performance
# ------------------------------------------------------------------------------

headline_outputs <- c(
  "disposition_postacute",
  "disposition_death_hospice",
  "icu_any",
  "icu_ge8",
  "vent_any",
  "vent_ge8",
  "hlos_ge20",
  "recorded_tl_operation"
)

headline_labels <- c(
  "disposition_postacute" =
    "Post-acute disposition",
  "disposition_death_hospice" =
    "Death / hospice",
  "icu_any" =
    "Any ICU",
  "icu_ge8" =
    "ICU ≥8 d",
  "vent_any" =
    "Any ventilation",
  "vent_ge8" =
    "Ventilation ≥8 d",
  "hlos_ge20" =
    "Hospital LOS ≥20 d",
  "recorded_tl_operation" =
    "Recorded thoracolumbar operation"
)

table2 <- copy(
  classification_2024_ci[
    output_id %in%
      headline_outputs
  ]
)

table2[
  ,
  output_order :=
    match(
      output_id,
      headline_outputs
    )
]

setorder(
  table2,
  output_order
)

table2_display <- table2[
  ,
  .(
    Outcome =
      unname(
        headline_labels[
          output_id
        ]
      ),
    `2024 N` =
      fmt_int(n),
    Prevalence =
      fmt_pct(prevalence),
    `AUROC (95% CI)` =
      fmt_ci3(
        auroc,
        auroc_lower,
        auroc_upper
      ),
    `AUPRC (95% CI)` =
      fmt_ci3(
        auprc,
        auprc_lower,
        auprc_upper
      ),
    `PR lift` =
      paste0(
        sprintf(
          "%.1f",
          pr_lift
        ),
        "x"
      ),
    Brier =
      fmt3(brier),
    `Calibration intercept` =
      fmt_signed3(
        calibration_intercept
      ),
    `Calibration slope` =
      fmt3(
        calibration_slope
      )
  )
]

fwrite(
  table2_display,
  file.path(
    SOURCE_DIR,
    "Table_2_Probability_Performance.csv"
  )
)

ft2 <- make_pub_ft(
  table2_display,
  note = paste(
    "Confidence intervals are percentile intervals from 1,000 patient-level bootstrap resamples of fixed 2024 predictions; models were not refit.",
    "PR lift is AUPRC divided by outcome prevalence.",
    "No recalibration was learned in the evaluation cohort.",
    "Recorded thoracolumbar operation denotes treatment receipt during the index hospitalization, not operative indication, spinal stability, or neurologic prognosis."
  )
)

# ------------------------------------------------------------------------------
# 6. MAIN TABLE 3 — Continuous-duration performance
# ------------------------------------------------------------------------------

duration_labels <- c(
  "hospital_los" =
    "Hospital LOS",
  "icu_los_conditional" =
    "ICU LOS among ICU users",
  "ventilator_days_conditional" =
    "Ventilator duration among ventilated patients"
)

duration_2024_extra <- duration_temporal[
  evaluation_year ==
    2024L,
  .(
    outcome,
    improvement_vs_constant_median_pct
  )
]

table3 <- merge(
  duration_2024_ci,
  duration_2024_extra,
  by =
    "outcome",
  all.x =
    TRUE
)

table3[
  ,
  duration_order :=
    match(
      outcome,
      names(
        duration_labels
      )
    )
]

setorder(
  table3,
  duration_order
)

table3_display <- table3[
  ,
  .(
    Outcome =
      unname(
        duration_labels[
          outcome
        ]
      ),
    `2024 N` =
      fmt_int(N),
    `Median prediction MAE, d (95% CI)` =
      fmt_days_ci2(
        mae_estimate,
        mae_lower,
        mae_upper
      ),
    `MAE improvement vs development-median benchmark, %` =
      sprintf(
        "%.1f",
        improvement_vs_constant_median_pct
      ),
    `Median bias, d` =
      sprintf(
        "%+.2f",
        median_bias_estimate
      ),
    `Q10–Q90 coverage (95% CI)` =
      paste0(
        scales::percent(
          coverage_estimate,
          accuracy =
            0.1
        ),
        " (",
        scales::percent(
          coverage_lower,
          accuracy =
            0.1
        ),
        "–",
        scales::percent(
          coverage_upper,
          accuracy =
            0.1
        ),
        ")"
      ),
    `Median Q10–Q90 width, d (95% CI)` =
      fmt_days_ci2(
        median_width_estimate,
        median_width_lower,
        median_width_upper
      )
  )
]

fwrite(
  table3_display,
  file.path(
    SOURCE_DIR,
    "Table_3_Duration_Performance.csv"
  )
)

ft3 <- make_pub_ft(
  table3_display,
  note = paste(
    "The development-median benchmark predicts the corresponding development-cohort median duration for every patient.",
    "Q10/Q50/Q90 models were fit on log1p-transformed duration and back-transformed to days with noncrossing row-wise rearrangement.",
    "The displayed Q10–Q90 range is an estimated 10th–90th percentile range, not a guaranteed 80% prediction interval.",
    "Discrete-day quantile evaluation uses a 1e-6-day numerical equality tolerance."
  )
)

# ------------------------------------------------------------------------------
# 7. eTABLE 1 — Cohort construction / analytic populations
# ------------------------------------------------------------------------------

flow_n <- function(
    step_name
) {
  d <- cohort_flow[
    step ==
      step_name
  ]

  pooled <- d[
    tolower(
      as.character(
        admission_year
      )
    ) ==
      "pooled",
    n
  ]

  if (length(pooled) == 1L) {
    return(
      as.numeric(
        pooled
      )
    )
  }

  sum(
    d[
      as.character(
        admission_year
      ) %in%
        as.character(
          2020:2024
        ),
      n
    ],
    na.rm =
      TRUE
  )
}

n_source_puf <- flow_n(
  "All PUF_TRAUMA records"
)

n_age_eligible <- flow_n(
  "Recorded age 18-89"
)

n_narrow <- flow_n(
  "Adult + narrow TL phenotype"
)

n_transfer_in <- flow_n(
  "Transfer-in = Yes"
)

n_transfer_unknown <- flow_n(
  "Transfer status unknown"
)

n_direct <- flow_n(
  "Known direct presentation"
)

n_transfer_out <- flow_n(
  "Future transfer-out among direct"
)

n_ama <- flow_n(
  "Future AMA among direct"
)

n_custody <- flow_n(
  "Future custody/ED-other among direct"
)

n_primary <- flow_n(
  "PRIMARY complete-trajectory cohort"
)

etable1a <- data.table(
  Stage = c(
    "Source ACS TQP PUF records",
    "Recorded age 18-89 years",
    "Adult narrow thoracolumbar phenotype",
    "Transfer-in excluded",
    "Unknown transfer status excluded",
    "Known direct presentations",
    "Subsequent transfer-out excluded",
    "Subsequent discharge against medical advice excluded",
    "Subsequent custody / ED-other trajectory loss excluded",
    "Retained complete index-hospital trajectory"
  ),
  N_num = c(
    n_source_puf,
    n_age_eligible,
    n_narrow,
    n_transfer_in,
    n_transfer_unknown,
    n_direct,
    n_transfer_out,
    n_ama,
    n_custody,
    n_primary
  ),
  Den_num = c(
    NA_real_,
    n_source_puf,
    n_age_eligible,
    n_narrow,
    n_narrow,
    n_narrow,
    n_direct,
    n_direct,
    n_direct,
    n_direct
  )
)

etable1a[
  ,
  `% of denominator` :=
    ifelse(
      is.na(
        Den_num
      ),
      "—",
      fmt_pct(
        N_num /
          Den_num
      )
    )
]

etable1a[
  ,
  N :=
    fmt_int(
      N_num
    )
]

etable1a[
  ,
  Denominator :=
    ifelse(
      is.na(
        Den_num
      ),
      "—",
      fmt_int(
        Den_num
      )
    )
]

etable1a <- etable1a[
  ,
  .(
    Stage,
    N,
    Denominator,
    `% of denominator`
  )
]

status_steps <- c(
  "PRIMARY complete-trajectory cohort",
  "Future transfer-out among direct",
  "Future AMA among direct",
  "Future custody/ED-other among direct"
)

status_labels <- c(
  "PRIMARY complete-trajectory cohort" =
    "Retained complete trajectory",
  "Future transfer-out among direct" =
    "Future transfer-out",
  "Future AMA among direct" =
    "Future AMA",
  "Future custody/ED-other among direct" =
    "Future custody / ED-other"
)

etable1b <- dcast(
  cohort_flow[
    step %in%
      status_steps &
      as.character(
        admission_year
      ) %in%
        as.character(
          2020:2024
        ),
    .(
      `Index-hospital trajectory status` =
        unname(
          status_labels[
            step
          ]
        ),
      admission_year =
        as.character(
          admission_year
        ),
      N =
        n
    )
  ],
  `Index-hospital trajectory status` ~
    admission_year,
  value.var =
    "N",
  fill =
    0
)

etable1b[
  ,
  Total :=
    rowSums(
      .SD
    ),
  .SDcols =
    intersect(
      as.character(
        2020:2024
      ),
      names(
        etable1b
      )
    )
]

for (
  v in setdiff(
    names(
      etable1b
    ),
    "Index-hospital trajectory status"
  )
) {
  etable1b[
    ,
    (v) :=
      fmt_int(
        get(v)
      )
  ]
}

class24 <- classification_temporal[
  evaluation_year ==
    2024L
]

get_class24 <- function(id) {
  class24[
    output_id ==
      id
  ]
}

duration24 <- duration_temporal[
  evaluation_year ==
    2024L
]

endpoint_definitions <- data.table(
  `Outcome / model` = c(
    "Disposition",
    "Hospital LOS trajectory",
    "ICU trajectory",
    "Mechanical ventilation trajectory",
    "Recorded thoracolumbar operation",
    "Continuous hospital LOS",
    "Conditional ICU LOS",
    "Conditional ventilator duration"
  ),
  `Operational definition` = c(
    "Home/home health; post-acute facility; death/hospice",
    "≤7 d; 8–19 d; ≥20 d",
    "No ICU; ICU 1–7 d; ICU ≥8 d",
    "No ventilation; ventilation 1–7 d; ventilation ≥8 d",
    "Therapeutic thoracolumbar fusion/fixation/decompression/augmentation/excision/revision code map; treatment receipt",
    "Q10/Q50/Q90 hospital days among patients with known hospital LOS",
    "Q10/Q50/Q90 ICU days among ICU users",
    "Q10/Q50/Q90 ventilator days among ventilated patients"
  ),
  `2024 N` = c(
    fmt_int(
      get_class24(
        "disposition_home"
      )$n[1L]
    ),
    fmt_int(
      get_class24(
        "hlos_ge20"
      )$n[1L]
    ),
    fmt_int(
      get_class24(
        "icu_any"
      )$n[1L]
    ),
    fmt_int(
      get_class24(
        "vent_any"
      )$n[1L]
    ),
    fmt_int(
      get_class24(
        "recorded_tl_operation"
      )$n[1L]
    ),
    fmt_int(
      duration24[
        outcome ==
          "hospital_los",
        N
      ][1L]
    ),
    fmt_int(
      duration24[
        outcome ==
          "icu_los_conditional",
        N
      ][1L]
    ),
    fmt_int(
      duration24[
        outcome ==
          "ventilator_days_conditional",
        N
      ][1L]
    )
  ),
  `2024 prevalence / analytic population` = c(
    paste0(
      "Home ",
      fmt_pct(
        get_class24(
          "disposition_home"
        )$prevalence[1L]
      ),
      "; post-acute ",
      fmt_pct(
        get_class24(
          "disposition_postacute"
        )$prevalence[1L]
      ),
      "; death/hospice ",
      fmt_pct(
        get_class24(
          "disposition_death_hospice"
        )$prevalence[1L]
      )
    ),
    paste0(
      "≥20 d ",
      fmt_pct(
        get_class24(
          "hlos_ge20"
        )$prevalence[1L]
      )
    ),
    paste0(
      "Any ICU ",
      fmt_pct(
        get_class24(
          "icu_any"
        )$prevalence[1L]
      ),
      "; ≥8 d ",
      fmt_pct(
        get_class24(
          "icu_ge8"
        )$prevalence[1L]
      )
    ),
    paste0(
      "Any ventilation ",
      fmt_pct(
        get_class24(
          "vent_any"
        )$prevalence[1L]
      ),
      "; ≥8 d ",
      fmt_pct(
        get_class24(
          "vent_ge8"
        )$prevalence[1L]
      )
    ),
    fmt_pct(
      get_class24(
        "recorded_tl_operation"
      )$prevalence[1L]
    ),
    "Patients with known hospital LOS",
    "ICU users",
    "Ventilated patients"
  )
)

# ------------------------------------------------------------------------------
# 8. eTABLE 2 — Predictor dictionary / missingness
# ------------------------------------------------------------------------------

predictor_rows <- lapply(
  final_context,
  function(v) {
    x <- dt[[v]]

    x24 <- dt[
      admission_year ==
        2024L
    ][[v]]

    data.table(
      Domain =
        predictor_domain(v),
      Predictor =
        pretty_name(v),
      `Analysis variable` =
        v,
      `Type / source` =
        paste0(
          predictor_type(v),
          "; ",
          predictor_source(v)
        ),
      `Overall missing / unknown` =
        paste0(
          fmt_int(
            missing_count(x)
          ),
          " (",
          fmt_pct(
            missing_count(x) /
              length(x)
          ),
          ")"
        ),
      `2024 missing / unknown` =
        paste0(
          fmt_int(
            missing_count(x24)
          ),
          " (",
          fmt_pct(
            missing_count(x24) /
              length(x24)
          ),
          ")"
        ),
      `Endpoint-use exception` =
        if (
          v %in%
            pg$social_context
        ) {
          "Disposition and categorical hospital-LOS models only"
        } else {
          ""
        }
    )
  }
)

etable2 <- rbindlist(
  predictor_rows,
  fill =
    TRUE
)

pupil_note <- ""

if (
  "pupil_clean" %in%
    etable2$`Analysis variable`
) {
  prow <- etable2[
    `Analysis variable` ==
      "pupil_clean"
  ][1L]

  pupil_note <- paste0(
    " Pupillary response was missing/unknown in ",
    prow$`Overall missing / unknown`,
    " overall and ",
    prow$`2024 missing / unknown`,
    " in 2024. Predictor-family adjudication and missingness-stratified robustness analyses are reported separately in the supplement."
  )
}


tripod_predictor_dictionary <- rbindlist(
  lapply(
    final_context,
    function(v) {
      data.table(
        `Display name` = pretty_name(v),
        `Analysis variable` = v,
        Domain = predictor_domain(v),
        Type = predictor_type(v),
        Units = predictor_units(v),
        `Observed categories / range` = implementation_levels(v),
        `Source / derivation` = predictor_source(v),
        `Prediction-time availability` =
          "Available after initial trauma-center evaluation and diagnostic workup",
        `Missing / unknown handling` =
          if (predictor_type(v) == "Categorical") {
            "Explicit __UNKNOWN__; unseen application-time levels map to __OTHER__"
          } else {
            "Native XGBoost missing-value routing"
          },
        `Final model use` =
          if (v %in% pg$social_context) {
            "Disposition and categorical hospital-LOS models only"
          } else {
            "All final clinical models"
          }
      )
    }
  ),
  fill = TRUE
)

# ------------------------------------------------------------------------------
# 9. eTABLE 3 — Predictor adjudication
# ------------------------------------------------------------------------------

stage2b <- if (
  !is.na(
    files$stage2b
  )
) {
  fread(
    files$stage2b
  )
} else {
  data.table()
}

stage2c <- if (
  !is.na(
    files$stage2c
  )
) {
  fread(
    files$stage2c
  )
} else {
  data.table()
}

range_text <- function(
    d,
    metric =
      "median_delta_auroc"
) {
  if (nrow(d) == 0L) {
    return(
      "Not available"
    )
  }

  vals <- d[[metric]]

  paste0(
    sprintf(
      "%+.4f",
      min(
        vals,
        na.rm = TRUE
      )
    ),
    " to ",
    sprintf(
      "%+.4f",
      max(
        vals,
        na.rm = TRUE
      )
    )
  )
}

variant_rows <- list(
  list(
    variant = "minus_associated_injury",
    family = "Associated injury phenotypes",
    construct = "Concomitant cranial, chest, abdominal/pelvic, and extremity injury phenotypes",
    rationale = "Polytrauma context available after diagnostic workup without AIS/ISS summary scores.",
    decision = "Retain"
  ),
  list(
    variant = "minus_pmhx",
    family = "Selected pre-existing conditions",
    construct = "Prespecified ACS pre-existing-condition fields",
    rationale = "Baseline comorbidity and physiologic reserve known at prediction time.",
    decision = "Retain"
  ),
  list(
    variant = "minus_gcs_qualifiers",
    family = "GCS examination qualifiers",
    construct = "Intubated; sedated/paralyzed; eye obstruction; unknown/unassessable",
    rationale = "Provides context for altered or unassessable GCS components.",
    decision = "Retain"
  ),
  list(
    variant = "minus_spine_soft_tissue",
    family = "Thoracolumbar soft-tissue injury",
    construct = "Ligament/dislocation/sprain flag within the qualifying cohort",
    rationale = "Adds treatment-pattern information beyond coarse fracture/neural flags.",
    decision = "Retain"
  ),
  list(
    variant = "plus_chest_detail",
    family = "Detailed chest injury",
    construct = "Rib, flail, pleural, pulmonary, diaphragm, thoracic vascular/cardiac detail",
    rationale = "Prediction-time critical-care severity information after diagnostic workup.",
    decision = "Add"
  ),
  list(
    variant = "plus_pupil",
    family = "Pupillary response",
    construct = "Admission pupillary response",
    rationale = "Early neurologic-severity information with acceptable timing.",
    decision = "Add"
  ),
  list(
    variant = "plus_airway_state",
    family = "Airway / respiratory-assistance state",
    construct = "Current respiratory-assistance / airway state",
    rationale = "Treatment-dependent near prediction time and potentially circular for ventilation outcomes.",
    decision = "Exclude from primary models"
  )
)

etable3a <- rbindlist(
  lapply(
    variant_rows,
    function(x) {
      d <- stage2b[
        variant ==
          x$variant
      ]

      data.table(
        `Predictor family` =
          x$family,
        `Variables / construct` =
          x$construct,
        `Clinical and timing rationale` =
          x$rationale,
        `Temporal performance effect` =
          paste0(
            "Median ΔAUROC across endpoints: ",
            range_text(d)
          ),
        `Final decision` =
          x$decision
      )
    }
  ),
  fill =
    TRUE
)

etable3b <- data.table(
  `Predictor-set comparison` = c(
    "Coarse vs granular thoracolumbar representation",
    "Thoracolumbar soft-tissue flag",
    "Race / ethnicity / payer for disposition and hospital LOS"
  ),
  `Temporal performance evidence` = c(
    if (
      nrow(stage2c) >
        0L
    ) {
      paste0(
        "Granular vs final clinical predictor set median ΔAUROC across endpoints: ",
        range_text(
          stage2c[
            variant ==
              "granular_plus_chest_pupil"
          ]
        ),
        "; differences were negligible."
      )
    } else {
      "Finalist comparison unavailable."
    },

    if (
      nrow(stage2c) >
        0L
    ) {
      paste0(
        "Removing the flag reduced operation AUROC by ",
        sprintf(
          "%.4f",
          abs(
            stage2c[
              variant ==
                "coarse_plus_chest_pupil_no_soft" &
                endpoint ==
                  "operative_spine_intervention",
              median_delta_auroc
            ][1L]
          )
        ),
        " and AUPRC by ",
        sprintf(
          "%.4f",
          abs(
            stage2c[
              variant ==
                "coarse_plus_chest_pupil_no_soft" &
                endpoint ==
                  "operative_spine_intervention",
              median_delta_auprc
            ][1L]
          )
        ),
        "."
      )
    } else {
      "Finalist comparison unavailable."
    },

    if (
      nrow(stage2c) >
        0L
    ) {
      paste0(
        "Context augmentation improved median AUROC by ",
        sprintf(
          "%+.4f",
          stage2c[
            variant ==
              "coarse_plus_chest_pupil_social" &
              endpoint ==
                "discharge4",
            median_delta_auroc
          ][1L]
        ),
        " for disposition and ",
        sprintf(
          "%+.4f",
          stage2c[
            variant ==
              "coarse_plus_chest_pupil_social" &
              endpoint ==
                "hlos3",
            median_delta_auroc
          ][1L]
        ),
        " for hospital LOS."
      )
    } else {
      "Finalist comparison unavailable."
    }
  ),
  `Final decision` = c(
    "Use the coarse representation for parsimony",
    "Retain in the clinical core",
    "Use only for disposition and categorical hospital-LOS models"
  )
)

# ------------------------------------------------------------------------------
# 10. eTABLE 4 — Complete temporal classification performance
# ------------------------------------------------------------------------------

etable4_display <- classification_temporal[
  ,
  .(
    Domain =
      classification_domain(
        output_id
      ),
    `Evaluation year` =
      as.character(
        evaluation_year
      ),
    Outcome =
      pretty_symbols(
        output_label
      ),
    N =
      fmt_int(n),
    Prevalence =
      fmt_pct(prevalence),
    AUROC =
      fmt3(auroc),
    AUPRC =
      fmt3(auprc),
    `PR lift` =
      paste0(
        sprintf(
          "%.2f",
          pr_lift
        ),
        "x"
      ),
    Brier =
      fmt3(brier),
    `Log loss` =
      fmt3(logloss),
    `Calibration intercept` =
      fmt_signed3(
        calibration_intercept
      ),
    `Calibration slope` =
      fmt3(
        calibration_slope
      )
  )
]

# ------------------------------------------------------------------------------
# 11. eTABLE 5 — Temporal durations + threshold rationale
# ------------------------------------------------------------------------------

etable5a <- duration_temporal[
  ,
  .(
    `Evaluation year` =
      as.character(
        evaluation_year
      ),
    Outcome =
      unname(
        duration_labels[
          outcome
        ]
      ),
    N =
      fmt_int(N),
    `MAE, d` =
      sprintf(
        "%.2f",
        MAE_days
      ),
    `MAE improvement vs development median, %` =
      sprintf(
        "%.1f",
        improvement_vs_constant_median_pct
      ),
    `Median bias, d` =
      sprintf(
        "%+.2f",
        median_bias_days
      ),
    `Mean pinball loss` =
      sprintf(
        "%.3f",
        mean_pinball
      ),
    `Q10–Q90 coverage` =
      fmt_pct(
        q10_q90_coverage_tolerance
      ),
    `Median Q10–Q90 width, d` =
      sprintf(
        "%.2f",
        median_q10_q90_width
      )
  )
]

dev_dt <- dt[
  admission_year %in%
    2020:2023
]

eval_dt <- dt[
  admission_year ==
    2024L
]

threshold_rows <- list()

add_threshold <- function(
    outcome,
    threshold,
    dev_x,
    eval_x
) {
  dev_x <- as.numeric(
    dev_x
  )

  eval_x <- as.numeric(
    eval_x
  )

  dev_x <- dev_x[
    is.finite(
      dev_x
    )
  ]

  eval_x <- eval_x[
    is.finite(
      eval_x
    )
  ]

  threshold_rows[[length(threshold_rows) + 1L]] <<- data.table(
    Outcome =
      outcome,
    Threshold =
      paste0(
        "≥",
        threshold,
        " d"
      ),
    `Development N` =
      length(dev_x),
    `Development prevalence` =
      mean(
        dev_x >=
          threshold
      ),
    `Approx. development percentile` =
      mean(
        dev_x <
          threshold
      ),
    `2024 prevalence` =
      mean(
        eval_x >=
          threshold
      )
  )
}

if (
  "hospital_los_days" %in%
    names(dt)
) {
  add_threshold(
    "Hospital LOS",
    20,
    dev_dt$hospital_los_days,
    eval_dt$hospital_los_days
  )
}

if (
  "icu_los_days" %in%
    names(dt)
) {
  add_threshold(
    "ICU LOS",
    8,
    dev_dt$icu_los_days,
    eval_dt$icu_los_days
  )
}

vent_var <- first_existing(
  c(
    "ventilator_days_clean",
    "ventilator_days"
  )
)

if (!is.na(vent_var)) {
  add_threshold(
    "Ventilator duration",
    8,
    dev_dt[[vent_var]],
    eval_dt[[vent_var]]
  )
}

etable5b <- rbindlist(
  threshold_rows,
  fill =
    TRUE
)

if (nrow(etable5b) > 0L) {
  etable5b[
    ,
    `Development N` :=
      fmt_int(
        `Development N`
      )
  ]

  etable5b[
    ,
    `Development prevalence` :=
      fmt_pct(
        `Development prevalence`
      )
  ]

  etable5b[
    ,
    `Approx. development percentile` :=
      fmt_pct(
        `Approx. development percentile`
      )
  ]

  etable5b[
    ,
    `2024 prevalence` :=
      fmt_pct(
        `2024 prevalence`
      )
  ]
}

# ------------------------------------------------------------------------------
# 12. eTABLE 6 — Model development and evaluation specification
# ------------------------------------------------------------------------------

class_configs <- if (
  !is.na(
    files$class_configs
  )
) {
  fread(
    files$class_configs
  )
} else {
  data.table()
}

duration_configs <- if (
  !is.na(
    files$duration_configs
  )
) {
  fread(
    files$duration_configs
  )
} else {
  data.table()
}

build_config_display <- function(
    d,
    model_family
) {
  if (nrow(d) == 0L) {
    return(
      data.table()
    )
  }

  endpoint_raw <- if (
    model_family ==
      "Classification"
  ) {
    col_or_na(
      d,
      c(
        "endpoint",
        "endpoint_id"
      ),
      ""
    )
  } else {
    col_or_na(
      d,
      c(
        "outcome",
        "endpoint",
        "endpoint_id"
      ),
      ""
    )
  }

  data.table(
    `Model family` =
      model_family,
    Output =
      humanize_model_endpoint(
        endpoint_raw
      ),
    `Evaluation year` =
      as.character(
        col_or_na(
          d,
          c(
            "evaluation_year",
            "eval_year"
          ),
          ""
        )
      ),
    `Predictor set` =
      humanize_predictor_policy(
        col_or_na(
          d,
          c(
            "predictor_policy",
            "variant"
          ),
          if (
            model_family ==
              "Classification"
          ) {
            "final_core"
          } else {
            "final_core"
          }
        )
      ),
    `No. predictors` =
      as.character(
        col_or_na(
          d,
          c(
            "n_predictors",
            "predictor_count"
          ),
          ""
        )
      ),
    `Max depth` =
      as.character(
        col_or_na(
          d,
          "max_depth",
          ""
        )
      ),
    `Min child weight` =
      as.character(
        col_or_na(
          d,
          "min_child_weight",
          ""
        )
      ),
    Subsample =
      as.character(
        col_or_na(
          d,
          "subsample",
          ""
        )
      ),
    `Column subsample` =
      as.character(
        col_or_na(
          d,
          "colsample_bytree",
          ""
        )
      ),
    `L2 penalty` =
      as.character(
        col_or_na(
          d,
          "lambda",
          ""
        )
      ),
    `Selected boosting rounds` =
      as.character(
        col_or_na(
          d,
          c(
            "selected_rounds",
            "best_iteration"
          ),
          ""
        )
      )
  )
}

etable6a <- rbindlist(
  list(
    build_config_display(
      class_configs,
      "Classification"
    ),
    build_config_display(
      duration_configs,
      "Quantile duration"
    )
  ),
  fill =
    TRUE
)

# Drop columns that are completely blank in the actual artifacts.
if (nrow(etable6a) > 0L) {
  blank_cols <- names(
    etable6a
  )[
    vapply(
      etable6a,
      function(x) {
        all(
          is.na(x) |
            trimws(
              as.character(x)
            ) ==
              ""
        )
      },
      logical(1)
    )
  ]

  blank_cols <- setdiff(
    blank_cols,
    c(
      "Model family",
      "Output"
    )
  )

  if (length(blank_cols) > 0L) {
    etable6a[
      ,
      (blank_cols) :=
        NULL
    ]
  }
}

etable6b <- data.table(
  Component = c(
    "Prediction-time population",
    "Primary XGBoost learning rate",
    "Structural search",
    "Rolling temporal evaluation",
    "Numeric missing data",
    "Categorical missing / unseen levels",
    "Probability calibration assessment",
    "Duration quantiles",
    "2024 confidence intervals",
    "Validation status"
  ),
  Specification = c(
    "Adults 18–89 years with the narrow acute thoracolumbar trauma phenotype, direct presentation, and an observable index-hospital trajectory",
    "eta = 0.05",
    "24 prespecified structural configurations evaluated during final model development",
    "Train 2020 / tune 2021 / evaluate 2022; train 2020–2021 / tune 2022 / evaluate 2023; train 2020–2022 / tune 2023 / evaluate 2024",
    "XGBoost native missing-value routing",
    "Explicit Unknown category; unseen application-time levels mapped to Other",
    "Calibration intercept/slope and flexible calibration curves; no evaluation-year recalibration",
    "Q10/Q50/Q90 on log1p(days), back-transformed with feasible support and noncrossing rearrangement",
    "1,000 patient-level bootstrap resamples of fixed predictions; models not refit",
    "Forward temporal evaluation is complete. Independent external/prospective validation is required before routine clinical deployment."
  )
)

# ------------------------------------------------------------------------------
# 13. eTABLE 7 — Model and predictor sensitivity
# ------------------------------------------------------------------------------

ridge <- if (
  !is.na(
    files$ridge
  )
) {
  fread(
    files$ridge
  )
} else {
  data.table()
}

etable7a <- data.table()

if (nrow(ridge) > 0L) {
  etable7a <- ridge[
    evaluation_year ==
      2024L &
      metric_level %in%
        c(
          "macro",
          "binary"
        ),
    .(
      Outcome =
        humanize_model_endpoint(
          endpoint
        ),
      `XGBoost AUROC` =
        fmt3(
          xgb_auroc
        ),
      `Ridge AUROC` =
        fmt3(
          ridge_auroc
        ),
      `ΔAUROC` =
        fmt_signed3(
          delta_auroc_xgb_minus_ridge
        ),
      `XGBoost AUPRC` =
        fmt3(
          xgb_auprc
        ),
      `Ridge AUPRC` =
        fmt3(
          ridge_auprc
        ),
      `XGBoost Brier` =
        fmt3(
          xgb_brier
        ),
      `Ridge Brier` =
        fmt3(
          ridge_brier
        )
    )
  ]
}

variant_public <- c(
  "minus_associated_injury" =
    "Remove associated-injury phenotypes",
  "minus_pmhx" =
    "Remove selected pre-existing conditions",
  "minus_gcs_qualifiers" =
    "Remove GCS examination qualifiers",
  "minus_spine_soft_tissue" =
    "Remove thoracolumbar soft-tissue flag",
  "plus_chest_detail" =
    "Add detailed chest-injury phenotypes",
  "plus_pupil" =
    "Add pupillary response",
  "plus_airway_state" =
    "Add airway / respiratory-assistance state",
  "plus_social_context" =
    "Add race / ethnicity / payer"
)

variant_decision <- c(
  "minus_associated_injury" =
    "Retain associated-injury phenotypes",
  "minus_pmhx" =
    "Retain selected pre-existing conditions",
  "minus_gcs_qualifiers" =
    "Retain GCS qualifiers",
  "minus_spine_soft_tissue" =
    "Retain thoracolumbar soft-tissue flag",
  "plus_chest_detail" =
    "Include detailed chest injury",
  "plus_pupil" =
    "Include pupillary response",
  "plus_airway_state" =
    "Exclude from primary models",
  "plus_social_context" =
    "Use only for disposition and hospital-LOS models"
)

etable7b <- data.table()

if (nrow(stage2b) > 0L) {
  etable7b <- stage2b[
    variant %in%
      names(
        variant_public
      ),
    .(
      `Predictor-family test` =
        unname(
          variant_public[
            variant[1L]
          ]
        ),
      `Median ΔAUROC across endpoints` =
        sprintf(
          "%+.4f",
          median(
            median_delta_auroc,
            na.rm = TRUE
          )
        ),
      `ΔAUROC range across endpoints` =
        paste0(
          sprintf(
            "%+.4f",
            min(
              median_delta_auroc,
              na.rm = TRUE
            )
          ),
          " to ",
          sprintf(
            "%+.4f",
            max(
              median_delta_auroc,
              na.rm = TRUE
            )
          )
        ),
      `Median ΔAUPRC across endpoints` =
        sprintf(
          "%+.4f",
          median(
            median_delta_auprc,
            na.rm = TRUE
          )
        ),
      `Final decision` =
        unname(
          variant_decision[
            variant[1L]
          ]
        )
    ),
    by =
      variant
  ][
    ,
    variant :=
      NULL
  ]
}

# ------------------------------------------------------------------------------
# 14. eTABLE 8 — Robustness
# ------------------------------------------------------------------------------

subgroup <- if (
  !is.na(
    files$subgroup
  )
) {
  fread(
    files$subgroup
  )
} else {
  data.table()
}

early_death <- if (
  !is.na(
    files$early_death
  )
) {
  fread(
    files$early_death
  )
} else {
  data.table()
}

airway <- if (
  !is.na(
    files$airway
  )
) {
  fread(
    files$airway
  )
} else {
  data.table()
}

broad <- if (
  !is.na(
    files$broad
  )
) {
  fread(
    files$broad
  )
} else {
  data.table()
}

composite_2024 <- if (
  !is.na(
    files$composite_2024
  )
) {
  fread(
    files$composite_2024
  )
} else {
  data.table()
}

key_domains <- c(
  "GCS severity",
  "Age",
  "Sex",
  "Clinical input missingness",
  "TL region",
  "SCI",
  "Neural pattern",
  "Concomitant cervical injury",
  "Concomitant TBI"
)

key_outputs <- c(
  "disposition_postacute",
  "disposition_death_hospice",
  "icu_ge8",
  "vent_ge8",
  "recorded_tl_operation"
)

subgroup_full_public <- data.table()
etable8a <- data.table()

if (nrow(subgroup) > 0L) {
  overall_auc <- classification_temporal[
    evaluation_year ==
      2024L,
    .(
      output_id,
      overall_AUROC =
        auroc
    )
  ]

  subgroup_full_public <- merge(
    subgroup[
      evaluation_year ==
        2024L &
        subgroup_domain %in%
          key_domains &
        output_id %in%
          key_outputs
    ],
    overall_auc,
    by =
      "output_id",
    all.x =
      TRUE
  )

  subgroup_full_public[
    ,
    delta_AUROC :=
      AUROC -
        overall_AUROC
  ]

  etable8a <- subgroup_full_public[
    ,
    {
      worst_idx <- which.min(
        delta_AUROC
      )

      list(
        `Strata evaluated` =
          uniqueN(
            subgroup_level
          ),
        `AUROC range` =
          paste0(
            sprintf(
              "%.3f",
              min(
                AUROC,
                na.rm = TRUE
              )
            ),
            "–",
            sprintf(
              "%.3f",
              max(
                AUROC,
                na.rm = TRUE
              )
            )
          ),
        `Largest AUROC decrease vs overall` =
          sprintf(
            "%+.3f",
            min(
              delta_AUROC,
              na.rm = TRUE
            )
          ),
        `Most affected subgroup / output` =
          paste0(
            humanize_subgroup_level(
              subgroup_domain[
                worst_idx
              ],
              subgroup_level[
                worst_idx
              ]
            ),
            " / ",
            pretty_symbols(
              fcase(
                output_id[
                  worst_idx
                ] ==
                  "recorded_tl_operation",
                "Recorded thoracolumbar operation",
                default =
                  output_label[
                    worst_idx
                  ]
              )
            )
          ),
        `Max |calibration intercept|` =
          sprintf(
            "%.3f",
            max(
              abs(
                calibration_intercept
              ),
              na.rm = TRUE
            )
          ),
        `Max |calibration slope−1|` =
          sprintf(
            "%.3f",
            max(
              abs(
                calibration_slope -
                  1
              ),
              na.rm = TRUE
            )
          )
      )
    },
    by =
      .(
        `Subgroup domain` =
          subgroup_domain
      )
  ]
}

etable8b <- data.table()

if (nrow(airway) > 0L) {
  etable8b <- airway[
    evaluation_year ==
      2024L &
      output_id %in%
        c(
          "vent_any",
          "vent_ge8"
        ),
    .(
      `Baseline airway state` =
        airway_state,
      Outcome =
        pretty_symbols(
          output_label
        ),
      N =
        fmt_int(N),
      Prevalence =
        fmt_pct(prevalence),
      AUROC =
        fmt3(AUROC),
      AUPRC =
        fmt3(AUPRC),
      `Calibration intercept` =
        fmt_signed3(
          calibration_intercept
        ),
      `Calibration slope` =
        fmt3(
          calibration_slope
        )
    )
  ]
}

etable8c <- data.table()

if (nrow(early_death) > 0L) {
  etable8c <- early_death[
    evaluation_year ==
      2024L &
      sensitivity_population !=
        "Overall",
    .(
      Section =
        "Early-death exclusion",
      Sensitivity =
        sensitivity_population,
      Outcome =
        pretty_symbols(
          output_label
        ),
      N =
        fmt_int(N),
      Prevalence =
        fmt_pct(prevalence),
      AUROC =
        fmt3(AUROC),
      AUPRC =
        fmt3(AUPRC),
      `Calibration intercept` =
        fmt_signed3(
          calibration_intercept
        ),
      `Calibration slope` =
        fmt3(
          calibration_slope
        )
    )
  ]
}

if (nrow(composite_2024) > 0L) {
  comp_rows <- composite_2024[
    ,
    .(
      Section =
        "Death-or-prolonged composite",
      Sensitivity =
        "Composite outcome",
      Outcome =
        pretty_symbols(
          endpoint_label
        ),
      N =
        fmt_int(N),
      Prevalence =
        fmt_pct(prevalence),
      AUROC =
        fmt3(AUROC),
      AUPRC =
        fmt3(AUPRC),
      `Calibration intercept` =
        fmt_signed3(
          calibration_intercept
        ),
      `Calibration slope` =
        fmt3(
          calibration_slope
        )
    )
  ]

  etable8c <- rbindlist(
    list(
      etable8c,
      comp_rows
    ),
    fill =
      TRUE
  )
}

etable8d <- data.table()

if (nrow(broad) > 0L) {
  desired_broad_outputs <- c(
    "Post-acute facility",
    "Death/hospice",
    ">=20 days",
    "ICU >=8 days",
    "Ventilation >=8 days",
    "Recorded thoracolumbar operative treatment"
  )

  etable8d <- broad[
    output_label %in%
      desired_broad_outputs,
    .(
      Population =
        fcase(
          population ==
            "Broad union",
          "Broad phenotype union",

          population ==
            "Broad-only sensitivity",
          "Broad-only phenotype",

          default =
            population
        ),
      Outcome =
        pretty_symbols(
          fcase(
            output_label ==
              ">=20 days",
            "Hospital LOS ≥20 d",

            output_label ==
              "Recorded thoracolumbar operative treatment",
            "Recorded thoracolumbar operation",

            default =
              output_label
          )
        ),
      N =
        fmt_int(N),
      Prevalence =
        fmt_pct(prevalence),
      AUROC =
        fmt3(AUROC),
      AUPRC =
        fmt3(AUPRC),
      Brier =
        fmt3(Brier),
      `Calibration intercept` =
        fmt_signed3(
          calibration_intercept
        ),
      `Calibration slope` =
        fmt3(
          calibration_slope
        )
    )
  ]
}

# ------------------------------------------------------------------------------
# 15. eTABLE 9 — Operative endpoint audit
# ------------------------------------------------------------------------------

map_dt <- if (
  !is.na(
    files$operative_map
  ) &&
    file.exists(
      files$operative_map
    )
) {
  fread(
    files$operative_map
  )
} else {
  data.table()
}

matched_dt <- if (
  !is.na(
    files$matched_procedures
  ) &&
    file.exists(
      files$matched_procedures
    )
) {
  as.data.table(
    readRDS(
      files$matched_procedures
    )
  )
} else {
  data.table()
}

cat(
  "\neTable 9A code-map audit:\n",
  "  operative_map: ",
  ifelse(
    is.na(files$operative_map),
    "NOT FOUND",
    files$operative_map
  ),
  "\n",
  "  map rows: ",
  nrow(map_dt),
  "\n",
  "  map columns: ",
  ifelse(
    ncol(map_dt) == 0L,
    "<none>",
    paste(
      names(map_dt),
      collapse = ", "
    )
  ),
  "\n",
  "  matched procedure rows: ",
  nrow(matched_dt),
  "\n",
  sep = ""
)

empty_operation_code_counts <- function() {
  data.table(
    Category =
      character(),
    `Unique ICD-10-PCS codes` =
      integer()
  )
}

operation_code_counts <- function(d) {
  if (
    is.null(d) ||
      nrow(d) ==
        0L
  ) {
    return(
      empty_operation_code_counts()
    )
  }

  category_col <- intersect(
    c(
      "spine_operation_category",
      "category",
      "operation_category",
      "Category"
    ),
    names(d)
  )

  code_col <- intersect(
    c(
      "procedure_code",
      "icd10pcs",
      "icd10_pcs",
      "icd10pcs_code",
      "pcs_code",
      "code",
      "ICDPROCEDURECODE"
    ),
    names(d)
  )

  if (
    length(category_col) <
      1L ||
      length(code_col) <
        1L
  ) {
    return(
      empty_operation_code_counts()
    )
  }

  out <- d[
    !is.na(
      get(
        code_col[1L]
      )
    ) &
      trimws(
        as.character(
          get(
            code_col[1L]
          )
        )
      ) !=
        "" &
      !is.na(
        get(
          category_col[1L]
        )
      ) &
      trimws(
        as.character(
          get(
            category_col[1L]
          )
        )
      ) !=
        "",
    .(
      `Unique ICD-10-PCS codes` =
        uniqueN(
          as.character(
            get(
              code_col[1L]
            )
          )
        )
    ),
    by =
      .(
        Category =
          as.character(
            get(
              category_col[1L]
            )
          )
      )
  ]

  if (nrow(out) == 0L) {
    return(
      empty_operation_code_counts()
    )
  }

  out
}

category_code_counts <- operation_code_counts(
  map_dt
)

if (
  nrow(category_code_counts) ==
    0L
) {
  category_code_counts <- operation_code_counts(
    matched_dt
  )
}

operative_pooled <- operative_categories[
  ,
  .(
    `Procedure rows` =
      sum(
        procedure_rows
      ),
    `Patients with category` =
      sum(
        patients
      )
  ),
  by =
    .(
      Category =
        as.character(
          spine_operation_category
        )
    )
]

# A missing optional code-map artifact should never terminate the entire
# manuscript-table build. If code-level data were found, merge the unique-code
# counts. Otherwise, emit the operation-category table without that column.
if (
  nrow(category_code_counts) >
    0L
) {
  etable9a <- merge(
    operative_pooled,
    category_code_counts,
    by =
      "Category",
    all.x =
      TRUE,
    sort =
      FALSE
  )
} else {
  etable9a <- copy(
    operative_pooled
  )

  warning(
    paste(
      "Could not locate a usable ICD-10-PCS code-level artifact for eTable 9A.",
      "The table will be produced without the 'Unique ICD-10-PCS codes' column."
    ),
    call. =
      FALSE
  )
}

etable9a[
  ,
  `Procedure rows` :=
    fmt_int(
      `Procedure rows`
    )
]

etable9a[
  ,
  `Patients with category` :=
    fmt_int(
      `Patients with category`
    )
]

# Never display an all-blank audit column.
if (
  "Unique ICD-10-PCS codes" %in%
    names(
      etable9a
    ) &&
    all(
      is.na(
        etable9a$`Unique ICD-10-PCS codes`
      )
    )
) {
  etable9a[
    ,
    `Unique ICD-10-PCS codes` :=
      NULL
  ]
}

etable9a_title <- if (
  "Unique ICD-10-PCS codes" %in%
    names(
      etable9a
    )
) {
  "Therapeutic thoracolumbar operation code-map categories"
} else {
  "Therapeutic thoracolumbar operation categories"
}

etable9b <- operative_rates[
  ,
  .(
    Year =
      as.character(
        admission_year
      ),
    N =
      fmt_int(n),
    `Operative events` =
      fmt_int(
        operative_events
      ),
    `Operative rate` =
      paste0(
        sprintf(
          "%.1f",
          operative_percent
        ),
        "%"
      )
  )
]

etable9c <- operative_timing[
  ,
  .(
    Year =
      as.character(
        admission_year
      ),
    `Matched procedure rows` =
      fmt_int(
        matched_procedure_rows
      ),
    `Matched patients` =
      fmt_int(
        matched_patients
      ),
    `Procedure day known` =
      fmt_pct(
        procedure_start_days_known_n /
          matched_procedure_rows
      ),
    `Procedure hour known` =
      fmt_pct(
        procedure_start_hrs_known_n /
          matched_procedure_rows
      ),
    `Negative timing values` =
      paste0(
        fmt_int(
          procedure_start_days_negative_n
        ),
        " day; ",
        fmt_int(
          procedure_start_hrs_negative_n
        ),
        " hour"
      )
  )
]


# ------------------------------------------------------------------------------
# 15B. eTABLE 10 — Temporal analytic sample sizes and outcome counts
# ------------------------------------------------------------------------------

dt[
  ,
  discharge3_tripod :=
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
      levels =
        c(
          "Home/home health",
          "Post-acute facility",
          "Death/hospice"
        )
    )
]

tripod_folds <- list(
  `2022` =
    list(
      train =
        2020L,
      tune =
        2021L,
      refit =
        2020:2021,
      evaluation =
        2022L
    ),
  `2023` =
    list(
      train =
        2020:2021,
      tune =
        2022L,
      refit =
        2020:2022,
      evaluation =
        2023L
    ),
  `2024` =
    list(
      train =
        2020:2022,
      tune =
        2023L,
      refit =
        2020:2023,
      evaluation =
        2024L
    )
)

class_specs <- list(
  `Disposition trajectory` =
    list(
      var =
        "discharge3_tripod",
      type =
        "factor"
    ),
  `Hospital LOS trajectory` =
    list(
      var =
        "hlos3",
      type =
        "factor"
    ),
  `ICU trajectory` =
    list(
      var =
        "icu3",
      type =
        "factor"
    ),
  `Ventilation trajectory` =
    list(
      var =
        "vent3",
      type =
        "factor"
    ),
  `Recorded thoracolumbar operation` =
    list(
      var =
        "operative_spine_intervention",
      type =
        "binary"
    )
)

vent_duration_var <- if (
  "ventilator_days_clean" %in%
    names(dt)
) {
  "ventilator_days_clean"
} else {
  "ventilator_days"
}

duration_specs <- list(
  `Hospital LOS quantiles` =
    list(
      var =
        "hospital_los_days",
      eligible =
        function(d) {
          rep(
            TRUE,
            nrow(d)
          )
        }
    ),
  `ICU LOS quantiles among ICU users` =
    list(
      var =
        "icu_los_days",
      eligible =
        function(d) {
          !is.na(
            d$icu3
          ) &
            as.character(
              d$icu3
            ) !=
              "No ICU"
        }
    ),
  `Ventilator-duration quantiles among ventilated patients` =
    list(
      var =
        vent_duration_var,
      eligible =
        function(d) {
          !is.na(
            d$vent3
          ) &
            as.character(
              d$vent3
            ) !=
              "No ventilation"
        }
    )
)

partition_year_label <- function(yrs) {
  yrs <- sort(
    unique(
      as.integer(yrs)
    )
  )

  if (length(yrs) == 1L) {
    return(
      as.character(yrs)
    )
  }

  paste0(
    min(yrs),
    "–",
    max(yrs)
  )
}

factor_partition_rows <- function(
    d,
    var,
    model_type,
    outcome,
    evaluation_year,
    partition,
    partition_years
) {
  x <- droplevels(
    factor(
      d[[var]]
    )
  )

  x <- x[
    !is.na(x)
  ]

  if (length(x) == 0L) {
    return(
      data.table(
        `Model type` =
          model_type,
        Outcome =
          outcome,
        `Evaluation year` =
          as.character(
            evaluation_year
          ),
        Partition =
          partition,
        `Partition years` =
          partition_year_label(
            partition_years
          ),
        `Analytic N` =
          0L,
        `Outcome level` =
          NA_character_,
        `Outcome level N` =
          NA_integer_
      )
    )
  }

  tab <- table(x)

  data.table(
    `Model type` =
      model_type,
    Outcome =
      outcome,
    `Evaluation year` =
      as.character(
        evaluation_year
      ),
    Partition =
      partition,
    `Partition years` =
      partition_year_label(
        partition_years
      ),
    `Analytic N` =
      length(x),
    `Outcome level` =
      pretty_symbols(
        names(tab)
      ),
    `Outcome level N` =
      as.integer(tab)
  )
}

binary_partition_rows <- function(
    d,
    var,
    model_type,
    outcome,
    evaluation_year,
    partition,
    partition_years
) {
  y <- suppressWarnings(
    as.integer(
      d[[var]]
    )
  )

  y <- y[
    !is.na(y)
  ]

  if (length(y) == 0L) {
    return(
      data.table(
        `Model type` =
          model_type,
        Outcome =
          outcome,
        `Evaluation year` =
          as.character(
            evaluation_year
          ),
        Partition =
          partition,
        `Partition years` =
          partition_year_label(
            partition_years
          ),
        `Analytic N` =
          0L,
        `Outcome level` =
          c(
            "Event",
            "Non-event"
          ),
        `Outcome level N` =
          c(
            0L,
            0L
          )
      )
    )
  }

  data.table(
    `Model type` =
      model_type,
    Outcome =
      outcome,
    `Evaluation year` =
      as.character(
        evaluation_year
      ),
    Partition =
      partition,
    `Partition years` =
      partition_year_label(
        partition_years
      ),
    `Analytic N` =
      length(y),
    `Outcome level` =
      c(
        "Event",
        "Non-event"
      ),
    `Outcome level N` =
      c(
        sum(
          y ==
            1L
        ),
        sum(
          y ==
            0L
        )
      )
  )
}

duration_partition_rows <- function(
    d,
    var,
    eligible,
    outcome,
    evaluation_year,
    partition,
    partition_years
) {
  y <- suppressWarnings(
    as.numeric(
      d[[var]]
    )
  )

  ok <-
    is.finite(y) &
    eligible(d)

  n_ok <- sum(ok)

  data.table(
    `Model type` =
      "Continuous duration",
    Outcome =
      outcome,
    `Evaluation year` =
      as.character(
        evaluation_year
      ),
    Partition =
      partition,
    `Partition years` =
      partition_year_label(
        partition_years
      ),
    `Analytic N` =
      n_ok,
    `Outcome level` =
      "Eligible duration records",
    `Outcome level N` =
      n_ok
  )
}

tripod_full_rows <- list()

for (eval_year in names(tripod_folds)) {
  fold <- tripod_folds[[eval_year]]

  partitions <- list(
    Train =
      fold$train,
    Tune =
      fold$tune,
    Refit =
      fold$refit,
    Evaluation =
      fold$evaluation
  )

  for (nm in names(class_specs)) {
    sp <- class_specs[[nm]]

    for (partition_name in names(partitions)) {
      yrs <- partitions[[partition_name]]

      dpart <- dt[
        admission_year %in%
          yrs
      ]

      row_block <- if (
        sp$type ==
          "factor"
      ) {
        factor_partition_rows(
          dpart,
          sp$var,
          "Classification",
          nm,
          eval_year,
          partition_name,
          yrs
        )
      } else {
        binary_partition_rows(
          dpart,
          sp$var,
          "Classification",
          nm,
          eval_year,
          partition_name,
          yrs
        )
      }

      tripod_full_rows[
        [
          length(tripod_full_rows) +
            1L
        ]
      ] <- row_block
    }
  }

  for (nm in names(duration_specs)) {
    sp <- duration_specs[[nm]]

    for (partition_name in names(partitions)) {
      yrs <- partitions[[partition_name]]

      dpart <- dt[
        admission_year %in%
          yrs
      ]

      tripod_full_rows[
        [
          length(tripod_full_rows) +
            1L
        ]
      ] <- duration_partition_rows(
        dpart,
        sp$var,
        sp$eligible,
        nm,
        eval_year,
        partition_name,
        yrs
      )
    }
  }
}

# Normalized machine-readable source table:
# one row per model / temporal fold / partition / outcome level.
etable10_full <- rbindlist(
  tripod_full_rows,
  fill =
    TRUE
)

# Compact reviewer-facing display:
# retain N for every partition but show detailed class/event counts only for
# the forward evaluation cohort. Complete train/tune/refit class/event counts
# remain available in etable10_full and are written to source_data below.
etable10_n <- unique(
  etable10_full[
    ,
    .(
      `Model type`,
      Outcome,
      `Evaluation year`,
      Partition,
      `Analytic N`
    )
  ]
)

etable10_n <- dcast(
  etable10_n,
  `Model type` +
    Outcome +
    `Evaluation year` ~
    Partition,
  value.var =
    "Analytic N"
)

evaluation_counts <- etable10_full[
  Partition ==
    "Evaluation" &
  `Model type` ==
    "Classification" &
  !is.na(
    `Outcome level`
  ),
  .(
    `Evaluation outcome counts` =
      paste(
        paste0(
          `Outcome level`,
          " ",
          fmt_int(
            `Outcome level N`
          )
        ),
        collapse =
          "; "
      )
  ),
  by =
    .(
      `Model type`,
      Outcome,
      `Evaluation year`
    )
]

etable10 <- merge(
  etable10_n,
  evaluation_counts,
  by =
    c(
      "Model type",
      "Outcome",
      "Evaluation year"
    ),
  all.x =
    TRUE,
  sort =
    FALSE
)

etable10[
  is.na(
    `Evaluation outcome counts`
  ),
  `Evaluation outcome counts` :=
    "—"
]

outcome_order <- c(
  names(class_specs),
  names(duration_specs)
)

etable10[
  ,
  `.model_order` :=
    match(
      `Model type`,
      c(
        "Classification",
        "Continuous duration"
      )
    )
]

etable10[
  ,
  `.outcome_order` :=
    match(
      Outcome,
      outcome_order
    )
]

etable10[
  ,
  `.year_order` :=
    as.integer(
      `Evaluation year`
    )
]

setorder(
  etable10,
  `.model_order`,
  `.year_order`,
  `.outcome_order`
)

for (
  v in intersect(
    c(
      "Train",
      "Tune",
      "Refit",
      "Evaluation"
    ),
    names(etable10)
  )
) {
  etable10[
    ,
    (v) :=
      fmt_int(
        get(v)
      )
  ]
}

setnames(
  etable10,
  old =
    intersect(
      c(
        "Train",
        "Tune",
        "Refit",
        "Evaluation"
      ),
      names(etable10)
    ),
  new =
    paste0(
      intersect(
        c(
          "Train",
          "Tune",
          "Refit",
          "Evaluation"
        ),
        names(etable10)
      ),
      " N"
    )
)

etable10[
  ,
  c(
    ".model_order",
    ".outcome_order",
    ".year_order"
  ) :=
    NULL
]

# ------------------------------------------------------------------------------
# 16. Write machine-readable source data
# ------------------------------------------------------------------------------

source_tables <- list(
  eTable_1A_Cohort_Flow =
    etable1a,
  eTable_1B_Observation_Status =
    etable1b,
  eTable_1C_Endpoint_Definitions =
    endpoint_definitions,
  eTable_2_Predictor_Dictionary =
    etable2,
  TRIPOD_AI_Predictor_Implementation_Dictionary =
    tripod_predictor_dictionary,
  eTable_3A_Family_Adjudication =
    etable3a,
  eTable_3B_Predictor_Set_Decisions =
    etable3b,
  eTable_4_Complete_Temporal_Classification =
    etable4_display,
  eTable_5A_Temporal_Duration =
    etable5a,
  eTable_5B_Threshold_Justification =
    etable5b,
  eTable_6A_Selected_Configs =
    etable6a,
  eTable_6B_Development_Evaluation_Specification =
    etable6b,
  eTable_7A_XGBoost_vs_Ridge =
    etable7a,
  eTable_7B_Predictor_Sensitivity_Summary =
    etable7b,
  eTable_8A_Subgroup_Robustness_Summary =
    etable8a,
  eTable_8A_FULL_Subgroup_Robustness =
    subgroup_full_public,
  eTable_8B_Airway =
    etable8b,
  eTable_8C_Early_Death_and_Composite =
    etable8c,
  eTable_8D_Broad_Phenotype =
    etable8d,
  eTable_9A_Operation_Categories =
    etable9a,
  eTable_9B_Annual_Operation_Rate =
    etable9b,
  eTable_9C_Operation_Timing_QC =
    etable9c,
  eTable_10_Temporal_Analytic_Sample_Sizes =
    etable10,
  TRIPOD_AI_eTable_10_FULL_Temporal_Analytic_Counts =
    etable10_full
)

for (nm in names(source_tables)) {
  fwrite(
    source_tables[[nm]],
    file.path(
      SOURCE_DIR,
      paste0(
        nm,
        ".csv"
      )
    )
  )
}

# ------------------------------------------------------------------------------
# 17. Build polished DOCX files
# ------------------------------------------------------------------------------

main_docx <- file.path(
  OUT_DIR,
  "Tables.docx"
)

main_doc <- read_docx()

main_doc <- body_set_default_section(
  main_doc,
  landscape_section
)

main_doc <- add_table_to_doc(
  main_doc,
  "Table 1.",
  "Cohort characteristics and temporal case-mix comparison",
  ft1
)

main_doc <- add_table_to_doc(
  main_doc,
  "Table 2.",
  "Final probability and trajectory performance in the 2024 forward temporal evaluation",
  ft2,
  page_break_before =
    TRUE
)

main_doc <- add_table_to_doc(
  main_doc,
  "Table 3.",
  "Continuous-duration performance in the 2024 forward temporal evaluation",
  ft3,
  page_break_before =
    TRUE
)

print(
  main_doc,
  target =
    main_docx
)

ft_etable2 <- make_pub_ft(
  etable2,
  section_col =
    "Domain",
  note = paste0(
    "All predictors were available at the intended prediction time. Numeric missingness uses native XGBoost routing; categorical missingness uses an explicit Unknown level, and unseen application-time levels map to Other. ",
    "The clinical core contains 54 predictors; race, ethnicity, and payer are added only to disposition and categorical hospital-LOS models. ",
    pupil_note
  )
)

ft_etable4 <- make_pub_ft(
  etable4_display,
  section_col =
    "Domain",
  note =
    "All values are from fixed forward temporal predictions. Multiclass outcomes are reported as class-specific one-vs-rest probabilities; no evaluation-year recalibration was learned."
)

ft_etable6a <- make_pub_ft(
  etable6a,
  section_col =
    if (
      "Model family" %in%
        names(
          etable6a
        )
    ) {
      "Model family"
    } else {
      NULL
    }
)

ft_etable8c <- make_pub_ft(
  etable8c,
  section_col =
    if (
      "Section" %in%
        names(
          etable8c
        )
    ) {
      "Section"
    } else {
      NULL
    }
)

supp_specs <- list(
  list(
    prefix =
      "eTable 1A.",
    title =
      "Complete cohort construction",
    ft =
      make_pub_ft(
        etable1a,
        note = paste(
          "Transfer-in and unknown transfer status are excluded before the intended-use direct-presentation cohort is defined.",
          "Subsequent transfer-out, discharge against medical advice, and custody/ED-other dispositions are excluded because the complete index-hospital trajectory is not observable."
        )
      )
  ),

  list(
    prefix =
      "eTable 1B.",
    title =
      "Index-hospital trajectory completeness by admission year",
    ft =
      make_pub_ft(
        etable1b
      )
  ),

  list(
    prefix =
      "eTable 1C.",
    title =
      "Endpoint definitions and 2024 analytic populations",
    ft =
      make_pub_ft(
        endpoint_definitions,
        note =
          "Recorded thoracolumbar operation is a treatment-receipt endpoint, not an operative-indication, stability, or neurologic-prognosis endpoint."
      )
  ),

  list(
    prefix =
      "eTable 2.",
    title =
      "Final predictor dictionary and missingness",
    ft =
      ft_etable2
  ),

  list(
    prefix =
      "eTable 3A.",
    title =
      "Predictor-family adjudication",
    ft =
      make_pub_ft(
        etable3a,
        note =
          "For removal analyses, a negative ΔAUROC supports retention; for addition analyses, a positive ΔAUROC supports inclusion. Final decisions also considered timing, calibration, burden, and conceptual redundancy."
      )
  ),

  list(
    prefix =
      "eTable 3B.",
    title =
      "Final predictor-set decisions",
    ft =
      make_pub_ft(
        etable3b
      )
  ),

  list(
    prefix =
      "eTable 4.",
    title =
      "Complete rolling temporal classification performance",
    ft =
      ft_etable4
  ),

  list(
    prefix =
      "eTable 5A.",
    title =
      "Complete rolling temporal continuous-duration performance",
    ft =
      make_pub_ft(
        etable5a,
        note =
          "Mean pinball loss is the average quantile loss across Q10, Q50, and Q90. Lower values indicate smaller quantile-prediction error."
      )
  ),

  list(
    prefix =
      "eTable 5B.",
    title =
      "Development-only trajectory-threshold justification",
    ft =
      make_pub_ft(
        etable5b,
        note =
          "Thresholds were fixed from development-period distributional review before final 2024 reporting. Ventilator duration ≥8 days denotes extended mechanical ventilation and is not presented as a universal prolonged-mechanical-ventilation definition."
      )
  ),

  list(
    prefix =
      "eTable 6A.",
    title =
      "Fold-specific selected structural configurations",
    ft =
      ft_etable6a
  ),

  list(
    prefix =
      "eTable 6B.",
    title =
      "Model development and evaluation specification",
    ft =
      make_pub_ft(
        etable6b
      )
  ),

  list(
    prefix =
      "eTable 7A.",
    title =
      "XGBoost versus ridge comparator in 2024",
    ft =
      make_pub_ft(
        etable7a,
        note =
          "For multiclass trajectories, comparator values are macro-averaged across classes; binary recorded-operation values are reported directly."
      )
  ),

  list(
    prefix =
      "eTable 7B.",
    title =
      "Predictor-family sensitivity summary",
    ft =
      make_pub_ft(
        etable7b
      )
  ),

  list(
    prefix =
      "eTable 8A.",
    title =
      "Summary of 2024 clinical subgroup robustness",
    ft =
      make_pub_ft(
        etable8a,
        section_col =
          if (
            "Subgroup domain" %in%
              names(
                etable8a
              )
          ) {
            "Subgroup domain"
          } else {
            NULL
          },
        note =
          "This table summarizes the complete prespecified subgroup grid. The full subgroup-by-output results are provided as machine-readable source data."
      )
  ),

  list(
    prefix =
      "eTable 8B.",
    title =
      "Ventilation performance by baseline airway state",
    ft =
      make_pub_ft(
        etable8b
      )
  ),

  list(
    prefix =
      "eTable 8C.",
    title =
      "Early-death and death-or-prolonged-resource sensitivities",
    ft =
      ft_etable8c
  ),

  list(
    prefix =
      "eTable 8D.",
    title =
      "Broad-phenotype sensitivity analysis",
    ft =
      make_pub_ft(
        etable8d,
        note =
          "The locked narrow-phenotype models were applied outside the intended-use phenotype without retraining or recalibration."
      )
  ),

  list(
    prefix =
      "eTable 9A.",
    title =
      etable9a_title,
    ft =
      make_pub_ft(
        etable9a,
        note =
          "Patients may contribute to more than one operation category when multiple therapeutic procedures are recorded during the index hospitalization."
      )
  ),

  list(
    prefix =
      "eTable 9B.",
    title =
      "Recorded thoracolumbar operation rate by year",
    ft =
      make_pub_ft(
        etable9b
      )
  ),

  list(
    prefix =
      "eTable 9C.",
    title =
      "Operative procedure timing quality control",
    ft =
      make_pub_ft(
        etable9c,
        note = paste(
          "The endpoint is derived from therapeutic thoracolumbar ICD-10-PCS procedures.",
          "Diagnostic, monitoring, injection/administration, inspection, removal-only, and other nontherapeutic codes are excluded from the curated map."
        )
      )
  ),

  list(
    prefix =
      "eTable 10.",
    title =
      "Temporal analytic sample sizes and evaluation outcome counts",
    ft =
      make_pub_ft(
        etable10,
        section_col =
          "Model type",
        note = paste(
          "Train, tune, and refit partitions contain only years preceding the corresponding forward temporal evaluation year.",
          "For compact presentation, the rendered table shows analytic N for every partition and exact class/event counts for the forward evaluation partition.",
          "Exact class/event counts for train, tune, refit, and evaluation partitions are retained in the accompanying machine-readable source data."
        )
      )
  )
)

supp_docx <- file.path(
  OUT_DIR,
  "Supplemental.docx"
)

supp_doc <- read_docx()

supp_doc <- body_set_default_section(
  supp_doc,
  landscape_section
)

for (i in seq_along(supp_specs)) {
  s <- supp_specs[[i]]

  supp_doc <- add_table_to_doc(
    supp_doc,
    s$prefix,
    s$title,
    s$ft,
    page_break_before =
      i >
        1L
  )
}

print(
  supp_doc,
  target =
    supp_docx
)


table_manifest <- data.table(
  exhibit = c(
    "Table 1", "Table 2", "Table 3",
    "eTable 1", "eTable 2", "eTable 3", "eTable 4", "eTable 5",
    "eTable 6", "eTable 7", "eTable 8", "eTable 9", "eTable 10"
  ),
  purpose = c(
    "Cohort characteristics and temporal case-mix comparison",
    "Final probability / trajectory performance",
    "Continuous-duration performance",
    "Cohort flow, observation status, endpoint definitions",
    "Final predictor dictionary and missingness",
    "Predictor adjudication and final predictor-set decisions",
    "Complete temporal classification performance",
    "Temporal duration performance and threshold justification",
    "Model development and evaluation specification",
    "XGBoost-versus-ridge and predictor-family sensitivity",
    "Subgroup, airway, early-death, composite, and broad-phenotype robustness",
    "Recorded thoracolumbar operative-treatment endpoint audit",
    "Temporal analytic sample sizes and evaluation outcome counts"
  )
)
fwrite(table_manifest, file.path(OUT_DIR, "TABLE_MANIFEST.csv"))

writeLines(
  c(
    "SPINE-TRACT PUBLIC / REVIEWER-FACING TABLE PIPELINE COMPLETE",
    "",
    paste0("Main tables: ", main_docx),
    paste0("Supplement: ", supp_docx),
    paste0("Machine-readable source data: ", SOURCE_DIR),
    "",
    "No model was fit, tuned, or recalibrated by this script."
  ),
  file.path(OUT_DIR, "TABLE_PIPELINE_SUMMARY.txt")
)

cat(
  "\n============================================================\n",
  "SPINE-TRACT FINAL TABLES COMPLETE\n",
  "Main: ", main_docx,
  "\nSupplement: ", supp_docx,
  "\n============================================================\n",
  sep = ""
)
