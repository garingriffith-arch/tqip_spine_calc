# ==============================================================================
# 31_build_SPINE_TRACT_final_figures_VERTICAL_FINAL_v4.R
#
# SPINE-TRACT PUBLIC / REVIEWER-FACING FIGURE PIPELINE
#
# Visual language intentionally mirrors the final TBI-TRACT paper while the
# content is adapted to acute thoracolumbar trauma.
#
# MAIN
#   Figure 1  Prediction framework / rolling temporal design / final architecture
#   Figure 2  Temporal performance fingerprint
#   Figure 3  2024 probability + continuous-duration calibration
#
# SUPPLEMENT
#   eFigure 1 Complete cohort flow + year-specific observation status
#   eFigure 2 Development-to-2024 case-mix drift love plot
#   eFigure 3 Calibration of all probability outputs in 2024
#   eFigure 4 Temporal duration-quantile calibration with discrete-day brackets
#   eFigure 5 Predictor-family adjudication heatmap
#   eFigure 6 Spine-specific subgroup robustness + contextual calibration
#
# Reviewer-facing principles:
#   - locked repaired-DX predictions only
#   - no ROC-curve wall or argmax confusion-matrix wall
#   - 2024 estimate + 95% CI foregrounded; pale bar = 2022-2024 range
#   - AUPRC displayed relative to prevalence
#   - flexible calibration curves use individual fixed predictions
#   - grouped observed proportions complement flexible calibration curves
#   - duration Q10-Q90 = estimated 10th-90th percentile range
#   - discrete-day quantiles use P(Y<Q) to P(Y<=Q) brackets
#   - predictor adjudication replaces generic feature-importance storytelling
#
# OUTPUT
#   outputs/SPINE_TRACT_FINAL_FIGURES/
#       PNG  300 dpi review copies
#       TIFF 600 dpi LZW submission copies
#
# This script does NOT fit, tune, or recalibrate the clinical prediction models.
# Evaluation-only calibration curves summarize fixed 2024 predictions.
# ==============================================================================

rm(list = ls())
gc()

required_packages <- c(
  "data.table",
  "ggplot2",
  "patchwork",
  "scales",
  "grid",
  "Matrix"
)

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
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
  library(Matrix)
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

figure_dir <- file.path(PROJECT_DIR, "outputs", "SPINE_TRACT_FINAL_FIGURES")
source_dir <- file.path(figure_dir, "source_data")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

search_roots <- c(
  file.path(PROJECT_DIR, "outputs"),
  file.path(PROJECT_DIR, "data", "derived")
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

  classification_temporal = find_file(
    "01_CLASSIFICATION_TEMPORAL_POINT_ESTIMATES.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),
  classification_ci = find_file(
    "04_CLASSIFICATION_2024_BOOTSTRAP_CI.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),
  duration_temporal = find_file(
    "07_DURATION_TEMPORAL_POINT_ESTIMATES.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
  ),
  duration_ci = find_file(
    "09_DURATION_2024_BOOTSTRAP_CI.csv",
    "stage5_manuscript_outputs_REPAIRED_DX"
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
  subgroup = find_file(
    "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv",
    "stage4_robustness_REPAIRED_DX",
    required = FALSE
  )
)

artifact_manifest <- data.table(
  artifact = names(files),
  path = unlist(files, use.names = FALSE)
)
artifact_manifest[, exists := file.exists(path)]
fwrite(artifact_manifest, file.path(figure_dir, "FIGURE_SOURCE_ARTIFACT_MANIFEST.csv"))

# ------------------------------------------------------------------------------
# 2. Locked data / final predictor policy
# ------------------------------------------------------------------------------

dt <- as.data.table(readRDS(files$primary))
spec <- readRDS(files$predictor_spec)
dt[, inc_key := as.character(inc_key)]
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

temporal <- fread(files$classification_temporal)
headline_ci <- fread(files$classification_ci)
duration_year <- fread(files$duration_temporal)
duration_ci <- fread(files$duration_ci)
cohort_flow <- fread(files$cohort_flow)

# ------------------------------------------------------------------------------
# 3. Visual system — same family as TBI-TRACT
# ------------------------------------------------------------------------------

BASE_FONT <- "Arial"
NAVY <- "#17365D"
BLUE <- "#4472C4"
GREEN <- "#548235"
PURPLE <- "#8064A2"
ORANGE <- "#C55A11"
GOLD <- "#BF9000"
RED <- "#C00000"
MID_GRAY <- "#7F7F7F"
LIGHT_GRAY <- "#D9D9D9"
VERY_LIGHT_GRAY <- "#F2F2F2"
BLACK <- "#1F1F1F"

domain_palette <- c(
  "Disposition" = BLUE,
  "Hospital LOS" = GOLD,
  "ICU" = GREEN,
  "Ventilation" = PURPLE,
  "Spine treatment" = ORANGE
)
domain_shapes <- c(
  "Disposition" = 16,
  "Hospital LOS" = 15,
  "ICU" = 17,
  "Ventilation" = 18,
  "Spine treatment" = 8
)

theme_tract <- function(base_size = 9) {
  theme_minimal(base_family = BASE_FONT, base_size = base_size) +
    theme(
      plot.title = element_text(
        face = "bold", color = NAVY, size = base_size + 2, hjust = 0
      ),
      plot.subtitle = element_text(color = MID_GRAY, size = base_size),
      plot.caption = element_text(color = MID_GRAY, size = base_size - 1, hjust = 0),
      axis.title = element_text(color = BLACK, face = "bold"),
      axis.text = element_text(color = BLACK),
      strip.text = element_text(face = "bold", color = NAVY),
      strip.background = element_rect(fill = VERY_LIGHT_GRAY, color = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.title = element_text(face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(8, 12, 8, 8)
    )
}

save_figure <- function(plot, stem, width, height) {
  png_file <- file.path(figure_dir, paste0(stem, ".png"))
  tiff_file <- file.path(figure_dir, paste0(stem, ".tiff"))

  ggsave(
    filename = png_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 300,
    bg = "white"
  )
  ggsave(
    filename = tiff_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = 600,
    device = "tiff",
    compression = "lzw",
    bg = "white"
  )
  invisible(list(png = png_file, tiff = tiff_file))
}

pretty_predictor <- function(x) {
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
    "gcsq_unknown_recovered" = "GCS qualifier: unknown",
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
    "dx_thoracic_vertebral_fracture" = "Thoracic fracture",
    "dx_lumbar_vertebral_fracture" = "Lumbar fracture",
    "dx_tl_spinal_cord_injury_any" = "Thoracolumbar spinal cord injury",
    "dx_tl_root_cauda_any" = "Root / cauda injury",
    "dx_tl_soft_tissue_any" = "Thoracolumbar ligament / sprain",
    "dx_cervical_spine_injury" = "Cervical spine injury",
    "dx_tbi_any" = "Traumatic brain injury",
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
    z <- tools::toTitleCase(z)
    z <- gsub("\\bCopd\\b", "COPD", z)
    z <- gsub("\\bChf\\b", "Congestive heart failure", z)
    z <- gsub("\\bTbi\\b", "TBI", z)
    z <- gsub("\\bTl\\b", "Thoracolumbar", z)
    z <- gsub("\\bGcs\\b", "GCS", z)
    out[miss] <- z
  }

  out
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
    z[
      z ==
        "SEX_1"
    ] <- "Male"

    z[
      z ==
        "SEX_2"
    ] <- "Female"

    z[
      z ==
        "SEX_3"
    ] <- "Non-binary"
  }

  z
}

humanize_context_level <- function(
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

  if (domain == "Ethnicity") {
    z[
      z %in%
        c(
          "ETH_1",
          "1"
        )
    ] <- "Hispanic or Latino"

    z[
      z %in%
        c(
          "ETH_2",
          "2"
        )
    ] <- "Not Hispanic or Latino"
  }

  if (domain == "Payer") {
    payer_map <- c(
      "PAYER_1" =
        "Medicaid",
      "PAYER_2" =
        "Not billed",
      "PAYER_3" =
        "Self-pay",
      "PAYER_4" =
        "Private/commercial",
      "PAYER_6" =
        "Medicare",
      "PAYER_7" =
        "Other government",
      "PAYER_10" =
        "Other",
      "1" =
        "Medicaid",
      "2" =
        "Not billed",
      "3" =
        "Self-pay",
      "4" =
        "Private/commercial",
      "6" =
        "Medicare",
      "7" =
        "Other government",
      "10" =
        "Other"
    )

    mapped <- unname(
      payer_map[z]
    )

    z[
      !is.na(
        mapped
      )
    ] <- mapped[
      !is.na(
        mapped
      )
    ]
  }

  z
}

# ------------------------------------------------------------------------------
# 4. Canonical target definitions and fixed-prediction loader
# ------------------------------------------------------------------------------

target_order <- c(
  "disposition_postacute",
  "disposition_death_hospice",
  "icu_any",
  "icu_ge8",
  "vent_any",
  "vent_ge8",
  "hlos_ge20",
  "recorded_tl_operation"
)

target_labels <- c(
  "disposition_postacute" = "Post-acute disposition",
  "disposition_death_hospice" = "Death / hospice",
  "icu_any" = "Any ICU",
  "icu_ge8" = "ICU \u22658 d",
  "vent_any" = "Any ventilation",
  "vent_ge8" = "Ventilation \u22658 d",
  "hlos_ge20" = "Hospital LOS \u226520 d",
  "recorded_tl_operation" = "Recorded thoracolumbar operation"
)

domain_for_target <- function(x) {
  fcase(
    x %in% c("disposition_postacute", "disposition_death_hospice"), "Disposition",
    x == "hlos_ge20", "Hospital LOS",
    x %in% c("icu_any", "icu_ge8"), "ICU",
    x %in% c("vent_any", "vent_ge8"), "Ventilation",
    x == "recorded_tl_operation", "Spine treatment",
    default = "Other"
  )
}

all_output_specs <- data.table(
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
    "Home / home health",
    "Post-acute disposition",
    "Death / hospice",
    "Hospital LOS \u22647 d",
    "Hospital LOS 8-19 d",
    "Hospital LOS \u226520 d",
    "No ICU",
    "ICU 1-7 d",
    "ICU \u22658 d",
    "Any ICU",
    "No ventilation",
    "Ventilation 1-7 d",
    "Ventilation \u22658 d",
    "Any ventilation",
    "Recorded TL operation"
  ),
  endpoint = c(
    rep("discharge3_final", 3),
    rep("hlos3_final", 3),
    rep("icu3_final", 4),
    rep("vent3_final", 4),
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
  )
)

find_prediction_file <- function(endpoint, year) {
  find_file(
    paste0(endpoint, "__eval", year, "__predictions.csv"),
    "stage3_final_classification_REPAIRED_DX"
  )
}

extract_binary_output <- function(year, output_id_value) {
  sp <- all_output_specs[output_id == output_id_value]
  if (nrow(sp) != 1L) {
    stop("Could not uniquely resolve output: ", output_id_value, call. = FALSE)
  }

  pred <- fread(find_prediction_file(sp$endpoint, year))
  pred[, inc_key := as.character(inc_key)]

  ep <- sp$endpoint
  target <- sp$target

  if (ep == "operative_spine_intervention") {
    return(data.table(
      admission_year = as.integer(pred$admission_year),
      inc_key = pred$inc_key,
      output_id = output_id_value,
      output_label = sp$output_label,
      y = as.integer(pred$y),
      p = as.numeric(pred$pred)
    ))
  }

  obs <- as.character(pred$observed_class)

  if (ep == "discharge3_final") {
    map <- c(
      "Home/home health" = "p_Home.home.health",
      "Post-acute facility" = "p_Post.acute.facility",
      "Death/hospice" = "p_Death.hospice"
    )
    p <- as.numeric(pred[[unname(map[target])]])
    y <- as.integer(obs == target)

  } else if (ep == "hlos3_final") {
    map <- c(
      "<=7 days" = "p_X..7.days",
      "8-19 days" = "p_X8.19.days",
      ">=20 days" = "p_X..20.days"
    )
    p <- as.numeric(pred[[unname(map[target])]])
    y <- as.integer(obs == target)

  } else if (ep == "icu3_final") {
    if (target == "ANY") {
      p <- as.numeric(pred$p_ICU.1.7.days) + as.numeric(pred$p_ICU...8.days)
      y <- as.integer(obs != "No ICU")
    } else {
      map <- c(
        "No ICU" = "p_No.ICU",
        "ICU 1-7 days" = "p_ICU.1.7.days",
        "ICU >=8 days" = "p_ICU...8.days"
      )
      p <- as.numeric(pred[[unname(map[target])]])
      y <- as.integer(obs == target)
    }

  } else if (ep == "vent3_final") {
    if (target == "ANY") {
      p <- as.numeric(pred$p_Ventilation.1.7.days) +
        as.numeric(pred$p_Ventilation...8.days)
      y <- as.integer(obs != "No ventilation")
    } else {
      map <- c(
        "No ventilation" = "p_No.ventilation",
        "Ventilation 1-7 days" = "p_Ventilation.1.7.days",
        "Ventilation >=8 days" = "p_Ventilation...8.days"
      )
      p <- as.numeric(pred[[unname(map[target])]])
      y <- as.integer(obs == target)
    }
  } else {
    stop("Unknown endpoint: ", ep, call. = FALSE)
  }

  data.table(
    admission_year = as.integer(pred$admission_year),
    inc_key = pred$inc_key,
    output_id = output_id_value,
    output_label = sp$output_label,
    y = y,
    p = p
  )
}

# ------------------------------------------------------------------------------
# 5. Evaluation-only flexible calibration
# ------------------------------------------------------------------------------

clamp_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

fit_calibration_curve <- function(d, grid_n = 150L, n_groups = 10L) {
  d <- d[is.finite(p) & !is.na(y)]
  if (nrow(d) < 200L || uniqueN(d$y) < 2L) {
    stop("Insufficient data for flexible calibration curve.", call. = FALSE)
  }

  d[, lp := qlogis(clamp_prob(p))]
  fit <- glm(
    y ~ splines::ns(lp, df = 4),
    data = d,
    family = binomial()
  )

  p_lo <- quantile(d$p, 0.005, na.rm = TRUE)
  p_hi <- quantile(d$p, 0.995, na.rm = TRUE)
  grid_p <- seq(p_lo, p_hi, length.out = grid_n)
  newd <- data.frame(lp = qlogis(clamp_prob(grid_p)))
  pr <- predict(fit, newdata = newd, type = "link", se.fit = TRUE)

  curve <- data.table(
    predicted_probability = grid_p,
    calibrated_probability = plogis(pr$fit),
    lower_95 = plogis(pr$fit - 1.96 * pr$se.fit),
    upper_95 = plogis(pr$fit + 1.96 * pr$se.fit)
  )

  d[, rank_fraction := frank(p, ties.method = "average") / .N]
  d[, group := pmin(n_groups, pmax(1L, ceiling(rank_fraction * n_groups)))]
  anchors <- d[
    ,
    .(
      N = .N,
      mean_predicted = mean(p),
      observed = mean(y)
    ),
    by = group
  ]

  list(curve = curve, anchors = anchors)
}

# ------------------------------------------------------------------------------
# 6. FIGURE 1 — Framework / temporal design / architecture
# ------------------------------------------------------------------------------

box_layer <- function(
    xmin, xmax, ymin, ymax, label,
    fill = "white", color = NAVY, text_color = NAVY,
    size = 3.0, fontface = "plain"
) {
  list(
    annotate(
      "rect",
      xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
      fill = fill, color = color, linewidth = 0.7
    ),
    annotate(
      "text",
      x = (xmin + xmax) / 2,
      y = (ymin + ymax) / 2,
      label = label,
      color = text_color,
      family = BASE_FONT,
      fontface = fontface,
      size = size,
      lineheight = 0.95
    )
  )
}

arrow_layer <- function(x, y, xend, yend, color = MID_GRAY) {
  annotate(
    "segment",
    x = x, y = y, xend = xend, yend = yend,
    color = color, linewidth = 0.7,
    arrow = arrow(length = unit(0.12, "in"), type = "closed")
  )
}

fig1_a <- ggplot() +
  box_layer(
    0.02, 0.19, 0.36, 0.64,
    "Trauma-center\narrival",
    fill = VERY_LIGHT_GRAY,
    fontface = "bold"
  ) +
  box_layer(
    0.25, 0.47, 0.29, 0.71,
    "Initial evaluation\n+ neurologic examination\n+ physiology\n+ diagnostic workup",
    fill = "#EAF0F8",
    fontface = "bold",
    size = 2.8
  ) +
  box_layer(
    0.53, 0.72, 0.27, 0.73,
    "SPINE-TRACT\nprediction point",
    fill = "#D9E2F3",
    size = 3.3,
    fontface = "bold"
  ) +
  box_layer(
    0.78, 0.98, 0.29, 0.71,
    "Future inpatient\nresource + disposition\n+ treatment trajectory",
    fill = "#FCE4D6",
    color = ORANGE,
    text_color = BLACK,
    fontface = "bold"
  ) +
  arrow_layer(0.19, 0.50, 0.25, 0.50) +
  arrow_layer(0.47, 0.50, 0.53, 0.50) +
  arrow_layer(0.72, 0.50, 0.78, 0.50) +
  annotate(
    "text",
    x = 0.625,
    y = 0.13,
    label = "Prediction occurs after the receiving trauma-center diagnostic workup but before subsequent inpatient utilization and operation.",
    family = BASE_FONT,
    color = MID_GRAY,
    size = 2.65
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = BASE_FONT) +
  labs(title = "A. Intended prediction time") +
  theme(
    plot.title = element_text(face = "bold", color = NAVY, size = 11, hjust = 0)
  )

fig1_b <- ggplot() +
  box_layer(0.02, 0.20, 0.68, 0.88, "Train\n2020", fill = "#EAF0F8", fontface = "bold") +
  box_layer(0.28, 0.46, 0.68, 0.88, "Tune\n2021", fill = "#FFF2CC", color = GOLD, text_color = BLACK, fontface = "bold") +
  box_layer(0.54, 0.72, 0.68, 0.88, "Evaluate\n2022", fill = "#E2F0D9", color = GREEN, text_color = BLACK, fontface = "bold") +
  arrow_layer(0.20, 0.78, 0.28, 0.78) +
  arrow_layer(0.46, 0.78, 0.54, 0.78) +

  box_layer(0.02, 0.20, 0.39, 0.59, "Train\n2020-2021", fill = "#EAF0F8", fontface = "bold") +
  box_layer(0.28, 0.46, 0.39, 0.59, "Tune\n2022", fill = "#FFF2CC", color = GOLD, text_color = BLACK, fontface = "bold") +
  box_layer(0.54, 0.72, 0.39, 0.59, "Evaluate\n2023", fill = "#E2F0D9", color = GREEN, text_color = BLACK, fontface = "bold") +
  arrow_layer(0.20, 0.49, 0.28, 0.49) +
  arrow_layer(0.46, 0.49, 0.54, 0.49) +

  box_layer(0.02, 0.20, 0.10, 0.30, "Train\n2020-2022", fill = "#EAF0F8", fontface = "bold") +
  box_layer(0.28, 0.46, 0.10, 0.30, "Tune\n2023", fill = "#FFF2CC", color = GOLD, text_color = BLACK, fontface = "bold") +
  box_layer(0.54, 0.72, 0.10, 0.30, "Evaluate\n2024", fill = "#E2F0D9", color = GREEN, text_color = BLACK, fontface = "bold") +
  arrow_layer(0.20, 0.20, 0.28, 0.20) +
  arrow_layer(0.46, 0.20, 0.54, 0.20) +

  box_layer(
    0.80, 0.98, 0.31, 0.67,
    "Model architecture frozen\nafter temporal evaluation\n\nNo evaluation-year\nrecalibration\n\nIndependent external\nvalidation next",
    fill = "#D9E2F3",
    fontface = "bold",
    size = 2.7
  ) +
  arrow_layer(0.72, 0.78, 0.80, 0.67, color = NAVY) +
  arrow_layer(0.72, 0.49, 0.80, 0.49, color = NAVY) +
  arrow_layer(0.72, 0.20, 0.80, 0.31, color = NAVY) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = BASE_FONT) +
  labs(title = "B. Rolling-origin temporal evaluation") +
  theme(
    plot.title = element_text(face = "bold", color = NAVY, size = 11, hjust = 0)
  )

fig1_c <- ggplot() +
  box_layer(
    0.02, 0.30, 0.33, 0.67,
    "Final clinical core\n\n54 predictors\n\nDirect-presentation intended use",
    fill = "#D9E2F3",
    size = 2.8,
    fontface = "bold"
  ) +
  box_layer(
    0.38, 0.62, 0.56, 0.90,
    "Disposition and hospital-LOS models\n\nDisposition trajectory\nHospital LOS trajectory\n\n+ race + ethnicity + payer\n(57 total inputs)",
    fill = "#EAF0F8",
    color = BLUE,
    text_color = BLACK,
    size = 2.7,
    fontface = "bold"
  ) +
  box_layer(
    0.38, 0.62, 0.10, 0.44,
    "ICU, ventilation, and operation models\n\nICU trajectory\nVentilation trajectory\nRecorded thoracolumbar operation",
    fill = "#E2F0D9",
    color = GREEN,
    text_color = BLACK,
    size = 2.7,
    fontface = "bold"
  ) +
  arrow_layer(0.30, 0.53, 0.38, 0.73, color = BLUE) +
  arrow_layer(0.30, 0.47, 0.38, 0.27, color = GREEN) +
  box_layer(
    0.70, 0.98, 0.56, 0.90,
    "Probability trajectories\n\nDisposition: 3 classes\nHLOS: \u22647 / 8-19 / \u226520 d\nICU: none / 1-7 / \u22658 d\nVentilation: none / 1-7 / \u22658 d\nRecorded operation: probability",
    fill = "white",
    color = ORANGE,
    text_color = BLACK,
    size = 2.5
  ) +
  box_layer(
    0.70, 0.98, 0.10, 0.44,
    "Continuous trajectories\n\nHospital LOS\nICU LOS among ICU users\nVentilator duration among ventilated patients\n\nMedian + estimated Q10-Q90 range",
    fill = "white",
    color = PURPLE,
    text_color = BLACK,
    size = 2.5
  ) +
  arrow_layer(0.62, 0.73, 0.70, 0.73, color = ORANGE) +
  arrow_layer(0.62, 0.27, 0.70, 0.27, color = PURPLE) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = BASE_FONT) +
  labs(title = "C. Final architecture and endpoint-specific predictor policy") +
  theme(
    plot.title = element_text(face = "bold", color = NAVY, size = 11, hjust = 0)
  )

figure_1 <-
  fig1_a /
  fig1_b /
  fig1_c +
  plot_layout(heights = c(0.72, 1.10, 1.05)) +
  plot_annotation(
    title = "Figure 1. SPINE-TRACT prediction framework and temporal evaluation",
    theme = theme(
      plot.title = element_text(
        family = BASE_FONT,
        face = "bold",
        color = NAVY,
        size = 14,
        hjust = 0
      )
    )
  )

save_figure(
  figure_1,
  "Figure_1_SPINE_TRACT_Framework",
  width = 10.5,
  height = 10.0
)

# ------------------------------------------------------------------------------
# 7. FIGURE 2 — Temporal performance fingerprint
# ------------------------------------------------------------------------------

figure2_data <- temporal[output_id %in% target_order]
figure2_data[, domain := domain_for_target(output_id)]
figure2_data[
  ,
  display_target := factor(
    unname(target_labels[output_id]),
    levels = rev(unname(target_labels[target_order]))
  )
]

summary2 <- figure2_data[
  ,
  .(
    AUROC_temporal_min = min(auroc),
    AUROC_temporal_max = max(auroc),
    PR_lift_temporal_min = min(pr_lift),
    PR_lift_temporal_max = max(pr_lift),
    AUROC_2024 = auroc[evaluation_year == 2024L][1L],
    PR_lift_2024 = pr_lift[evaluation_year == 2024L][1L],
    prevalence_2024 = prevalence[evaluation_year == 2024L][1L],
    domain = domain[1L]
  ),
  by = .(output_id, display_target)
]

ci2 <- headline_ci[output_id %in% target_order]
summary2 <- merge(
  summary2,
  ci2[
    ,
    .(
      output_id,
      auroc_lower,
      auroc_upper,
      auprc_lower,
      auprc_upper,
      prevalence
    )
  ],
  by = "output_id",
  all.x = TRUE
)
summary2[
  ,
  `:=`(
    PR_lift_low = auprc_lower / prevalence,
    PR_lift_high = auprc_upper / prevalence,
    domain = factor(domain, levels = names(domain_palette))
  )
]
fwrite(summary2, file.path(source_dir, "Figure_2_Temporal_Performance_Data.csv"))

pr_lift_upper <- max(
  c(summary2$PR_lift_temporal_max, summary2$PR_lift_high, summary2$PR_lift_2024),
  na.rm = TRUE
) * 1.05

p2a <- ggplot(summary2, aes(y = display_target, color = domain, shape = domain)) +
  geom_segment(
    aes(
      x = AUROC_temporal_min,
      xend = AUROC_temporal_max,
      yend = display_target
    ),
    linewidth = 2.5,
    alpha = 0.22
  ) +
  geom_segment(
    aes(x = auroc_lower, xend = auroc_upper, yend = display_target),
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = AUROC_2024),
    size = 2.8,
    fill = "white",
    stroke = 1.0
  ) +
  scale_color_manual(values = domain_palette) +
  scale_shape_manual(values = domain_shapes) +
  scale_x_continuous(
    limits = c(0.70, 0.97),
    breaks = seq(0.70, 0.95, by = 0.05)
  ) +
  labs(
    title = "A. Discrimination",
    subtitle = "2024 estimate (95% CI)\nPale bar = 2022-2024 range",
    x = "AUROC",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_tract(10) +
  guides(shape = "none") +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal"
  )

p2b <- ggplot(summary2, aes(y = display_target, color = domain, shape = domain)) +
  geom_segment(
    aes(
      x = PR_lift_temporal_min,
      xend = PR_lift_temporal_max,
      yend = display_target
    ),
    linewidth = 2.5,
    alpha = 0.22
  ) +
  geom_segment(
    aes(x = PR_lift_low, xend = PR_lift_high, yend = display_target),
    linewidth = 0.8
  ) +
  geom_point(
    aes(x = PR_lift_2024),
    size = 2.8,
    fill = "white",
    stroke = 1.0
  ) +
  geom_vline(xintercept = 1, color = MID_GRAY, linetype = 2, linewidth = 0.5) +
  scale_color_manual(values = domain_palette) +
  scale_shape_manual(values = domain_shapes) +
  scale_x_continuous(
    limits = c(1, pr_lift_upper),
    labels = function(x) paste0(x, "x")
  ) +
  labs(
    title = "B. Precision-recall performance",
    subtitle = "AUPRC divided by prevalence\n1x = no-skill baseline",
    x = "Precision-recall lift",
    y = NULL,
    color = NULL,
    shape = NULL
  ) +
  theme_tract(10) +
  guides(shape = "none") +
  theme(
    axis.text.y = element_blank(),
    legend.position = "none"
  )

figure_2 <-
  p2a + p2b +
    plot_layout(
      widths = c(1, 1)
    )

figure_2 <- figure_2 +
  plot_annotation(
    title = "Figure 2. Temporal performance of SPINE-TRACT outputs",
    theme = theme(
      plot.title = element_text(
        family = BASE_FONT,
        face = "bold",
        color = NAVY,
        size = 14,
        hjust = 0
      )
    )
  )

save_figure(
  figure_2,
  "Figure_2_Temporal_Performance_Fingerprint",
  width = 11.0,
  height = 6.5
)

# ------------------------------------------------------------------------------
# 8. FIGURE 3A — Curated 2024 probability calibration
# ------------------------------------------------------------------------------

calibration_targets <- c(
  "disposition_postacute",
  "disposition_death_hospice",
  "icu_ge8",
  "vent_ge8",
  "recorded_tl_operation"
)
calibration_curves <- list()
calibration_anchor_rows <- list()

for (id in calibration_targets) {
  d <- extract_binary_output(2024L, id)
  cal <- fit_calibration_curve(d)

  curve <- cal$curve
  curve[, `:=`(
    output_id = id,
    output_label = unname(target_labels[id])
  )]

  anchors <- cal$anchors
  anchors[, `:=`(
    output_id = id,
    output_label = unname(target_labels[id])
  )]

  calibration_curves[[length(calibration_curves) + 1L]] <- curve
  calibration_anchor_rows[[length(calibration_anchor_rows) + 1L]] <- anchors
}

calibration <- rbindlist(calibration_curves, fill = TRUE)
calibration_anchors <- rbindlist(calibration_anchor_rows, fill = TRUE)
fwrite(
  calibration,
  file.path(source_dir, "Figure_3A_Flexible_Calibration_Curves.csv")
)
fwrite(
  calibration_anchors,
  file.path(source_dir, "Figure_3A_Calibration_Anchors.csv")
)

make_calibration_panel <- function(output_id_value) {
  curve <- calibration[output_id == output_id_value][order(predicted_probability)]
  points <- calibration_anchors[output_id == output_id_value][order(mean_predicted)]
  metric <- headline_ci[output_id == output_id_value]

  upper <- max(
    c(
      curve$predicted_probability,
      curve$calibrated_probability,
      curve$upper_95,
      points$mean_predicted,
      points$observed
    ),
    na.rm = TRUE
  )
  upper <- min(1, max(0.12, upper * 1.05))
  domain <- domain_for_target(output_id_value)

  ggplot() +
    geom_abline(
      intercept = 0,
      slope = 1,
      color = MID_GRAY,
      linetype = 2,
      linewidth = 0.6
    ) +
    geom_ribbon(
      data = curve,
      aes(
        x = predicted_probability,
        ymin = lower_95,
        ymax = upper_95
      ),
      fill = domain_palette[domain],
      alpha = 0.13
    ) +
    geom_line(
      data = curve,
      aes(
        x = predicted_probability,
        y = calibrated_probability
      ),
      color = domain_palette[domain],
      linewidth = 0.9
    ) +
    geom_point(
      data = points,
      aes(
        x = mean_predicted,
        y = observed
      ),
      shape = 21,
      fill = "white",
      color = domain_palette[domain],
      size = 2.4,
      stroke = 0.7
    ) +
    annotate(
      "label",
      x = upper * 0.04,
      y = upper * 0.96,
      hjust = 0,
      vjust = 1,
      size = 2.25,
      family = BASE_FONT,
      label = paste0(
        "Intercept ",
        sprintf("%+.2f", metric$calibration_intercept[1L]),
        "\nSlope ",
        sprintf("%.2f", metric$calibration_slope[1L])
      ),
      label.size = 0.15,
      fill = alpha("white", 0.85),
      color = BLACK
    ) +
    scale_x_continuous(
      limits = c(0, upper),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, upper),
      labels = label_percent(accuracy = 1),
      expand = expansion(mult = c(0, 0.02))
    ) +

    coord_equal() +
    labs(
      title = unname(target_labels[output_id_value]),
      x = "Predicted",
      y = "Observed"
    ) +
    theme_tract(7.7) +
    theme(
      plot.title = element_text(size = 9, face = "bold", color = NAVY),
      legend.position = "none"
    )
}

cal_panels <- lapply(
  calibration_targets,
  make_calibration_panel
)

figure_3a_top <- wrap_plots(
  cal_panels[1:3],
  ncol = 3
)

figure_3a_bottom <-
  plot_spacer() +
  cal_panels[[4L]] +
  cal_panels[[5L]] +
  plot_spacer() +
  plot_layout(
    widths = c(0.5, 1, 1, 0.5)
  )

figure_3a <-
  figure_3a_top /
  figure_3a_bottom +
  plot_layout(
    heights = c(1, 1)
  )

# ------------------------------------------------------------------------------
# 9. FIGURE 3B — 2024 continuous-duration calibration
# ------------------------------------------------------------------------------

duration_order <- c(
  "hospital_los",
  "icu_los_conditional",
  "ventilator_days_conditional"
)
duration_labels <- c(
  "hospital_los" = "Hospital LOS",
  "icu_los_conditional" = "ICU LOS conditional on ICU use",
  "ventilator_days_conditional" = "Ventilator duration conditional on ventilation"
)
duration_colors <- c(
  "hospital_los" = GOLD,
  "icu_los_conditional" = GREEN,
  "ventilator_days_conditional" = PURPLE
)

find_duration_prediction <- function(outcome, year) {
  find_file(
    paste0(
      outcome,
      "__eval",
      year,
      "__Q10_Q50_Q90_predictions.csv"
    ),
    "stage3_duration_quantiles_REPAIRED_DX"
  )
}

make_duration_bins <- function(outcome, year, n_groups = 10L) {
  d <- fread(find_duration_prediction(outcome, year))
  d[, rank_fraction := frank(q50_median_days, ties.method = "average") / .N]
  d[, group := pmin(n_groups, pmax(1L, ceiling(rank_fraction * n_groups)))]

  out <- d[
    ,
    .(
      N = as.integer(.N),
      predicted_q10 = as.numeric(median(as.numeric(q10_days), na.rm = TRUE)),
      predicted_median = as.numeric(median(as.numeric(q50_median_days), na.rm = TRUE)),
      predicted_q90 = as.numeric(median(as.numeric(q90_days), na.rm = TRUE)),
      observed_q10 = as.numeric(
        quantile(as.numeric(observed_days), probs = 0.10, type = 1, na.rm = TRUE)
      ),
      observed_median = as.numeric(
        median(as.numeric(observed_days), na.rm = TRUE)
      ),
      observed_q90 = as.numeric(
        quantile(as.numeric(observed_days), probs = 0.90, type = 1, na.rm = TRUE)
      )
    ),
    by = group
  ]
  out[, `:=`(
    outcome = outcome,
    duration_label = unname(duration_labels[outcome]),
    evaluation_year = year
  )]
  out
}

duration_calibration <- rbindlist(
  lapply(duration_order, function(id) make_duration_bins(id, 2024L)),
  fill = TRUE
)
fwrite(
  duration_calibration,
  file.path(source_dir, "Figure_3B_Duration_Calibration_Bins.csv")
)

make_duration_panel <- function(duration_id_value) {
  d <- duration_calibration[outcome == duration_id_value][order(predicted_median)]

  ggplot(d, aes(x = predicted_median, y = observed_median)) +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = "dashed",
      linewidth = 0.45
    ) +
    geom_errorbar(
      aes(ymin = observed_q10, ymax = observed_q90),
      width = 0,
      linewidth = 0.45,
      color = duration_colors[duration_id_value]
    ) +
    geom_segment(
      aes(
        x = predicted_q10,
        xend = predicted_q90,
        y = observed_median,
        yend = observed_median
      ),
      linewidth = 0.45,
      color = duration_colors[duration_id_value]
    ) +
    geom_point(
      size = 1.8,
      color = duration_colors[duration_id_value]
    ) +
    labs(
      title = unname(duration_labels[duration_id_value]),
      x = "Predicted duration, days",
      y = "Observed duration, days"
    ) +
    coord_equal() +
    theme_classic(base_family = BASE_FONT) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.title = element_text(size = 8),
      axis.text = element_text(size = 7)
    )
}

duration_panels <- lapply(duration_order, make_duration_panel)
figure_3b <- wrap_plots(duration_panels, nrow = 1)

section_header <- function(label, subtitle = NULL) {
  ggplot() +
    annotate(
      "text",
      x = 0,
      y = 0.70,
      label = label,
      hjust = 0,
      family = BASE_FONT,
      fontface = "bold",
      color = NAVY,
      size = 4.2
    ) +
    {
      if (!is.null(subtitle)) {
        annotate(
          "text",
          x = 0,
          y = 0.18,
          label = subtitle,
          hjust = 0,
          family = BASE_FONT,
          color = MID_GRAY,
          size = 3.0
        )
      }
    } +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void(base_family = BASE_FONT) +
    theme(plot.margin = margin(0, 8, 0, 8))
}

figure_3 <-
  section_header(
    "A. Probability calibration in the 2024 temporal cohort",
    "Flexible calibration curves with 95% confidence bands and grouped observed proportions"
  ) /
  figure_3a /
  section_header(
    "B. Continuous resource-duration calibration in 2024",
    "Horizontal bars = predicted Q10–Q90; vertical bars = observed group Q10–Q90"
  ) /
  figure_3b +
  plot_layout(heights = c(0.14, 1.45, 0.14, 0.90)) +
  plot_annotation(
    title = "Figure 3. Calibration of SPINE-TRACT probability and continuous-duration predictions",
    theme = theme(
      plot.title = element_text(
        family = BASE_FONT,
        face = "bold",
        color = NAVY,
        size = 14,
        hjust = 0
      )
    )
  )

save_figure(
  figure_3,
  "Figure_3_Calibration",
  width = 12.0,
  height = 11.8
)

# ------------------------------------------------------------------------------
# 10. eFIGURE 1 — Complete cohort flow + year-specific observation status
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

n_source_puf <- flow_n("All PUF_TRAUMA records")
n_age_eligible <- flow_n("Recorded age 18-89")
n_narrow <- flow_n("Adult + narrow TL phenotype")
n_transfer_in <- flow_n("Transfer-in = Yes")
n_transfer_unknown <- flow_n("Transfer status unknown")
n_direct <- flow_n("Known direct presentation")
n_transfer_out <- flow_n("Future transfer-out among direct")
n_ama <- flow_n("Future AMA among direct")
n_custody <- flow_n("Future custody/ED-other among direct")
n_primary <- flow_n("PRIMARY complete-trajectory cohort")

ef1_a <- ggplot() +
  box_layer(
    0.28, 0.72, 0.82, 0.96,
    paste0(
      "Source ACS TQP PUF\n2020-2024\nN = ",
      comma(n_source_puf)
    ),
    fill = "#D9E2F3",
    fontface = "bold",
    size = 3.0
  ) +
  box_layer(
    0.28, 0.72, 0.63, 0.77,
    paste0(
      "Recorded age 18-89 years\nN = ",
      comma(n_age_eligible)
    ),
    fill = "#D9E2F3",
    fontface = "bold",
    size = 3.0
  ) +
  box_layer(
    0.28, 0.72, 0.44, 0.58,
    paste0(
      "Adult narrow thoracolumbar phenotype\nN = ",
      comma(n_narrow)
    ),
    fill = "#D9E2F3",
    fontface = "bold",
    size = 3.0
  ) +
  box_layer(
    0.28, 0.72, 0.25, 0.39,
    paste0(
      "Known direct presentations\nN = ",
      comma(n_direct)
    ),
    fill = "#D9E2F3",
    fontface = "bold",
    size = 3.0
  ) +
  box_layer(
    0.28, 0.72, 0.06, 0.20,
    paste0(
      "Retained complete index-hospital trajectory\nN = ",
      comma(n_primary)
    ),
    fill = "#D9E2F3",
    fontface = "bold",
    size = 3.0
  ) +
  box_layer(
    0.76, 0.99, 0.42, 0.61,
    paste0(
      "Excluded before direct presentation\n\nTransfer-in: ",
      comma(n_transfer_in),
      "\nUnknown transfer: ",
      comma(n_transfer_unknown)
    ),
    fill = "#FCE4D6",
    color = RED,
    text_color = BLACK,
    size = 2.45
  ) +
  box_layer(
    0.01, 0.24, 0.18, 0.42,
    paste0(
      "Excluded after direct presentation\n\nTransfer-out: ",
      comma(n_transfer_out),
      "\nAMA: ",
      comma(n_ama),
      "\nCustody / ED-other: ",
      comma(n_custody)
    ),
    fill = "#FCE4D6",
    color = RED,
    text_color = BLACK,
    size = 2.35
  ) +
  arrow_layer(0.50, 0.82, 0.50, 0.77, color = NAVY) +
  arrow_layer(0.50, 0.63, 0.50, 0.58, color = NAVY) +
  arrow_layer(0.50, 0.44, 0.50, 0.39, color = NAVY) +
  arrow_layer(0.50, 0.25, 0.50, 0.20, color = NAVY) +
  arrow_layer(0.72, 0.51, 0.76, 0.51, color = RED) +
  arrow_layer(0.28, 0.32, 0.24, 0.32, color = RED) +
  coord_cartesian(
    xlim = c(0, 1),
    ylim = c(0, 1),
    clip = "off"
  ) +
  theme_void(
    base_family = BASE_FONT
  )

efigure_1 <-
  ef1_a +
    plot_annotation(
      title = "eFigure 1. Cohort selection for the intended-use population",
      theme = theme(
        plot.title = element_text(
          family = BASE_FONT,
          face = "bold",
          color = NAVY,
          size = 14
        )
      )
    )

save_figure(
  efigure_1,
  "eFigure_1_Cohort_Flow",
  width = 10.5,
  height = 7.0
)


# ------------------------------------------------------------------------------
# 11. eFIGURE 2 — Development-to-2024 case-mix drift
# ------------------------------------------------------------------------------

smd_continuous <- function(x1, x2) {
  x1 <- as.numeric(x1)
  x2 <- as.numeric(x2)
  denom <- sqrt((var(x1, na.rm = TRUE) + var(x2, na.rm = TRUE)) / 2)
  if (!is.finite(denom) || denom == 0) return(0)
  (mean(x2, na.rm = TRUE) - mean(x1, na.rm = TRUE)) / denom
}

smd_binary <- function(p1, p2) {
  denom <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (!is.finite(denom) || denom == 0) return(0)
  (p2 - p1) / denom
}

max_level_smd <- function(x1, x2) {
  z1 <- trimws(as.character(x1))
  z2 <- trimws(as.character(x2))
  z1[is.na(z1) | z1 == ""] <- "__UNKNOWN__"
  z2[is.na(z2) | z2 == ""] <- "__UNKNOWN__"
  lev <- union(unique(z1), unique(z2))
  vals <- vapply(
    lev,
    function(lv) {
      abs(
        smd_binary(
          mean(z1 == lv),
          mean(z2 == lv)
        )
      )
    },
    numeric(1)
  )
  max(vals, na.rm = TRUE)
}

dev <- dt[admission_year %in% 2020:2023]
eval24 <- dt[admission_year == 2024L]

drift_rows <- lapply(final_context, function(v) {
  x1 <- dev[[v]]
  x2 <- eval24[[v]]

  if (is.numeric(x1) || is.integer(x1)) {
    vals <- unique(x1[!is.na(x1)])
    if (length(vals) <= 2L && all(vals %in% c(0, 1))) {
      smd <- abs(
        smd_binary(
          mean(x1 == 1L, na.rm = TRUE),
          mean(x2 == 1L, na.rm = TRUE)
        )
      )
      type <- "Binary"
    } else {
      smd <- abs(smd_continuous(x1, x2))
      type <- "Numeric"
    }
  } else {
    smd <- max_level_smd(x1, x2)
    type <- "Categorical: maximum level-specific SMD"
  }

  data.table(
    predictor = v,
    display_predictor = pretty_predictor(v),
    type = type,
    abs_smd = smd
  )
})

drift <- rbindlist(drift_rows, fill = TRUE)
setorder(drift, -abs_smd)
fwrite(drift, file.path(source_dir, "eFigure_2_Complete_Case_Mix_SMD.csv"))

drift_plot <- drift[seq_len(min(30L, .N))]
drift_plot[
  ,
  display_predictor := factor(
    display_predictor,
    levels = rev(display_predictor)
  )
]

efigure_2 <- ggplot(
  drift_plot,
  aes(x = abs_smd, y = display_predictor)
) +
  geom_vline(
    xintercept = 0.10,
    linetype = 2,
    color = RED,
    linewidth = 0.7
  ) +
  geom_segment(
    aes(x = 0, xend = abs_smd, yend = display_predictor),
    color = LIGHT_GRAY,
    linewidth = 0.7
  ) +
  geom_point(
    aes(fill = abs_smd >= 0.10),
    shape = 21,
    size = 2.7,
    color = NAVY
  ) +
  scale_fill_manual(
    values = c(`FALSE` = "white", `TRUE` = NAVY),
    guide = "none"
  ) +
  labs(
    title = "eFigure 2. Predictors with the largest development-to-2024 case-mix differences",
    subtitle = "Thirty final-model inputs with the largest absolute standardized mean differences; dashed line denotes |SMD| = 0.10",
    x = "Absolute standardized mean difference",
    y = NULL
  ) +
  theme_tract(9)

save_figure(
  efigure_2,
  "eFigure_2_Case_Mix_Drift_Love_Plot",
  width = 10.5,
  height = 8.5
)

# ------------------------------------------------------------------------------
# 12. eFIGURE 3 — Calibration of all probability outputs in 2024
# ------------------------------------------------------------------------------

all_calibration_curves <- list()
all_calibration_anchors <- list()

for (id in all_output_specs$output_id) {
  d <- extract_binary_output(2024L, id)
  cal <- fit_calibration_curve(d)

  curve <- cal$curve
  curve[, `:=`(
    output_id = id,
    output_label = all_output_specs[output_id == id, output_label]
  )]

  anchors <- cal$anchors
  anchors[, `:=`(
    output_id = id,
    output_label = all_output_specs[output_id == id, output_label]
  )]

  all_calibration_curves[[length(all_calibration_curves) + 1L]] <- curve
  all_calibration_anchors[[length(all_calibration_anchors) + 1L]] <- anchors
}

all_calibration <- rbindlist(all_calibration_curves, fill = TRUE)
all_anchors <- rbindlist(all_calibration_anchors, fill = TRUE)
fwrite(
  all_calibration,
  file.path(source_dir, "eFigure_3_Complete_2024_Calibration_Curves.csv")
)

make_complete_cal_panel <- function(id) {
  curve <- all_calibration[output_id == id][order(predicted_probability)]
  anchors <- all_anchors[output_id == id]

  upper <- max(
    c(
      curve$predicted_probability,
      curve$calibrated_probability,
      curve$upper_95,
      anchors$mean_predicted,
      anchors$observed
    ),
    na.rm = TRUE
  )
  upper <- min(1, max(0.12, upper * 1.05))

  ggplot() +
    geom_abline(
      intercept = 0,
      slope = 1,
      linetype = 2,
      color = MID_GRAY,
      linewidth = 0.45
    ) +
    geom_ribbon(
      data = curve,
      aes(
        x = predicted_probability,
        ymin = lower_95,
        ymax = upper_95
      ),
      fill = BLUE,
      alpha = 0.10
    ) +
    geom_line(
      data = curve,
      aes(
        x = predicted_probability,
        y = calibrated_probability
      ),
      color = NAVY,
      linewidth = 0.75
    ) +
    geom_point(
      data = anchors,
      aes(
        x = mean_predicted,
        y = observed
      ),
      shape = 21,
      fill = "white",
      color = NAVY,
      size = 1.4,
      stroke = 0.5
    ) +
    scale_x_continuous(
      limits = c(0, upper),
      labels = label_percent(accuracy = 1)
    ) +
    scale_y_continuous(
      limits = c(0, upper),
      labels = label_percent(accuracy = 1)
    ) +
    coord_equal() +
    labs(
      title = all_output_specs[output_id == id, output_label],
      x = "Predicted",
      y = "Observed"
    ) +
    theme_tract(6.8) +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      legend.position = "none"
    )
}

complete_panels <- lapply(all_output_specs$output_id, make_complete_cal_panel)

efigure_3 <- wrap_plots(complete_panels, ncol = 4) +
  plot_annotation(
    title = "eFigure 3. Calibration of all probability outputs in 2024",
    subtitle = "Flexible calibration curves with grouped observed proportions for all trajectory classes, derived probabilities, and recorded thoracolumbar operation",
    theme = theme(
      plot.title = element_text(
        family = BASE_FONT,
        face = "bold",
        color = NAVY,
        size = 14
      ),
      plot.subtitle = element_text(
        family = BASE_FONT,
        color = MID_GRAY,
        size = 9
      )
    )
  )

save_figure(
  efigure_3,
  "eFigure_3_Complete_2024_Calibration_Grid",
  width = 12.0,
  height = 12.8
)

# ------------------------------------------------------------------------------
# 13. eFIGURE 4 — Temporal duration quantile calibration
# ------------------------------------------------------------------------------

DURATION_TOL <- 1e-6
quantile_rows <- list()
coverage_rows <- list()

for (yr in 2022:2024) {
  for (id in duration_order) {
    d <- fread(find_duration_prediction(id, yr))

    y <- as.numeric(d$observed_days)
    q10 <- as.numeric(d$q10_days)
    q50 <- as.numeric(d$q50_median_days)
    q90 <- as.numeric(d$q90_days)

    q_defs <- list(
      Q10 = list(q = q10, nominal = 0.10),
      Q50 = list(q = q50, nominal = 0.50),
      Q90 = list(q = q90, nominal = 0.90)
    )

    for (qn in names(q_defs)) {
      qv <- q_defs[[qn]]$q
      quantile_rows[[length(quantile_rows) + 1L]] <- data.table(
        outcome = id,
        duration_label = unname(duration_labels[id]),
        evaluation_year = yr,
        quantile = qn,
        nominal = q_defs[[qn]]$nominal,
        empirical_lower = mean(y < qv - DURATION_TOL),
        empirical_upper = mean(y <= qv + DURATION_TOL)
      )
    }

    coverage_rows[[length(coverage_rows) + 1L]] <- data.table(
      outcome = id,
      duration_label = unname(duration_labels[id]),
      evaluation_year = yr,
      coverage = mean(
        y + DURATION_TOL >= q10 &
        y - DURATION_TOL <= q90
      )
    )
  }
}

quantile_reliability <- rbindlist(quantile_rows)
coverage_time <- rbindlist(coverage_rows)

fwrite(
  quantile_reliability,
  file.path(source_dir, "eFigure_4_Temporal_Quantile_Reliability.csv")
)
fwrite(
  coverage_time,
  file.path(source_dir, "eFigure_4_Temporal_Q10_Q90_Coverage.csv")
)

quantile_reliability[, midpoint := (empirical_lower + empirical_upper) / 2]

ef4_a <- ggplot(
  quantile_reliability,
  aes(
    x = nominal,
    y = midpoint,
    ymin = empirical_lower,
    ymax = empirical_upper,
    color = factor(evaluation_year),
    group = factor(evaluation_year)
  )
) +
  geom_abline(
    intercept = 0,
    slope = 1,
    linetype = 2,
    color = MID_GRAY,
    linewidth = 0.5
  ) +
  geom_errorbar(width = 0.015, linewidth = 0.6) +
  geom_point(size = 2.1) +
  facet_wrap(~ duration_label, nrow = 1) +
  scale_x_continuous(
    breaks = c(0.10, 0.50, 0.90),
    labels = c("Q10", "Q50", "Q90"),
    limits = c(0.05, 0.95)
  ) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1)
  ) +
  labs(
    title = "A. Quantile reliability",
    subtitle = "Vertical interval spans P(Y < predicted quantile) to P(Y ≤ predicted quantile), appropriate for discrete-day outcomes",
    x = "Predicted quantile",
    y = "Empirical cumulative probability",
    color = "Evaluation year"
  ) +
  theme_tract(9)

ef4_b <- ggplot(
  coverage_time,
  aes(
    x = evaluation_year,
    y = coverage,
    color = duration_label,
    group = duration_label
  )
) +
  geom_hline(
    yintercept = 0.80,
    linetype = 2,
    color = MID_GRAY,
    linewidth = 0.6
  ) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.2) +
  scale_x_continuous(breaks = 2022:2024) +
  scale_y_continuous(labels = label_percent(accuracy = 1)) +
  labs(
    title = "B. Empirical Q10–Q90 coverage",
    subtitle = "Dashed line = nominal 80% reference",
    x = "Evaluation year",
    y = "Empirical coverage",
    color = NULL
  ) +
  theme_tract(9)

efigure_4 <-
  ef4_a /
  ef4_b +
  plot_layout(heights = c(1.15, 0.85)) +
  plot_annotation(
    title = "eFigure 4. Temporal calibration of continuous duration quantiles",
    theme = theme(
      plot.title = element_text(
        family = BASE_FONT,
        face = "bold",
        color = NAVY,
        size = 14
      )
    )
  )

save_figure(
  efigure_4,
  "eFigure_4_Duration_Quantile_Calibration_Over_Time",
  width = 11.0,
  height = 8.8
)

# ------------------------------------------------------------------------------
# 14. eFIGURE 5 — Predictor-family adjudication
# ------------------------------------------------------------------------------

stage2b <- if (!is.na(files$stage2b)) fread(files$stage2b) else data.table()
stage2c <- if (!is.na(files$stage2c)) fread(files$stage2c) else data.table()

endpoint_labels <- c(
  "discharge4" = "Disposition",
  "hlos3" = "HLOS",
  "icu3" = "ICU",
  "vent3" = "Ventilation",
  "operative_spine_intervention" = "TL operation"
)

variant_labels_b <- c(
  "minus_associated_injury" = "Remove associated injury",
  "minus_pmhx" = "Remove selected pre-existing conditions",
  "minus_gcs_qualifiers" = "Remove GCS qualifiers",
  "minus_spine_soft_tissue" = "Remove TL soft-tissue flag",
  "plus_chest_detail" = "Add detailed chest injury",
  "plus_pupil" = "Add pupillary response",
  "plus_airway_state" = "Add airway state",
  "plus_social_context" = "Add race / ethnicity / payer"
)

heat_b <- stage2b[variant %in% names(variant_labels_b)]
if (nrow(heat_b) > 0L) {
  heat_b[, `:=`(
    Variant = factor(
      unname(variant_labels_b[variant]),
      levels = rev(unname(variant_labels_b))
    ),
    Endpoint = factor(
      unname(endpoint_labels[endpoint]),
      levels = unname(endpoint_labels)
    )
  )]
}

variant_labels_c <- c(
  "granular_plus_chest_pupil" = "Use granular thoracolumbar representation",
  "coarse_plus_chest_pupil_social" = "Add race / ethnicity / payer"
)

heat_c <- stage2c[variant %in% names(variant_labels_c)]
if (nrow(heat_c) > 0L) {
  heat_c[, `:=`(
    Variant = factor(
      unname(variant_labels_c[variant]),
      levels = rev(unname(variant_labels_c))
    ),
    Endpoint = factor(
      unname(endpoint_labels[endpoint]),
      levels = unname(endpoint_labels)
    )
  )]
}

fwrite(heat_b, file.path(source_dir, "eFigure_5_Family_Screening_Data.csv"))
fwrite(heat_c, file.path(source_dir, "eFigure_5_Finalist_Comparison_Data.csv"))

if (nrow(heat_b) > 0L && nrow(heat_c) > 0L) {
  heat_scale <- max(
    abs(c(heat_b$median_delta_auroc, heat_c$median_delta_auroc)),
    na.rm = TRUE
  )

  make_heat <- function(d, title) {
    ggplot(
      d,
      aes(
        x = Endpoint,
        y = Variant,
        fill = median_delta_auroc
      )
    ) +
      geom_tile(color = "white", linewidth = 0.4) +
      geom_text(
        aes(label = sprintf("%+.3f", median_delta_auroc)),
        size = 2.7,
        family = BASE_FONT
      ) +
      scale_fill_gradient2(
        low = RED,
        mid = "white",
        high = BLUE,
        midpoint = 0,
        limits = c(-heat_scale, heat_scale),
        name = "Delta AUROC"
      ) +
      labs(
        title = title,
        x = NULL,
        y = NULL
      ) +
      theme_tract(9) +
      theme(
        axis.text.x = element_text(angle = 35, hjust = 1)
      )
  }

  ef5_a <- make_heat(
    heat_b,
    "A. Predictor-family screening vs clinical reference predictor set"
  )
  ef5_b <- make_heat(
    heat_c,
    "B. Final predictor-set comparisons"
  )

  efigure_5 <-
    (
      ef5_a /
      ef5_b +
      plot_layout(heights = c(1.55, 0.85), guides = "collect")
    ) &
      theme(legend.position = "bottom")

  efigure_5 <-
    efigure_5 +
    plot_annotation(
      title = "eFigure 5. Temporal predictor-family adjudication",
      subtitle = "Positive ΔAUROC indicates improved discrimination relative to the stated reference; final decisions also considered timing, calibration, burden, and redundancy.",
      theme = theme(
        plot.title = element_text(
          family = BASE_FONT,
          face = "bold",
          color = NAVY,
          size = 14
        ),
        plot.subtitle = element_text(
          family = BASE_FONT,
          color = MID_GRAY,
          size = 9
        )
      )
    )

  save_figure(
    efigure_5,
    "eFigure_5_Predictor_Adjudication_Heatmap",
    width = 11.0,
    height = 9.2
  )
}

# ------------------------------------------------------------------------------
# 15. eFIGURE 6 — Spine-specific subgroup robustness + contextual calibration
# ------------------------------------------------------------------------------

if (
  is.na(
    files$subgroup
  )
) {
  warning(
    "Subgroup robustness file not found; eFigure 6 will be skipped."
  )
} else {
  subgroup <- fread(
    files$subgroup
  )

  overall_2024 <- temporal[
    evaluation_year ==
      2024L,
    .(
      output_id,
      overall_auroc =
        auroc
    )
  ]

  subgroup <- merge(
    subgroup,
    overall_2024,
    by =
      "output_id",
    all.x =
      TRUE
  )

  subgroup[
    ,
    delta_auroc :=
      AUROC -
        overall_auroc
  ]

  heat_outputs <- c(
    "disposition_postacute",
    "disposition_death_hospice",
    "icu_ge8",
    "vent_ge8",
    "recorded_tl_operation"
  )

  heat_domains <- c(
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

  # Display all prespecified clinical subgroup levels meeting the minimum-N
  # reporting rule rather than selecting the largest observed differences.
  sub_heat <- subgroup[
    evaluation_year ==
      2024L &
      output_id %in%
        heat_outputs &
      subgroup_domain %in%
        heat_domains &
      N >=
        500
  ]

  sub_heat[
    ,
    subgroup_display :=
      mapply(
        humanize_subgroup_level,
        subgroup_domain,
        subgroup_level,
        USE.NAMES =
          FALSE
      )
  ]

  sub_heat[
    ,
    subgroup_domain :=
      factor(
        subgroup_domain,
        levels =
          heat_domains
      )
  ]

  sub_heat[
    ,
    display_output :=
      factor(
        unname(
          target_labels[
            output_id
          ]
        ),
        levels =
          unname(
            target_labels[
              heat_outputs
            ]
          )
      )
  ]

  fwrite(
    sub_heat,
    file.path(
      source_dir,
      "eFigure_6_Subgroup_Robustness_Data.csv"
    )
  )

  limit6 <- max(
    abs(
      sub_heat$delta_auroc
    ),
    na.rm =
      TRUE
  )

  ef6_a <- ggplot(
    sub_heat,
    aes(
      x =
        display_output,
      y =
        subgroup_display,
      fill =
        delta_auroc
    )
  ) +
    geom_tile(
      color =
        "white",
      linewidth =
        0.35
    ) +
    geom_text(
      aes(
        label =
          sprintf(
            "%+.02f",
            delta_auroc
          )
      ),
      size =
        2.35,
      family =
        BASE_FONT
    ) +
    scale_fill_gradient2(
      low =
        RED,
      mid =
        "white",
      high =
        BLUE,
      midpoint =
        0,
      limits =
        c(
          -limit6,
          limit6
        ),
      name =
        "\u0394AUROC\nvs overall"
    ) +
    facet_grid(
      subgroup_domain ~ .,
      scales =
        "free_y",
      space =
        "free_y",
      switch =
        "y"
    ) +
    labs(
      title =
        "A. 2024 subgroup discrimination relative to overall performance",
      x =
        NULL,
      y =
        NULL
    ) +
    theme_tract(
      8.2
    ) +
    theme(
      axis.text.x =
        element_text(
          angle =
            35,
          hjust =
            1
        ),
      strip.placement =
        "outside",
      strip.text.y.left =
        element_text(
          angle =
            0,
          hjust =
            1
        )
    )

  context <- subgroup[
    evaluation_year ==
      2024L &
      subgroup_domain %in%
        c(
          "Race",
          "Ethnicity",
          "Payer"
        ) &
      output_id %in%
        c(
          "disposition_postacute",
          "disposition_death_hospice",
          "hlos_ge20"
        ) &
      N >=
        500
  ]

  context[
    ,
    context_level :=
      mapply(
        humanize_context_level,
        subgroup_domain,
        subgroup_level,
        USE.NAMES =
          FALSE
      )
  ]

  context[
    ,
    subgroup_domain :=
      factor(
        subgroup_domain,
        levels =
          c(
            "Race",
            "Ethnicity",
            "Payer"
          )
      )
  ]

  context[
    ,
    display_output :=
      factor(
        unname(
          target_labels[
            output_id
          ]
        ),
        levels =
          unname(
            target_labels[
              c(
                "disposition_postacute",
                "disposition_death_hospice",
                "hlos_ge20"
              )
            ]
          )
      )
  ]

  fwrite(
    context,
    file.path(
      source_dir,
      "eFigure_6_Contextual_Subgroup_Calibration_Data.csv"
    )
  )

  context_limit <- max(
    abs(
      context$calibration_intercept
    ),
    na.rm =
      TRUE
  )

  ef6_b <- ggplot(
    context,
    aes(
      x =
        display_output,
      y =
        context_level,
      fill =
        calibration_intercept
    )
  ) +
    geom_tile(
      color =
        "white",
      linewidth =
        0.35
    ) +
    geom_text(
      aes(
        label =
          sprintf(
            "%.2f",
            calibration_slope
          )
      ),
      size =
        2.35,
      family =
        BASE_FONT
    ) +
    scale_fill_gradient2(
      low =
        RED,
      mid =
        "white",
      high =
        BLUE,
      midpoint =
        0,
      limits =
        c(
          -context_limit,
          context_limit
        ),
      name =
        "Calibration\nintercept"
    ) +
    facet_grid(
      subgroup_domain ~ .,
      scales =
        "free_y",
      space =
        "free_y",
      switch =
        "y"
    ) +
    labs(
      title =
        "B. Contextual subgroup calibration in 2024",
      subtitle =
        "Tile color = calibration intercept; printed value = calibration slope",
      x =
        NULL,
      y =
        NULL
    ) +
    theme_tract(
      8.2
    ) +
    theme(
      axis.text.x =
        element_text(
          angle =
            35,
          hjust =
            1
        ),
      strip.placement =
        "outside",
      strip.text.y.left =
        element_text(
          angle =
            0,
          hjust =
            1
        )
    )

  efigure_6 <-
    ef6_a /
      ef6_b +
      plot_layout(
        heights =
          c(
            1.25,
            0.95
          ),
        guides =
          "keep"
      ) +
      plot_annotation(
        title =
          "eFigure 6. Subgroup robustness and contextual calibration",
        caption =
          "Subgroup analyses are descriptive robustness assessments. Panel B displays all eligible race, ethnicity, and payer groups meeting the prespecified minimum-N rule.",
        theme =
          theme(
            plot.title =
              element_text(
                family =
                  BASE_FONT,
                face =
                  "bold",
                color =
                  NAVY,
                size =
                  14
              ),
            plot.caption =
              element_text(
                family =
                  BASE_FONT,
                color =
                  MID_GRAY,
                size =
                  8,
                hjust =
                  0
              )
          )
      )

  save_figure(
    efigure_6,
    "eFigure_6_Subgroup_and_Context_Calibration",
    width =
      12.0,
    height =
      15.0
  )
}


# ------------------------------------------------------------------------------
# 16. Figure manifest
# ------------------------------------------------------------------------------

manifest <- data.table(
  exhibit = c(
    "Figure 1",
    "Figure 2",
    "Figure 3",
    "eFigure 1",
    "eFigure 2",
    "eFigure 3",
    "eFigure 4",
    "eFigure 5",
    "eFigure 6"
  ),
  file_stem = c(
    "Figure_1_SPINE_TRACT_Framework",
    "Figure_2_Temporal_Performance_Fingerprint",
    "Figure_3_Calibration",
    "eFigure_1_Cohort_Flow",
    "eFigure_2_Case_Mix_Drift_Love_Plot",
    "eFigure_3_Complete_2024_Calibration_Grid",
    "eFigure_4_Duration_Quantile_Calibration_Over_Time",
    "eFigure_5_Predictor_Adjudication_Heatmap",
    "eFigure_6_Subgroup_and_Context_Calibration"
  ),
  purpose = c(
    "Prediction time, temporal design, and final architecture",
    "Temporal AUROC and precision-recall lift",
    "Probability and continuous-duration calibration",
    "Cohort selection for the intended-use population",
    "Development-to-2024 case-mix drift",
    "Calibration of all 2024 probability outputs",
    "Temporal discrete-duration quantile calibration",
    "Predictor-family adjudication and finalist comparison",
    "Spine-specific subgroup robustness and contextual calibration"
  )
)

fwrite(manifest, file.path(figure_dir, "FIGURE_MANIFEST.csv"))

writeLines(
  c(
    "SPINE-TRACT FINAL FIGURE PIPELINE COMPLETE",
    "",
    paste0("Output directory: ", figure_dir),
    "",
    "Main:",
    "  Figure 1  framework / rolling temporal design / architecture",
    "  Figure 2  temporal performance fingerprint",
    "  Figure 3  probability + continuous-duration calibration",
    "",
    "Supplement:",
    "  eFigure 1 complete cohort flow",
    "  eFigure 2 case-mix drift",
    "  eFigure 3 complete 2024 calibration grid",
    "  eFigure 4 duration quantile calibration over time",
    "  eFigure 5 predictor-family adjudication",
    "  eFigure 6 subgroup/context calibration",
    "",
    "All submission TIFFs: 600 dpi, LZW compression.",
    "PNG copies are generated at 300 dpi for review.",
    "No clinical prediction model was fit or recalibrated by this script."
  ),
  file.path(figure_dir, "FIGURE_PIPELINE_SUMMARY.txt")
)

cat(
  "\n============================================================\n",
  "SPINE-TRACT FINAL FIGURES COMPLETE\n",
  "Output: ", figure_dir,
  "\n============================================================\n",
  sep = ""
)
