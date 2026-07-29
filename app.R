# ============================================================
# SPINE-TRACT: Spine Trauma Resource and Acute Care Trajectory Calculator
# Shiny app for shinyapps.io deployment
#
# Required files:
#   data/model_bundle.rds
#   data/predictor_metadata.rds
#
# Optional files:
#   manifest.json
#   www/ohsu_logo.png
# ============================================================

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(data.table)
  library(Matrix)
  library(xgboost)
  library(scales)
  library(stringr)
})

`%||%` <- function(x, y) if (is.null(x)) y else x

# ------------------------------------------------------------
# Load model objects
# ------------------------------------------------------------

bundle_path <- file.path("data", "model_bundle.rds")
metadata_path <- file.path("data", "predictor_metadata.rds")

if (!file.exists(bundle_path)) {
  stop(
    paste0(
      "Could not find data/model_bundle.rds.\n\n",
      "Place the locked SPINE-TRACT model bundle at: ",
      normalizePath(bundle_path, winslash = "/", mustWork = FALSE)
    ),
    call. = FALSE
  )
}

if (!file.exists(metadata_path)) {
  stop(
    paste0(
      "Could not find data/predictor_metadata.rds.\n\n",
      "Place the SPINE-TRACT predictor metadata at: ",
      normalizePath(metadata_path, winslash = "/", mustWork = FALSE)
    ),
    call. = FALSE
  )
}

bundle <- readRDS(bundle_path)
meta_obj <- readRDS(metadata_path)
metadata <- as.data.table(meta_obj$metadata %||% meta_obj)

required_metadata_cols <- c(
  "variable", "label", "group", "input_type", "model_type", "choices",
  "default", "min", "max", "app_order"
)
for (nm in required_metadata_cols) {
  if (!nm %in% names(metadata)) metadata[, (nm) := NA_character_]
}

metadata[, variable := as.character(variable)]
metadata <- metadata[!is.na(variable) & nzchar(variable)]
metadata[, label := fifelse(is.na(label) | !nzchar(label), variable, as.character(label))]
metadata[, group := fifelse(is.na(group) | !nzchar(group), "Other inputs", as.character(group))]
metadata[, input_type := fifelse(is.na(input_type) | !nzchar(input_type), "auto", as.character(input_type))]
metadata[, choices := as.character(choices)]
metadata[, default := as.character(default)]
metadata[, app_order_num := suppressWarnings(as.numeric(app_order))]
metadata[is.na(app_order_num), app_order_num := 9999]

# ------------------------------------------------------------
# Spine-facing labels and grouping polish
# ------------------------------------------------------------

label_map <- c(
  age = "Age, years",
  sex_clean = "Sex",
  race_clean = "Race",
  ethnicity_clean = "Ethnicity",
  insurance_clean = "Insurance",
  transfer_clean = "Interfacility transfer",
  mechanism_clean = "Mechanism of injury",
  gcs_eye_clean = "GCS eye",
  gcs_motor_clean = "GCS motor",
  gcs_verbal_clean = "GCS verbal",
  gcs_total_aug = "Total GCS",
  pupil_clean = "Pupillary response",
  sbp_clean = "Systolic blood pressure, mmHg",
  pulse_clean = "Heart rate, beats/min",
  rr_clean = "Respiratory rate, breaths/min",
  spo2_clean = "Oxygen saturation, %",
  respiratoryassistance_clean = "Respiratory assistance on arrival",
  bleeding_disorder = "Bleeding disorder / coagulopathy",
  diabetes = "Diabetes",
  copd = "COPD",
  hypertension = "Hypertension",
  current_smoker = "Current smoker",
  dx_thoracic_vertebral_fracture = "Thoracic vertebral fracture",
  dx_lumbar_vertebral_fracture = "Lumbar vertebral fracture",
  dx_both_tl_fracture_levels = "Both thoracic and lumbar fracture levels",
  dx_tl_neural_injury = "Thoracolumbar neural injury",
  dx_tl_spinal_cord_injury = "Thoracolumbar spinal cord injury",
  dx_cauda_equina_or_nerve_root = "Cauda equina or nerve-root injury",
  dx_fracture_with_tl_neural_injury = "Fracture with thoracolumbar neural injury",
  dx_cervical_spine_injury = "Concomitant cervical spine injury",
  dx_tbi = "Concomitant TBI",
  dx_intracranial_hemorrhage_any = "Concomitant intracranial hemorrhage",
  dx_thoracic_injury = "Thoracic injury",
  dx_rib_fracture = "Rib fracture",
  dx_multiple_rib_fractures = "Multiple rib fractures",
  dx_flail_chest = "Flail chest",
  dx_pulmonary_contusion_laceration = "Pulmonary contusion/laceration",
  dx_pneumothorax_hemothorax = "Pneumothorax/hemothorax",
  dx_major_thoracic_organ_vascular = "Major thoracic organ/vascular injury",
  dx_abdominal_pelvic_injury = "Abdominal/pelvic injury",
  dx_upper_extremity_injury = "Upper-extremity injury",
  dx_lower_extremity_injury = "Lower-extremity injury",
  operative_spine_intervention = "Operative spine intervention"
)

ui_group_map <- c(
  age = "Demographics",
  sex_clean = "Demographics",
  race_clean = "Demographics",
  ethnicity_clean = "Demographics",
  insurance_clean = "Demographics",
  transfer_clean = "Transfer and mechanism",
  mechanism_clean = "Transfer and mechanism",
  gcs_eye_clean = "Neurologic status",
  gcs_motor_clean = "Neurologic status",
  gcs_verbal_clean = "Neurologic status",
  gcs_total_aug = "Neurologic status",
  pupil_clean = "Neurologic status",
  sbp_clean = "Vital signs and respiratory support",
  pulse_clean = "Vital signs and respiratory support",
  rr_clean = "Vital signs and respiratory support",
  spo2_clean = "Vital signs and respiratory support",
  respiratoryassistance_clean = "Vital signs and respiratory support",
  bleeding_disorder = "Comorbidities",
  diabetes = "Comorbidities",
  copd = "Comorbidities",
  hypertension = "Comorbidities",
  current_smoker = "Comorbidities"
)

metadata[variable %in% names(label_map), label := unname(label_map[variable])]
metadata[variable %in% names(ui_group_map), group := unname(ui_group_map[variable])]
metadata[grepl("thoracic|lumbar|spine|cauda|nerve_root", variable, ignore.case = TRUE), group := "Thoracolumbar spine phenotype"]
metadata[grepl("rib|flail|pulmonary|pneumo|hemo", variable, ignore.case = TRUE), group := "Thoracic and pulmonary injuries"]
metadata[grepl("tbi|intracranial|abdomen|pelvic|extremity|polyregion|fracture_region", variable, ignore.case = TRUE) &
           !grepl("thoracic|lumbar|spine|cauda|nerve_root|rib|pulmonary|pneumo|hemo", variable, ignore.case = TRUE),
         group := "Concomitant injury profile"]
metadata[grepl("preexisting|comorbid|smoker|copd|diabetes|hypertension|bleeding", variable, ignore.case = TRUE), group := "Comorbidities"]

# ------------------------------------------------------------
# Endpoint/model helpers
# ------------------------------------------------------------

clip_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

apply_binary_recalibration <- function(p, obj) {
  cal <- obj$calibration %||% obj$recalibration %||% NULL
  if (is.null(cal)) return(clip_prob(p))
  intercept <- suppressWarnings(as.numeric(cal$intercept %||% cal$alpha %||% NA_real_))
  slope <- suppressWarnings(as.numeric(cal$slope %||% cal$beta %||% NA_real_))
  eps <- suppressWarnings(as.numeric(cal$clip_eps %||% 1e-6))
  if (!is.finite(intercept) || !is.finite(slope)) return(clip_prob(p))
  if (!is.finite(eps) || eps <= 0 || eps >= 0.5) eps <- 1e-6
  as.numeric(plogis(intercept + slope * qlogis(clip_prob(p, eps))))
}

find_xgb <- function(obj, max_depth = 10) {
  if (max_depth < 0 || is.null(obj)) return(NULL)
  if (inherits(obj, "xgb.Booster")) return(obj)
  if (is.list(obj)) {
    for (nm in names(obj)) {
      hit <- find_xgb(obj[[nm]], max_depth - 1)
      if (!is.null(hit)) return(hit)
    }
  }
  NULL
}

get_model <- function(obj) {
  model <- find_xgb(obj)
  if (is.null(model)) stop("Could not find xgb.Booster in endpoint object.", call. = FALSE)
  model
}

pick_endpoint_obj <- function(candidates) {
  nms <- names(bundle)
  for (cand in candidates) {
    if (cand %in% nms) return(bundle[[cand]])
  }
  for (cand in candidates) {
    hit <- nms[tolower(nms) == tolower(cand)][1]
    if (!is.na(hit)) return(bundle[[hit]])
  }
  NULL
}

endpoint_specs <- list(
  discharge = list(
    label = "Discharge disposition",
    type = "multiclass",
    candidates = c("discharge", "discharge_disposition", "discharge_model", "discharge_multiclass")
  ),
  icu_admission = list(
    label = "ICU admission",
    type = "binary",
    candidates = c("icu_admission", "icu", "icu_admit")
  ),
  mechanical_ventilation = list(
    label = "Mechanical ventilation",
    type = "binary",
    candidates = c("mechanical_ventilation", "mech_vent", "ventilation")
  ),
  hlos_ge20 = list(
    label = "Hospital LOS ≥20 days",
    type = "binary",
    candidates = c("hlos_ge20", "hospital_los_ge20", "hospital_los_20")
  ),
  icu_los_ge8 = list(
    label = "ICU LOS ≥8 days",
    type = "binary",
    candidates = c("icu_los_ge8", "icu_los_ge8_full", "icu_los_8")
  ),
  vent_days_ge8 = list(
    label = "Ventilator duration ≥8 days",
    type = "binary",
    candidates = c("vent_days_ge8", "vent_days_ge8_full", "ventilator_days_ge8", "ventilator_duration_ge8")
  ),
  operative_spine_intervention = list(
    label = "Operative spine intervention",
    type = "binary",
    candidates = c("operative_spine_intervention", "spine_operation", "operative_spine", "surgery")
  )
)

available_endpoints <- names(Filter(function(spec) !is.null(pick_endpoint_obj(spec$candidates)), endpoint_specs))
if (length(available_endpoints) == 0) {
  stop("No recognized SPINE-TRACT endpoint objects were found in model_bundle.rds.", call. = FALSE)
}

get_endpoint_obj <- function(endpoint) {
  obj <- pick_endpoint_obj(endpoint_specs[[endpoint]]$candidates)
  if (is.null(obj)) stop("Endpoint not found in model bundle: ", endpoint, call. = FALSE)
  obj
}

get_predictors <- function(endpoint, obj) {
  preds <- obj$predictors %||% obj$raw_predictors %||% NULL
  if (!is.null(preds) && length(preds) > 0) return(as.character(preds))
  if (endpoint == "mechanical_ventilation") {
    preds <- meta_obj$predictors_no_resp %||% meta_obj$predictors_full %||% NULL
  } else {
    preds <- meta_obj$predictors_full %||% NULL
  }
  if (!is.null(preds) && length(preds) > 0) return(as.character(preds))
  unique(metadata$variable)
}

get_feature_names <- function(obj, model) {
  fn <- obj$feature_names %||% model$feature_names %||% NULL
  if (!is.null(fn) && length(fn) > 0) return(as.character(fn))
  NULL
}

# ------------------------------------------------------------
# Input/default helpers
# ------------------------------------------------------------

input_id <- function(v) paste0("var__", v)

split_choices <- function(x, fallback = character(0)) {
  if (is.null(x) || length(x) == 0 || is.na(x[1]) || !nzchar(as.character(x[1]))) return(fallback)
  out <- unlist(strsplit(as.character(x[1]), "\\|\\|"))
  out <- out[!is.na(out) & nzchar(out)]
  if (length(out) == 0) fallback else out
}

safe_numeric <- function(x, default = 0) {
  if (is.null(x) || length(x) == 0) return(default)
  y <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(y)) return(default)
  y
}

safe_yesno <- function(x) {
  if (is.null(x) || length(x) == 0) return(0L)
  as.integer(as.character(x[[1]]) %in% c("1", "Yes", "yes", "TRUE", "true", TRUE))
}

row_for_var <- function(v) {
  r <- metadata[variable == v][1]
  if (nrow(r) == 0) {
    data.table(
      variable = v,
      label = ifelse(v %in% names(label_map), unname(label_map[v]), v),
      group = "Other inputs",
      input_type = "auto",
      choices = NA_character_,
      default = NA_character_,
      min = NA_character_,
      max = NA_character_,
      app_order_num = 9999
    )
  } else {
    r
  }
}

should_show_input <- function(v) {
  r <- row_for_var(v)
  it <- tolower(as.character(r$input_type %||% "auto"))
  !it %in% c("derived", "hidden", "exclude", "none")
}

make_default_value <- function(r) {
  if (!is.na(r$default) && nzchar(r$default)) return(as.character(r$default))
  choices <- split_choices(r$choices)
  if (length(choices) > 0) return(choices[1])
  it <- tolower(as.character(r$input_type))
  if (grepl("age", r$variable)) return("65")
  if (grepl("gcs_total", r$variable)) return("15")
  if (grepl("sbp", r$variable)) return("120")
  if (grepl("pulse", r$variable)) return("80")
  if (grepl("rr", r$variable)) return("16")
  if (grepl("spo2", r$variable)) return("98")
  if (it %in% c("yesno", "binary", "logical")) return("0")
  "0"
}

make_input_control <- function(r) {
  v <- r$variable
  lab <- r$label
  it <- tolower(as.character(r$input_type))
  choices <- split_choices(r$choices)
  default <- make_default_value(r)

  if (length(choices) > 0 && !(it %in% c("numeric", "integer", "number"))) {
    if (all(choices %in% c("0", "1")) || it %in% c("yesno", "binary", "logical")) {
      return(selectInput(input_id(v), lab, choices = c("No" = "0", "Yes" = "1"), selected = default))
    }
    return(selectInput(input_id(v), lab, choices = stats::setNames(choices, choices), selected = default))
  }

  if (it %in% c("yesno", "binary", "logical") || grepl("^dx_|fracture|injury|copd|diabetes|smoker|hypertension|bleeding", v, ignore.case = TRUE)) {
    return(selectInput(input_id(v), lab, choices = c("No" = "0", "Yes" = "1"), selected = default))
  }

  min_val <- suppressWarnings(as.numeric(r$min))
  max_val <- suppressWarnings(as.numeric(r$max))
  if (!is.finite(min_val)) min_val <- NA_real_
  if (!is.finite(max_val)) max_val <- NA_real_
  value <- suppressWarnings(as.numeric(default))
  if (!is.finite(value)) value <- 0

  numericInput(input_id(v), lab, value = value, min = min_val, max = max_val)
}

get_input_value <- function(input, v) {
  r <- row_for_var(v)
  raw <- input[[input_id(v)]]
  if (is.null(raw)) raw <- make_default_value(r)
  it <- tolower(as.character(r$input_type))
  choices <- split_choices(r$choices)

  if (it %in% c("numeric", "integer", "number") || length(choices) == 0 && suppressWarnings(!is.na(as.numeric(raw)))) {
    return(safe_numeric(raw, safe_numeric(make_default_value(r), 0)))
  }
  if (all(choices %in% c("0", "1")) || it %in% c("yesno", "binary", "logical")) {
    return(safe_yesno(raw))
  }
  as.character(raw[[1]])
}

first_choice <- function(v, choices, fallback) {
  hit <- choices[grepl(v, choices, ignore.case = TRUE)][1]
  ifelse(is.na(hit), fallback, hit)
}

# ------------------------------------------------------------
# Derived variable handling
# ------------------------------------------------------------

derive_variable <- function(v, input, r) {
  choices <- split_choices(r$choices, c("0", "1"))

  age <- safe_numeric(input[[input_id("age")]], 65)
  gcs <- safe_numeric(input[[input_id("gcs_total_aug")]], 15)
  sbp <- safe_numeric(input[[input_id("sbp_clean")]], 120)
  pulse <- safe_numeric(input[[input_id("pulse_clean")]], 80)
  rr <- safe_numeric(input[[input_id("rr_clean")]], 16)
  spo2 <- safe_numeric(input[[input_id("spo2_clean")]], 98)

  if (v %in% c("age_group_aug", "age_group")) {
    if (age < 40) return(first_choice("18|39|young", choices, choices[1]))
    if (age < 65) return(first_choice("40|64|adult", choices, choices[1]))
    if (age < 75) return(first_choice("65|74", choices, choices[min(length(choices), 2)]))
    return(first_choice("75|80|older|elder", choices, choices[length(choices)]))
  }

  if (v %in% c("gcs_severity_aug", "gcs_severity")) {
    if (gcs >= 13) return(first_choice("mild|13", choices, choices[1]))
    if (gcs >= 9) return(first_choice("moderate|9", choices, choices[1]))
    return(first_choice("severe|3", choices, choices[length(choices)]))
  }

  if (v %in% c("hypotension_sbp90_aug", "hypotension_sbp90")) return(as.integer(sbp < 90))
  if (v %in% c("hypoxia_spo2_90_aug", "hypoxia_spo2_90")) return(as.integer(spo2 <= 90))
  if (v %in% c("tachycardia_120_aug", "tachycardia_120")) return(as.integer(pulse >= 120))
  if (v %in% c("abnormal_rr_aug", "abnormal_rr")) return(as.integer(rr < 10 | rr > 29))

  yesno_names <- metadata[tolower(input_type) %in% c("yesno", "binary", "logical"), variable]
  comorb_names <- yesno_names[grepl("bleeding|diabetes|copd|hypertension|smoker|comorbid|preexisting", yesno_names, ignore.case = TRUE)]
  if (v %in% c("n_preexisting_conditions", "preexisting_condition_count")) {
    return(sum(vapply(comorb_names, function(z) safe_yesno(input[[input_id(z)]]), integer(1)), na.rm = TRUE))
  }

  dx_names <- metadata[grepl("^dx_", variable), variable]
  if (v %in% c("n_unique_dx_codes", "n_traumatic_dx_codes", "n_total_dx_codes")) {
    return(sum(vapply(dx_names, function(z) safe_yesno(input[[input_id(z)]]), integer(1)), na.rm = TRUE))
  }

  if (v %in% c("dx_tl_fracture_any", "tl_fracture_any")) {
    return(as.integer(safe_yesno(input[[input_id("dx_thoracic_vertebral_fracture")]]) == 1L || safe_yesno(input[[input_id("dx_lumbar_vertebral_fracture")]]) == 1L))
  }

  if (v %in% c("dx_both_tl_fracture_levels", "both_tl_fracture_levels")) {
    return(as.integer(safe_yesno(input[[input_id("dx_thoracic_vertebral_fracture")]]) == 1L && safe_yesno(input[[input_id("dx_lumbar_vertebral_fracture")]]) == 1L))
  }

  if (v %in% c("dx_fracture_with_tl_neural_injury", "fracture_with_tl_neural_injury")) {
    fx <- safe_yesno(input[[input_id("dx_thoracic_vertebral_fracture")]]) == 1L || safe_yesno(input[[input_id("dx_lumbar_vertebral_fracture")]]) == 1L
    neural <- safe_yesno(input[[input_id("dx_tl_neural_injury")]]) == 1L || safe_yesno(input[[input_id("dx_tl_spinal_cord_injury")]]) == 1L || safe_yesno(input[[input_id("dx_cauda_equina_or_nerve_root")]]) == 1L
    return(as.integer(fx && neural))
  }

  default <- make_default_value(r)
  if (all(choices %in% c("0", "1"))) return(safe_yesno(default))
  if (suppressWarnings(!is.na(as.numeric(default)))) return(safe_numeric(default, 0))
  as.character(default)
}

value_for_predictor <- function(input, v) {
  r <- row_for_var(v)
  it <- tolower(as.character(r$input_type))
  if (it %in% c("derived", "hidden", "exclude", "none") || !should_show_input(v)) {
    return(derive_variable(v, input, r))
  }
  get_input_value(input, v)
}

make_prediction_row <- function(input, predictors) {
  vals <- lapply(predictors, function(v) value_for_predictor(input, v))
  names(vals) <- predictors
  as.data.frame(vals, stringsAsFactors = FALSE, check.names = FALSE)
}

coerce_prediction_row <- function(df, predictors) {
  for (v in predictors) {
    r <- row_for_var(v)
    it <- tolower(as.character(r$input_type))
    choices <- split_choices(r$choices)
    if (it %in% c("numeric", "integer", "number") || suppressWarnings(!is.na(as.numeric(df[[v]][1])) && length(choices) == 0)) {
      df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
    } else if (all(choices %in% c("0", "1")) || it %in% c("yesno", "binary", "logical")) {
      df[[v]] <- suppressWarnings(as.numeric(df[[v]]))
    } else if (length(choices) > 0) {
      df[[v]] <- factor(as.character(df[[v]]), levels = choices)
    } else {
      df[[v]] <- as.character(df[[v]])
    }
  }
  df
}

make_design_matrix <- function(df, predictors, feature_names = NULL) {
  df <- coerce_prediction_row(df, predictors)
  mm <- model.matrix(stats::as.formula(paste("~", paste(sprintf("`%s`", predictors), collapse = " + "), "- 1")), data = df)
  if (!is.null(feature_names)) {
    aligned <- matrix(0, nrow = nrow(mm), ncol = length(feature_names))
    colnames(aligned) <- feature_names
    common <- intersect(colnames(mm), feature_names)
    if (length(common) > 0) aligned[, common] <- mm[, common, drop = FALSE]
    mm <- aligned
  }
  Matrix::Matrix(mm, sparse = TRUE)
}

predict_endpoint <- function(endpoint, input) {
  spec <- endpoint_specs[[endpoint]]
  obj <- get_endpoint_obj(endpoint)
  model <- get_model(obj)
  predictors <- get_predictors(endpoint, obj)
  predictors <- predictors[predictors %in% unique(c(metadata$variable, predictors))]
  df <- make_prediction_row(input, predictors)
  feature_names <- get_feature_names(obj, model)
  x <- make_design_matrix(df, predictors, feature_names)
  dmat <- xgb.DMatrix(x)
  raw <- as.numeric(predict(model, dmat))

  if (spec$type == "multiclass") {
    levels <- obj$levels %||% obj$outcome_levels %||% meta_obj$discharge_levels %||% c("Home/home health", "Post-acute facility", "Death/hospice")
    if (length(raw) != length(levels)) {
      levels <- paste0("Class ", seq_along(raw))
    }
    p <- raw / sum(raw)
    return(data.table(endpoint = spec$label, class = levels, probability = p, type = "multiclass"))
  }

  p <- apply_binary_recalibration(raw[1], obj)
  data.table(endpoint = spec$label, class = spec$label, probability = p, type = "binary")
}

# ------------------------------------------------------------
# UI helpers
# ------------------------------------------------------------

endpoint_choices <- stats::setNames(
  available_endpoints,
  vapply(available_endpoints, function(e) endpoint_specs[[e]]$label, character(1))
)

all_predictors <- unique(unlist(lapply(available_endpoints, function(e) get_predictors(e, get_endpoint_obj(e)))))
visible_vars <- unique(metadata[variable %in% all_predictors & !tolower(input_type) %in% c("derived", "hidden", "exclude", "none")][order(app_order_num), variable])
if (length(visible_vars) == 0) visible_vars <- unique(metadata[!tolower(input_type) %in% c("derived", "hidden", "exclude", "none")][order(app_order_num), variable])

logo_ui <- if (file.exists(file.path("www", "ohsu_logo.png"))) {
  img(src = "ohsu_logo.png", class = "ohsu-logo")
} else {
  div("OHSU", class = "ohsu-logo-fallback")
}

risk_band <- function(p) {
  ifelse(p < 0.10, "Lower", ifelse(p < 0.30, "Moderate", "Higher"))
}

ui <- page_fluid(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    primary = "#1f4e79",
    bg = "#f4f7fb",
    fg = "#243447"
  ),
  tags$head(
    tags$meta(name = "viewport", content = "width=device-width, initial-scale=1"),
    tags$style(HTML("\n      :root { --page-max: 1340px; --card-radius: 24px; --shadow-soft: 0 8px 28px rgba(31, 52, 73, 0.07); --border-soft: #e7edf5; --text-main: #243447; --text-muted: #5b6b7f; --accent: #1f4e79; }\n      body { background: #f4f7fb; }\n      .app-container { max-width: var(--page-max); margin: 0 auto; padding: 24px 22px 36px 22px; }\n      .app-header { background: #ffffff; border-radius: 28px; padding: clamp(18px, 2.2vw, 30px); margin-bottom: 24px; box-shadow: var(--shadow-soft); border: 1px solid var(--border-soft); }\n      .header-grid { display: grid; grid-template-columns: minmax(70px, 96px) 1fr; gap: 20px; align-items: center; }\n      .logo-wrap { display: flex; align-items: center; justify-content: center; }\n      .ohsu-logo { width: clamp(58px, 6vw, 92px); height: auto; display: block; }\n      .ohsu-logo-fallback { width: 88px; height: 88px; border-radius: 22px; display: flex; align-items: center; justify-content: center; background: #1f4e79; color: white; font-weight: 900; letter-spacing: 0.06em; }\n      .header-title { margin: 0 0 8px 0; font-weight: 800; line-height: 1.04; font-size: clamp(1.9rem, 3.4vw, 3.2rem); color: var(--text-main); max-width: 1000px; }\n      .header-subtitle { color: var(--text-muted); margin: 0 0 3px 0; font-size: 1.05rem; }\n      .input-card, .metric-card, .plot-card, .detail-card { background: #ffffff; border: 1px solid var(--border-soft) !important; border-radius: var(--card-radius) !important; box-shadow: var(--shadow-soft); }\n      .input-card .card-body { padding: 18px; }\n      .metric-card .card-body, .plot-card .card-body, .detail-card .card-body { padding: 22px; }\n      .sticky-panel { position: sticky; top: 24px; max-height: calc(100vh - 48px); overflow-y: auto; padding-right: 4px; scrollbar-width: thin; }\n      .section-title { font-weight: 800; color: var(--text-main); margin-bottom: 14px; line-height: 1.06; font-size: clamp(1.55rem, 2vw, 2rem); }\n      .form-label { font-weight: 650; color: #2f4257; margin-bottom: 5px; font-size: 0.97rem; }\n      .shiny-input-container { margin-bottom: 10px; }\n      .form-control, .form-select { border-radius: 14px !important; border: 1px solid #d4dde8 !important; min-height: 44px; box-shadow: none !important; }\n      .btn-primary { background-color: #245789 !important; border-color: #245789 !important; border-radius: 14px !important; font-weight: 750; min-height: 46px; margin-top: 6px; }\n      .metric-grid { display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 14px; margin-bottom: 18px; }\n      .metric-value { font-size: clamp(1.5rem, 2.2vw, 2.15rem); line-height: 1; font-weight: 800; color: var(--accent); margin-bottom: 10px; }\n      .metric-label { font-size: 0.92rem; color: var(--text-muted); line-height: 1.35; }\n      .note { color: #5b6b7f; font-size: 0.94rem; line-height: 1.45; }\n      .group-title { font-weight: 800; margin-top: 12px; margin-bottom: 8px; color: #243447; border-bottom: 1px solid #e7edf5; padding-bottom: 6px; }\n      @media (max-width: 1199px) { .metric-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); } .sticky-panel { position: static; max-height: none; overflow-y: visible; padding-right: 0; } }\n      @media (max-width: 767px) { .app-container { padding: 18px 14px 28px 14px; } .header-grid, .metric-grid { grid-template-columns: 1fr; } .header-grid { text-align: center; } }\n    "))
  ),
  div(
    class = "app-container",
    div(
      class = "app-header",
      div(
        class = "header-grid",
        div(class = "logo-wrap", logo_ui),
        div(
          h1("SPINE-TRACT", class = "header-title"),
          p("Spine Trauma Resource and Acute Care Trajectory Calculator", class = "header-subtitle"),
          p("Admission-era estimates after acute thoracolumbar spine trauma", class = "header-subtitle")
        )
      )
    ),
    layout_columns(
      col_widths = c(4, 8),
      card(
        class = "input-card sticky-panel",
        card_body(
          div(class = "section-title", "Patient inputs"),
          selectInput("endpoint", "Endpoint", choices = endpoint_choices, selected = available_endpoints[1]),
          uiOutput("dynamic_inputs"),
          actionButton("calculate", "Calculate", class = "btn-primary w-100")
        )
      ),
      div(
        card(
          class = "metric-card",
          card_body(
            div(class = "section-title", "Estimated trajectory"),
            uiOutput("prediction_cards")
          )
        ),
        div(style = "height: 18px;"),
        card(
          class = "detail-card",
          card_body(
            h4("Interpretation"),
            p(class = "note", "SPINE-TRACT estimates observed resource-utilization trajectories. Outputs should support early resource planning and structured communication, not determine whether ICU admission, mechanical ventilation, prolonged hospitalization, or operative intervention is clinically indicated."),
            p(class = "note", paste0("Model years: ", paste(meta_obj$model_years %||% "2020–2024", collapse = ", "), ". Development cohort N: ", format(meta_obj$model_n %||% NA_integer_, big.mark = ","), "."))
          )
        )
      )
    )
  )
)

server <- function(input, output, session) {
  output$dynamic_inputs <- renderUI({
    endpoint <- input$endpoint %||% available_endpoints[1]
    obj <- get_endpoint_obj(endpoint)
    predictors <- get_predictors(endpoint, obj)
    vars <- unique(metadata[variable %in% predictors & variable %in% visible_vars][order(app_order_num), variable])
    if (length(vars) == 0) vars <- visible_vars

    md <- rbindlist(lapply(vars, row_for_var), fill = TRUE)
    md <- md[order(app_order_num, group, label)]
    groups <- unique(md$group)

    tagList(lapply(groups, function(g) {
      rows <- md[group == g]
      tagList(
        div(class = "group-title", g),
        lapply(seq_len(nrow(rows)), function(i) make_input_control(rows[i]))
      )
    }))
  })

  prediction <- eventReactive(input$calculate, {
    req(input$endpoint)
    predict_endpoint(input$endpoint, input)
  }, ignoreInit = FALSE)

  output$prediction_cards <- renderUI({
    res <- prediction()
    req(res)

    if (unique(res$type) == "multiclass") {
      res <- res[order(-probability)]
      div(
        class = "metric-grid",
        lapply(seq_len(nrow(res)), function(i) {
          div(
            class = "metric-card",
            div(class = "card-body",
                div(class = "metric-value", scales::percent(res$probability[i], accuracy = 0.1)),
                div(class = "metric-label", res$class[i]))
          )
        })
      )
    } else {
      p <- res$probability[1]
      div(
        class = "metric-grid",
        div(class = "metric-card", div(class = "card-body", div(class = "metric-value", scales::percent(p, accuracy = 0.1)), div(class = "metric-label", res$endpoint[1]))),
        div(class = "metric-card", div(class = "card-body", div(class = "metric-value", risk_band(p)), div(class = "metric-label", "Descriptive risk band")))
      )
    }
  })
}

shinyApp(ui, server)
