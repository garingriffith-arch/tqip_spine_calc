# SPINE-TRACT | Trauma Resource and Acute Care Trajectory Calculator
# Source specification: 2020–2024 repaired-DX deployment, Stage 21.
# UI mirrors TBI-TRACT. This application does not provide independent validation.

suppressPackageStartupMessages({
  library(shiny)
  library(bslib)
  library(data.table)
  library(Matrix)
  library(xgboost)
  library(reticulate)
})

`%||%` <- function(x, y) if (is.null(x) || !length(x)) y else x

artifact_root <- Sys.getenv("SPINE_TRACT_APP_ARTIFACT_DIR", unset =
  file.path("data", "deployment_2020_2024"))
model_dir <- file.path(artifact_root, "models")
encoder_dir <- file.path(artifact_root, "encoders")
model_ids <- c(
  "discharge3_final", "hlos3_final", "icu3_final", "vent3_final",
  "operative_spine_intervention", "hospital_los",
  "icu_los_conditional", "ventilator_days_conditional"
)
duration_ids <- c("hospital_los", "icu_los_conditional", "ventilator_days_conditional")
context_ids <- c("discharge3_final", "hlos3_final")

model_paths <- setNames(file.path(model_dir, paste0(model_ids, ".json")), model_ids)
encoder_paths <- setNames(file.path(encoder_dir,
                                   paste0(model_ids, "_encoder.rds")), model_ids)
all_paths <- c(model_paths, encoder_paths)
if (any(!file.exists(all_paths))) stop(
  "Missing deployment artifacts (do not substitute the legacy model bundle):\n",
  paste(all_paths[!file.exists(all_paths)], collapse = "\n"), call. = FALSE
)

read_entry <- function(id) {
  md <- readRDS(encoder_paths[[id]])
  required <- c("endpoint_id", "type", "predictors", "encoder", "feature_names",
                "predictor_policy", "training_n", "training_years", "selected_rounds")
  if (length(setdiff(required, names(md)))) stop(
    id, " metadata missing: ", paste(setdiff(required, names(md)), collapse = ", "))
  if (!identical(as.character(md$endpoint_id), id)) stop(id, " endpoint metadata mismatch")
  if (!identical(as.integer(md$training_years), 2020:2024)) stop(id, " training-year mismatch")
  want_policy <- if (id %in% context_ids) "final_context" else "final_core"
  if (!identical(md$predictor_policy, want_policy)) stop(id, " predictor policy mismatch")
  expected <- if (want_policy == "final_context") 57L else 54L
  if (length(md$predictors) != expected) stop(id, " expected ", expected, " predictors")
  if (length(md$feature_names) < expected || anyDuplicated(md$feature_names)) {
    stop(id, " invalid encoded feature names")
  }
  if (!all(md$predictors %in% names(md$encoder$per_predictor))) {
    stop(id, " missing per-predictor encoding specifications")
  }
  if (id %in% duration_ids && !identical(md$type, "joint_quantile")) {
    stop(id, " should have joint_quantile objective")
  }
  if (id %in% duration_ids && !identical(as.numeric(md$quantile_alphas), c(.1,.5,.9))) {
    stop(id, " unexpected quantile alphas")
  }
  md
}
metadata <- setNames(lapply(model_ids, read_entry), model_ids)

# R XGBoost classifiers: trained with sparse Matrix::Matrix inputs.
class_models <- setNames(lapply(setdiff(model_ids, duration_ids), function(id) {
  xgboost::xgb.load(model_paths[[id]])
}), setdiff(model_ids, duration_ids))

# Python XGBoost v3 joint-quantile boosters. Load once per R process.
# Posit Connect Cloud provisions Python from requirements.txt/manifest.json.
# Allow reticulate to use RETICULATE_PYTHON (if set by host) or its
# environment discovery. Never select a workstation-specific Conda env here.
py_xgb <- tryCatch(reticulate::import("xgboost", convert = FALSE),
                   error = function(e) stop(
                     "Python XGBoost is required for three joint-quantile models. ",
                     "Install xgboost>=2 in the reticulate Python environment. Details: ",
                     conditionMessage(e), call. = FALSE))
py_ver <- as.character(reticulate::py_to_r(py_xgb$`__version__`))
if (as.integer(strsplit(py_ver, "\\.")[[1L]][1L]) < 2L) {
  stop("Python XGBoost >=2 required for joint-quantile inference; found ", py_ver)
}
duration_models <- setNames(lapply(duration_ids, function(id) {
  py_xgb$Booster(model_file = normalizePath(model_paths[[id]], winslash = "/"))
}), duration_ids)

used <- unique(unlist(lapply(metadata, `[[`, "predictors"), use.names = FALSE))
if (length(used) != 57L) stop("Unexpected total SPINE-TRACT predictor count: ", length(used))
context_vars <- setdiff(metadata$discharge3_final$predictors,
                        metadata$icu3_final$predictors)
if (!setequal(context_vars, c("race_clean", "ethnicity_clean", "payer_clean"))) {
  stop("Unexpected social/contextual predictor set: ", paste(context_vars, collapse = ", "))
}
for (id in setdiff(model_ids, context_ids)) {
  if (length(intersect(context_vars, metadata[[id]]$predictors))) {
    stop("Social/contextual variables leaked into clinical-only endpoint: ", id)
  }
}

is_numeric <- function(v) {
  any(vapply(metadata, function(m) {
    identical(m$encoder$per_predictor[[v]]$type, "numeric")
  }, logical(1)))
}
numeric_predictors <- used[vapply(used, is_numeric, logical(1))]
categorical_predictors <- setdiff(used, numeric_predictors)
binary_predictors <- intersect(used, unique(c(
  grep("^(dx_|pmhx_)", used, value = TRUE),
  grep("(gcsq_|^comorbidity_)", used, value = TRUE)
)))
# Only known numeric binary flags get Yes/No/Unknown menus.
binary_predictors <- intersect(binary_predictors, numeric_predictors)
continuous_predictors <- setdiff(numeric_predictors, binary_predictors)
get_levels <- function(v) {
  for (m in metadata) {
    d <- m$encoder$per_predictor[[v]]
    if (!is.null(d) && identical(d$type, "categorical")) return(as.character(d$levels))
  }
  character()
}

known_labels <- c(
  age_years = "Age, years", sex_clean = "Sex", race_clean = "Race",
  ethnicity_clean = "Ethnicity", payer_clean = "Payer / insurance",
  mechanism_clean = "Mechanism of injury", gcs_eye_clean = "GCS eye",
  gcs_verbal_clean = "GCS verbal", gcs_motor_clean = "GCS motor",
  gcsq_intubated_recovered = "Intubated during GCS assessment",
  gcsq_sedated_paralyzed_recovered = "Sedated / paralyzed during GCS assessment",
  gcsq_eye_obstruction_recovered = "GCS eye score affected by obstruction",
  gcsq_unknown_recovered = "GCS assessment unknown",
  pupil_clean = "Pupillary response", sbp_clean = "Systolic blood pressure, mm Hg",
  pulse_clean = "Heart rate, beats/min", rr_clean = "Respiratory rate, breaths/min",
  spo2_clean = "Oxygen saturation, %", temperature_c_recovered = "Temperature, °C",
  respiratory_assistance_clean = "Respiratory assistance at initial assessment",
  dx_thoracic_vertebral_fracture = "Thoracic vertebral fracture",
  dx_lumbar_vertebral_fracture = "Lumbar vertebral fracture",
  dx_tl_spinal_cord_injury_any = "Thoracolumbar spinal cord injury",
  dx_tl_root_cauda_any = "Thoracolumbar nerve root / cauda equina injury",
  dx_tl_soft_tissue_any = "Thoracolumbar soft-tissue injury",
  dx_cervical_spine_injury = "Concomitant cervical spine injury",
  dx_tbi_any = "Concomitant traumatic brain injury",
  dx_rib_fracture = "Rib fracture", dx_multiple_rib_fractures = "Multiple rib fractures",
  dx_flail_chest = "Flail chest", dx_pulmonary_contusion_laceration = "Pulmonary contusion / laceration",
  dx_pneumothorax_hemothorax = "Pneumothorax / hemothorax",
  dx_abdominal_pelvic_injury = "Abdominal / pelvic injury",
  dx_upper_extremity_injury = "Upper-extremity injury",
  dx_lower_extremity_injury = "Lower-extremity injury"
)
label_for <- function(v) {
  if (v %in% names(known_labels)) return(unname(known_labels[[v]]))
  z <- gsub("^(pmhx_|dx_)", "", v)
  z <- gsub("_clean$|_recovered$|_any$", "", z)
  z <- gsub("_", " ", z)
  paste0(toupper(substring(z, 1L, 1L)), substring(z, 2L))
}

# Some category IDs are registry-specific. Preserve the true encoded values;
# only use an explicit verified display map for unambiguous categories.
pretty_level <- function(v, z) {
  z <- trimws(as.character(z))
  if (z %in% c("__UNKNOWN__", "Unknown", "NA", "")) return("Unknown / not recorded")
  if (z == "__OTHER__") return("Other / not listed")
  if (v == "sex_clean") {
    if (z == "SEX_1") return("Male")
    if (z == "SEX_2") return("Female")
    if (z == "SEX_3") return("Nonbinary / other (registry category)")
  }
  if (v == "ethnicity_clean") {
    if (z == "ETH_1") return("Hispanic or Latino")
    if (z == "ETH_2") return("Not Hispanic or Latino")
  }
  if (v == "payer_clean") {
    # Verified against the previously reviewed SPINE-TRACT eTable 11 category definitions.
    payer_labels <- c(PAYER_1="Medicaid", PAYER_3="Self-pay",
      PAYER_4="Private / commercial insurance", PAYER_6="Medicare",
      PAYER_7="Other government", PAYER_10="Other insurance")
    if (z %in% names(payer_labels)) return(unname(payer_labels[[z]]))
    if (grepl("^PAYER_", z)) return(paste0("Registry payer code ", sub("^PAYER_", "", z)))
  }
  if (grepl("^(SEX|ETH)_[[:digit:]]+$", z)) return(z)
  if (grepl("^Transport/MVC$", z)) return("Motor vehicle / transport-related")
  gsub("_", " ", z, fixed = TRUE)
}

choices_for <- function(v) {
  lev <- unique(get_levels(v)); if (!length(lev)) lev <- "__UNKNOWN__"
  lev <- setdiff(lev, "__OTHER__") # encoder fallback, not a clinical answer
  labels <- vapply(lev, function(z) pretty_level(v, z), character(1))
  keep <- !duplicated(labels)
  lev <- lev[keep]; labels <- labels[keep]
  unk <- grepl("Unknown|not recorded", labels, ignore.case = TRUE)
  ord <- order(unk, tolower(labels))
  stats::setNames(lev[ord], labels[ord])
}
input_id <- function(v) paste0("var__", v)

numeric_limits <- list(
  age_years = list(min=18, max=89, step=1, note="Study population: age 18–89 years"),
  sbp_clean = list(min=0, max=300, step=1, note="Interface range: 0–300 mm Hg"),
  pulse_clean = list(min=0, max=300, step=1, note="Interface range: 0–300 beats/min"),
  rr_clean = list(min=0, max=100, step=1, note="Interface range: 0–100 breaths/min"),
  spo2_clean = list(min=0, max=100, step=1, note="Interface range: 0–100%"),
  temperature_c_recovered = list(min=25, max=45, step=.1, note="Interface range: 25–45 °C")
)
gcs_choices <- list(
  gcs_eye_clean = c("1 - None"="1", "2 - To pain"="2", "3 - To speech"="3",
                    "4 - Spontaneous"="4", "Unknown / not recorded"=""),
  gcs_verbal_clean = c("1 - None"="1", "2 - Incomprehensible sounds"="2",
                       "3 - Inappropriate words"="3", "4 - Confused"="4",
                       "5 - Oriented"="5", "Unknown / not recorded"=""),
  gcs_motor_clean = c("1 - None"="1", "2 - Extension"="2", "3 - Flexion"="3",
                      "4 - Withdraws"="4", "5 - Localizes"="5",
                      "6 - Obeys commands"="6", "Unknown / not recorded"="")
)

# Groups are generated from the actual locked metadata, not legacy UI metadata.
choose_group <- function(v) {
  if (v %in% c("age_years","sex_clean","race_clean","ethnicity_clean","payer_clean"))
    return("Demographics and context")
  if (grepl("mechanism|gcs|pupil", v, ignore.case = TRUE))
    return("Mechanism and neurologic status")
  if (grepl("sbp|pulse|respir|spo2|oxygen|temperature|heart_rate|rr_clean", v, ignore.case=TRUE))
    return("Admission physiology and respiratory support")
  if (grepl("^dx_", v) && grepl("thoracic|lumbar|tl_|cauda|nerve_root|spine|vertebr|cord", v))
    return("Thoracolumbar spine phenotype")
  if (grepl("^dx_", v) && grepl("rib|flail|chest|pulmon|pneum|hemothorax", v))
    return("Thoracic and pulmonary injuries")
  if (grepl("^dx_", v)) return("Concomitant injury pattern")
  if (grepl("^pmhx_|diabetes|copd|smoker|hypertension|renal|chf|cirrhosis",v))
    return("Pre-existing conditions")
  "Other initial evaluation inputs"
}
group_order <- c("Demographics and context", "Mechanism and neurologic status",
                 "Admission physiology and respiratory support", "Thoracolumbar spine phenotype",
                 "Thoracic and pulmonary injuries", "Concomitant injury pattern",
                 "Pre-existing conditions", "Other initial evaluation inputs")
groups <- split(used, factor(vapply(used, choose_group, character(1)), levels=group_order))
groups <- groups[lengths(groups)>0]

# A preview is an explicitly reduced-input analysis, not separately validated.
quick_common <- intersect(used, c(
  "age_years", "sex_clean", "gcs_eye_clean", "gcs_verbal_clean", "gcs_motor_clean",
  "gcsq_intubated_recovered", "gcsq_sedated_paralyzed_recovered", "pupil_clean",
  "sbp_clean", "pulse_clean", "rr_clean", "spo2_clean", "mechanism_clean",
  "respiratory_assistance_clean", "dx_thoracic_vertebral_fracture",
  "dx_lumbar_vertebral_fracture", "dx_tl_spinal_cord_injury_any",
  "dx_tl_root_cauda_any", "dx_tl_soft_tissue_any", "dx_cervical_spine_injury",
  "dx_tbi_any", "pmhx_diabetes", "pmhx_anticoagulant_therapy"
))
quick_sets <- list(
  disposition = intersect(unique(c(quick_common, context_vars, "pmhx_dementia", "pmhx_chf")), used),
  hlos = intersect(unique(c(quick_common, context_vars, "dx_pulmonary_contusion_laceration")), used),
  icu = intersect(unique(c(quick_common, "dx_rib_fracture", "dx_pneumothorax_hemothorax")), used),
  ventilation = intersect(unique(c(quick_common, "dx_pulmonary_contusion_laceration")), used),
  operation = intersect(unique(c(quick_common, "dx_tl_spinal_cord_injury_any", "dx_tl_root_cauda_any")), used)
)
quick_choices <- c("Discharge disposition"="disposition",
                   "Hospital length of stay"="hlos", "ICU trajectory"="icu",
                   "Mechanical ventilation trajectory"="ventilation",
                   "Recorded thoracolumbar operation"="operation")

make_control <- function(v, cached = NULL) {
  id <- input_id(v); lab <- label_for(v)
  cached_chr <- if (length(cached)) as.character(cached[[1L]]) else NA_character_
  if (v %in% names(gcs_choices)) {
    vals <- unname(gcs_choices[[v]])
    selected <- if (!is.na(cached_chr) && cached_chr %in% vals) cached_chr else ""
    return(selectInput(id,lab,gcs_choices[[v]],selected=selected,selectize=FALSE))
  }
  if (v %in% binary_predictors) {
    opts <- c("No"="0", "Yes"="1", "Unknown / not recorded"="")
    selected <- if (!is.na(cached_chr) && cached_chr %in% unname(opts)) cached_chr else ""
    return(selectInput(id, lab, opts, selected=selected, selectize=FALSE))
  }
  if (v %in% continuous_predictors) {
    num <- suppressWarnings(as.numeric(cached_chr))
    val <- if (length(num)==1L && is.finite(num)) num else NA_real_
    lim <- numeric_limits[[v]]
    if (!is.null(lim)) return(tagList(numericInput(id,lab,value=val,min=lim$min,
                        max=lim$max,step=lim$step),div(class="input-hint",lim$note)))
    return(numericInput(id,lab,value=val))
  }
  opts <- choices_for(v)
  vals <- unname(opts)
  fallback <- if ("__UNKNOWN__" %in% vals) "__UNKNOWN__" else vals[[1L]]
  selected <- if (!is.na(cached_chr) && cached_chr %in% vals) cached_chr else fallback
  selectInput(id,lab,opts,selected=selected,selectize=FALSE)
}

validate_numeric_inputs <- function(input, allowed) {
  errs <- character()
  for (v in intersect(names(numeric_limits), allowed)) {
    raw <- input[[input_id(v)]]
    if (is.null(raw) || !length(raw) || is.na(raw[[1L]]) || identical(raw[[1L]], "")) next
    z <- suppressWarnings(as.numeric(raw[[1L]])); lim <- numeric_limits[[v]]
    if (!is.finite(z) || z < lim$min || z > lim$max)
      errs <- c(errs, paste0(label_for(v), " must be between ", lim$min, " and ", lim$max, "."))
  }
  errs
}
collect_values <- function(input, allowed) {
  vals <- setNames(vector("list", length(used)), used)
  for (v in used) {
    if (!(v %in% allowed)) {
      vals[[v]] <- if (v %in% categorical_predictors) "__UNKNOWN__" else NA_real_
      next
    }
    x <- input[[input_id(v)]]
    if (v %in% numeric_predictors) {
      z <- suppressWarnings(as.numeric(if (length(x)) x[[1L]] else NA_real_))
      vals[[v]] <- if (length(z)==1L && is.finite(z)) z else NA_real_
    } else {
      z <- if (length(x)) as.character(x[[1L]]) else ""
      vals[[v]] <- if (is.na(z) || !nzchar(z)) "__UNKNOWN__" else z
    }
  }
  vals
}

# Byte-for-byte feature order and preprocessing parity with Stage 21.
encode_values <- function(meta, vals, sparse = FALSE) {
  features <- as.character(meta$feature_names)
  X <- matrix(0, nrow=1L, ncol=length(features), dimnames=list(NULL,features))
  for (v in meta$predictors) {
    info <- meta$encoder$per_predictor[[v]]
    if (is.null(info)) stop("Missing encoder specification: ", v)
    if (identical(info$type, "numeric")) {
      if (!(v %in% features)) stop("Missing numeric encoded feature: ", v)
      X[1L, v] <- suppressWarnings(as.numeric(vals[[v]] %||% NA_real_))
    } else {
      levels <- as.character(info$levels)
      z <- trimws(as.character(vals[[v]] %||% "__UNKNOWN__"))
      if (is.na(z) || !nzchar(z) || z %in% c("NA","NaN","<NA>")) z <- "__UNKNOWN__"
      if (!(z %in% levels)) z <- "__OTHER__"
      if (!(z %in% levels)) stop("Encoder missing __OTHER__ fallback: ", v)
      all_names <- paste0(v, "__", make.names(levels, unique=TRUE))
      if (!all(all_names %in% features)) stop("Categorical feature mismatch for ", v)
      X[1L, all_names[[match(z,levels)]]] <- 1
    }
  }
  if (sparse) return(Matrix::Matrix(X, sparse=TRUE))
  storage.mode(X) <- "double"
  X
}

raw_class <- function(id, vals) {
  md <- metadata[[id]]
  X <- encode_values(md,vals,sparse=TRUE) # mirrors Stage 21 classification sparse input
  dm <- xgboost::xgb.DMatrix(X, missing=NA_real_)
  as.numeric(stats::predict(class_models[[id]],dm))
}
raw_duration <- function(id, vals) {
  md <- metadata[[id]]
  X <- encode_values(md,vals,sparse=FALSE) # mirrors Stage 21 duration dense input
  dm <- py_xgb$DMatrix(X, missing=NA_real_)
  z <- duration_models[[id]]$predict(dm)
  as.numeric(reticulate::py_to_r(z))
}
pred_multi <- function(id, vals) {
  p <- raw_class(id,vals)
  lev <- as.character(metadata[[id]]$outcome_levels)
  if (length(lev)!=3L || length(p)!=3L || any(!is.finite(p)) ||
      any(p < -1e-6) || any(p > 1+1e-6)) stop(id, " invalid 3-class output")
  p <- pmax(0,pmin(1,p))
  if (abs(sum(p)-1)>1e-4) stop(id, " multiclass probabilities fail normalization QC")
  data.table(class=lev,probability=p/sum(p))
}
pred_binary <- function(id, vals) {
  p <- raw_class(id,vals)
  if (length(p)!=1L || !is.finite(p) || p < -1e-6 || p > 1+1e-6)
    stop("Invalid operation probability")
  pmax(0,pmin(1,p))
}
pred_quantiles <- function(id, vals) {
  p <- raw_duration(id,vals)
  if (length(p)!=3L || any(!is.finite(p))) stop(id, " invalid Q10/Q50/Q90")
  bound <- as.numeric(metadata[[id]]$lower_bound %||% 0)
  z <- sort(pmax(expm1(p),bound))
  setNames(z,c("q10","q50","q90"))
}

fmt_prob <- function(p) {
  x <- 100*as.numeric(p)
  if (!is.finite(x)) return("—")
  if (x > 0 && x < .1) return("<0.1%")
  if (x < 5 || x > 95) sprintf("%.1f%%",x) else sprintf("%.0f%%",x)
}
fmt_days <- function(x) {
  x <- as.numeric(x)
  if (!is.finite(x)) return("—")
  if (x < 10) sprintf("%.1f d",x) else sprintf("%.0f d",x)
}
pretty_class <- function(x) {
  map <- c("Home/home health"="Home / home health",
           "Post-acute facility"="Post-acute facility", "Death/hospice"="Death / hospice",
           "<=7 days"="≤7 days", "8-19 days"="8–19 days", ">=20 days"="≥20 days",
           "No ICU"="No ICU", "ICU 1-7 days"="ICU 1–7 days", "ICU >=8 days"="ICU ≥8 days",
           "No ventilation"="No ventilation", "Ventilation 1-7 days"="Ventilation 1–7 days",
           "Ventilation >=8 days"="Ventilation ≥8 days")
  if (x %in% names(map)) unname(map[[x]]) else x
}
probability_card <- function(name, p, cls) {
  div(class=paste("result-card",cls),
      div(class="result-top",span(class="result-name",name),
          span(class="result-pct",fmt_prob(p))),
      div(class="bar", div(class="fill",style=sprintf("width:%.5f%%",100*p))))
}
prob_cards <- function(d, cls) div(class="result-grid",lapply(seq_len(nrow(d)),function(i)
  probability_card(pretty_class(d$class[[i]]),d$probability[[i]],cls)))
duration_card <- function(name,q,cls,note=NULL) {
  div(class=paste("duration-card",cls),
      div(class="duration-label",name),
      div(class="duration-main",fmt_days(q[["q50"]])),
      div(class="duration-range",paste0("Estimated 10th–90th percentile range: ",
                                        fmt_days(q[["q10"]])," – ",fmt_days(q[["q90"]]))),
      if (!is.null(note)) div(class="small-note",note))
}
section_ui <- function(title,note,body) div(class="result-section",
    div(class="result-heading",h3(title),span(note)),body)
traj_duration <- function(d,q,cls,title,note=NULL) tagList(prob_cards(d,cls),
   duration_card(title,q,cls,note))
operation_section <- function(p) section_ui("Recorded thoracolumbar operation",
  "Probability of recorded treatment receipt",div(class="result-grid",
   probability_card("Recorded thoracolumbar operative treatment",p,"result-neuro")))

# Same visual design system as the TBI-TRACT calculator.
css <- "
body{background:#f4f7fb;color:#243447}.app{max-width:1450px;margin:auto;padding:22px}.hero,.card{background:white;border:1px solid #e3ebf3;border-radius:22px;box-shadow:0 8px 26px rgba(31,52,73,.06)}
.hero{padding:22px 26px;margin-bottom:20px}.hero-grid{display:grid;grid-template-columns:84px 1fr;gap:18px;align-items:center}.logo{width:78px}.title{font-weight:850;font-size:clamp(1.8rem,3vw,2.8rem);margin:0}.subtitle{color:#5b6d7f;margin:4px 0 0}
.sticky{position:sticky;top:18px}.scroll{max-height:calc(100vh - 355px);overflow-y:auto;padding-right:5px}.mode{background:#edf4fb;border-radius:15px;padding:12px;margin-bottom:12px}.btn-primary{border-radius:13px;font-weight:750;min-height:44px}.section-title{font-weight:800}.result-section{padding:4px 0 19px;border-bottom:1px solid #edf1f5;margin-bottom:17px}.result-heading{display:flex;justify-content:space-between;align-items:baseline;gap:12px}.result-heading h3{font-size:1.05rem;font-weight:800}.result-heading span{font-size:.82rem;color:#6b7a8c;font-weight:650}.result-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:11px}.result-card,.duration-card{border:1px solid #dbe8f5;border-radius:17px;padding:13px;background:#f9fbfd}.result-top{display:flex;justify-content:space-between;gap:10px}.result-name{font-weight:700;color:#425466}.result-pct{font-size:1.7rem;font-weight:900}.bar{height:11px;background:#e7edf5;border-radius:99px;overflow:hidden;margin-top:10px}.fill{height:100%;border-radius:99px}.result-disposition .fill{background:#1f4e79}.result-hospital .fill{background:#2f7d32}.result-icu .fill{background:#b45f06}.result-vent .fill{background:#187b80}.result-neuro .fill{background:#6f42c1}.duration-card{margin-top:11px;background:white}.duration-main{font-size:2rem;font-weight:900}.duration-label,.duration-range{font-weight:700;color:#425466}.small-note,.footnote,.quick-note,.input-hint{font-size:.84rem;color:#6b7a8c}.input-hint{margin-top:-10px;margin-bottom:10px}.notice{background:#fff8e6;border:1px solid #f1d99a;border-radius:14px;padding:11px 13px;margin-bottom:13px}.details{margin-top:18px}.detail-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:18px 28px}.detail-grid h3{font-size:1rem;font-weight:800}@media(max-width:1199px){.sticky{position:static}.scroll{max-height:none}}@media(max-width:900px){.result-grid,.detail-grid{grid-template-columns:1fr}.hero-grid{grid-template-columns:1fr}.app{padding:12px}}
"

ui <- bslib::page_fluid(
  theme=bslib::bs_theme(version=5,bootswatch="flatly",primary="#1f4e79",
                        bg="#f4f7fb",fg="#243447"),
  tags$head(tags$style(HTML(css))),
  div(class="app",
    div(class="hero",div(class="hero-grid",
      div(if (file.exists("www/ohsu_logo.png")) img(src="ohsu_logo.png",class="logo") else strong("OHSU")),
      div(h1("SPINE-TRACT: Trauma Resource and Acute Care Trajectory Calculator",class="title"),
          p("Early inpatient trajectory estimates for adults with acute thoracolumbar spine trauma after initial trauma-center evaluation and diagnostic workup",class="subtitle")))),
    bslib::layout_columns(col_widths=c(4,8),
      div(class="sticky",bslib::card(bslib::card_body(
        h2("Calculator mode",class="section-title"),
        div(class="mode",radioButtons("mode",NULL,
           c("Complete calculator (recommended)"="complete","Quick endpoint preview"="quick"),"complete"),
           conditionalPanel("input.mode == 'quick'",
             selectInput("endpoint","Endpoint of interest",quick_choices,"operation",selectize=FALSE),
             div(class="quick-note",paste("Quick preview uses the same deployed endpoint model; hidden",
                    "predictors are passed as missing/unknown. Reduced-input predictions were not separately validated.")))),
        h2("Initial evaluation",class="section-title"),uiOutput("mode_note"),
        div(class="scroll",uiOutput("inputs")),
        actionButton("calc","Calculate trajectory",class="btn-primary w-100")
      ))),
      bslib::card(bslib::card_body(h2("Predicted inpatient trajectory",class="section-title"),
            uiOutput("results"),div(class="footnote",uiOutput("result_note"))))
    ),
    bslib::card(class="details",bslib::card_body(
      h2("Clinical context, performance, and intended use",class="section-title"),
      div(class="detail-grid",
        div(h3("Who this model represents"),p(paste(
          "Developed from 424,017 direct-presenting adults aged 18–89 years with acute qualifying",
          "thoracolumbar spine trauma in ACS TQP PUF 2020–2024. The 2024 temporal evaluation",
          "cohort contained 90,649 patients before endpoint-specific exclusions. Patients",
          "transferred in or without an observable index-hospital trajectory were excluded."))),
        div(h3("When to use it"),p(paste(
          "Use after the initial trauma-center evaluation and diagnostic workup when initial",
          "clinical findings, physiology, and relevant imaging-derived diagnostic injury",
          "phenotypes are available. Registry-derived diagnoses are not individually timestamped;",
          "some treatment or ICU/ventilation exposure can already be present at assessment."))),
        div(h3("What it estimates"),p(paste(
          "Three-class discharge disposition; categorized hospital LOS (≤7, 8–19, ≥20 days);",
          "three-class ICU and invasive mechanical-ventilation trajectories; unconditional",
          "hospital LOS and conditional ICU/ventilator-duration quantiles; and likelihood",
          "of recorded thoracolumbar operative treatment receipt."))),
        div(h3("2024 temporal performance"),p(paste(
          "AUROC: post-acute facility 0.776; death/hospice 0.925; any ICU 0.874;",
          "ICU ≥8 days 0.898; any ventilation 0.936; ventilation ≥8 days 0.932;",
          "hospital LOS ≥20 days 0.854; recorded thoracolumbar operation 0.761."))),
        div(h3("Duration performance"),p(paste(
          "2024 median-prediction MAE was 4.82 days for hospital LOS, 4.21 days",
          "for conditional ICU LOS, and 5.46 days for conditional ventilation duration.",
          "Observed coverage of estimated Q10–Q90 ranges was 79.5%, 82.0%, and 79.5%, respectively."))),
        div(h3("How to interpret it"),p(paste(
          "Predictions are estimates, not guarantees. Q10–Q90 forecasts are estimated quantile",
          "ranges, not guaranteed 80% prediction intervals. Race, ethnicity, and payer",
          "are used only for disposition and categorized hospital LOS and can encode inequities;",
          "model calibration alone does not establish equitable use."))),
        div(h3("Validation status"),p(paste(
          "Published performance reflects development-informed rolling-origin temporal evaluation",
          "during 2022–2024, not independent external validation. Deployment models were then fit",
          "on eligible 2020–2024 records, with internal CV used only to choose boosting rounds.",
          "External and prospective validation remain pending."))),
        div(h3("Important clinical limitation"),p(paste(
          "The registry lacks CT/MRI fracture morphology, posterior ligamentous complex integrity,",
          "TLICS/AO classifications and ASIA/AIS neurologic examinations. Operation reflects",
          "recorded treatment receipt—not operative indication or appropriateness. SPINE-TRACT",
          "must not determine spinal stability, operative selection, neurologic prognosis,",
          "or withholding of care. Intended for research, counseling and resource planning.")))
      )
    ))
  )
)

server <- function(input,output,session) {
  mode <- reactive(input$mode %||% "complete")
  endpoint <- reactive(input$endpoint %||% "operation")
  input_cache <- reactiveValues()
  lapply(used,function(v) local({
    var <- v; id <- input_id(var)
    observeEvent(input[[id]], {
      value <- input[[id]]
      if (!is.null(value) && length(value)) input_cache[[var]] <- value
    },ignoreNULL=TRUE,priority=100)
  }))
  allowed <- reactive(if (mode()=="complete") used else quick_sets[[endpoint()]])
  output$mode_note <- renderUI({
    if (mode()=="complete") {
      div(class="quick-note",paste0(
        "Complete mode displays all ",length(allowed())," inputs used across the final models.",
        " Unavailable information may be left unknown/blank. Input values are preserved when switching modes."))
    } else div(class="quick-note",paste0(
      "Quick preview shows ",length(allowed())," focused inputs; all other inputs are treated as",
      " missing/unknown. Entered values are preserved for this session."))
  })
  output$inputs <- renderUI({
    a <- allowed()
    panels <- lapply(names(groups),function(g) {
      vars <- intersect(groups[[g]],a)
      if (!length(vars)) return(NULL)
      ctrl <- lapply(vars,function(v) make_control(v,isolate(input_cache[[v]])))
      do.call(bslib::accordion_panel,c(list(title=g),ctrl))
    })
    panels <- panels[!vapply(panels,is.null,logical(1))]
    do.call(bslib::accordion,c(list(id="acc",open=FALSE),panels))
  })
  results_data <- eventReactive(input$calc,{
    errors <- validate_numeric_inputs(input,allowed())
    validate(need(!length(errors),paste(errors,collapse=" ")))
    vals <- collect_values(input,allowed()); ep <- endpoint()
    disp <- function() pred_multi("discharge3_final",vals)
    hlos <- function() list(trajectory=pred_multi("hlos3_final",vals),
                             duration=pred_quantiles("hospital_los",vals))
    icu <- function() list(trajectory=pred_multi("icu3_final",vals),
                            duration=pred_quantiles("icu_los_conditional",vals))
    vent <- function() list(trajectory=pred_multi("vent3_final",vals),
                             duration=pred_quantiles("ventilator_days_conditional",vals))
    operation <- function() pred_binary("operative_spine_intervention",vals)
    airway_positive <- isTRUE(vals$gcsq_intubated_recovered==1) ||
      (is.character(vals$respiratory_assistance_clean) &&
       grepl("assist",vals$respiratory_assistance_clean,ignore.case=TRUE))
    if (mode()=="quick") {
      if (ep=="disposition") return(list(mode="quick",ep=ep,disposition=disp()))
      if (ep=="hlos") return(c(list(mode="quick",ep=ep),hlos()))
      if (ep=="icu") return(c(list(mode="quick",ep=ep),icu()))
      if (ep=="ventilation") return(c(list(mode="quick",ep=ep,airway=airway_positive),vent()))
      return(list(mode="quick",ep=ep,operation=operation()))
    }
    list(mode="complete",ep="all",disposition=disp(),hlos=hlos(),icu=icu(),
         ventilation=vent(),operation=operation(),airway=airway_positive)
  },ignoreInit=FALSE)
  output$results <- renderUI({
    r <- results_data(); req(r)
    airway_note <- if (isTRUE(r$airway)) div(class="notice",paste(
      "Respiratory assistance or intubation is already present at initial assessment.",
      "Interpret any-ventilation probability as a full-hospitalization trajectory, not a",
      "purely future incidence estimate.")) else NULL
    if (r$mode=="quick") {
      if (r$ep=="disposition") return(section_ui("Discharge disposition",
        "Mutually exclusive probabilities",prob_cards(r$disposition,"result-disposition")))
      if (r$ep=="hlos") return(section_ui("Hospital length of stay",
        "Trajectory + continuous forecast",traj_duration(r$trajectory,r$duration,
        "result-hospital","Predicted median hospital LOS")))
      if (r$ep=="icu") return(section_ui("ICU trajectory",
        "Trajectory + conditional duration",traj_duration(r$trajectory,r$duration,
        "result-icu","If ICU care occurs: predicted median ICU LOS",
        "Duration is conditional on ICU use.")))
      if (r$ep=="ventilation") return(tagList(airway_note,section_ui(
        "Mechanical ventilation trajectory","Trajectory + conditional duration",
        traj_duration(r$trajectory,r$duration,"result-vent",
          "If ventilation occurs: predicted median ventilator duration",
          "Duration is conditional on invasive mechanical ventilation."))))
      return(operation_section(r$operation))
    }
    tagList(airway_note,
      section_ui("Discharge disposition","Mutually exclusive probabilities",
        prob_cards(r$disposition,"result-disposition")),
      section_ui("Hospital length of stay","Trajectory + continuous forecast",
        traj_duration(r$hlos$trajectory,r$hlos$duration,
          "result-hospital","Predicted median hospital LOS")),
      section_ui("ICU trajectory","Trajectory + conditional duration",
        traj_duration(r$icu$trajectory,r$icu$duration,"result-icu",
          "If ICU care occurs: predicted median ICU LOS","Duration is conditional on ICU use.")),
      section_ui("Mechanical ventilation trajectory","Trajectory + conditional duration",
        traj_duration(r$ventilation$trajectory,r$ventilation$duration,"result-vent",
          "If ventilation occurs: predicted median ventilator duration",
          "Duration is conditional on invasive mechanical ventilation.")),
      operation_section(r$operation)
    )
  })
  output$result_note <- renderUI({
    if (mode()=="quick") {
      span(paste(
        "Quick preview uses the same deployed endpoint model but omitted information is",
        "missing/unknown. Use Complete mode for the full-input model workflow."))
    } else span(paste(
      "Within each multiclass endpoint, probabilities sum to 100%.",
      "ICU and ventilation duration estimates are conditional on resource use.",
      "Ranges are estimated 10th–90th percentile intervals."))
  })
}
shinyApp(ui,server)
