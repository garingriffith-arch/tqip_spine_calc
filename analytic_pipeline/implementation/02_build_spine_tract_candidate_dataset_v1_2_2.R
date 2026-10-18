# ==============================================================================
# 02_build_spine_tract_candidate_dataset.R
# SPINE-TRACT Stage 2A — v1.2
# Fixes:
#   v1.1: split dependent data.table composite assignments.
#   v1.2: remove scalar endpoint labels from data.table by= lists and add
#         those labels after grouping.
#
# Build the FROZEN candidate analysis dataset and prespecified predictor families
# after Gate 1 v2.1.
#
# This script:
#   1) loads the Gate 1 v2.1 patient-level cohort;
#   2) freezes the primary complete-trajectory cohort (expected N = 424,017);
#   3) keeps the broad-phenotype complete-trajectory cohort for sensitivity;
#   4) rebuilds richer, admission-era diagnosis phenotypes from ICD-10-CM;
#   5) derives mechanism from the primary external-cause code;
#   6) derives selected baseline PMHx from PUF_PREEXISTINGCONDITIONS;
#   7) rebuilds the therapeutic thoracolumbar operation endpoint for the union
#      of primary + broad-sensitivity patients;
#   8) creates the locked endpoint architecture;
#   9) writes candidate predictor families and adjudication rules;
#  10) saves a frozen dataset + MD5 fingerprint for reproducibility.
#
# IMPORTANT:
#   - No model is fit here.
#   - No complete-case predictor restriction is applied.
#   - Numeric NA remains NA for XGBoost native missing routing.
#   - Categorical missingness is represented as "__UNKNOWN__".
#   - AIS/ISS and downstream treatment/complication fields are not predictors.
#
# RUN AFTER:
#   01_spine_tract_gate1_audit_v2_1.R
#
# HARDWARE:
#   Targets ~80% of logical CPU threads for data.table.
# ==============================================================================

rm(list = ls())
gc()

required_pkgs <- c("data.table")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0L) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(data.table)
})

# ------------------------------------------------------------------------------
# 0. Configuration
# ------------------------------------------------------------------------------

YEARS <- 2020:2024

RAW_ROOT <- Sys.getenv(
  "SPINE_TRACT_RAW_ROOT",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/Griffith_G-selected"
)

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

R_DIR <- file.path(PROJECT_DIR, "R")
REFERENCE_DIR <- file.path(PROJECT_DIR, "reference")
STAGE2_DIR <- file.path(PROJECT_DIR, "outputs", "stage2_predictor_adjudication")
DATA_DIR <- file.path(PROJECT_DIR, "data", "derived", "stage2_predictor_adjudication")

for (d in c(REFERENCE_DIR, STAGE2_DIR, DATA_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

EXPECTED_PRIMARY_N <- 424017L
EXPECTED_BROAD_N <- 440380L

TARGET_CPU_FRACTION <- 0.80
LOGICAL_CORES <- parallel::detectCores(logical = TRUE)
if (is.na(LOGICAL_CORES) || LOGICAL_CORES < 1L) LOGICAL_CORES <- 1L
N_THREADS <- max(1L, min(LOGICAL_CORES, as.integer(round(LOGICAL_CORES * TARGET_CPU_FRACTION))))
data.table::setDTthreads(N_THREADS)

Sys.setenv(
  OMP_NUM_THREADS = as.character(N_THREADS),
  OPENBLAS_NUM_THREADS = as.character(N_THREADS),
  MKL_NUM_THREADS = as.character(N_THREADS),
  NUMEXPR_NUM_THREADS = as.character(N_THREADS)
)

cat("\n============================================================\n")
cat("SPINE-TRACT STAGE 2A: FROZEN CANDIDATE DATASET\n")
cat("============================================================\n")
cat("Raw root:       ", RAW_ROOT, "\n", sep = "")
cat("Project:        ", PROJECT_DIR, "\n", sep = "")
cat("Logical cores:  ", LOGICAL_CORES, "\n", sep = "")
cat("DT threads:     ", getDTthreads(), "\n", sep = "")

if (!dir.exists(RAW_ROOT)) stop("RAW_ROOT does not exist: ", RAW_ROOT, call. = FALSE)
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR, call. = FALSE)

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

norm_key <- function(x) {
  z <- trimws(as.character(x))
  z <- gsub("\\.0$", "", z)
  z[z %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  z
}

norm_code <- function(x) {
  z <- toupper(as.character(x))
  z <- gsub("[^A-Z0-9]", "", z)
  z[z %in% c("", "NA", "NULL")] <- NA_character_
  z
}

safe_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
norm_text <- function(x) tolower(trimws(as.character(x)))

is_yes_flag <- function(x) {
  y <- norm_text(x)
  yn <- safe_num(x)
  (!is.na(yn) & yn == 1) |
    y %in% c("yes", "y", "true", "t", "present", "selected", "checked")
}

is_no_flag <- function(x) {
  y <- norm_text(x)
  yn <- safe_num(x)
  (!is.na(yn) & yn %in% c(0, 2)) |
    y %in% c("no", "n", "false", "f", "not selected")
}

clean_cat <- function(x, prefix = NULL) {
  z <- trimws(as.character(x))
  z[is.na(z) | z == "" | toupper(z) %in% c("NA", "N/A", "NULL", "UNKNOWN")] <- "__UNKNOWN__"
  if (!is.null(prefix)) {
    idx <- z != "__UNKNOWN__"
    z[idx] <- paste0(prefix, z[idx])
  }
  z
}

starts_any <- function(x, prefixes) {
  x2 <- x
  x2[is.na(x2)] <- ""
  out <- rep(FALSE, length(x2))
  for (p in prefixes) out <- out | startsWith(x2, p)
  out
}

annual_csv_dir <- function(year) {
  exact <- file.path(
    RAW_ROOT,
    sprintf("PUF AY %d", year),
    sprintf("PUF AY %d", year),
    "CSV"
  )
  if (dir.exists(exact)) return(exact)

  year_root <- file.path(RAW_ROOT, sprintf("PUF AY %d", year))
  if (!dir.exists(year_root)) stop("Annual folder not found: ", year_root, call. = FALSE)

  csvs <- list.files(
    year_root, pattern = "\\.csv$", recursive = TRUE,
    full.names = TRUE, ignore.case = TRUE
  )
  dirs <- unique(dirname(csvs))
  if (length(dirs) == 1L) return(dirs)

  stop(
    "Could not uniquely resolve CSV directory for ", year,
    ". Candidates: ", paste(dirs, collapse = " | "),
    call. = FALSE
  )
}

find_annual_file <- function(year, basename_target, required = TRUE) {
  d <- annual_csv_dir(year)
  direct <- file.path(d, basename_target)
  if (file.exists(direct)) return(normalizePath(direct, winslash = "/", mustWork = TRUE))

  hits <- list.files(
    d,
    pattern = paste0("^", gsub("\\.", "\\\\.", basename_target), "$"),
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) == 1L) return(normalizePath(hits, winslash = "/", mustWork = TRUE))

  if (required) {
    stop(
      "Expected exactly one ", basename_target, " for ", year,
      "; found ", length(hits), call. = FALSE
    )
  }
  NA_character_
}

read_header <- function(path) {
  names(fread(path, nrows = 0, showProgress = FALSE, check.names = FALSE))
}

canon <- function(x) {
  z <- toupper(as.character(x))
  z <- gsub("[^A-Z0-9]+", "_", z)
  z <- gsub("^_+|_+$", "", z)
  z
}

detect_col <- function(nms, exact = character(), regex = character(), exclude = character()) {
  cn <- canon(nms)
  for (e in exact) {
    idx <- which(cn == canon(e))
    if (length(idx) > 0L) return(nms[idx[1L]])
  }

  keep <- rep(FALSE, length(nms))
  for (r in regex) keep <- keep | grepl(r, cn, ignore.case = TRUE, perl = TRUE)
  if (length(exclude) > 0L) {
    for (r in exclude) keep <- keep & !grepl(r, cn, ignore.case = TRUE, perl = TRUE)
  }
  idx <- which(keep)
  if (length(idx) == 0L) NA_character_ else nms[idx[1L]]
}

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(hit[1L], winslash = "/", mustWork = TRUE)
}

# External-cause mechanism categories intentionally mirror the final TBI-TRACT
# conceptual groups while avoiding dependence on a year-specific numeric lookup.
mechanism_from_ecode <- function(x) {
  code <- norm_code(x)
  out <- rep("Other", length(code))
  out[is.na(code)] <- "__UNKNOWN__"

  transport <- startsWith(code, "V")
  fall <- starts_any(code, sprintf("W%02d", 0:19))

  firearm_explosive_prefix <- c(
    sprintf("W%02d", 32:40),
    sprintf("X%02d", 72:75),
    sprintf("X%02d", 93:96),
    sprintf("Y%02d", 22:25)
  )
  firearm_explosive <- starts_any(code, firearm_explosive_prefix)

  assault_prefix <- c(sprintf("X%02d", 85:99), sprintf("Y%02d", 0:9))
  assault <- starts_any(code, assault_prefix)

  # Priority matters because firearm assault should remain firearm/explosive.
  out[transport] <- "Transport"
  out[fall] <- "Fall"
  out[assault] <- "Assault"
  out[firearm_explosive] <- "Firearm/explosive"
  out
}

build_race_category <- function(d) {
  needed <- c(
    "race_american_indian", "race_asian", "race_black",
    "race_pacific_islander", "race_other", "race_white",
    "race_na", "race_uk"
  )
  if (!all(needed %in% names(d))) return(rep("__UNKNOWN__", nrow(d)))

  ai <- is_yes_flag(d$race_american_indian)
  asn <- is_yes_flag(d$race_asian)
  blk <- is_yes_flag(d$race_black)
  pac <- is_yes_flag(d$race_pacific_islander)
  oth <- is_yes_flag(d$race_other)
  wht <- is_yes_flag(d$race_white)
  nul <- is_yes_flag(d$race_na) | is_yes_flag(d$race_uk)

  mat <- cbind(ai, asn, blk, pac, oth, wht)
  nsel <- rowSums(mat, na.rm = TRUE)

  out <- rep("__UNKNOWN__", nrow(d))
  out[nul] <- "__UNKNOWN__"
  out[!nul & nsel > 1L] <- "Multiple"
  out[!nul & nsel == 1L & ai] <- "American Indian/Alaska Native"
  out[!nul & nsel == 1L & asn] <- "Asian"
  out[!nul & nsel == 1L & blk] <- "Black/African American"
  out[!nul & nsel == 1L & pac] <- "Native Hawaiian/Pacific Islander"
  out[!nul & nsel == 1L & oth] <- "Other"
  out[!nul & nsel == 1L & wht] <- "White"
  out
}

# ------------------------------------------------------------------------------
# 2. Locate and load Gate 1 v2.1 patient-level data
# ------------------------------------------------------------------------------

gate1_candidates <- c(
  file.path(
    PROJECT_DIR, "data", "derived", "gate1_audit_v2_1",
    "spine_gate1_v2_patient_level_2020_2024.rds"
  ),
  file.path(
    PROJECT_DIR, "data", "derived", "gate1_audit_v2_1",
    "spine_gate1_v2_1_patient_level_2020_2024.rds"
  )
)

gate1_file <- first_existing(gate1_candidates)

if (is.na(gate1_file)) {
  recursive_hits <- list.files(
    file.path(PROJECT_DIR, "data", "derived"),
    pattern = "spine_gate1.*patient_level.*2020_2024\\.rds$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(recursive_hits) > 0L) {
    info <- file.info(recursive_hits)
    gate1_file <- normalizePath(
      recursive_hits[order(info$mtime, decreasing = TRUE)][1L],
      winslash = "/",
      mustWork = TRUE
    )
  }
}

if (is.na(gate1_file) || !file.exists(gate1_file)) {
  stop(
    "Could not locate the Gate 1 v2.1 patient-level RDS. ",
    "Run 01_spine_tract_gate1_audit_v2_1.R first.",
    call. = FALSE
  )
}

cat("Gate 1 source:\n  ", gate1_file, "\n", sep = "")
dt_all <- as.data.table(readRDS(gate1_file))

required_gate1 <- c(
  "inc_key", "admission_year",
  "primary_complete_trajectory", "broad_complete_trajectory",
  "age_years", "hospital_los_days", "icu_los_days",
  "ventilator_days_clean", "discharge_4class"
)
missing_gate1 <- setdiff(required_gate1, names(dt_all))
if (length(missing_gate1) > 0L) {
  stop("Gate 1 RDS is missing: ", paste(missing_gate1, collapse = ", "), call. = FALSE)
}

dt_all[, inc_key := norm_key(inc_key)]
dt <- dt_all[broad_complete_trajectory == TRUE]
rm(dt_all)
gc()

if (anyDuplicated(dt[, paste(admission_year, inc_key, sep = "|")])) {
  stop("Duplicate patient-year keys in broad Gate 1 cohort.", call. = FALSE)
}

primary_n <- dt[primary_complete_trajectory == TRUE, .N]
broad_n <- nrow(dt)

if (primary_n != EXPECTED_PRIMARY_N) {
  stop(
    "Primary Gate 1 cohort changed: expected ", EXPECTED_PRIMARY_N,
    ", observed ", primary_n,
    ". Do not proceed until the cohort discrepancy is explained.",
    call. = FALSE
  )
}

if (broad_n != EXPECTED_BROAD_N) {
  stop(
    "Broad Gate 1 cohort changed: expected ", EXPECTED_BROAD_N,
    ", observed ", broad_n,
    ". Do not proceed until the cohort discrepancy is explained.",
    call. = FALSE
  )
}

cat("Primary N: ", format(primary_n, big.mark = ","), "\n", sep = "")
cat("Broad N:   ", format(broad_n, big.mark = ","), "\n", sep = "")

# Stable categorical representations.
dt[, sex_clean := clean_cat(sex, "SEX_")]
dt[, ethnicity_clean := clean_cat(ethnicity, "ETH_")]
dt[, payer_clean := clean_cat(payer, "PAYER_")]
dt[, race_clean := build_race_category(dt)]
dt[, mechanism_clean := mechanism_from_ecode(primary_ecode)]

if ("pupil" %in% names(dt)) {
  dt[, pupil_clean := clean_cat(pupil, "PUPIL_")]
} else {
  dt[, pupil_clean := "__UNKNOWN__"]
}

if ("respiratory_assistance" %in% names(dt)) {
  dt[, respiratory_assistance_clean := clean_cat(respiratory_assistance, "RESP_")]
} else {
  dt[, respiratory_assistance_clean := "__UNKNOWN__"]
}

if (!"supplemental_oxygen_clean" %in% names(dt)) {
  dt[, supplemental_oxygen_clean := "__UNKNOWN__"]
}

# Standardize numeric naming.
if ("hr_clean" %in% names(dt) && !"pulse_clean" %in% names(dt)) {
  setnames(dt, "hr_clean", "pulse_clean")
}

# ------------------------------------------------------------------------------
# 3. Rebuild richer admission-era ICD injury phenotypes
# ------------------------------------------------------------------------------

cat("\nDeriving spine-specific and associated-injury phenotypes...\n")

dx_flags_list <- list()
dx_file_audit <- list()

for (yr in YEARS) {
  f <- find_annual_file(yr, "PUF_ICDDIAGNOSIS.csv", required = TRUE)
  hdr <- read_header(f)

  id_col <- detect_col(
    hdr,
    exact = c("INC_KEY", "INCKEY"),
    regex = c("^INC.*KEY$", "TRAUMA.*KEY", "PATIENT.*ID")
  )
  code_col <- detect_col(
    hdr,
    exact = c("ICDDIAGNOSISCODE", "ICD10CM", "ICDDIAGNOSIS"),
    regex = c("ICD.*DIAG", "DIAG.*CODE", "ICD10.*CM"),
    exclude = c("AIS", "DESC", "LOOKUP", "VERSION")
  )

  if (is.na(id_col) || is.na(code_col)) {
    stop("Could not identify diagnosis columns for ", yr, call. = FALSE)
  }

  dx_file_audit[[as.character(yr)]] <- data.table(
    admission_year = yr,
    path = f,
    id_col = id_col,
    code_col = code_col
  )

  x <- fread(
    f,
    select = c(id_col, code_col),
    integer64 = "character",
    showProgress = TRUE,
    nThread = N_THREADS
  )
  setnames(x, c(id_col, code_col), c("inc_key", "code"))
  x[, inc_key := norm_key(inc_key)]
  x[, code := norm_code(code)]

  ids <- dt[admission_year == yr, unique(inc_key)]
  x <- x[inc_key %chin% ids & !is.na(code)]

  x[, `:=`(
    # Thoracolumbar primary phenotype and modifiers
    dx_thoracic_vertebral_fracture_tmp = starts_any(code, c("S220", "S221")),
    dx_lumbar_vertebral_fracture_tmp = starts_any(code, c("S320")),
    dx_thoracic_spinal_cord_injury_tmp = starts_any(code, c("S240", "S241")),
    dx_lumbar_spinal_cord_injury_tmp = starts_any(code, c("S340", "S341")),
    dx_thoracic_nerve_root_injury_tmp = starts_any(code, c("S242")),
    dx_lumbar_sacral_nerve_root_injury_tmp = starts_any(code, c("S342")),
    dx_cauda_equina_injury_tmp = starts_any(code, c("S343")),
    dx_lumbosacral_plexus_injury_tmp = starts_any(code, c("S344")),
    dx_thoracic_ligament_dislocation_sprain_tmp =
      starts_any(code, c("S230", "S231", "S232", "S233", "S234")),
    dx_lumbar_ligament_dislocation_sprain_tmp =
      starts_any(code, c("S330", "S331", "S332", "S333", "S334", "S335")),
    dx_cervical_spine_injury_tmp = starts_any(code, c("S12", "S13", "S14")),

    # Head / neck associated injury
    dx_tbi_any_tmp = startsWith(code, "S06"),
    dx_intracranial_hemorrhage_any_tmp =
      starts_any(code, c("S063", "S064", "S065", "S066")),
    dx_skull_fracture_any_tmp = startsWith(code, "S02"),
    dx_neck_vascular_injury_tmp = startsWith(code, "S15"),

    # Chest detail
    dx_rib_fracture_any_tmp = starts_any(code, c("S223", "S224", "S225")),
    dx_multiple_rib_fractures_tmp = startsWith(code, "S224"),
    dx_flail_chest_tmp = startsWith(code, "S225"),
    dx_sternum_fracture_tmp = startsWith(code, "S222"),
    dx_pulmonary_contusion_or_laceration_tmp = startsWith(code, "S273"),
    dx_traumatic_pneumothorax_tmp = startsWith(code, "S270"),
    dx_traumatic_hemothorax_tmp = startsWith(code, "S271"),
    dx_traumatic_hemopneumothorax_tmp = startsWith(code, "S272"),
    dx_diaphragm_injury_tmp = startsWith(code, "S278"),
    dx_thoracic_vascular_injury_tmp = startsWith(code, "S25"),
    dx_cardiac_injury_tmp = startsWith(code, "S26"),

    # Other associated injury
    dx_sacral_coccyx_fracture_tmp = starts_any(code, c("S321", "S322")),
    dx_pelvic_fracture_tmp = starts_any(code, c("S323", "S324", "S325", "S327", "S328")),
    dx_abdominal_visceral_vascular_injury_tmp =
      starts_any(code, c("S30", "S31", "S35", "S36", "S37", "S38", "S39")),
    dx_upper_extremity_injury_tmp = grepl("^S[456][0-9]", code),
    dx_lower_extremity_injury_tmp = grepl("^S[789][0-9]", code)
  )]

  tmp_cols <- grep("_tmp$", names(x), value = TRUE)

  sx <- x[, lapply(.SD, function(v) as.integer(any(v, na.rm = TRUE))),
          by = inc_key, .SDcols = tmp_cols]
  setnames(sx, tmp_cols, sub("_tmp$", "", tmp_cols))

  # data.table evaluates the right-hand side of a multi-column := before
  # assigning any of the new columns. Therefore, create the pleural composite
  # first, then use it in the broader chest composite in a second statement.
  sx[, `:=`(
    dx_tl_spinal_cord_injury_any =
      as.integer(dx_thoracic_spinal_cord_injury == 1L |
                   dx_lumbar_spinal_cord_injury == 1L),
    dx_tl_root_cauda_any =
      as.integer(dx_thoracic_nerve_root_injury == 1L |
                   dx_lumbar_sacral_nerve_root_injury == 1L |
                   dx_cauda_equina_injury == 1L),
    dx_tl_soft_tissue_any =
      as.integer(dx_thoracic_ligament_dislocation_sprain == 1L |
                   dx_lumbar_ligament_dislocation_sprain == 1L),
    dx_pleural_air_blood_injury =
      as.integer(dx_traumatic_pneumothorax == 1L |
                   dx_traumatic_hemothorax == 1L |
                   dx_traumatic_hemopneumothorax == 1L)
  )]

  sx[, dx_nonspinal_chest_injury :=
       as.integer(dx_rib_fracture_any == 1L |
                    dx_sternum_fracture == 1L |
                    dx_pulmonary_contusion_or_laceration == 1L |
                    dx_pleural_air_blood_injury == 1L |
                    dx_diaphragm_injury == 1L |
                    dx_thoracic_vascular_injury == 1L |
                    dx_cardiac_injury == 1L)]

  sx[, admission_year := yr]
  dx_flags_list[[as.character(yr)]] <- sx

  rm(x, sx)
  gc()
}

dx_file_audit <- rbindlist(dx_file_audit, fill = TRUE)
fwrite(dx_file_audit, file.path(STAGE2_DIR, "01_diagnosis_source_audit.csv"))

dx_flags <- rbindlist(dx_flags_list, fill = TRUE)
pre_n <- nrow(dt)
dt <- merge(dt, dx_flags, by = c("admission_year", "inc_key"), all.x = TRUE, sort = FALSE)
if (nrow(dt) != pre_n) stop("Diagnosis merge changed row count.", call. = FALSE)

derived_dx_cols <- setdiff(names(dx_flags), c("admission_year", "inc_key"))
for (v in derived_dx_cols) {
  set(dt, which(is.na(dt[[v]])), v, 0L)
}

# ------------------------------------------------------------------------------
# 4. Selected baseline PMHx
# ------------------------------------------------------------------------------

cat("\nDeriving selected baseline PMHx...\n")

pmhx_map <- data.table(
  pmhx_code = c("4", "31", "23", "11", "19", "8", "15", "26", "7", "9", "25", "24"),
  predictor = c(
    "pmhx_bleeding_disorder",
    "pmhx_anticoagulant_therapy",
    "pmhx_copd",
    "pmhx_diabetes",
    "pmhx_hypertension",
    "pmhx_current_smoker",
    "pmhx_functional_dependence",
    "pmhx_dementia",
    "pmhx_chf",
    "pmhx_chronic_renal_failure",
    "pmhx_cirrhosis",
    "pmhx_steroid_use"
  ),
  label = c(
    "Bleeding disorder",
    "Anticoagulant therapy",
    "COPD",
    "Diabetes mellitus",
    "Hypertension",
    "Current smoker",
    "Functionally dependent health status",
    "Dementia",
    "Congestive heart failure",
    "Chronic renal failure",
    "Cirrhosis",
    "Steroid use"
  )
)
fwrite(pmhx_map, file.path(STAGE2_DIR, "02_selected_pmhx_code_map.csv"))

pmhx_by_year <- list()
pmhx_source_audit <- list()

for (yr in YEARS) {
  f <- find_annual_file(yr, "PUF_PREEXISTINGCONDITIONS.csv", required = TRUE)
  hdr <- read_header(f)

  id_col <- detect_col(hdr, exact = c("INC_KEY", "INCKEY"), regex = c("^INC.*KEY$"))
  cond_col <- detect_col(
    hdr,
    exact = c("PREEXISTINGCONDITION"),
    regex = c("^PREEXISTING.*CONDITION$"),
    exclude = c("ANSWER", "BIU")
  )
  ans_col <- detect_col(
    hdr,
    exact = c("PREEXISTINGCONDITIONANSWER"),
    regex = c("^PREEXISTING.*CONDITION.*ANSWER$"),
    exclude = c("BIU")
  )
  biu_col <- detect_col(
    hdr,
    exact = c("PREEXISTINGCONDITIONANSWER_BIU"),
    regex = c("^PREEXISTING.*CONDITION.*ANSWER.*BIU$")
  )

  if (anyNA(c(id_col, cond_col, ans_col))) {
    stop("Could not detect PMHx columns for ", yr, call. = FALSE)
  }

  pmhx_source_audit[[as.character(yr)]] <- data.table(
    admission_year = yr,
    path = f,
    id_col = id_col,
    condition_col = cond_col,
    answer_col = ans_col,
    answer_biu_col = biu_col
  )

  select_cols <- unique(na.omit(c(id_col, cond_col, ans_col, biu_col)))
  p <- fread(
    f,
    select = select_cols,
    integer64 = "character",
    showProgress = TRUE,
    nThread = N_THREADS
  )

  setnames(p, id_col, "inc_key")
  setnames(p, cond_col, "pmhx_code")
  setnames(p, ans_col, "answer")
  if (!is.na(biu_col) && biu_col %in% names(p)) setnames(p, biu_col, "answer_biu")
  if (!"answer_biu" %in% names(p)) p[, answer_biu := NA_character_]

  p[, inc_key := norm_key(inc_key)]
  p[, pmhx_code := norm_key(pmhx_code)]

  ids <- dt[admission_year == yr, unique(inc_key)]
  p <- p[inc_key %chin% ids & pmhx_code %in% pmhx_map$pmhx_code]

  p[, answer_num := suppressWarnings(as.integer(gsub("\\.0$", "", as.character(answer))))]
  p[, answer_text := toupper(trimws(as.character(answer)))]
  p[, biu_num := suppressWarnings(as.integer(gsub("\\.0$", "", as.character(answer_biu))))]

  p[, pmhx_state := fcase(
    answer_num == 1L | answer_text %in% c("YES", "Y", "TRUE", "PRESENT"), 1L,
    answer_num == 2L | answer_text %in% c("NO", "N", "FALSE", "ABSENT"), 0L,
    biu_num == 2L, NA_integer_,
    default = NA_integer_
  )]

  # If duplicate rows exist for the same condition, a positive answer wins,
  # explicit No is next, and otherwise the condition remains unknown.
  agg <- p[, .(
    state = if (any(pmhx_state == 1L, na.rm = TRUE)) {
      1L
    } else if (any(pmhx_state == 0L, na.rm = TRUE)) {
      0L
    } else {
      NA_integer_
    }
  ), by = .(inc_key, pmhx_code)]

  agg <- merge(agg, pmhx_map[, .(pmhx_code, predictor)], by = "pmhx_code", all.x = TRUE)
  wide <- dcast(
    agg,
    inc_key ~ predictor,
    value.var = "state",
    fun.aggregate = function(z) {
      if (any(z == 1L, na.rm = TRUE)) 1L
      else if (any(z == 0L, na.rm = TRUE)) 0L
      else NA_integer_
    }
  )
  wide[, admission_year := yr]
  pmhx_by_year[[as.character(yr)]] <- wide

  rm(p, agg, wide)
  gc()
}

pmhx_source_audit <- rbindlist(pmhx_source_audit, fill = TRUE)
fwrite(pmhx_source_audit, file.path(STAGE2_DIR, "03_pmhx_source_audit.csv"))

pmhx_wide <- rbindlist(pmhx_by_year, fill = TRUE)
pre_n <- nrow(dt)
dt <- merge(dt, pmhx_wide, by = c("admission_year", "inc_key"), all.x = TRUE, sort = FALSE)
if (nrow(dt) != pre_n) stop("PMHx merge changed row count.", call. = FALSE)

pmhx_predictors <- pmhx_map$predictor
for (v in pmhx_predictors) {
  if (!v %in% names(dt)) dt[, (v) := NA_integer_]
}

# PMHx availability audit: do NOT replace unknown history with zero.
pmhx_missingness <- rbindlist(lapply(pmhx_predictors, function(v) {
  dt[primary_complete_trajectory == TRUE,
     .(
       N = .N,
       positive_n = sum(get(v) == 1L, na.rm = TRUE),
       negative_n = sum(get(v) == 0L, na.rm = TRUE),
       unknown_n = sum(is.na(get(v))),
       unknown_pct = 100 * mean(is.na(get(v)))
     ),
     by = admission_year][, predictor := v]
}))
setcolorder(pmhx_missingness, c(
  "predictor", "admission_year", "N", "positive_n", "negative_n",
  "unknown_n", "unknown_pct"
))
fwrite(pmhx_missingness, file.path(STAGE2_DIR, "04_pmhx_missingness_prevalence_by_year.csv"))

# ------------------------------------------------------------------------------
# 5. Rebuild operative endpoint for primary + broad sensitivity union
# ------------------------------------------------------------------------------

cat("\nRebuilding operative endpoint for primary + broad union...\n")

operative_map_file <- first_existing(c(
  file.path(REFERENCE_DIR, "spine_operative_icd10pcs_map.csv"),
  file.path(R_DIR, "spine_operative_icd10pcs_map.csv"),
  file.path(PROJECT_DIR, "spine_operative_icd10pcs_map.csv")
))

if (is.na(operative_map_file)) {
  stop(
    "spine_operative_icd10pcs_map.csv not found in reference/, R/, or project root.",
    call. = FALSE
  )
}

op_map <- fread(operative_map_file, showProgress = FALSE)
setnames(op_map, names(op_map), tolower(names(op_map)))
if (!all(c("procedure_code", "spine_operation_category") %in% names(op_map))) {
  stop("Operative map requires procedure_code and spine_operation_category.", call. = FALSE)
}
op_map[, procedure_code := norm_code(procedure_code)]
op_map <- unique(op_map[!is.na(procedure_code)])

op_rows <- list()
op_source_audit <- list()

for (yr in YEARS) {
  f <- find_annual_file(yr, "PUF_ICDPROCEDURE.csv", required = TRUE)
  hdr <- read_header(f)

  id_col <- detect_col(
    hdr,
    exact = c("INC_KEY", "INCKEY"),
    regex = c("^INC.*KEY$")
  )
  code_col <- detect_col(
    hdr,
    exact = c("ICDPROCEDURECODE", "ICD10PCS", "ICDPROCEDURE"),
    regex = c("ICD.*PROCEDURE", "PROCEDURE.*CODE", "PCS")
  )
  day_col <- detect_col(
    hdr,
    exact = c("HOSPITALPROCEDURESTARTDAYS"),
    regex = c("HOSP.*PROC.*START.*DAY")
  )

  if (is.na(id_col) || is.na(code_col)) {
    stop("Could not detect procedure columns for ", yr, call. = FALSE)
  }

  op_source_audit[[as.character(yr)]] <- data.table(
    admission_year = yr, path = f, id_col = id_col,
    code_col = code_col, start_day_col = day_col
  )

  sel <- unique(na.omit(c(id_col, code_col, day_col)))
  p <- fread(
    f,
    select = sel,
    integer64 = "character",
    showProgress = TRUE,
    nThread = N_THREADS
  )
  setnames(p, id_col, "inc_key")
  setnames(p, code_col, "procedure_code")
  if (!is.na(day_col) && day_col %in% names(p)) setnames(p, day_col, "procedure_start_days")
  if (!"procedure_start_days" %in% names(p)) p[, procedure_start_days := NA_real_]

  p[, inc_key := norm_key(inc_key)]
  p[, procedure_code := norm_code(procedure_code)]
  p[, procedure_start_days_num := safe_num(procedure_start_days)]

  ids <- dt[admission_year == yr, unique(inc_key)]
  p <- p[
    inc_key %chin% ids &
      procedure_code %in% op_map$procedure_code &
      (is.na(procedure_start_days_num) | procedure_start_days_num >= 0)
  ]

  p <- merge(p, op_map, by = "procedure_code", all.x = TRUE, sort = FALSE)
  p[, admission_year := yr]
  op_rows[[as.character(yr)]] <- p

  rm(p)
  gc()
}

op_source_audit <- rbindlist(op_source_audit, fill = TRUE)
fwrite(op_source_audit, file.path(STAGE2_DIR, "05_operative_source_audit.csv"))

op_long <- rbindlist(op_rows, fill = TRUE)
op_patient <- unique(op_long[, .(admission_year, inc_key)])
op_patient[, operative_spine_intervention := 1L]

pre_n <- nrow(dt)
dt <- merge(dt, op_patient, by = c("admission_year", "inc_key"), all.x = TRUE, sort = FALSE)
if (nrow(dt) != pre_n) stop("Operative endpoint merge changed row count.", call. = FALSE)
dt[is.na(operative_spine_intervention), operative_spine_intervention := 0L]

op_rate <- dt[, .(
  N = .N,
  events = sum(operative_spine_intervention == 1L),
  prevalence_pct = 100 * mean(operative_spine_intervention == 1L)
), by = .(admission_year, cohort = fifelse(
  primary_complete_trajectory == TRUE, "Primary",
  "Broad-only sensitivity"
))]
fwrite(op_rate, file.path(STAGE2_DIR, "06_operative_endpoint_rate_union.csv"))

# ------------------------------------------------------------------------------
# 6. Locked endpoint architecture
# ------------------------------------------------------------------------------

dt[, discharge4 := factor(
  discharge_4class,
  levels = c(
    "Home/home health",
    "Inpatient rehabilitation",
    "Other post-acute facility",
    "Death/hospice"
  )
)]

dt[, hlos3 := factor(
  fifelse(
    is.na(hospital_los_days), NA_character_,
    fifelse(
      hospital_los_days <= 7, "<=7 days",
      fifelse(hospital_los_days < 20, "8-19 days", ">=20 days")
    )
  ),
  levels = c("<=7 days", "8-19 days", ">=20 days")
)]

dt[, icu3 := factor(
  fifelse(
    is.na(icu_los_days), NA_character_,
    fifelse(
      icu_los_days <= 0, "No ICU",
      fifelse(icu_los_days < 8, "ICU 1-7 days", "ICU >=8 days")
    )
  ),
  levels = c("No ICU", "ICU 1-7 days", "ICU >=8 days")
)]

dt[, vent3 := factor(
  fifelse(
    is.na(ventilator_days_clean), NA_character_,
    fifelse(
      ventilator_days_clean <= 0, "No ventilation",
      fifelse(
        ventilator_days_clean < 8,
        "Ventilation 1-7 days",
        "Ventilation >=8 days"
      )
    )
  ),
  levels = c(
    "No ventilation",
    "Ventilation 1-7 days",
    "Ventilation >=8 days"
  )
)]

endpoint_manifest <- data.table(
  endpoint = c(
    "discharge4", "hlos3", "icu3", "vent3",
    "operative_spine_intervention",
    "hospital_los_days", "icu_los_days", "ventilator_days_clean"
  ),
  model_role = c(
    "Primary multiclass", "Primary multiclass", "Primary multiclass",
    "Primary multiclass", "Primary binary",
    "Quantile duration", "Conditional quantile duration", "Conditional quantile duration"
  ),
  definition = c(
    "Home/home health | inpatient rehabilitation | other post-acute facility | death/hospice",
    "<=7 d | 8-19 d | >=20 d",
    "No ICU | ICU 1-7 d | ICU >=8 d",
    "No ventilation | ventilation 1-7 d | ventilation >=8 d",
    "Any mapped therapeutic thoracic/lumbar/thoracolumbar operation during index hospitalization",
    "Observed index-hospital length of stay; Q10/Q50/Q90",
    "Observed ICU LOS among ICU users; Q10/Q50/Q90",
    "Observed ventilator days among ventilated patients; Q10/Q50/Q90"
  )
)
fwrite(endpoint_manifest, file.path(STAGE2_DIR, "07_locked_endpoint_manifest.csv"))

# ------------------------------------------------------------------------------
# 7. Prespecified predictor families
# ------------------------------------------------------------------------------

predictor_groups <- list(
  demographics = c(
    "age_years", "sex_clean"
  ),

  gcs_components = c(
    "gcs_eye_clean", "gcs_verbal_clean", "gcs_motor_clean"
  ),

  gcs_qualifiers = c(
    "gcsq_intubated_recovered",
    "gcsq_sedated_paralyzed_recovered",
    "gcsq_eye_obstruction_recovered",
    "gcsq_unknown_recovered"
  ),

  physiology = c(
    "sbp_clean", "pulse_clean", "rr_clean", "spo2_clean", "temperature_c_clean"
  ),

  supplemental_oxygen = c(
    "supplemental_oxygen_clean"
  ),

  mechanism = c(
    "mechanism_clean"
  ),

  spine_granular = c(
    "dx_thoracic_vertebral_fracture",
    "dx_lumbar_vertebral_fracture",
    "dx_thoracic_spinal_cord_injury",
    "dx_lumbar_spinal_cord_injury",
    "dx_thoracic_nerve_root_injury",
    "dx_lumbar_sacral_nerve_root_injury",
    "dx_cauda_equina_injury",
    "dx_thoracic_ligament_dislocation_sprain",
    "dx_lumbar_ligament_dislocation_sprain",
    "dx_cervical_spine_injury"
  ),

  spine_coarse = c(
    "dx_thoracic_vertebral_fracture",
    "dx_lumbar_vertebral_fracture",
    "dx_tl_spinal_cord_injury_any",
    "dx_tl_root_cauda_any",
    "dx_tl_soft_tissue_any",
    "dx_cervical_spine_injury"
  ),

  associated_injury = c(
    "dx_tbi_any",
    "dx_intracranial_hemorrhage_any",
    "dx_skull_fracture_any",
    "dx_neck_vascular_injury",
    "dx_nonspinal_chest_injury",
    "dx_abdominal_visceral_vascular_injury",
    "dx_sacral_coccyx_fracture",
    "dx_pelvic_fracture",
    "dx_upper_extremity_injury",
    "dx_lower_extremity_injury"
  ),

  chest_detail = c(
    "dx_rib_fracture_any",
    "dx_multiple_rib_fractures",
    "dx_flail_chest",
    "dx_sternum_fracture",
    "dx_pulmonary_contusion_or_laceration",
    "dx_pleural_air_blood_injury",
    "dx_diaphragm_injury",
    "dx_thoracic_vascular_injury",
    "dx_cardiac_injury"
  ),

  pmhx = pmhx_predictors,

  pupil = c(
    "pupil_clean"
  ),

  airway_state = c(
    "respiratory_assistance_clean"
  ),

  social_context = c(
    "race_clean", "ethnicity_clean", "payer_clean"
  )
)

# Ensure every listed field exists before freezing.
all_group_predictors <- unique(unlist(predictor_groups, use.names = FALSE))
missing_predictors <- setdiff(all_group_predictors, names(dt))
if (length(missing_predictors) > 0L) {
  stop(
    "Candidate predictor fields missing from frozen dataset: ",
    paste(missing_predictors, collapse = ", "),
    call. = FALSE
  )
}

clinical_backbone <- unique(c(
  predictor_groups$demographics,
  predictor_groups$gcs_components,
  predictor_groups$gcs_qualifiers,
  predictor_groups$physiology,
  predictor_groups$supplemental_oxygen,
  predictor_groups$mechanism,
  predictor_groups$spine_granular,
  predictor_groups$associated_injury,
  predictor_groups$pmhx
))

predictor_sets <- list(
  clinical_backbone = clinical_backbone,

  minus_gcs_qualifiers = setdiff(
    clinical_backbone,
    predictor_groups$gcs_qualifiers
  ),

  minus_pmhx = setdiff(
    clinical_backbone,
    predictor_groups$pmhx
  ),

  minus_associated_injury = setdiff(
    clinical_backbone,
    predictor_groups$associated_injury
  ),

  coarse_spine = unique(c(
    setdiff(clinical_backbone, predictor_groups$spine_granular),
    predictor_groups$spine_coarse
  )),

  minus_spine_soft_tissue = setdiff(
    clinical_backbone,
    c(
      "dx_thoracic_ligament_dislocation_sprain",
      "dx_lumbar_ligament_dislocation_sprain"
    )
  ),

  plus_chest_detail = unique(c(
    clinical_backbone,
    predictor_groups$chest_detail
  )),

  plus_pupil = unique(c(
    clinical_backbone,
    predictor_groups$pupil
  )),

  plus_airway_state = unique(c(
    clinical_backbone,
    predictor_groups$airway_state
  )),

  plus_social_context = unique(c(
    clinical_backbone,
    predictor_groups$social_context
  ))
)

saveRDS(
  list(
    predictor_groups = predictor_groups,
    predictor_sets = predictor_sets,
    clinical_backbone = clinical_backbone
  ),
  file.path(DATA_DIR, "spine_tract_predictor_sets.rds"),
  compress = FALSE
)

predictor_set_table <- rbindlist(lapply(names(predictor_sets), function(set_name) {
  data.table(
    predictor_set = set_name,
    predictor = predictor_sets[[set_name]]
  )
}))
fwrite(predictor_set_table, file.path(STAGE2_DIR, "08_predictor_sets.csv"))

group_table <- rbindlist(lapply(names(predictor_groups), function(group_name) {
  data.table(
    predictor_group = group_name,
    predictor = predictor_groups[[group_name]]
  )
}))
fwrite(group_table, file.path(STAGE2_DIR, "09_predictor_groups.csv"))

adjudication_plan <- data.table(
  family = c(
    "Age + sex",
    "GCS components",
    "GCS qualifiers",
    "Admission physiology",
    "Supplemental oxygen",
    "Mechanism",
    "Granular spine phenotype",
    "Associated injury phenotype",
    "Selected baseline PMHx",
    "Detailed chest/pulmonary injury",
    "Pupillary response",
    "Respiratory assistance/current airway state",
    "Race/ethnicity/payer",
    "GCS total/severity",
    "AIS/ISS/registry severity",
    "Aggregate injury/comorbidity counts",
    "Transfer status",
    "Lumbosacral plexus injury"
  ),
  initial_role = c(
    "Prespecified retain",
    "Prespecified retain",
    "Empirically stress-test removal",
    "Prespecified retain",
    "Prespecified retain",
    "Prespecified retain",
    "Empirically compare granular vs coarse",
    "Empirically stress-test removal",
    "Empirically stress-test removal",
    "Optional augmentation",
    "Optional augmentation",
    "High-scrutiny optional augmentation; never for ventilation",
    "Test only for discharge + categorical HLOS",
    "Exclude redundant",
    "Exclude",
    "Exclude redundant / coding-intensity prone",
    "Exclude; transfer-in/unknown are outside primary population",
    "Audit-only; does not define narrow TL cohort"
  ),
  rationale = c(
    "Baseline demographics known at prediction time.",
    "Admission neurologic examination without redundant total score.",
    "Contextualizes intubated/sedated/obstructed examination; assess whether burden is justified in spine trauma.",
    "Early physiology plausibly available at intended prediction time.",
    "Early respiratory severity without using current invasive ventilation as an outcome proxy.",
    "Injury context known from primary external-cause code.",
    "Spine-specific anatomy/neurologic phenotype is central to the SPINE-TRACT use case.",
    "Polytrauma materially affects ICU, ventilation, LOS, surgery, and disposition.",
    "Admission-known reserve/frailty modifiers; preserve unknown PMHx rather than forcing No.",
    "May improve ICU/ventilation prediction beyond a broad chest flag but increases input burden and redundancy.",
    "TBI-oriented variable with questionable value for a spine-specific bedside tool; retain only if clearly additive.",
    "Can reveal treatment/support already underway and risks state leakage; test only on distal outcomes.",
    "Social/health-system context, not biological mechanism; evaluate only where placement/LOS processes plausibly mediate outcomes.",
    "Deterministic/redundant with GCS components.",
    "Retrospective registry severity summaries are not acceptable near-admission bedside inputs.",
    "Adds coding intensity and duplicates individual phenotypes.",
    "Prediction population is restricted to direct presentations with complete trajectory.",
    "S34.4 is anatomically distinct from cord/root/cauda injury and is not a qualifying primary phenotype."
  )
)
fwrite(adjudication_plan, file.path(STAGE2_DIR, "10_predictor_adjudication_plan.csv"))

# ------------------------------------------------------------------------------
# 8. QC: missingness, prevalence, endpoint counts
# ------------------------------------------------------------------------------

predictor_type <- function(x) {
  if (is.numeric(x) || is.integer(x)) "numeric" else "categorical"
}

predictor_qc <- rbindlist(lapply(unique(c(
  clinical_backbone,
  predictor_groups$chest_detail,
  predictor_groups$pupil,
  predictor_groups$airway_state,
  predictor_groups$social_context
)), function(v) {
  x <- dt[primary_complete_trajectory == TRUE][[v]]
  type <- predictor_type(x)

  if (type == "numeric") {
    data.table(
      predictor = v,
      type = type,
      n = length(x),
      missing_n = sum(is.na(x)),
      missing_pct = 100 * mean(is.na(x)),
      unique_nonmissing = uniqueN(x[!is.na(x)]),
      positive_pct = if (all(na.omit(unique(x)) %in% c(0, 1))) {
        100 * mean(x == 1, na.rm = TRUE)
      } else {
        NA_real_
      }
    )
  } else {
    xx <- as.character(x)
    miss <- is.na(xx) | xx == "" | xx == "__UNKNOWN__"
    data.table(
      predictor = v,
      type = type,
      n = length(x),
      missing_n = sum(miss),
      missing_pct = 100 * mean(miss),
      unique_nonmissing = uniqueN(xx[!miss]),
      positive_pct = NA_real_
    )
  }
}))
fwrite(predictor_qc, file.path(STAGE2_DIR, "11_candidate_predictor_QC.csv"))

endpoint_counts <- rbindlist(list(
  dt[primary_complete_trajectory == TRUE,
     .N,
     by = .(admission_year, level = as.character(discharge4))
  ][, endpoint := "discharge4"],

  dt[primary_complete_trajectory == TRUE,
     .N,
     by = .(admission_year, level = as.character(hlos3))
  ][, endpoint := "hlos3"],

  dt[primary_complete_trajectory == TRUE,
     .N,
     by = .(admission_year, level = as.character(icu3))
  ][, endpoint := "icu3"],

  dt[primary_complete_trajectory == TRUE,
     .N,
     by = .(admission_year, level = as.character(vent3))
  ][, endpoint := "vent3"],

  dt[primary_complete_trajectory == TRUE,
     .N,
     by = .(
       admission_year,
       level = as.character(operative_spine_intervention)
     )
  ][, endpoint := "operative_spine_intervention"]
), fill = TRUE)

setcolorder(endpoint_counts, c("admission_year", "endpoint", "level", "N"))
fwrite(endpoint_counts, file.path(STAGE2_DIR, "12_endpoint_counts_by_year.csv"))

# ------------------------------------------------------------------------------
# 9. Freeze datasets
# ------------------------------------------------------------------------------

# Keep the broad union with a flag so later sensitivity models use the same
# predictor derivation. Primary models always filter primary_complete_trajectory.
keep_cols <- unique(c(
  "admission_year", "inc_key",
  "primary_complete_trajectory", "broad_complete_trajectory",
  "discharge4", "hlos3", "icu3", "vent3",
  "operative_spine_intervention",
  "hospital_los_days", "icu_los_days", "ventilator_days_clean",
  "final_disposition", "final_disposition_source",
  "age_years", "gcs_total_clean",
  all_group_predictors
))
keep_cols <- intersect(keep_cols, names(dt))

frozen_union <- dt[, ..keep_cols]
setorder(frozen_union, admission_year, inc_key)

primary <- frozen_union[primary_complete_trajectory == TRUE]

if (nrow(primary) != EXPECTED_PRIMARY_N) {
  stop("Primary row count changed immediately before freeze.", call. = FALSE)
}
if (nrow(frozen_union) != EXPECTED_BROAD_N) {
  stop("Broad-union row count changed immediately before freeze.", call. = FALSE)
}

union_file <- file.path(DATA_DIR, "spine_tract_frozen_candidate_union_2020_2024.rds")
primary_file <- file.path(DATA_DIR, "spine_tract_frozen_primary_2020_2024.rds")

saveRDS(frozen_union, union_file, compress = FALSE)
saveRDS(primary, primary_file, compress = FALSE)

fingerprint <- data.table(
  artifact = c(
    "Frozen primary dataset",
    "Frozen primary+broad union",
    "Predictor-set specification"
  ),
  path = c(
    primary_file,
    union_file,
    file.path(DATA_DIR, "spine_tract_predictor_sets.rds")
  )
)
fingerprint[, size_mb := round(file.info(path)$size / 1024^2, 2)]
fingerprint[, md5 := unname(tools::md5sum(path))]
fingerprint[, created_at := as.character(Sys.time())]
fwrite(fingerprint, file.path(STAGE2_DIR, "13_frozen_dataset_fingerprint.csv"))

source_manifest <- data.table(
  source = c(
    "Gate 1 patient-level RDS",
    "Operative code map"
  ),
  path = c(
    gate1_file,
    operative_map_file
  )
)
fwrite(source_manifest, file.path(STAGE2_DIR, "14_stage2_source_manifest.csv"))

readme <- c(
  "SPINE-TRACT Stage 2A complete.",
  "",
  paste0("Primary frozen N: ", format(nrow(primary), big.mark = ",")),
  paste0("Broad-union frozen N: ", format(nrow(frozen_union), big.mark = ",")),
  "",
  "Primary candidate architecture:",
  "- 4-class discharge",
  "- HLOS <=7 / 8-19 / >=20 days + later Q10/Q50/Q90",
  "- ICU no / 1-7 / >=8 days + later conditional Q10/Q50/Q90",
  "- Ventilation no / 1-7 / >=8 days + later conditional Q10/Q50/Q90",
  "- Recorded therapeutic thoracolumbar operation",
  "",
  "Predictor policy:",
  "- Clinical backbone is prespecified; family-level ablations are performed in Stage 2B.",
  "- GCS total is retained for description/QC only, not the primary predictor set.",
  "- AIS/ISS and aggregate counts are not primary predictors.",
  "- Respiratory assistance is a high-scrutiny optional augmentation and is never tested for the ventilation endpoint.",
  "- Race/ethnicity/payer augmentation is tested only for discharge and categorical HLOS.",
  "- Numeric missingness is retained as NA.",
  "- Categorical missingness is explicit __UNKNOWN__.",
  "",
  "Next:",
  "Run 03_run_spine_tract_predictor_adjudication_gpu.R"
)
writeLines(readme, file.path(STAGE2_DIR, "README_STAGE2A.txt"))

cat("\n============================================================\n")
cat("STAGE 2A COMPLETE\n")
cat("============================================================\n")
cat("Frozen primary:\n  ", primary_file, "\n", sep = "")
cat("Frozen union:\n  ", union_file, "\n", sep = "")
cat("Predictor spec:\n  ", file.path(DATA_DIR, "spine_tract_predictor_sets.rds"), "\n", sep = "")
cat("Next: source R/03_run_spine_tract_predictor_adjudication_gpu.R\n")
