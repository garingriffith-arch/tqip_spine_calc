
# ==============================================================================
# 01_spine_tract_gate1_audit_v2_1_1.R
# SPINE-TRACT Gate 1 v2.1: cohort / endpoint / data-integrity audit
#
# PURPOSE
#   - NO model fitting.
#   - NO random train/test split.
#   - NO threshold selection based on model performance.
#   - Build and audit the 2020-2024 direct-presentation thoracolumbar trauma cohort.
#   - Exclude transfer-ins from the primary direct-presentation population.
#   - Exclude transfer-outs from trajectory analyses.
#   - Test whether a 4-class discharge outcome is feasible.
#   - Audit HLOS 14/20/28-day thresholds.
#   - Audit ICU / ventilator 8-day thresholds and duration distributions.
#   - Audit spine phenotype definitions, including S24.2 / S34.2 / S34.3
#     separately from S34.4 (lumbosacral plexus injury).
#   - Audit the manually reviewed thoracolumbar ICD-10-PCS operative map.
#
# HARDWARE TARGET
#   Your PC: i9-14900KF / 64 GB RAM / RTX 4090.
#   This script targets ~80% of available logical CPU threads for data.table.
#   The Gate 1 workload is CSV parsing / grouping / QC, so CUDA does not provide
#   a meaningful advantage here. The RTX 4090 should be used in later XGBoost
#   model-training scripts (device = "cuda", tree_method = "hist").
#
# EXPECTED RAW DATA STRUCTURE EXAMPLE
#   C:/Users/garin/OneDrive/OHSU/Research/TQIP/Griffith_G-selected/
#     PUF AY 2024/PUF AY 2024/CSV/PUF_TRAUMA.csv
#     PUF AY 2024/PUF AY 2024/CSV/PUF_ICDDIAGNOSIS.csv
#     PUF AY 2024/PUF AY 2024/CSV/PUF_ICDPROCEDURE.csv
#
# PROJECT
#   C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization
#
# REQUIRED PACKAGE
#   data.table
#
# RUN
#   source("C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization/R/01_spine_tract_gate1_audit_v2_1_1.R")
# ==============================================================================

rm(list = setdiff(ls(), character(0)))

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package 'data.table' is required. Install it with install.packages('data.table').", call. = FALSE)
}
suppressPackageStartupMessages(library(data.table))

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
OUTPUT_DIR <- file.path(PROJECT_DIR, "outputs", "gate1_audit_v2_1")
DERIVED_DIR <- file.path(PROJECT_DIR, "data", "derived", "gate1_audit_v2_1")

dir.create(REFERENCE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DERIVED_DIR, recursive = TRUE, showWarnings = FALSE)

OPERATIVE_MAP_CANDIDATES <- c(
  file.path(REFERENCE_DIR, "spine_operative_icd10pcs_map.csv"),
  file.path(R_DIR, "spine_operative_icd10pcs_map.csv"),
  file.path(PROJECT_DIR, "spine_operative_icd10pcs_map.csv")
)
OPERATIVE_MAP_FILE <- OPERATIVE_MAP_CANDIDATES[file.exists(OPERATIVE_MAP_CANDIDATES)][1]
if (length(OPERATIVE_MAP_FILE) == 0L || is.na(OPERATIVE_MAP_FILE)) {
  OPERATIVE_MAP_FILE <- file.path(REFERENCE_DIR, "spine_operative_icd10pcs_map.csv")
}

# Use about 80% of logical CPU threads. On a 14900KF (32 logical threads)
# this should generally resolve to ~26 data.table threads.
TARGET_CPU_FRACTION <- 0.80
LOGICAL_CORES <- parallel::detectCores(logical = TRUE)
if (is.na(LOGICAL_CORES) || LOGICAL_CORES < 1L) LOGICAL_CORES <- 1L
N_THREADS <- max(1L, min(LOGICAL_CORES, as.integer(round(LOGICAL_CORES * TARGET_CPU_FRACTION))))

data.table::setDTthreads(N_THREADS)

# Avoid accidental BLAS/OpenMP oversubscription in any helper code.
Sys.setenv(
  OMP_NUM_THREADS = as.character(N_THREADS),
  OPENBLAS_NUM_THREADS = as.character(N_THREADS),
  MKL_NUM_THREADS = as.character(N_THREADS),
  NUMEXPR_NUM_THREADS = as.character(N_THREADS)
)

message("SPINE-TRACT Gate 1 audit v2.1")
message("Raw root:      ", RAW_ROOT)
message("Project dir:   ", PROJECT_DIR)
message("Output dir:    ", OUTPUT_DIR)
message("Logical cores: ", LOGICAL_CORES)
message("DT threads:    ", data.table::getDTthreads(), " (~", round(100 * N_THREADS / LOGICAL_CORES), "% of logical threads)")
message("NOTE: GPU is intentionally not used in Gate 1; this stage is primarily disk / RAM / CPU parsing and aggregation.")

if (!dir.exists(RAW_ROOT)) stop("RAW_ROOT does not exist: ", RAW_ROOT, call. = FALSE)
if (!dir.exists(PROJECT_DIR)) stop("PROJECT_DIR does not exist: ", PROJECT_DIR, call. = FALSE)

# ------------------------------------------------------------------------------
# 1. Helpers
# ------------------------------------------------------------------------------

norm_key <- function(x) {
  z <- as.character(x)
  z <- trimws(z)
  z <- gsub("\\.0$", "", z)
  z[z %in% c("", "NA", "NaN", "NULL", "null")] <- NA_character_
  z
}

canon_name <- function(x) {
  x <- toupper(as.character(x))
  x <- gsub("[^A-Z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

norm_icd <- function(x) {
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
    y %in% c("yes", "y", "true", "t", "present", "selected")
}

clean_numeric <- function(x, lower = -Inf, upper = Inf) {
  z <- safe_num(x)
  z[is.infinite(z)] <- NA_real_
  # Common registry sentinel values; BIU-specific handling happens separately.
  z[z %in% c(-99, -88, -77, -9, -8, -7, 9999, 99999)] <- NA_real_
  z[z < lower | z > upper] <- NA_real_
  z
}

is_missing_like <- function(x) {
  z <- toupper(trimws(as.character(x)))
  is.na(x) | z %in% c(
    "", "NA", "N/A", "NULL", "UNKNOWN", "UNK",
    "NOT KNOWN", "NOT KNOWN/NOT RECORDED", "NOT RECORDED",
    "NOT APPLICABLE"
  )
}

starts_any <- function(x, prefixes) {
  out <- rep(FALSE, length(x))
  x[is.na(x)] <- ""
  for (p in prefixes) out <- out | startsWith(x, p)
  out
}

first_existing <- function(header, exact = character(), regex = character(), exclude = character()) {
  if (length(header) == 0) return(NA_character_)
  can <- canon_name(header)

  for (e in exact) {
    idx <- which(can == canon_name(e))
    if (length(idx) > 0L) return(header[idx[1L]])
  }

  keep <- rep(FALSE, length(header))
  for (r in regex) keep <- keep | grepl(r, can, ignore.case = TRUE, perl = TRUE)
  if (length(exclude) > 0L) {
    for (r in exclude) keep <- keep & !grepl(r, can, ignore.case = TRUE, perl = TRUE)
  }
  idx <- which(keep)
  if (length(idx) == 0L) NA_character_ else header[idx[1L]]
}

find_biu_for <- function(header, base_col, extra_regex = character()) {
  if (is.na(base_col) || !nzchar(base_col)) return(NA_character_)
  can <- canon_name(header)
  base <- canon_name(base_col)

  exact <- c(
    paste0(base, "_BIU"),
    paste0(base, "BIU"),
    paste0("BIU_", base)
  )
  idx <- which(can %in% exact)
  if (length(idx) > 0L) return(header[idx[1L]])

  patterns <- c(
    paste0("^", gsub("_", ".*", base), ".*BIU$"),
    paste0("^BIU.*", gsub("_", ".*", base), "$"),
    extra_regex
  )
  first_existing(header, regex = patterns)
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
  if (!dir.exists(year_root)) {
    stop("Annual folder not found for ", year, ": ", year_root, call. = FALSE)
  }

  csv_dirs <- unique(dirname(list.files(
    year_root,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )))
  if (length(csv_dirs) == 1L) return(csv_dirs)

  stop(
    "Could not uniquely resolve the CSV directory for ", year,
    ". Candidate directories: ", paste(csv_dirs, collapse = " | "),
    call. = FALSE
  )
}

find_unique_annual_file <- function(year, basename_target) {
  d <- annual_csv_dir(year)

  direct <- file.path(d, basename_target)
  if (file.exists(direct)) return(normalizePath(direct, winslash = "/", mustWork = TRUE))

  hits <- list.files(
    d,
    pattern = paste0("^", gsub("\\.", "\\\\.", basename_target), "$"),
    full.names = TRUE,
    recursive = FALSE,
    ignore.case = TRUE
  )
  if (length(hits) != 1L) {
    stop(
      "Expected exactly one ", basename_target, " for ", year,
      " in ", d, "; found ", length(hits),
      if (length(hits) > 0L) paste0(": ", paste(hits, collapse = " | ")) else "",
      call. = FALSE
    )
  }
  normalizePath(hits, winslash = "/", mustWork = TRUE)
}

read_header <- function(path) {
  names(fread(path, nrows = 0, check.names = FALSE, showProgress = FALSE))
}

decode_transfer <- function(x) {
  z_chr <- toupper(trimws(as.character(x)))
  z_num <- suppressWarnings(as.integer(gsub("\\.0$", "", z_chr)))
  out <- rep(NA_character_, length(x))
  out[z_num == 1L | z_chr %in% c("YES", "Y", "TRUE")] <- "Yes"
  out[z_num == 2L | z_num == 0L | z_chr %in% c("NO", "N", "FALSE")] <- "No"
  out
}

decode_discharge <- function(x) {
  z <- suppressWarnings(as.integer(gsub("\\.0$", "", as.character(x))))
  out <- rep(NA_character_, length(z))
  out[z %in% c(3L, 6L)] <- "Home/home health"
  out[z == 11L] <- "Inpatient rehabilitation"
  out[z %in% c(2L, 7L, 12L, 13L, 14L)] <- "Other post-acute facility"
  out[z %in% c(5L, 8L)] <- "Death/hospice"
  out[z == 1L] <- "Transfer to short-term acute hospital"
  out[z == 4L] <- "AMA/discontinued care"
  out[z == 10L] <- "Court/law enforcement"
  out
}

decode_ed_discharge <- function(x) {
  z <- suppressWarnings(as.integer(gsub("\\.0$", "", as.character(x))))
  out <- rep(NA_character_, length(z))
  out[z %in% c(4L, 9L)] <- "Home/home health"
  out[z == 5L] <- "Death/hospice"
  out[z == 10L] <- "AMA/discontinued care"
  out[z == 11L] <- "Transfer to short-term acute hospital"
  out[z == 6L] <- "ED other/institutional/custody"
  # Defensive future-facing mapping; not expected for AY 2020-2024.
  out[z == 13L] <- "Death/hospice"
  out
}

repair_final_disposition <- function(hospital_raw, ed_raw) {
  hosp <- decode_discharge(hospital_raw)
  ed <- decode_ed_discharge(ed_raw)
  ed_code <- suppressWarnings(as.integer(gsub("\\.0$", "", as.character(ed_raw))))

  final <- hosp
  source <- fifelse(!is.na(hosp), "Hospital discharge disposition", NA_character_)

  use_ed <- is.na(final) & !is.na(ed)
  final[use_ed] <- ed[use_ed]
  source[use_ed] <- "ED discharge disposition"

  status <- rep("Retain", length(final))
  status[is.na(final)] <- "Unresolved after disposition repair"
  status[final == "Transfer to short-term acute hospital"] <- "Exclude: short-term acute-care transfer"
  status[final == "AMA/discontinued care"] <- "Exclude: AMA"
  status[final %in% c("Court/law enforcement", "ED other/institutional/custody")] <-
    "Exclude: custody/ED-other institutional"

  admitted_unresolved <- is.na(final) & ed_code %in% c(1L, 2L, 3L, 7L, 8L)
  status[admitted_unresolved] <- "Unresolved: admitted ED disposition but hospital disposition missing"

  data.table(
    final_disposition = final,
    final_disposition_source = source,
    trajectory_status = status
  )
}

duration_from_value_biu <- function(value, biu = NULL, upper = 575) {
  val <- clean_numeric(value, lower = 0, upper = upper)

  if (is.null(biu)) {
    return(data.table(
      duration = val,
      duration_status = fifelse(!is.na(val), "recorded_numeric", "unknown_or_unresolved")
    ))
  }

  b_chr <- toupper(trimws(as.character(biu)))
  b_num <- suppressWarnings(as.integer(gsub("\\.0$", "", b_chr)))

  structural_zero <- b_num == 1L | b_chr %in% c("NOT APPLICABLE", "NA")
  unknown <- b_num == 2L | b_chr %in% c(
    "NOT KNOWN/NOT RECORDED", "NOT KNOWN", "NOT RECORDED", "UNKNOWN", "NK/NR"
  )

  out <- val
  out[structural_zero %in% TRUE] <- 0
  out[unknown %in% TRUE] <- NA_real_

  status <- rep("unknown_or_unresolved", length(out))
  status[!is.na(val) & !(structural_zero %in% TRUE) & !(unknown %in% TRUE)] <- "recorded_numeric"
  status[structural_zero %in% TRUE] <- "structural_zero_not_applicable"
  status[unknown %in% TRUE] <- "unknown_not_recorded"

  data.table(duration = out, duration_status = status)
}

smd_cont <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) != 2L || length(x) == 0L) return(NA_real_)
  lev <- sort(unique(g))
  a <- x[g == lev[1]]; b <- x[g == lev[2]]
  if (length(a) < 2L || length(b) < 2L) return(NA_real_)
  den <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (is.na(den) || den == 0) return(0)
  (mean(b) - mean(a)) / den
}

smd_binary <- function(x, g) {
  ok <- !is.na(x) & !is.na(g)
  x <- x[ok]; g <- g[ok]
  if (length(unique(g)) != 2L || length(x) == 0L) return(NA_real_)
  lev <- sort(unique(g))
  p1 <- mean(x[g == lev[1]] == 1)
  p2 <- mean(x[g == lev[2]] == 1)
  den <- sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / 2)
  if (is.na(den) || den == 0) return(0)
  (p2 - p1) / den
}

q_safe <- function(x, probs) {
  if (all(is.na(x))) return(rep(NA_real_, length(probs)))
  as.numeric(stats::quantile(x, probs = probs, na.rm = TRUE, names = FALSE, type = 7))
}

top_values <- function(dt, year, variable, n = 20L) {
  if (!variable %in% names(dt)) return(NULL)
  tmp <- dt[, .N, by = .(value = as.character(get(variable)))][order(-N)]
  tmp[, `:=`(year = year, variable = variable)]
  setcolorder(tmp, c("year", "variable", "value", "N"))
  head(tmp, n)
}

# ------------------------------------------------------------------------------
# 2. Hardware audit
# ------------------------------------------------------------------------------

gpu_info <- tryCatch(
  system2(
    "nvidia-smi",
    c(
      "--query-gpu=name,memory.total,driver_version",
      "--format=csv,noheader"
    ),
    stdout = TRUE,
    stderr = TRUE
  ),
  error = function(e) paste("nvidia-smi unavailable:", conditionMessage(e))
)

hardware_audit <- data.table(
  item = c(
    "logical_cores",
    "target_cpu_fraction",
    "data_table_threads",
    "gpu_query",
    "R_version",
    "platform"
  ),
  value = c(
    as.character(LOGICAL_CORES),
    as.character(TARGET_CPU_FRACTION),
    as.character(data.table::getDTthreads()),
    paste(gpu_info, collapse = " | "),
    R.version.string,
    R.version$platform
  )
)
fwrite(hardware_audit, file.path(OUTPUT_DIR, "00_hardware_audit.csv"))

# ------------------------------------------------------------------------------
# 3. Resolve annual files and freeze one source file per required type per year
# ------------------------------------------------------------------------------

manifest <- rbindlist(lapply(YEARS, function(y) {
  files <- c(
    trauma = find_unique_annual_file(y, "PUF_TRAUMA.csv"),
    diagnosis = find_unique_annual_file(y, "PUF_ICDDIAGNOSIS.csv"),
    procedure = find_unique_annual_file(y, "PUF_ICDPROCEDURE.csv")
  )

  data.table(
    admission_year = y,
    file_type = names(files),
    path = unname(files),
    size_mb = round(file.info(unname(files))$size / 1024^2, 1),
    modified = as.character(file.info(unname(files))$mtime)
  )
}))
fwrite(manifest, file.path(OUTPUT_DIR, "01_source_file_manifest.csv"))

message("Resolved annual files successfully.")

# ------------------------------------------------------------------------------
# 4. Build annual diagnosis phenotype flags
# ------------------------------------------------------------------------------

diagnosis_field_map <- list()
diagnosis_flags_by_year <- list()
diagnosis_code_audit <- list()

for (y in YEARS) {
  message("\n[", y, "] Diagnosis phenotype audit")
  f <- manifest[admission_year == y & file_type == "diagnosis", path]
  hdr <- read_header(f)

  id_col <- first_existing(
    hdr,
    exact = c("INC_KEY", "INCKEY"),
    regex = c("^INC.*KEY$", "TRAUMA.*KEY", "PATIENT.*ID")
  )
  code_col <- first_existing(
    hdr,
    exact = c("ICD10CM", "ICD10DIAGNOSIS", "ICDDIAGNOSIS", "DIAGNOSISCODE"),
    regex = c("ICD.*DIAG", "DIAG.*CODE", "ICD10.*CM"),
    exclude = c("AIS", "LOOKUP", "DESCRIPTION", "VERSION")
  )

  if (is.na(id_col) || is.na(code_col)) {
    stop("Could not identify diagnosis key/code columns in ", f, call. = FALSE)
  }

  diagnosis_field_map[[as.character(y)]] <- data.table(
    admission_year = y,
    construct = c("inc_key", "icd10cm_diagnosis"),
    source_column = c(id_col, code_col)
  )

  d <- fread(
    f,
    select = c(id_col, code_col),
    integer64 = "character",
    showProgress = TRUE,
    nThread = N_THREADS
  )
  setnames(d, c(id_col, code_col), c("inc_key", "dx_code_raw"))
  d[, inc_key := norm_key(inc_key)]
  d[, dx_code := norm_icd(dx_code_raw)]
  d <- d[!is.na(inc_key) & !is.na(dx_code) & nchar(dx_code) >= 3L]

  # Exact spine phenotype concepts.
  d[, `:=`(
    dx_thoracic_vertebral_fracture = as.integer(starts_any(dx_code, c("S220", "S221"))),
    dx_lumbar_vertebral_fracture = as.integer(starts_any(dx_code, c("S320"))),

    dx_thoracic_spinal_cord_injury = as.integer(starts_any(dx_code, c("S240", "S241"))),
    dx_thoracic_nerve_root_injury = as.integer(starts_any(dx_code, c("S242"))),

    dx_lumbar_spinal_cord_injury = as.integer(starts_any(dx_code, c("S340", "S341"))),
    dx_lumbar_sacral_nerve_root_injury = as.integer(starts_any(dx_code, c("S342"))),
    dx_cauda_equina_injury = as.integer(starts_any(dx_code, c("S343"))),
    dx_lumbosacral_plexus_injury = as.integer(starts_any(dx_code, c("S344"))),

    dx_thoracic_ligament_dislocation_sprain = as.integer(starts_any(dx_code, c("S230", "S231", "S232", "S233", "S234"))),
    dx_lumbar_ligament_dislocation_sprain = as.integer(starts_any(dx_code, c("S330", "S331", "S332", "S333", "S334", "S335"))),

    dx_cervical_spine_injury = as.integer(starts_any(dx_code, c("S12", "S13", "S14"))),
    dx_tbi_any = as.integer(starts_any(dx_code, c("S06"))),
    dx_thoracic_injury_any = as.integer(starts_any(dx_code, sprintf("S%02d", 20:29))),
    dx_abdominal_pelvic_injury_any = as.integer(starts_any(dx_code, sprintf("S%02d", 30:39))),
    dx_upper_extremity_injury_any = as.integer(starts_any(dx_code, sprintf("S%02d", 40:69))),
    dx_lower_extremity_injury_any = as.integer(starts_any(dx_code, sprintf("S%02d", 70:99)))
  )]

  flag_cols <- grep("^dx_", names(d), value = TRUE)

  patient_flags <- d[, lapply(.SD, max, na.rm = TRUE), by = inc_key, .SDcols = flag_cols]

  patient_flags[, dx_tl_fracture_any := as.integer(
    dx_thoracic_vertebral_fracture == 1L | dx_lumbar_vertebral_fracture == 1L
  )]
  patient_flags[, dx_tl_spinal_cord_injury_any := as.integer(
    dx_thoracic_spinal_cord_injury == 1L | dx_lumbar_spinal_cord_injury == 1L
  )]
  patient_flags[, dx_tl_root_or_cauda_any := as.integer(
    dx_thoracic_nerve_root_injury == 1L |
      dx_lumbar_sacral_nerve_root_injury == 1L |
      dx_cauda_equina_injury == 1L
  )]
  patient_flags[, dx_tl_neural_primary_any := as.integer(
    dx_tl_spinal_cord_injury_any == 1L | dx_tl_root_or_cauda_any == 1L
  )]
  patient_flags[, dx_tl_soft_tissue_any := as.integer(
    dx_thoracic_ligament_dislocation_sprain == 1L |
      dx_lumbar_ligament_dislocation_sprain == 1L
  )]

  # Important correction: S34.4 lumbosacral plexus injury is audited separately
  # and does NOT qualify the narrow primary thoracolumbar spine cohort by itself.
  patient_flags[, spine_narrow_dx := as.integer(
    dx_tl_fracture_any == 1L | dx_tl_neural_primary_any == 1L
  )]
  patient_flags[, spine_broad_dx := as.integer(
    spine_narrow_dx == 1L | dx_tl_soft_tissue_any == 1L
  )]

  patient_flags[, fracture_region_pattern := fifelse(
    dx_thoracic_vertebral_fracture == 1L & dx_lumbar_vertebral_fracture == 1L, "Thoracic + lumbar",
    fifelse(dx_thoracic_vertebral_fracture == 1L, "Thoracic only",
            fifelse(dx_lumbar_vertebral_fracture == 1L, "Lumbar only", "No TL vertebral fracture"))
  )]

  patient_flags[, admission_year := y]
  diagnosis_flags_by_year[[as.character(y)]] <- patient_flags

  # Audit exact observed S24 / S34 diagnoses; useful for the root/cauda/plexus decision.
  diagnosis_code_audit[[as.character(y)]] <- d[
    starts_any(dx_code, c("S24", "S34")),
    .N,
    by = dx_code
  ][order(dx_code)][, admission_year := y]

  rm(d, patient_flags)
  invisible(gc())
}

diagnosis_field_map <- rbindlist(diagnosis_field_map, fill = TRUE)
fwrite(diagnosis_field_map, file.path(OUTPUT_DIR, "02_diagnosis_field_map.csv"))

diagnosis_code_audit <- rbindlist(diagnosis_code_audit, fill = TRUE)
setcolorder(diagnosis_code_audit, c("admission_year", "dx_code", "N"))
fwrite(diagnosis_code_audit, file.path(OUTPUT_DIR, "03_s24_s34_exact_code_audit.csv"))

# ------------------------------------------------------------------------------
# 5. Build annual trauma base and apply dictionary-aware endpoint handling
# ------------------------------------------------------------------------------

trauma_field_map <- list()
cohort_by_year <- list()
raw_value_audits <- list()

for (y in YEARS) {
  message("\n[", y, "] Trauma / endpoint audit")
  f <- manifest[admission_year == y & file_type == "trauma", path]
  hdr <- read_header(f)

  fields <- list(
    inc_key = first_existing(hdr, exact = c("INC_KEY", "INCKEY"), regex = c("^INC.*KEY$", "TRAUMA.*KEY", "PATIENT.*ID")),
    age = first_existing(hdr, exact = c("AGEYEARS", "AGE"), regex = c("^AGE.*YEAR", "^AGE$")),
    sex = first_existing(hdr, exact = c("SEX"), regex = c("^SEX$", "GENDER")),
    ethnicity = first_existing(hdr, exact = c("ETHNICITY"), regex = c("ETHNIC", "HISPANIC")),
    payer = first_existing(hdr, exact = c("PRIMARYMETHODPAYMENT", "PRIMARYPAYMENT"), regex = c("PAYMENT", "PAYER", "PAYOR", "INSURANCE")),

    # TQP PUF race is flattened into Y/N fields rather than one categorical RACE column.
    race_american_indian = first_existing(hdr, exact = c("AMERICANINDIAN")),
    race_asian = first_existing(hdr, exact = c("ASIAN")),
    race_black = first_existing(hdr, exact = c("BLACK")),
    race_pacific_islander = first_existing(hdr, exact = c("PACIFICISLANDER")),
    race_other = first_existing(hdr, exact = c("RACEOTHER")),
    race_white = first_existing(hdr, exact = c("WHITE")),
    race_na = first_existing(hdr, exact = c("RACE_NA")),
    race_uk = first_existing(hdr, exact = c("RACE_UK")),

    transfer = first_existing(hdr, exact = c("INTERFACILITYTRANSFER", "TRANSFER", "INTER_FACILITY_TRANSFER"), regex = c("INTER.*FACIL.*TRANSFER", "^TRANSFER$")),
    transfer_biu = first_existing(hdr, exact = c("INTERFACILITYTRANSFER_BIU"), regex = c("INTER.*FACIL.*TRANSFER.*BIU")),

    hospital_discharge = first_existing(
      hdr,
      exact = c("HOSPDISCHARGEDISPOSITION", "HOSPITALDISCHARGEDISPOSITION", "FINALDISCHARGEDISPOSITION"),
      regex = c("HOSP.*DISCH.*DISP", "FINAL.*DISCH.*DISP"),
      exclude = c("^ED")
    ),
    ed_discharge = first_existing(
      hdr,
      exact = c("EDDISCHARGEDISPOSITION"),
      regex = c("^ED.*DISCH.*DISP"),
      exclude = c("BIU")
    ),
    ed_discharge_biu = first_existing(
      hdr,
      exact = c("EDDISCHARGEDISPOSITION_BIU"),
      regex = c("^ED.*DISCH.*DISP.*BIU")
    ),

    hospital_los = first_existing(
      hdr,
      exact = c(
        "FINALDISCHARGEDAYS", "FINAL_DISCHARGE_DAYS", "FINALDISCHDAYS",
        "HOSPDISCHARGEDAYS", "HOSP_DISCHARGE_DAYS", "HOSPDISCHDAYS",
        "HOSPITALDISCHARGEDAYS", "HOSPITAL_DISCHARGE_DAYS",
        "HOSPDAYS", "HOSPITALDAYS", "HOSPLOS", "HOSPITALLOS", "TOTALHOSPITALDAYS", "HLOS", "TOTALLOS"
      ),
      regex = c("FINAL.*DISCH.*DAYS", "HOSP.*DISCH.*DAYS", "HOSP.*LOS", "TOTAL.*HOSP.*DAYS"),
      exclude = c("ARRIVAL", "ADMISSION", "^ED", "BIU")
    ),
    icu_los = first_existing(
      hdr,
      exact = c("TOTALICULOS", "TOTALICUDAYS", "ICUDAYS"),
      regex = c("TOTAL.*ICU.*LOS", "ICU.*LOS", "ICU.*DAYS"),
      exclude = c("BIU")
    ),
    vent_days = first_existing(
      hdr,
      exact = c("TOTALVENTDAYS", "TOTALVENTILATORDAYS"),
      regex = c("TOTAL.*VENT.*DAYS", "VENTILATOR.*DAYS", "VENT.*DAYS"),
      exclude = c("BIU")
    ),

    gcs_eye = first_existing(hdr, exact = c("GCSEYE"), regex = c("^GCS.*EYE$"), exclude = c("FIELD", "40", "BIU")),
    gcs_verbal = first_existing(hdr, exact = c("GCSVERBAL"), regex = c("^GCS.*VERBAL$"), exclude = c("FIELD", "40", "BIU")),
    gcs_motor = first_existing(hdr, exact = c("GCSMOTOR"), regex = c("^GCS.*MOTOR$"), exclude = c("FIELD", "40", "BIU", "HIGHEST")),
    gcs_total = first_existing(hdr, exact = c("TOTALGCS", "GCSTOTAL"), regex = c("^TOTAL.*GCS$", "^GCS.*TOTAL$"), exclude = c("FIELD", "HIGHEST", "BIU")),

    gcsq_sedated = first_existing(hdr, exact = c("GCSQ_SEDATEDPARALYZED")),
    gcsq_eye_obstruction = first_existing(hdr, exact = c("GCSQ_EYEOBSTRUCTION")),
    gcsq_intubated = first_existing(hdr, exact = c("GCSQ_INTUBATED")),
    gcsq_valid = first_existing(hdr, exact = c("GCSQ_VALID")),
    gcsq_na = first_existing(hdr, exact = c("GCSQ_NA")),
    gcsq_uk = first_existing(hdr, exact = c("GCSQ_UK")),

    sbp = first_existing(hdr, exact = c("SBP"), regex = c("^SBP$", "INITIAL.*SYSTOLIC"), exclude = c("FIELD", "LOWEST", "BIU")),
    hr = first_existing(hdr, exact = c("PULSERATE", "HR"), regex = c("^PULSE.*RATE$", "HEART.*RATE", "^HR$"), exclude = c("FIELD", "BIU")),
    rr = first_existing(hdr, exact = c("RESPIRATORYRATE", "RR"), regex = c("^RESPIRATORY.*RATE$", "^RR$"), exclude = c("FIELD", "ASSIST", "BIU")),
    spo2 = first_existing(hdr, exact = c("PULSEOXIMETRY", "SPO2"), regex = c("PULSE.*OX", "OXYGEN.*SAT", "SPO2"), exclude = c("FIELD", "BIU")),
    temperature = first_existing(hdr, exact = c("TEMPERATURE"), regex = c("^TEMPERATURE$"), exclude = c("BIU")),
    respiratory_assistance = first_existing(hdr, exact = c("RESPIRATORYASSISTANCE"), regex = c("^RESPIRATORY.*ASSISTANCE$"), exclude = c("BIU")),
    supplemental_oxygen = first_existing(hdr, exact = c("SUPPLEMENTALOXYGEN"), regex = c("^SUPPLEMENTAL.*OXYGEN$"), exclude = c("BIU")),
    primary_ecode = first_existing(hdr, exact = c("PRIMARYECODEICD10"), regex = c("^PRIMARY.*ECODE.*ICD10$")),
    pupil = first_existing(hdr, exact = c("TBIPUPILLARYRESPONSE"), regex = c("PUPILLARY.*RESPONSE"), exclude = c("BIU", "FIELD"))
  )

  fields$icu_los_biu <- find_biu_for(
    hdr, fields$icu_los,
    extra_regex = c("ICU.*LOS.*BIU", "BIU.*ICU.*LOS")
  )
  fields$vent_days_biu <- find_biu_for(
    hdr, fields$vent_days,
    extra_regex = c("VENT.*DAYS.*BIU", "BIU.*VENT.*DAYS")
  )

  required <- c("inc_key", "age", "transfer", "hospital_discharge", "ed_discharge", "hospital_los", "icu_los", "vent_days")
  missing_required <- required[vapply(fields[required], function(z) is.na(z) || !nzchar(z), logical(1))]
  if (length(missing_required) > 0L) {
    stop(
      "Missing required PUF_TRAUMA fields for ", y, ": ",
      paste(missing_required, collapse = ", "),
      ". Review the header and 04_trauma_field_map.csv after adjusting detection patterns.",
      call. = FALSE
    )
  }

  trauma_field_map[[as.character(y)]] <- data.table(
    admission_year = y,
    construct = names(fields),
    source_column = unlist(fields, use.names = FALSE)
  )

  select_cols <- unique(na.omit(unlist(fields, use.names = FALSE)))
  d <- fread(
    f,
    select = select_cols,
    integer64 = "character",
    showProgress = TRUE,
    nThread = N_THREADS
  )

  # Rename selected columns to stable internal names.
  for (nm in names(fields)) {
    src <- fields[[nm]]
    if (!is.na(src) && src %in% names(d)) setnames(d, src, nm)
  }

  d[, inc_key := norm_key(inc_key)]
  if (anyDuplicated(d$inc_key)) {
    dup <- d[duplicated(inc_key) | duplicated(inc_key, fromLast = TRUE)]
    fwrite(dup, file.path(OUTPUT_DIR, sprintf("FATAL_duplicate_trauma_keys_%d.csv", y)))
    stop("Duplicate INC_KEY rows in PUF_TRAUMA for ", y, ". See FATAL_duplicate_trauma_keys_", y, ".csv", call. = FALSE)
  }

  d[, admission_year := y]
  d[, age_years := clean_numeric(age, lower = 0, upper = 120)]
  d[, transfer_status := decode_transfer(transfer)]
  d[, hospital_discharge_detail := decode_discharge(hospital_discharge)]

  disp <- repair_final_disposition(d$hospital_discharge, d$ed_discharge)
  d[, `:=`(
    final_disposition = disp$final_disposition,
    final_disposition_source = disp$final_disposition_source,
    trajectory_status = disp$trajectory_status
  )]
  rm(disp)

  # Recovered GCS qualifier indicators using the same flattened PUF logic as
  # the finalized TBI-TRACT repair pipeline.
  gcs_valid <- if ("gcsq_valid" %in% names(d)) is_yes_flag(d$gcsq_valid) else rep(FALSE, nrow(d))
  gcs_intubated <- if ("gcsq_intubated" %in% names(d)) is_yes_flag(d$gcsq_intubated) else rep(FALSE, nrow(d))
  gcs_sedated <- if ("gcsq_sedated" %in% names(d)) is_yes_flag(d$gcsq_sedated) else rep(FALSE, nrow(d))
  gcs_eye_obs <- if ("gcsq_eye_obstruction" %in% names(d)) is_yes_flag(d$gcsq_eye_obstruction) else rep(FALSE, nrow(d))
  gcs_uk <- if ("gcsq_uk" %in% names(d)) is_yes_flag(d$gcsq_uk) else rep(FALSE, nrow(d))
  gcs_na <- if ("gcsq_na" %in% names(d)) is_yes_flag(d$gcsq_na) else rep(FALSE, nrow(d))
  gcs_any_known <- gcs_valid | gcs_intubated | gcs_sedated | gcs_eye_obs

  d[, `:=`(
    gcsq_intubated_recovered = as.integer(gcs_intubated),
    gcsq_sedated_paralyzed_recovered = as.integer(gcs_sedated),
    gcsq_eye_obstruction_recovered = as.integer(gcs_eye_obs),
    gcsq_valid_recovered = as.integer(gcs_valid),
    gcsq_unknown_recovered = as.integer(gcs_uk | gcs_na | !gcs_any_known),
    gcsq_conflict_recovered = as.integer((gcs_uk | gcs_na) & gcs_any_known)
  )]
  rm(gcs_valid, gcs_intubated, gcs_sedated, gcs_eye_obs, gcs_uk, gcs_na, gcs_any_known)

  # Admission-era fields recovered for later predictor adjudication.
  if ("temperature" %in% names(d)) {
    d[, temperature_c_clean := clean_numeric(temperature, lower = 10, upper = 45)]
  }
  if ("supplemental_oxygen" %in% names(d)) {
    so2n <- safe_num(d$supplemental_oxygen)
    so2t <- norm_text(d$supplemental_oxygen)
    d[, supplemental_oxygen_clean := fifelse(
      (!is.na(so2n) & so2n == 1) | so2t %in% c("no supplemental oxygen", "no"),
      "No supplemental oxygen",
      fifelse(
        (!is.na(so2n) & so2n == 2) | so2t %in% c("supplemental oxygen", "yes"),
        "Supplemental oxygen",
        "Unknown/not recorded"
      )
    )]
    rm(so2n, so2t)
  }

  # Dictionary-aware ICU / ventilator interpretation.
  icu_dec <- duration_from_value_biu(
    d$icu_los,
    if ("icu_los_biu" %in% names(d)) d$icu_los_biu else NULL
  )
  d[, `:=`(
    icu_los_days = icu_dec$duration,
    icu_los_status = icu_dec$duration_status
  )]
  rm(icu_dec)

  vent_dec <- duration_from_value_biu(
    d$vent_days,
    if ("vent_days_biu" %in% names(d)) d$vent_days_biu else NULL
  )
  d[, `:=`(
    ventilator_days_clean = vent_dec$duration,
    vent_days_status = vent_dec$duration_status
  )]
  rm(vent_dec)

  d[, hospital_los_days := clean_numeric(hospital_los, lower = 0, upper = 3650)]

  # Candidate admission-era numeric fields; no imputation in Gate 1.
  for (nm in c("gcs_eye", "gcs_verbal", "gcs_motor", "gcs_total", "sbp", "hr", "rr", "spo2")) {
    if (nm %in% names(d)) {
      bounds <- switch(
        nm,
        gcs_eye = c(1, 4),
        gcs_verbal = c(1, 5),
        gcs_motor = c(1, 6),
        gcs_total = c(3, 15),
        sbp = c(1, 350),
        hr = c(1, 300),
        rr = c(1, 100),
        spo2 = c(1, 100),
        c(-Inf, Inf)
      )
      d[, (paste0(nm, "_clean")) := clean_numeric(get(nm), bounds[1], bounds[2])]
    }
  }

  # Merge diagnosis phenotype.
  dx <- diagnosis_flags_by_year[[as.character(y)]]
  d <- merge(d, dx, by = c("inc_key", "admission_year"), all.x = TRUE, sort = FALSE)

  dx_cols <- grep("^dx_|^spine_", names(d), value = TRUE)
  dx_cols <- setdiff(dx_cols, c("spine_narrow_dx", "spine_broad_dx"))
  for (cc in unique(c(dx_cols, "spine_narrow_dx", "spine_broad_dx"))) {
    if (cc %in% names(d)) set(d, which(is.na(d[[cc]])), cc, 0L)
  }

  # Cohort flags.
  # AGEYEARS is the released PUF age field; values >89 are suppressed, so the
  # recorded-age adult cohort is effectively 18-89.
  d[, adult_recorded_age := !is.na(age_years) & age_years >= 18]
  d[, adult_narrow_spine := adult_recorded_age & spine_narrow_dx == 1L]
  d[, adult_broad_spine := adult_recorded_age & spine_broad_dx == 1L]

  # Primary presentation rule: known direct presentation only.
  d[, direct_presentation_known := transfer_status == "No"]
  d[, transfer_in := transfer_status == "Yes"]
  d[, transfer_status_unknown := is.na(transfer_status)]

  # Downstream trajectory status is based on REPAIRED final disposition:
  # hospital disposition where present, otherwise direct-ED disposition.
  d[, transfer_out := final_disposition == "Transfer to short-term acute hospital"]
  d[, future_ama := final_disposition == "AMA/discontinued care"]
  d[, future_custody_ed_other := final_disposition %in% c(
    "Court/law enforcement", "ED other/institutional/custody"
  )]
  d[, disposition_unresolved := grepl("^Unresolved", trajectory_status)]

  # TBI-TRACT-aligned complete index-center trajectory:
  #   - direct presentation
  #   - not later transferred to another acute-care hospital
  #   - not AMA
  #   - not custody/ED-other institutional
  #   - interpretable final disposition
  d[, primary_complete_trajectory := adult_narrow_spine &
       direct_presentation_known &
       trajectory_status == "Retain"]

  d[, broad_complete_trajectory := adult_broad_spine &
       direct_presentation_known &
       trajectory_status == "Retain"]

  d[, discharge_4class := fifelse(
    final_disposition %in% c(
      "Home/home health",
      "Inpatient rehabilitation",
      "Other post-acute facility",
      "Death/hospice"
    ),
    final_disposition,
    NA_character_
  )]

  # Binary / trajectory audit endpoints.
  d[, icu_admission := fifelse(
    icu_los_status == "unknown_not_recorded" | is.na(icu_los_days),
    NA_integer_,
    as.integer(icu_los_days > 0)
  )]
  d[, mechanical_ventilation := fifelse(
    vent_days_status == "unknown_not_recorded" | is.na(ventilator_days_clean),
    NA_integer_,
    as.integer(ventilator_days_clean > 0)
  )]

  d[, hlos_ge14 := fifelse(is.na(hospital_los_days), NA_integer_, as.integer(hospital_los_days >= 14))]
  d[, hlos_ge20 := fifelse(is.na(hospital_los_days), NA_integer_, as.integer(hospital_los_days >= 20))]
  d[, hlos_ge28 := fifelse(is.na(hospital_los_days), NA_integer_, as.integer(hospital_los_days >= 28))]
  d[, icu_los_ge8 := fifelse(is.na(icu_los_days), NA_integer_, as.integer(icu_los_days >= 8))]
  d[, vent_days_ge8 := fifelse(is.na(ventilator_days_clean), NA_integer_, as.integer(ventilator_days_clean >= 8))]

  # Impossible-value QC.
  d[, icu_gt_hlos := !is.na(icu_los_days) & !is.na(hospital_los_days) & icu_los_days > hospital_los_days]
  d[, vent_gt_hlos := !is.na(ventilator_days_clean) & !is.na(hospital_los_days) & ventilator_days_clean > hospital_los_days]

  cohort_by_year[[as.character(y)]] <- d

  # Raw-value audits for critical registry fields.
  for (v in intersect(
    c("age", "transfer", "hospital_discharge", "ed_discharge", "icu_los", "icu_los_biu", "vent_days", "vent_days_biu", "gcsq_intubated", "gcsq_sedated", "gcsq_eye_obstruction", "gcsq_valid", "gcsq_na", "gcsq_uk"),
    names(d)
  )) {
    raw_value_audits[[paste(y, v, sep = "_")]] <- top_values(d, y, v, 25L)
  }

  rm(d, dx)
  invisible(gc())
}

trauma_field_map <- rbindlist(trauma_field_map, fill = TRUE)
fwrite(trauma_field_map, file.path(OUTPUT_DIR, "04_trauma_field_map.csv"))

raw_value_audits <- rbindlist(raw_value_audits, fill = TRUE)
fwrite(raw_value_audits, file.path(OUTPUT_DIR, "05_critical_raw_value_audit.csv"))

cohort <- rbindlist(cohort_by_year, fill = TRUE, use.names = TRUE)

# Preserve a fast local derived file. compress=FALSE favors speed on your hardware.
saveRDS(
  cohort,
  file.path(DERIVED_DIR, "spine_gate1_v2_patient_level_2020_2024.rds"),
  compress = FALSE
)

# ------------------------------------------------------------------------------
# 6. Cohort flow
# ------------------------------------------------------------------------------

flow_rows <- list()

for (y in c(YEARS, NA_integer_)) {
  d <- if (is.na(y)) cohort else cohort[admission_year == y]
  label <- if (is.na(y)) "Pooled" else as.character(y)

  counts <- c(
    "All PUF_TRAUMA records" = nrow(d),
    "Recorded age 18-89" = d[adult_recorded_age == TRUE, .N],
    "Adult + narrow TL phenotype" = d[adult_narrow_spine == TRUE, .N],
    "Transfer-in = Yes" = d[adult_narrow_spine == TRUE & transfer_in == TRUE, .N],
    "Transfer status unknown" = d[adult_narrow_spine == TRUE & transfer_status_unknown == TRUE, .N],
    "Known direct presentation" = d[adult_narrow_spine == TRUE & direct_presentation_known == TRUE, .N],
    "Future transfer-out among direct" = d[
      adult_narrow_spine == TRUE & direct_presentation_known == TRUE & transfer_out == TRUE, .N
    ],
    "Future AMA among direct" = d[
      adult_narrow_spine == TRUE & direct_presentation_known == TRUE & future_ama == TRUE, .N
    ],
    "Future custody/ED-other among direct" = d[
      adult_narrow_spine == TRUE & direct_presentation_known == TRUE & future_custody_ed_other == TRUE, .N
    ],
    "Unresolved final disposition among direct" = d[
      adult_narrow_spine == TRUE & direct_presentation_known == TRUE & disposition_unresolved == TRUE, .N
    ],
    "PRIMARY complete-trajectory cohort" = d[primary_complete_trajectory == TRUE, .N],
    "4-class discharge eligible" = d[
      primary_complete_trajectory == TRUE & !is.na(discharge_4class), .N
    ],
    "Broad-phenotype complete-trajectory sensitivity cohort" = d[
      broad_complete_trajectory == TRUE, .N
    ]
  )

  flow_rows[[label]] <- data.table(
    admission_year = label,
    step = names(counts),
    n = as.integer(counts)
  )
}

cohort_flow <- rbindlist(flow_rows)
fwrite(cohort_flow, file.path(OUTPUT_DIR, "06_cohort_flow_by_year.csv"))

# ------------------------------------------------------------------------------
# 7. Spine phenotype audit
# ------------------------------------------------------------------------------

phenotype_vars <- c(
  "dx_thoracic_vertebral_fracture",
  "dx_lumbar_vertebral_fracture",
  "dx_thoracic_spinal_cord_injury",
  "dx_thoracic_nerve_root_injury",
  "dx_lumbar_spinal_cord_injury",
  "dx_lumbar_sacral_nerve_root_injury",
  "dx_cauda_equina_injury",
  "dx_lumbosacral_plexus_injury",
  "dx_tl_fracture_any",
  "dx_tl_spinal_cord_injury_any",
  "dx_tl_root_or_cauda_any",
  "dx_tl_neural_primary_any",
  "dx_tl_soft_tissue_any",
  "spine_narrow_dx",
  "spine_broad_dx",
  "dx_cervical_spine_injury",
  "dx_tbi_any"
)

phenotype_counts <- rbindlist(lapply(c(YEARS, NA_integer_), function(y) {
  d <- if (is.na(y)) cohort else cohort[admission_year == y]
  yy <- if (is.na(y)) "Pooled" else as.character(y)

  rbindlist(lapply(phenotype_vars, function(v) {
    data.table(
      admission_year = yy,
      phenotype = v,
      n_all_adult = d[adult_recorded_age == TRUE, sum(get(v) == 1L, na.rm = TRUE)],
      n_primary_complete_trajectory = d[primary_complete_trajectory == TRUE, sum(get(v) == 1L, na.rm = TRUE)],
      n_broad_complete_trajectory = d[broad_complete_trajectory == TRUE, sum(get(v) == 1L, na.rm = TRUE)]
    )
  }))
}))
fwrite(phenotype_counts, file.path(OUTPUT_DIR, "07_spine_phenotype_counts.csv"))

fracture_pattern <- cohort[
  primary_complete_trajectory == TRUE,
  .N,
  by = .(admission_year, fracture_region_pattern)
][, percent := 100 * N / sum(N), by = admission_year]
fwrite(fracture_pattern, file.path(OUTPUT_DIR, "08_fracture_region_pattern_by_year.csv"))

# ------------------------------------------------------------------------------
# 8. Four-class discharge feasibility audit
# ------------------------------------------------------------------------------

discharge_counts <- cohort[
  primary_complete_trajectory == TRUE,
  .N,
  by = .(admission_year, final_disposition)
][order(admission_year, -N)]
discharge_counts[, percent := 100 * N / sum(N), by = admission_year]
fwrite(discharge_counts, file.path(OUTPUT_DIR, "09_detailed_discharge_by_year.csv"))

discharge4 <- cohort[
  primary_complete_trajectory == TRUE & !is.na(discharge_4class),
  .N,
  by = .(admission_year, discharge_4class)
][order(admission_year, discharge_4class)]
discharge4[, percent := 100 * N / sum(N), by = admission_year]

pooled4 <- cohort[
  primary_complete_trajectory == TRUE & !is.na(discharge_4class),
  .N,
  by = discharge_4class
]
pooled4[, `:=`(
  admission_year = 9999L,
  percent = 100 * N / sum(N)
)]
discharge4_out <- rbind(
  discharge4,
  pooled4[, .(admission_year, discharge_4class, N, percent)],
  fill = TRUE
)
discharge4_out[admission_year == 9999L, admission_year := NA_integer_]
fwrite(discharge4_out, file.path(OUTPUT_DIR, "10_discharge_4class_feasibility.csv"))

discharge_gate <- discharge4[
  ,
  .(
    total_eligible = sum(N),
    smallest_class_n = min(N),
    smallest_class_percent = min(percent),
    n_classes_present = uniqueN(discharge_4class)
  ),
  by = admission_year
]
fwrite(discharge_gate, file.path(OUTPUT_DIR, "11_discharge_4class_gate_summary.csv"))

# ------------------------------------------------------------------------------
# 9. Hospital LOS threshold and duration audit
# ------------------------------------------------------------------------------

hlos_summary <- rbindlist(lapply(c(YEARS, NA_integer_), function(y) {
  d <- if (is.na(y)) cohort[primary_complete_trajectory == TRUE] else
    cohort[admission_year == y & primary_complete_trajectory == TRUE]
  yy <- if (is.na(y)) "Pooled" else as.character(y)

  qs <- q_safe(d$hospital_los_days, c(.10, .25, .50, .75, .90, .95, .99))
  data.table(
    admission_year = yy,
    n = nrow(d),
    n_hlos_known = sum(!is.na(d$hospital_los_days)),
    pct_hlos_missing = 100 * mean(is.na(d$hospital_los_days)),
    q10 = qs[1], q25 = qs[2], median = qs[3], q75 = qs[4],
    q90 = qs[5], q95 = qs[6], q99 = qs[7],
    pct_ge14 = 100 * mean(d$hlos_ge14 == 1L, na.rm = TRUE),
    pct_ge20 = 100 * mean(d$hlos_ge20 == 1L, na.rm = TRUE),
    pct_ge28 = 100 * mean(d$hlos_ge28 == 1L, na.rm = TRUE)
  )
}))
fwrite(hlos_summary, file.path(OUTPUT_DIR, "12_hlos_threshold_audit.csv"))

# ------------------------------------------------------------------------------
# 10. ICU / ventilation audit
# ------------------------------------------------------------------------------

resource_summary <- rbindlist(lapply(c(YEARS, NA_integer_), function(y) {
  d <- if (is.na(y)) cohort[primary_complete_trajectory == TRUE] else
    cohort[admission_year == y & primary_complete_trajectory == TRUE]
  yy <- if (is.na(y)) "Pooled" else as.character(y)

  icu_q <- q_safe(d[icu_admission == 1L, icu_los_days], c(.10, .50, .90))
  vent_q <- q_safe(d[mechanical_ventilation == 1L, ventilator_days_clean], c(.10, .50, .90))

  data.table(
    admission_year = yy,
    n = nrow(d),

    icu_known_n = sum(!is.na(d$icu_admission)),
    icu_admission_pct = 100 * mean(d$icu_admission == 1L, na.rm = TRUE),
    icu_ge8_pct_full = 100 * mean(d$icu_los_ge8 == 1L, na.rm = TRUE),
    icu_q10_if_used = icu_q[1],
    icu_q50_if_used = icu_q[2],
    icu_q90_if_used = icu_q[3],

    vent_known_n = sum(!is.na(d$mechanical_ventilation)),
    mechanical_ventilation_pct = 100 * mean(d$mechanical_ventilation == 1L, na.rm = TRUE),
    vent_ge8_pct_full = 100 * mean(d$vent_days_ge8 == 1L, na.rm = TRUE),
    vent_q10_if_used = vent_q[1],
    vent_q50_if_used = vent_q[2],
    vent_q90_if_used = vent_q[3],

    icu_gt_hlos_n = sum(d$icu_gt_hlos, na.rm = TRUE),
    vent_gt_hlos_n = sum(d$vent_gt_hlos, na.rm = TRUE)
  )
}))
fwrite(resource_summary, file.path(OUTPUT_DIR, "13_icu_ventilator_audit.csv"))

duration_status <- rbindlist(list(
  cohort[primary_complete_trajectory == TRUE, .N, by = .(admission_year, resource_status = icu_los_status)][
    , resource := "ICU LOS"
  ],
  cohort[primary_complete_trajectory == TRUE, .N, by = .(admission_year, resource_status = vent_days_status)][
    , resource := "Ventilator days"
  ]
), fill = TRUE)
setcolorder(duration_status, c("admission_year", "resource", "resource_status", "N"))
fwrite(duration_status, file.path(OUTPUT_DIR, "14_duration_BIU_status_by_year.csv"))

# ------------------------------------------------------------------------------
# 11. Candidate predictor missingness / raw-value audit
# ------------------------------------------------------------------------------

candidate_predictor_cols <- intersect(
  c(
    "age_years", "sex", "ethnicity", "payer",
    "gcs_eye_clean", "gcs_verbal_clean", "gcs_motor_clean", "gcs_total_clean",
    "gcsq_intubated_recovered", "gcsq_sedated_paralyzed_recovered",
    "gcsq_eye_obstruction_recovered", "gcsq_valid_recovered",
    "gcsq_unknown_recovered", "gcsq_conflict_recovered",
    "sbp_clean", "hr_clean", "rr_clean", "spo2_clean", "temperature_c_clean",
    "respiratory_assistance", "supplemental_oxygen_clean", "primary_ecode",
    "race_american_indian", "race_asian", "race_black", "race_pacific_islander",
    "race_other", "race_white", "race_na", "race_uk",
    "dx_thoracic_vertebral_fracture", "dx_lumbar_vertebral_fracture",
    "dx_thoracic_spinal_cord_injury", "dx_lumbar_spinal_cord_injury",
    "dx_thoracic_nerve_root_injury", "dx_lumbar_sacral_nerve_root_injury",
    "dx_cauda_equina_injury", "dx_cervical_spine_injury", "dx_tbi_any"
  ),
  names(cohort)
)

predictor_missingness <- rbindlist(lapply(YEARS, function(y) {
  d <- cohort[admission_year == y & primary_complete_trajectory == TRUE]
  rbindlist(lapply(candidate_predictor_cols, function(v) {
    x <- d[[v]]
    miss <- if (is.numeric(x) || is.integer(x)) is.na(x) else is_missing_like(x)
    data.table(
      admission_year = y,
      variable = v,
      n = length(x),
      missing_n = sum(miss),
      missing_percent = 100 * mean(miss)
    )
  }))
}))
fwrite(predictor_missingness, file.path(OUTPUT_DIR, "15_predictor_missingness_by_year.csv"))

predictor_value_audit <- rbindlist(lapply(YEARS, function(y) {
  d <- cohort[admission_year == y & primary_complete_trajectory == TRUE]
  rbindlist(lapply(
    intersect(c("sex", "ethnicity", "payer", "respiratory_assistance", "supplemental_oxygen_clean", "primary_ecode"), names(d)),
    function(v) top_values(d, y, v, 20L)
  ), fill = TRUE)
}), fill = TRUE)
fwrite(predictor_value_audit, file.path(OUTPUT_DIR, "16_categorical_predictor_top_values.csv"))

# ------------------------------------------------------------------------------
# 12. Selection audit: retained direct/no-transfer-out vs transfer-excluded
# ------------------------------------------------------------------------------

selection_base <- cohort[adult_narrow_spine == TRUE]
selection_base[, selection_group := fifelse(
  primary_complete_trajectory == TRUE,
  "Retained_complete_trajectory",
  "Excluded_transfer_or_trajectory_loss"
)]

selection_vars_cont <- intersect(
  c("age_years", "gcs_total_clean", "gcs_motor_clean", "sbp_clean", "hr_clean", "rr_clean", "spo2_clean"),
  names(selection_base)
)
selection_vars_bin <- intersect(
  c(
    "dx_thoracic_vertebral_fracture", "dx_lumbar_vertebral_fracture",
    "dx_thoracic_spinal_cord_injury", "dx_lumbar_spinal_cord_injury",
    "dx_thoracic_nerve_root_injury", "dx_lumbar_sacral_nerve_root_injury",
    "dx_cauda_equina_injury", "dx_cervical_spine_injury", "dx_tbi_any"
  ),
  names(selection_base)
)

selection_smd <- rbindlist(list(
  rbindlist(lapply(selection_vars_cont, function(v) {
    data.table(
      variable = v,
      type = "continuous",
      smd = smd_cont(selection_base[[v]], selection_base$selection_group)
    )
  })),
  rbindlist(lapply(selection_vars_bin, function(v) {
    data.table(
      variable = v,
      type = "binary",
      smd = smd_binary(selection_base[[v]], selection_base$selection_group)
    )
  }))
), fill = TRUE)
selection_smd[, abs_smd := abs(smd)]
selection_smd <- selection_smd[order(-abs_smd)]
fwrite(selection_smd, file.path(OUTPUT_DIR, "17_selection_SMD_retained_vs_transfer_excluded.csv"))

# ------------------------------------------------------------------------------
# 13. Procedure endpoint / timing audit
# ------------------------------------------------------------------------------

if (!file.exists(OPERATIVE_MAP_FILE)) {
  stop(
    "Required operative map not found. Expected one of:\n  ",
    paste(OPERATIVE_MAP_CANDIDATES, collapse = "\n  "),
    "\nPlace spine_operative_icd10pcs_map.csv in the project/reference folder and rerun. ",
    "Gate 1 v2 will not silently skip the operative endpoint audit.",
    call. = FALSE
  )
} else {
  op_map <- fread(OPERATIVE_MAP_FILE, showProgress = FALSE)
  setnames(op_map, names(op_map), canon_name(names(op_map)))

  if (!all(c("PROCEDURE_CODE", "SPINE_OPERATION_CATEGORY") %in% names(op_map))) {
    stop("Operative map must contain procedure_code and spine_operation_category columns.", call. = FALSE)
  }

  op_map[, PROCEDURE_CODE := norm_icd(PROCEDURE_CODE)]
  op_map <- unique(op_map[!is.na(PROCEDURE_CODE)])

  if (anyDuplicated(op_map$PROCEDURE_CODE)) {
    dup_map <- op_map[duplicated(PROCEDURE_CODE) | duplicated(PROCEDURE_CODE, fromLast = TRUE)]
    fwrite(dup_map, file.path(OUTPUT_DIR, "FATAL_duplicate_operative_map_codes.csv"))
    stop("Duplicate codes in operative map. See FATAL_duplicate_operative_map_codes.csv", call. = FALSE)
  }

  proc_field_map <- list()
  proc_matches <- list()
  proc_code_counts <- list()
  proc_timing_status <- list()

  for (y in YEARS) {
    message("\n[", y, "] Procedure endpoint audit")
    f <- manifest[admission_year == y & file_type == "procedure", path]
    hdr <- read_header(f)

    id_col <- first_existing(
      hdr,
      exact = c("INC_KEY", "INCKEY"),
      regex = c("^INC.*KEY$", "TRAUMA.*KEY", "PATIENT.*ID")
    )
    code_col <- first_existing(
      hdr,
      exact = c("ICD10PCS", "ICDPROCEDURE", "ICD10PROCEDURE"),
      regex = c("ICD.*PROCEDURE", "ICD10.*PCS", "PROCEDURE.*CODE", "PCS")
    )
    start_day_col <- first_existing(
      hdr,
      exact = c("HOSPITALPROCEDURESTARTDAYS", "HOSPITALPROCEDURESTARTDAY", "PROCEDURESTARTDAYS"),
      regex = c("HOSP.*PROC.*START.*DAY")
    )
    start_hrs_col <- first_existing(
      hdr,
      exact = c("HOSPITALPROCEDURESTARTHRS", "HOSPITALPROCEDURESTARTHOURS", "PROCEDURESTARTHRS"),
      regex = c("HOSP.*PROC.*START.*HRS", "HOSP.*PROC.*START.*HOUR")
    )
    start_date_col <- first_existing(
      hdr,
      exact = c("HOSPITALPROCEDURESTARTDATE", "PROCEDURESTARTDATE"),
      regex = c("HOSP.*PROC.*START.*DATE")
    )
    start_time_col <- first_existing(
      hdr,
      exact = c("HOSPITALPROCEDURESTARTTIME", "PROCEDURESTARTTIME"),
      regex = c("HOSP.*PROC.*START.*TIME")
    )

    if (is.na(id_col) || is.na(code_col)) {
      stop("Could not detect procedure key/code columns for ", y, call. = FALSE)
    }

    proc_field_map[[as.character(y)]] <- data.table(
      admission_year = y,
      construct = c("inc_key", "procedure_code", "procedure_start_days", "procedure_start_hrs", "procedure_start_date", "procedure_start_time"),
      source_column = c(id_col, code_col, start_day_col, start_hrs_col, start_date_col, start_time_col)
    )

    sel <- unique(na.omit(c(id_col, code_col, start_day_col, start_hrs_col, start_date_col, start_time_col)))
    p <- fread(
      f,
      select = sel,
      integer64 = "character",
      showProgress = TRUE,
      nThread = N_THREADS
    )
    setnames(p, id_col, "inc_key")
    setnames(p, code_col, "procedure_code")
    if (!is.na(start_day_col) && start_day_col %in% names(p)) setnames(p, start_day_col, "procedure_start_days")
    if (!is.na(start_hrs_col) && start_hrs_col %in% names(p)) setnames(p, start_hrs_col, "procedure_start_hrs")
    if (!is.na(start_date_col) && start_date_col %in% names(p)) setnames(p, start_date_col, "procedure_start_date")
    if (!is.na(start_time_col) && start_time_col %in% names(p)) setnames(p, start_time_col, "procedure_start_time")

    p[, inc_key := norm_key(inc_key)]
    p[, procedure_code := norm_icd(procedure_code)]
    p <- p[!is.na(inc_key) & !is.na(procedure_code)]

    # Restrict the endpoint audit to the primary transfer-clean cohort.
    cohort_keys <- cohort[admission_year == y & primary_complete_trajectory == TRUE, unique(inc_key)]
    p <- p[inc_key %chin% cohort_keys]

    proc_code_counts[[as.character(y)]] <- p[
      procedure_code %in% op_map$PROCEDURE_CODE,
      .N,
      by = procedure_code
    ][order(-N)][, admission_year := y]

    pm <- merge(
      p[procedure_code %in% op_map$PROCEDURE_CODE],
      op_map,
      by.x = "procedure_code",
      by.y = "PROCEDURE_CODE",
      all.x = TRUE,
      sort = FALSE
    )
    pm[, admission_year := y]
    proc_matches[[as.character(y)]] <- pm

    proc_timing_status[[as.character(y)]] <- data.table(
      admission_year = y,
      matched_procedure_rows = nrow(pm),
      matched_patients = uniqueN(pm$inc_key),
      procedure_start_days_field_present = "procedure_start_days" %in% names(pm),
      procedure_start_days_known_n = if ("procedure_start_days" %in% names(pm)) sum(!is_missing_like(pm$procedure_start_days)) else 0L,
      procedure_start_days_negative_n = if ("procedure_start_days" %in% names(pm)) sum(safe_num(pm$procedure_start_days) < 0, na.rm = TRUE) else 0L,
      procedure_start_hrs_field_present = "procedure_start_hrs" %in% names(pm),
      procedure_start_hrs_known_n = if ("procedure_start_hrs" %in% names(pm)) sum(!is_missing_like(pm$procedure_start_hrs)) else 0L,
      procedure_start_hrs_negative_n = if ("procedure_start_hrs" %in% names(pm)) sum(safe_num(pm$procedure_start_hrs) < 0, na.rm = TRUE) else 0L,
      procedure_start_date_field_present = "procedure_start_date" %in% names(pm),
      procedure_start_date_known_n = if ("procedure_start_date" %in% names(pm)) sum(!is_missing_like(pm$procedure_start_date)) else 0L,
      procedure_start_time_field_present = "procedure_start_time" %in% names(pm),
      procedure_start_time_known_n = if ("procedure_start_time" %in% names(pm)) sum(!is_missing_like(pm$procedure_start_time)) else 0L
    )

    rm(p, pm)
    invisible(gc())
  }

  proc_field_map <- rbindlist(proc_field_map, fill = TRUE)
  fwrite(proc_field_map, file.path(OUTPUT_DIR, "18_procedure_field_map.csv"))

  proc_code_counts <- rbindlist(proc_code_counts, fill = TRUE)
  fwrite(proc_code_counts, file.path(OUTPUT_DIR, "19_operative_code_counts_by_year.csv"))

  proc_matches <- rbindlist(proc_matches, fill = TRUE)
  saveRDS(
    proc_matches,
    file.path(DERIVED_DIR, "spine_gate1_v2_matched_operative_procedures_2020_2024.rds"),
    compress = FALSE
  )

  proc_category_counts <- proc_matches[
    ,
    .(
      procedure_rows = .N,
      patients = uniqueN(inc_key)
    ),
    by = .(admission_year, spine_operation_category = SPINE_OPERATION_CATEGORY)
  ][order(admission_year, -patients)]
  fwrite(proc_category_counts, file.path(OUTPUT_DIR, "20_operative_category_counts_by_year.csv"))

  proc_patient <- unique(proc_matches[, .(inc_key, admission_year)])
  proc_patient[, operative_spine_intervention := 1L]

  op_event_summary <- merge(
    cohort[primary_complete_trajectory == TRUE, .(inc_key, admission_year)],
    proc_patient,
    by = c("inc_key", "admission_year"),
    all.x = TRUE
  )
  op_event_summary[is.na(operative_spine_intervention), operative_spine_intervention := 0L]

  op_event_summary <- op_event_summary[
    ,
    .(
      n = .N,
      operative_events = sum(operative_spine_intervention),
      operative_percent = 100 * mean(operative_spine_intervention)
    ),
    by = admission_year
  ]
  fwrite(op_event_summary, file.path(OUTPUT_DIR, "21_operative_event_rate_by_year.csv"))

  proc_timing_status <- rbindlist(proc_timing_status, fill = TRUE)
  fwrite(proc_timing_status, file.path(OUTPUT_DIR, "22_procedure_timing_completeness.csv"))
}

# ------------------------------------------------------------------------------
# 14. Age-release / upper-bound empirical audit
# ------------------------------------------------------------------------------

age_audit <- cohort[
  ,
  .(
    n = .N,
    age_known_n = sum(!is.na(age_years)),
    min_age = suppressWarnings(min(age_years, na.rm = TRUE)),
    median_age = suppressWarnings(stats::median(age_years, na.rm = TRUE)),
    max_age = suppressWarnings(max(age_years, na.rm = TRUE)),
    n_age_90plus_released = sum(age_years >= 90, na.rm = TRUE),
    n_age_89 = sum(age_years == 89, na.rm = TRUE)
  ),
  by = admission_year
]
age_audit[!is.finite(min_age), min_age := NA_real_]
age_audit[!is.finite(median_age), median_age := NA_real_]
age_audit[!is.finite(max_age), max_age := NA_real_]
fwrite(age_audit, file.path(OUTPUT_DIR, "23_age_release_empirical_audit.csv"))

# ------------------------------------------------------------------------------
# 15. Gate 1 summary
# ------------------------------------------------------------------------------

summary_lines <- c(
  "SPINE-TRACT Gate 1 audit v2.1 completed.",
  paste0("Years: ", paste(YEARS, collapse = ", ")),
  paste0("Raw root: ", RAW_ROOT),
  paste0("CPU threads used by data.table: ", data.table::getDTthreads(), " / ", LOGICAL_CORES),
  "",
  "Key v2 corrections relative to the first audit run:",
  "- Repairs final disposition using hospital discharge disposition first, then direct-ED disposition when hospital disposition is absent.",
  "- Restores direct-ED home/home-health and death patients, mirroring the finalized TBI-TRACT trajectory logic.",
  "- Excludes future acute-care transfer, AMA, custody/ED-other institutional, and unresolved final trajectories from the primary cohort.",
  "- Uses TOTALGCS plus flattened GCSQ_* fields where available.",
  "- Audits flattened race fields rather than searching for a nonexistent single RACE column.",
  "- Produces a true broad-phenotype complete-trajectory sensitivity cohort.",
  "- Requires the operative ICD-10-PCS map and will not silently skip procedure auditing.",
  "",
  "Locked cohort logic in this audit:",
  "- Recorded PUF age >=18; AGEYEARS values >89 are suppressed in the PUF, so this is the released 18-89 cohort.",
  "- Narrow TL phenotype includes thoracic/lumbar vertebral fracture, thoracic/lumbosacral SCI, thoracic root injury, lumbar/sacral root injury, or cauda equina injury.",
  "- S34.4 lumbosacral plexus injury is audited separately and does not qualify the narrow cohort by itself.",
  "- Transfer-in = Yes and unknown transfer status are excluded.",
  "- Subsequent transfer-out, AMA, custody/ED-other institutional disposition, and unresolved final disposition are excluded from the complete trajectory cohort.",
  "- Four-class discharge is evaluated descriptively only.",
  "- HLOS 14/20/28 and ICU/ventilator 8-day thresholds are audited descriptively only.",
  "- BIU Not Applicable for ICU/ventilator duration is treated as structural zero; NK/NR remains missing.",
  "- No model fitting occurs in Gate 1.",
  "",
  "Review these first:",
  "06_cohort_flow_by_year.csv",
  "07_spine_phenotype_counts.csv",
  "10_discharge_4class_feasibility.csv",
  "11_discharge_4class_gate_summary.csv",
  "12_hlos_threshold_audit.csv",
  "13_icu_ventilator_audit.csv",
  "14_duration_BIU_status_by_year.csv",
  "15_predictor_missingness_by_year.csv",
  "17_selection_SMD_retained_vs_transfer_excluded.csv",
  "21_operative_event_rate_by_year.csv",
  "22_procedure_timing_completeness.csv",
  "23_age_release_empirical_audit.csv"
)

writeLines(summary_lines, file.path(OUTPUT_DIR, "README_GATE1_RESULTS.txt"))

message("\n============================================================")
message("SPINE-TRACT Gate 1 audit v2.1 complete.")
message("Outputs: ", OUTPUT_DIR)
message("Derived: ", DERIVED_DIR)
message("Start with: ", file.path(OUTPUT_DIR, "README_GATE1_RESULTS.txt"))
message("============================================================")
