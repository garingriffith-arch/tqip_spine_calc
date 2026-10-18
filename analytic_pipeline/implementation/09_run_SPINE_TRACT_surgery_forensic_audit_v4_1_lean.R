# ==============================================================================
# 09_run_SPINE_TRACT_surgery_forensic_audit_v4_1_lean.R
# SPINE-TRACT Stage 3E v4.1 — LEAN operative-treatment forensic audit
#
# PURPOSE
#   Explain the decline in recorded thoracolumbar operative-treatment AUROC
#   WITHOUT rebuilding the historical pipeline from raw TQIP files.
#
# THIS VERSION USES ONLY EXISTING DERIVED ARTIFACTS:
#   1) current Gate 1 patient-level RDS (contains ALL narrow TL patients,
#      including transfer-ins);
#   2) Gate 1 matched operative-procedure RDS;
#   3) current frozen Stage 2 candidate dataset + predictor specification;
#   4) current Stage 3A selected surgery hyperparameters;
#   5) packaged historical operative ICD-10-PCS code map for direct comparison.
#
# NO RAW CSV REBUILD.
#
# PRIMARY QUESTIONS
#   A. Did the old vs current operation code map change?
#   B. How much signal is carried by transfer-in status / transfer population?
#   C. How much discrimination can be recovered in CURRENT direct presenters by
#      adding back defensible-but-redundant/contextual admission-era features?
#
# IMPORTANT
#   This is forensic only. It does not alter the locked analysis.
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

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_surgery_forensic_audit_lean"
)
PRED_DIR <- file.path(OUT_DIR, "predictions")
dir.create(PRED_DIR, recursive = TRUE, showWarnings = FALSE)

SEED <- 20261001L
MAX_ROUNDS <- 4000L
EARLY_STOP <- 100L

logical_cores <- parallel::detectCores(logical = TRUE)
if (is.na(logical_cores) || logical_cores < 1L) logical_cores <- 1L
CPU_THREADS <- max(
  1L,
  min(
    logical_cores,
    as.integer(round(logical_cores * 0.80))
  )
)
data.table::setDTthreads(CPU_THREADS)

first_existing <- function(x) {
  hit <- unique(x[file.exists(x)])
  if (length(hit) == 0L) return(NA_character_)
  normalizePath(hit[1L], winslash = "/", mustWork = TRUE)
}

latest_recursive <- function(root, pattern) {
  if (!dir.exists(root)) return(NA_character_)
  hits <- list.files(
    root,
    pattern = pattern,
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )
  if (length(hits) == 0L) return(NA_character_)
  info <- file.info(hits)
  normalizePath(
    hits[order(info$mtime, decreasing = TRUE)][1L],
    winslash = "/",
    mustWork = TRUE
  )
}

# ------------------------------------------------------------------------------
# 1. Locate CURRENT derived artifacts only
# ------------------------------------------------------------------------------

gate1_file <- first_existing(c(
  file.path(
    PROJECT_DIR,
    "data",
    "derived",
    "gate1_audit_v2_1",
    "spine_gate1_v2_patient_level_2020_2024.rds"
  ),
  file.path(
    PROJECT_DIR,
    "data",
    "derived",
    "gate1_audit_v2_1",
    "spine_gate1_v2_1_patient_level_2020_2024.rds"
  ),
  file.path(
    PROJECT_DIR,
    "data",
    "derived",
    "gate1_audit_v2_1",
    "spine_gate1_patient_level_2020_2024.rds"
  )
))
if (is.na(gate1_file)) {
  gate1_file <- latest_recursive(
    file.path(PROJECT_DIR, "data", "derived"),
    "spine_gate1.*patient_level.*2020_2024\\.rds$"
  )
}

proc_file <- latest_recursive(
  file.path(PROJECT_DIR, "data", "derived"),
  "spine_gate1.*matched.*operative.*2020_2024\\.rds$"
)

current_data_file <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication",
  "spine_tract_frozen_primary_2020_2024.rds"
)
current_spec_file <- file.path(
  PROJECT_DIR,
  "data",
  "derived",
  "stage2_predictor_adjudication",
  "spine_tract_predictor_sets.rds"
)
current_cfg_file <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification",
  "rolling_origin_gpu",
  "14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv"
)

historical_map_file <- file.path(
  PROJECT_DIR,
  "reference",
  "old_spine_operative_icd10pcs_map_legacy.csv"
)
current_map_file <- file.path(
  PROJECT_DIR,
  "reference",
  "spine_operative_icd10pcs_map.csv"
)

needed <- c(
  gate1_file,
  proc_file,
  current_data_file,
  current_spec_file,
  current_cfg_file,
  historical_map_file,
  current_map_file
)

bad <- needed[
  is.na(needed) |
    !file.exists(needed)
]
if (length(bad) > 0L) {
  stop(
    "Missing required derived artifact(s):\n",
    paste(bad, collapse = "\n"),
    call. = FALSE
  )
}

manifest <- data.table(
  artifact = c(
    "Gate1 patient-level",
    "Gate1 matched operative procedures",
    "Current frozen Stage2 dataset",
    "Current predictor spec",
    "Current selected Stage3A configs",
    "Historical operative map",
    "Current operative map"
  ),
  path = needed
)
manifest[, size_mb := round(
  file.info(path)$size / 1024^2,
  2
)]
fwrite(
  manifest,
  file.path(OUT_DIR, "00_SOURCE_MANIFEST.csv")
)
print(manifest)

# ------------------------------------------------------------------------------
# 2. First answer: did the operative ICD-10-PCS map change?
# ------------------------------------------------------------------------------

old_map <- fread(historical_map_file)
new_map <- fread(current_map_file)

setnames(
  old_map,
  names(old_map),
  tolower(names(old_map))
)
setnames(
  new_map,
  names(new_map),
  tolower(names(new_map))
)

required_map_cols <- c(
  "procedure_code",
  "spine_operation_category"
)
if (
  !all(required_map_cols %in% names(old_map)) ||
    !all(required_map_cols %in% names(new_map))
) {
  stop("Unexpected operative-map columns.", call. = FALSE)
}

old_map[, procedure_code := toupper(
  gsub("[^A-Z0-9]", "", procedure_code)
)]
new_map[, procedure_code := toupper(
  gsub("[^A-Z0-9]", "", procedure_code)
)]

old_map <- unique(old_map[
  ,
  ..required_map_cols
])
new_map <- unique(new_map[
  ,
  ..required_map_cols
])

only_old <- fsetdiff(old_map, new_map)
only_new <- fsetdiff(new_map, old_map)

map_compare <- data.table(
  historical_rows = nrow(old_map),
  current_rows = nrow(new_map),
  identical_rows =
    nrow(only_old) == 0L &&
    nrow(only_new) == 0L,
  historical_only_rows = nrow(only_old),
  current_only_rows = nrow(only_new)
)

fwrite(
  map_compare,
  file.path(
    OUT_DIR,
    "01_OPERATIVE_MAP_COMPARISON.csv"
  )
)
fwrite(
  only_old,
  file.path(
    OUT_DIR,
    "01a_OPERATIVE_MAP_HISTORICAL_ONLY.csv"
  )
)
fwrite(
  only_new,
  file.path(
    OUT_DIR,
    "01b_OPERATIVE_MAP_CURRENT_ONLY.csv"
  )
)

cat("\nOperative code-map comparison:\n")
print(map_compare)

# ------------------------------------------------------------------------------
# 3. Load current frozen direct-presentation dataset
# ------------------------------------------------------------------------------

cur <- as.data.table(
  readRDS(current_data_file)
)
spec <- readRDS(current_spec_file)
cfg <- fread(current_cfg_file)

if (nrow(cur) != 424017L) {
  stop(
    "Current frozen primary N changed: ",
    nrow(cur),
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
    "Expected 53 final_core predictors; found ",
    length(final_core),
    ".",
    call. = FALSE
  )
}

# Add old-style redundant admission-era constructs where available.
cur[, oldstyle_gcs_total := fifelse(
  !is.na(gcs_eye_clean) &
    !is.na(gcs_verbal_clean) &
    !is.na(gcs_motor_clean),
  gcs_eye_clean +
    gcs_verbal_clean +
    gcs_motor_clean,
  NA_real_
)]

if ("gcs_total_clean" %in% names(cur)) {
  cur[
    !is.na(gcs_total_clean),
    oldstyle_gcs_total :=
      as.numeric(gcs_total_clean)
  ]
}

cur[, oldstyle_gcs_severity := fifelse(
  is.na(oldstyle_gcs_total),
  "__UNKNOWN__",
  fifelse(
    oldstyle_gcs_total <= 8,
    "Severe_3_8",
    fifelse(
      oldstyle_gcs_total <= 12,
      "Moderate_9_12",
      "Mild_13_15"
    )
  )
)]

cur[, oldstyle_hypotension_sbp90 := fifelse(
  is.na(sbp_clean),
  NA_integer_,
  as.integer(sbp_clean <= 90)
)]
cur[, oldstyle_hypoxia_spo2_90 := fifelse(
  is.na(spo2_clean),
  NA_integer_,
  as.integer(spo2_clean <= 90)
)]
cur[, oldstyle_tachycardia_hr120 := fifelse(
  is.na(pulse_clean),
  NA_integer_,
  as.integer(pulse_clean >= 120)
)]
cur[, oldstyle_abnormal_rr := fifelse(
  is.na(rr_clean),
  NA_integer_,
  as.integer(
    rr_clean < 10 |
      rr_clean > 29
  )
)]

oldstyle_missing_source <- intersect(
  c(
    "gcs_eye_clean",
    "gcs_verbal_clean",
    "gcs_motor_clean",
    "oldstyle_gcs_total",
    "sbp_clean",
    "pulse_clean",
    "rr_clean",
    "spo2_clean"
  ),
  names(cur)
)
for (v in oldstyle_missing_source) {
  cur[
    ,
    (paste0(v, "_oldstyle_missing")) :=
      as.integer(is.na(get(v)))
  ]
}

# Maximal CURRENT admission-era candidate set: every prespecified candidate
# group, plus the old-style redundant GCS / physiology constructs.
max_current <- unique(unlist(
  pg,
  use.names = FALSE
))
max_current <- intersect(
  max_current,
  names(cur)
)

oldstyle_extras <- grep(
  "^oldstyle_|_oldstyle_missing$",
  names(cur),
  value = TRUE
)

max_current <- unique(c(
  max_current,
  oldstyle_extras
))

# Intermediate set: final_core + social/airway + old-style redundancies.
plus_redundant <- unique(c(
  final_core,
  intersect(
    c(
      pg$social_context,
      pg$airway_state,
      oldstyle_extras
    ),
    names(cur)
  )
))

predictor_sets <- list(
  current_final_core = final_core,
  current_plus_redundant =
    plus_redundant,
  current_maximal_admission =
    max_current
)

fwrite(
  rbindlist(lapply(
    names(predictor_sets),
    function(nm) {
      data.table(
        predictor_set = nm,
        predictor =
          predictor_sets[[nm]]
      )
    }
  )),
  file.path(
    OUT_DIR,
    "02_CURRENT_DIRECT_PREDICTOR_SETS.csv"
  )
)

# ------------------------------------------------------------------------------
# 4. Lightweight Gate1 transfer dataset
# ------------------------------------------------------------------------------

cat("\nLoading existing Gate1 RDS (no raw-file rebuild)...\n")
gate <- as.data.table(
  readRDS(gate1_file)
)

if (!"adult_narrow_spine" %in% names(gate)) {
  stop(
    "Gate1 RDS does not contain adult_narrow_spine.",
    call. = FALSE
  )
}

gate <- gate[
  adult_narrow_spine == TRUE &
    admission_year %in% 2020:2024
]

# Keep only variables needed for this forensic transfer analysis before merging.
keep_gate <- intersect(
  c(
    "inc_key",
    "admission_year",
    "age_years",
    "sex",
    "race",
    "ethnicity",
    "payer",
    "transfer_in",
    "direct_presentation_known",
    "transfer_status_unknown",
    "transfer_status",
    "gcs_eye_clean",
    "gcs_verbal_clean",
    "gcs_motor_clean",
    "gcs_total_clean",
    "gcs_qualifier",
    "sbp_clean",
    "hr_clean",
    "rr_clean",
    "spo2_clean",
    "dx_thoracic_vertebral_fracture",
    "dx_lumbar_vertebral_fracture",
    "dx_thoracic_spinal_cord_injury",
    "dx_lumbar_spinal_cord_injury",
    "dx_thoracic_nerve_root_injury",
    "dx_lumbar_sacral_nerve_root_injury",
    "dx_cauda_equina_injury",
    "dx_cervical_spine_injury",
    "dx_tbi_any"
  ),
  names(gate)
)
gate <- gate[
  ,
  ..keep_gate
]

proc <- as.data.table(
  readRDS(proc_file)
)
proc[, inc_key := as.character(inc_key)]
gate[, inc_key := as.character(inc_key)]

op <- unique(
  proc[
    ,
    .(
      inc_key,
      admission_year =
        as.integer(admission_year)
    )
  ]
)
op[, operative_spine_intervention := 1L]

gate <- merge(
  gate,
  op,
  by = c(
    "inc_key",
    "admission_year"
  ),
  all.x = TRUE,
  sort = FALSE
)
gate[
  is.na(operative_spine_intervention),
  operative_spine_intervention := 0L
]

rm(proc, op)
gc()

# Normalize categorical values without trying to reinterpret their clinical meaning.
clean_cat <- function(x) {
  z <- trimws(as.character(x))
  z[
    is.na(z) |
      z == ""
  ] <- "__UNKNOWN__"
  z
}

for (v in intersect(
  c(
    "sex",
    "race",
    "ethnicity",
    "payer",
    "gcs_qualifier",
    "transfer_status"
  ),
  names(gate)
)) {
  gate[
    ,
    (paste0(v, "_forensic")) :=
      clean_cat(get(v))
  ]
}

gate[, oldstyle_gcs_total := gcs_total_clean]
gate[
  is.na(oldstyle_gcs_total) &
    !is.na(gcs_eye_clean) &
    !is.na(gcs_verbal_clean) &
    !is.na(gcs_motor_clean),
  oldstyle_gcs_total :=
    gcs_eye_clean +
    gcs_verbal_clean +
    gcs_motor_clean
]

gate[, oldstyle_gcs_severity := fifelse(
  is.na(oldstyle_gcs_total),
  "__UNKNOWN__",
  fifelse(
    oldstyle_gcs_total <= 8,
    "Severe_3_8",
    fifelse(
      oldstyle_gcs_total <= 12,
      "Moderate_9_12",
      "Mild_13_15"
    )
  )
)]

gate[, oldstyle_hypotension_sbp90 := fifelse(
  is.na(sbp_clean),
  NA_integer_,
  as.integer(sbp_clean <= 90)
)]
gate[, oldstyle_hypoxia_spo2_90 := fifelse(
  is.na(spo2_clean),
  NA_integer_,
  as.integer(spo2_clean <= 90)
)]
gate[, oldstyle_tachycardia_hr120 := fifelse(
  is.na(hr_clean),
  NA_integer_,
  as.integer(hr_clean >= 120)
)]
gate[, oldstyle_abnormal_rr := fifelse(
  is.na(rr_clean),
  NA_integer_,
  as.integer(rr_clean < 10 | rr_clean > 29)
)]

gate_basic <- intersect(
  c(
    "age_years",
    "sex_forensic",
    "race_forensic",
    "ethnicity_forensic",
    "payer_forensic",
    "gcs_eye_clean",
    "gcs_verbal_clean",
    "gcs_motor_clean",
    "oldstyle_gcs_total",
    "oldstyle_gcs_severity",
    "gcs_qualifier_forensic",
    "sbp_clean",
    "hr_clean",
    "rr_clean",
    "spo2_clean",
    "oldstyle_hypotension_sbp90",
    "oldstyle_hypoxia_spo2_90",
    "oldstyle_tachycardia_hr120",
    "oldstyle_abnormal_rr",
    "dx_thoracic_vertebral_fracture",
    "dx_lumbar_vertebral_fracture",
    "dx_thoracic_spinal_cord_injury",
    "dx_lumbar_spinal_cord_injury",
    "dx_thoracic_nerve_root_injury",
    "dx_lumbar_sacral_nerve_root_injury",
    "dx_cauda_equina_injury",
    "dx_cervical_spine_injury",
    "dx_tbi_any"
  ),
  names(gate)
)

gate_with_transfer <- unique(c(
  gate_basic,
  intersect(
    c(
      "transfer_in",
      "transfer_status_forensic"
    ),
    names(gate)
  )
))

# Cohort summary itself can already tell us whether transfer-in strongly enriches
# operative treatment.
gate[, transfer_group := fifelse(
  transfer_in == TRUE,
  "Transfer-in",
  fifelse(
    direct_presentation_known == TRUE,
    "Direct",
    "Unknown"
  )
)]

transfer_summary <- gate[
  ,
  .(
    N = .N,
    operative_events =
      sum(operative_spine_intervention),
    operative_prevalence =
      mean(operative_spine_intervention)
  ),
  by = .(
    admission_year,
    transfer_group
  )
][
  order(
    admission_year,
    transfer_group
  )
]

fwrite(
  transfer_summary,
  file.path(
    OUT_DIR,
    "03_TRANSFER_OPERATION_PREVALENCE.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. Modeling helpers
# ------------------------------------------------------------------------------

is_num <- function(x) {
  is.numeric(x) ||
    is.integer(x)
}

fit_encoder <- function(d, predictors) {
  enc <- vector(
    "list",
    length(predictors)
  )
  names(enc) <- predictors

  for (v in predictors) {
    x <- d[[v]]

    if (is_num(x)) {
      enc[[v]] <- list(
        type = "numeric"
      )
    } else {
      z <- clean_cat(x)
      lev <- unique(c(
        sort(unique(z)),
        "__UNKNOWN__",
        "__OTHER__"
      ))
      enc[[v]] <- list(
        type = "categorical",
        levels = lev
      )
    }
  }

  enc
}

encode_with <- function(
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

    if (info$type == "numeric") {
      m <- Matrix(
        as.numeric(x),
        ncol = 1L,
        sparse = TRUE
      )
      colnames(m) <- v
      blocks[[i]] <- m
    } else {
      z <- clean_cat(x)
      z[
        !z %in% info$levels
      ] <- "__OTHER__"
      f <- factor(
        z,
        levels = info$levels
      )
      m <- sparseMatrix(
        i = seq_along(f),
        j = as.integer(f),
        x = 1,
        dims = c(
          length(f),
          length(info$levels)
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

  do.call(cbind, blocks)
}

fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  n1 <- as.numeric(
    sum(y == 1L)
  )
  n0 <- as.numeric(
    sum(y == 0L)
  )
  if (n1 == 0 || n0 == 0) {
    return(NA_real_)
  }

  r <- rank(
    p,
    ties.method = "average"
  )
  (
    sum(r[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}

fast_auprc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
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
  ys <- y[ord]
  precision <- cumsum(
    ys == 1L
  ) / seq_along(ys)
  mean(
    precision[ys == 1L]
  )
}

gpu_params <- function() {
  if (
    packageVersion("xgboost") <
      package_version("2.0.0")
  ) {
    list(
      tree_method = "gpu_hist",
      predictor = "gpu_predictor",
      gpu_id = 0L
    )
  } else {
    list(
      tree_method = "hist",
      device = "cuda:0"
    )
  }
}

cfg24 <- cfg[
  endpoint ==
    "operative_spine_intervention" &
    evaluation_year == 2024L
]
if (nrow(cfg24) != 1L) {
  stop(
    "Could not resolve current 2024 surgery structural configuration.",
    call. = FALSE
  )
}

params <- c(
  list(
    objective = "binary:logistic",
    eval_metric = "logloss",
    max_depth =
      as.integer(cfg24$max_depth),
    eta =
      as.numeric(cfg24$eta),
    min_child_weight =
      as.numeric(
        cfg24$min_child_weight
      ),
    subsample =
      as.numeric(cfg24$subsample),
    colsample_bytree =
      as.numeric(
        cfg24$colsample_bytree
      ),
    lambda = 1,
    alpha = 0,
    nthread = CPU_THREADS,
    verbosity = 0
  ),
  gpu_params()
)

train_es <- function(
    dtrain,
    dtune
) {
  if (
    packageVersion("xgboost") <
      package_version("2.0.0")
  ) {
    xgb.train(
      params = params,
      data = dtrain,
      nrounds = MAX_ROUNDS,
      watchlist = list(
        train = dtrain,
        tune = dtune
      ),
      early_stopping_rounds =
        EARLY_STOP,
      maximize = FALSE,
      verbose = 0
    )
  } else {
    xgb.train(
      params = params,
      data = dtrain,
      nrounds = MAX_ROUNDS,
      evals = list(
        train = dtrain,
        tune = dtune
      ),
      early_stopping_rounds =
        EARLY_STOP,
      maximize = FALSE,
      verbose = 0
    )
  }
}

extract_best <- function(m) {
  b <- tryCatch(
    m$best_iteration,
    error = function(e) NULL
  )
  if (
    !is.null(b) &&
      length(b) == 1L &&
      is.finite(b)
  ) {
    return(as.integer(b))
  }

  a <- tryCatch(
    xgb.attr(m, "best_iteration"),
    error = function(e) NULL
  )
  if (
    !is.null(a) &&
      length(a) == 1L
  ) {
    return(
      as.integer(a) + 1L
    )
  }

  MAX_ROUNDS
}

fit_temporal_2024 <- function(
    d,
    predictors,
    experiment
) {
  predictors <- intersect(
    unique(predictors),
    names(d)
  )

  x <- d[
    admission_year %in% 2020:2024 &
      !is.na(
        operative_spine_intervention
      )
  ]

  train <- x[
    admission_year %in% 2020:2022
  ]
  tune <- x[
    admission_year == 2023L
  ]
  dev <- x[
    admission_year %in% 2020:2023
  ]
  eva <- x[
    admission_year == 2024L
  ]

  enc_train <- fit_encoder(
    train,
    predictors
  )

  X_train <- encode_with(
    train,
    predictors,
    enc_train
  )
  X_tune <- encode_with(
    tune,
    predictors,
    enc_train
  )

  dtrain <- xgb.DMatrix(
    X_train,
    label =
      train$operative_spine_intervention,
    missing = NA_real_
  )
  dtune <- xgb.DMatrix(
    X_tune,
    label =
      tune$operative_spine_intervention,
    missing = NA_real_
  )

  set.seed(SEED)
  es <- train_es(
    dtrain,
    dtune
  )
  bi <- extract_best(es)

  rm(
    train,
    tune,
    X_train,
    X_tune,
    dtrain,
    dtune,
    es,
    enc_train
  )
  gc()

  enc_dev <- fit_encoder(
    dev,
    predictors
  )
  X_dev <- encode_with(
    dev,
    predictors,
    enc_dev
  )
  X_eval <- encode_with(
    eva,
    predictors,
    enc_dev
  )

  ddev <- xgb.DMatrix(
    X_dev,
    label =
      dev$operative_spine_intervention,
    missing = NA_real_
  )
  deval <- xgb.DMatrix(
    X_eval,
    label =
      eva$operative_spine_intervention,
    missing = NA_real_
  )

  set.seed(SEED + 1000L)
  fit <- xgb.train(
    params = params,
    data = ddev,
    nrounds = bi,
    verbose = 0
  )

  p <- as.numeric(
    predict(
      fit,
      deval
    )
  )
  y <- as.integer(
    eva$operative_spine_intervention
  )

  pred <- data.table(
    inc_key = eva$inc_key,
    admission_year =
      eva$admission_year,
    y = y,
    pred = p
  )
  fwrite(
    pred,
    file.path(
      PRED_DIR,
      paste0(
        experiment,
        "__2024.csv"
      )
    )
  )

  prev <- mean(y)
  auprc <- fast_auprc(y, p)

  out <- data.table(
    experiment = experiment,
    n_predictors =
      length(predictors),
    n_2024 = length(y),
    events_2024 = sum(y),
    prevalence_2024 = prev,
    auroc_2024 =
      fast_auc(y, p),
    auprc_2024 = auprc,
    pr_lift_2024 =
      auprc / prev,
    selected_rounds = bi
  )

  rm(
    dev,
    eva,
    X_dev,
    X_eval,
    ddev,
    deval,
    fit,
    pred,
    enc_dev
  )
  gc()

  out
}

# ------------------------------------------------------------------------------
# 6. Run only the models that answer the forensic question
# ------------------------------------------------------------------------------

results <- list()

run_case <- function(
    d,
    predictors,
    label
) {
  cat(
    "\n============================================================\n",
    label,
    "\n============================================================\n",
    sep = ""
  )

  t0 <- Sys.time()
  z <- fit_temporal_2024(
    d,
    predictors,
    label
  )
  z[, runtime_minutes := as.numeric(
    difftime(
      Sys.time(),
      t0,
      units = "mins"
    )
  )]
  results[[length(results) + 1L]] <<- z

  fwrite(
    rbindlist(
      results,
      fill = TRUE
    ),
    file.path(
      OUT_DIR,
      "04_FORENSIC_MODEL_RESULTS_CHECKPOINT.csv"
    )
  )

  print(z)
  invisible(z)
}

# Current direct target: verify baseline then add back removed/redundant families.
run_case(
  cur,
  final_core,
  "A_current_direct_final_core"
)

run_case(
  cur,
  plus_redundant,
  "B_current_direct_plus_redundant"
)

run_case(
  cur,
  max_current,
  "C_current_direct_maximal_admission"
)

# Transfer decomposition using only variables already in Gate1.
known_transfer <- gate[
  transfer_group %in%
    c("Direct", "Transfer-in")
]

run_case(
  known_transfer,
  gate_with_transfer,
  "D_all_known_transfer_with_transfer_indicator"
)

run_case(
  known_transfer,
  gate_basic,
  "E_all_known_transfer_without_transfer_indicator"
)

run_case(
  gate[
    transfer_group == "Direct"
  ],
  gate_basic,
  "F_gate1_direct_only_basic_oldstyle"
)

run_case(
  gate[
    transfer_group == "Transfer-in"
  ],
  gate_basic,
  "G_gate1_transfer_in_only_basic_oldstyle"
)

final <- rbindlist(
  results,
  fill = TRUE
)

# Add known historical reference rows for direct visual comparison only.
historical_reference <- data.table(
  experiment = c(
    "REFERENCE_old_random_validation",
    "REFERENCE_old_2024_temporal"
  ),
  n_predictors = NA_integer_,
  n_2024 = NA_integer_,
  events_2024 = NA_integer_,
  prevalence_2024 = NA_real_,
  auroc_2024 = c(
    0.751,
    0.736
  ),
  auprc_2024 = NA_real_,
  pr_lift_2024 = NA_real_,
  selected_rounds = NA_integer_,
  runtime_minutes = NA_real_
)

final_with_reference <- rbind(
  final,
  historical_reference,
  fill = TRUE
)

fwrite(
  final_with_reference,
  file.path(
    OUT_DIR,
    "05_FINAL_FORENSIC_MODEL_RESULTS.csv"
  )
)

# Direct contrasts that answer the question.
auc <- function(id) {
  z <- final[
    experiment == id,
    auroc_2024
  ]
  if (length(z) == 0L) return(NA_real_)
  z[1L]
}

contrasts <- data.table(
  contrast = c(
    "Add redundant/contextual current features to final_core",
    "Maximal current admission-era set vs final_core",
    "Transfer indicator contribution in mixed direct+transfer population",
    "Direct-only vs transfer-in-only discrimination under same basic feature set"
  ),
  reference_auroc = c(
    auc(
      "A_current_direct_final_core"
    ),
    auc(
      "A_current_direct_final_core"
    ),
    auc(
      "E_all_known_transfer_without_transfer_indicator"
    ),
    auc(
      "F_gate1_direct_only_basic_oldstyle"
    )
  ),
  comparator_auroc = c(
    auc(
      "B_current_direct_plus_redundant"
    ),
    auc(
      "C_current_direct_maximal_admission"
    ),
    auc(
      "D_all_known_transfer_with_transfer_indicator"
    ),
    auc(
      "G_gate1_transfer_in_only_basic_oldstyle"
    )
  )
)
contrasts[
  ,
  delta_auroc :=
    comparator_auroc -
    reference_auroc
]

fwrite(
  contrasts,
  file.path(
    OUT_DIR,
    "06_FORENSIC_AUROC_CONTRASTS.csv"
  )
)

cat("\n============================================================\n")
cat("LEAN SURGERY FORENSIC AUDIT COMPLETE\n")
cat("============================================================\n")
print(final_with_reference)
cat("\nContrasts:\n")
print(contrasts)
cat("\nTransfer prevalence:\n")
print(
  transfer_summary[
    admission_year == 2024L
  ]
)

writeLines(
  c(
    "SPINE-TRACT SURGERY FORENSIC AUDIT v4 LEAN",
    "",
    "No raw TQIP rebuild was performed.",
    "",
    "Interpretation:",
    "1. If 01_OPERATIVE_MAP_COMPARISON shows identical_rows=TRUE, the operative code map did not cause the AUROC decline.",
    "2. A vs B/C quantifies how much current predictor cleanup reduced discrimination in direct presenters.",
    "3. D vs E quantifies the transfer-indicator signal in a mixed population.",
    "4. F vs G shows whether direct and transfer-in patients are intrinsically different prediction problems.",
    "5. 03_TRANSFER_OPERATION_PREVALENCE shows whether transfer-ins are strongly enriched for operation.",
    "",
    "Do not automatically promote the maximal forensic feature set into the final calculator."
  ),
  file.path(
    OUT_DIR,
    "README_RESULTS.txt"
  )
)
