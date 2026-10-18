# ==============================================================================
# 02c_repair_SPINE_TRACT_frozen_dx_phenotypes_v1.R
# SPINE-TRACT Stage 2A diagnosis-phenotype repair
#
# PURPOSE
#   Repair the canonical frozen Stage 2A datasets after QC showed that several
#   diagnosis-derived predictors were stored as 100% NA despite the cohort being
#   defined by those diagnoses.
#
# IMPORTANT
#   - Reads ONLY the five annual PUF_ICDDIAGNOSIS.csv files.
#   - Does NOT rerun Gate 1, PMHx, procedure endpoint derivation, or any model.
#   - Uses the exact ICD-10-CM definitions from Stage 2A v1.2.
#   - Backs up the current canonical frozen RDS files before replacing them.
#   - Hard-checks the repaired primary cohort against known Gate 1 phenotype
#     counts before writing canonical replacements.
#
# AFTER THIS SCRIPT
#   Do NOT reuse old model performance as final. The prior models were trained
#   with several key diagnosis predictors entirely missing. Targeted
#   re-adjudication/refitting is required after this repair is confirmed.
# ==============================================================================

rm(list = ls())
gc()

required_pkgs <- c("data.table")
missing_pkgs <- required_pkgs[
  !vapply(
    required_pkgs,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]
if (length(missing_pkgs) > 0L) {
  stop(
    "Missing package(s): ",
    paste(
      missing_pkgs,
      collapse = ", "
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

YEARS <- 2020:2024

RAW_ROOT <- Sys.getenv(
  "SPINE_TRACT_RAW_ROOT",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/Griffith_G-selected"
)

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

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage2_dx_repair_v1"
)

dir.create(
  OUT_DIR,
  recursive = TRUE,
  showWarnings = FALSE
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

needed <- c(
  PRIMARY_FILE,
  UNION_FILE,
  SPEC_FILE
)

bad <- needed[
  !file.exists(needed)
]

if (length(bad) > 0L) {
  stop(
    "Missing required Stage 2 artifact(s):\n",
    paste(
      bad,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

EXPECTED_PRIMARY_N <- 424017L
EXPECTED_UNION_N <- 440380L

TARGET_CPU_FRACTION <- 0.80
LOGICAL_CORES <- parallel::detectCores(
  logical = TRUE
)

if (
  is.na(LOGICAL_CORES) ||
    LOGICAL_CORES < 1L
) {
  LOGICAL_CORES <- 1L
}

N_THREADS <- max(
  1L,
  min(
    LOGICAL_CORES,
    as.integer(
      round(
        LOGICAL_CORES *
          TARGET_CPU_FRACTION
      )
    )
  )
)

data.table::setDTthreads(
  N_THREADS
)

# ------------------------------------------------------------------------------
# 1. Helpers copied from Stage 2A v1.2
# ------------------------------------------------------------------------------

norm_key <- function(x) {
  z <- trimws(
    as.character(x)
  )

  z <- gsub(
    "\\.0$",
    "",
    z
  )

  z[
    z %in% c(
      "",
      "NA",
      "NaN",
      "NULL",
      "null"
    )
  ] <- NA_character_

  z
}

norm_code <- function(x) {
  z <- toupper(
    as.character(x)
  )

  z <- gsub(
    "[^A-Z0-9]",
    "",
    z
  )

  z[
    z %in% c(
      "",
      "NA",
      "NULL"
    )
  ] <- NA_character_

  z
}

starts_any <- function(
    x,
    prefixes
) {
  x2 <- x
  x2[
    is.na(x2)
  ] <- ""

  out <- rep(
    FALSE,
    length(x2)
  )

  for (p in prefixes) {
    out <- out |
      startsWith(
        x2,
        p
      )
  }

  out
}

annual_csv_dir <- function(year) {
  exact <- file.path(
    RAW_ROOT,
    sprintf(
      "PUF AY %d",
      year
    ),
    sprintf(
      "PUF AY %d",
      year
    ),
    "CSV"
  )

  if (dir.exists(exact)) {
    return(exact)
  }

  year_root <- file.path(
    RAW_ROOT,
    sprintf(
      "PUF AY %d",
      year
    )
  )

  if (!dir.exists(year_root)) {
    stop(
      "Annual folder not found: ",
      year_root,
      call. = FALSE
    )
  }

  csvs <- list.files(
    year_root,
    pattern = "\\.csv$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  dirs <- unique(
    dirname(csvs)
  )

  if (length(dirs) == 1L) {
    return(dirs)
  }

  stop(
    "Could not uniquely resolve CSV directory for ",
    year,
    ". Candidates: ",
    paste(
      dirs,
      collapse = " | "
    ),
    call. = FALSE
  )
}

find_annual_file <- function(
    year,
    basename_target
) {
  d <- annual_csv_dir(year)

  direct <- file.path(
    d,
    basename_target
  )

  if (file.exists(direct)) {
    return(
      normalizePath(
        direct,
        winslash = "/",
        mustWork = TRUE
      )
    )
  }

  hits <- list.files(
    d,
    pattern = paste0(
      "^",
      gsub(
        "\\.",
        "\\\\.",
        basename_target
      ),
      "$"
    ),
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(hits) != 1L) {
    stop(
      "Expected exactly one ",
      basename_target,
      " for ",
      year,
      "; found ",
      length(hits),
      call. = FALSE
    )
  }

  normalizePath(
    hits[1L],
    winslash = "/",
    mustWork = TRUE
  )
}

read_header <- function(path) {
  names(
    fread(
      path,
      nrows = 0,
      showProgress = FALSE,
      check.names = FALSE
    )
  )
}

canon <- function(x) {
  z <- toupper(
    as.character(x)
  )

  z <- gsub(
    "[^A-Z0-9]+",
    "_",
    z
  )

  gsub(
    "^_+|_+$",
    "",
    z
  )
}

detect_col <- function(
    nms,
    exact = character(),
    regex = character(),
    exclude = character()
) {
  cn <- canon(nms)

  for (e in exact) {
    idx <- which(
      cn == canon(e)
    )

    if (length(idx) > 0L) {
      return(
        nms[
          idx[1L]
        ]
      )
    }
  }

  keep <- rep(
    FALSE,
    length(nms)
  )

  for (r in regex) {
    keep <- keep |
      grepl(
        r,
        cn,
        ignore.case = TRUE,
        perl = TRUE
      )
  }

  if (length(exclude) > 0L) {
    for (r in exclude) {
      keep <- keep &
        !grepl(
          r,
          cn,
          ignore.case = TRUE,
          perl = TRUE
        )
    }
  }

  idx <- which(keep)

  if (length(idx) == 0L) {
    NA_character_
  } else {
    nms[
      idx[1L]
    ]
  }
}

# ------------------------------------------------------------------------------
# 2. Load canonical frozen union and audit BEFORE repair
# ------------------------------------------------------------------------------

union <- as.data.table(
  readRDS(
    UNION_FILE
  )
)

union[
  ,
  inc_key :=
    norm_key(
      inc_key
    )
]

union[
  ,
  admission_year :=
    as.integer(
      admission_year
    )
]

if (nrow(union) != EXPECTED_UNION_N) {
  stop(
    "Expected frozen union N=",
    EXPECTED_UNION_N,
    "; observed ",
    nrow(union),
    ".",
    call. = FALSE
  )
}

if (
  sum(
    union$
      primary_complete_trajectory ==
      TRUE,
    na.rm = TRUE
  ) !=
    EXPECTED_PRIMARY_N
) {
  stop(
    "Expected primary flag N=",
    EXPECTED_PRIMARY_N,
    " within frozen union.",
    call. = FALSE
  )
}

# Save a manifest of the exact input artifacts before mutation.
source_manifest <- data.table(
  artifact = c(
    "Canonical frozen primary before repair",
    "Canonical frozen union before repair",
    "Predictor specification"
  ),
  path = c(
    PRIMARY_FILE,
    UNION_FILE,
    SPEC_FILE
  ),
  size_mb = round(
    file.info(
      c(
        PRIMARY_FILE,
        UNION_FILE,
        SPEC_FILE
      )
    )$size /
      1024^2,
    2
  ),
  md5 = unname(
    tools::md5sum(
      c(
        PRIMARY_FILE,
        UNION_FILE,
        SPEC_FILE
      )
    )
  )
)

fwrite(
  source_manifest,
  file.path(
    OUT_DIR,
    "00_SOURCE_ARTIFACTS_BEFORE_REPAIR.csv"
  )
)

# ------------------------------------------------------------------------------
# 3. Rebuild exact current diagnosis phenotypes for all union patients
# ------------------------------------------------------------------------------

dx_flags_list <- list()
dx_file_audit <- list()

for (yr in YEARS) {
  cat(
    "\nRebuilding diagnosis flags for ",
    yr,
    "...\n",
    sep = ""
  )

  f <- find_annual_file(
    yr,
    "PUF_ICDDIAGNOSIS.csv"
  )

  hdr <- read_header(f)

  id_col <- detect_col(
    hdr,
    exact = c(
      "INC_KEY",
      "INCKEY"
    ),
    regex = c(
      "^INC.*KEY$",
      "TRAUMA.*KEY",
      "PATIENT.*ID"
    )
  )

  code_col <- detect_col(
    hdr,
    exact = c(
      "ICDDIAGNOSISCODE",
      "ICD10CM",
      "ICDDIAGNOSIS"
    ),
    regex = c(
      "ICD.*DIAG",
      "DIAG.*CODE",
      "ICD10.*CM"
    ),
    exclude = c(
      "AIS",
      "DESC",
      "LOOKUP",
      "VERSION"
    )
  )

  if (
    is.na(id_col) ||
      is.na(code_col)
  ) {
    stop(
      "Could not identify diagnosis columns for ",
      yr,
      call. = FALSE
    )
  }

  dx_file_audit[[length(dx_file_audit) + 1L]] <- data.table(
    admission_year = yr,
    path = f,
    id_col = id_col,
    code_col = code_col
  )

  x <- fread(
    f,
    select = c(
      id_col,
      code_col
    ),
    integer64 = "character",
    showProgress = TRUE,
    nThread = N_THREADS
  )

  setnames(
    x,
    c(
      id_col,
      code_col
    ),
    c(
      "inc_key",
      "code"
    )
  )

  x[
    ,
    inc_key :=
      norm_key(
        inc_key
      )
  ]

  x[
    ,
    code :=
      norm_code(
        code
      )
  ]

  ids <- union[
    admission_year == yr,
    unique(inc_key)
  ]

  x <- x[
    inc_key %chin%
      ids &
      !is.na(code)
  ]

  x[, `:=`(
    # Thoracolumbar primary phenotype and modifiers
    dx_thoracic_vertebral_fracture_tmp =
      starts_any(
        code,
        c(
          "S220",
          "S221"
        )
      ),
    dx_lumbar_vertebral_fracture_tmp =
      starts_any(
        code,
        c("S320")
      ),
    dx_thoracic_spinal_cord_injury_tmp =
      starts_any(
        code,
        c(
          "S240",
          "S241"
        )
      ),
    dx_lumbar_spinal_cord_injury_tmp =
      starts_any(
        code,
        c(
          "S340",
          "S341"
        )
      ),
    dx_thoracic_nerve_root_injury_tmp =
      starts_any(
        code,
        c("S242")
      ),
    dx_lumbar_sacral_nerve_root_injury_tmp =
      starts_any(
        code,
        c("S342")
      ),
    dx_cauda_equina_injury_tmp =
      starts_any(
        code,
        c("S343")
      ),
    dx_lumbosacral_plexus_injury_tmp =
      starts_any(
        code,
        c("S344")
      ),
    dx_thoracic_ligament_dislocation_sprain_tmp =
      starts_any(
        code,
        c(
          "S230",
          "S231",
          "S232",
          "S233",
          "S234"
        )
      ),
    dx_lumbar_ligament_dislocation_sprain_tmp =
      starts_any(
        code,
        c(
          "S330",
          "S331",
          "S332",
          "S333",
          "S334",
          "S335"
        )
      ),
    dx_cervical_spine_injury_tmp =
      starts_any(
        code,
        c(
          "S12",
          "S13",
          "S14"
        )
      ),

    # Head / neck associated injury
    dx_tbi_any_tmp =
      startsWith(
        code,
        "S06"
      ),
    dx_intracranial_hemorrhage_any_tmp =
      starts_any(
        code,
        c(
          "S063",
          "S064",
          "S065",
          "S066"
        )
      ),
    dx_skull_fracture_any_tmp =
      startsWith(
        code,
        "S02"
      ),
    dx_neck_vascular_injury_tmp =
      startsWith(
        code,
        "S15"
      ),

    # Chest detail
    dx_rib_fracture_any_tmp =
      starts_any(
        code,
        c(
          "S223",
          "S224",
          "S225"
        )
      ),
    dx_multiple_rib_fractures_tmp =
      startsWith(
        code,
        "S224"
      ),
    dx_flail_chest_tmp =
      startsWith(
        code,
        "S225"
      ),
    dx_sternum_fracture_tmp =
      startsWith(
        code,
        "S222"
      ),
    dx_pulmonary_contusion_or_laceration_tmp =
      startsWith(
        code,
        "S273"
      ),
    dx_traumatic_pneumothorax_tmp =
      startsWith(
        code,
        "S270"
      ),
    dx_traumatic_hemothorax_tmp =
      startsWith(
        code,
        "S271"
      ),
    dx_traumatic_hemopneumothorax_tmp =
      startsWith(
        code,
        "S272"
      ),
    dx_diaphragm_injury_tmp =
      startsWith(
        code,
        "S278"
      ),
    dx_thoracic_vascular_injury_tmp =
      startsWith(
        code,
        "S25"
      ),
    dx_cardiac_injury_tmp =
      startsWith(
        code,
        "S26"
      ),

    # Other associated injury
    dx_sacral_coccyx_fracture_tmp =
      starts_any(
        code,
        c(
          "S321",
          "S322"
        )
      ),
    dx_pelvic_fracture_tmp =
      starts_any(
        code,
        c(
          "S323",
          "S324",
          "S325",
          "S327",
          "S328"
        )
      ),
    dx_abdominal_visceral_vascular_injury_tmp =
      starts_any(
        code,
        c(
          "S30",
          "S31",
          "S35",
          "S36",
          "S37",
          "S38",
          "S39"
        )
      ),
    dx_upper_extremity_injury_tmp =
      grepl(
        "^S[456][0-9]",
        code
      ),
    dx_lower_extremity_injury_tmp =
      grepl(
        "^S[789][0-9]",
        code
      )
  )]

  tmp_cols <- grep(
    "_tmp$",
    names(x),
    value = TRUE
  )

  sx <- x[
    ,
    lapply(
      .SD,
      function(v) {
        as.integer(
          any(
            v,
            na.rm = TRUE
          )
        )
      }
    ),
    by = inc_key,
    .SDcols = tmp_cols
  ]

  setnames(
    sx,
    tmp_cols,
    sub(
      "_tmp$",
      "",
      tmp_cols
    )
  )

  sx[, `:=`(
    dx_tl_spinal_cord_injury_any =
      as.integer(
        dx_thoracic_spinal_cord_injury ==
          1L |
          dx_lumbar_spinal_cord_injury ==
            1L
      ),
    dx_tl_root_cauda_any =
      as.integer(
        dx_thoracic_nerve_root_injury ==
          1L |
          dx_lumbar_sacral_nerve_root_injury ==
            1L |
          dx_cauda_equina_injury ==
            1L
      ),
    dx_tl_soft_tissue_any =
      as.integer(
        dx_thoracic_ligament_dislocation_sprain ==
          1L |
          dx_lumbar_ligament_dislocation_sprain ==
            1L
      ),
    dx_pleural_air_blood_injury =
      as.integer(
        dx_traumatic_pneumothorax ==
          1L |
          dx_traumatic_hemothorax ==
            1L |
          dx_traumatic_hemopneumothorax ==
            1L
      )
  )]

  sx[
    ,
    dx_nonspinal_chest_injury :=
      as.integer(
        dx_rib_fracture_any ==
          1L |
          dx_sternum_fracture ==
            1L |
          dx_pulmonary_contusion_or_laceration ==
            1L |
          dx_pleural_air_blood_injury ==
            1L |
          dx_diaphragm_injury ==
            1L |
          dx_thoracic_vascular_injury ==
            1L |
          dx_cardiac_injury ==
            1L
      )
  ]

  sx[
    ,
    admission_year := yr
  ]

  dx_flags_list[[length(dx_flags_list) + 1L]] <- sx

  rm(
    x,
    sx
  )
  gc()
}

dx_file_audit <- rbindlist(
  dx_file_audit,
  fill = TRUE
)

fwrite(
  dx_file_audit,
  file.path(
    OUT_DIR,
    "01_DIAGNOSIS_SOURCE_AUDIT.csv"
  )
)

dx_flags <- rbindlist(
  dx_flags_list,
  fill = TRUE
)

derived_dx_cols <- setdiff(
  names(dx_flags),
  c(
    "admission_year",
    "inc_key"
  )
)

# ------------------------------------------------------------------------------
# 4. Before / after audit and safe replacement
# ------------------------------------------------------------------------------

before_audit <- rbindlist(
  lapply(
    derived_dx_cols,
    function(v) {
      if (!v %in% names(union)) {
        return(
          data.table(
            variable = v,
            before_missing_n =
              nrow(union),
            before_positive_n =
              NA_integer_
          )
        )
      }

      data.table(
        variable = v,
        before_missing_n =
          sum(
            is.na(
              union[[v]]
            )
          ),
        before_positive_n =
          sum(
            union[[v]] ==
              1L,
            na.rm = TRUE
          )
      )
    }
  ),
  fill = TRUE
)

# Remove the stale diagnosis-derived columns so the merge cannot suffix them.
drop_existing <- intersect(
  derived_dx_cols,
  names(union)
)

if (length(drop_existing) > 0L) {
  union[
    ,
    (drop_existing) := NULL
  ]
}

pre_n <- nrow(union)

union <- merge(
  union,
  dx_flags,
  by = c(
    "admission_year",
    "inc_key"
  ),
  all.x = TRUE,
  sort = FALSE
)

if (nrow(union) != pre_n) {
  stop(
    "Diagnosis repair merge changed union row count.",
    call. = FALSE
  )
}

# Patients with no matching diagnosis row get 0 for diagnosis-derived flags.
for (v in derived_dx_cols) {
  set(
    union,
    which(
      is.na(
        union[[v]]
      )
    ),
    v,
    0L
  )
}

setorder(
  union,
  admission_year,
  inc_key
)

primary <- union[
  primary_complete_trajectory ==
    TRUE
]

if (nrow(primary) != EXPECTED_PRIMARY_N) {
  stop(
    "Primary N changed after repair.",
    call. = FALSE
  )
}

if (nrow(union) != EXPECTED_UNION_N) {
  stop(
    "Union N changed after repair.",
    call. = FALSE
  )
}

after_audit <- rbindlist(
  lapply(
    derived_dx_cols,
    function(v) {
      data.table(
        variable = v,
        after_missing_n =
          sum(
            is.na(
              union[[v]]
            )
          ),
        after_positive_n =
          sum(
            union[[v]] ==
              1L,
            na.rm = TRUE
          ),
        primary_positive_n =
          sum(
            primary[[v]] ==
              1L,
            na.rm = TRUE
          )
      )
    }
  ),
  fill = TRUE
)

repair_audit <- merge(
  before_audit,
  after_audit,
  by = "variable",
  all = TRUE,
  sort = FALSE
)

fwrite(
  repair_audit,
  file.path(
    OUT_DIR,
    "02_DX_BEFORE_AFTER_REPAIR_AUDIT.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. Hard Gate 1 concordance checks for primary cohort
# ------------------------------------------------------------------------------

expected_primary_counts <- data.table(
  variable = c(
    "dx_thoracic_vertebral_fracture",
    "dx_lumbar_vertebral_fracture",
    "dx_thoracic_spinal_cord_injury",
    "dx_lumbar_spinal_cord_injury",
    "dx_thoracic_nerve_root_injury",
    "dx_lumbar_sacral_nerve_root_injury",
    "dx_cauda_equina_injury",
    "dx_tl_spinal_cord_injury_any",
    "dx_tl_root_cauda_any",
    "dx_cervical_spine_injury",
    "dx_tbi_any"
  ),
  expected_n = c(
    218002L,
    259812L,
    14713L,
    7825L,
    278L,
    896L,
    1104L,
    21661L,
    2215L,
    65376L,
    92755L
  )
)

expected_primary_counts[
  ,
  observed_n :=
    vapply(
      variable,
      function(v) {
        sum(
          primary[[v]] ==
            1L,
          na.rm = TRUE
        )
      },
      integer(1)
    )
]

expected_primary_counts[
  ,
  concordant :=
    observed_n ==
      expected_n
]

fwrite(
  expected_primary_counts,
  file.path(
    OUT_DIR,
    "03_GATE1_PRIMARY_PHENOTYPE_CONCORDANCE.csv"
  )
)

if (
  !all(
    expected_primary_counts$
      concordant
  )
) {
  print(
    expected_primary_counts
  )

  stop(
    "Rebuilt diagnosis phenotypes do not reproduce known Gate 1 primary counts. Canonical RDS files were NOT replaced. Review 03_GATE1_PRIMARY_PHENOTYPE_CONCORDANCE.csv.",
    call. = FALSE
  )
}

# Every primary patient must still meet the narrow phenotype.
primary[
  ,
  narrow_rederived :=
    dx_thoracic_vertebral_fracture ==
      1L |
      dx_lumbar_vertebral_fracture ==
        1L |
      dx_thoracic_spinal_cord_injury ==
        1L |
      dx_lumbar_spinal_cord_injury ==
        1L |
      dx_tl_root_cauda_any ==
        1L
]

if (
  any(
    !primary$
      narrow_rederived
  )
) {
  stop(
    "At least one primary patient no longer meets the rederived narrow phenotype. Canonical RDS files were NOT replaced.",
    call. = FALSE
  )
}

primary[
  ,
  narrow_rederived := NULL
]

# ------------------------------------------------------------------------------
# 6. Backup current canonical files, then write repaired canonical artifacts
# ------------------------------------------------------------------------------

stamp <- format(
  Sys.time(),
  "%Y%m%d_%H%M%S"
)

backup_dir <- file.path(
  DATA_DIR,
  paste0(
    "pre_dx_repair_backup_",
    stamp
  )
)

dir.create(
  backup_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

backup_primary <- file.path(
  backup_dir,
  basename(
    PRIMARY_FILE
  )
)

backup_union <- file.path(
  backup_dir,
  basename(
    UNION_FILE
  )
)

ok1 <- file.copy(
  PRIMARY_FILE,
  backup_primary,
  overwrite = FALSE
)

ok2 <- file.copy(
  UNION_FILE,
  backup_union,
  overwrite = FALSE
)

if (
  !isTRUE(ok1) ||
    !isTRUE(ok2)
) {
  stop(
    "Could not create canonical dataset backups. Repaired data were NOT written.",
    call. = FALSE
  )
}

# Write to temporary files first, then atomically replace.
tmp_primary <- paste0(
  PRIMARY_FILE,
  ".repair_tmp"
)

tmp_union <- paste0(
  UNION_FILE,
  ".repair_tmp"
)

saveRDS(
  primary,
  tmp_primary,
  compress = FALSE
)

saveRDS(
  union,
  tmp_union,
  compress = FALSE
)

# Re-read and verify before replacing canonical files.
verify_primary <- as.data.table(
  readRDS(
    tmp_primary
  )
)

verify_union <- as.data.table(
  readRDS(
    tmp_union
  )
)

if (
  nrow(verify_primary) !=
    EXPECTED_PRIMARY_N ||
    nrow(verify_union) !=
      EXPECTED_UNION_N
) {
  stop(
    "Temporary repaired RDS verification failed. Canonical files were NOT replaced.",
    call. = FALSE
  )
}

unlink(
  PRIMARY_FILE
)
unlink(
  UNION_FILE
)

ok_primary <- file.rename(
  tmp_primary,
  PRIMARY_FILE
)

ok_union <- file.rename(
  tmp_union,
  UNION_FILE
)

if (
  !isTRUE(ok_primary) ||
    !isTRUE(ok_union)
) {
  stop(
    "Canonical replacement failed after backup. Restore from: ",
    backup_dir,
    call. = FALSE
  )
}

# ------------------------------------------------------------------------------
# 7. Final fingerprint and QC
# ------------------------------------------------------------------------------

final_fingerprint <- data.table(
  artifact = c(
    "Repaired frozen primary dataset",
    "Repaired frozen primary+broad union",
    "Predictor-set specification"
  ),
  path = c(
    PRIMARY_FILE,
    UNION_FILE,
    SPEC_FILE
  )
)

final_fingerprint[
  ,
  size_mb :=
    round(
      file.info(path)$size /
        1024^2,
      2
    )
]

final_fingerprint[
  ,
  md5 :=
    unname(
      tools::md5sum(path)
    )
]

final_fingerprint[
  ,
  repaired_at :=
    as.character(
      Sys.time()
    )
]

fwrite(
  final_fingerprint,
  file.path(
    OUT_DIR,
    "04_REPAIRED_FROZEN_DATASET_FINGERPRINT.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT Stage 2A diagnosis repair complete.",
    "",
    paste0(
      "Primary N: ",
      format(
        nrow(primary),
        big.mark = ","
      )
    ),
    paste0(
      "Union N: ",
      format(
        nrow(union),
        big.mark = ","
      )
    ),
    "",
    paste0(
      "Backup directory: ",
      backup_dir
    ),
    "",
    "The canonical frozen Stage 2A primary and union RDS files now contain rederived ICD-10-CM diagnosis phenotypes using the Stage 2A v1.2 definitions.",
    "",
    "IMPORTANT: Prior Stage 2B/2C/3 model results are no longer final because several diagnosis predictors had been entirely missing. Do not delete the old outputs; they are useful for before/after comparison.",
    "",
    "Next step: targeted corrected predictor adjudication and rolling-origin refit. Do not rerun the entire historical pipeline."
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
  "STAGE 2A DX REPAIR COMPLETE\n"
)
cat(
  "============================================================\n"
)
cat(
  "Canonical primary and union RDS files were repaired.\n"
)
cat(
  "Backup:\n  ",
  backup_dir,
  "\n",
  sep = ""
)
cat(
  "\nGate 1 phenotype concordance:\n"
)
print(
  expected_primary_counts
)
cat(
  "\nReturn this output folder for review:\n  ",
  OUT_DIR,
  "\n",
  sep = ""
)
