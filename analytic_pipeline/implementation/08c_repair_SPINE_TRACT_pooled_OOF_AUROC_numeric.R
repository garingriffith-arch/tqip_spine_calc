# ==============================================================================
# 08c_repair_SPINE_TRACT_pooled_OOF_AUROC_numeric.R
# Repair ONLY the pooled 2022-2024 AUROC values in Stage 3D.
#
# Why:
#   In the original fast_auc(), n1 and n0 were integers. With ~265k pooled
#   observations, n1*n0 can exceed R's 32-bit integer limit, producing NA/NaN.
#   The year-specific 2024 AUROCs and bootstrap CIs are unaffected.
#
# NO MODEL REFIT. NO BOOTSTRAP RERUN.
# ==============================================================================

suppressPackageStartupMessages(library(data.table))

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

CLASS_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification",
  "rolling_origin_gpu"
)

OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_manuscript_outputs",
  "frozen_predictions"
)

OUT_FILE <- file.path(
  OUT_DIR,
  "02_CLASSIFICATION_POOLED_OOF_POINT_ESTIMATES.csv"
)

fast_auc_safe <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  n1 <- as.numeric(sum(y == 1L))
  n0 <- as.numeric(sum(y == 0L))

  if (n1 == 0 || n0 == 0) return(NA_real_)

  r <- rank(p, ties.method = "average")

  (
    sum(r[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}

load_output <- function(year, id) {
  pred_dir <- file.path(CLASS_DIR, "predictions")

  if (id %in% c("home", "postacute", "death")) {
    d <- fread(file.path(
      pred_dir,
      paste0("discharge3_final__eval", year, "__predictions.csv")
    ))

    if (id == "home") {
      return(data.table(
        y = as.integer(d$observed_class == "Home/home health"),
        p = d$p_Home.home.health
      ))
    }
    if (id == "postacute") {
      return(data.table(
        y = as.integer(d$observed_class == "Post-acute facility"),
        p = d$p_Post.acute.facility
      ))
    }
    return(data.table(
      y = as.integer(d$observed_class == "Death/hospice"),
      p = d$p_Death.hospice
    ))
  }

  if (id == "hlos20") {
    d <- fread(file.path(
      pred_dir,
      paste0("hlos3_final__eval", year, "__predictions.csv")
    ))
    return(data.table(
      y = as.integer(d$observed_class == ">=20 days"),
      p = d$p_X..20.days
    ))
  }

  if (id %in% c("icu_any", "icu8")) {
    d <- fread(file.path(
      pred_dir,
      paste0("icu3_final__eval", year, "__predictions.csv")
    ))
    if (id == "icu_any") {
      return(data.table(
        y = as.integer(d$observed_class != "No ICU"),
        p = d$p_ICU.1.7.days + d$p_ICU...8.days
      ))
    }
    return(data.table(
      y = as.integer(d$observed_class == "ICU >=8 days"),
      p = d$p_ICU...8.days
    ))
  }

  if (id %in% c("vent_any", "vent8")) {
    d <- fread(file.path(
      pred_dir,
      paste0("vent3_final__eval", year, "__predictions.csv")
    ))
    if (id == "vent_any") {
      return(data.table(
        y = as.integer(d$observed_class != "No ventilation"),
        p = d$p_Ventilation.1.7.days + d$p_Ventilation...8.days
      ))
    }
    return(data.table(
      y = as.integer(d$observed_class == "Ventilation >=8 days"),
      p = d$p_Ventilation...8.days
    ))
  }

  if (id == "surgery") {
    d <- fread(file.path(
      pred_dir,
      paste0(
        "operative_spine_intervention__eval",
        year,
        "__predictions.csv"
      )
    ))
    return(data.table(
      y = as.integer(d$y),
      p = as.numeric(d$pred)
    ))
  }

  stop("Unknown id: ", id, call. = FALSE)
}

map <- data.table(
  id = c(
    "home", "postacute", "death", "hlos20",
    "icu_any", "icu8", "vent_any", "vent8", "surgery"
  ),
  output_id = c(
    "disposition_home",
    "disposition_postacute",
    "disposition_death_hospice",
    "hlos_ge20",
    "icu_any",
    "icu_ge8",
    "vent_any",
    "vent_ge8",
    "recorded_tl_operation"
  )
)

out <- fread(OUT_FILE)

# If every original pooled AUROC was NaN, fread may infer this column as logical.
# Force numeric storage before assigning repaired numeric AUROC values.
out[, auroc := as.numeric(auroc)]

for (i in seq_len(nrow(map))) {
  pooled <- rbindlist(lapply(
    2022:2024,
    function(y) load_output(y, map$id[i])
  ))

  auc <- fast_auc_safe(
    pooled$y,
    pooled$p
  )

  out[
    output_id == map$output_id[i],
    auroc := auc
  ]
}

if (any(!is.finite(out$auroc))) {
  stop(
    "At least one pooled AUROC remains non-finite after repair.",
    call. = FALSE
  )
}

fwrite(out, OUT_FILE)

cat("\nPOOLED OOF AUROC REPAIR COMPLETE\n")
cat("No model refit. No bootstrap rerun.\n")
cat("Updated:\n  ", OUT_FILE, "\n\n", sep = "")
print(
  out[, .(
    output_label,
    n,
    prevalence,
    auroc,
    auprc,
    calibration_slope
  )]
)
