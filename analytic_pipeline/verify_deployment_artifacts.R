# Verify published final deployment artifacts against Stage 21 manifest.
# No ACS source data are accessed, no model fitting or inference is performed.
# Run from the repository root:
#   Rscript analytic_pipeline/verify_deployment_artifacts.R

if (!requireNamespace("digest", quietly = TRUE)) {
  stop("Install the R package 'digest' to verify SHA-256 checksums.", call. = FALSE)
}
base <- file.path("data", "deployment_2020_2024")
manifest_path <- file.path(base, "artifact_manifest.csv")
if (!file.exists(manifest_path)) stop("Artifact manifest not found.", call. = FALSE)
a <- read.csv(manifest_path, stringsAsFactors = FALSE, check.names = FALSE)
expected_ids <- c("discharge3_final", "hlos3_final", "icu3_final",
                  "vent3_final", "operative_spine_intervention",
                  "hospital_los", "icu_los_conditional",
                  "ventilator_days_conditional")
need <- c("endpoint_id", "model_sha256", "meta_sha256",
          "selected_rounds", "training_n", "training_years",
          "internal_cv_folds", "model_size_bytes",
          "max_abs_logit_or_quantile_difference")
if (!all(need %in% names(a)) ||
    nrow(a) != 8L || anyDuplicated(a$endpoint_id) ||
    !setequal(a$endpoint_id, expected_ids)) {
  stop("Unexpected deployment artifact manifest.", call. = FALSE)
}
if (any(a$training_years != "2020-2024") ||
    any(a$internal_cv_folds != 3L) ||
    any(!is.finite(a$selected_rounds) | a$selected_rounds <= 0) ||
    any(!is.finite(a$training_n) | a$training_n <= 0) ||
    any(!is.finite(a$max_abs_logit_or_quantile_difference) |
        a$max_abs_logit_or_quantile_difference > 1e-6)) {
  stop("Deployment manifest training or parity QC failed.", call. = FALSE)
}
for (i in seq_len(nrow(a))) {
  id <- a$endpoint_id[i]
  files <- c(model = file.path(base, "models", paste0(id, ".json")),
             metadata = file.path(base, "encoders", paste0(id, "_encoder.rds")))
  if (any(!file.exists(files))) {
    stop("Missing deployment file(s) for ", id, ": ",
         paste(files[!file.exists(files)], collapse = ", "), call. = FALSE)
  }
  if (!identical(as.numeric(file.info(files[["model"]])$size),
                 as.numeric(a$model_size_bytes[i]))) {
    stop("Model size mismatch for ", id, call. = FALSE)
  }
  hashes <- vapply(files,
                   function(f) digest::digest(file = f, algo = "sha256",
                                              serialize = FALSE),
                   character(1))
  if (!identical(unname(hashes), c(a$model_sha256[i], a$meta_sha256[i]))) {
    stop("Deployment SHA-256 mismatch for ", id, call. = FALSE)
  }
}
cat("PASS: eight model JSONs, eight encoder RDS files, and all checksums match.\n")
cat("Training years, internal CV and Stage 21 serialization parity values verified.\n")
cat("Not checked: clinical performance, individual inference results or external validation.\n")
