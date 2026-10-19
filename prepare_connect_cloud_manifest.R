# Generate a Posit Connect Cloud manifest for the mixed R/Python Shiny app.
#
# Run this FROM A FRESH R SESSION with getwd() set to the GitHub clone root.
# Does not refit models or read any patient-level files. Only inspects the
# eight already-distributed JSON model files and eight RDS encoder filenames.
#
# Posit references:
#   https://rstudio.github.io/rsconnect/reference/writeManifest.html
#   https://forum.posit.co/t/publishing-content-that-uses-r-and-python-to-connect-cloud/210165
#
# Important: Do NOT put any local RETICULATE_PYTHON path into app_core.R,
# .Renviron, or a Git commit. The manifest contains version/dependency
# information; the host chooses its installed Python runtime.

required_r_packages <- c(
  "shiny", "bslib", "data.table", "Matrix", "xgboost", "reticulate",
  "rsconnect", "jsonlite"
)
missing_packages <- required_r_packages[
  !vapply(required_r_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop(
    "Install missing R packages, restart R, then re-run: ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
ids <- c(
  "discharge3_final",
  "hlos3_final",
  "icu3_final",
  "vent3_final",
  "operative_spine_intervention",
  "hospital_los",
  "icu_los_conditional",
  "ventilator_days_conditional"
)
models <- file.path("data", "deployment_2020_2024", "models", paste0(ids, ".json"))
encoders <- file.path("data", "deployment_2020_2024", "encoders",
                      paste0(ids, "_encoder.rds"))
app_files <- c(
  "app.R", "app_core.R", "requirements.txt", "www/ohsu_logo.png",
  "data/deployment_2020_2024/artifact_manifest.csv", models, encoders
)
missing_files <- app_files[!file.exists(file.path(root, app_files))]
if (length(missing_files)) {
  stop("Required deployment files not found in Git clone:\n",
       paste(missing_files, collapse = "\n"), call. = FALSE)
}
if (length(models) != 8L || length(encoders) != 8L) {
  stop("Incorrect model inventory", call. = FALSE)
}

# For LOCAL manifest generation only. Posit will supply its own interpreter.
py <- Sys.getenv("RETICULATE_PYTHON", unset = "")
if (!nzchar(py)) {
  py <- reticulate::conda_python("tbi-tract-xgb-gpu")
}
py <- normalizePath(py, winslash = "/", mustWork = TRUE)
if (!file.exists(py)) stop("Python executable not found: ", py)
Sys.setenv(RETICULATE_PYTHON = py)

# Fail early if the Conda environment does not have the intended XGBoost build.
version_str <- system2(
  py, args = c("-c", shQuote("import xgboost; print(xgboost.__version__)")),
  stdout = TRUE, stderr = TRUE
)
exit_code <- attr(version_str, "status")
if (!is.null(exit_code) && exit_code != 0L) {
  stop("Python XGBoost verification failed:\n",
       paste(version_str, collapse = "\n"), call. = FALSE)
}
version_str <- trimws(tail(version_str, 1L))
cat("Python for manifest:", py, "\n")
cat("Python XGBoost:", version_str, "\n")
if (!identical(version_str, "3.2.0")) {
  stop("Expected the Stage 21 Python XGBoost version 3.2.0; got ",
       version_str, call. = FALSE)
}

cat("Generating new manifest.json ...\n")
# appFiles is intentional: prevents accidentally including source ACS PUFs,
# analysis caches, or other unreviewed files in deployment metadata.
rsconnect::writeManifest(
  appDir = root,
  appFiles = app_files,
  appMode = "shiny",
  python = py,
  forceGeneratePythonEnvironment = FALSE
)

manifest <- jsonlite::read_json(file.path(root, "manifest.json"))
r_pkgs <- names(manifest$packages)
required_app_r <- c("shiny", "bslib", "data.table", "Matrix", "xgboost", "reticulate")
missing_r <- setdiff(required_app_r, r_pkgs)
if (length(missing_r)) {
  stop("Manifest missing R runtime packages: ",
       paste(missing_r, collapse = ", "),
       "; do NOT push this manifest.", call. = FALSE)
}
if (is.null(manifest$python) ||
    is.null(manifest$python$version) ||
    is.null(manifest$python$package_manager$package_file)) {
  stop(
    "Manifest lacks Python environment metadata or requirements path. ",
    "Do NOT push it; inspect rsconnect::writeManifest() logs.",
    call. = FALSE
  )
}
if (!identical(
  as.character(manifest$python$package_manager$package_file),
  "requirements.txt"
)) {
  stop("Manifest references unexpected Python requirements file; do NOT push.",
       call. = FALSE)
}
req <- readLines(file.path(root, "requirements.txt"), warn = FALSE)
if (!any(grepl("^\\s*xgboost==3\\.2\\.0\\s*$", req))) {
  stop("requirements.txt is not pinned to xgboost==3.2.0.",
       call. = FALSE)
}
cat("\nPASS: mixed-language manifest generated.\n",
    "R: ", manifest$platform, "\n",
    "Python: ", manifest$python$version, "\n",
    "R packages include: ", paste(required_app_r, collapse = ", "), "\n",
    "Python requirements: requirements.txt (xgboost==3.2.0)\n",
    "Push ONLY manifest.json to the integration branch, then deploy from GitHub.\n",
    sep = "")
