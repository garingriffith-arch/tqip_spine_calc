# Stage 21: all-2020–2024 deployment fitting (historical executed-source candidate).
# This fit is distinct from rolling-origin temporal evaluation.
# Do not rerun unless explicitly rebuilding deployment models.
project_dir <- Sys.getenv("SPINE_TRACT_PROJECT_DIR", unset = "")
if (!nzchar(project_dir) || !dir.exists(project_dir)) {
  stop("Set SPINE_TRACT_PROJECT_DIR to the authorized local SPINE-TRACT workspace.", call. = FALSE)
}
implementation <- file.path("analytic_pipeline", "implementation", "21_build_SPINE_TRACT_deployment_2020_2024_FIXED_20261019.R")
if (!file.exists(implementation)) stop("Missing implementation: ", implementation, call. = FALSE)
source(implementation, local = FALSE, echo = FALSE)
