# Stage 31: manuscript figure generation from the frozen local analysis outputs.
# Latest supplied correction; execution-to-approved-output parity is not established.
project_dir <- Sys.getenv("SPINE_TRACT_PROJECT_DIR", unset = "")
if (!nzchar(project_dir) || !dir.exists(project_dir)) {
  stop("Set SPINE_TRACT_PROJECT_DIR to the authorized local SPINE-TRACT workspace.", call. = FALSE)
}
implementation <- file.path("analytic_pipeline", "implementation", "31_build_SPINE_TRACT_figures_v4_SYMMETRIC_FIG2.R")
if (!file.exists(implementation)) stop("Missing implementation: ", implementation, call. = FALSE)
source(implementation, local = FALSE, echo = FALSE)
