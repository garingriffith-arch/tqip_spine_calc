# Stage 21: locked 2020–2024 all-years deployment fitting.
# Requires licensed repaired-DX cohort and local Stage 3 tuning results.
# Internal CV selects boosting rounds only; it is not independent validation.
implementation <- file.path("analytic_pipeline", "implementation",
                            "21_build_SPINE_TRACT_deployment_2020_2024_CVFIX_v2.R")
if (!file.exists(implementation)) {
  stop("Missing implementation: ", implementation, call. = FALSE)
}
source(implementation, local = FALSE, echo = FALSE)
