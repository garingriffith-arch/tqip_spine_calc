# Stage 32: 2024 fixed-prediction subgroup uncertainty for eTable 11.
# Requires licensed repaired-DX cohort and evaluation prediction caches.
# Does not refit any prediction model. Set SPINE_TRACT_ETABLE11_B=1000.
implementation <- file.path("analytic_pipeline", "implementation",
                            "32_build_SPINE_TRACT_eTable11.R")
if (!file.exists(implementation)) {
  stop("Missing implementation: ", implementation, call. = FALSE)
}
source(implementation, local = FALSE, echo = FALSE)
