# SPINE-TRACT reporting entry point: manuscript tables and supplement
#
# Run this from the repository root after setting SPINE_TRACT_PROJECT_DIR to
# the local ACS analysis workspace that contains the locked repaired-DX inputs.
# This script builds Word tables and source CSVs only. It does not fit models.
#
# Example:
#   Sys.setenv(SPINE_TRACT_PROJECT_DIR = "D:/Research/TQIP Spine Utilization")
#   source("analytic_pipeline/05_reporting/30_build_tables.R")

implementation <- file.path(
  "analytic_pipeline", "implementation",
  "30_build_SPINE_TRACT_final_tables_PUBLIC_REVIEWER_v4_1.R"
)
if (!file.exists(implementation)) {
  stop("Missing implementation: ", implementation, call. = FALSE)
}
source(implementation, local = FALSE, echo = FALSE)
