# SPINE-TRACT local environment template
#
# Copy to R/00_config.R for local notes, and customize paths.
# The archived table script uses the SPINE_TRACT_PROJECT_DIR environment
# variable directly; sourcing this example does not configure other steps.

Sys.setenv(
  SPINE_TRACT_PROJECT_DIR = "PATH/TO/YOUR/LOCAL/TQIP Spine Utilization"
)

# ACS TQP PUF raw files and derived patient-level artifacts are stored locally
# under the authorized study environment; they are not part of this repository.
# Reproduction requires the corresponding source PUFs for 2020–2024.
