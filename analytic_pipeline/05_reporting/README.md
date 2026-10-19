# SPINE-TRACT manuscript reporting

## Available archived code

- `16_build_prediction_reporting.R` generates fixed-prediction temporal performance, calibration and bootstrap outputs.
- `17_build_final_manuscript_outputs.R` assembles the manuscript-facing model summaries.
- `18_repair_pooled_auroc.R` and `19_numeric_auroc_qc.R` implement targeted pooled-AUROC numerical corrections without refitting models.
- `30_build_tables.R` builds three main tables and eTables 1–10 from the repaired diagnosis cohort and locked upstream aggregates.
- `32_build_TRIPOD_AI_support.R` builds eTable 11 fixed-prediction 2024 subgroup uncertainty: AUROC/AUPRC/Brier 1,000 patient-level bootstrap percentile intervals and logistic-calibration Wald intervals. It excludes subgroups with fewer than 500 eligible records or fewer than 25 events/non-events; it does not refit any classifier.

The eTable 11 code was executed and reviewed locally, but its licensed record-level prediction inputs and checkpoint RDS files are not distributed in the public repository.

## Unavailable source

**The final revised code for Figure 1–3 and eFigures 1–6 has not been archived.** Do not claim those figures can be regenerated from this repository until the final executed plotting source is supplied and verified. See [pipeline manifest](../PIPELINE_MANIFEST.csv).

## Execution

Set `SPINE_TRACT_PROJECT_DIR` to the licensed workspace. Run numbered entrypoints from the repository root after upstream repaired-DX datasets and model prediction caches exist. Do not substitute all-years deployment model predictions for the historical 2024 fold-specific evaluation predictions. Reporting is not a new model-fitting or validation stage.
