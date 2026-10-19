# SPINE-TRACT manuscript reporting

## Available archived code

- `16_build_prediction_reporting.R` generates fixed-prediction temporal performance, calibration and bootstrap outputs.
- `17_build_final_manuscript_outputs.R` assembles the manuscript-facing model summaries.
- `18_repair_pooled_auroc.R` and `19_numeric_auroc_qc.R` implement targeted pooled-AUROC numerical corrections without refitting models.
- `30_build_tables.R` builds three main tables and eTables 1–10 from the repaired diagnosis cohort and locked upstream aggregates.
- `32_build_TRIPOD_AI_support.R` builds eTable 11 fixed-prediction 2024 subgroup uncertainty: AUROC/AUPRC/Brier 1,000 patient-level bootstrap percentile intervals and logistic-calibration Wald intervals. It excludes subgroups with fewer than 500 eligible records or fewer than 25 events/non-events; it does not refit any classifier.

The eTable 11 code was executed and reviewed locally, but its licensed record-level prediction inputs and checkpoint RDS files are not distributed in the public repository.

## Figure source and limits

`31_build_figures.R` references `../implementation/31_build_SPINE_TRACT_figures_v4_SYMMETRIC_FIG2.R`, the latest supplied corrected plot builder. The original source bytes are preserved. It writes Figures 1–3 and eFigures 1–6 as 300-dpi PNG and 600-dpi LZW TIFF files, plus figure source-data CSVs and manifests. No prediction model is fit.

**Output parity remains unverified.** The available figure ZIPs predate the last supplied Figure 2 correction; the repository does not contain the fully approved final rendered figures as a verified comparison set. Do not claim numerical or pixel-level reproduction merely because the source is archived.

To run, set `SPINE_TRACT_PROJECT_DIR` to the authorized local workspace; historical source files preserve their original fallback paths for provenance, while the entrypoint explicitly requires this environment variable. Required upstream files include the repaired frozen primary RDS, predictor-set RDS, Gate 1 cohort-flow CSV, Stage 5 frozen temporal classification/duration and 2024 bootstrap CI CSVs, Stage 2B/2C predictor adjudication CSVs, and Stage 4 subgroup robustness CSV. The script audits its discovered paths in `FIGURE_SOURCE_ARTIFACT_MANIFEST.csv`. These patient-level/intermediate inputs are **not** in a clean GitHub clone.

## Execution

Set `SPINE_TRACT_PROJECT_DIR` to the licensed workspace. Run numbered entrypoints from the repository root after upstream repaired-DX datasets and model prediction caches exist. Do not substitute all-years deployment model predictions for the historical 2024 fold-specific evaluation predictions. Reporting is not a new model-fitting or validation stage.
