# SPINE-TRACT analytic pipeline

This directory archives the executable source for the final development-informed rolling-origin 2020–2024 SPINE-TRACT manuscript analyses. The directory organization follows the TBI-TRACT pipeline: stable numbered entrypoints plus preserved, versioned implementations.

## Study scope

ACS Trauma Quality Programs Participant Use Files, 2020–2024; adults 18–89 years presenting directly with an acute qualifying thoracolumbar trauma phenotype and an observable index-hospital trajectory (N=424,017). Evaluation years are 2022, 2023 and 2024. The 2024 evaluation cohort contained 90,649 patients prior to outcome-specific exclusions.

Classification targets: three-class disposition, hospital LOS, ICU and ventilation trajectories, and binary **recorded therapeutic thoracolumbar operation** (treatment receipt, not clinical indication). Duration targets: Q10/Q50/Q90 hospital LOS and conditional ICU/ventilator duration. Core and contextual policies comprise 54 and 57 predictors, respectively.

Evaluation is **internal, development-informed, forward temporal** rather than independent external validation. No external/prospective validation has been performed.

## Running

Review [`PIPELINE_MANIFEST.csv`](PIPELINE_MANIFEST.csv) and [`SOURCE_ARCHIVE_MANIFEST.csv`](SOURCE_ARCHIVE_MANIFEST.csv) before executing. From the repository root, `Rscript analytic_pipeline/verify_archive.R` checks the 20 numbered entrypoints and implementation files without accessing ACS data or fitting models. Set `SPINE_TRACT_PROJECT_DIR` to the authorized local ACS analysis workspace. Numbered entrypoints are in `01_cohort/`, `02_model_development/`, `03_sensitivity/`, and `05_reporting/`. The `implementation/` scripts are preserved without cosmetic rewriting because they are the analysis provenance sources. Several stages require frozen prior outputs, the repaired diagnosis cohort, a locally built GPU-enabled XGBoost version, and local path configuration. The archive is **not** a no-data, one-command reproduction.

Workflow: cohort/endpoint audit → Stage 2 candidate construction and required diagnosis repair → repaired predictor screening/finalist lock → temporal classification/duration fits → ridge comparator and robustness → reporting, numerical QC, and publication tables. Numeric duration ties and pooled AUROC have targeted postprocessing corrections; **do not rerun expensive fitting solely to apply those corrections**.

## Availability limitations

The source package includes all executable final-analysis stages recovered and reconciled from the locked code archives, with verified Git blob SHA-1 identifiers in SOURCE_ARCHIVE_MANIFEST.csv. However, the **latest locally revised Figure 1–3/eFigure 1–6 code** and the final verified **B=1000 TRIPOD+AI subgroup reporting script** are being finalized and are not included as purportedly final in this snapshot. The independent final 2020–2024 deployment fit is pending.

No ACS patient-level PUF files, derived patient-level datasets, individual prediction caches or record-level identifiers are distributed. A curated therapeutic ICD-10-PCS code map may be shared as non-patient-level metadata under the existing repository `reference/` directory. No usable facility keys are available, so centre-level heterogeneity was not estimated.

Operation prediction is treatment receipt, not operative indication. The PUF does not include CT/MRI morphology, posterior ligamentous complex integrity, TLICS/AO classifications, or ASIA/AIS neurological examinations. The tool must not be used for stability assessment, operative selection or neurological prognosis.
