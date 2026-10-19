# SPINE-TRACT analysis pipeline

This directory contains the numbered analysis workflow and source archive for the repaired-diagnosis 2020–2024 SPINE-TRACT study. **A clean clone verifies the source files and deployment models, but cannot reproduce the patient-level analyses without licensed ACS TQP PUF records and the specified locally generated intermediates.**

## Methods and data

- Primary cohort: 424,017 direct-presenting adults aged 18–89 years (2020–2024) after the repaired diagnosis-phenotype cohort construction.
- Evaluation: rolling-origin internal temporal evaluation in 2022, 2023, and 2024, with model-architecture decisions informed by evaluation. **Not independent external validation.**
- Predictors: 54 clinical variables for the ICU, ventilation, operation and duration models; 57 contextual variables for disposition and categorical hospital LOS.
- Model families: multiclass XGBoost for four trajectory outcomes; binary XGBoost for recorded therapeutic thoracolumbar operation; joint Q10/Q50/Q90 XGBoost for three continuous-duration endpoints.
- Missing data: native numeric missing-value routing; categorical `__UNKNOWN__` and `__OTHER__` encoding, using training-derived levels.
- Recorded operation predicts treatment **receipt**, not necessity or appropriateness.

## Requirements

From the repository root, first check source integrity:

```sh
Rscript analytic_pipeline/verify_archive.R
Rscript analytic_pipeline/verify_deployment_artifacts.R
```

These tests require Git; deployment SHA-256 checks also require the R package `digest`. They do not read the licensed PUF.

For analysis reruns, set `SPINE_TRACT_PROJECT_DIR` to the authorized local research workspace and prepare the upstream PUF fields. Most historical sources use the repaired frozen cohort at:

```text
<SPINE_TRACT_PROJECT_DIR>/data/derived/stage2_predictor_adjudication/
  spine_tract_frozen_primary_2020_2024.rds
  spine_tract_predictor_sets.rds
```

Some stages additionally require the source ACS annual PUFs, phenotype mappings, prior-stage prediction caches, Python/CUDA XGBoost, or saved grid-search outputs. Their preflight checks and the [pipeline manifest](PIPELINE_MANIFEST.csv) describe these dependencies. This archive does not promise a clean-machine, single-command reanalysis.

## Numbered workflow

| Steps | Location | Purpose |
| --- | --- | --- |
| 01–03 | [01_cohort](01_cohort/) | Source audit; candidate construction; **mandatory repaired diagnosis-phenotype step** |
| 04–08 | [02_model_development](02_model_development/) | Predictor adjudication, locked 54/57 architecture, rolling-origin classification and joint duration fits |
| 09–11 | [02_model_development](02_model_development/) | Duration calibration postprocessing and ridge-comparator sensitivity |
| 12–15 | [03_sensitivity](03_sensitivity/) | Fixed-prediction, composite, phenotype and operation-endpoint robustness |
| 16–20 | [05_reporting](05_reporting/) | Discrimination, calibration, numerical QC and main/eTables 1–10 |
| 21 | [04_deployment](04_deployment/) | Final full-cohort deployment fits; internal CV selects boosting rounds only |
| 32 | [05_reporting](05_reporting/) | eTable 11 fixed-prediction 2024 subgroup bootstrap (B=1,000) |

**Unarchived:** The final code used to generate Figure 1–3 and eFigures 1–6 (pipeline step 31). It is listed as pending rather than reconstructed without evidence. The figures can be examined in the existing manuscript package, but the last plotting revision cannot presently be reproduced from this GitHub source alone. **Stage 21 source fidelity:** the repository archives a corrected CV-parser deployment builder, but the original executable recorded an MD5 of `64cc76848eeeef02fec044f7345ee464` in the artifact manifest. The exact script with that checksum has not been supplied, so byte-for-byte agreement with the executed fit has not been verified.

The deployment models are final all-years fitted artifacts and must not be substituted for the 2022–2024 fold-specific models when recreating temporal performance. The source [PIPELINE_MANIFEST.csv](PIPELINE_MANIFEST.csv) maps entrypoints to versions; [SOURCE_ARCHIVE_MANIFEST.csv](SOURCE_ARCHIVE_MANIFEST.csv) records Git blob SHA-1 identifiers. Neither manuscript performance nor external validity can be inferred from deployment cross-validation used for choosing boosting rounds.

## Restricted inputs and outputs

This repository omits raw ACS PUF records, derived patient-level cohorts, individual predictions, record identifiers, and analysis caches. It contains an aggregate deployment artifact manifest and native trained model files; no patient-level records are required for the Shiny app. Review ACS/institutional terms before redistribution. There are no usable facility identifiers in the manuscript cohort, so centre-level heterogeneity was not evaluated.

## Known corrections preserved

The archived workflow explicitly includes the diagnosis-derived phenotype repair, the duration calibration and pooled-AUROC numerical corrections, and the corrected Python XGBoost joint-quantile CV result parser. These are historical analysis steps, not optional cosmetic fixes. Repeating model fitting is **not** necessary merely to apply postprocessing corrections.
