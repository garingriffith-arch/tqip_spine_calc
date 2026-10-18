# SPINE-TRACT

**Trauma Resource and Acute Care Trajectory Calculator for Adults With Acute Thoracolumbar Spine Trauma**

SPINE-TRACT is a multi-endpoint prediction-model research framework that characterizes observed inpatient trajectories after acute thoracolumbar spine trauma. This repository contains a **legacy Shiny application** and an evolving, separate **manuscript-facing reproducibility archive**. The former should not be confused with the models evaluated in the current revised study.

## Revised manuscript study at a glance

| Component | Final manuscript specification |
|---|---|
| Data source | American College of Surgeons (ACS) Trauma Quality Programs Participant Use Files, 2020–2024 |
| Intended cohort | 424,017 directly presenting adults aged 18–89 years with a narrow qualifying thoracolumbar trauma phenotype and observable index-hospital trajectory |
| 2024 evaluation cohort | 90,649 patients before outcome-specific missingness restrictions |
| Prediction time | After initial trauma-center assessment, stabilization, and diagnostic workup, but before subsequent inpatient resource utilization or operation |
| Evaluation | Rolling-origin internal temporal evaluation in 2022, 2023, and 2024 |
| Predictors | 54 clinical predictors; disposition/categorical hospital-LOS models additionally use race, ethnicity, and payer (57 total) |
| Probabilities | 3-class disposition, hospital LOS, ICU, and ventilation trajectories; recorded thoracolumbar operative treatment |
| Duration estimates | Q10/Q50/Q90 hospital LOS; ICU LOS conditional on ICU use; ventilation duration conditional on ventilation |
| Comparator | Ridge-based classification models trained under matched predictor and evaluation policies |
| Validation | Internal temporal evaluation; independent external and prospective validation pending |

**Selected 2024 discrimination:** AUROC 0.776 for post-acute facility disposition, 0.925 for death/hospice, 0.874 for any ICU, 0.898 for ICU ≥8 days, 0.936 for any ventilation, 0.932 for ventilation ≥8 days, 0.854 for hospital LOS ≥20 days, and 0.761 for recorded thoracolumbar operative treatment.

**2024 continuous-duration accuracy:** median-prediction mean absolute error of 4.82 days for hospital LOS, 4.21 days for ICU LOS among ICU users, and 5.46 days for ventilation duration among ventilated patients. Displayed Q10–Q90 values are *estimated 10th–90th percentile ranges*, not guaranteed coverage intervals.

**Important interpretation:** The 2024 cohort was used in development-informed forward temporal evaluation, not as an untouched independent external holdout. This model estimates **treatment receipt**, not the clinical indication for spine surgery. It does not incorporate CT/MRI fracture morphology, posterior ligamentous complex integrity, TLICS/AO Spine classification, or ASIA/AIS examination and is not intended to determine spinal stability or neurological prognosis.

## Repository layout

- `app.R` — **legacy** Shiny app from an earlier SPINE-TRACT model iteration, currently untouched.
- `data/model_bundle.rds` and `data/predictor_metadata.rds` — associated **legacy** application artifacts; **not final revised-manuscript models**.
- `analytic_pipeline/` — revised-study research code, stages, and documentation, structured in parallel with [TBI-TRACT](https://github.com/garingriffith-arch/tqip_tbi_calc).
- `analytic_pipeline/implementation/` — preserved executable sources; currently includes the reviewer-facing table builder v4.1.
- `analytic_pipeline/05_reporting/` — numbered reporting entry point and documentation.
- `R/00_config.example.R` — example of configuring a local, authorized data workspace.
- `CODE_AVAILABILITY.md` — current public-code scope and archival limitations.

For exact archive status, see the [pipeline README](analytic_pipeline/README.md) and [manifest](analytic_pipeline/PIPELINE_MANIFEST.csv). The repository is **not yet a complete clean-clone reproduction** of the revised study; additional finalized stage sources will be added only after provenance/QA checks.

## Running the archived manuscript table builder

From the repository root, after obtaining authorized access to the source data and constructing the locked local analysis outputs:

```r
install.packages(c("data.table", "officer", "flextable", "scales"))
Sys.setenv(SPINE_TRACT_PROJECT_DIR = "PATH/TO/YOUR/LOCAL/TQIP Spine Utilization")
source("analytic_pipeline/05_reporting/30_build_tables.R")
```

It writes `Tables.docx`, `Supplemental.docx`, and associated source CSVs in the local analysis workspace's `outputs/SPINE_TRACT_PUBLICATION_TABLES/` directory. It does **not** train, tune, or recalibrate clinical models.

## Legacy app

The repository's existing `app.R` accepts an older binary resource-model bundle. It has **not** been migrated to the final three-class/quantile trajectory architecture. Do not use its output as the current manuscript model or treat the displayed estimated risks as independently validated predictions.

Following final model freeze and independent validation, deployment code and metadata can be updated separately without changing the already reported 2022–2024 temporal estimates.

## Data availability and model availability

ACS patient-level TQP PUF files, derived patient-level datasets, and prediction caches are **not distributed**. Reproduction requires authorized access and compliance with ACS data-use terms.

The final revised deployment model objects have not been archived in this repository. No facility-level heterogeneity analysis was possible because usable facility keys were unavailable.

This project is intended for research, counseling, and anticipatory resource characterization—not autonomous treatment allocation, operative selection, or neurological prognostication.

## Related project

[TBI-TRACT: traumatic brain injury trajectory calculator and analytic pipeline](https://github.com/garingriffith-arch/tqip_tbi_calc).

*Documentation updated October 2026.*
