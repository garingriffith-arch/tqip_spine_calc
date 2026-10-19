# SPINE-TRACT

**Trauma Resource and Acute Care Trajectory Calculator for Adults With Acute Thoracolumbar Spine Trauma**

**[Open the SPINE-TRACT calculator](https://01a154f9-e4de-f916-8a37-2663ed8ed9d1.share.connect.posit.cloud)**

SPINE-TRACT estimates inpatient trajectories and recorded resource use after the initial assessment of adults with acute thoracolumbar spine trauma. This repository contains the Shiny application, its trained deployment models and preprocessing encoders, and the archived R analysis pipeline.

## Cohort and evaluation

The study uses the American College of Surgeons Trauma Quality Programs Participant Use Files, 2020–2024. The primary cohort includes **424,017 direct-presenting adults aged 18–89 years** with a qualifying acute thoracolumbar injury and an observable index-hospital trajectory.

Models were developed and evaluated using rolling-origin temporal splits with evaluation years **2022, 2023, and 2024**. The 2024 evaluation cohort included **90,649 patients** before endpoint-specific exclusions. The analysis represents **development-informed internal temporal evaluation**, not independent external or prospective validation.

| 2024 classification output | AUROC |
| --- | ---: |
| Post-acute facility disposition | 0.776 |
| Death or hospice disposition | 0.925 |
| Any ICU use | 0.874 |
| ICU stay ≥8 days | 0.898 |
| Any mechanical ventilation | 0.936 |
| Mechanical ventilation ≥8 days | 0.932 |
| Hospital LOS ≥20 days | 0.854 |
| Recorded thoracolumbar operation | 0.761 |

For 2024, median-prediction mean absolute error was **4.82 days** for hospital length of stay (LOS), **4.21 days** for ICU LOS conditional on ICU use, and **5.46 days** for ventilator days conditional on ventilation. The model-estimated Q10–Q90 ranges are not guaranteed 80% prediction intervals.

## Models

The locked predictor policy consists of **54 clinical predictors**, with race, ethnicity, and payer added only for discharge disposition and categorical hospital LOS (**57 contextual predictors**). The deployment build fits **eight XGBoost models** on the eligible 2020–2024 cohort:

- Five classifiers: three-class discharge disposition; three-class hospital LOS, ICU, and ventilation trajectories; and binary recorded thoracolumbar operation.
- Three joint-quantile duration models: Q10/Q50/Q90 for hospital LOS, ICU LOS conditional on ICU use, and ventilator days conditional on ventilation.

The all-years deployment models are distinct from the earlier rolling-origin evaluation models. Three duration models use Python XGBoost 3.2.0 for joint-quantile inference; the classifiers use R XGBoost.

## Repository layout

| Path | Contents |
| --- | --- |
| `app.R`, `app_core.R` | Shiny entrypoint, model inference, clinical inputs, output presentation |
| `data/deployment_2020_2024/` | Eight model JSONs, eight encoder RDSs, SHA-256 artifact manifest |
| `test_app_contract.R` | Structural tests for all eight inference interfaces |
| `analytic_pipeline/` | Numbered source entrypoints, archived implementations, provenance and integrity checks |
| `README_DEPLOYMENT.md` | Local app testing and Posit Connect Cloud publishing |
| `www/ohsu_logo.png` | Application logo |

The [analysis pipeline documentation](analytic_pipeline/README.md) gives execution order, inputs, outcome definitions, model-fit dependencies and the archive's remaining limits.

From the repository root:

```sh
Rscript analytic_pipeline/verify_archive.R
Rscript analytic_pipeline/verify_deployment_artifacts.R
```

These checks verify archived source integrity and deployment artifact checksums; they do **not** run the analyses or establish external validation.

## Clinical interpretation

The operative endpoint predicts **recorded therapeutic thoracolumbar operation**, not surgical indication or treatment appropriateness. ACS registry data do not encode fracture morphology from CT/MRI, posterior ligamentous complex integrity, TLICS/AO classification, or ASIA/AIS neurological examination. SPINE-TRACT does not assess spinal stability, select surgery, or predict neurological recovery.

Outputs are research-derived estimates for counseling and anticipatory resource planning, not stand-alone treatment recommendations. Social-context predictors can reflect unequal care pathways and should not be interpreted as biological effects or used to ration care. External and prospective validation remain outstanding.

## Data and code availability

The licensed ACS PUF source files, patient-level derived records, encounter identifiers and individual prediction caches are **not redistributed**. The public deployment model objects and preprocessing metadata are supplied under `data/deployment_2020_2024/`; sharing is subject to applicable ACS and institutional terms.

The first 20 analysis sources, the completed all-years deployment source (step 21), and the fixed-prediction eTable 11 source (step 32) are archived. The final figure-generation script for Figure 1–3 and eFigures 1–6 is **not yet archived**. See [code availability](CODE_AVAILABILITY.md) for exact limits. A clean clone cannot reproduce the patient-level analyses without licensed source data and required intermediate files.
