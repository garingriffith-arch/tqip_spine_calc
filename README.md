# SPINE-TRACT

**Trauma Resource and Acute Care Trajectory Calculator for Adults With Acute Thoracolumbar Spine Trauma**

SPINE-TRACT is a multi-output prediction-model research project evaluating early inpatient trajectories after acute thoracolumbar spine trauma. The [analysis pipeline](analytic_pipeline/README.md) is organized to parallel the TBI-TRACT publication repository.

The revised analysis uses ACS Trauma Quality Programs Participant Use Files (2020–2024) and a narrow cohort of 424,017 direct-presenting adults aged 18–89 years with an observable index-hospital trajectory. Evaluation used temporally ordered 2022–2024 cohorts (2024 N=90,649 before endpoint-specific missingness). This is **development-informed internal forward temporal evaluation**, not external validation.

### Selected 2024 discrimination

| Output | AUROC |
|---|---:|
| Post-acute disposition | 0.776 |
| Death/hospice | 0.925 |
| Any ICU | 0.874 |
| ICU ≥8 days | 0.898 |
| Any ventilation | 0.936 |
| Ventilation ≥8 days | 0.932 |
| Hospital LOS ≥20 days | 0.854 |
| Recorded thoracolumbar operation | 0.761 |

The median-prediction mean absolute error in 2024 was 4.82 days for hospital LOS, 4.21 days for ICU LOS conditional on ICU use, and 5.46 days for ventilator duration conditional on ventilation. Estimated Q10–Q90 ranges are not guaranteed 80% prediction intervals.

### Repository organization

- `analytic_pipeline/01_cohort/`: cohort construction, procedure phenotype, diagnosis QC
- `analytic_pipeline/02_model_development/`: temporal predictor adjudication, classification, duration quantiles, ridge comparator
- `analytic_pipeline/03_sensitivity/`: clinical, airway, phenotype and composite sensitivities
- `analytic_pipeline/05_reporting/`: fixed-prediction performance, calibration, main/supplemental tables
- `analytic_pipeline/implementation/`: preserved versioned source code
- `analytic_pipeline/04_deployment/`: final deployment artifacts and policy documentation (all-years fit **pending**)
- `app.R` and `data/model_bundle.rds`: **legacy** calculator—do not conflate with revised manuscript models

Details and run order: [`analytic_pipeline/README.md`](analytic_pipeline/README.md). File-level provenance and checksums: [`analytic_pipeline/SOURCE_ARCHIVE_MANIFEST.csv`](analytic_pipeline/SOURCE_ARCHIVE_MANIFEST.csv).

### Interpretation and availability

SPINE-TRACT predicts **recorded operative treatment receipt**, not operative indication. The PUF does not provide CT/MRI fracture morphology, posterior ligamentous complex integrity, TLICS/AO Spine classifications, or ASIA/AIS neurologic examination. It is not designed to assess spinal stability, choose surgery, or predict neurological recovery. Independent external and prospective validation are pending.

Patient-level ACS TQP PUF data, derived patient-level records, and prediction caches are not included. Access to those source data is separately governed by ACS terms. All-years deployment models and revised calculator integration are forthcoming and will be distinguished explicitly from the legacy Shiny app.
