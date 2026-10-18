# SPINE-TRACT analytic pipeline

This directory is the **manuscript/reproducibility archive**, separate from the previously published Shiny app in the repository root. The organization mirrors the TBI-TRACT research repository while preserving SPINE-TRACT-specific cohort definitions, endpoints, and limitations.

## Locked study design

- **Source:** ACS TQP PUF, admission years 2020–2024.
- **Intended population:** direct-presenting adults aged 18–89 years with qualifying acute thoracolumbar vertebral fracture, thoracolumbar spinal cord injury, or root/cauda injury, and an observable index-hospital trajectory.
- **Primary cohort:** 424,017 patients; evaluation year 2024 contained 90,649 patients.
- **Prediction time:** after initial trauma-center examination, physiologic stabilization, and diagnostic workup, before later resource receipt/disposition or operative treatment.
- **Temporal evaluation:** train on 2020, tune on 2021, refit on 2020–2021 and evaluate 2022; then move each step forward for 2023 and 2024.
- **Model family:** XGBoost; ridge comparator for probability/trajectory outcomes.
- **Predictor policy:** final clinical core has 54 predictors; disposition and categorical hospital-LOS models add race, ethnicity, and payer (57). The other classification and duration models use the clinical core.
- **Probability targets:** three-class disposition, hospital length-of-stay, ICU, and ventilation trajectories, plus recorded thoracolumbar operative treatment.
- **Continuous targets:** Q10/Q50/Q90 hospital LOS and conditional ICU/ventilator duration.
- **Validation:** internal forward temporal evaluation only; independent external and prospective evaluation are still needed.

The surgical endpoint measures **recorded treatment receipt**, not operative indication, spinal instability, suitability for surgery, or neurological outcome. The registry does not capture CT/MRI fracture morphology, posterior ligamentous complex integrity, TLICS/AO classification, or ASIA/AIS examination. SPINE-TRACT must not be used to infer these.

## Folder layout

| Folder | Scope | Archive status |
|---|---|---|
| `01_cohort/` | Diagnosis/procedure phenotype, transfer exclusions, and dataset reconciliation | **647-code therapeutic operation map archived**; executable cohort source pending |
| `02_model_development/` | Predictor adjudication, 2022–2024 rolling folds, ridge comparator, duration quantiles | Source transfer pending |
| `03_sensitivity/` | Airway, broad phenotype, subgroups, early death, composite outcomes | Source transfer pending |
| `04_deployment/` | Final all-years deployment model fit, encoders, and deployment checks | Not finalized; legacy app is separate |
| `05_reporting/` | Tables, figures, TRIPOD+AI reporting, and model-independent QC | **Table builder and temporal-metric source/pooled-AUROC correction archived**; final figure source pending |
| `implementation/` | Historical versioned executable sources | See individual file names |

Read `PIPELINE_MANIFEST.md` before treating the public archive as a complete end-to-end pipeline.

## Available executable source

`05_reporting/30_build_tables.R` sources the final reviewer-facing Word table pipeline from `implementation/30_build_SPINE_TRACT_final_tables_PUBLIC_REVIEWER_v4_1.R`.

The implementation accepts `SPINE_TRACT_PROJECT_DIR` and discovers precomputed **local** analytic data and reporting artifacts beneath that directory. It outputs `Tables.docx`, `Supplemental.docx`, and machine-readable table sources to `outputs/SPINE_TRACT_PUBLICATION_TABLES/`. It performs no prediction-model fitting or recalibration.

**Requirements:** R; `data.table`, `officer`, `flextable`, and `scales`; the frozen repaired diagnosis cohort and predictor specification, Gate 1 flow/operation summaries, Stage 5 classification/duration metrics/CIs, and optional Stage 2–4 robustness/comparator artifacts. Not all of these are publicly distributable.

## Provenance and next steps

A vetted non-patient-level therapeutic operation map is available at [`reference/spine_operative_icd10pcs_map.csv`](reference/spine_operative_icd10pcs_map.csv), with 647 unique ICD-10-PCS codes across six categories. The Stage 5 fixed-prediction reporting implementation and the numerical pooled AUROC correction are archived in `implementation/`; do not use the original pooled AUROC output without that correction.\n\nThe diagnostic phenotype was corrected before final model development; final analyses consistently use the repaired 54/57-predictor architecture. Old random-split SPINE-TRACT and early application models must not be mixed with these estimates.

The repository should eventually include checked versions of the canonical cohort, predictor adjudication, final classification/duration, ridge, robustness, reporting and figure-generation scripts; software/seed manifests; aggregate QC outputs; and, after the deployment work is frozen, the all-years model objects and encoder metadata.

Do not commit patient-level data or intermediate identifiable/record-level files. There are **no usable facility keys**, so facility-specific performance heterogeneity is not estimated.

Last documentation update: October 2026.
