# Deployment-model inventory and status

SPINE-TRACT's revised manuscript analyses are **not** the same as the repository-root legacy Shiny calculator. The existing `app.R`, `data/model_bundle.rds`, and `data/predictor_metadata.rds` belong to an earlier generation and must not be relabeled or redeployed as the final revised model.

## Final revised deployment models — not yet fitted

The locked 2020–2024 architecture calls for **14 fitted XGBoost objects**:

| Model | Count | Predictor policy |
|---|---:|---|
| Three-class discharge disposition (home/home health; post-acute facility; death/hospice) | 1 | 57 contextual |
| Three-class hospital LOS (≤7; 8–19; ≥20 days) | 1 | 57 contextual |
| Three-class ICU trajectory (none; 1–7; ≥8 days) | 1 | 54 clinical |
| Three-class mechanical-ventilation trajectory (none; 1–7; ≥8 days) | 1 | 54 clinical |
| Recorded therapeutic thoracolumbar operation (binary treatment receipt) | 1 | 54 clinical |
| Hospital LOS Q10, Q50, Q90 | 3 | 54 clinical |
| ICU LOS Q10, Q50, Q90 conditional on ICU use | 3 | 54 clinical |
| Ventilator days Q10, Q50, Q90 conditional on ventilation | 3 | 54 clinical |
| **Total** | **14** | |

The 57-predictor policy is the final 54-predictor clinical core **plus race, ethnicity, and payer** for the disposition and categorical hospital-LOS trajectories only. The contextual variables are not used by ICU, ventilation, operative-treatment, or continuous-duration models.

Derived probabilities such as “any ICU use” and “ICU ≥8 days” come from the ICU multiclass model. They do **not** require separately fitted binary deployment models.

## Existing evaluated model files: not the all-years deployment fit

The repaired-diagnosis rolling-origin classification models were saved locally under:

`outputs/stage3_final_classification_REPAIRED_DX/rolling_origin_gpu/models/`

The repaired-diagnosis duration-quantile models were saved locally under:

`outputs/stage3_duration_quantiles_REPAIRED_DX/rolling_origin_gpu/models/`

These are year-specific forward-evaluation fits for 2022, 2023, and 2024, **not** final models fitted using all five years. Their encoders, selected configuration and evaluation year must be preserved for reproducibility, but they must not be silently moved into the live calculator.

For the authoritative fold configurations see `14_FINAL_SELECTED_CONFIGS_BY_FOLD.csv` in the classification output folder and `04_SELECTED_DURATION_CONFIGS.csv` in the duration output folder. The endpoint-matched ridge comparator is for manuscript benchmarking, **not** deployment.

## Artifact destination after the all-years refit

Recommended namespaced location:

- `data/deployment_2020_2024/models/` — five classifier model objects plus nine quantile models (or five endpoint bundles preserving all 14 fitted models)
- `data/deployment_2020_2024/encoders/` — endpoint-specific preprocessing and levels/column-order metadata
- `data/deployment_2020_2024/manifest.csv` — script provenance, training N, predictor-set policy, selected structural settings, all-years round-selection method, versions, and SHA-256 checksums

Retain the trained model objects **with** their complete encoders and inference/preprocessing function; a bare XGBoost booster is insufficient to reproduce patient-level inference. Do not use the existing `data/models/` destination until intentionally migrating and testing the legacy application.

The final refit must reuse the locked predictor/endpoint policy and prevent leakage; changes to the deployment-fitting procedure must not retroactively alter the manuscript's 2022–2024 forward-temporal performance estimates.

The all-years final refit and accompanying inference validation are **pending**. External and prospective validation are required before routine clinical use. SPINE-TRACT does not use CT/MRI fracture morphology, posterior ligamentous complex status, TLICS/AO classification, or ASIA/AIS neurological examination. Recorded operation predicts treatment receipt, not operative indication.
