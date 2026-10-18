# Deployment model status

A final, all-2020–2024 deployment fit using the locked revised SPINE-TRACT predictor policies **has not yet been run or archived**. The `models/` directories created by the rolling-origin 2022–2024 scripts contain **fold-specific evaluation models** and must not be copied to the Shiny app as final deployment models.

The existing repository-root `app.R`, `data/model_bundle.rds`, and `data/predictor_metadata.rds` belong to a superseded model generation. Preserve them as legacy, not the revised manuscript architecture.

## Final deployment artifacts to produce

The deployment refit should generate five classifier model objects (three-class disposition, three-class hospital LOS, three-class ICU, three-class ventilation, and binary recorded thoracolumbar operation) and nine quantile regressors (Q10/Q50/Q90 separately for hospital LOS, ICU LOS conditional on ICU use, and ventilator days conditional on ventilation). **Total: 14 fitted XGBoost models**, normally saved as one RDS per classifier and one RDS per quantile, or as five endpoint bundles preserving all underlying model objects.

For classification, use the final **57-input contextual predictor policy** for disposition and categorical hospital LOS, and the **54-input clinical core** for ICU, ventilation, and recorded thoracolumbar operation. All nine duration models use the 54-input clinical core. For each model preserve the encoder, feature order, categorical levels, numeric missing routing, selected structural grid configuration, chosen all-years boosting rounds, outcome labels, training population and exact script/seed/software versions.

Recommended final artifact location after refit and inference QA: `data/deployment_2020_2024/models/` and `data/deployment_2020_2024/encoders/`, with an immutable manifest and SHA-256 checksums. **Do not** use `data/models/` for revised objects until the legacy app has been intentionally migrated and verified.

Independent external/prospective validation is still required before routine clinical use. Do not describe internal temporal evaluation as external validation.
