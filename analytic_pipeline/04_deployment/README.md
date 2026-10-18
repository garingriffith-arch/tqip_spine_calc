# Deployment is distinct from manuscript temporal evaluation

The root `app.R` and `data/model_bundle.rds` represent a **previous SPINE-TRACT model generation**. They are retained to avoid disrupting an existing application, but they do not implement the new 54/57-predictor hybrid trajectory and Q10/Q50/Q90 framework.

The manuscript-facing 2022–2024 rolling evaluations are **not** deployment models. A later targeted all-2020–2024 fitting step must produce the definitive model objects and encoders, with a tested inference function, before replacing or releasing the new calculator.

No final deployment fitting script or all-years model objects are being asserted to exist in this GitHub archive. Independent external/prospective validation remains pending. The tool is not intended to determine operative indications, stability, neurological recovery, or whether ICU/ventilation is required.
