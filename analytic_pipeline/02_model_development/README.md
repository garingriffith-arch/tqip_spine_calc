# Temporal model development

The final SPINE-TRACT probability framework uses the following **rolling-origin internal temporal evaluation**:

| Evaluation | Initial training | Hyperparameter tuning | Refit | Evaluation cohort |
|---|---|---|---|---|
| 2022 | 2020 | 2021 | 2020–2021 | 2022 |
| 2023 | 2020–2021 | 2022 | 2020–2022 | 2023 |
| 2024 | 2020–2022 | 2023 | 2020–2023 | 2024 |

Final architecture:

- 54 predictor clinical core; 57 contextual inputs for disposition and categorical hospital LOS only;
- three-class disposition, hospital LOS, ICU and ventilator trajectories;
- binary recorded thoracolumbar operative treatment;
- Q10/Q50/Q90 continuous durations: hospital LOS in patients with observed LOS, ICU LOS conditional on ICU use, and ventilation duration conditional on ventilation;
- ridge comparator under the same forward temporal policy.

**Archive status:** Numbered entrypoints 04–11 and their preserved implementations are archived. Reproduction requires the licensed PUFs, repaired frozen cohort and the documented local GPU/XGBoost environment. Published performance estimates derive from the frozen revised analysis, never the legacy Shiny bundle.

The evaluation years were also consulted during predictor architecture development; they are not independent external validation cohorts.
