# Run from the app folder after installing R dependencies and Python XGBoost >=2.
# Checks structural prediction contracts; it does not validate clinical accuracy.
options(warn = 1)
env <- new.env(parent = globalenv())
source("app_core.R", local = env)
stopifnot(length(env$model_ids) == 8L, length(env$used) == 57L)
vals <- setNames(vector("list", length(env$used)), env$used)
for (v in env$used) {
  vals[[v]] <- if (v %in% env$categorical_predictors) "__UNKNOWN__" else NA_real_
}
for (v in c("age_years", "sbp_clean", "pulse_clean", "rr_clean", "spo2_clean",
            "gcs_eye_clean", "gcs_verbal_clean", "gcs_motor_clean")) {
  if (v %in% names(vals)) vals[[v]] <- switch(v,
    age_years = 50, sbp_clean = 125, pulse_clean = 85, rr_clean = 18,
    spo2_clean = 97, gcs_eye_clean = 4, gcs_verbal_clean = 5, gcs_motor_clean = 6)
}
for (id in c("discharge3_final", "hlos3_final", "icu3_final", "vent3_final")) {
  d <- env$pred_multi(id, vals)
  stopifnot(nrow(d) == 3L, all(is.finite(d$probability)),
            all(d$probability >= 0), abs(sum(d$probability) - 1) < 1e-6)
  cat("PASS multiclass:", id, "\n")
}
op <- env$pred_binary("operative_spine_intervention", vals)
stopifnot(length(op) == 1L, is.finite(op), op >= 0, op <= 1)
cat("PASS operation probability\n")
for (id in env$duration_ids) {
  q <- env$pred_quantiles(id, vals)
  stopifnot(identical(names(q), c("q10", "q50", "q90")),
            all(is.finite(q)), all(diff(q) >= 0),
            all(q >= as.numeric(env$metadata[[id]]$lower_bound)))
  cat("PASS joint quantiles:", id, "\n")
}
# A fixed one-row feature-order contract for every endpoint.
for (id in env$model_ids) {
  m <- env$metadata[[id]]
  x <- env$encode_values(m, vals, sparse = !(id %in% env$duration_ids))
  stopifnot(ncol(x) == length(m$feature_names),
            identical(colnames(x), as.character(m$feature_names)))
}
cat("\nPASS: all eight model interfaces and feature-order checks.\n")
