# ==============================================================================
# 10r_build_SPINE_TRACT_manuscript_outputs_REPAIRED_DX_v1_1.R
# SPINE-TRACT manuscript outputs — repaired-DX locked prediction caches
#
# PURPOSE
#   Create manuscript-facing output data from ALREADY-FROZEN Stage 3A
#   classification predictions and Stage 3C duration predictions.
#
# NO MODEL IS REFIT.
#
# OUTPUTS
#   - pooled rolling-origin 2022-2024 point estimates;
#   - 2024 headline classification performance with B=1000 bootstrap CIs;
#   - 2024 duration performance with B=1000 bootstrap CIs;
#   - flexible logistic-spline calibration curves for 2024 classification;
#   - 20-bin observed-vs-predicted calibration data;
#   - temporal 2022/2023/2024 point-estimate tables;
#   - canonical manifest / QC.
#
# PRIMARY DESIGN
#   Direct presentations only.
#   Disposition = Home/home health / Post-acute facility / Death/hospice.
#
# IMPORTANT
#   Bootstrap resamples patients within the fixed evaluation-year prediction
#   cache. This quantifies uncertainty of performance on the fixed locked model;
#   it is NOT model-refit stability.
# ==============================================================================

rm(list = ls())
gc()

required <- c("data.table")
missing <- required[
  !vapply(required, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing) > 0L) {
  stop(
    "Missing package(s): ", paste(missing, collapse = ", "),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(data.table)
})

PROJECT_DIR <- Sys.getenv(
  "SPINE_TRACT_PROJECT_DIR",
  unset = "C:/Users/garin/OneDrive/OHSU/Research/TQIP/TQIP Spine Utilization"
)

CLASS_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_final_classification_REPAIRED_DX",
  "rolling_origin_gpu"
)
DURATION_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_duration_quantiles_REPAIRED_DX",
  "rolling_origin_gpu"
)
RIDGE_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage3_ridge_comparator_REPAIRED_DX_FINAL"
)
ROBUST_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage4_robustness_REPAIRED_DX",
  "fixed_predictions_v1"
)
OUT_DIR <- file.path(
  PROJECT_DIR,
  "outputs",
  "stage5_manuscript_outputs_REPAIRED_DX",
  "frozen_predictions_v1"
)

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

B <- 1000L
SEED <- 20260929L
DURATION_TOL <- 1e-6
set.seed(SEED)

LOGICAL_CORES <- parallel::detectCores(logical = TRUE)
if (is.na(LOGICAL_CORES) || LOGICAL_CORES < 2L) LOGICAL_CORES <- 2L

N_WORKERS <- suppressWarnings(
  as.integer(Sys.getenv("SPINE_TRACT_BOOTSTRAP_WORKERS", unset = "20"))
)
if (is.na(N_WORKERS) || N_WORKERS < 1L) N_WORKERS <- 1L
N_WORKERS <- min(N_WORKERS, max(1L, LOGICAL_CORES - 2L))

BOOT_CL <- parallel::makeCluster(N_WORKERS, type = "PSOCK")
on.exit(parallel::stopCluster(BOOT_CL), add = TRUE)

cat(
  "\nBootstrap workers: ",
  N_WORKERS, "/", LOGICAL_CORES,
  " logical cores\n",
  sep = ""
)

# ------------------------------------------------------------------------------
# 1. Metric helpers
# ------------------------------------------------------------------------------

clamp_prob <- function(p, eps = 1e-7) {
  pmin(pmax(p, eps), 1 - eps)
}

fast_auc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  n1 <- as.double(sum(y == 1L))
  n0 <- as.double(sum(y == 0L))

  if (n1 == 0L || n0 == 0L) return(NA_real_)

  r <- rank(p, ties.method = "average")
  (
    sum(r[y == 1L]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}

fast_auprc <- function(y, p) {
  ok <- !is.na(y) & !is.na(p)
  y <- as.integer(y[ok])
  p <- as.numeric(p[ok])

  npos <- sum(y == 1L)
  if (npos == 0L) return(NA_real_)

  ord <- order(p, decreasing = TRUE)
  ys <- y[ord]
  precision <- cumsum(ys == 1L) / seq_along(ys)

  mean(precision[ys == 1L])
}

binary_point_metrics <- function(y, p) {
  y <- as.integer(y)
  p <- as.numeric(p)
  pp <- clamp_prob(p)
  prev <- mean(y == 1L)

  fit_slope <- tryCatch(
    suppressWarnings(
      glm(
        y ~ qlogis(pp),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  fit_int <- tryCatch(
    suppressWarnings(
      glm(
        y ~ offset(qlogis(pp)),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  data.table(
    n = length(y),
    prevalence = prev,
    auroc = fast_auc(y, p),
    auprc = fast_auprc(y, p),
    pr_lift = if (prev > 0) fast_auprc(y, p) / prev else NA_real_,
    brier = mean((y - p)^2),
    logloss = -mean(
      y * log(pp) +
        (1 - y) * log(1 - pp)
    ),
    calibration_intercept = if (is.null(fit_int)) {
      NA_real_
    } else {
      unname(coef(fit_int)[1L])
    },
    calibration_slope = if (is.null(fit_slope)) {
      NA_real_
    } else {
      unname(coef(fit_slope)[2L])
    }
  )
}

ci <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) {
    return(c(lower = NA_real_, upper = NA_real_))
  }
  unname(
    quantile(
      x,
      probs = c(0.025, 0.975),
      na.rm = TRUE,
      type = 6
    )
  )
}

# Patient-level ordinary nonparametric bootstrap.
# Parallelized over a persistent Windows PSOCK cluster.
bootstrap_binary <- function(y, p, B, seed, cl) {
  y <- as.integer(y)
  p <- as.numeric(p)
  n <- length(y)

  chunks <- split(
    seq_len(B),
    rep(seq_along(cl), length.out = B)
  )

  parallel::clusterSetRNGStream(cl, iseed = seed)

  pieces <- parallel::parLapply(
    cl,
    chunks,
    function(ids, y, p, n) {
      clamp_local <- function(x, eps = 1e-7) {
        pmin(pmax(x, eps), 1 - eps)
      }

      auc_local <- function(yy, pp) {
        n1 <- sum(yy == 1L)
        n0 <- sum(yy == 0L)
        if (n1 == 0L || n0 == 0L) return(NA_real_)
        r <- rank(pp, ties.method = "average")
        (
          sum(r[yy == 1L]) -
            n1 * (n1 + 1) / 2
        ) / (n1 * n0)
      }

      auprc_local <- function(yy, pp) {
        npos <- sum(yy == 1L)
        if (npos == 0L) return(NA_real_)
        ord <- order(pp, decreasing = TRUE)
        ys <- yy[ord]
        precision <- cumsum(ys == 1L) / seq_along(ys)
        mean(precision[ys == 1L])
      }

      out <- matrix(
        NA_real_,
        nrow = length(ids),
        ncol = 4L,
        dimnames = list(
          NULL,
          c("auroc", "auprc", "brier", "logloss")
        )
      )

      for (j in seq_along(ids)) {
        idx <- sample.int(n, size = n, replace = TRUE)
        yy <- y[idx]
        pp <- p[idx]
        cpp <- clamp_local(pp)

        out[j, "auroc"] <- auc_local(yy, pp)
        out[j, "auprc"] <- auprc_local(yy, pp)
        out[j, "brier"] <- mean((yy - pp)^2)
        out[j, "logloss"] <- -mean(
          yy * log(cpp) +
            (1 - yy) * log(1 - cpp)
        )
      }

      out
    },
    y = y,
    p = p,
    n = n
  )

  as.data.table(do.call(rbind, pieces))
}

# ------------------------------------------------------------------------------
# 2. Canonical classification-output loader
# ------------------------------------------------------------------------------

load_binary_output <- function(
    year,
    id,
    label
) {
  pred_dir <- file.path(
    CLASS_DIR,
    "predictions"
  )

  if (id == "home") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "discharge3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )
    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "disposition_home",
        output_label = "Home/home health",
        y = as.integer(
          d$observed_class == "Home/home health"
        ),
        p = d$p_Home.home.health
      )
    )
  }

  if (id == "postacute") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "discharge3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )
    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "disposition_postacute",
        output_label = "Post-acute facility",
        y = as.integer(
          d$observed_class == "Post-acute facility"
        ),
        p = d$p_Post.acute.facility
      )
    )
  }

  if (id == "death") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "discharge3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )
    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "disposition_death_hospice",
        output_label = "Death/hospice",
        y = as.integer(
          d$observed_class == "Death/hospice"
        ),
        p = d$p_Death.hospice
      )
    )
  }

  if (id == "hlos20") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "hlos3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )
    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "hlos_ge20",
        output_label = "Hospital LOS >=20 days",
        y = as.integer(
          d$observed_class == ">=20 days"
        ),
        p = d$p_X..20.days
      )
    )
  }

  if (id %in% c("icu_any", "icu8")) {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "icu3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    if (id == "icu_any") {
      return(
        data.table(
          admission_year = d$admission_year,
          inc_key = d$inc_key,
          output_id = "icu_any",
          output_label = "Any ICU use",
          y = as.integer(
            d$observed_class != "No ICU"
          ),
          p = d$p_ICU.1.7.days +
            d$p_ICU...8.days
        )
      )
    }

    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "icu_ge8",
        output_label = "ICU >=8 days",
        y = as.integer(
          d$observed_class == "ICU >=8 days"
        ),
        p = d$p_ICU...8.days
      )
    )
  }

  if (id %in% c("vent_any", "vent8")) {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "vent3_final__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    if (id == "vent_any") {
      return(
        data.table(
          admission_year = d$admission_year,
          inc_key = d$inc_key,
          output_id = "vent_any",
          output_label = "Any invasive mechanical ventilation",
          y = as.integer(
            d$observed_class != "No ventilation"
          ),
          p = d$p_Ventilation.1.7.days +
            d$p_Ventilation...8.days
        )
      )
    }

    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "vent_ge8",
        output_label = "Ventilation >=8 days",
        y = as.integer(
          d$observed_class == "Ventilation >=8 days"
        ),
        p = d$p_Ventilation...8.days
      )
    )
  }

  if (id == "surgery") {
    d <- fread(
      file.path(
        pred_dir,
        paste0(
          "operative_spine_intervention__eval",
          year,
          "__predictions.csv"
        )
      )
    )

    return(
      data.table(
        admission_year = d$admission_year,
        inc_key = d$inc_key,
        output_id = "recorded_tl_operation",
        output_label =
          "Recorded thoracolumbar operative treatment",
        y = as.integer(d$y),
        p = as.numeric(d$pred)
      )
    )
  }

  stop("Unknown output id: ", id, call. = FALSE)
}

output_specs <- data.table(
  id = c(
    "home",
    "postacute",
    "death",
    "hlos20",
    "icu_any",
    "icu8",
    "vent_any",
    "vent8",
    "surgery"
  ),
  output_id = c(
    "disposition_home",
    "disposition_postacute",
    "disposition_death_hospice",
    "hlos_ge20",
    "icu_any",
    "icu_ge8",
    "vent_any",
    "vent_ge8",
    "recorded_tl_operation"
  )
)

# ------------------------------------------------------------------------------
# 3. Temporal point estimates and pooled OOF point estimates
# ------------------------------------------------------------------------------

temporal_rows <- list()
pooled_rows <- list()

for (i in seq_len(nrow(output_specs))) {
  id <- output_specs$id[i]

  per_year <- list()

  for (year in 2022:2024) {
    d <- load_binary_output(
      year,
      id,
      id
    )

    m <- binary_point_metrics(
      d$y,
      d$p
    )
    m[, `:=`(
      output_id = d$output_id[1L],
      output_label = d$output_label[1L],
      evaluation_year = year
    )]

    temporal_rows[[length(temporal_rows) + 1L]] <- m
    per_year[[length(per_year) + 1L]] <- d
  }

  pooled <- rbindlist(
    per_year,
    use.names = TRUE
  )

  pm <- binary_point_metrics(
    pooled$y,
    pooled$p
  )
  pm[, `:=`(
    output_id = pooled$output_id[1L],
    output_label = pooled$output_label[1L],
    evaluation_year = "Pooled 2022-2024"
  )]

  pooled_rows[[length(pooled_rows) + 1L]] <- pm
}

temporal_metrics <- rbindlist(
  temporal_rows,
  fill = TRUE
)
pooled_metrics <- rbindlist(
  pooled_rows,
  fill = TRUE
)

setcolorder(
  temporal_metrics,
  c(
    "output_id",
    "output_label",
    "evaluation_year",
    setdiff(
      names(temporal_metrics),
      c(
        "output_id",
        "output_label",
        "evaluation_year"
      )
    )
  )
)

fwrite(
  temporal_metrics,
  file.path(
    OUT_DIR,
    "01_CLASSIFICATION_TEMPORAL_POINT_ESTIMATES.csv"
  )
)
fwrite(
  pooled_metrics,
  file.path(
    OUT_DIR,
    "02_CLASSIFICATION_POOLED_OOF_POINT_ESTIMATES.csv"
  )
)

# ------------------------------------------------------------------------------
# 4. 2024 classification bootstrap CIs
# ------------------------------------------------------------------------------

bootstrap_summary <- list()

for (i in seq_len(nrow(output_specs))) {
  id <- output_specs$id[i]

  cat(
    "\nBootstrap ",
    i,
    "/",
    nrow(output_specs),
    ": ",
    id,
    "\n",
    sep = ""
  )

  d <- load_binary_output(
    2024L,
    id,
    id
  )

  point <- binary_point_metrics(
    d$y,
    d$p
  )

  boot <- bootstrap_binary(
    d$y,
    d$p,
    B = B,
    seed = SEED + i,
    cl = BOOT_CL
  )

  auc_ci <- ci(boot$auroc)
  pr_ci <- ci(boot$auprc)
  brier_ci <- ci(boot$brier)
  ll_ci <- ci(boot$logloss)

  row <- data.table(
    output_id = d$output_id[1L],
    output_label = d$output_label[1L],
    evaluation_year = 2024L,
    B = B,
    n = point$n,
    prevalence = point$prevalence,
    auroc = point$auroc,
    auroc_lower = auc_ci[1L],
    auroc_upper = auc_ci[2L],
    auprc = point$auprc,
    auprc_lower = pr_ci[1L],
    auprc_upper = pr_ci[2L],
    pr_lift = point$pr_lift,
    brier = point$brier,
    brier_lower = brier_ci[1L],
    brier_upper = brier_ci[2L],
    logloss = point$logloss,
    logloss_lower = ll_ci[1L],
    logloss_upper = ll_ci[2L],
    calibration_intercept =
      point$calibration_intercept,
    calibration_slope =
      point$calibration_slope
  )

  bootstrap_summary[[length(bootstrap_summary) + 1L]] <- row

  fwrite(
    rbindlist(
      bootstrap_summary,
      fill = TRUE
    ),
    file.path(
      OUT_DIR,
      "03_CLASSIFICATION_2024_BOOTSTRAP_CI_CHECKPOINT.csv"
    )
  )
}

classification_2024_ci <- rbindlist(
  bootstrap_summary,
  fill = TRUE
)

fwrite(
  classification_2024_ci,
  file.path(
    OUT_DIR,
    "04_CLASSIFICATION_2024_BOOTSTRAP_CI.csv"
  )
)

# ------------------------------------------------------------------------------
# 5. 2024 calibration data: 20-bin + flexible logistic spline
# ------------------------------------------------------------------------------

calibration_bins <- list()
calibration_spline <- list()

for (i in seq_len(nrow(output_specs))) {
  id <- output_specs$id[i]
  d <- load_binary_output(
    2024L,
    id,
    id
  )

  d[, bin := frank(
    p,
    ties.method = "average"
  )]
  d[, bin := pmax(
    1L,
    pmin(
      20L,
      ceiling(20 * bin / .N)
    )
  )]

  bins <- d[
    ,
    .(
      n = .N,
      mean_predicted = mean(p),
      observed_rate = mean(y)
    ),
    by = bin
  ]
  bins[, `:=`(
    output_id = d$output_id[1L],
    output_label = d$output_label[1L],
    evaluation_year = 2024L
  )]

  calibration_bins[[length(calibration_bins) + 1L]] <- bins

  pp <- clamp_prob(d$p)
  lp <- qlogis(pp)

  spline_fit <- tryCatch(
    suppressWarnings(
      glm(
        d$y ~ splines::ns(lp, df = 4L),
        family = binomial()
      )
    ),
    error = function(e) NULL
  )

  if (!is.null(spline_fit)) {
    grid_p <- unique(
      as.numeric(
        quantile(
          d$p,
          probs = seq(
            0.01,
            0.99,
            length.out = 99L
          ),
          na.rm = TRUE,
          type = 8
        )
      )
    )

    grid_lp <- qlogis(
      clamp_prob(grid_p)
    )

    pred <- predict(
      spline_fit,
      newdata = data.frame(
        lp = grid_lp
      ),
      type = "response",
      se.fit = FALSE
    )

    curve <- data.table(
      output_id = d$output_id[1L],
      output_label = d$output_label[1L],
      evaluation_year = 2024L,
      predicted_probability = grid_p,
      spline_observed_probability =
        as.numeric(pred)
    )

    calibration_spline[[length(calibration_spline) + 1L]] <- curve
  }
}

fwrite(
  rbindlist(
    calibration_bins,
    fill = TRUE
  ),
  file.path(
    OUT_DIR,
    "05_CLASSIFICATION_2024_CALIBRATION_20BIN.csv"
  )
)

fwrite(
  rbindlist(
    calibration_spline,
    fill = TRUE
  ),
  file.path(
    OUT_DIR,
    "06_CLASSIFICATION_2024_CALIBRATION_SPLINE.csv"
  )
)

# ------------------------------------------------------------------------------
# 6. Duration point estimates + 2024 bootstrap CIs
# ------------------------------------------------------------------------------

duration_metrics_candidates <- c(
  file.path(
    DURATION_DIR,
    "05b_FINAL_DURATION_TEMPORAL_METRICS_DISCRETE_TOLERANCE.csv"
  ),
  file.path(
    DURATION_DIR,
    "05_FINAL_DURATION_TEMPORAL_METRICS.csv"
  )
)

duration_metrics_file <-
  duration_metrics_candidates[
    file.exists(
      duration_metrics_candidates
    )
  ][1L]

if (
  length(duration_metrics_file) == 0L ||
    is.na(duration_metrics_file)
) {
  stop(
    "Stage 3C repaired-DX duration metrics not found.",
    call. = FALSE
  )
}

duration_temporal <- fread(
  duration_metrics_file
)

fwrite(
  duration_temporal,
  file.path(
    OUT_DIR,
    "07_DURATION_TEMPORAL_POINT_ESTIMATES.csv"
  )
)

duration_specs <- data.table(
  outcome = c(
    "hospital_los",
    "icu_los_conditional",
    "ventilator_days_conditional"
  )
)

duration_bootstrap_rows <- list()

duration_metric_one <- function(d) {
  y <- as.numeric(d$observed_days)
  q10 <- as.numeric(d$q10_days)
  q50 <- as.numeric(d$q50_median_days)
  q90 <- as.numeric(d$q90_days)

  pinball <- function(y, q, alpha) {
    e <- y - q
    mean(
      ifelse(
        e >= 0,
        alpha * e,
        (alpha - 1) * e
      )
    )
  }

  c(
    mae = mean(abs(y - q50)),
    median_bias = median(q50 - y),
    mean_pinball = mean(c(
      pinball(y, q10, 0.10),
      pinball(y, q50, 0.50),
      pinball(y, q90, 0.90)
    )),
    p_y_lt_q10 =
      mean(y < q10 - DURATION_TOL),
    p_y_le_q10 =
      mean(y <= q10 + DURATION_TOL),
    p_y_lt_q50 =
      mean(y < q50 - DURATION_TOL),
    p_y_le_q50 =
      mean(y <= q50 + DURATION_TOL),
    p_y_lt_q90 =
      mean(y < q90 - DURATION_TOL),
    p_y_le_q90 =
      mean(y <= q90 + DURATION_TOL),
    coverage =
      mean(
        y + DURATION_TOL >= q10 &
          y - DURATION_TOL <= q90
      ),
    median_width =
      median(q90 - q10)
  )
}

bootstrap_duration <- function(d, B, seed, cl) {
  y <- as.numeric(d$observed_days)
  q10 <- as.numeric(d$q10_days)
  q50 <- as.numeric(d$q50_median_days)
  q90 <- as.numeric(d$q90_days)
  n <- length(y)
  tol <- DURATION_TOL

  chunks <- split(
    seq_len(B),
    rep(seq_along(cl), length.out = B)
  )

  parallel::clusterSetRNGStream(
    cl,
    iseed = seed
  )

  pieces <- parallel::parLapply(
    cl,
    chunks,
    function(
        ids,
        y,
        q10,
        q50,
        q90,
        n,
        tol
    ) {
      pinball <- function(y, q, alpha) {
        e <- y - q

        mean(
          ifelse(
            e >= 0,
            alpha * e,
            (alpha - 1) * e
          )
        )
      }

      one_metric <- function(idx) {
        yy <- y[idx]
        a <- q10[idx]
        b <- q50[idx]
        c <- q90[idx]

        c(
          mae =
            mean(
              abs(
                yy - b
              )
            ),
          median_bias =
            median(
              b - yy
            ),
          mean_pinball =
            mean(c(
              pinball(
                yy,
                a,
                0.10
              ),
              pinball(
                yy,
                b,
                0.50
              ),
              pinball(
                yy,
                c,
                0.90
              )
            )),
          p_y_lt_q10 =
            mean(
              yy < a - tol
            ),
          p_y_le_q10 =
            mean(
              yy <= a + tol
            ),
          p_y_lt_q50 =
            mean(
              yy < b - tol
            ),
          p_y_le_q50 =
            mean(
              yy <= b + tol
            ),
          p_y_lt_q90 =
            mean(
              yy < c - tol
            ),
          p_y_le_q90 =
            mean(
              yy <= c + tol
            ),
          coverage =
            mean(
              yy + tol >= a &
                yy - tol <= c
            ),
          median_width =
            median(
              c - a
            )
        )
      }

      out <- matrix(
        NA_real_,
        nrow =
          length(ids),
        ncol =
          12L,
        dimnames = list(
          NULL,
          c(
            "mae",
            "median_bias",
            "mean_pinball",
            "p_y_lt_q10",
            "p_y_le_q10",
            "p_y_lt_q50",
            "p_y_le_q50",
            "p_y_lt_q90",
            "p_y_le_q90",
            "coverage",
            "median_width",
            "unused_placeholder"
          )
        )
      )

      # Drop the placeholder column after filling the 11 real metrics.
      for (j in seq_along(ids)) {
        z <- one_metric(
          sample.int(
            n,
            size = n,
            replace = TRUE
          )
        )

        out[
          j,
          seq_along(z)
        ] <- z
      }

      out[
        ,
        colnames(out) !=
          "unused_placeholder",
        drop = FALSE
      ]
    },
    y = y,
    q10 = q10,
    q50 = q50,
    q90 = q90,
    n = n,
    tol = tol
  )

  as.data.table(
    do.call(
      rbind,
      pieces
    )
  )
}

for (i in seq_len(nrow(duration_specs))) {
  id <- duration_specs$outcome[i]

  pfile <- file.path(
    DURATION_DIR,
    "predictions",
    paste0(
      id,
      "__eval2024__Q10_Q50_Q90_predictions.csv"
    )
  )

  d <- fread(pfile)
  point <- duration_metric_one(d)

  boot <- bootstrap_duration(
    d = d,
    B = B,
    seed = SEED + 100L + i,
    cl = BOOT_CL
  )

  row <- data.table(
    outcome = id,
    evaluation_year = 2024L,
    B = B,
    N = nrow(d)
  )

  for (nm in names(point)) {
    q <- ci(boot[[nm]])
    row[
      ,
      paste0(nm, "_estimate") := point[nm]
    ]
    row[
      ,
      paste0(nm, "_lower") := q[1L]
    ]
    row[
      ,
      paste0(nm, "_upper") := q[2L]
    ]
  }

  duration_bootstrap_rows[[length(duration_bootstrap_rows) + 1L]] <- row

  fwrite(
    rbindlist(
      duration_bootstrap_rows,
      fill = TRUE
    ),
    file.path(
      OUT_DIR,
      "08_DURATION_2024_BOOTSTRAP_CI_CHECKPOINT.csv"
    )
  )
}

duration_2024_ci <- rbindlist(
  duration_bootstrap_rows,
  fill = TRUE
)

fwrite(
  duration_2024_ci,
  file.path(
    OUT_DIR,
    "09_DURATION_2024_BOOTSTRAP_CI.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Locked comparator / robustness source copies
# ------------------------------------------------------------------------------

ridge_source <- file.path(
  RIDGE_DIR,
  "06_FINAL_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv"
)

if (file.exists(ridge_source)) {
  file.copy(
    ridge_source,
    file.path(
      OUT_DIR,
      "10_FINAL_XGBOOST_VS_RIDGE_HEADLINE_COMPARISON.csv"
    ),
    overwrite = TRUE
  )
}

robust_source <- file.path(
  ROBUST_DIR,
  "06_2024_SUBGROUP_RANGE_SUMMARY.csv"
)

if (file.exists(robust_source)) {
  file.copy(
    robust_source,
    file.path(
      OUT_DIR,
      "11_2024_SUBGROUP_RANGE_SUMMARY.csv"
    ),
    overwrite = TRUE
  )
}

source_manifest <- data.table(
  source = c(
    "Classification prediction directory",
    "Duration prediction directory",
    "Final Ridge comparator",
    "Robustness suite"
  ),
  path = c(
    CLASS_DIR,
    DURATION_DIR,
    ridge_source,
    robust_source
  ),
  exists = c(
    dir.exists(CLASS_DIR),
    dir.exists(DURATION_DIR),
    file.exists(ridge_source),
    file.exists(robust_source)
  )
)

fwrite(
  source_manifest,
  file.path(
    OUT_DIR,
    "12_SOURCE_ARTIFACT_MANIFEST.csv"
  )
)

# ------------------------------------------------------------------------------
# 7. Manifest / QC
# ------------------------------------------------------------------------------

manifest <- data.table(
  item = c(
    "Primary population",
    "Disposition",
    "Classification evaluation years",
    "Duration evaluation years",
    "Bootstrap replicates",
    "Bootstrap interpretation",
    "Classification calibration curve",
    "Duration range terminology",
    "Discrete duration comparison tolerance"
  ),
  specification = c(
    "Direct presentations only; repaired frozen primary cohort",
    "Home/home health; post-acute facility; death/hospice",
    "2022, 2023, 2024 repaired-DX rolling-origin fixed prediction caches",
    "2022, 2023, 2024 repaired-DX rolling-origin fixed prediction caches",
    as.character(B),
    "Patient-level nonparametric resampling of fixed evaluation predictions; not model-refit stability",
    "2024 20-bin observed/predicted plus 4-df logistic natural-spline curve",
    "Estimated 10th-90th percentile range; not a guaranteed 80% prediction interval",
    paste0(DURATION_TOL, " day; used only for numerical equality in discrete-day duration calibration")
  )
)

fwrite(
  manifest,
  file.path(
    OUT_DIR,
    "13_MANUSCRIPT_OUTPUT_MANIFEST.csv"
  )
)

writeLines(
  c(
    "SPINE-TRACT REPAIRED-DX MANUSCRIPT PREDICTION OUTPUTS COMPLETE",
    "",
    "NO MODEL WAS REFIT.",
    "",
    paste0(
      "Bootstrap B=",
      B,
      " using fixed 2024 prediction caches."
    ),
    "",
    "Review first:",
    "  04_CLASSIFICATION_2024_BOOTSTRAP_CI.csv",
    "  02_CLASSIFICATION_POOLED_OOF_POINT_ESTIMATES.csv",
    "  06_CLASSIFICATION_2024_CALIBRATION_SPLINE.csv",
    "  09_DURATION_2024_BOOTSTRAP_CI.csv",
    "  07_DURATION_TEMPORAL_POINT_ESTIMATES.csv"
  ),
  file.path(
    OUT_DIR,
    "README_STAGE3D_RESULTS.txt"
  )
)

cat("\n============================================================\n")
cat("REPAIRED-DX MANUSCRIPT OUTPUTS COMPLETE — NO MODEL REFIT\n")
cat("============================================================\n")
cat(
  "Output: ",
  OUT_DIR,
  "\n",
  sep = ""
)
