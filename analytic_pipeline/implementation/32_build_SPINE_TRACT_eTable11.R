# =============================================================================
# SPINE-TRACT / Stage 32: TRIPOD+AI subgroup performance uncertainty (eTable 11)
# Locked repaired-DX models; 2024 FIXED PREDICTIONS ONLY. NO MODEL REFITS.
#
# Matches the subgroup derivations and classification prediction columns in:
# 08a_run_SPINE_TRACT_fixed_prediction_robustness_REPAIRED_DX_v1.R
# and the key endpoints/domains summarized in publication eTable 8A.
#
# Statistical methods:
#  * AUROC, average precision (AUPRC), Brier: patient-level nonparametric
#    percentile bootstrap (B=1000 by default) of fixed prediction/outcome pairs.
#  * Calibration-in-the-large: glm(y ~ offset(logit(p))), binomial.
#  * Calibration slope: glm(y ~ logit(p)), binomial.
#  * Calibration 95% CIs: normal/Wald coefficient +/- 1.96 SE (NOT bootstrap).
#  * Missing outcome / nonfinite prediction pairs excluded endpoint-specifically.
#  * Groups with N<500 or <25 cases/noncases are not estimable (logged).
#  * Checkpoints keyed by B, seed, code-method version, and input-file MD5.
#
# Required: data.table. Optional for .docx: officer, flextable.
# Set SPINE_TRACT_PROJECT_DIR before running (the licensed local workspace).
# Optional: SPINE_TRACT_ETABLE11_B (default 1000),
#           SPINE_TRACT_ETABLE11_WORKERS (default 4),
#           SPINE_TRACT_ETABLE11_SEED (default 20261018).
# Use B=1000 for manuscript publication. B<1000 is a development test only.
# =============================================================================

suppressPackageStartupMessages(library(data.table))

METHOD_VERSION <- "2026-10-18_v1"
PROJECT_DIR <- Sys.getenv("SPINE_TRACT_PROJECT_DIR", unset = "")
if (!nzchar(PROJECT_DIR) || !dir.exists(PROJECT_DIR)) {
  stop("Set SPINE_TRACT_PROJECT_DIR to the existing licensed SPINE-TRACT project folder.")
}
B <- as.integer(Sys.getenv("SPINE_TRACT_ETABLE11_B", "1000"))
WORKERS <- as.integer(Sys.getenv("SPINE_TRACT_ETABLE11_WORKERS", "4"))
SEED <- as.integer(Sys.getenv("SPINE_TRACT_ETABLE11_SEED", "20261018"))
if (is.na(B) || B < 20L) stop("Invalid bootstrap B (must be >=20).")
if (is.na(WORKERS) || WORKERS < 1L) stop("Invalid worker count.")
if (is.na(SEED) || SEED < 1L) stop("Invalid seed.")
ncores <- parallel::detectCores(logical = TRUE)
if (is.na(ncores) || ncores < 1L) ncores <- WORKERS
WORKERS <- min(WORKERS, ncores, B)

DATA_DIR <- file.path(PROJECT_DIR, "data", "derived", "stage2_predictor_adjudication")
PRIMARY_FILE <- file.path(DATA_DIR, "spine_tract_frozen_primary_2020_2024.rds")
SPEC_FILE <- file.path(DATA_DIR, "spine_tract_predictor_sets.rds")
CLASS_DIR <- file.path(PROJECT_DIR, "outputs", "stage3_final_classification_REPAIRED_DX", "rolling_origin_gpu", "predictions")
ROBUST_FILE <- file.path(PROJECT_DIR, "outputs", "stage4_robustness_REPAIRED_DX", "fixed_predictions_v1", "03_CLASSIFICATION_SUBGROUP_ROBUSTNESS.csv")

output_specs <- list(
  disposition_postacute = list(endpoint="discharge3_final", label="Post-acute facility", target="Post-acute facility", col="p_Post.acute.facility"),
  disposition_death_hospice = list(endpoint="discharge3_final", label="Death/hospice", target="Death/hospice", col="p_Death.hospice"),
  hlos_ge20 = list(endpoint="hlos3_final", label="Hospital LOS >=20 days", target=">=20 days", col="p_X..20.days"),
  icu_ge8 = list(endpoint="icu3_final", label="ICU >=8 days", target="ICU >=8 days", col="p_ICU...8.days"),
  vent_ge8 = list(endpoint="vent3_final", label="Ventilation >=8 days", target="Ventilation >=8 days", col="p_Ventilation...8.days"),
  recorded_tl_operation = list(endpoint="operative_spine_intervention", label="Recorded thoracolumbar operation", target="EVENT", col="pred")
)
# Main eTable 8A clinical domains and outcomes; plus selected system-mediated
# output by contextual group. Individual-level ALL results are also saved.
clinical_domains <- c("Age", "Sex", "GCS severity", "Clinical input missingness", "TL region", "SCI", "Neural pattern", "Concomitant cervical injury", "Concomitant TBI")
clinical_outputs <- c("disposition_postacute", "disposition_death_hospice", "icu_ge8", "vent_ge8", "recorded_tl_operation")
context_domains <- c("Race", "Ethnicity", "Payer")
context_outputs <- c("disposition_postacute", "hlos_ge20")

pred_files <- setNames(file.path(CLASS_DIR, paste0(unique(vapply(output_specs, `[[`, "", "endpoint")), "__eval2024__predictions.csv")),
                       unique(vapply(output_specs, `[[`, "", "endpoint")))
required <- c(PRIMARY_FILE, SPEC_FILE, unname(pred_files), ROBUST_FILE)
if (!all(file.exists(required))) {
  stop("Missing locally generated prerequisite(s):\n", paste(required[!file.exists(required)], collapse="\n"),
       "\nRun the final classification and Stage 08a robustness outputs before eTable 11.")
}

OUT_DIR <- file.path(PROJECT_DIR, "outputs", "SPINE_TRACT_PUBLICATION_TABLES", "eTable11_TRIPOD_AI")
CHECK_DIR <- file.path(OUT_DIR, sprintf("checkpoints_%s_B%d_seed%d", METHOD_VERSION, B, SEED))
dir.create(CHECK_DIR, recursive=TRUE, showWarnings=FALSE)

# MD5 verifies that a resumed run uses the SAME frozen cohort and predictions.
input_hash <- unname(tools::md5sum(required))
manifest <- list(method=METHOD_VERSION, B=B, seed=SEED,
                 input_files=normalizePath(required, winslash="/", mustWork=TRUE),
                 md5=input_hash,
                 CI_methods=c("AUROC/AUPRC/Brier: fixed-prediction percentile bootstrap",
                              "Calibration intercept/slope: GLM Wald"))
manifest_path <- file.path(CHECK_DIR, "CHECKPOINT_MANIFEST.rds")
if (file.exists(manifest_path)) {
  old <- readRDS(manifest_path)
  if (!identical(old, manifest)) {
    stop("Checkpoint manifest mismatch (B/seed/input version/files changed). Choose a fresh checkpoint directory; existing results were NOT reused.")
  }
} else {
  saveRDS(manifest, manifest_path)
}
cat("eTable11 | B=", B, " | workers=", WORKERS, " | checkpoint=", CHECK_DIR, "\n", sep="")

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
norm_chr <- function(x) {
  z <- trimws(as.character(x))
  z[is.na(z) | z=="" | z %in% c("NA", "NaN", "<NA>")] <- "__UNKNOWN__"
  z
}

# Derive strata using the definitions in the repaired-DX Stage 08a pipeline.
raw <- as.data.table(readRDS(PRIMARY_FILE))
stopifnot(nrow(raw)==424017L)
stopifnot(all(c("inc_key", "admission_year") %in% names(raw)))
raw[, admission_year := as.integer(admission_year)]
d <- raw[admission_year==2024L]
rm(raw); invisible(gc())
if (nrow(d)!=90649L) stop("2024 primary count changed: ",nrow(d)," (expected 90649).")
d[, inc_key := as.character(inc_key)]
if (anyNA(d$inc_key) || anyDuplicated(d$inc_key)) stop("Missing or duplicate 2024 inc_key.")

sp <- readRDS(SPEC_FILE)
pg <- sp$predictor_groups
clinical_backbone <- sp$predictor_sets$clinical_backbone
final_core <- unique(c(setdiff(clinical_backbone, pg$spine_granular), pg$spine_coarse, pg$chest_detail, pg$pupil))
if (length(final_core)!=54L || !"dx_tl_soft_tissue_any" %in% final_core) stop("Final repaired-DX 54-predictor policy mismatch.")
if (!all(c(final_core,"age_years","sex_clean","race_clean","ethnicity_clean","payer_clean") %in% names(d))) stop("Missing source subgroup/predictor fields.")

ge <- num(d$gcs_eye_clean); gv <- num(d$gcs_verbal_clean); gm <- num(d$gcs_motor_clean)
gt <- rep(NA_real_,nrow(d)); gok <- !is.na(ge)&!is.na(gv)&!is.na(gm)
gt[gok] <- ge[gok]+gv[gok]+gm[gok]
d[, gcs_severity_sens := fcase(!is.na(gt) & gt >= 13 & gt <= 15, "Mild GCS 13-15",
                              !is.na(gt) & gt >= 9 & gt <= 12, "Moderate GCS 9-12",
                              !is.na(gt) & gt >= 3 & gt <= 8, "Severe GCS 3-8",
                              default="GCS unknown/unclassifiable")]
agev <- num(d$age_years)
d[, age_group_sens := fcase(!is.na(agev)&agev>=18&agev<=39,"18-39",
                           !is.na(agev)&agev>=40&agev<=64,"40-64",
                           !is.na(agev)&agev>=65&agev<=79,"65-79",
                           !is.na(agev)&agev>=80&agev<=89,"80-89",default="Age unknown")]
miss_count <- integer(nrow(d))
for (v in final_core) {
  x <- d[[v]]
  if (is.numeric(x) || is.integer(x)) {
    miss_count <- miss_count + as.integer(is.na(x) | !is.finite(as.numeric(x)))
  } else {
    miss_count <- miss_count + as.integer(norm_chr(x) %in% c("__UNKNOWN__", "__OTHER__"))
  }
}
d[, missingness_sens := fcase(miss_count==0L,"0 missing/unknown",miss_count==1L,"1 missing/unknown",
                              miss_count==2L,"2 missing/unknown",miss_count>=3L,">=3 missing/unknown",default="Unknown")]
thoracic <- d$dx_thoracic_vertebral_fracture==1L | d$dx_thoracic_spinal_cord_injury==1L | d$dx_thoracic_nerve_root_injury==1L
lumbar <- d$dx_lumbar_vertebral_fracture==1L | d$dx_lumbar_spinal_cord_injury==1L | d$dx_lumbar_sacral_nerve_root_injury==1L | d$dx_cauda_equina_injury==1L
d[, spine_region_sens := fcase(thoracic & lumbar,"Thoracic + lumbar",thoracic & !lumbar,"Thoracic only",
                              !thoracic & lumbar,"Lumbar only",default="Region unclassifiable")]
fracture <- d$dx_thoracic_vertebral_fracture==1L | d$dx_lumbar_vertebral_fracture==1L
neural <- d$dx_tl_spinal_cord_injury_any==1L | d$dx_tl_root_cauda_any==1L
d[, fracture_neural_sens := fcase(fracture & !neural,"Fracture only",neural,"Neural injury present",default="Other qualifying phenotype")]
d[, sci_sens := fifelse(dx_tl_spinal_cord_injury_any==1L,"SCI present","No TL SCI")]
d[, neural_pattern_sens := fcase(dx_tl_spinal_cord_injury_any==1L & dx_tl_root_cauda_any==1L,"Cord + root/cauda",
                                dx_tl_spinal_cord_injury_any==1L,"Cord injury",
                                dx_tl_root_cauda_any==1L,"Root/cauda injury",default="No neural injury")]
d[, cervical_sens := fifelse(dx_cervical_spine_injury==1L,"Concomitant cervical injury","No cervical injury")]
d[, tbi_sens := fifelse(dx_tbi_any==1L,"Concomitant TBI","No TBI")]
subgroup_vars <- c("Age"="age_group_sens","Sex"="sex_clean","GCS severity"="gcs_severity_sens",
  "Clinical input missingness"="missingness_sens","TL region"="spine_region_sens",
  "SCI"="sci_sens","Neural pattern"="neural_pattern_sens",
  "Concomitant cervical injury"="cervical_sens","Concomitant TBI"="tbi_sens",
  "Race"="race_clean","Ethnicity"="ethnicity_clean","Payer"="payer_clean")
subgroups <- copy(d[, c("inc_key", unname(subgroup_vars)), with=FALSE])
for (v in unname(subgroup_vars)) set(subgroups, j=v, value=norm_chr(subgroups[[v]]))
rm(d,sp,pg,clinical_backbone,ge,gv,gm,gt,gok,agev,miss_count,thoracic,lumbar,fracture,neural); invisible(gc())

# Compute identical metrics to the Stage 08a robustness analysis.
fast_auc <- function(y,p) {
  n1 <- sum(y==1L); n0 <- sum(y==0L)
  if (!n1 || !n0) return(NA_real_)
  r <- rank(p, ties.method="average")
  (sum(r[y==1L]) - n1*(n1+1)/2)/(n1*n0)
}
fast_auprc <- function(y,p) {
  if (!any(y==1L)) return(NA_real_)
  yy <- y[order(p, decreasing=TRUE)]
  mean((cumsum(yy==1L)/seq_along(yy))[yy==1L])
}
wald_calibration <- function(y,p) {
  lp <- qlogis(pmin(pmax(p,1e-7),1-1e-7))
  fi <- tryCatch(suppressWarnings(glm(y ~ offset(lp),family=binomial())),error=function(e) NULL)
  fs <- tryCatch(suppressWarnings(glm(y ~ lp,family=binomial())),error=function(e) NULL)
  pull <- function(f,coef_name) {
    if (is.null(f) || !isTRUE(f$converged)) return(c(est=NA_real_,low=NA_real_,high=NA_real_))
    z <- coef(summary(f))
    if (!coef_name %in% rownames(z)) return(c(est=NA_real_,low=NA_real_,high=NA_real_))
    est <- unname(z[coef_name,"Estimate"]); se <- unname(z[coef_name,"Std. Error"])
    if (!is.finite(est) || !is.finite(se)) return(c(est=NA_real_,low=NA_real_,high=NA_real_))
    c(est=est,low=est-qnorm(.975)*se,high=est+qnorm(.975)*se)
  }
  list(intercept=pull(fi,"(Intercept)"),slope=pull(fs,"lp"))
}

bootstrap_chunk <- function(reps,y,p,seeds) {
  n <- length(y)
  z <- matrix(NA_real_,nrow=length(reps),ncol=3L,dimnames=list(NULL,c("AUROC","AUPRC","Brier")))
  for (i in seq_along(reps)) {
    set.seed(seeds[reps[i]]) # independent replicate seeds; worker count invariant
    k <- sample.int(n,size=n,replace=TRUE)
    yy <- y[k]; pp <- p[k]
    if (sum(yy==1L)==0L || sum(yy==0L)==0L) next
    z[i,] <- c(fast_auc(yy,pp),fast_auprc(yy,pp),mean((yy-pp)^2))
  }
  z
}

cl <- NULL
if (WORKERS>1L) {
  cl <- parallel::makeCluster(WORKERS, type="PSOCK")
  parallel::clusterExport(cl,c("fast_auc","fast_auprc","bootstrap_chunk"),envir=environment())
}
# Stop cluster even when script errors; run inside a top-level function below.
bootstrap_ci <- function(y,p,seed) {
  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max,B,replace=FALSE)
  chunks <- split(seq_len(B),rep(seq_len(WORKERS),length.out=B))
  pieces <- if (is.null(cl)) {
    lapply(chunks,bootstrap_chunk,y=y,p=p,seeds=seeds)
  } else {
    parallel::parLapply(cl,chunks,bootstrap_chunk,y=y,p=p,seeds=seeds)
  }
  bm <- do.call(rbind,pieces)
  valid <- colSums(is.finite(bm))
  low <- high <- setNames(rep(NA_real_,3L),colnames(bm))
  for (m in colnames(bm)) {
    if (valid[m]>=ceiling(0.95*B)) {
      q <- quantile(bm[,m],c(.025,.975),na.rm=TRUE,names=FALSE,type=7)
      low[m] <- q[1L]; high[m] <- q[2L]
    }
  }
  list(low=low,high=high,valid=valid)
}

MIN_N <- 500L; MIN_EVENTS <- 25L
all_rows <- list(); qc_rows <- list(); ijob <- 0L
old <- fread(ROBUST_FILE)
old <- old[evaluation_year==2024L]

# Helper: group name enters file key with a checksum, not untrusted raw filename.
key_id <- function(s) {
  f <- tempfile("spine_etable11_key_")
  on.exit(unlink(f), add=TRUE)
  writeLines(enc2utf8(s), f, useBytes=TRUE)
  unname(tools::md5sum(f))
}
get_y_p <- function(pred,spec) {
  required <- c("inc_key",spec$col)
  if (spec$endpoint!="operative_spine_intervention") required <- c(required,"observed_class") else required <- c(required,"y")
  if (!all(required %in% names(pred))) stop("Missing prediction columns: ",paste(setdiff(required,names(pred)),collapse=","))
  obs <- if (spec$endpoint=="operative_spine_intervention") as.integer(pred$y) else as.character(pred$observed_class)
  y <- if (spec$endpoint=="operative_spine_intervention") obs else as.integer(obs==spec$target)
  data.table(inc_key=as.character(pred$inc_key), y=y, p=as.numeric(pred[[spec$col]]))
}

run <- function() {
  for (out_id in names(output_specs)) {
    spc <- output_specs[[out_id]]
    pred <- fread(pred_files[[spc$endpoint]])
    if ("admission_year" %in% names(pred) && any(num(pred$admission_year)!=2024L,na.rm=TRUE)) stop("Wrong prediction year: ",spc$endpoint)
    yy <- get_y_p(pred,spc)
    if (anyDuplicated(yy$inc_key) || anyNA(yy$inc_key)) stop("Missing/duplicate prediction key for ",out_id)
    merged <- merge(yy,subgroups,by="inc_key",all.x=TRUE,sort=FALSE)
    if (nrow(merged)!=nrow(yy) || anyNA(merged[[unname(subgroup_vars[1L])]])) stop("Prediction-to-cohort join mismatch for ",out_id)
    if (any(is.finite(merged$p)&(merged$p < -1e-9 | merged$p > 1+1e-9))) stop("Probabilities outside [0,1]: ",out_id)
    # No silent clipping of non-probability predictions.
    if (any(merged$p[is.finite(merged$p)]<0 | merged$p[is.finite(merged$p)]>1)) stop("Invalid probabilities: ",out_id)
    domains <- if (out_id %in% clinical_outputs) clinical_domains else character()
    if (out_id %in% context_outputs) domains <- c(domains,context_domains)
    if (!length(domains)) next
    for (domain in domains) {
      v <- unname(subgroup_vars[[domain]])
      for (lev in sort(unique(merged[[v]]))) {
        # Stage 08a applies minimum counts before metric-specific complete pairs.
        z_raw <- merged[get(v)==lev]
        raw_n <- nrow(z_raw)
        raw_events <- sum(z_raw$y==1L, na.rm=TRUE)
        raw_nonevents <- sum(z_raw$y==0L, na.rm=TRUE)
        if (raw_n<MIN_N || raw_events<MIN_EVENTS || raw_nonevents<MIN_EVENTS) {
          qc_rows[[length(qc_rows)+1L]] <<- data.table(output_id=out_id,subgroup_domain=domain,subgroup_level=lev,N=raw_n,events=raw_events,non_events=raw_nonevents,status="excluded_minimum_counts")
          next
        }
        z <- z_raw[!is.na(y) & y %in% c(0L,1L) & is.finite(p)]
        n <- nrow(z); ne <- sum(z$y==1L); nn <- n-ne
        if (ne==0L || nn==0L) {
          qc_rows[[length(qc_rows)+1L]] <<- data.table(output_id=out_id,subgroup_domain=domain,subgroup_level=lev,N=n,events=ne,non_events=nn,status="excluded_no_binary_variation")
          next
        }
        k <- paste(out_id,domain,lev,sep="|")
        ck <- file.path(CHECK_DIR,paste0("row_",key_id(k),".rds"))
        if (file.exists(ck)) {
          obj <- readRDS(ck)
          if (!identical(obj$key,k) || !identical(obj$B,B) || !identical(obj$method,METHOD_VERSION)) stop("Checkpoint key/method mismatch: ",ck)
          row <- obj$row
          if (row$N != n || row$events != ne) stop("Checkpoint subgroup count mismatch: ",k)
          cat("RESUME ",k,"\n",sep="")
        } else {
          ijob <<- ijob+1L
          seed_group <- as.integer((SEED + sum(utf8ToInt(k)*seq_along(utf8ToInt(k)))) %% 2147483646L)+1L
          cat("BOOT  ",k," | N=",n," events=",ne," B=",B,"\n",sep="")
          y <- z$y; p <- z$p
          cals <- wald_calibration(y,p)
          boot <- bootstrap_ci(y,p,seed=seed_group)
          row <- data.table(
            evaluation_year=2024L,output_id=out_id,output_label=spc$label,
            subgroup_domain=domain,subgroup_level=lev,
            N=n,events=ne,non_events=nn,prevalence=ne/n,
            AUROC=fast_auc(y,p),AUROC_low=boot$low["AUROC"],AUROC_high=boot$high["AUROC"],
            AUPRC=fast_auprc(y,p),AUPRC_low=boot$low["AUPRC"],AUPRC_high=boot$high["AUPRC"],
            Brier=mean((y-p)^2),Brier_low=boot$low["Brier"],Brier_high=boot$high["Brier"],
            calibration_intercept=cals$intercept["est"],calibration_intercept_low=cals$intercept["low"],calibration_intercept_high=cals$intercept["high"],
            calibration_slope=cals$slope["est"],calibration_slope_low=cals$slope["low"],calibration_slope_high=cals$slope["high"],
            bootstrap_B=B,bootstrap_valid_AUROC=boot$valid["AUROC"],bootstrap_valid_AUPRC=boot$valid["AUPRC"],bootstrap_valid_Brier=boot$valid["Brier"],
            classification_ci="patient bootstrap percentile; fixed predictions", calibration_ci="binomial GLM Wald"
          )
          saveRDS(list(key=k,B=B,method=METHOD_VERSION,row=row),ck)
        }
        # Point-estimate parity against locked existing subgroup output.
        prev <- old[output_id==out_id & subgroup_domain==domain & subgroup_level==lev]
        if (nrow(prev)!=1L) stop("No unique 2024 Stage 08a reference for ",k)
        tests <- c(N=row$N-prev$N,events=row$events-prev$events,
                   AUROC=row$AUROC-prev$AUROC,AUPRC=row$AUPRC-prev$AUPRC,
                   Brier=row$Brier-prev$Brier,
                   cal_intercept=row$calibration_intercept-prev$calibration_intercept,
                   cal_slope=row$calibration_slope-prev$calibration_slope)
        if (any(!is.finite(tests)) || any(abs(tests)>1e-6)) {
          stop("LOCKED POINT-ESTIMATE PARITY FAILURE for ",k,". Differences: ",paste(names(tests),format(tests,digits=4),collapse="; "))
        }
        all_rows[[length(all_rows)+1L]] <<- row
      }
    }
    rm(pred,yy,merged); invisible(gc())
  }
}
# tryCatch ensures cleanup even after errors.
tryCatch(run(),finally={if (!is.null(cl)) parallel::stopCluster(cl)})

results <- rbindlist(all_rows,fill=TRUE)
qc <- rbindlist(qc_rows,fill=TRUE)
if (!nrow(results)) stop("No eligible subgroups were produced.")
setorder(results,subgroup_domain,subgroup_level,output_id)
if (B!=1000L) warning("B=",B," is a development run, NOT a manuscript-ready B=1000 analysis.")
# Explicitly mark weak estimates; no misleading intervals are filled in.
results[, complete_CI := is.finite(AUROC_low)&is.finite(AUPRC_low)&is.finite(Brier_low)&
                           is.finite(calibration_intercept_low)&is.finite(calibration_slope_low)]
if (any(!results$complete_CI)) {
  warning(sum(!results$complete_CI)," subgroup(s) have unavailable confidence limits. See numeric output; do not interpret NA as no uncertainty.")
}

fwrite(results,file.path(OUT_DIR,sprintf("eTable11_ALL_NUMERIC_B%d.csv",B)))
if (nrow(qc)) fwrite(qc,file.path(OUT_DIR,"eTable11_EXCLUDED_GROUPS_QC.csv"))

# A = clinical domains from eTable 8A; B = contextual subgroup uncertainty
A <- results[subgroup_domain %in% clinical_domains & output_id %in% clinical_outputs]
C <- results[subgroup_domain %in% context_domains & output_id %in% context_outputs]
fmt_ci <- function(est,lo,hi) {
  ifelse(is.finite(est)&is.finite(lo)&is.finite(hi),
         sprintf("%.3f (%.3f–%.3f)",est,lo,hi),
         ifelse(is.finite(est),sprintf("%.3f (NE)",est),"NE"))
}
public_table <- function(x) {
  x[,.(
    `Subgroup domain`=subgroup_domain,
    `Subgroup`=subgroup_level,
    `Outcome`=output_label,
    `N`=format(N,big.mark=",",trim=TRUE,scientific=FALSE),
    `Events`=format(events,big.mark=",",trim=TRUE,scientific=FALSE),
    `AUROC (95% CI)`=fmt_ci(AUROC,AUROC_low,AUROC_high),
    `AUPRC (95% CI)`=fmt_ci(AUPRC,AUPRC_low,AUPRC_high),
    `Brier (95% CI)`=fmt_ci(Brier,Brier_low,Brier_high),
    `Calibration intercept (95% CI)`=fmt_ci(calibration_intercept,calibration_intercept_low,calibration_intercept_high),
    `Calibration slope (95% CI)`=fmt_ci(calibration_slope,calibration_slope_low,calibration_slope_high)
  )]
}
Ap <- public_table(A); Cp <- public_table(C)
fwrite(Ap,file.path(OUT_DIR,sprintf("eTable11A_CLINICAL_B%d.csv",B)))
fwrite(Cp,file.path(OUT_DIR,sprintf("eTable11B_CONTEXTUAL_B%d.csv",B)))

if (requireNamespace("officer",quietly=TRUE) && requireNamespace("flextable",quietly=TRUE)) {
  doc <- officer::read_docx()
  landscape <- officer::prop_section(
    page_size=officer::page_size(orient="landscape",width=11,height=8.5),
    page_margins=officer::page_mar(top=.55,bottom=.55,left=.55,right=.55,gutter=0),
    type="continuous")
  doc <- officer::body_set_default_section(doc,landscape)
  add_panel <- function(doc,title,tbl,cols) {
    doc <- officer::body_add_par(doc,title,style="heading 1")
    ft <- flextable::flextable(as.data.frame(tbl[, ..cols]))
    ft <- flextable::theme_booktabs(ft)
    ft <- flextable::fontsize(ft,size=7,part="all")
    ft <- flextable::bold(ft,part="header")
    ft <- flextable::padding(ft,padding=2,part="all")
    ft <- flextable::width(ft,j=1:5,width=c(1.35,1.40,1.70,.50,.55))
    if (length(cols)==8L) ft <- flextable::width(ft,j=6:8,width=c(1.42,1.42,1.42))
    if (length(cols)==7L) ft <- flextable::width(ft,j=6:7,width=c(2.05,2.05))
    ft <- flextable::set_table_properties(ft,layout="fixed",opts_word=list(repeat_headers=TRUE,split=FALSE))
    flextable::body_add_flextable(doc,ft)
  }
  discr_cols <- c("Subgroup domain","Subgroup","Outcome","N","Events",
                  "AUROC (95% CI)","AUPRC (95% CI)","Brier (95% CI)")
  calibr_cols <- c("Subgroup domain","Subgroup","Outcome","N","Events",
                   "Calibration intercept (95% CI)","Calibration slope (95% CI)")
  doc <- add_panel(doc,"eTable 11A. Clinical subgroup discrimination and scoring, 2024",Ap,discr_cols)
  doc <- add_panel(doc,"eTable 11B. Clinical subgroup calibration, 2024",Ap,calibr_cols)
  doc <- add_panel(doc,"eTable 11C. Contextual subgroup discrimination and scoring, 2024",Cp,discr_cols)
  doc <- add_panel(doc,"eTable 11D. Contextual subgroup calibration, 2024",Cp,calibr_cols)
  doc <- officer::body_add_par(doc,
    paste0("Notes: 2024 repaired-DX direct-presentation cohort. AUROC, AUPRC, and Brier intervals: ",B,
           " fixed-prediction, patient-level percentile bootstrap resamples; calibration intercept and slope ",
           "intervals: model-based binomial GLM Wald. Models not refit; estimates do not capture ",
           "model-selection uncertainty. Eligibility: at least 500 patients, 25 events and 25 non-events. ",
           "NE = not estimable. Descriptive subgroup metrics do not establish algorithmic fairness."))
  print(doc,target=file.path(OUT_DIR,sprintf("eTable11_B%d.docx",B)))
} else {
  message("DOCX not built (officer/flextable unavailable); publication CSVs were saved.")
}

writeLines(c(
  "SPINE-TRACT eTable 11 statistical methods and QC",
  paste("Method version:",METHOD_VERSION),paste("Bootstrap replicates requested:",B),
  paste("2024 cohort:",90649L),paste("Subgroup rows:",nrow(results)),
  paste("Excluded count-limited strata:",nrow(qc)),
  paste("Missing one or more CI:",sum(!results$complete_CI)),
  "AUROC/AUPRC/Brier: fixed-prediction patient bootstrap percentile intervals.",
  "Calibration intercept/slope: logistic GLM model-based Wald intervals.",
  "Prediction-model fit/architecture was not repeated in bootstrap replicates.",
  "All eligible results checked against Stage 08a exact point estimates, tolerance 1e-6.",
  paste("Input MD5 hashes:",paste(basename(required),input_hash,sep="=",collapse="; "))
),file.path(OUT_DIR,sprintf("eTable11_QC_README_B%d.txt",B)))
cat("COMPLETE. Outputs: ",OUT_DIR,"\n",sep="")