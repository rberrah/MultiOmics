#!/usr/bin/env Rscript
# FINAL STAGE: only now open NKI ER labels, compare frozen MAINZ-only scores
# with scores offset-calibrated on the SEPARATE TRANSBIG study.
# Important disclosure: NKI has already been analyzed in earlier PR40, so
# this experiment is retrospective methodological development, not fresh
# independent confirmatory validation.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript scripts/evaluate_transbig_calibration_pmx.R outdir")
out<-args[[1L]]
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required")
source(file.path("multiomics-engine","external_validation.R"))
truth<-read.csv(file.path(out,"external-truth.csv"),stringsAsFactors=FALSE)
src<-jsonlite::read_json(file.path(out,"source-contract.json"),simplifyVector=TRUE)
contract<-jsonlite::read_json(file.path(out,"transbig-calibration-contract.json"),simplifyVector=TRUE)
if(nrow(truth)!=337L||nrow(truth)!=src$nki_n||
   anyDuplicated(truth$subject_id)||anyNA(truth$subject_id)||
   !setequal(unique(truth$outcome),c("ER_positive","ER_negative"))||
   contract$eligible_calibration_participants<150L)
 stop("Unexpected NKI external outcomes or calibration provenance")
truth$outcome<-as.integer(truth$outcome=="ER_positive")
p_train<-as.numeric(unlist(src$mainz_er_counts)[["ER_positive"]])/src$mainz_n
if(!is.finite(p_train)||p_train<=0||p_train>=1)
 stop("Training prevalence not estimable.")
clip<-function(p)pmax(1e-6,pmin(1-1e-6,p))
log_loss<-function(y,p)-mean(y*log(clip(p))+(1-y)*log(1-clip(p)))
reliability<-function(y,p){
  bins<-cut(p,breaks=seq(0,1,by=.2),include.lowest=TRUE,right=TRUE)
  totals<-as.integer(table(bins))
  empirical<-vapply(levels(bins),function(k){
    ix<-which(as.character(bins)==k)
    if(!length(ix))return(NA_real_)
    mean(y[ix])
  },numeric(1))
  means<-vapply(levels(bins),function(k){
    ix<-which(as.character(bins)==k)
    if(!length(ix))return(NA_real_)
    mean(p[ix])
  },numeric(1))
  list(
    boundaries=c(0,.2,.4,.6,.8,1),
    counts=as.list(totals),
    observed_fraction=as.list(empirical),
    mean_prediction=as.list(means),
    ece=sum(totals/length(p)*abs(empirical-means),na.rm=TRUE),
    note="Five fixed-width probability bins, descriptive only; empty bins not used in ECE."
  )
}
evaluate_arm<-function(arm){
 dir<-file.path(out,arm)
 raw<-read.csv(file.path(dir,"external-predictions.csv"),stringsAsFactors=FALSE)
 cal<-read.csv(file.path(dir,"nki-frozen-transbig-calibrated-scores.csv"),stringsAsFactors=FALSE)
 calibrator<-jsonlite::read_json(file.path(dir,"transbig-only-calibrator.json"),simplifyVector=TRUE)
 if(nrow(raw)!=nrow(truth)||nrow(cal)!=nrow(truth)||
    anyDuplicated(raw$subject_id)||anyDuplicated(cal$subject_id)||
    !setequal(truth$subject_id,raw$subject_id)||
    !setequal(truth$subject_id,cal$subject_id)||
    !identical(calibrator$fit_cohort,"TRANSBIG")||
    !identical(calibrator$score_cohort,"NKI")||
    calibrator$fit_n!=contract$eligible_calibration_participants||
    calibrator$slope_fixed!=1L ||
    !isTRUE(calibrator$nki_labels_ever_opened_in_calibration==FALSE))
  stop("Post-training calibration/NKI coverage or study isolation mismatch")
 d<-merge(merge(truth,raw,by="subject_id"),cal,
   by="subject_id",suffixes=c("_raw","_cal"))
 if(nrow(d)!=nrow(truth)||anyNA(d)||
    any(!is.finite(d$prediction_raw))||
    any(!is.finite(d$prediction_cal))||
    any(d$prediction_raw<0|d$prediction_raw>1)||
    any(d$prediction_cal<0|d$prediction_cal>1))
  stop("Missing or impossible evaluation probabilities")
 assert_replayed<-stats::plogis(stats::qlogis(clip(d$prediction_raw))+
                                 calibrator$intercept)
 if(max(abs(assert_replayed-d$prediction_cal))>1e-12)
  stop("Saved calibrator fails exact numeric replay on frozen NKI scores")
 uncal<-validate_external_predictions(
   data.frame(subject_id=d$subject_id,outcome=d$outcome,
      prediction=d$prediction_raw),
   outcome_type="binary",independent_cohort=FALSE,
   cohort_label="Retrospective NKI (previously inspected in PR40)",
   bootstrap_repetitions=300L,seed=20261011L)
 recal<-validate_external_predictions(
   data.frame(subject_id=d$subject_id,outcome=d$outcome,
      prediction=d$prediction_cal),
   outcome_type="binary",independent_cohort=FALSE,
   cohort_label="Retrospective NKI calibrated using TRANSBIG only",
   bootstrap_repetitions=300L,seed=20261011L)
 if(!identical(uncal$status,"independent_status_not_asserted")||
    !identical(recal$status,"independent_status_not_asserted")||
    !isTRUE(uncal$cohort_coverage$complete)||
    !isTRUE(recal$cohort_coverage$complete))
  stop("Retrospective repeated NKI test incorrectly classified as fresh external validation")
 baseline<-mean((d$outcome-p_train)^2)
 # Stratified paired bootstrap estimates the same NKI population's metric
 # differences, without new model fitting or choosing a recalibrator.
 pos<-which(d$outcome==1L);neg<-which(d$outcome==0L)
 set.seed(20261011L)
 changes<-replicate(500L,{
   k<-c(sample(pos,length(pos),replace=TRUE),sample(neg,length(neg),replace=TRUE))
   yy<-d$outcome[k];a<-d$prediction_raw[k];b<-d$prediction_cal[k]
   c(
     brier=mean((yy-b)^2)-mean((yy-a)^2),
     log_loss=log_loss(yy,b)-log_loss(yy,a),
     auc=validation_binary_auc(yy,b)-validation_binary_auc(yy,a)
   )
 })
 delta_ci<-function(k)as.list(unname(stats::quantile(changes[k,],c(.025,.5,.975))))
 list(
  n=nrow(d),calibrator=calibrator,
  uncalibrated=uncal$metrics,transbig_intercept_calibrated=recal$metrics,
  brier_skill_reference=baseline,
  uncalibrated_brier_skill=1-uncal$metrics$brier/baseline,
  calibrated_brier_skill=1-recal$metrics$brier/baseline,
  paired_change_calibrated_minus_uncalibrated=list(
    brier=mean((d$outcome-d$prediction_cal)^2)-
          mean((d$outcome-d$prediction_raw)^2),
    brier_bootstrap_95_ci=delta_ci("brier"),
    log_loss=log_loss(d$outcome,d$prediction_cal)-
             log_loss(d$outcome,d$prediction_raw),
    log_loss_bootstrap_95_ci=delta_ci("log_loss"),
    auc=validation_binary_auc(d$outcome,d$prediction_cal)-
        validation_binary_auc(d$outcome,d$prediction_raw),
    auc_bootstrap_95_ci=delta_ci("auc"),
    bootstrap_repetitions=500L,seed=20261011L
  ),
  uncalibrated_reliability=reliability(d$outcome,d$prediction_raw),
  calibrated_reliability=reliability(d$outcome,d$prediction_cal),
  statistical_note="Stratified paired bootstrap describes uncertainty conditional on this previously-inspected NKI cohort, not calibration transportability to a new center.",
  evaluation_has_already_been_exposed_in_pr40=TRUE,
  scientific_certification=FALSE
 )
}
primary<-evaluate_arm("primary_rank")
sensitivity<-evaluate_arm("sensitivity_raw")
report<-list(
 experiment="Retrospective three-cohort MAINZ training, TRANSBIG one-parameter calibration, NKI previously-inspected Agilent evaluation",
 mainz_training_n=src$mainz_n,
 transbig_calibration_n=contract$eligible_calibration_participants,
 nki_final_n=src$nki_n,
 common_entrez_gene_count=src$common_gene_count,
 primary=primary,sensitivity=sensitivity,
 compliance=list(training_feature_model_frozen_on_mainz=TRUE,
  calibration_parameter_fitted_only_on_transbig=TRUE,
  nki_labels_not_used_to_fit_calibrator=TRUE,
  calibration_methods_selected_after_prior_nki_research=TRUE,
  cannot_claim_fresh_independent_confirmatory_validation=TRUE,
  calibration_cohort_different_from_final_measurement_platform=TRUE),
 clinical_certification=FALSE,
 interpretation="A negative or positive calibration gain is exploratory because the NKI outcome results had been viewed before this method was selected. Final prospective/independent validation must use a new fully held-out external dataset.",
 publicationReady=FALSE
)
jsonlite::write_json(report,file.path(out,"three-cohort-retrospective-calibration-report.json"),
 pretty=TRUE,auto_unbox=TRUE,na="null",null="null")
for(arm in c("primary_rank","sensitivity_raw")){
  e<-if(arm=="primary_rank")primary else sensitivity
  delta<-e$paired_change_calibrated_minus_uncalibrated
  limits<-unlist(delta$brier_bootstrap_95_ci,use.names=FALSE)
  cat(sprintf("THREE-COHORT CALIBRATION DELTA: arm=%s n=%d TRANSBIG_intercept=%.6f Brier_change=%+.6f 95%%_bootstrap_CI=[%+.6f,%+.6f] log_loss_change=%+.6f baseline_Brier=%.6f calibration_slope_before=%.4f calibration_slope_after=%.4f\n",
   arm,e$n,e$calibrator$intercept,delta$brier,
   limits[[1L]],limits[[3L]],delta$log_loss,e$brier_skill_reference,
   e$uncalibrated$calibration_slope,
   e$transbig_intercept_calibrated$calibration_slope))
}
cat(sprintf("PMx THREE-COHORT CALIBRATION RETROSPECTIVE PASS: TRANSBIG n=%d, NKI n=%d; rank Brier before=%.4f after=%.4f AUC before=%.4f after=%.4f; raw Brier before=%.4f after=%.4f; NO NEW CONFIRMATORY CLAIM\n",
 contract$eligible_calibration_participants,src$nki_n,
 primary$uncalibrated$brier,primary$transbig_intercept_calibrated$brier,
 primary$uncalibrated$auc,primary$transbig_intercept_calibrated$auc,
 sensitivity$uncalibrated$brier,sensitivity$transbig_intercept_calibrated$brier))
