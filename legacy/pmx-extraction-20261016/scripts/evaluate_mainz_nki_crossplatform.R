#!/usr/bin/env Rscript
# Evaluate both PREDECLARED cross-platform PMx models ONLY after JS freezes
# all probabilities and ten fully permuted MAINZ native models per arm.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)
 stop("Usage: Rscript scripts/evaluate_mainz_nki_crossplatform.R outdir")
dir<-args[[1L]]
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required.")
source(file.path("multiomics-engine","external_validation.R"))
src<-jsonlite::read_json(file.path(dir,"source-contract.json"),simplifyVector=TRUE)
truth<-utils::read.csv(file.path(dir,"external-truth.csv"),stringsAsFactors=FALSE)
if(nrow(truth)!=src$nki_n||anyDuplicated(truth$subject_id)||
   anyNA(truth$subject_id)||!setequal(unique(truth$outcome),c("ER_positive","ER_negative")))
 stop("NKI clinical labels are unavailable or duplicate.")
if(src$common_gene_count<5000L||
   src$raw_sample_id_overlap!=0L||src$suspicious_profile_pairs!=0L)
 stop("Cross-platform annotation/identity audit failed.")
truth$outcome<-as.integer(truth$outcome=="ER_positive")
train_counts<-unlist(src$mainz_er_counts)
training_prevalence<-as.numeric(train_counts[["ER_positive"]])/src$mainz_n
if(!is.finite(training_prevalence)||training_prevalence<=0||training_prevalence>=1)
 stop("Invalid frozen training prevalence.")
evaluate_arm<-function(arm){
 path<-file.path(dir,arm)
 p<-utils::read.csv(file.path(path,"external-predictions.csv"),stringsAsFactors=FALSE)
 model<-jsonlite::read_json(file.path(path,"frozen-pmx-model.json"),simplifyVector=TRUE)
 if(nrow(p)!=src$nki_n||anyDuplicated(p$subject_id)||
    !setequal(p$subject_id,truth$subject_id)||
    anyNA(p$prediction)||any(!is.finite(p$prediction))||
    any(p$prediction<0|p$prediction>1))
  stop(paste("Invalid frozen NKI probabilities for",arm))
 if(length(model$trainingSubjectIds)!=src$mainz_n||
    any(model$trainingSubjectIds %in% truth$subject_id))
  stop("External patient has been used for PMx training or sample audit failed.")
 if(model$fullGenePanel!=src$common_gene_count||
    model$fullTrainingSampleCount!=src$mainz_n)
  stop("Frozen model gene panel differs from the source mapping.")
 eval<-merge(truth,p,by="subject_id",sort=FALSE)
 if(nrow(eval)!=src$nki_n)stop("Lost a clinical outcome during patient matching.")
 v<-validate_external_predictions(eval,
   outcome_type="binary",prediction_kind="probability",
   independent_cohort=TRUE,
   cohort_label=paste0("NKI Agilent two-channel (external platform): ",arm),
   bootstrap_repetitions=500L,seed=20261009L)
 if(!identical(v$status,"external_validation")||
    !isTRUE(v$cohort_coverage$complete)||
    !is.finite(v$metrics$auc)||!is.finite(v$metrics$brier)||
    !is.finite(v$metrics$log_loss))
  stop(paste("NKI external validation metrics incomplete in",arm))
 null<-utils::read.csv(file.path(path,"native-train-null-predictions.csv"),
   stringsAsFactors=FALSE,check.names=FALSE)
 cols<-grep("^perm_[0-9]+$",names(null),value=TRUE)
 if(length(cols)!=10L||nrow(null)!=nrow(truth)||
    anyDuplicated(null$subject_id)||!setequal(null$subject_id,truth$subject_id))
  stop("Full PMx training-label null refit is missing.")
 joined<-merge(truth,null,by="subject_id",sort=FALSE)
 neg_auc<-vapply(cols,function(col){
   score<-joined[[col]]
   if(anyNA(score)||any(!is.finite(score))||any(score<0|score>1))
    stop("PMx negative model produced invalid scores.")
   validation_binary_auc(joined$outcome,score)
 },numeric(1))
 if(any(!is.finite(neg_auc)))stop("Invalid negative refit AUC")
 set.seed(20261009L)
 test_perm<-replicate(300L,
    validation_binary_auc(sample(eval$outcome),eval$prediction))
 if(any(!is.finite(test_perm)))stop("Invalid shuffled NKI clinical outcomes.")
 baseline_brier<-mean((eval$outcome-training_prevalence)^2)
 brier_skill<-1-v$metrics$brier/baseline_brier
 list(
   metrics=v$metrics,bootstrap=v$bootstrap,full_coverage=v$cohort_coverage,
   n=src$nki_n,selected_lambda=model$lambda,
   features_selected=if(is.data.frame(model$requiredFeatures))
      nrow(model$requiredFeatures) else length(model$requiredFeatures),
   observed_er_positive_fraction=mean(eval$outcome),
   mainz_training_er_positive_fraction=training_prevalence,
   mean_predicted_er_positive_fraction=mean(eval$prediction),
   training_prevalence_only_brier=baseline_brier,
   brier_skill_vs_frozen_training_prevalence=brier_skill,
   native_training_label_null=list(n=10L,mean_auc=mean(neg_auc),values=as.list(neg_auc)),
   external_label_null=list(n=300L,mean_auc=mean(test_perm)),
   frozen_training_sha256=model$trainingJsonSha256,
   independent_cohort=TRUE,independent_platform=TRUE,
   interpretation="Research only. Different published cohorts and different assay technologies; no clinical validation."
 )
}
primary<-evaluate_arm("primary_rank")
sensitivity<-evaluate_arm("sensitivity_raw")
p<-read.csv(file.path(dir,"primary_rank","external-predictions.csv"))
q<-read.csv(file.path(dir,"sensitivity_raw","external-predictions.csv"))
paired<-merge(merge(truth,p,by="subject_id"),q,by="subject_id",
   suffixes=c("_rank","_raw"))
if(nrow(paired)!=src$nki_n)stop("Paired external model comparison lost people.")
set.seed(20261010L)
pairs<-replicate(500L,{
 idx<-sample(seq_len(nrow(paired)),replace=TRUE)
 yy<-paired$outcome[idx]
 if(length(unique(yy))<2L)return(c(NA_real_,NA_real_))
 c(
  validation_binary_auc(yy,paired$prediction_rank[idx])-
   validation_binary_auc(yy,paired$prediction_raw[idx]),
  mean((yy-paired$prediction_rank[idx])^2)-
   mean((yy-paired$prediction_raw[idx])^2)
 )
})
paired_ci<-function(x){
 x<-x[is.finite(x)]
 if(length(x)<400L)stop("Too many invalid paired bootstrap resamples.")
 unname(stats::quantile(x,c(.025,.5,.975)))
}
report<-list(
 study="MAINZ Affymetrix HG-U133A to NKI Agilent/Rosetta independent-study external assay transfer",
 source_manifest=src,
 primary_analysis_predeclared="Per-patient within-mapped-gene percentile ranks, without fitting on NKI",
 sensitivity_predeclared="Original published expression values, no cross-platform batch normalization",
 primary_rank=primary,sensitivity_raw=sensitivity,
 paired_bootstrap_delta_rank_minus_raw=list(n=500L,
  auc_ci_95=as.list(paired_ci(pairs[1L,])),
  brier_ci_95=as.list(paired_ci(pairs[2L,]))),
 all_performance_results_retained=TRUE,
 clinical_certification=FALSE,
 scientific_limitations=c(
  "NKI and MAINZ are independent published cohorts measured using different assay platforms, not an interlaboratory prospective clinical validation.",
  "Entrez gene mapping excludes ambiguous IDs; probes are collapsed by median without NKI labels. Source annotation versions and gene assignments may have errors.",
  "Mainz measures absolute single-channel intensity while NKI uses dual-channel log-ratios. Rank-based comparison reduces numerical scale differences but changes the biological meaning of predictors.",
  "External preprocessing is already present upstream; this experiment does not validate the preprocessing from raw CEL, TIFF or IDAT files.",
  "Only ER receptor status is predicted; ER status is not a survival prediction or clinical decision.",
  "Both the primary rank transfer and raw sensitivity are always reported. No algorithm or calibration is chosen after looking at NKI labels.",
  "Sample identity can never be fully authenticated from anonymous published data; extreme gene-rank profile checks cannot rule out all duplicate patients.",
  "Ten full native training-label permutations are software negative controls, not a calibrated hypothesis test or guaranteed protection against dataset shifts.",
  "Successful cross-platform classification does not establish clinical validity, multivariate omics integration or robustness to all laboratories."
 )
)
jsonlite::write_json(report,file.path(dir,"mainz-nki-crossplatform-pmx-report.json"),
 auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
cat(sprintf("PMx MAINZ->NKI CROSS-PLATFORM PASS: genes=%d train=%d external=%d; rank AUC=%.4f Brier=%.4f nullAUC=%.4f; raw AUC=%.4f Brier=%.4f nullAUC=%.4f; rank Brier skill=%.4f\n",
 src$common_gene_count,src$mainz_n,src$nki_n,
 primary$metrics$auc,primary$metrics$brier,
 primary$native_training_label_null$mean_auc,
 sensitivity$metrics$auc,sensitivity$metrics$brier,
 sensitivity$native_training_label_null$mean_auc,
 primary$brier_skill_vs_frozen_training_prevalence))
