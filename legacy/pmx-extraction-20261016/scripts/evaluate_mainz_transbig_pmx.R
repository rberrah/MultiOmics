#!/usr/bin/env Rscript
# EVALUATION ONLY: after frozen native PMx predictions were created on
# study TRANSBIG without revealing a single evaluation phenotype.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript scripts/evaluate_mainz_transbig_pmx.R output_dir")
dest<-args[[1L]]
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required.")
source(file.path("multiomics-engine","external_validation.R"))
manifest<-jsonlite::read_json(file.path(dest,"source-contract.json"),simplifyVector=TRUE)
model<-jsonlite::read_json(file.path(dest,"frozen-pmx-model.json"),simplifyVector=TRUE)
prediction<-read.csv(file.path(dest,"transbig-frozen-predictions.csv"),stringsAsFactors=FALSE)
truth<-read.csv(file.path(dest,"holdout-truth.csv"),stringsAsFactors=FALSE)
if(nrow(prediction)!=manifest$external_n||nrow(truth)!=manifest$external_n ||
   anyDuplicated(prediction$subject_id)||anyDuplicated(truth$subject_id) ||
   !setequal(prediction$subject_id,truth$subject_id))
  stop("External clinical sample identities mismatch or are duplicated.")
if(length(model$trainingSubjectIds)!=manifest$train_n ||
   any(model$trainingSubjectIds%in%truth$subject_id))
  stop("The published training and external evaluation subjects overlap.")
if(manifest$expression_duplicate_audit$suspected_pairs_above_0995!=0L||
   manifest$cross_study_name_duplicate_count!=0L)
  stop("Possible cross-study duplicates have not been excluded.")
if(manifest$all_shared_probes!=22283L ||
   manifest$global_probes_removed!=0L)
  stop("Unexpected probe panel or global filtering.")
if(!setequal(unique(truth$outcome),c("ER_positive","ER_negative")))
  stop("ER receptor phenotype not interpretable.")
evaluation<-merge(truth,prediction,by="subject_id",sort=FALSE)
if(nrow(evaluation)!=manifest$external_n||any(!is.finite(evaluation$prediction)))
  stop("At least one external patient has no prediction.")
evaluation$outcome<-as.integer(evaluation$outcome=="ER_positive")
result<-validate_external_predictions(
  evaluation,outcome_type="binary",prediction_kind="probability",
  independent_cohort=TRUE,
  cohort_label="TRANSBIG (different study from MAINZ, same Affymetrix platform)",
  bootstrap_repetitions=500L,seed=20261009L
)
if(!identical(result$status,"external_validation")||
   !isTRUE(result$cohort_coverage$complete)||
   !is.finite(result$metrics$auc)||
   !is.finite(result$metrics$brier)||
   !is.finite(result$metrics$log_loss))
  stop("External evaluation has incomplete or nonfinite metrics.")

null<-read.csv(file.path(dest,"train-label-null-scores.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
cols<-grep("^perm_[0-9]+$",names(null),value=TRUE)
if(length(cols)!=10L||nrow(null)!=manifest$external_n||
   anyDuplicated(null$subject_id)||!setequal(null$subject_id,truth$subject_id))
  stop("External negative-control score panel is incomplete.")
matched<-merge(evaluation[,c("subject_id","outcome")],null,
  by="subject_id",sort=FALSE)
train_null_auc<-vapply(cols,function(key){
  p<-matched[[key]]
  if(any(!is.finite(p))||any(p<0|p>1))stop("Negative refit generated invalid probabilities.")
  validation_binary_auc(matched$outcome,p)
},numeric(1))
if(any(!is.finite(train_null_auc)))
  stop("Native PMx training-label null metrics not estimable.")
set.seed(20261009L)
test_label_null<-replicate(300L,
  validation_binary_auc(sample(evaluation$outcome),evaluation$prediction))
if(any(!is.finite(test_label_null)))
  stop("Null permutations of external outcomes failed.")
# Report negative results without filtering by AUC. Training-label null AUC
# could be nontrivially above chance in an unbalanced or shifted test cohort.
# It is a diagnostic, not a significance test or nominal false positive rate.
metrics<-result$metrics
report<-list(
  benchmark="MAINZ_TO_TRANSBIG_INDEPENDENT_STUDY_FROZEN_PMX_2026",
  original_train_study=manifest$train_study,
  independent_test_study=manifest$evaluation_study,
  train_n=manifest$train_n,external_n=manifest$external_n,
  train_class_counts=manifest$train_class_counts,
  external_class_counts=manifest$external_class_counts,
  probes_considered=manifest$all_shared_probes,
  feature_set_policy="All shared HG-U133A probes, exact IDs; PMx selects 12 from training only",
  features_selected=if(is.data.frame(model$requiredFeatures))
    nrow(model$requiredFeatures) else length(model$requiredFeatures),
  lambda=model$lambda,training_input_sha256=model$trainingInputSha256,
  independent_study_asserted=TRUE,external_platform_identical=TRUE,
  duplicate_exclusion_audit=manifest$expression_duplicate_audit,
  sample_ids_common=manifest$cross_study_name_duplicate_count,
  all_external_patients_scored=result$cohort_coverage,
  metrics=metrics,bootstrap=result$bootstrap,
  training_label_null=list(n=length(train_null_auc),mean_auc=mean(train_null_auc),
    all_aucs=as.list(train_null_auc)),
  external_label_null=list(n=length(test_label_null),mean_auc=mean(test_label_null)),
  scientific_certification=FALSE,
  interpretation="External cohort discrimination and calibration of a frozen experimental PMx ridge predictor; different studies but same assay technology",
  limitations=c(
    "Distinct clinical studies MAINZ and TRANSBIG, but the public preprocessed assays share a platform and may share undocumented batch effects.",
    "Cross-study subject uniqueness is screened by sample ID and extreme expression-profile correlation; anonymous re-arrayed duplicates cannot be definitively excluded.",
    "Both published eSets were already normalized upstream; this is not a raw-data preprocessing or cross-platform transportability validation.",
    "ER clinical phenotype labelling and receptor methods may differ between source studies.",
    "External study prevalence and enrollment criteria differ; AUC does not establish calibrated clinical risk or treatment utility.",
    "Small set of null training-label fits is a software control only; it does not establish FDR, sensitivity or specificity.",
    "One transcriptomic microarray only; no validation of multiomics, DIABLO/MOFA2, or RNA-seq count analysis.",
    "No clinical decisions can be based on this research benchmark."
  )
)
jsonlite::write_json(report,file.path(dest,"mainz-transbig-pmx-report.json"),
  pretty=TRUE,auto_unbox=TRUE,na="null",null="null")
cat(sprintf("PMx CROSS-STUDY FROZEN EXTERNAL PASS: MAINZ n=%d TRANSBIG n=%d probes=%d AUC=%.4f Brier=%.4f nullTrainAUC=%.4f nullTestAUC=%.4f\n",
  manifest$train_n,manifest$external_n,manifest$all_shared_probes,
  metrics$auc,metrics$brier,mean(train_null_auc),mean(test_label_null)))
