#!/usr/bin/env Rscript
# Blind native PMx predictions were generated and frozen before this
# step reads a single heldout label.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript scripts/evaluate_all_unselected_pmx.R output_dir")
dest<-args[[1L]]
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required")
source(file.path("multiomics-engine","external_validation.R"))
src<-jsonlite::read_json(file.path(dest,"source-contract.json"),simplifyVector=TRUE)
model<-jsonlite::read_json(file.path(dest,"frozen-pmx-model.json"),simplifyVector=TRUE)
pred<-utils::read.csv(file.path(dest,"frozen-pmx-predictions.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
truth<-utils::read.csv(file.path(dest,"holdout-truth.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(pred)!=src$heldout_n||nrow(truth)!=src$heldout_n||
   src$measured_features!=12625L||src$removed_features!=0L)
  stop("Wrong heldout cohort or measured probe denominator.")
if(anyDuplicated(pred$subject_id)||anyDuplicated(truth$subject_id)||
   !setequal(pred$subject_id,truth$subject_id))
  stop("Predictions and eligible holdout patient IDs differ.")
if(any(model$trainingSubjectIds %in% truth$subject_id)||
   length(model$trainingSubjectIds)!=src$train_n)
  stop("Training patient reused as heldout.")
if(anyNA(pred$prediction)||any(!is.finite(pred$prediction))||
   any(pred$prediction<0|pred$prediction>1))
  stop("Some subjects have no valid PMx heldout probability.")
if(!identical(sort(unique(truth$outcome)),sort(c("BCR/ABL","NEG"))))
  stop("Unexpected phenotypic truth group.")
joined<-merge(truth,pred,by="subject_id",sort=FALSE)
if(nrow(joined)!=nrow(truth))stop("Holdout evaluation dropped subjects.")
joined$outcome<-as.integer(joined$outcome=="BCR/ABL")
validation<-validate_external_predictions(joined,
  outcome_type="binary",prediction_kind="probability",
  independent_cohort=FALSE,
  cohort_label="ALL internal same-dataset stratified holdout",
  bootstrap_repetitions=500L,seed=20261009L)
if(!identical(validation$status,"independent_status_not_asserted") ||
   !identical(validation$evaluation_status,"ok")||
   !isTRUE(validation$cohort_coverage$complete)||
   !is.finite(validation$metrics$auc)||
   !is.finite(validation$metrics$brier))
  stop("Heldout assessment was incomplete or non-estimable.")
# No requirement on positive predictive performance; low accuracy remains
# a valid *negative result* and must be archived, not discarded.
null<-utils::read.csv(file.path(dest,"pmx-train-label-null-predictions.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
null_cols<-grep("^perm_[0-9]+$",names(null),value=TRUE)
if(length(null_cols)!=10L||nrow(null)!=nrow(truth)||
   anyDuplicated(null$subject_id)||!setequal(null$subject_id,truth$subject_id))
  stop("Not all predeclared training-label permutation models scored the heldout cohort.")
z<-merge(joined[,c("subject_id","outcome")],null,
  by="subject_id",sort=FALSE)
train_null_auc<-vapply(null_cols,function(col){
  p<-z[[col]]
  if(anyNA(p)||any(!is.finite(p))||any(p<0|p>1))
    stop("Native PMx permutation predictions contain invalid values.")
  validation_binary_auc(z$outcome,p)
},numeric(1))
if(any(!is.finite(train_null_auc)))
  stop("Some high-dimensional negative model AUC was not estimable.")
set.seed(20261009L)
test_null_auc<-replicate(300L,
  validation_binary_auc(sample(joined$outcome),joined$prediction))
if(any(!is.finite(test_null_auc)))
  stop("Holdout label permutation has non-estimable AUC.")
report<-list(
  study="Bioconductor ALL high-dimensional native PMx benchmark",
  source="bioconductor-source/ALL immutable GitHub SHA in workflow",
  train_n=src$train_n,heldout_n=src$heldout_n,
  train_class_counts=src$training_class_counts,
  heldout_class_counts=src$heldout_class_counts,
  features_measured=src$measured_features,
  features_global_prefiltered=src$removed_features,
  model_features=if(is.data.frame(model$requiredFeatures))
    nrow(model$requiredFeatures) else length(model$requiredFeatures),
  training_input_sha256=model$trainingInputSha256,
  selected_lambda=model$lambda,
  full_heldout_coverage=validation$cohort_coverage,
  auc=validation$metrics$auc,brier=validation$metrics$brier,
  log_loss=validation$metrics$log_loss,metrics=validation$metrics,
  bootstrap=validation$bootstrap,
  training_label_permutation_auc=list(n=length(train_null_auc),
    mean=mean(train_null_auc),values=as.list(train_null_auc)),
  heldout_label_permutation_auc=list(n=length(test_null_auc),
    mean=mean(test_null_auc)),
  status="same_source_holdout_no_external_validation_claim",
  method="Native PMx ridgeGLM; feature selection and tuning in training only",
  limitations=c(
    "This is a fixed-seed within-source holdout; not an independent institution, assay platform or prospective validation.",
    "12,625 measured probes were included without significance-based prior feature selection, but source microarray preprocessing was already performed before splitting.",
    "The clinical phenotypes BCR/ABL and NEG can be associated with other disease and treatment characteristics; classification does not establish causality.",
    "The small holdout cohort and 10 null refits give imprecise estimates; negative results are retained.",
    "A successful model benchmark does not validate transcriptomic count models, multi-omics or clinical decisions.",
    "The model is an experimental binary ridge predictor, not an approved medical device."
  ),
  publicationReady=FALSE
)
jsonlite::write_json(report,file.path(dest,"all-native-heldout-report.json"),
  auto_unbox=TRUE,pretty=TRUE,na="null",null="null")
cat(sprintf("ALL PMx NATIVE HOLDOUT PASS: train=%d heldout=%d probes=%d AUC=%.4f Brier=%.4f trainNullAUC=%.4f labelNullAUC=%.4f\n",
  src$train_n,src$heldout_n,src$measured_features,
  validation$metrics$auc,validation$metrics$brier,
  mean(train_null_auc),mean(test_null_auc)))
