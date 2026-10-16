#!/usr/bin/env Rscript
# Evaluate a REAL PMx frozen JS model on a held-out publisher-defined TCGA
# partition. The previous JS trainer NEVER reads holdout-truth.csv.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript scripts/evaluate_tcga_native_pmx.R exported_dir")
if(!requireNamespace("jsonlite",quietly=TRUE))stop("Install jsonlite")
source(file.path("multiomics-engine","external_validation.R"))
dir<-args[[1L]]
p<-read.csv(file.path(dir,"frozen-pmx-predictions.csv"),stringsAsFactors=FALSE)
truth<-read.csv(file.path(dir,"holdout-truth.csv"),stringsAsFactors=FALSE)
model<-jsonlite::read_json(file.path(dir,"frozen-pmx-model.json"),simplifyVector=TRUE)
src<-jsonlite::read_json(file.path(dir,"source-contract.json"),simplifyVector=TRUE)
if(nrow(p)!=70L||nrow(truth)!=70L)stop("All original heldout people must be represented.")
if(anyNA(p$subject_id)||anyDuplicated(p$subject_id) ||
  anyNA(truth$subject_id)||anyDuplicated(truth$subject_id))
  stop("Heldout patient identity is incomplete or duplicated.")
if(!setequal(p$subject_id,truth$subject_id))
  stop("Frozen predictions do not exactly match the heldout source subjects.")
if(anyNA(p$prediction)||any(!is.finite(p$prediction)) ||
  any(p$prediction<0|p$prediction>1))
  stop("PMx model generated invalid frozen predictions.")
if(length(model$trainingSubjectIds)!=105L ||
  any(model$trainingSubjectIds %in% p$subject_id))
  stop("Training/evaluation subject partitions overlap or size changed.")
eligible<-truth$outcome %in% c("Her2","LumA")
if(sum(eligible)!=49L||length(unique(truth$outcome[eligible]))!=2L)
  stop("Binary target eligibility has changed.")
evaluation<-merge(truth[eligible,,drop=FALSE],p,by="subject_id",sort=FALSE)
evaluation$outcome<-as.integer(evaluation$outcome=="Her2")
if(nrow(evaluation)!=49L)stop("Lost a heldout eligible patient on ID merge.")
res<-validate_external_predictions(evaluation,
  outcome_type="binary",prediction_kind="probability",
  independent_cohort=TRUE,
  cohort_label="Publisher-defined TCGA within-source holdout",
  bootstrap_repetitions=300L,seed=20261009L)
if(!identical(res$status,"external_validation") ||
   !identical(res$evaluation_status,"ok") ||
   res$cohort_coverage$excluded_rows!=0L)
  stop("The original PMx scoring failed frozen heldout evaluator.")
if(!is.finite(res$metrics$auc)||!is.finite(res$metrics$brier) ||
   !is.finite(res$metrics$log_loss))
  stop("Heldout discrimination/calibration metrics must be finite.")
# Predeclared negative reference using 200 shuffled outcomes: the actual
# frozen predictions and feature coefficients are NEVER changed.
set.seed(20261009L)
null_auc<-replicate(200L,
  validation_binary_auc(sample(evaluation$outcome),evaluation$prediction))
if(any(!is.finite(null_auc))||abs(mean(null_auc)-.5)>.08)
  stop("Shuffled holdout scores do not resemble the null as expected.")

# Independent check of PMx NATIVE train-label permutations. The JS training
# program had no access to this R file's holdout class labels.
null_scores<-read.csv(file.path(dir,"pmx-train-label-null-predictions.csv"),
  stringsAsFactors=FALSE,check.names=FALSE)
if(nrow(null_scores)!=70L||anyDuplicated(null_scores$subject_id)||
   !setequal(null_scores$subject_id,p$subject_id))
  stop("Native training-label null predictions have incorrect subject coverage.")
null_columns<-grep("^perm_[0-9]+$",names(null_scores),value=TRUE)
if(length(null_columns)!=30L)
  stop("Native training-label null controls must contain 30 refitted PMx models.")
null_matched<-merge(evaluation[,c("subject_id","outcome")],null_scores,
  by="subject_id",sort=FALSE)
if(nrow(null_matched)!=49L)stop("Null controls lost an eligible test person.")
train_null_auc<-vapply(null_columns,function(col){
  prob<-null_matched[[col]]
  if(any(!is.finite(prob))||any(prob<0|prob>1))
    stop("Native PMx training-label null model generated invalid probabilities.")
  validation_binary_auc(null_matched$outcome,prob)
},numeric(1))
if(any(!is.finite(train_null_auc)) || abs(mean(train_null_auc)-.5)>.20)
  stop("Native PMx training-label permutations show unexpected residual discrimination.")

report<-list(
  benchmark="pmx_native_frozen_predictor_source_heldout_v1",
  publisherTrain=src$fullTrain,publisherTest=src$fullTest,
  pmxTrain=length(model$trainingSubjectIds),
  pmxScoredTest=nrow(p),binaryEvaluatedTest=nrow(evaluation),
  modelName=model$engine,
  featuresUsed=if(is.data.frame(model$requiredFeatures)) nrow(model$requiredFeatures) else length(model$requiredFeatures),
  selectedLambda=model$lambda,
  fittedOnHoldout=FALSE,
  sourcePreprocessingBeforeSplit=src$sourcePreprocessing,
  heldoutProteomicsAvailable=FALSE,
  sourceSameCohort=TRUE,
  trainingInputSha256=model$trainingInputSha256,
  metrics=res$metrics,
  bootstrap=res$bootstrap,
  heldoutPermutation=list(replicates=length(null_auc),meanAuc=mean(null_auc)),
  nativeTrainingLabelPermutation=list(replicates=length(train_null_auc),meanAuc=mean(train_null_auc)),
  guardrails=res$cohort_coverage,
  scientificCertification=FALSE,
  verdict="software-pipeline-heldout-test_only_not_clinically_validated",
  caveats=c(
    "Real PMx feature selection, scaling and ridgeGLM are trained on 105 subjects, not refitted on test.",
    "This is an mRNA-only test; 70 source-defined heldout patients include 49 evaluable Her2/LumA cases.",
    "Source input features were already normalized and preselected upstream by mixOmics; possible upstream leakage is not excluded.",
    "Both partitions originate in one TCGA cohort, not a different hospital, site, lab or time.",
    "AUC may be artificially perfect for this educational reference subset.",
    "This is not evidence of MOFA2/DIABLO integration, external deployment or a clinically validated biomarker."
  ))
jsonlite::write_json(report,file.path(dir,"pmx-native-heldout-report.json"),
  auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
cat(sprintf("REAL PMx MODEL HOLDOUT PASS: trained=%d scored=%d evaluated=%d AUC=%.4f Brier=%.4f permAUC=%.4f\n",
  report$pmxTrain,report$pmxScoredTest,report$binaryEvaluatedTest,
  res$metrics$auc,res$metrics$brier,mean(null_auc)))
cat(sprintf("NATIVE PMx REFIT LABEL-NULL PASS: 30 refits mean AUC=%.4f\n",mean(train_null_auc)))
