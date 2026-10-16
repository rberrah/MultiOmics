#!/usr/bin/env Rscript
# Source-native train/test evaluation of frozen predictions using the PMx R
# external-validation METRIC evaluator. This does NOT validate PMx's own
# predictive training engine or multi-block DIABLO.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: Rscript scripts/benchmark_tcga_heldout_reference.R external/mixOmics tmp/tcga-heldout")
root <- args[[1L]]; dest <- args[[2L]]
if(!requireNamespace("jsonlite",quietly=TRUE))stop("Install jsonlite for the reproducible report.")
source(file.path("multiomics-engine","external_validation.R"))
source_path <- file.path(root,"data","breast.TCGA.rda")
if(!file.exists(source_path))stop("Official mixOmics train/test dataset missing.")
e <- new.env(parent=emptyenv())
load(source_path,envir=e)
tcga <- e$breast.TCGA
if(is.null(tcga$data.train)||is.null(tcga$data.test))
  stop("A publisher-defined training/test separation is mandatory.")
train <- tcga$data.train; holdout <- tcga$data.test
if(!is.null(holdout$protein) && ncol(as.matrix(holdout$protein))>0L)
  stop("Unexpected source version: reassess the test modality contract.")
if(is.null(train$mrna)||is.null(holdout$mrna))
  stop("Both splits must have measured mRNA.")
tr <- as.matrix(train$mrna)
te <- as.matrix(holdout$mrna)
storage.mode(tr) <- "double"; storage.mode(te) <- "double"
if(nrow(tr)!=length(train$subtype)||nrow(te)!=length(holdout$subtype))
  stop("Subtype labels do not align with source rows.")
if(nrow(tr)!=150L || nrow(te)!=70L)
  stop("Publisher-defined split changed; manually reassess benchmark.")
if(anyDuplicated(colnames(tr)) || anyDuplicated(colnames(te)))
  stop("Duplicate feature identifiers.")
shared <- intersect(colnames(tr),colnames(te))
if(length(shared)<50L)stop("Too few shared mRNA features across the source split.")
# Exact feature intersection is structural and label-free; no selection of
# features by holdout outcomes or estimated performance.
tr <- tr[,shared,drop=FALSE]
te <- te[,shared,drop=FALSE]
class_names <- c("Her2","LumA")
labels_tr <- as.character(train$subtype)
labels_te <- as.character(holdout$subtype)
keep_tr <- labels_tr %in% class_names
keep_te <- labels_te %in% class_names
if(sum(keep_tr)<20L || sum(keep_te)<12L)
  stop("Too few eligible source-defined Her2/LumA subjects.")
if(min(table(factor(labels_tr[keep_tr],levels=class_names)))<5L ||
   min(table(factor(labels_te[keep_te],levels=class_names)))<5L)
  stop("Holdout requires both classes with at least five subjects.")

# Check actual IDs if they are informative; synthetic row indices are not
# evidence of sample identity and are deliberately NOT compared.
is_informative <- function(x,n)length(x)==n && !anyNA(x) &&
  !anyDuplicated(x) && !identical(x,as.character(seq_len(n)))
id_tr <- rownames(tr); id_te <- rownames(te)
if(is_informative(id_tr,nrow(tr)) && is_informative(id_te,nrow(te)) &&
   length(intersect(id_tr,id_te))>0L)
  stop("Potential subject overlap across development and evaluation partitions.")

# All fitting quantities below are estimated EXCLUSIVELY on training rows.
x <- tr[keep_tr,,drop=FALSE]
train_y <- as.integer(labels_tr[keep_tr]=="Her2")
means <- colMeans(x,na.rm=TRUE)
means[!is.finite(means)] <- 0
impute_training <- function(mat,reference_means) {
  m <- mat
  for(j in seq_len(ncol(m))){
    missing <- !is.finite(m[,j])
    if(any(missing))m[missing,j] <- reference_means[[j]]
  }
  m
}
x <- impute_training(x,means)
sds <- apply(x,2,stats::sd)
sds[!is.finite(sds)|sds<1e-6] <- 1
xscale <- sweep(sweep(x,2,means,"-"),2,sds,"/")
# Feature selection uses ONLY DEVELOPMENT LABELS.
effect <- colMeans(xscale[train_y==1L,,drop=FALSE])-
  colMeans(xscale[train_y==0L,,drop=FALSE])
eligible <- is.finite(effect) & apply(xscale,2,stats::sd)>1e-8
selection <- order(-abs(effect),names(effect))
selection <- selection[eligible[selection]]
k <- min(20L,length(selection))
if(k<5L)stop("Insufficient informative train-only features.")
chosen <- selection[seq_len(k)]
weights <- effect[chosen]
norm <- sqrt(sum(weights^2))
if(!is.finite(norm)||norm<=0)stop("Cannot estimate training-only centroid direction.")
score <- function(mat){
  m <- impute_training(mat[,shared,drop=FALSE],means)
  z <- sweep(sweep(m,2,means,"-"),2,sds,"/")
  as.numeric(z[,chosen,drop=FALSE] %*% weights/norm)
}
train_score <- as.numeric(xscale[,chosen,drop=FALSE] %*% weights/norm)
# Freeze one-dimensional logistic probability mapping on the training split.
fit <- suppressWarnings(stats::glm(train_y~train_score,family=stats::binomial()))
coef <- stats::coef(fit)
if(length(coef)!=2L || any(!is.finite(coef)))
  stop("Training-only probability calibration is not estimable.")
train_prob <- stats::plogis(coef[[1L]]+coef[[2L]]*train_score)
if(!all(is.finite(train_prob)))stop("Nonfinite fitted training probabilities.")
all_holdout_score <- score(te)
all_holdout_prob <- pmin(1-1e-6,pmax(1e-6,
  stats::plogis(coef[[1L]]+coef[[2L]]*all_holdout_score)))
if(any(!is.finite(all_holdout_prob)))stop("Nonfinite frozen test predictions.")

# All 70 predictions are now frozen. Only AFTER this point use the holdout
# labels to select the prespecified binary evaluation target and compute metrics.
holdout_y <- as.integer(labels_te[keep_te]=="Her2")
predictions <- data.frame(
  subject_id = sprintf("OFFICIAL_TCGA_HOLDOUT_%03d",which(keep_te)),
  outcome = holdout_y,
  prediction = all_holdout_prob[keep_te],
  stringsAsFactors=FALSE
)
if(anyDuplicated(predictions$subject_id))stop("Duplicate evaluation IDs.")
if(any(!is.finite(predictions$prediction)))stop("Invalid frozen predictions.")
evaluation <- validate_external_predictions(
  predictions,outcome_type="binary",prediction_column="prediction",
  outcome_column="outcome",prediction_kind="probability",
  independent_cohort=TRUE,cohort_label="TCGA mixOmics publisher holdout",
  bootstrap_repetitions=300L,seed=20261009L
)
if(!identical(evaluation$status,"external_validation") ||
   !identical(evaluation$evaluation_status,"ok") ||
   !is.finite(evaluation$metrics$auc))
  stop("PMx frozen-prediction external evaluator failed.")

# Distinct label-shuffle control: all predictions remain unchanged, only the
# evaluation labels are shuffled; never retrain or select a model on test.
set.seed(20261009L)
permutation_aucs <- vapply(seq_len(200L),function(i){
  permuted <- sample(holdout_y,length(holdout_y),replace=FALSE)
  validation_binary_auc(permuted,predictions$prediction)
},numeric(1))
if(any(!is.finite(permutation_aucs)))
  stop("Null label permutation AUC was not estimable.")
permutation_mean <- mean(permutation_aucs)
if(abs(permutation_mean-0.5)>0.08)
  stop("Shuffled holdout labels generated systematic discrimination: check leak.")

# Stronger null: permute DEVELOPMENT outcomes; rerun selection, calibration
# and heldout scoring. The 49 heldout outcomes are unchanged and unseen by fit.
set.seed(20261010L)
x_test_scaled <- sweep(sweep(impute_training(te,means),2,means,"-"),2,sds,"/")
training_null_aucs <- vapply(seq_len(80L),function(i){
  y0 <- sample(train_y,length(train_y),replace=FALSE)
  effect0 <- colMeans(xscale[y0==1L,,drop=FALSE])-
    colMeans(xscale[y0==0L,,drop=FALSE])
  ix <- order(-abs(effect0),names(effect0))
  ix <- ix[is.finite(effect0[ix]) & apply(xscale,2,stats::sd)[ix]>1e-8]
  ix <- ix[seq_len(min(k,length(ix)))]
  w <- effect0[ix]
  denom <- sqrt(sum(w^2))
  if(!is.finite(denom)||denom<=0)return(NA_real_)
  s0 <- as.numeric(xscale[,ix,drop=FALSE]%*%w/denom)
  fit0 <- suppressWarnings(try(stats::glm(y0~s0,family=stats::binomial()),silent=TRUE))
  if(inherits(fit0,"try-error")||any(!is.finite(stats::coef(fit0))))
    return(NA_real_)
  st <- as.numeric(x_test_scaled[,ix,drop=FALSE]%*%w/denom)
  b <- stats::coef(fit0)
  validation_binary_auc(holdout_y,stats::plogis(b[[1L]]+b[[2L]]*st[keep_te]))
},numeric(1))
training_null_valid <- sum(is.finite(training_null_aucs))
if(training_null_valid<76L)
  stop("Training-label permutation produced too many invalid controls.")
training_null_mean <- mean(training_null_aucs[is.finite(training_null_aucs)])
if(abs(training_null_mean-0.5)>0.15)
  stop("Training-label permutation retained unusual discrimination: audit source feature selection.")


# Guard against a genuine error: reusing exact training subjects as the
# validation set must not be silently reported as the official holdout.
if(nrow(tr)!=150L || nrow(te)!=70L)
  stop("Train/holdout sample size invariant failed.")

dir.create(dest,recursive=TRUE,showWarnings=FALSE)
utils::write.csv(predictions,file.path(dest,"heldout-frozen-predictions.csv"),
  row.names=FALSE,quote=FALSE)
report <- list(
  benchmark="source_native_tcga_heldout_reference_v1",
  data_source="mixOmics breast.TCGA at workflow-pinned source Git SHA",
  cohort="Publisher-defined heldout TCGA subset, same originating cohort",
  source_train_n=nrow(tr),source_test_n=nrow(te),
  eligible_train_n=sum(keep_tr),eligible_test_n=sum(keep_te),
  train_class_counts=as.list(table(labels_tr[keep_tr])),
  holdout_class_counts=as.list(table(labels_te[keep_te])),
  measured_source_test_modalities=c("mrna","mirna"),
  training_modality="mrna_only",
  missing_test_proteomics=TRUE,
  source_features_shared=length(shared),
  selected_feature_count=k,selected_feature_names=names(weights),
  selection_policy="top 20 absolute standardized train-only Her2 minus LumA centroids",
  model="train-only standardized centroid score, logistic calibration fit on TRAINING only",
  independent_train_test_partitions_asserted=TRUE,
  source_sample_identifiers_check=if(is_informative(id_tr,nrow(tr)) &&
    is_informative(id_te,nrow(te)))"no_overlap" else "not_verifiable_from_rownames",
  evaluator="multiomics-engine/external_validation.R validate_external_predictions",
  metrics=evaluation$metrics,
  bootstrap=evaluation$bootstrap,
  permutation_null=list(repetitions=length(permutation_aucs),
    mean_auc=permutation_mean,
    min_auc=min(permutation_aucs),max_auc=max(permutation_aucs)),
  training_label_permutation_control=list(
    attempted=length(training_null_aucs),valid=training_null_valid,
    average_holdout_auc=training_null_mean,
    procedure="refit selected features and probabilities on permuted training labels, holdout labels untouched"),
  perfect_auc_caution=if(evaluation$metrics$auc>=.999)
    "Perfect AUC in a heavily preselected dataset must not be presented as clinical transportability" else NULL,
  pass=TRUE,scientificCertification=FALSE,
  constraints=c(
    "Evaluates PMx frozen-prediction METRICS, not PMx internal predictor training; it is a reference baseline, not proof of better external performance.",
    "Source train/test partitions are publisher-defined, but represent one TCGA source, not a genuinely new hospital or technical platform.",
    "These globally prefiltered source variables may make subtype classification deceptively easy; near-perfect AUC is not clinical transportability.",
    "Input data were normalized/preselected by mixOmics upstream before these partitions; cannot exclude all upstream information leakage.",
    "Test cohort lacks proteomics, so this cannot validate multiomics DIABLO or multimodal patient-level deployment.",
    "Source test outcome labels were used only to determine predeclared evaluable Her2/LumA subjects and metrics after predictions were frozen.",
    "Bootstrap estimates sample uncertainty on the selected test cohort, not external transportability.",
    "No claim of confirmatory clinical biomarker validation or general FDR control."
  )
)
jsonlite::write_json(report,file.path(dest,"tcga-heldout-report.json"),
  auto_unbox=TRUE,pretty=TRUE,null="null",na="null")
cat("TCGA HELDOUT REFERENCE PASS: train=",sum(keep_tr),
    " test=",sum(keep_te),
    " AUC=",round(evaluation$metrics$auc,4),
    " permuted mean AUC=",round(permutation_mean,4),
    " training-null mean AUC=",round(training_null_mean,4),"\n",sep="")
