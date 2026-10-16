#!/usr/bin/env Rscript
# Real biological cohort, independent from the Nutrimouse animal experiment.
# TCGA breast.TCGA from pinned mixOmics commit. This is an INTERNAL held-out
# validation on publicly processed and preselected samples, NOT an independent
# prospective cohort or a blinded multicentre generalization test.
#
# Goal: evaluate actual mixOmics multiblock model alongside the two mono-omic
# baselines on IDENTICAL held-out subjects, with a prespecified negative control.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop("Usage: Rscript tcga_two_view_research_validation.R MIXOMICS_ROOT OUTPUT_DIR")
root <- args[[1L]]; out <- args[[2L]]
source("multiomics-engine/advanced_methods.R")
if(!requireNamespace("mixOmics",quietly=TRUE) ||
   !requireNamespace("jsonlite",quietly=TRUE))
  stop("Genuine mixOmics and jsonlite required. No browser/ML surrogate.")
file <- file.path(root,"data","breast.TCGA.rda")
if(!file.exists(file))stop("Pinned breast.TCGA biological source unavailable")
study <- new.env(parent=emptyenv())
load(file,envir=study)
if(!exists("breast.TCGA",envir=study,inherits=FALSE))
 stop("Invalid breast.TCGA experimental source")
data <- study$breast.TCGA$data.train
if(!is.list(data)|| !all(c("mrna","protein","subtype") %in% names(data)))
 stop("TCGA biological input schema changed")
eligible <- as.character(data$subtype)%in%c("Her2","LumA")
omics <- list(transcriptomics=as.matrix(data$mrna[eligible,,drop=FALSE]),
              proteomics=as.matrix(data$protein[eligible,,drop=FALSE]))
classes <- as.character(data$subtype[eligible])
ids <- sprintf("TCGA_SUBJECT_%03d",seq_along(classes))
if(length(classes)<32L||min(table(classes))<12L ||
  nrow(omics$transcriptomics)!=length(ids) ||
  nrow(omics$proteomics)!=length(ids) ||
  any(vapply(omics,function(x)any(!is.finite(x)),logical(1))))
 stop("Insufficient or invalid TCGA matched patients")
for(v in names(omics)){
 if(is.null(colnames(omics[[v]]))||anyDuplicated(colnames(omics[[v]])))
  stop("Ambiguous source feature IDs")
 rownames(omics[[v]])<-ids
 storage.mode(omics[[v]])<-"double"
}
y <- setNames(classes,ids)
seeds <- as.integer(c(20261013L,20261014L,20261015L,20261016L,20261017L))
# All hyperparameters and target are declared before running. 1 component
# and 10 features PER VIEW is a method-specific budget, not an equal
# total-dimensionality comparison. The control cannot certify added value.
keepX <- 10L; ncomp <- 1L; test_fraction <- 0.25
dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(file.exists(file.path(out,"predictions.csv")))
 stop("Refusing to overwrite previous biological validation evidence")
result_rows <- list();splits <- list();z<-0L
models <- c("multiomics_DIABLO","transcriptomics_sPLSDA",
            "proteomics_sPLSDA","scrambled_protein_DIABLO")
for(seed in seeds){
  split <- prepare_diablo_holdout_split(
   omics,y,seed=seed,test_fraction=test_fraction)
  if(length(intersect(split$trainIds,split$testIds)) ||
    !all(vapply(split$recipes,function(x)
      identical(x$fittedOn,"training_subjects_only"),logical(1))))
    stop("Data leakage: overlapping subjects or pre-trained preprocessing")
  setup <- matrix(0.1,nrow=2,ncol=2,
     dimnames=list(names(omics),names(omics)))
  diag(setup) <- 0
  keep <- setNames(lapply(split$train,function(x)min(keepX,ncol(x))),
    names(omics))
  set.seed(seed)
  fit <- mixOmics::block.splsda(
    X=split$train,Y=split$y_train,ncomp=ncomp,
    keepX=keep,design=setup,scale=FALSE)
  vote <- predict(fit,newdata=split$test,dist="max.dist")$WeightedVote[["max.dist"]]
  if(is.null(vote))stop("Real DIABLO weighted votes missing")
  multi <- as.character(as.matrix(vote)[,1L])
  pred <- list(multiomics_DIABLO=multi)
  for(view in names(omics)){
    k <- min(keepX,ncol(split$train[[view]]))
    set.seed(seed)
    m <- mixOmics::splsda(X=split$train[[view]],
         Y=split$y_train,ncomp=1L,keepX=k,scale=FALSE)
    vote <- predict(m,newdata=split$test[[view]],
                    dist="max.dist")$class[["max.dist"]]
    if(is.null(vote))stop("Real mono-omic class output missing")
    pred[[paste0(view,"_sPLSDA")]] <- as.character(as.matrix(vote)[,1L])
  }
  # Negative control: scramble protein rows separately in TRAIN and TEST,
  # preserve the RNA data, labels, and the original patient holdout.
  # This destroys subject correspondence for the protein view.
  negTrain<-split$train
  negTest<-split$test
  set.seed(seed+100000L)
  train_shuffle<-sample(seq_len(nrow(negTrain$proteomics)))
  test_shuffle<-sample(seq_len(nrow(negTest$proteomics)))
  negTrain$proteomics <- negTrain$proteomics[train_shuffle,,drop=FALSE]
  negTest$proteomics <- negTest$proteomics[test_shuffle,,drop=FALSE]
  rownames(negTrain$proteomics)<-split$trainIds
  rownames(negTest$proteomics)<-split$testIds
  set.seed(seed)
  neg <- mixOmics::block.splsda(
    X=negTrain,Y=split$y_train,ncomp=1L,
    keepX=keep,design=setup,scale=FALSE)
  negvote <- predict(neg,newdata=negTest,dist="max.dist")$WeightedVote[["max.dist"]]
  if(is.null(negvote))stop("Scrambled-pair control predictions missing")
  pred$scrambled_protein_DIABLO <- as.character(as.matrix(negvote)[,1L])
  for(method in models){
    classes_pred <- pred[[method]]
    truth <- as.character(split$y_test)
    if(length(classes_pred)!=length(truth)||anyNA(classes_pred)||
       any(!classes_pred%in%levels(split$y_test)))stop("Incomplete predictions")
    z<-z+1L
    result_rows[[z]] <- data.frame(
      seed=seed,subject_id=split$testIds,
      method=method,true_label=truth,predicted=classes_pred,
      correct=as.integer(truth==classes_pred),
      class_target="Her2_vs_LumA",
      training_only_preprocessing=TRUE,stringsAsFactors=FALSE)
  }
  splits[[as.character(seed)]]<-list(
    n_train=length(split$trainIds),n_test=length(split$testIds),
    train_ids=split$trainIds,test_ids=split$testIds,
    feature_counts=vapply(split$train,ncol,integer(1)),
    train_only_preprocessing=TRUE)
  cat("TCGA genuine biological matched holdout:",seed,"subjects",length(split$testIds),"\n")
}
rows <- do.call(rbind,result_rows)
stopifnot(nrow(rows)==length(seeds)*length(models)*splits[[1L]]$n_test ||
  nrow(rows)>=length(seeds)*length(models)*6L)
if(anyDuplicated(paste(rows$seed,rows$subject_id,rows$method)))
 stop("Duplicate holdout observation")
for(seed in seeds){
 part<-rows[rows$seed==seed,,drop=FALSE]
 keys<-split(part,part$method)
 ref<-keys[[models[[1L]]]]
 for(method in models){
   other<-keys[[method]]
   stopifnot(setequal(other$subject_id,ref$subject_id))
   idx<-match(ref$subject_id,other$subject_id)
   stopifnot(all(ref$true_label==other$true_label[idx]))
 }
}
utils::write.csv(rows,file.path(out,"predictions.csv"),row.names=FALSE)
folds<-do.call(rbind,lapply(split(rows,list(rows$seed,rows$method),drop=TRUE),
  function(x){
    recalls<-vapply(c("Her2","LumA"),function(g)
      mean(x$predicted[x$true_label==g]==g),numeric(1))
    data.frame(seed=x$seed[[1]],method=x$method[[1]],
      test_n=nrow(x),BER=mean(1-recalls),
      accuracy=mean(x$correct),stringsAsFactors=FALSE)
  }))
utils::write.csv(folds,file.path(out,"fold_metrics.csv"),row.names=FALSE)
averages<-aggregate(cbind(BER,accuracy)~method,folds,mean)
averages$substantive_improvement_demonstrated <- FALSE
averages$external_cohort_replicated <- FALSE
utils::write.csv(averages,file.path(out,"descriptive_summary.csv"),row.names=FALSE)
manifest <- list(
  schema="pmx-omics-tcga-matched-holdout/v1",
  dataset="public mixOmics breast.TCGA$data.train",origin="TCGA processed source",
  source_repository="mixOmics-org/mixOmics",
  source_commit="ef3e760526623d9e91e945e4b50be48d764efc40",
  source_sample_selection="Her2 and LumA from already processed data.train",
  n_subjects=length(ids),methods=models,
  number_of_repeated_splits=length(seeds),
  seeds=as.list(seeds),train_test_splits=splits,
  reference_engine="mixOmics",
  mixomics_version=as.character(utils::packageVersion("mixOmics")),
  r_version=as.character(getRversion()),
  feature_budget="10 per view, 1 component, design 0.1, no parameter tuning",
  mandatory_sensitivity="break protein subject correspondence in train and test",
  limitations=c(
   "Same cohort is reused across repeated splits; no independent external replication",
   "Source data were previously processed and may have been selected upstream",
   "The multi-omic arm uses 20 selected features in total versus 10 for each single-omic arm",
   "Model tuning is fixed, but published source data may include upstream transformations",
   "No confidence interval treating overlapping splits as independent observations",
   "No evidence of improvement unless explicitly supported by the observed results",
   "No biological biomarker validation and no clinical decision-support claim"),
  independent_external_validation=FALSE,
  confirmatory_certification=FALSE)
jsonlite::write_json(manifest,file.path(out,"manifest.json"),
                     pretty=TRUE,auto_unbox=TRUE)
print(averages,row.names=FALSE)
cat("REAL TCGA ANALYSES COMPLETE: paired subjects and honest negative controls; NOT external cohort replication\n")
