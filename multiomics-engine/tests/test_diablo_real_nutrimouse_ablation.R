#!/usr/bin/env Rscript
# REAL Nutrimouse ablation: is two-omic DIABLO better than RNA or lipid alone?
# This test uses the EXACT same externally untouched animals, seeds, labels,
# and fold-local recipes as the existing real-data DIABLO benchmark.
# It invokes real mixOmics::splsda (not a PMx native surrogate) on each
# individual view, 1 component / 10 features with NO tuning.
#
# Critical: two folds reuse some animals; fold observations are NOT
# independent. No bootstrap CI or superiority p-value is claimed.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
 stop("Usage: Rscript multiomics-engine/tests/test_diablo_real_nutrimouse_ablation.R external/mixOmics benchmark_dir")
root<-args[[1L]]
out<-args[[2L]]
if(!requireNamespace("mixOmics",quietly=TRUE)||
   !requireNamespace("jsonlite",quietly=TRUE))
 stop("Actual mixOmics and jsonlite required. No substitute allowed.")
source("multiomics-engine/advanced_methods.R")
env<-new.env(parent=emptyenv())
load(file.path(root,"data","nutrimouse.rda"),envir=env)
if(!exists("nutrimouse",envir=env,inherits=FALSE))
 stop("Pinned Nutrimouse study absent")
d<-env$nutrimouse
if(!is.list(d)||!all(c("gene","lipid","genotype","diet")%in%names(d)))
 stop("Pinned Nutrimouse data schema changed")
rna<-as.matrix(d$gene)
lipid<-as.matrix(d$lipid)
storage.mode(rna)<-"double"
storage.mode(lipid)<-"double"
n<-nrow(rna)
if(n!=40L||nrow(lipid)!=n||ncol(rna)!=120L||ncol(lipid)!=21L||
   any(!is.finite(rna))||any(!is.finite(lipid))||
   length(d$genotype)!=n||length(d$diet)!=n)
 stop("Unexpected experimental assay dimensions or missing observations")
if(is.null(colnames(rna))||anyDuplicated(colnames(rna))||
   is.null(colnames(lipid))||anyDuplicated(colnames(lipid)))
 stop("Feature IDs are ambiguous")
ids<-sprintf("NUTRIMOUSE_%03d",seq_len(n))
rownames(rna)<-ids; rownames(lipid)<-ids
y<-setNames(as.character(d$genotype),ids)
diet<-setNames(as.character(d$diet),ids)
if(anyNA(y)||length(unique(y))!=2L||min(table(y))<12L)
 stop("Public genotype classification contract has changed")
blocks<-list(transcriptomics=rna,metabolomics=lipid)
seeds<-c(20261011L,20261012L)
observed_file<-file.path(out,"observed","diablo_outer_holdout.csv")
if(!file.exists(observed_file))
 stop("First execute original paired DIABLO benchmark on identical source and seeds")
multi<-utils::read.csv(observed_file,stringsAsFactors=FALSE)
expected_n<-sum(vapply(seeds,function(seed){
 split<-prepare_diablo_holdout_split(blocks,y,seed=seed)
 length(split$testIds)
},integer(1)))
if(nrow(multi)!=expected_n||
   anyDuplicated(paste(multi$seed,multi$sampleId))||
   !identical(sort(unique(multi$seed)),sort(seeds))||
   anyNA(multi$actual)||anyNA(multi$predicted)||
   any(!multi$sampleId%in%ids)||
   any(multi$actual!=unname(y[multi$sampleId]))||
   any(!multi$predicted%in%unique(y)))
 stop("Existing DIABLO predictions do not match the complete public holdout contract")
all_arms<-list()
summaries<-list()
for(view in names(blocks)){
 one_view_rows<-list()
 for(seed in seeds){
   # EXACT same holdout ids and TRAIN-only variance filtering / imputation /
   # centering/scaling fitted by the canonical DIABLO reference engine.
   split<-prepare_diablo_holdout_split(blocks,y,seed=seed)
   train<-split$train[[view]]
   test<-split$test[[view]]
   k<-min(10L,ncol(train))
   stopifnot(identical(rownames(train),split$trainIds),
             identical(rownames(test),split$testIds),
             length(intersect(split$trainIds,split$testIds))==0L,
             identical(colnames(train),colnames(test)),
             identical(split$recipes[[view]]$fittedOn,"training_subjects_only"))
   set.seed(as.integer(seed))
   model<-mixOmics::splsda(
     X=train,Y=split$y_train,ncomp=1L,keepX=k,scale=FALSE
   )
   pred<-predict(model,newdata=test,dist="max.dist")
   vote<-pred$class[["max.dist"]]
   if(is.null(vote)||nrow(as.matrix(vote))!=length(split$testIds))
     stop(paste("Actual one-view mixOmics predictions missing:",view,seed))
   classes<-as.character(as.matrix(vote)[,1L])
   truth<-as.character(split$y_test)
   if(length(classes)!=length(truth)||anyNA(classes)||
      any(!classes%in%levels(split$y_test)))
     stop(paste("One-omic classifier produced invalid categories:",view,seed))
   one_view_rows[[length(one_view_rows)+1L]]<-data.frame(
     seed=as.integer(seed),sampleId=split$testIds,
     actual=truth,predicted=classes,view=view,
     correct=as.integer(truth==classes),stringsAsFactors=FALSE
   )
 }
 results<-do.call(rbind,one_view_rows)
 if(nrow(results)!=expected_n ||
    anyDuplicated(paste(results$seed,results$sampleId)))
   stop(paste("One-omic holdout predictions missing or duplicated:",view))
 actual_keys<-paste(multi$seed,multi$sampleId)
 view_keys<-paste(results$seed,results$sampleId)
 if(!setequal(actual_keys,view_keys))
   stop(paste("Single omics comparison did not use exactly identical outer folds:",view))
 ix<-match(actual_keys,view_keys)
 if(anyNA(ix)||any(multi$actual!=results$actual[ix]))
   stop("Outcome mismatch between DIABLO and ablated view")
 fold_metrics<-lapply(seeds,function(seed){
   z<-results[results$seed==seed,,drop=FALSE]
   recall<-vapply(sort(unique(y)),function(label)
     mean(z$predicted[z$actual==label]==label),numeric(1))
   list(seed=seed,test_n=nrow(z),BER=mean(1-recall),
        accuracy=mean(z$correct),
        class_recall=as.list(recall))
 })
 ber<-mean(vapply(fold_metrics,function(x)x$BER,numeric(1)))
 summaries[[view]]<-list(method="actual_mixOmics_sPLSDA",
    view=view,features_considered=ncol(blocks[[view]]),
    keepX=10L,ncomp=1L,scale=FALSE,distance="max.dist",
    full_training_only_preprocessing=TRUE,
    fold_metrics=fold_metrics,mean_BER=ber)
 all_arms[[view]]<-results
 utils::write.csv(results,file.path(out,paste0("single_",view,"_heldout.csv")),
    row.names=FALSE)
 cat(sprintf("REAL NUTRIMOUSE SINGLE-OMICS PASS: view=%s heldout=%d BER=%.4f\n",
   view,nrow(results),ber))
}
single<-do.call(rbind,all_arms)
multi$method<-"diablo_paired"
single$method<-paste0("splsda_",single$view)
single$view<-NULL
full<-rbind(
 multi[,c("seed","sampleId","actual","predicted","correct","method")],
 single[,c("seed","sampleId","actual","predicted","correct","method")]
)
utils::write.csv(full,file.path(out,"paired_omics_ablation_predictions.csv"),row.names=FALSE)
# Compare independent animal IDs within each fold; the signs and numerical
# differences are descriptive, not p-values from duplicated fold observations.
paired<-data.frame(seed=multi$seed,sampleId=multi$sampleId,
 actual=multi$actual,
 diablo_correct=as.integer(multi$actual==multi$predicted))
for(view in names(all_arms)){
 q<-all_arms[[view]]
 key<-paste(q$seed,q$sampleId)
 ix<-match(paste(paired$seed,paired$sampleId),key)
 stopifnot(!anyNA(ix),all(paired$actual==q$actual[ix]))
 paired[[paste0(view,"_correct")]]<-q$correct[ix]
 paired[[paste0("diablo_rescue_vs_",view)]]<-
   as.integer(paired$diablo_correct==1L & q$correct[ix]==0L)
 paired[[paste0("diablo_harm_vs_",view)]]<-
   as.integer(paired$diablo_correct==0L & q$correct[ix]==1L)
}
utils::write.csv(paired,
 file.path(out,"paired_correctness_animal_and_seed.csv"),row.names=FALSE)
multi_metrics<-lapply(seeds,function(seed){
 z<-multi[multi$seed==seed,,drop=FALSE]
 rec<-vapply(sort(unique(y)),function(label)
   mean(z$predicted[z$actual==label]==label),numeric(1))
 list(seed=seed,n_test=nrow(z),BER=mean(1-rec))
})
multi_ber<-mean(vapply(multi_metrics,function(x)x$BER,numeric(1)))
comparisons<-lapply(names(summaries),function(view){
 mono<-summaries[[view]]
 list(view=view,mean_BER_diablo=multi_ber,
  mean_BER_single_omic=mono$mean_BER,
  delta_BER_diablo_minus_single=multi_ber-mono$mean_BER,
  paired_correctness_rescued=sum(paired[[paste0("diablo_rescue_vs_",view)]]),
  paired_correctness_harmed=sum(paired[[paste0("diablo_harm_vs_",view)]]),
  descriptive_only=TRUE)
})
names(comparisons)<-names(summaries)
audit<-list(
 dataset="mixOmics Nutrimouse, real paired transcriptome and lipids",
 source="mixOmics-org/mixOmics@ef3e760526623d9e91e945e4b50be48d764efc40",
 actual_software="mixOmics::block.splsda vs mixOmics::splsda",
 reference_version=as.character(utils::packageVersion("mixOmics")),
 study_scope="Research-only within-cohort comparison; animal study, not external clinical validation",
 n_biological_animals=n,
 n_split_seeds=length(seeds),
 repeated_holdout_predictions=expected_n,
 unique_holdout_animals=length(unique(multi$sampleId)),
 same_samples_and_fold_local_recipes_verified=TRUE,
 fixed_parameters=list(keepX=10L,ncomp=1L,distance="max.dist",
    tuning_performed=FALSE,train_only_feature_imputation_and_scaling=TRUE),
 diablo_reference=list(BER=multi_ber,fold_metrics=multi_metrics),
 single_omic_references=summaries,
 diablo_vs_single_view=comparisons,
 common_subject_pairing_is_required_but_synergy_not_assumed=TRUE,
 independent_repeat_warning="Different splits can reuse mice; do not compute a nominal superiority CI or p-value treating repeated predictions as independent.",
 comparison_caveat="DIABLO block WeightedVote[max.dist] versus single-view sPLS-DA class[max.dist] are default reference prediction summaries and not necessarily identical decision rules.",
 no_calibrated_probabilities_or_human_generalization=TRUE,
 significant_benefit_claimed=FALSE,
 benchmark_success_means_execution_and_complete_results_not_accuracy=TRUE,
 caveats=c(
  "No new independent laboratory or external experimental cohort",
  "Previously inspected Nutrimouse genotype outcomes, retrospective method comparison",
  "Only two prespecified splits across 40 animals, finite small precision",
  "No tuning and no comparison to more expressive or nonlinear reference methods",
  "The 120 transcriptomic genes were PRESELECTED from approximately 30,000 candidates by original Nutrimouse study; therefore this benchmark is not a whole-transcriptome unbiased-feature test",
  "Sample selection and published source preprocessing were not independently audited",
  "Genotype and diet may associate; DIABLO and single-view fits do not adjust genotype effects for diet",
  "Difference in balanced error alone is NOT a test of incremental multi-omic predictive utility"
 )
)
jsonlite::write_json(audit,
 file.path(out,"real_nutrimouse_diablo_vs_single_omics.json"),
 pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
cat(sprintf("REAL NUTRIMOUSE SAME-FOLD INCREMENTAL VALUE PASS: DIABLO BER %.4f; transcriptome only %.4f; lipids only %.4f; comparisons descriptive, NO multiomics synergy claim\n",
 multi_ber,summaries$transcriptomics$mean_BER,summaries$metabolomics$mean_BER))
