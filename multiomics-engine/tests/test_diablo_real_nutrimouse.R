#!/usr/bin/env Rscript
# REAL public two-omic validation of the ACTUAL mixOmics::block.splsda
# used by PMx. Nutrimouse contains transcriptome + liver fatty-acid/lipid
# profiles on the SAME 40 animals and experimental genotype/diet metadata.
#
# Primary study: genotype classification; 2 predeclared independent-subject
# holdouts, no hyperparameter tuning (fixed keepX=10 and ncomp=1).
# Comparator: 8 FULL genotype-label permutations (re-fit DIABLO each time).
# Additional falsification: deliberately BREAK the lipid-to-animal pairing
# and refit on exactly the same genotype labels and outer split seeds.
# A diet-only baseline uses TRAIN diet->genotype association, never test labels.
#
# This is a public REAL-DATA software/biology stress test, not external-cohort
# validation, NOT validation of clinical biomarkers, and not formal p-value
# control with only eight permutations. Preserve poor performance honestly.

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
  stop("Usage: Rscript multiomics-engine/tests/test_diablo_real_nutrimouse.R external/mixOmics outputdir")
root <- args[[1L]]
dest <- args[[2L]]
if(!requireNamespace("mixOmics",quietly=TRUE))
  stop("Real mixOmics required. No surrogate or invented DIABLO model.")
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("jsonlite required for machine-readable audit.")
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/validation_benchmarks.R")

e <- new.env(parent=emptyenv())
load(file.path(root,"data","nutrimouse.rda"),envir=e)
if(!exists("nutrimouse",envir=e,inherits=FALSE))
  stop("Immutable official public Nutrimouse object absent.")
data <- e$nutrimouse
if(!is.list(data)||!all(c("gene","lipid","genotype","diet")%in%names(data)))
  stop("Nutrimouse assay blocks and experimental phenotypes absent.")
gene <- as.matrix(data$gene)
lipid <- as.matrix(data$lipid)
storage.mode(gene)<-"double"
storage.mode(lipid)<-"double"
n <- nrow(gene)
if(n!=40L||nrow(lipid)!=n||ncol(gene)<50L||ncol(lipid)<10L)
  stop("Public Nutrimouse measured multi-omic dimensions changed.")
if(length(data$genotype)!=n||length(data$diet)!=n)
  stop("Nutrimouse genotype/diet phenotype lengths differ from assay rows.")
if(any(!is.finite(gene))||any(!is.finite(lipid)))
  stop("Nutrimouse missing values require an explicit missingness audit.")
if(anyDuplicated(colnames(gene))||anyDuplicated(colnames(lipid))||
   is.null(colnames(gene))||is.null(colnames(lipid)))
  stop("Real assay feature names missing/duplicated.")
y <- as.character(data$genotype)
diet <- as.character(data$diet)
if(anyNA(y)||anyNA(diet)||any(!nzchar(y))||any(!nzchar(diet))||
   length(unique(y))!=2L||min(table(y))<12L)
  stop("Observed genotype labels do not meet predeclared binary holdout contract.")
if(length(unique(diet))<3L)stop("Nutrimouse diet design missing.")
ids <- sprintf("NUTRIMOUSE_%03d",seq_len(n))
rownames(gene)<-ids
rownames(lipid)<-ids
names(y)<-ids
names(diet)<-ids
blocks <- list(transcriptomics=gene,metabolomics=lipid)
if(!identical(rownames(blocks$transcriptomics),
             rownames(blocks$metabolomics)))
  stop("Subject pairing failed before actual DIABLO training.")

# Locked before fitting. Repeated holdouts REUSE animals across seeds:
# do not treat their 20 test predictions as independent sampling units.
seeds <- c(20261011L,20261012L)
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
observed <- run_diablo_outer_holdout(
  blocks,y,file.path(dest,"observed"),seeds=seeds,tune=FALSE
)
if(!identical(observed$status,"ok") ||
   nrow(observed$predictions)!=sum(vapply(observed$summary$folds,
         function(x)x$heldOutSubjects,integer(1))) ||
   any(!is.finite(observed$summary$meanBER)))
  stop("The actual DIABLO model did not score all predeclared held-out animals.")
if(anyDuplicated(paste(observed$predictions$seed,observed$predictions$sampleId))||
   any(!observed$predictions$sampleId%in%ids))
  stop("Duplicated/foreign heldout animals in DIABLO predictions.")
# Explicit class- and diet-balanced audit of the two original factors.
cross_tab <- table(
  genotype=factor(y,levels=sort(unique(y))),
  diet=factor(diet,levels=sort(unique(diet)))
)
# Must report diet/genotype association rather than falsely calling diet
# randomized or clinically adjusted merely because class balance exists.

# A diet-only predictor trained independently within EACH outer split;
# the baseline tests whether genotype is predictable using diet alone.
baseline <- list()
for(seed in seeds){
  split<-prepare_diablo_holdout_split(blocks,y,seed=seed)
  tr<-split$trainIds
  te<-split$testIds
  dtrain<-diet[tr];dtest<-diet[te]
  train_freq<-table(factor(y[tr],levels=sort(unique(y))))
  majority_global<-names(train_freq)[which.max(train_freq)]
  predictions<-vapply(seq_along(te),function(i){
    local<-y[tr][dtrain==dtest[[i]]]
    if(!length(local))return(majority_global)
    freq<-table(factor(local,levels=sort(unique(y))))
    names(freq)[which.max(freq)]
  },character(1))
  actual<-y[te]
  ber<-mean(vapply(sort(unique(y)),function(label){
    1-mean(predictions[actual==label]==label)
  },numeric(1)))
  baseline[[length(baseline)+1L]]<-data.frame(
    seed=seed,sampleId=te,actual=unname(actual),diet=unname(diet[te]),
    dietOnlyPrediction=predictions,correct=as.integer(predictions==actual),
    stringsAsFactors=FALSE
  )
  if(!is.finite(ber))stop("Non-estimable training-only diet baseline.")
}
diet_baseline<-do.call(rbind,baseline)
diet_ber<-vapply(seeds,function(seed){
 r<-diet_baseline[diet_baseline$seed==seed,,drop=FALSE]
 mean(vapply(sort(unique(y)),function(label)
  1-mean(r$dietOnlyPrediction[r$actual==label]==label),numeric(1)))
},numeric(1))
utils::write.csv(diet_baseline,
  file.path(dest,"diet_only_fold_local_predictions.csv"),row.names=FALSE)

# Negative control 1: permute the genotype at SUBJECT level; retrain all
# fold-specific preprocessing, feature selection and DIABLO coefficients.
# Eight null draws allow no nominal p-value below 1/9; this is a smoke test.
null<-run_diablo_permutation_benchmark(
  blocks,y,file.path(dest,"full_label_permutations"),
  n_permutations=8L,permutation_seed=20261011L,
  split_seeds=seeds,tune=FALSE
)
if(!identical(null$status,"exploratory_null_distribution_not_formal_calibration")||
   length(null$nullBERValues)!=8L ||
   any(!is.finite(null$nullBERValues)) ||
   !isTRUE(all.equal(null$observedMeanBER,observed$summary$meanBER,
                     tolerance=1e-12)) ||
   null$minimumPossibleP!=1/9)
  stop("Full model negative-control or observed DIABLO benchmark mismatch.")

# Negative control 2: deliberately break cross-omic SUBJECT correspondence
# while leaving every lipid vector and every genotype label intact.
# Do not reconstruct a single missing sample from heldout outcomes.
set.seed(20261011L)
permuted_ids<-sample(ids,size=n,replace=FALSE)
if(sum(permuted_ids==ids)>3L)
  stop("Deliberately mispaired lipid block accidentally kept too many samples.")
mispaired<-blocks
mispaired$metabolomics<-lipid[permuted_ids,,drop=FALSE]
rownames(mispaired$metabolomics)<-ids
unpaired<-run_diablo_outer_holdout(
  mispaired,y,file.path(dest,"broken_pairing"),
  seeds=seeds,tune=FALSE
)
if(!identical(unpaired$status,"ok") ||
   any(!is.finite(c(unpaired$summary$meanBER,unpaired$summary$maxBER))) ||
   !identical(sort(unique(unpaired$predictions$seed)),sort(seeds)))
  stop("Broken-pairing falsification did not re-fit both DIABLO fold models.")
# Do not REQUIRE broken pairing to lower accuracy: one view might be
# sufficiently predictive, and a fixed random permutation can be lucky.
# The whole paired-vs-unpaired comparison is retained without threshold.

# Evaluate out-of-fold genotype accuracy by actual diet group. Same animal
# can appear in more than one seed and is NOT an independent replicate.
joined<-merge(observed$predictions,
 data.frame(sampleId=ids,diet=unname(diet),stringsAsFactors=FALSE),
 by="sampleId",sort=FALSE)
if(nrow(joined)!=nrow(observed$predictions))
  stop("Diet-conditional evaluation dropped a real heldout prediction.")
diet_audit<-lapply(sort(unique(diet)),function(k){
 z<-joined[joined$diet==k,,drop=FALSE]
 list(diet=k,evaluations=nrow(z),unique_animals=length(unique(z$sampleId)),
      accuracy=if(nrow(z))mean(z$actual==z$predicted) else NA_real_)
})

audit<-list(
  study="Bioconductor mixOmics Nutrimouse, real paired liver transcriptome and lipid assays",
  pinned_source="mixOmics-org/mixOmics@ef3e760526623d9e91e945e4b50be48d764efc40",
  package_version=as.character(utils::packageVersion("mixOmics")),
  true_reference_method="mixOmics::block.splsda via PMx run_diablo_outer_holdout",
  n_animals=n,
  assays=list(transcriptomics_features=ncol(gene),lipid_features=ncol(lipid)),
  genotype_levels=as.list(table(y)),
  diet_levels=as.list(table(diet)),
  genotype_by_diet=unclass(as.data.frame(cross_tab)),
  split_seeds=as.integer(seeds),
  fixed_model_specification=list(ncomp=1L,keepX_per_view=10L,
      outer_holdout_fraction=.25,training_only_preprocessing=TRUE,
      tuning=FALSE),
  observed=list(mean_heldout_BER=observed$summary$meanBER,
    fold_metrics=observed$summary$folds,
    repeated_holdout_evaluations=nrow(observed$predictions),
    heldout_unique_animals=length(unique(observed$predictions$sampleId)),
    signature_jaccard=observed$summary$signatureStability),
  negative_control_full_training_label_permutations=list(
    n=null$nPermutations,null_mean_BER=null$nullMeanBER,
    null_BERs=as.list(null$nullBERValues),
    empirical_lower_tail_P_with_plus_one_correction=null$lowerTailBEREmpiricalP,
    smallest_possible_P=null$minimumPossibleP,
    NOT_formal_FDR_calibration=TRUE
  ),
  negative_control_lipid_subject_pairing_destroyed=list(
    fixed_random_seed=20261011L,
    paired_animal_matches=sum(permuted_ids==ids),
    mean_heldout_BER=unpaired$summary$meanBER,
    fold_metrics=unpaired$summary$folds,
    no_expected_direction_or_significance_claim=TRUE
  ),
  baseline_diet_only_train_fold_majority=list(
    per_seed_BER=as.list(diet_ber),mean_BER=mean(diet_ber),
    no_holdout_labels_for_fitting=TRUE),
  heldout_diet_strata=diet_audit,
  study_certification="real_study_subject_holdout_not_independent_external_validation",
  scientific_scope=c(
    "Nutrimouse is a single small published controlled animal experiment with 40 subjects, not clinical data.",
    "Diet and genotype can be associated; report both and a training-only diet baseline.",
    "Repeated heldout folds reuse biological animals, so pooled observations are not independent replicates.",
    "Source preprocessing and original acquisition QC are upstream and not independently verified.",
    "No tuning or clinical biomarker selection: fixed keepX 10 per block, one component, two predeclared splits.",
    "Eight full label permutations cannot support formal p < 0.05 or false-discovery guarantees.",
    "Deliberately unpairing one block is a falsification check, not a calibrated hypothesis test.",
    "MOFA2 unsupervised and other multi-omics integration methods are not validated by a DIABLO result.",
    "All model outcomes, including chance or worse performance, are retained; CI success measures execution and audit completeness, NOT high classification accuracy."
  ),
  publication_ready=FALSE,clinical_usable=FALSE
)
jsonlite::write_json(audit,file.path(dest,"real_nutrimouse_diablo_audit.json"),
  auto_unbox=TRUE,pretty=TRUE,na="null",digits=10L)
cat(sprintf(
  "REAL PUBLIC 2-OMICS DIABLO BENCHMARK PASS (no performance threshold): animals=%d genes=%d lipids=%d heldout_BER=%.4f nullLabelBER=%.4f brokenPairBER=%.4f dietOnlyBER=%.4f nNull=%d\n",
  n,ncol(gene),ncol(lipid),observed$summary$meanBER,
  null$nullMeanBER,unpaired$summary$meanBER,mean(diet_ber),
  null$nPermutations
))
