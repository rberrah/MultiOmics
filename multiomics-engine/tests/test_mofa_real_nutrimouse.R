#!/usr/bin/env Rscript
# Actual MOFA2 / mofapy2 on paired, REAL public Nutrimouse omics.
# Nutrimouse RNA panel contains 120 AUTHOR-PRESELECTED genes of ~30k.
# Distinct genotype & diet are experimental biological factors, NOT batch.
# No measured technical batch labels are available: technical separation
# CANNOT be validated on this study. Preserve all negative findings.
#
# Fit: 2 prespecified random seeds; exact 40 paired animals, 3 latent factors.
# Exploratory post-fit audit: genotype conditional on diet, and diet
# conditional on genotype, 199 restricted null permutations, max over all
# factors (no selecting factor after seeing labels and quoting naïve p).
# Sensitivity: one MOFA2 fit with lipid-to-animal pairing deliberately broken;
# a separate label-free per-feature z-score sensitivity on the original pairing.
# No factor is discarded for weak recovery; no superiority threshold.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)
 stop("Usage: Rscript multiomics-engine/tests/test_mofa_real_nutrimouse.R external/mixOmics output_dir")
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/mofa_real_metadata_audit.R")
if(!requireNamespace("MOFA2",quietly=TRUE)||
   !requireNamespace("jsonlite",quietly=TRUE))
 stop("This is a real MOFA2/mofapy2 benchmark: install reference packages.")
e<-new.env(parent=emptyenv())
load(file.path(args[[1L]],"data","nutrimouse.rda"),envir=e)
if(!exists("nutrimouse",envir=e,inherits=FALSE))
 stop("Immutable official Nutrimouse object not found")
study<-e$nutrimouse
if(!is.list(study)||!all(c("gene","lipid","genotype","diet")%in%names(study)))
 stop("Nutrimouse data schema changed: cannot interpret cohorts")
rna<-as.matrix(study$gene)
lipid<-as.matrix(study$lipid)
storage.mode(rna)<-"double"
storage.mode(lipid)<-"double"
if(!identical(dim(rna),c(40L,120L))||
   !identical(dim(lipid),c(40L,21L))||
   any(!is.finite(rna))||any(!is.finite(lipid))||
   is.null(colnames(rna))||is.null(colnames(lipid))||
   anyDuplicated(colnames(rna))||anyDuplicated(colnames(lipid)))
 stop("Original real Nutrimouse assays changed or have invalid features")
ids<-sprintf("NUTRIMOUSE_%03d",seq_len(40L))
rownames(rna)<-ids; rownames(lipid)<-ids
geno<-setNames(as.character(study$genotype),ids)
diet<-setNames(as.character(study$diet),ids)
if(length(geno)!=40L||length(diet)!=40L||
   anyNA(geno)||anyNA(diet)||
   length(unique(geno))!=2L||length(unique(diet))!=5L)
 stop("Nutrimouse biological experiment design changed")
design_check<-model.matrix(~factor(geno)+factor(diet))
if(qr(design_check)$rank!=ncol(design_check))
 stop("Genotype/diet are inseparable, no conditional association claims")
blocks<-list(transcriptomics=rna,lipidomics=lipid)
preflight<-mofa2_input_preflight(blocks,min_samples=32L)
if(!identical(preflight$status,"eligible_for_exploration")||
   preflight$nSharedSubjects!=40L)
 stop("Real MOFA2 pairing/preflight not eligible")
output<-args[[2L]]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
seeds<-c(20261011L,20261012L)
records<-list()
score_sets<-list()
fit_and_audit<-function(experimental_blocks,seed,scenario){
  path<-file.path(output,paste0(scenario,"_seed_",seed))
  # Keep warning evidence (including MOFA's factor-versus-average-level
  # normalization warning) in the archived report, not only in CI stdout.
  fit_warnings<-character()
  fit<-withCallingHandlers(run_mofa2_blocks(
    experimental_blocks,path,factors=3L,seed=seed,
    convergence_mode="fast",min_samples=32L,
    max_features_per_block=NULL
  ), warning=function(w){
    fit_warnings<<-unique(c(fit_warnings,conditionMessage(w)))
  })
  if(!identical(fit$summary$matrixOrientation,
      "feature_rows_subject_columns_for_MOFA2")||
     !identical(fit$summary$implementationContract,
      "mofa2-feature-rows-sample-columns-v2")||
     !isTRUE(fit$summary$subjectIdentityPreserved)||
     fit$summary$samples!=40L)
    stop("Real MOFA2 changed biological sample identity or orientation")
  factor_df<-MOFA2::get_factors(fit$model,as.data.frame=TRUE)
  scores<-mofa_real_factor_matrix(factor_df,ids)
  audit<-mofa_real_metadata_audit(scores,geno,diet,
     n_permutations=199L,seed=20261011L)
  view_qc<-mofa_real_view_level_qc(scores,experimental_blocks,threshold=0.80)
  if(!all(c("genotype_given_diet","diet_given_genotype")%in%
          names(audit$empirical_p_max_across_factors)))
    stop("No genotype/diet conditional association audit available")
  # Simultaneous two biological endpoints, adjusted descriptively.
  # pmin takes attributes from its FIRST argument. pmin(1, named_vector)
  # silently removes the genotype/diet names, corrupting downstream pairing
  # comparisons. Keep and verify experimental-effect identities explicitly.
  bonferroni<-pmin(2*audit$empirical_p_max_across_factors,1)
  names(bonferroni)<-names(audit$empirical_p_max_across_factors)
  if(!identical(sort(names(bonferroni)),
                sort(c("genotype_given_diet","diet_given_genotype"))))
    stop("Bonferroni-adjusted experimental effects lost their identities")
  summary<-list(
    scenario=scenario,seed=seed,
    n_biological_animals=40L,
    effective_factors=ncol(scores),
    factor_names=colnames(scores),
    inputs_preselected_by_publication=TRUE,
    omic_count=2L,
    source_transcriptome_features=120L,
    source_lipid_features=21L,
    genotype_levels=as.list(table(geno)),
    diet_levels=as.list(table(diet)),
    max_conditional_r2=as.list(audit$observed_max_partial_r2),
    conditional_max_factor_empirical_p=as.list(audit$empirical_p_max_across_factors),
    two_endpoint_bonferroni_adjusted_p=as.list(bonferroni),
    factor_vs_assay_mean_qc=view_qc,
    actual_r_fit_warnings=as.list(fit_warnings),
    normalization_issue_not_ruled_out=TRUE,
    permutation_count=audit$n_permutations,
    permutation_scope="Genotype permuted WITHIN diet, diet permuted WITHIN genotype, on fixed unsupervised factors, maximum over all factors for each test",
    model_backend="genuine Bioconductor MOFA2 via mofapy2 and HDF5, PMx run_mofa2_blocks",
    source_data_phenotypes_not_used_to_fit_mofa=TRUE,
    has_measured_technical_batch=FALSE,
    technical_batch_separation_validated=FALSE,
    interpretation="Exploratory retrospective biological associations, no proof of causation or cross-study reproducibility"
  )
  jsonlite::write_json(summary,file.path(path,"metadata_associations.json"),
    pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
  utils::write.csv(audit$observed_partial_r2,
    file.path(path,"factor_conditional_partial_R2.csv"),row.names=TRUE)
  utils::write.csv(audit$null_max_partial_r2,
    file.path(path,"restricted_metadata_nulls.csv"),row.names=FALSE)
  utils::write.csv(scores,file.path(path,"ordered_animal_factor_scores.csv"),
    row.names=TRUE)
  qc_rows<-do.call(rbind,lapply(names(view_qc$views),function(view){
    z<-view_qc$views[[view]]$correlations
    do.call(rbind,lapply(names(z),function(metric)
      data.frame(view=view,metric=metric,factor=names(z[[metric]]),
        correlation=as.numeric(z[[metric]]),stringsAsFactors=FALSE)))
  }))
  utils::write.csv(qc_rows,file.path(path,"factor_vs_view_mean_qc.csv"),row.names=FALSE)
  cat(sprintf("REAL MOFA2 NUTRIMOUSE FIT PASS: scenario=%s seed=%d n=40 factors=%d maxR2_genotype_given_diet=%.4f maxR2_diet_given_genotype=%.4f p_genotype=%.4f p_diet=%.4f\n",
    scenario,seed,ncol(scores),
    audit$observed_max_partial_r2[["genotype_given_diet"]],
    audit$observed_max_partial_r2[["diet_given_genotype"]],
    audit$empirical_p_max_across_factors[["genotype_given_diet"]],
    audit$empirical_p_max_across_factors[["diet_given_genotype"]]))
  list(summary=summary,scores=scores)
}
for(seed in seeds){
  z<-fit_and_audit(blocks,seed,"paired")
  records[[length(records)+1L]]<-z$summary
  score_sets[[as.character(seed)]]<-z$scores
}
stability<-mofa_real_between_seeds(
  score_sets[[as.character(seeds[[1L]])]],
  score_sets[[as.character(seeds[[2L]])]]
)
# Sensitivity: break biological pairing by reassigning whole lipid profiles.
# The "omics are now on the wrong animal" perturbation is NOT a legitimate
# platform input; score difference should be disclosed but never used to
# choose an omic integration method from results on the same mice.
set.seed(20261011L)
shuffled<-sample(ids,size=length(ids),replace=FALSE)
if(sum(shuffled==ids)>3L)
 stop("Nutrimouse deliberately broken lipid pairing retained too many true pairs")
mispaired<-blocks
mispaired$lipidomics<-lipid[match(shuffled,ids),,drop=FALSE]
rownames(mispaired$lipidomics)<-ids
perturbed<-fit_and_audit(mispaired,seeds[[1L]],"broken_lipid_subject_pairing")
records[[length(records)+1L]]<-perturbed$summary
# We prespecify a second pipeline with every feature individually centered
# and scaled (view-wide scaling is already enabled by the backend).
# This tests normalization sensitivity; it is NOT an independent cohort.
zscore_blocks<-mofa_real_zscore_features(blocks)
zscore_fit<-fit_and_audit(zscore_blocks,seeds[[1L]],
                         "label_free_feature_zscore_paired_sensitivity")
records[[length(records)+1L]]<-zscore_fit$summary
normalization_sensitivity<-list(
  method="Mean-zero / SD-one per assayed feature within each view; no phenotype labels accessed",
  different_preprocessing_from_reference=TRUE,
  predeclared_paired_seed=seeds[[1L]],
  paired_original_max_r2=records[[1L]]$max_conditional_r2,
  paired_zscore_max_r2=zscore_fit$summary$max_conditional_r2,
  original_view_level_qc=records[[1L]]$factor_vs_assay_mean_qc,
  zscore_view_level_qc=zscore_fit$summary$factor_vs_assay_mean_qc,
  normalization_fix_proven=FALSE,
  interpretation="Same 40 animals; changed pre-processing can change latent factors and associations, but cannot determine whether factor-to-mean correlation is biological or technical."
)
jsonlite::write_json(normalization_sensitivity,
  file.path(output,"feature_normalization_sensitivity.json"),
  pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
pairing_falsification<-mofa_real_pairing_falsification(
  records[[1L]],perturbed$summary)
jsonlite::write_json(pairing_falsification,
  file.path(output,"pairing_falsification_negative_result.json"),
  pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
cat(sprintf("REAL MOFA2 NEGATIVE FALSIFICATION: both_genotype_diet_associations_persist=%s; paired_R2=%.4f/%.4f; broken_R2=%.4f/%.4f; no_synergy_inference\n",
  pairing_falsification$both_associations_still_significant,
  records[[1L]]$max_conditional_r2$genotype_given_diet,
  records[[1L]]$max_conditional_r2$diet_given_genotype,
  perturbed$summary$max_conditional_r2$genotype_given_diet,
  perturbed$summary$max_conditional_r2$diet_given_genotype))
# Independent BASE-R single-omic PCA comparators, computed strictly without
# phenotype labels. Keep ALL four prespecified view/preprocessing baselines,
# not only a PC or arm with the strongest apparent genotype/diet association.
# Identical 40 animals, identical 3 post-fit tested factors and 199 stratified
# permutations as MOFA2. This is descriptive and not a paired predictive
# test, a model-agnostic synergy estimate or external confirmation.
pca_comparators<-list()
for(view in names(blocks)){
  for(preprocessing in c("raw_centered","per_feature_zscore")){
    standardize<-identical(preprocessing,"per_feature_zscore")
    baseline<-mofa_real_single_view_pca(blocks[[view]],
      n_components=3L,feature_standardization=standardize)
    baseline_audit<-mofa_real_metadata_audit(
      baseline$scores,geno,diet,n_permutations=199L,seed=20261011L)
    key<-paste(view,preprocessing,sep="__")
    pvals<-baseline_audit$empirical_p_max_across_factors
    if(!identical(sort(names(pvals)),
                  sort(c("genotype_given_diet","diet_given_genotype"))))
      stop("Single-omic comparator outcomes lost their metadata identities")
    pca_comparators[[key]]<-list(
      view=view,preprocessing=preprocessing,
      animals=nrow(baseline$scores),features=baseline$measured_feature_count,
      components=baseline$n_components,
      variance_explained_by_PC=baseline$explained_variance_fraction,
      max_partial_r2=as.list(baseline_audit$observed_max_partial_r2),
      empirical_p_max_across_three_PCs=as.list(pvals),
      bonferroni_adjusted_across_four_views_and_two_endpoints=
        as.list(stats::setNames(pmin(1,8*pvals),names(pvals))),
      permutations=199L,preprocessing_label_blind=TRUE,
      source_feature_selection_by_original_authors=TRUE,
      interpretation="Exploratory retrospective single-omic PCA. The four arms and two phenotype endpoints are not independent tests; adjusted p-values are a conservative post-fit descriptive screen, not external validation or evidence of integration synergy."
    )
    utils::write.csv(baseline$scores,
      file.path(output,paste0("single_omic_pca_",key,"_scores.csv")),
      row.names=TRUE)
    cat(sprintf("NUTRIMOUSE SINGLE-OMIC PCA BASELINE: view=%s preprocessing=%s geno_cond_R2=%.4f diet_cond_R2=%.4f p_geno=%.4f p_diet=%.4f (4 arms + 2 endpoints, NOT a synergy test)\n",
      view,preprocessing,
      baseline_audit$observed_max_partial_r2[["genotype_given_diet"]],
      baseline_audit$observed_max_partial_r2[["diet_given_genotype"]],
      pvals[["genotype_given_diet"]],pvals[["diet_given_genotype"]]))
  }
}
if(length(pca_comparators)!=4L)
  stop("Every prespecified single-view comparator must be retained")
pca_report<-list(
  scenario="Four unsupervised single-omic PCA reference analyses on real paired Nutrimouse",
  views=names(blocks),preprocessing=c("raw_centered","per_feature_zscore"),
  subjects=40L,per_arm_factors=3L,per_arm_postfit_permutations=199L,
  any_phenotype_used_for_preprocessing_or_pca=FALSE,
  full_original_unfiltered_transcriptome=FALSE,
  both_views_assayed_on_same_original_animals=TRUE,
  comparator_kind="single-omic PCA, not single-view MOFA2 or single-view prediction",
  comparisons=pca_comparators,
  caveat="Descriptive associations cannot show incremental value of multi-omic integration. Comparability of partial R2 is not a matched cross-validated predictive performance comparison; both feature selection by original authors and repeated inspection of the same cohort prevent confirmatory claims."
)
jsonlite::write_json(pca_report,
  file.path(output,"single_omic_pca_comparison.json"),
  pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
report<-list(
 study="Real Nutrimouse, true paired hepatic transcriptomic and lipid assays",
 source="mixOmics-org/mixOmics@ef3e760526623d9e91e945e4b50be48d764efc40",
 reference_package="Bioconductor MOFA2",
 package_version=as.character(utils::packageVersion("MOFA2")),
 exact_model_engine="PMx run_mofa2_blocks -> MOFA2::run_mofa(use_basilisk=TRUE)",
 total_biological_animals=40L,
 initializations=as.list(as.integer(seeds)),
 assayed_views=list(transcriptomics=120L,lipidomics=21L),
 gene_selection_caveat="Published transcriptomic panel 120 genes was preselected from roughly 30000 original genes in source study.",
 genotype_diet_not_technical_batch="Both are experimental biological/exposure labels; no machine batch IDs exist and technical batch is NOT assessed",
 model_fits=records,
 paired_two_seed_factor_correlations=stability,
 pairing_falsification_negative_result=pairing_falsification,
 single_omic_pca_comparators=pca_report,
 single_omic_comparison_cannot_establish_synergy=TRUE,
 per_fit_normalization_qc_included=TRUE,
 featurewise_normalization_sensitivity=normalization_sensitivity,
 featurewise_normalization_is_unsupervised_same_cohort=TRUE,
 normalization_artifact_not_excluded=TRUE,
 broken_pairing_matches=sum(shuffled==ids),
 broken_pairing_was_chosen_in_advance=TRUE,
 one_broken_pairing_seed_only=TRUE,
 no_technical_batch_validation=TRUE,
 no_external_cohort=TRUE,
 no_untouched_confirmatory_study=TRUE,
 no_supervised_biomarker_performance_claim=TRUE,
 caveats=c(
   "One retrospective 40-mouse study with author-preselected features; no unfiltered RNA-seq.",
   "Latent factor association is fitted and tested ON SAME experimental animals; label permutations are post-fit descriptive, not external confirmation.",
   "Genotype may be associated with diet; conditional partial-R2 and within-other-factor stratified nulls limit but cannot eliminate confounding.",
   "MOFA factors may rotate, flip sign, reorder or disappear across model initializations.",
   "Even with the lipid-animal matching broken, latent factors may remain associated with genotype/diet; these same-animal post-fit tests do NOT validate multimodal pairing or added value.",
   "A MOFA warning about latent factors tracking average sample expression must be archived with factor-vs-view mean QC; this may be biology, technical variation or normalization, not a validated biological mechanism.",
   "Intentionally mismatching lipids is a falsification stress test, not a formal null or proof of synergy.",
   "No technical batch metadata exists; cannot claim biological-vs-technical factor separation on this dataset.",
   "Four MOFA model runs (two paired seeds, one broken-pairing, one featurewise-zscore paired sensitivity) do not establish best factor number, generalizability, clinical validity or a validated biological mechanism.",
   "All associations, including null effects, are recorded, never filtered to winning seeds."
 ),
 publication_ready=FALSE,clinical_certification=FALSE
)
jsonlite::write_json(report,file.path(output,"real_nutrimouse_mofa2_report.json"),
  pretty=TRUE,auto_unbox=TRUE,na="null",digits=12L)
cat(sprintf("REAL MOFA2 NUTRIMOUSE ALL INPUTS AUDITED PASS: 40 mice, true paired seeds=2, broken pairing=1, feature_zscore_paired=1, no technical batch and NO CONFIRMATORY CLAIM\n"))
