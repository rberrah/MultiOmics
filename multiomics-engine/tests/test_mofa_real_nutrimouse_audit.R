#!/usr/bin/env Rscript
# Strict base-R tests: real MOFA2 audit doesn't need MOFA2 installed.
source("multiomics-engine/mofa_real_metadata_audit.R")
set.seed(20261011L)
n<-40L
ids<-sprintf("M%03d",seq_len(n))
diet<-setNames(rep(paste0("diet",1:5),each=8L),ids)
geno<-setNames(rep(rep(c("WT","KO"),each=4L),5L),ids)
x<-ifelse(geno=="KO",2,-2)+stats::rnorm(n,0,.1)
z<-seq_len(n)/40+stats::rnorm(n,0,.1)
f<-rbind(
  data.frame(sample=ids,factor="Factor2",value=as.numeric(-z)),
  data.frame(sample=ids,factor="Factor1",value=as.numeric(x))
)
mat<-mofa_real_factor_matrix(f,ids)
stopifnot(identical(dim(mat),c(40L,2L)),
          identical(rownames(mat),ids),
          all(mat[,2L]==x))
audit<-mofa_real_metadata_audit(mat,geno,diet,49L,20261011L)
stopifnot(audit$status=="retrospective_exploratory_post_fit_associations",
  audit$observed_max_partial_r2[["genotype_given_diet"]]>.90,
  all(audit$empirical_p_max_across_factors>=1/50),
  identical(dim(audit$null_max_partial_r2),c(49L,2L)),
  all(is.finite(audit$null_max_partial_r2)),
  all(audit$null_max_partial_r2>=0),
  all(audit$null_max_partial_r2<=1))
again<-mofa_real_metadata_audit(mat,geno,diet,49L,20261011L)
stopifnot(identical(audit$null_max_partial_r2,again$null_max_partial_r2))
set.seed(1)
pg<-mofa_stratum_permutation(as.character(geno),as.character(diet))
for(k in unique(diet)){
 ix<-which(diet==k)
 stopifnot(identical(sort(pg[ix]),sort(unname(geno[ix]))))
}
# Factor sign flips and factor order are irrelevant to reported max-partial
# R2. Best seed match compares full profiles and discloses nonunique links.
swapped<-mat[,2:1,drop=FALSE]
swapped[,1L]<- -swapped[,1L]
stab<-mofa_real_between_seeds(mat,swapped)
stopifnot(all(stab$first_seed_factor_best_abs_correlation>.999999),
          isTRUE(stab$all_factors_match_uniquely))
perturbed<-f[-1L,,drop=FALSE]
stopifnot(inherits(try(mofa_real_factor_matrix(perturbed,ids),silent=TRUE),"try-error"))
bad<-f;bad$sample[[1L]]<-bad$sample[[2L]]
stopifnot(inherits(try(mofa_real_factor_matrix(bad,ids),silent=TRUE),"try-error"))
bad<-f;bad$value[[1L]]<-Inf
stopifnot(inherits(try(mofa_real_factor_matrix(bad,ids),silent=TRUE),"try-error"))
bad_diet<-diet;bad_diet[[1L]]<-NA_character_
stopifnot(inherits(try(mofa_real_metadata_audit(mat,geno,bad_diet,19L),silent=TRUE),"try-error"))
bad_ids<-geno;names(bad_ids)[1L]<-"UNKNOWN"
stopifnot(inherits(try(mofa_real_metadata_audit(mat,bad_ids,diet,19L),silent=TRUE),"try-error"))
confounded<-setNames(ifelse(geno=="KO","high","low"),ids)
stopifnot(inherits(try(mofa_real_metadata_audit(mat,geno,confounded,19L),silent=TRUE),"try-error"))
stopifnot(inherits(try(mofa_real_metadata_audit(mat,geno,diet,8L),silent=TRUE),"try-error"))
# Deterministic QC detects a factor tracking per-animal assay-wide level,
# without calling this correlation a biological or technical mechanism.
input_signal<-matrix(stats::rnorm(n*12),nrow=n,ncol=12,
  dimnames=list(ids,paste0("rna",1:12)))
input_signal<-sweep(input_signal,1L,seq_len(n)*0.5,"+")
input_lipid<-matrix(stats::rnorm(n*6),nrow=n,ncol=6,
  dimnames=list(ids,paste0("lipid",1:6)))
qc_factor<-cbind(Level=seq_len(n)+stats::rnorm(n,0,0.001),
                 Independent=stats::rnorm(n))
rownames(qc_factor)<-ids
quality<-mofa_real_view_level_qc(qc_factor,
 list(transcriptomics=input_signal,lipidomics=input_lipid),
 threshold=0.80)
stopifnot(isTRUE(quality$any_view_level_warning),
 isTRUE(quality$views$transcriptomics$view_level_warning),
 quality$views$transcriptomics$max_absolute_correlations$raw_feature_mean>0.99,
 all(c("Level","Independent")%in%
   names(quality$views$lipidomics$correlations$raw_feature_mean)))
bad_lipid<-input_lipid[n:1,,drop=FALSE]
stopifnot(inherits(try(mofa_real_view_level_qc(qc_factor,
 list(transcriptomics=input_signal,lipidomics=bad_lipid)),
 silent=TRUE),"try-error"))
# Label-free z-score must preserve biological subject and feature identity,
# and refuse silent feature dropping / sample reordering.
scaled<-mofa_real_zscore_features(list(transcriptomics=input_signal,
                                      lipidomics=input_lipid))
stopifnot(identical(rownames(scaled$transcriptomics),ids),
          identical(colnames(scaled$transcriptomics),colnames(input_signal)),
          max(abs(colMeans(scaled$lipidomics)))<1e-10,
          max(abs(apply(scaled$lipidomics,2L,stats::sd)-1))<1e-10)
constant<-input_lipid; constant[,1L]<-0
stopifnot(inherits(try(mofa_real_zscore_features(list(
 transcriptomics=input_signal,lipidomics=constant)),
 silent=TRUE),"try-error"))
stopifnot(inherits(try(mofa_real_zscore_features(list(
 transcriptomics=input_signal,lipidomics=input_lipid[n:1,,drop=FALSE])),
 silent=TRUE),"try-error"))
# Paired, label-blind single-view PCA is deterministic, keeps exactly the
# measured feature and biological subject identities and never uses outcomes
# until the existing post-fit partial-R2/permutation audit.
pca_raw<-mofa_real_single_view_pca(input_signal,3L,FALSE)
pca_scaled<-mofa_real_single_view_pca(input_signal,3L,TRUE)
stopifnot(identical(dim(pca_raw$scores),c(n,3L)),
 identical(rownames(pca_scaled$scores),ids),
 identical(colnames(pca_raw$scores),paste0("PC",1:3)),
 isFALSE(pca_raw$labels_supplied_to_pca),
 max(abs(colMeans(pca_scaled$scores)))<1e-9,
 sum(unlist(pca_raw$explained_variance_fraction))<=1+1e-8,
 all(unlist(pca_scaled$explained_variance_fraction)>=0))
again_pca<-mofa_real_single_view_pca(input_signal,3L,TRUE)
stopifnot(isTRUE(all.equal(pca_scaled$scores,again_pca$scores,
                            tolerance=1e-10)))
pca_audit<-mofa_real_metadata_audit(pca_scaled$scores,geno,diet,
                                    n_permutations=19L,seed=20261011L)
stopifnot(length(pca_audit$observed_max_partial_r2)==2L,
          nrow(pca_audit$observed_partial_r2)==3L,
          all(pca_audit$empirical_p_max_across_factors>=1/20))
wrong_ids<-input_signal;rownames(wrong_ids)[1L]<-rownames(wrong_ids)[2L]
stopifnot(inherits(try(mofa_real_single_view_pca(wrong_ids),
                       silent=TRUE),"try-error"),
          inherits(try(mofa_real_single_view_pca(constant),
                       silent=TRUE),"try-error"))
wrong_metadata<-geno;names(wrong_metadata)[1L]<-"not-a-subject"
stopifnot(inherits(try(mofa_real_metadata_audit(
  pca_scaled$scores,wrong_metadata,diet,19L),silent=TRUE),"try-error"))
paired_summary<-list(
 max_conditional_r2=list(genotype_given_diet=0.5313,
                         diet_given_genotype=0.6846),
 two_endpoint_bonferroni_adjusted_p=list(genotype_given_diet=0.01,
                                          diet_given_genotype=0.01))
broken_summary<-list(
 max_conditional_r2=list(genotype_given_diet=0.5568,
                         diet_given_genotype=0.5032),
 two_endpoint_bonferroni_adjusted_p=list(genotype_given_diet=0.01,
                                          diet_given_genotype=0.01))
negative<-mofa_real_pairing_falsification(paired_summary,broken_summary)
stopifnot(identical(negative$status,
 "descriptive_negative_control_not_external_validation"),
 isTRUE(negative$both_associations_still_significant),
 identical(negative$valid_cross_omic_synergy_conclusion,FALSE),
 abs(negative$broken_minus_paired_partial_r2$genotype_given_diet-0.0255)<1e-10)
broken_summary$two_endpoint_bonferroni_adjusted_p$diet_given_genotype<-0.8
other<-mofa_real_pairing_falsification(paired_summary,broken_summary)
stopifnot(!other$both_associations_still_significant,
          isFALSE(other$associations_still_significant_after_breaking_pairing$diet_given_genotype))
broken_summary$max_conditional_r2$diet_given_genotype<-Inf
stopifnot(inherits(try(mofa_real_pairing_falsification(
 paired_summary,broken_summary),silent=TRUE),"try-error"))

cat("MOFA2 real metadata audit base-R contracts: paired ID, genotype/diet adjustment, permutation max factors, RNG replay, factor sign/order, confounding refusal and falsification QC PASS\n")
