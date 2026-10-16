# Pure base-R, post-fit association audit of UNSUPERVISED MOFA2 factors.
# Does not fit any model or use outcomes in factor generation.
# Do NOT treat experimental diet as batch: it is an experimental exposure.
mofa_real_factor_matrix <- function(factors, subject_ids) {
  if(!is.data.frame(factors) ||
     !all(c("sample","factor","value")%in%names(factors)))
    stop("MOFA2 must supply sample/factor/value long format")
  ids <- as.character(subject_ids)
  if(length(ids)<16L || anyNA(ids)||anyDuplicated(ids)||
     any(!nzchar(ids)))
    stop("MOFA2 audit requires distinct known biological subjects")
  sample<-as.character(factors$sample)
  factor<-as.character(factors$factor)
  value<-suppressWarnings(as.numeric(factors$value))
  if(anyNA(sample)||anyNA(factor)||any(!nzchar(sample))||
     any(!nzchar(factor))||any(!is.finite(value))||
     anyDuplicated(paste(sample,factor,sep="\r")))
    stop("Missing, nonfinite or duplicate MOFA2 factor scores")
  if(!setequal(sample,ids))
    stop("Factor sample sets differ from original assay biological subjects")
  names_f <- unique(factor)
  if(length(names_f)<1L||length(names_f)>10L)
    stop("MOFA2 factor count outside allowed 1..10")
  z<-vapply(names_f,function(name) {
     ix<-which(factor==name)
     if(length(ix)!=length(ids))
       stop("MOFA2 factor is missing some biological subjects")
     v<-value[ix][match(ids,sample[ix])]
     if(any(!is.finite(v)))stop("Unmatched factor score")
     v
  },numeric(length(ids)))
  if(is.null(dim(z)))z<-matrix(z,nrow=length(ids),ncol=1L)
  dimnames(z)<-list(ids,names_f)
  if(any(apply(z,2L,stats::sd)<=1e-8))
    stop("Degenerate MOFA2 factor scores: associations cannot be estimated")
  z
}

mofa_partial_r2 <- function(values,geno,diet) {
  if(length(values)!=length(geno)||length(values)!=length(diet)||
     any(!is.finite(values))||anyNA(geno)||anyNA(diet))
    stop("MOFA2 covariate association needs complete scores/genotype/diet")
  meta<-data.frame(geno=factor(geno),diet=factor(diet))
  if(nlevels(meta$geno)<2L||nlevels(meta$diet)<2L)
    stop("At least two observed levels per experimental variable")
  x_full<-stats::model.matrix(~geno+diet,meta)
  x_diet<-stats::model.matrix(~diet,meta)
  x_geno<-stats::model.matrix(~geno,meta)
  if(qr(x_full)$rank!=ncol(x_full))
    stop("Genotype and diet cannot be separated by the experimental design")
  rss<-function(X)sum(stats::lm.fit(X,values)$residuals^2)
  full<-rss(x_full)
  base_diet<-rss(x_diet)
  base_geno<-rss(x_geno)
  if(base_diet<=1e-12||base_geno<=1e-12)
    stop("Partial R2 is undefined with zero conditional residual variance")
  # Both are variance fractions CONDITIONAL on the other factor;
  # neither proves causality, independent biological signal nor an omic effect.
  c(
   genotype_given_diet=pmax(0,pmin(1,(base_diet-full)/base_diet)),
   diet_given_genotype=pmax(0,pmin(1,(base_geno-full)/base_geno))
  )
}

mofa_stratum_permutation <- function(labels,strata) {
  if(length(labels)!=length(strata)||anyNA(labels)||anyNA(strata))
    stop("Permutation requires complete labels and strata")
  result<-labels
  for(ix in split(seq_along(labels),factor(strata)))
    if(length(ix)>1L)result[ix]<-sample(labels[ix],length(ix),replace=FALSE)
  if(anyNA(result))stop("Permutation produced missing labels")
  result
}

mofa_real_metadata_audit <- function(scores,genotype,diet,
                                     n_permutations=199L,seed=20261011L) {
  if(!is.matrix(scores)||is.null(rownames(scores))||
     length(genotype)!=nrow(scores)||length(diet)!=nrow(scores))
    stop("MOFA2 post-fit scores and metadata cannot be paired")
  ids<-rownames(scores)
  if(anyDuplicated(ids)||anyNA(ids)||any(!is.finite(scores))||
     is.null(names(genotype))||is.null(names(diet))||
     anyDuplicated(names(genotype))||anyDuplicated(names(diet))||
     !setequal(names(genotype),ids)||!setequal(names(diet),ids))
    stop("Score/sample metadata identities must match without dropping anyone")
  if(length(n_permutations)!=1L||n_permutations<19L||
     n_permutations!=as.integer(n_permutations))
    stop("A prespecified >=19 metadata permutations are required")
  g<-as.character(genotype[ids])
  d<-as.character(diet[ids])
  factor_assoc<-function(gg,dd) {
     t(vapply(seq_len(ncol(scores)),function(j)
         mofa_partial_r2(scores[,j],gg,dd),numeric(2L)))
  }
  observed<-factor_assoc(g,d)
  colnames(observed)<-c("genotype_given_diet","diet_given_genotype")
  rownames(observed)<-colnames(scores)
  # Max across all MOFA factors BEFORE null comparison; no selecting the
  # best-looking factor and then quoting an unadjusted per-factor p-value.
  obs_max<-apply(observed,2L,max)
  set.seed(as.integer(seed))
  null<-matrix(NA_real_,nrow=n_permutations,ncol=2L,
    dimnames=list(NULL,colnames(observed)))
  for(i in seq_len(n_permutations)){
    perm_g<-mofa_stratum_permutation(g,d) # conditional on diet
    perm_d<-mofa_stratum_permutation(d,g) # conditional on genotype
    a<-factor_assoc(perm_g,d)
    b<-factor_assoc(g,perm_d)
    null[i,1L]<-max(a[,1L])
    null[i,2L]<-max(b[,2L])
  }
  p<-vapply(seq_len(ncol(null)),function(j)
     (1+sum(null[,j]>=obs_max[[j]]))/(n_permutations+1),numeric(1))
  names(p)<-colnames(null)
  list(
    status="retrospective_exploratory_post_fit_associations",
    observed_partial_r2=observed,
    observed_max_partial_r2=obs_max,
    null_max_partial_r2=null,
    empirical_p_max_across_factors=p,
    n_permutations=as.integer(n_permutations),
    min_possible_p=1/(n_permutations+1),
    design_adjustment="Genotype permutation within diet, diet permutation within genotype; partial R2 conditional on other experimental factor",
    caveat="Unsupervised factor scores were fitted without metadata; permutations are post-fit exploratory and do not provide independent biological validation, correct selected-source bias or ensure causal claims."
  )
}

mofa_real_between_seeds <- function(first,second) {
  if(!is.matrix(first)||!is.matrix(second)||
     !identical(rownames(first),rownames(second))||
     ncol(first)<1L||ncol(second)<1L)
    stop("MOFA2 seeds must contain exactly identical ordered biological subjects")
  correlations<-abs(stats::cor(first,second))
  if(any(!is.finite(correlations)))
    stop("Unidentifiable cross-seed latent factor correlations")
  # Do not silently assume factor order or sign is stable. Best matches
  # here are descriptive; duplicates are ALLOWED and disclosed explicitly,
  # never described as guaranteed one-to-one recovery.
  best<-apply(correlations,1L,max)
  matches<-colnames(correlations)[max.col(correlations,ties.method="first")]
  list(
    absolute_correlation=correlations,
    first_seed_factor_best_match=stats::setNames(matches,rownames(correlations)),
    first_seed_factor_best_abs_correlation=stats::setNames(best,rownames(correlations)),
    all_factors_match_uniquely=!anyDuplicated(matches),
    interpretation="Descriptive full-cohort factor reproducibility across two initializations; factors may rotate, permute or flip sign. Matching is not proof of stable biological signatures."
  )
}


# Measure whether MOFA scores track a view-wide intensity offset. This is a QC
# association on the SAME animals, NOT a correction, confounder test or proof
# of a technical artifact. Raw and feature-standardised row means are reported.
mofa_real_view_level_qc <- function(scores, blocks, threshold=0.80) {
  ids<-rownames(scores)
  if(!is.matrix(scores)||is.null(ids)||anyDuplicated(ids)||
     any(!is.finite(scores))||!is.list(blocks)||length(blocks)<2L||
     is.null(names(blocks))||anyDuplicated(names(blocks))||
     any(!nzchar(names(blocks)))||length(threshold)!=1L||
     !is.finite(threshold)||threshold<=0||threshold>=1)
    stop("Invalid MOFA scores, views or view-mean correlation threshold")
  results<-list()
  for(view in names(blocks)){
    X<-blocks[[view]]
    if(!is.matrix(X)||is.null(rownames(X))||
       !identical(rownames(X),ids)||ncol(X)<2L||
       any(!is.finite(X))||!is.numeric(X))
      stop("MOFA QC requires exactly ordered finite feature matrices")
    raw_mean<-rowMeans(X)
    sd_by_feature<-apply(X,2L,stats::sd)
    keep<-which(is.finite(sd_by_feature)&sd_by_feature>1e-12)
    if(length(keep)<2L)stop("MOFA QC needs at least two varying features per view")
    standardized<-sweep(sweep(X[,keep,drop=FALSE],2L,
                                   colMeans(X[,keep,drop=FALSE]),"-"),
                        2L,sd_by_feature[keep],"/")
    means<-list(raw_feature_mean=raw_mean,
                standardized_feature_mean=rowMeans(standardized))
    correlations<-lapply(means,function(v){
      if(stats::sd(v)<=1e-12)
        return(stats::setNames(rep(NA_real_,ncol(scores)),colnames(scores)))
      z<-vapply(seq_len(ncol(scores)),
                function(j)stats::cor(v,scores[,j]),numeric(1L))
      stats::setNames(z,colnames(scores))
    })
    maxima<-vapply(correlations,function(z)
      if(all(is.na(z)))NA_real_ else max(abs(z),na.rm=TRUE),numeric(1L))
    results[[view]]<-list(
      feature_count=ncol(X),varying_feature_count=length(keep),
      correlations=correlations,
      max_absolute_correlations=as.list(maxima),
      view_level_warning=any(is.finite(maxima)&maxima>=threshold),
      raw_feature_mean_sd=stats::sd(raw_mean),
      standardized_feature_mean_sd=stats::sd(means$standardized_feature_mean)
    )
  }
  list(
    threshold=threshold,views=results,
    any_view_level_warning=any(vapply(results,function(x)x$view_level_warning,logical(1L))),
    caveat="Correlation of factors with average sample signal can reflect biological differences, normalization or technical artifacts. This descriptive diagnostic cannot determine the cause, distinguish biology from batch, or rescue a mispaired study."
  )
}

# A negative falsification result is evidence, not an error to suppress.
# The fixed post-fit phenotype association test remains significant even when
# assays are deliberately matched to the WRONG animals in this benchmark.
mofa_real_pairing_falsification <- function(paired,broken) {
  keys<-c("genotype_given_diet","diet_given_genotype")
  unpack<-function(record,field){
    if(!is.list(record)||is.null(record[[field]]))
      stop("Paired and falsified summaries must contain identical QC fields")
    z<-unlist(record[[field]],use.names=TRUE)
    if(!identical(sort(names(z)),sort(keys))||
       length(z)!=2L||any(!is.finite(z)))
      stop("Invalid post-fit phenotype fields")
    z[keys]
  }
  p_r2<-unpack(paired,"max_conditional_r2")
  b_r2<-unpack(broken,"max_conditional_r2")
  p_p<-unpack(paired,"two_endpoint_bonferroni_adjusted_p")
  b_p<-unpack(broken,"two_endpoint_bonferroni_adjusted_p")
  if(any(p_r2<0|p_r2>1|b_r2<0|b_r2>1)||
     any(p_p<0|p_p>1|b_p<0|b_p>1))
    stop("Invalid partial R2 or multiplicity-adjusted p-values")
  persistence<-b_p<=0.05
  names(persistence)<-keys
  list(
    status="descriptive_negative_control_not_external_validation",
    observed_paired_partial_r2=as.list(p_r2),
    broken_pairing_partial_r2=as.list(b_r2),
    broken_minus_paired_partial_r2=as.list(b_r2-p_r2),
    broken_pairing_bonferroni_p=as.list(b_p),
    associations_still_significant_after_breaking_pairing=as.list(persistence),
    both_associations_still_significant=isTRUE(all(persistence)),
    valid_cross_omic_synergy_conclusion=FALSE,
    interpretation="Genotype/diet associations in latent scores can persist when lipids are assigned to the wrong animals. Even a low post-fit permutation p is not proof of correct assay pairing, multi-omic added value, independent replication or causality. One broken pairing is a sensitivity check, not a formal null distribution."
  )
}


# Unsupervised sensitivity ONLY. Standardise every measured feature across the
# same study subjects with no genotype/diet information. This is deliberately
# not used as an estimate of out-of-sample preprocessing quality.
mofa_real_zscore_features <- function(blocks) {
  if(!is.list(blocks)||length(blocks)<2L||is.null(names(blocks))||
     anyDuplicated(names(blocks))||any(!nzchar(names(blocks))))
    stop("Per-feature normalization needs at least two named assay views")
  ids<-rownames(blocks[[1L]])
  if(is.null(ids)||anyDuplicated(ids)||anyNA(ids))
    stop("Per-feature normalization requires exact biological IDs")
  out<-lapply(blocks,function(X){
    if(!is.matrix(X)||!is.numeric(X)||
       !identical(rownames(X),ids)||ncol(X)<2L||
       any(!is.finite(X))||anyDuplicated(colnames(X)))
      stop("Cannot normalize missing or mispaired assay values")
    s<-apply(X,2L,stats::sd)
    if(any(!is.finite(s))||any(s<=1e-12))
      stop("Feature normalization fails on constant/non-identifiable features")
    z<-sweep(sweep(X,2L,colMeans(X),"-"),2L,s,"/")
    if(any(!is.finite(z))||max(abs(colMeans(z)))>1e-9||
       max(abs(apply(z,2L,stats::sd)-1))>1e-9)
      stop("Feature-normalized assay has an invalid numeric mean/SD")
    z
  })
  names(out)<-names(blocks)
  out
}


# Matched single-omic UNSUPERVISED PCA baseline for the real Nutrimouse
# benchmark. Does not consume genotype/diet; labels enter only the identical
# post-fit descriptive partial-R2/permutation audit used for MOFA2.
# A PCA result is not a direct same-model comparison with MOFA2, and phenotype
# association of either method cannot establish genuine cross-omic synergy.
mofa_real_single_view_pca <- function(X, n_components=3L,
                                      feature_standardization=FALSE) {
  if(!is.matrix(X)||!is.numeric(X)||nrow(X)<16L||
     ncol(X)<n_components||is.null(rownames(X))||
     is.null(colnames(X))||anyDuplicated(rownames(X))||
     anyDuplicated(colnames(X))||anyNA(rownames(X))||
     anyNA(colnames(X))||any(!nzchar(rownames(X)))||
     any(!nzchar(colnames(X)))||any(!is.finite(X))||
     length(n_components)!=1L||!is.finite(n_components)||
     n_components!=as.integer(n_components)||
     n_components<1L||n_components>10L||
     length(feature_standardization)!=1L||
     !is.logical(feature_standardization)||
     is.na(feature_standardization))
    stop("Single-view PCA requires complete finite uniquely identified samples/features")
  s<-apply(X,2L,stats::sd)
  if(any(!is.finite(s))||any(s<=1e-12))
    stop("PCA baseline cannot silently filter invariant variables")
  fit<-stats::prcomp(X,center=TRUE,scale.=feature_standardization,
                     rank.=as.integer(n_components))
  scores<-fit$x[,seq_len(n_components),drop=FALSE]
  colnames(scores)<-paste0("PC",seq_len(n_components))
  if(!identical(rownames(scores),rownames(X))||
     any(!is.finite(scores))||
     any(apply(scores,2L,stats::sd)<=1e-8))
    stop("PCA returned degenerate scores or changed subject identities")
  total_variance<-if(feature_standardization) ncol(X) else sum(s^2)
  explained<-vapply(seq_len(n_components),function(j)
    stats::var(scores[,j])/total_variance,numeric(1L))
  names(explained)<-colnames(scores)
  if(any(!is.finite(explained))||any(explained<0)||
     sum(explained)>1+1e-8)
    stop("PCA explained-variance decomposition invalid")
  list(
    scores=scores,
    explained_variance_fraction=as.list(explained),
    measured_feature_count=ncol(X),
    n_components=as.integer(n_components),
    preprocessing=if(feature_standardization)
      "Centered and independently standardized EACH feature, label-blind"
    else "Centered each feature, original between-feature scales, label-blind",
    labels_supplied_to_pca=FALSE,
    caveat="PCA is an unsupervised single-omic comparator on the SAME study animals and author-preselected features. Its post-fit metadata association is exploratory, not a clinical test or proof of one-omic sufficiency."
  )
}
