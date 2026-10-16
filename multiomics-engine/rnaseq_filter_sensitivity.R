# Explicit, frozen sensitivity experiment. This code does NOT silently change
# the default production RNA-seq prefilter or estimate a good cutoff from p/q.
# All comparisons are paired on exactly the same synthetic biological study.
rnaseq_dilute_counts <- function(study, replicate) {
  if(!is.list(study)||!all(c("x","meta","truth","seed")%in%names(study)))
    stop("A complete pre-generated study with truth and sample IDs is required.")
  x <- as.matrix(study$x)
  if(!is.numeric(x)||is.null(dimnames(x))||
     !identical(rownames(x),rownames(study$meta))||
     !identical(colnames(x),names(study$truth))||
     anyNA(x)||any(!is.finite(x))||any(x<0)||
     any(x!=floor(x))||anyDuplicated(colnames(x)) ||
     length(replicate)!=1L||is.na(replicate)||
     replicate<1||replicate!=floor(replicate))
    stop("Invalid source raw-counts matrix or sample/feature correspondence.")
  # Independent label-blind binomial thinning changes gene detectability,
  # without changing the within-gene treatment-vs-control fold ratio.
  # Three prespecified expression strata, assigned by gene ID (not by q):
  # low abundance 30%, moderate abundance 30%, unchanged 40%.
  n_features <- ncol(x)
  p <- rep(c(.025,.20,1.0),times=c(floor(.30*n_features),
           floor(.30*n_features),
           n_features-2L*floor(.30*n_features)))
  local_seed <- as.integer(6200000L+as.integer(study$seed))
  if(is.na(local_seed) || local_seed<1L) stop("Invalid deterministic seed.")
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) initial_rng <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(had_seed) assign(".Random.seed",initial_rng,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(local_seed)
  diluted <- matrix(stats::rbinom(length(x),size=as.integer(x),
    prob=rep(p,each=nrow(x))),nrow=nrow(x),dimnames=dimnames(x))
  list(x=diluted,meta=study$meta,truth=study$truth,
       seed=study$seed,dilution_seed=local_seed,
       expression_retention=p,
       counts_diluted=TRUE)
}

rnaseq_independent_filter <- function(sample_by_gene, policy) {
  x <- as.matrix(sample_by_gene)
  if(!is.numeric(x)||!nrow(x)||!ncol(x)||
     anyNA(x)||any(!is.finite(x))||any(x<0)||any(x!=floor(x))||
     anyDuplicated(colnames(x))||
     !policy%in%c("reference","strict_prespecified"))
    stop("RNA-seq filter requires valid counts and a declared policy.")
  keep <- if(policy=="reference") rep(TRUE,ncol(x)) else
    colSums(x>=20L)>=ceiling(nrow(x)*.40)
  if(!any(keep)) stop("The predeclared filter removes all genes.")
  list(x=x[,keep,drop=FALSE],
       kept=colnames(x)[keep],discarded=colnames(x)[!keep],
       policy=policy,total_genes=ncol(x),
       retained_genes=sum(keep),
       filtering_blind_to_condition=TRUE)
}

rnaseq_filter_pair_summary <- function(records, reps=8L,bootstrap=1000L,
                                       seed=20261012L) {
  needed<-c("scenario","method","filter_policy","replicate","seed",
            "power","false_discovery_proportion","total_planted_features",
            "n_true","selected_features","external_retained_genes",
            "excluded_planted_genes")
  if(!is.data.frame(records)||!nrow(records)||
     !all(needed%in%names(records))||anyNA(records[,needed])||
     any(!records$filter_policy%in%c("reference","strict_prespecified"))||
     any(records$power<0|records$power>1)||
     any(records$false_discovery_proportion<0|
         records$false_discovery_proportion>1)||
     any(records$excluded_planted_genes<0|
         records$excluded_planted_genes>records$n_true)||
     length(bootstrap)!=1L||!is.finite(bootstrap)||bootstrap<100L)
    stop("Invalid paired preprocessing records.")
  groups <- split(records,list(records$scenario,records$method),drop=TRUE)
  had_seed <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had_seed) initial_rng <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(had_seed) assign(".Random.seed",initial_rng,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  out <- lapply(groups,function(group) {
    a<-group[group$filter_policy=="reference",,drop=FALSE]
    b<-group[group$filter_policy=="strict_prespecified",,drop=FALSE]
    if(nrow(a)<reps||nrow(a)!=nrow(b)||
       anyDuplicated(a$replicate)||anyDuplicated(b$replicate)||
       anyDuplicated(a$seed)||anyDuplicated(b$seed)||
       !setequal(a$replicate,b$replicate))
      stop("Unpaired or duplicated simulated studies.")
    a<-a[order(a$replicate),,drop=FALSE]
    b<-b[order(b$replicate),,drop=FALSE]
    if(!identical(a$seed,b$seed)||
       !identical(a$n_true,b$n_true)||
       !identical(a$total_planted_features,b$total_planted_features))
      stop("Pairing changed the true biological dataset or seed.")
    dp<-b$power-a$power
    df<-b$false_discovery_proportion-a$false_discovery_proportion
    n<-length(dp)
    draws<-replicate(bootstrap,sample.int(n,size=n,replace=TRUE),
                     simplify=FALSE)
    power_ci<-stats::quantile(vapply(draws,function(k) mean(dp[k]),
      numeric(1)),c(.025,.975),names=FALSE)
    fdr_ci<-stats::quantile(vapply(draws,function(k) mean(df[k]),
      numeric(1)),c(.025,.975),names=FALSE)
    data.frame(scenario=a$scenario[1L],method=a$method[1L],
       independent_paired_datasets=n,
       strict_minus_reference_power=mean(dp),
       paired_power_delta_ci95_low=power_ci[[1L]],
       paired_power_delta_ci95_high=power_ci[[2L]],
       strict_minus_reference_fdr=mean(df),
       paired_fdr_delta_ci95_low=fdr_ci[[1L]],
       paired_fdr_delta_ci95_high=fdr_ci[[2L]],
       mean_genes_removed_by_strict=mean(
           a$total_planted_features-b$external_retained_genes),
       mean_true_signals_removed_by_strict=mean(b$excluded_planted_genes),
       effect_direction="strict_minus_reference",
       certified=FALSE,stringsAsFactors=FALSE)
  })
  result<-do.call(rbind,out)
  rownames(result)<-NULL
  result[order(result$scenario,result$method),,drop=FALSE]
}
