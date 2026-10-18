# Prespecified diagnostic of original MultiOmics negative-binomial global-null
# GENERATOR. It mirrors the formula and scenario-specific design in the frozen
# 128-study benchmark, without importing/sourcing that executable test script.
# No change to default DESeq2, filtering or user-facing inference.
rnaseq_oracle_null_study <- function(scenario,replicate) {
  scenarios <- c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
  k <- match(scenario,scenarios)
  if(is.na(k)||length(replicate)!=1L||!is.finite(replicate)||
     replicate!=floor(replicate)||replicate<1L||replicate>9000L)
    stop("Invalid predeclared negative-binomial source study.")
  seed <- as.integer(780000L+k*10000L+as.integer(replicate))
  had <- exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  if(had) prior <- get(".Random.seed",envir=.GlobalEnv)
  on.exit(if(had)assign(".Random.seed",prior,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv),add=TRUE)
  set.seed(seed)
  genes <- sprintf("N%04d",seq_len(400L))
  means <- exp(seq(log(45),log(350),length.out=length(genes)))
  n_control <- if(k==3L)12L else 20L
  n_treated <- if(k==3L)28L else 20L
  n <- n_control+n_treated
  names <- sprintf("S%02d",seq_len(n))
  group <- factor(c(rep("control",n_control),rep("treated",n_treated)),
                  levels=c("control","treated"))
  if(k==3L) {
    batch <- c(sample(c(rep("B1",9L),rep("B2",3L))),
               sample(c(rep("B1",9L),rep("B2",19L))))
  } else {
    batch <- c(rep(c("B1","B2"),n_control/2L),
               rep(c("B1","B2"),n_treated/2L))
  }
  metadata <- data.frame(condition=group,
       batch=factor(batch,levels=c("B1","B2")),row.names=names)
  design <- stats::model.matrix(~batch+condition,metadata)
  stopifnot(qr(design)$rank==ncol(design),
            all(table(metadata$condition,metadata$batch)>0))
  high <- k==2L
  dispersion <- if(high) seq(.30,.70,length.out=length(genes))
                else seq(.08,.20,length.out=length(genes))
  depth <- if(high) exp(stats::runif(n,log(.3),log(3)))
           else stats::runif(n,.9,1.1)
  batch_effect <- exp(seq(-.50,.50,length.out=length(genes)))
  counts <- vapply(seq_len(n),function(i) {
    mu <- means*depth[[i]]*if(metadata$batch[[i]]=="B2")
      batch_effect else 1
    stats::rnbinom(length(genes),mu=mu,size=1/dispersion)
  },numeric(length(genes)))
  rownames(counts) <- genes
  colnames(counts) <- names
  if(!identical(dim(counts),c(400L,40L)) ||
     anyNA(counts)||any(counts<0)||any(counts!=floor(counts)))
    stop("Broken independent frozen null-generator input.")
  list(counts=t(counts),meta=metadata,depth=as.numeric(depth),
       dispersion=dispersion,seed=seed,scenario=scenario,
       replicate=as.integer(replicate),planted_condition_log2fc=0,
       design=design)
}

rnaseq_null_oracle_summary <- function(rows,replicates=32L) {
  needed <- c("scenario","replicate","seed","method","any_false_discovery",
              "finite_q","size_factor_log_rmse")
  if(!is.data.frame(rows)||!nrow(rows)||!all(needed%in%names(rows))||
     anyNA(rows[,needed])||
     any(!rows$method%in%c("DESeq2_default","DESeq2_oracle_true_depth"))||
     any(!rows$any_false_discovery%in%c(0L,1L)) ||
     any(!is.finite(rows$size_factor_log_rmse)))
    stop("Invalid paired native DESeq2 normalization comparison.")
  cells <- split(rows,list(rows$scenario,rows$method),drop=TRUE)
  scenarios <- c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
  if(length(cells)!=6L ||
     !setequal(unique(rows$scenario),scenarios) ||
     any(vapply(cells,nrow,integer(1))!=replicates))
    stop("Incomplete or post-selected normalization diagnostic.")
  for(s in scenarios) {
    z <- rows[rows$scenario==s,,drop=FALSE]
    a <- z[z$method=="DESeq2_default",]
    b <- z[z$method=="DESeq2_oracle_true_depth",]
    if(anyDuplicated(a$replicate)||anyDuplicated(b$replicate)||
       !setequal(a$replicate,b$replicate) ||
       !identical(a$seed[order(a$replicate)],b$seed[order(b$replicate)]))
      stop("Default and oracle factors were not tested on the same cohorts.")
  }
  out <- lapply(cells,function(z) {
    k <- sum(z$any_false_discovery)
    n <- nrow(z)
    ci <- stats::binom.test(k,n,conf.level=.95)$conf.int
    data.frame(scenario=z$scenario[[1L]],method=z$method[[1L]],
      independent_null_studies=n,
      studies_with_any_false_discovery=k,
      empirical_global_null_fdr=k/n,
      exact_binomial_ci_low=ci[[1L]],
      exact_binomial_ci_high=ci[[2L]],
      mean_size_factor_log_rmse=mean(z$size_factor_log_rmse),
      normalization_contributes_to_fdr_not_proven=TRUE,
      normalizations_are_synthetic_oracle_only=TRUE,
      scientific_calibration_certified=FALSE)
  })
  ans <- do.call(rbind,out)
  rownames(ans) <- NULL
  ans[order(ans$scenario,ans$method),,drop=FALSE]
}
