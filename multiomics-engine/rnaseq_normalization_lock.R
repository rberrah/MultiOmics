# Explicit RNA-seq scaling intervention for *benchmark-only* analyses.
# Native baseline = usual method-specific filter + native normalization.
# Strict frozen = strict external filter, native method's internal filter,
# but normalization fixed to the corresponding native baseline sample scales.
# Strict native = strict external filter + native normalization re-estimated.
# The path-dependent additive decomposition attributes all non-scale model
# changes to the first step and all downstream effects of changing sample
# scales within strict filtering to the second. It does NOT isolate a pure
# gene-filter-only causal effect.

rnaseq_reference_sample_scales <- function(method,fit) {
  scales <- switch(method,
    DESeq2=stats::setNames(
      as.numeric(DESeq2::sizeFactors(fit$model)),colnames(fit$model)),
    edgeR_QL=stats::setNames(
      fit$dge$samples$lib.size*fit$dge$samples$norm.factors,
      rownames(fit$dge$samples)),
    limma_voom=stats::setNames(
      fit$dge$samples$lib.size*fit$dge$samples$norm.factors,
      rownames(fit$dge$samples)),
    stop("Unrecognized RNA-seq method."))
  if(anyNA(scales)||any(!is.finite(scales))||any(scales<=0)||
     !length(scales)||is.null(names(scales))||
     anyDuplicated(names(scales)))
    stop("Native R fit had invalid sample scaling factors.")
  scales
}

rnaseq_locked_scale_fit <- function(method,counts,metadata,
                                    reference_scales) {
  method <- match.arg(method,c("DESeq2","edgeR_QL","limma_voom"))
  x <- as.matrix(counts)
  meta <- as.data.frame(metadata)
  if(!is.numeric(x)||!nrow(x)||!ncol(x)||
     is.null(rownames(x))||is.null(colnames(x))||
     !identical(rownames(x),rownames(meta))||
     anyDuplicated(rownames(x))||anyDuplicated(colnames(x))||
     anyNA(x)||any(!is.finite(x))||any(x<0)||any(x!=floor(x))||
     is.null(names(reference_scales))||
     !setequal(names(reference_scales),rownames(x))||
     anyNA(reference_scales)||any(!is.finite(reference_scales))||
     any(reference_scales<=0))
    stop("Scaling intervention requires aligned raw integer counts and reference sample scales.")
  design <- stats::model.matrix(~batch+condition,meta)
  if(qr(design)$rank!=ncol(design)||
     !("conditiontreated"%in%colnames(design)))
    stop("RNA-seq condition is not identified independently of batch.")
  locked <- unname(reference_scales[rownames(x)])
  names(locked) <- rownames(x)
  if(method=="DESeq2") {
    if(!requireNamespace("DESeq2",quietly=TRUE)) stop("DESeq2 not installed.")
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData=t(x),colData=meta,design=~batch+condition)
    keep <- rowSums(DESeq2::counts(dds)>=10L) >=
      max(2L,ceiling(ncol(dds)*.2))
    if(!any(keep)) stop("No DESeq2 genes survived internal filtering.")
    dds <- dds[keep,]
    DESeq2::sizeFactors(dds) <- locked[colnames(dds)]
    dds <- DESeq2::DESeq(dds,quiet=TRUE)
    res <- DESeq2::results(dds,
      contrast=c("condition","treated","control"),alpha=.05)
    tab <- data.frame(feature=rownames(res),
      q=as.numeric(res$padj),log2fc=as.numeric(res$log2FoldChange),
      stringsAsFactors=FALSE)
    effective <- stats::setNames(
      as.numeric(DESeq2::sizeFactors(dds)),colnames(dds))
  } else {
    if(!requireNamespace("edgeR",quietly=TRUE)||
       !requireNamespace("limma",quietly=TRUE)) stop("Required R reference engines not installed.")
    y <- edgeR::DGEList(counts=t(x))
    keep <- edgeR::filterByExpr(y,design=design)
    if(!any(keep)) stop("No edgeR genes survived internal filtering.")
    y <- y[keep,,keep.lib.sizes=FALSE]
    # lib.size is the post-filter library size. Fixed effective lib sizes
    # are imposed via offsets; do NOT recompute TMM under the frozen arm.
    y$samples$norm.factors <- locked[rownames(y$samples)]/
      y$samples$lib.size
    effective <- stats::setNames(
      y$samples$lib.size*y$samples$norm.factors,
      rownames(y$samples))
    if(method=="edgeR_QL") {
      y <- edgeR::estimateDisp(y,design=design,robust=TRUE)
      fit <- edgeR::glmQLFit(y,design=design,robust=TRUE)
      test <- edgeR::glmQLFTest(fit,coef="conditiontreated")
      tt <- edgeR::topTags(test,n=Inf,sort.by="PValue")$table
      tab <- data.frame(feature=rownames(tt),q=as.numeric(tt$FDR),
        log2fc=as.numeric(tt$logFC),stringsAsFactors=FALSE)
    } else {
      v <- limma::voom(y,design=design,plot=FALSE)
      fit <- limma::eBayes(limma::lmFit(v,design),
                           robust=TRUE,trend=FALSE)
      tt <- limma::topTable(fit,coef="conditiontreated",
                            number=Inf,sort.by="P")
      tab <- data.frame(feature=rownames(tt),
        q=as.numeric(tt[["adj.P.Val"]]),log2fc=as.numeric(tt$logFC),
        stringsAsFactors=FALSE)
    }
  }
  if(!identical(names(effective),names(locked))||
     !isTRUE(all.equal(unname(effective),unname(locked),
                        tolerance=1e-9)))
    stop("Model recomputed or changed the frozen normalization scales.")
  list(results=tab,effective_scales=effective,
       normalization_locked_verified=TRUE,
       internal_testing_retained_genes=nrow(tab))
}

rnaseq_scale_path_summary <- function(records,minimum=8L) {
  required <- c("scenario","method","replicate","seed","arm",
                "power","false_discovery_proportion","n_true")
  arms <- c("reference_native","strict_frozen",
            "strict_native")
  if(!is.data.frame(records)||!nrow(records)||
     !all(required%in%names(records))||anyNA(records[,required])||
     any(!records$arm%in%arms)||any(records$power<0|records$power>1)||
     any(records$false_discovery_proportion<0|
         records$false_discovery_proportion>1))
    stop("Malformed RNA-seq normalization decomposition dataset.")
  groups <- split(records,list(records$scenario,records$method),
                  drop=TRUE)
  result <- lapply(groups,function(group){
    cells <- lapply(arms,function(a)
      group[group$arm==a,,drop=FALSE])
    names(cells)<-arms
    n <- nrow(cells[[1L]])
    if(n<minimum||any(vapply(cells,nrow,integer(1))!=n))
      stop("Insufficient or unpaired studies in scale decomposition.")
    sorted <- lapply(cells,function(x){
      if(anyDuplicated(x$replicate)||anyDuplicated(x$seed))
        stop("Duplicate scientific replicates.")
      x[order(x$replicate),,drop=FALSE]
    })
    a <- sorted[[1L]];b <- sorted[[2L]];c <- sorted[[3L]]
    if(!identical(a$replicate,b$replicate)||
       !identical(a$replicate,c$replicate)||
       !identical(a$seed,b$seed)||!identical(a$seed,c$seed)||
       !identical(a$n_true,b$n_true)||!identical(a$n_true,c$n_true))
      stop("Decomposition compares different underlying studies.")
    row <- list(scenario=a$scenario[[1L]],method=a$method[[1L]],
                independent_datasets=n)
    for(metric in c("power","false_discovery_proportion")) {
      total <- c[[metric]]-a[[metric]]
      first <- b[[metric]]-a[[metric]]
      second<-c[[metric]]-b[[metric]]
      if(!isTRUE(all.equal(total,first+second,tolerance=1e-12)))
        stop("Path decomposition is not exactly additive.")
      key <- if(metric=="power") "power" else "fdr"
      row[[paste0(key,"_total_change")]]<-mean(total)
      row[[paste0(key,"_filter_and_refit_at_frozen_scales")]]<-mean(first)
      row[[paste0(key,"_renormalize_on_strict_filter")]]<-mean(second)
    }
    row$certified <- FALSE
    as.data.frame(row,stringsAsFactors=FALSE)
  })
  out <- do.call(rbind,result)
  rownames(out)<-NULL
  out[order(out$scenario,out$method),,drop=FALSE]
}


# Paired Monte Carlo uncertainty across INDEPENDENT simulated studies (NOT
# features/genes). Descriptive percentile bootstrap, not a statistical
# calibration or a confidence bound on the biological FDR of a user's cohort.
rnaseq_scale_path_paired_uncertainty <- function(records,minimum=8L,
                                                 bootstrap_repeats=2000L,
                                                 seed=20261015L) {
  if(length(bootstrap_repeats)!=1L || !is.finite(bootstrap_repeats) ||
     bootstrap_repeats!=floor(bootstrap_repeats) ||
     bootstrap_repeats<100L || bootstrap_repeats>50000L ||
     length(seed)!=1L || !is.finite(seed) || seed!=floor(seed) ||
     seed<0L || seed>2147483647L)
    stop("Require an explicit finite seed and 100..50000 bootstrap draws.")
  # Reuse exactly the reference paired-validation contract, including study
  # ID and planted-signal alignment, rather than accepting arbitrary rows.
  rnaseq_scale_path_summary(records,minimum=minimum)
  had_rng<-exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE)
  previous_rng<-if(had_rng)get(".Random.seed",envir=.GlobalEnv) else NULL
  on.exit({
    if(had_rng)assign(".Random.seed",previous_rng,envir=.GlobalEnv)
    else if(exists(".Random.seed",envir=.GlobalEnv,inherits=FALSE))
      rm(".Random.seed",envir=.GlobalEnv)
  },add=TRUE)
  set.seed(as.integer(seed))
  groups<-split(records,list(records$scenario,records$method),drop=TRUE)
  arms<-c("reference_native","strict_frozen","strict_native")
  contrasts<-list(
    total=c("strict_native","reference_native"),
    filter_and_refit_at_frozen_scales=c("strict_frozen","reference_native"),
    renormalize_on_strict_filter=c("strict_native","strict_frozen"))
  output<-list()
  for(group in groups) {
    aligned<-lapply(arms,function(arm) {
      x<-group[group$arm==arm,,drop=FALSE]
      x[order(x$replicate),,drop=FALSE]
    })
    names(aligned)<-arms
    n<-nrow(aligned[[1L]])
    # Each column is one bootstrap replicate of the paired DATASETS. Use
    # exactly the same study indexes for every arm and metric.
    index<-matrix(sample.int(n,size=n*bootstrap_repeats,replace=TRUE),
                  nrow=n,ncol=bootstrap_repeats)
    for(metric in c("power","false_discovery_proportion")) {
      for(label in names(contrasts)) {
        pair<-contrasts[[label]]
        delta<-aligned[[pair[[1L]]]][[metric]] -
               aligned[[pair[[2L]]]][[metric]]
        draws<-colMeans(matrix(delta[index],nrow=n,
                               ncol=bootstrap_repeats))
        ci<-as.numeric(stats::quantile(draws,c(.025,.975),
              names=FALSE,type=7))
        output[[length(output)+1L]]<-data.frame(
          scenario=aligned[[1L]]$scenario[[1L]],
          method=aligned[[1L]]$method[[1L]],
          metric=if(metric=="power") "power" else "empirical_mean_FDP",
          comparison=label,
          independent_studies=n,
          mean_paired_difference=mean(delta),
          percentile_95_lower=ci[[1L]],
          percentile_95_upper=ci[[2L]],
          bootstrap_repeats=as.integer(bootstrap_repeats),
          bootstrap_seed=as.integer(seed),
          interval_scope="paired_iid_simulated_studies_descriptive_only",
          clinical_certification=FALSE,
          nominal_fdr_certified=FALSE,
          stringsAsFactors=FALSE)
      }
    }
  }
  do.call(rbind,output)
}
