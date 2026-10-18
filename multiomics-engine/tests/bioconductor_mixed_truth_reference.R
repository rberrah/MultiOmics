#!/usr/bin/env Rscript
# Real mixed-truth negative-binomial simulation: 1%, 5%, 20% of features have
# prespecified nonzero treatment effects (positive and negative). Real DESeq2,
# edgeR QL and limma-voom fit the IDENTICAL datasets.
# Neither a single successful run nor the permissive regression gate below
# demonstrates 5% population FDR calibration or biological transferability.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/mixed_truth_reference_metrics.R")
source("multiomics-engine/reference_inference_audit.R")
for(pkg in c("DESeq2","edgeR","limma","statmod"))
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("Required REAL R package unavailable: ",pkg)

# A larger, predeclared research campaign uses fresh simulations rather
# than re-labelling the eight-study regression-screen seeds.
raw_offset <- Sys.getenv("MULTIOMICS_MIXED_REFERENCE_REPLICATE_OFFSET","0")
if(!grepl("^[0-9]+$",raw_offset))
  stop("MULTIOMICS_MIXED_REFERENCE_REPLICATE_OFFSET must be a non-negative integer.")
replicate_offset <- suppressWarnings(as.integer(raw_offset))
if(length(replicate_offset)!=1L || !is.finite(replicate_offset) ||
   replicate_offset<0L || replicate_offset>9000L)
  stop("Invalid mixed-truth replicate offset.")
raw_repeats <- Sys.getenv("MULTIOMICS_MIXED_REFERENCE_REPLICATES","12")
repeats <- suppressWarnings(as.integer(raw_repeats))
if(!grepl("^[0-9]+$",raw_repeats) || !is.finite(repeats) || repeats<8L)
  stop("At least eight independent mixed-truth datasets per scenario required.")
n_genes <- 400L
genes <- sprintf("G%04d",seq_len(n_genes))
base_expression <- exp(seq(log(70),log(500),length.out=n_genes))
dispersion <- seq(0.10,0.24,length.out=n_genes)
scenarios <- data.frame(
  scenario=c("one_pct","five_pct","twenty_pct"),
  n_alternatives=c(4L,20L,80L),replicates=rep(repeats,3L),
  planted_abs_log2fc=1.5,stringsAsFactors=FALSE)

make_study <- function(scenario_index,replicate) {
  seed <- 620000L+scenario_index*10000L+replicate
  set.seed(seed)
  ids <- sprintf("S%02d",seq_len(40L))
  group <- factor(rep(c("control","treated"),each=20L),
    levels=c("control","treated"))
  batch <- factor(rep(c("B1","B2"),times=20L),
    levels=c("B1","B2"))
  meta <- data.frame(subject_id=sprintf("P%02d",seq_len(40L)),
    condition=group,batch=batch,row.names=ids)
  n_signal <- scenarios$n_alternatives[[scenario_index]]
  chosen <- sample(seq_len(n_genes),size=n_signal,replace=FALSE)
  truth <- setNames(rep(0,n_genes),genes)
  truth[chosen] <- rep(c(1.5,-1.5),length.out=n_signal)
  stopifnot(sum(truth!=0)==n_signal,sum(truth>0)==sum(truth<0),
            !anyDuplicated(meta$subject_id),
            qr(stats::model.matrix(~batch+condition,meta))$rank==3L)
  depth <- exp(stats::runif(40L,log(.75),log(1.3)))
  batch_effect <- exp(seq(-.20,.20,length.out=n_genes))
  counts <- vapply(seq_len(40L),function(i) {
    mean_counts <- base_expression*depth[[i]]*
      if(batch[[i]]=="B2") batch_effect else 1
    if(group[[i]]=="treated") mean_counts <- mean_counts*2^unname(truth)
    stats::rnbinom(n_genes,mu=mean_counts,size=1/dispersion)
  },numeric(n_genes))
  rownames(counts) <- genes
  colnames(counts) <- ids
  list(x=t(counts),meta=meta,truth=truth,seed=seed)
}

out <- Sys.getenv("MULTIOMICS_VALIDATION_OUTPUT","")
if(nzchar(out)) dir.create(out,recursive=TRUE,showWarnings=FALSE)
root <- tempfile("multiomics_mixed_real_R_")
dir.create(root)
records <- list()
coverage_records <- list()
integrity_records <- list()

for(s in seq_len(nrow(scenarios))) {
  scenario <- scenarios$scenario[[s]]
  for(iteration in seq_len(repeats)) {
    rep <- as.integer(replicate_offset + iteration)
    study <- make_study(s,rep)
    prefix <- file.path(root,sprintf("%s_%03d",scenario,rep))
    d <- run_deseq2_counts(study$x,study$meta,paste0(prefix,"_d"),
      design_formula=~batch+condition,
      contrast=c("condition","treated","control"))
    e <- run_edger_ql_counts(study$x,study$meta,paste0(prefix,"_e"),
      design_formula=~batch+condition,coefficient="conditiontreated")
    v <- run_voom_counts(study$x,study$meta,paste0(prefix,"_v"),
      design_formula=~batch+condition,coefficient="conditiontreated")
    tables <- list(
      DESeq2=data.frame(feature=d$results$feature,q=d$results$padj,
                        log2fc=d$results$log2FoldChange),
      edgeR_QL=data.frame(feature=e$results$feature,q=e$results$FDR,
                         log2fc=e$results$logFC),
      limma_voom=data.frame(feature=v$results$feature,q=v$results[["adj.P.Val"]],
                           log2fc=v$results$logFC))
    native <- list(
      DESeq2=list(table=d$results,p="pvalue",q="padj"),
      edgeR_QL=list(table=e$results,p="PValue",q="FDR"),
      limma_voom=list(table=v$results,p="P.Value",q="adj.P.Val")
    )
    for(name in names(tables)) {
      entry <- native[[name]]
      integrity_records[[length(integrity_records)+1L]] <-
        audit_native_r_bh(entry$table,name,scenario,rep,
          feature_col="feature",p_col=entry$p,q_col=entry$q,
          expected_features=names(study$truth))
      record <- mixed_truth_feature_metrics(tables[[name]],study$truth,
        name,scenario,rep,study$seed)
      if(record$finite_q<300L ||
         record$n_true!=scenarios$n_alternatives[[s]])
        stop("Unexpected severe feature loss or truth mismatch in ",scenario," ",name)
      records[[length(records)+1L]] <- record
    }
    # Only DESeq2 exposes native per-feature Wald SE in the current production
    # adapter. Measure the unshrunk *approximate* 95% coverage, DO NOT imply
    # equivalent interval estimates exist for edgeR or limma-voom.
    dt <- d$results
    expected <- unname(study$truth[dt$feature])
    usable <- is.finite(dt$log2FoldChange) &
      is.finite(dt$lfcSE) & dt$lfcSE>0 & expected!=0
    coverage <- if(any(usable)) mean(
      abs(dt$log2FoldChange[usable]-expected[usable]) <=
        stats::qnorm(.975)*dt$lfcSE[usable]) else NA_real_
    coverage_records[[length(coverage_records)+1L]] <- data.frame(
      scenario=scenario,replicate=rep,seed=study$seed,
      method="DESeq2_Wald_normal_approx",ci_coverage=coverage,
      n_signal_estimates_with_se=sum(usable),
      n_planted_alternatives=sum(study$truth!=0),
      stringsAsFactors=FALSE)
    cat(sprintf("MIXED TRUTH %s rep %02d/%02d | DESeq2 FDR=%.3f power=%.3f\n",
        scenario,iteration,repeats,
        records[[length(records)-2L]]$false_discovery_proportion,
        records[[length(records)-2L]]$power))
  }
}
result <- do.call(rbind,records)
integrity <- do.call(rbind,integrity_records)
stopifnot(nrow(integrity)==nrow(result),
  all(integrity$native_bh_integrity_verified),
  all(!integrity$fdr_calibration_certified))
summary <- mixed_truth_study_summary(result,min_replicates=8L)
coverage <- do.call(rbind,coverage_records)
coverage_summary <- aggregate(
  coverage[,c("ci_coverage","n_signal_estimates_with_se")],
  by=list(scenario=coverage$scenario,method=coverage$method),
  FUN=function(x) mean(x,na.rm=TRUE))
if(nzchar(out)) {
  utils::write.csv(integrity,
    file.path(out,"bioconductor_mixed_truth_native_bh_integrity.csv"),
    row.names=FALSE)
  utils::write.csv(result,file.path(out,"bioconductor_mixed_truth_per_study.csv"),
                   row.names=FALSE)
  utils::write.csv(summary,file.path(out,"bioconductor_mixed_truth_summary.csv"),
                   row.names=FALSE)
  utils::write.csv(coverage,file.path(out,"bioconductor_mixed_truth_deseq2_ci.csv"),
                   row.names=FALSE)
  utils::write.csv(coverage_summary,
    file.path(out,"bioconductor_mixed_truth_deseq2_ci_summary.csv"),
    row.names=FALSE)
  utils::write.csv(scenarios,
    file.path(out,"bioconductor_mixed_truth_scenarios.csv"),row.names=FALSE)
  writeLines(c(
    paste("github_sha:",Sys.getenv("GITHUB_SHA","not_recorded")),
    paste("replicate_offset:",replicate_offset),
    paste("per_scenario_repetitions:",repeats),
    "prespecified_new_studies: distinct seeds from prior short screens when offset>0",
    "Design: 20 control / 20 treated, each independent biological subject",
    "Primary contrast: treated vs control with measured balanced batch",
    "Planted: 1%, 5%, 20% alternatives, plus and minus log2 FC 1.5",
    "FDR estimand: mean V/max(R,1) across independent simulations",
    "Power estimand: all planted alternatives, even if removed by filtering",
    "DESeq2 CI: unshrunk normal-approximation Wald, not cross-method coverage",
    "Independent BH replay: native raw p to finite q, preserving independent filtering",
    "BH integrity does NOT mean native model null p-values are calibrated",
    "Bootstrap: fixed seed, percentile resampling of whole datasets",
    "Simultaneous Hoeffding 95% bounds: independent studies, all scenario-method comparisons",
    "Not an external-validation or clinical-certification result",
    capture.output(sessionInfo())),
    file.path(out,"bioconductor_mixed_truth_environment.txt"))
}
print(summary,row.names=FALSE)
print(coverage_summary,row.names=FALSE)
unlink(root,recursive=TRUE)
# Large threshold is a software-regression alarm, not an alpha=5% test.
# Always archive negative calibration results and do not select a preferred
# method according to observed performance on this small simulated suite.
if(any(summary$estimated_fdr>0.35) || any(summary$mean_power<0.02))
  stop("Gross mixed-truth false-discovery inflation or loss of planted-signal recovery.")
cat("REAL BIOCONDUCTOR MIXED-TRUTH BENCHMARK: PASS (gross check only, NOT nominal FDR certification)\n")
