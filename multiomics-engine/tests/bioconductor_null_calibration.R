#!/usr/bin/env Rscript
# Real reference-package RNA-seq GLOBAL-NULL regression screen.
#
# Three prespecified designs use negative-binomial counts under no biological
# treatment effect. In an all-null experiment BH discoveries are all false;
# P(any BH discovery) = FDR = FWER. The binomial intervals quantify simulation
# uncertainty (NOT biological confidence intervals and NOT a calibration proof).
#
# Executed DESeq2, edgeR QL and limma-voom must remain separate from the
# fast browser/Node simulation benchmarks, and must not certify user studies.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/null_reference_metrics.R")
source("multiomics-engine/reference_inference_audit.R")
for (pkg in c("DESeq2","edgeR","limma","statmod")) {
  if (!requireNamespace(pkg,quietly=TRUE))
    stop("Missing real reference R package: ",pkg)
}
as_replicates <- function(name, default) {
  raw <- Sys.getenv(name,as.character(default))
  n <- suppressWarnings(as.integer(raw))
  if (!grepl("^[0-9]+$",raw) || length(n)!=1L || !is.finite(n) || n<8L)
    stop(name,": at least eight integer replicates are required.")
  n
}
base_repeats <- as_replicates("MULTIOMICS_NULL_REPLICATES",12L)
stress_repeats <- as_replicates("MULTIOMICS_NULL_STRESS_REPLICATES",8L)
cases <- data.frame(
  scenario=c("balanced_base","balanced_depth_dispersion","unequal_partial_batch"),
  n_control=c(20L,20L,12L), n_treated=c(20L,20L,28L),
  replicates=c(base_repeats,stress_repeats,stress_repeats),
  stringsAsFactors=FALSE)
genes <- sprintf("N%04d",seq_len(400L))
base_expression <- exp(seq(log(45),log(350),length.out=length(genes)))
cases$mechanism <- c(
  "balanced groups, mild depth variation and balanced measured batch effect",
  "balanced groups, 10-fold depth variation and high negative-binomial dispersion",
  "12/28 groups, partial batch-group confounding and measured batch effect")
stopifnot(!anyDuplicated(cases$scenario))

make_null_study <- function(scenario,replicate,scenario_index) {
  seed <- as.integer(780000L+scenario_index*10000L+replicate)
  set.seed(seed)
  n_control <- cases$n_control[[scenario_index]]
  n_treated <- cases$n_treated[[scenario_index]]
  n <- n_control+n_treated
  ids <- sprintf("S%02d",seq_len(n))
  group <- factor(c(rep("control",n_control),rep("treated",n_treated)),
                  levels=c("control","treated"))
  if (scenario=="unequal_partial_batch") {
    # Deliberately associated batch + group, but both groups occur in
    # both batches. The treatment coefficient remains identifiable.
    batches <- c(sample(c(rep("B1",9L),rep("B2",3L))),
                 sample(c(rep("B1",9L),rep("B2",19L))))
  } else {
    batches <- c(rep(c("B1","B2"),n_control/2),
                 rep(c("B1","B2"),n_treated/2))
  }
  meta <- data.frame(condition=group,
    batch=factor(batches,levels=c("B1","B2")),row.names=ids)
  design <- stats::model.matrix(~batch+condition,meta)
  stopifnot(qr(design)$rank==ncol(design),
    all(table(meta$condition,meta$batch)>0))
  is_high <- identical(scenario,"balanced_depth_dispersion")
  dispersion <- if (is_high) seq(0.30,0.70,length.out=length(genes))
                else seq(0.08,0.20,length.out=length(genes))
  depth <- if (is_high) exp(stats::runif(n,log(0.3),log(3.0)))
           else stats::runif(n,0.9,1.1)
  # Feature-specific nuisance batch effects; treatment has EXACTLY zero
  # planted effect. This tests adjustment, not a simple label shuffle.
  batch_effect <- exp(seq(-0.50,0.50,length.out=length(genes)))
  counts <- vapply(seq_len(n),function(i) {
    mean_counts <- base_expression*depth[[i]]*
      if(meta$batch[[i]]=="B2") batch_effect else 1
    stats::rnbinom(length(genes),mu=mean_counts,size=1/dispersion)
  },numeric(length(genes)))
  rownames(counts) <- genes
  colnames(counts) <- ids
  list(counts=t(counts),meta=meta,seed=seed,
       observed_depth_ratio=max(depth)/min(depth),
       max_dispersion=max(dispersion),n_control=n_control,n_treated=n_treated)
}

records <- list()
bh_integrity_records <- list()
root <- tempfile("multiomics_null_reference_")
dir.create(root)
output <- Sys.getenv("MULTIOMICS_VALIDATION_OUTPUT","")
if (nzchar(output)) dir.create(output,recursive=TRUE,showWarnings=FALSE)

for (s in seq_len(nrow(cases))) {
  scenario <- cases$scenario[[s]]
  for (replicate in seq_len(cases$replicates[[s]])) {
    study <- make_null_study(scenario,replicate,s)
    prefix <- file.path(root,sprintf("%s_rep%03d",scenario,replicate))
    x <- study$counts
    meta <- study$meta
    d <- run_deseq2_counts(x,meta,paste0(prefix,"_d"),
      design_formula=~batch+condition,contrast=c("condition","treated","control"))
    v <- run_voom_counts(x,meta,paste0(prefix,"_v"),
      design_formula=~batch+condition,coefficient="conditiontreated")
    e <- run_edger_ql_counts(x,meta,paste0(prefix,"_e"),
      design_formula=~batch+condition,coefficient="conditiontreated")
    native <- list(
      DESeq2=list(data=d$results,p="pvalue",q="padj"),
      limma_voom=list(data=v$results,p="P.Value",q="adj.P.Val"),
      edgeR_QL=list(data=e$results,p="PValue",q="FDR")
    )
    for (method in names(native)) {
      entry <- native[[method]]
      integrity <- audit_native_r_bh(entry$data,method,scenario,replicate,
        feature_col="feature",p_col=entry$p,q_col=entry$q,
        expected_features=genes)
      bh_integrity_records[[length(bh_integrity_records)+1L]] <- integrity
      q <- as.numeric(entry$data[[entry$q]])
      if(length(q)!=length(genes) || any(is.finite(q)&(q<0 | q>1)))
        stop("Invalid adjusted p-values: ",scenario," / ",method)
      usable <- sum(is.finite(q))
      if (usable < 300L)
        stop("Most q-values are unavailable: ",scenario," / ",method)
      selected <- sum(is.finite(q)&q<=0.05)
      records[[length(records)+1L]] <- data.frame(
        scenario=scenario,replicate=replicate,seed=study$seed,
        method=method,n_control=study$n_control,n_treated=study$n_treated,
        depth_ratio=study$observed_depth_ratio,
        max_dispersion=study$max_dispersion,
        features_tested=length(q),finite_q=usable,
        discoveries_q05=selected,
        false_discovery_fraction=as.integer(selected>0),
        stringsAsFactors=FALSE)
    }
    cat(sprintf("GLOBAL NULL %s replicate %02d/%02d complete\n",
      scenario,replicate,cases$replicates[[s]]))
  }
}
results <- do.call(rbind,records)
integrity <- do.call(rbind,bh_integrity_records)
stopifnot(nrow(integrity)==nrow(results),
  all(integrity$native_bh_integrity_verified),
  all(!integrity$fdr_calibration_certified))
summary <- null_reference_summary(results,min_replicates=8L)
if (nzchar(output)) {
  utils::write.csv(integrity,
    file.path(output,"bioconductor_null_native_bh_integrity.csv"),
    row.names=FALSE)
  utils::write.csv(results,file.path(output,"bioconductor_null_replicates.csv"),
                   row.names=FALSE)
  utils::write.csv(summary,file.path(output,"bioconductor_null_summary.csv"),
                   row.names=FALSE)
  utils::write.csv(cases,file.path(output,"bioconductor_null_scenarios.csv"),
                   row.names=FALSE)
  writeLines(c(
    paste("github_sha:",Sys.getenv("GITHUB_SHA","not_recorded")),
    paste("alpha_BH:",0.05),
    paste("gross_regression_threshold:",0.35),
    paste("reanalysis_scope: reference R global-null simulation only"),
    "native_BH_audit: independent base R replay on finite adjusted-q features",
    "DESeq2_independent_filtering: raw p can remain finite when adjusted q excluded",
    "integrity_not_calibration: exact BH matching says nothing about null p-value uniformity",
    paste("date_UTC:",format(Sys.time(),tz="UTC",usetz=TRUE)),
    "interpretation: no simulated cohort gives scientific certification",
    capture.output(sessionInfo())),
    file.path(output,"bioconductor_null_environment.txt"))
}
print(summary,row.names=FALSE)
# Deliberately insensitive: catches egregious regressions, NOT a claim that
# nominal FDR is achieved. At n=12, even zero discoveries yields a 95% exact
# upper bound above 0.20.
gross <- null_reference_gross_failures(summary,threshold=0.35)
unlink(root,recursive=TRUE)
if(nrow(gross))
  stop("Possible gross global-null false-discovery inflation: ",
       paste(paste(gross$scenario,gross$method,sep="/"),collapse=", "))
cat("REAL REFERENCE GLOBAL-NULL STRESS BENCHMARK: PASS (limited screen, NOT 5% FDR certification)\n")
