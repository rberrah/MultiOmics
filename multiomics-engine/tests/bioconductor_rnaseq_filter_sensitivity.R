#!/usr/bin/env Rscript
# Label-blind, paired raw-count filtering + normalization sensitivity audit.
# Uses real reference DESeq2 / edgeR QL / limma-voom implementations.
# This benchmark does NOT automatically change the production analysis policy.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/rnaseq_stress_designs.R")
source("multiomics-engine/rnaseq_filter_sensitivity.R")
source("multiomics-engine/mixed_truth_reference_metrics.R")
for(pkg in c("DESeq2","edgeR","limma","statmod"))
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("Missing required real reference package: ",pkg)
raw <- Sys.getenv("MULTIOMICS_FILTER_SENSITIVITY_REPLICATES","8")
reps <- suppressWarnings(as.integer(raw))
if(length(raw)!=1L||!grepl("^[0-9]+$",raw)||
   !is.finite(reps)||reps<8L||reps>500L)
  stop("Use 8 to 500 independent synthetic cohorts per scenario.")
grid <- rnaseq_stress_grid()
chosen <- c("small_cohort_6v6","unequal_batch_12v28","sparse_signals_20v20")
grid <- grid[match(chosen,grid$scenario),,drop=FALSE]
stopifnot(nrow(grid)==3L,!anyNA(grid$scenario))
policies <- c("reference","strict_prespecified")
out <- Sys.getenv("MULTIOMICS_VALIDATION_OUTPUT","")
if(nzchar(out)) dir.create(out,recursive=TRUE,showWarnings=FALSE)
root <- tempfile("rnaseq_prefilter_real_R_")
dir.create(root)
records <- list()
normalizations <- list()
design <- list()
for(s in seq_len(nrow(grid))) {
  scenario <- grid$scenario[[s]]
  for(rep in seq_len(reps)) {
    # Exactly the same simulated subjects, seed, raw counts and truth serve
    # both comparison arms. Thinning is independent of group membership.
    original <- rnaseq_stress_generate(grid[s,,drop=FALSE],rep)
    study <- rnaseq_dilute_counts(original,rep)
    norm_vectors <- list()
    for(policy in policies) {
      filtered <- rnaseq_independent_filter(study$x,policy)
      if(filtered$retained_genes<60L)
        stop("Prespecified filtering left too few features to benchmark.")
      paths <- file.path(root,paste0(scenario,"_",rep,"_",policy,
                      c("_d","_e","_v")))
      d <- run_deseq2_counts(filtered$x,study$meta,paths[[1L]],
        design_formula=~batch+condition,
        contrast=c("condition","treated","control"))
      e <- run_edger_ql_counts(filtered$x,study$meta,paths[[2L]],
        design_formula=~batch+condition,coefficient="conditiontreated")
      v <- run_voom_counts(filtered$x,study$meta,paths[[3L]],
        design_formula=~batch+condition,coefficient="conditiontreated")
      tabs <- list(
        DESeq2=data.frame(feature=d$results$feature,
          q=d$results$padj,log2fc=d$results$log2FoldChange),
        edgeR_QL=data.frame(feature=e$results$feature,
          q=e$results$FDR,log2fc=e$results$logFC),
        limma_voom=data.frame(feature=v$results$feature,
          q=v$results[["adj.P.Val"]],log2fc=v$results$logFC))
      excluded <- sum(study$truth[filtered$discarded]!=0)
      for(method in names(tabs)) {
        metric <- mixed_truth_feature_metrics(tabs[[method]],study$truth,
          method,scenario,rep,study$seed)
        metric$filter_policy <- policy
        metric$external_retained_genes <- filtered$retained_genes
        metric$excluded_planted_genes <- excluded
        metric$eligible_for_nominal_calibration <- FALSE
        if(metric$tested_features<60L || metric$finite_q<40L)
          stop("Unusable BH family after filtering: ",
               scenario," / ",method," / ",policy)
        records[[length(records)+1L]] <- metric
      }
      # Normalization is refitted after filtering in the current adapters:
      # treat the pair comparison as a FULL-PIPELINE sensitivity, not a
      # pure hypothesis-filter effect with fixed normalization.
      norm_vectors[[policy]] <- list(
        DESeq2=stats::setNames(as.numeric(DESeq2::sizeFactors(d$model)),
                              colnames(d$model)),
        edgeR_QL=stats::setNames(
          e$dge$samples$lib.size*e$dge$samples$norm.factors,
          rownames(e$dge$samples)),
        limma_voom=stats::setNames(
          v$dge$samples$lib.size*v$dge$samples$norm.factors,
          rownames(v$dge$samples)))
      unlink(paths,recursive=TRUE)
      rm(d,e,v,tabs,filtered)
    }
    for(method in c("DESeq2","edgeR_QL","limma_voom")) {
      reference <- norm_vectors$reference[[method]]
      strict <- norm_vectors$strict_prespecified[[method]]
      if(length(reference)!=nrow(study$x)||
         !identical(names(reference),names(strict))||
         any(!is.finite(reference))||any(!is.finite(strict))||
         any(reference<=0)||any(strict<=0))
        stop("Invalid paired reference normalization factors: ",method)
      delta <- log2(strict)-log2(reference)
      centered <- delta-mean(delta)
      normalizations[[length(normalizations)+1L]] <- data.frame(
        scenario=scenario,replicate=rep,seed=study$seed,
        method=method,
        max_abs_centered_log2_effective_scale_shift=max(abs(centered)),
        sd_centered_log2_effective_scale_shift=stats::sd(centered),
        normalization_reestimated_in_each_arm=TRUE,
        stringsAsFactors=FALSE)
    }
    design[[length(design)+1L]] <- data.frame(
      scenario=scenario,replicate=rep,seed=study$seed,
      n_control=original$n_control,n_treated=original$n_treated,
      counts_diluted=study$counts_diluted,
      n_source_genes=ncol(study$x),stringsAsFactors=FALSE)
    if(rep%%10L==0L||rep==reps)
      cat(sprintf("FILTER SENSITIVITY %s study %d/%d complete\n",
          scenario,rep,reps))
  }
}
records <- do.call(rbind,records)
paired <- rnaseq_filter_pair_summary(records,reps=reps,
  bootstrap=1000L,seed=20261012L)
by_policy <- records
by_policy$scenario <- paste(records$scenario,records$filter_policy,sep="__")
per_policy_summary <- mixed_truth_study_summary(by_policy,min_replicates=8L,
  bootstrap_repeats=1000L,seed=20261012L)
norm <- do.call(rbind,normalizations)
cohorts <- do.call(rbind,design)
if(nzchar(out)) {
  utils::write.csv(records,
    file.path(out,"rnaseq_filter_sensitivity_per_study.csv"),row.names=FALSE)
  utils::write.csv(paired,
    file.path(out,"rnaseq_filter_sensitivity_paired.csv"),row.names=FALSE)
  utils::write.csv(per_policy_summary,
    file.path(out,"rnaseq_filter_sensitivity_by_policy.csv"),row.names=FALSE)
  utils::write.csv(norm,
    file.path(out,"rnaseq_filter_sensitivity_normalization_shift.csv"),
    row.names=FALSE)
  utils::write.csv(cohorts,
    file.path(out,"rnaseq_filter_sensitivity_studies.csv"),row.names=FALSE)
  writeLines(c(
    paste("github_sha:",Sys.getenv("GITHUB_SHA","not_recorded")),
    paste("independent_studies_per_scenario:",reps),
    "Paired comparison: same patients, batch, depth and planted log2FC.",
    "Label-blind thinning: 30% genes 0.025x, 30% 0.2x, 40% 1x.",
    "Strict independent prefilter: count >=20 in >=40% samples.",
    "Reference: no additional external prefilter; each method's own filter.",
    "All arms recompute their native normalization after prefiltering.",
    "Therefore the paired perturbation includes normalization changes.",
    "No post-hoc filtering/policy selection based on p, q or performance.",
    "Bootstrap FDR intervals are descriptive: all-zero study FDPs yield degenerate zero intervals.",
    "Simultaneous Hoeffding 95% FDR bounds require independent iid simulated cohorts.",
    "No experimental biological, clinical or nominal-FDR certification.",
    paste("utc_timestamp:",format(Sys.time(),tz="UTC",usetz=TRUE)),
    capture.output(sessionInfo())),
    file.path(out,"rnaseq_filter_sensitivity_environment.txt"))
}
print(paired,row.names=FALSE)
print(per_policy_summary,row.names=FALSE)
print(aggregate(norm$max_abs_centered_log2_effective_scale_shift,
   list(scenario=norm$scenario,method=norm$method),mean),row.names=FALSE)
unlink(root,recursive=TRUE)
cat("REAL RNA-SEQ PAIRED PREFILTER / NORMALIZATION SENSITIVITY: EXECUTED; not a method-ranking or clinical certificate\n")
