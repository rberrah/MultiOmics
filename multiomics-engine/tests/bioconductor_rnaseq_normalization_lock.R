#!/usr/bin/env Rscript
# Real Bioconductor intervention: conditional path decomposition of a strict
# external prefilter's total power/FDR change. Not a causal filter-only test:
# dispersions, residual weights, BH family and downstream models are refitted.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/rnaseq_stress_designs.R")
source("multiomics-engine/rnaseq_filter_sensitivity.R")
source("multiomics-engine/rnaseq_normalization_lock.R")
source("multiomics-engine/mixed_truth_reference_metrics.R")
for(pkg in c("DESeq2","edgeR","limma","statmod"))
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("Required real reference R package not installed: ",pkg)
raw<-Sys.getenv("MULTIOMICS_NORM_LOCK_REPLICATES","8")
reps<-suppressWarnings(as.integer(raw))
if(length(raw)!=1L||!grepl("^[0-9]+$",raw)||
   !is.finite(reps)||reps<8L||reps>500L)
  stop("Set 8 to 500 independent synthetic biological datasets.")
grid<-rnaseq_stress_grid()
requested<-c("small_cohort_6v6","unequal_batch_12v28",
             "sparse_signals_20v20")
grid<-grid[match(requested,grid$scenario),,drop=FALSE]
if(nrow(grid)!=3L||anyNA(grid$scenario)) stop("Frozen scenarios missing.")
root<-tempfile("rnaseq_norm_lock_");dir.create(root)
out<-Sys.getenv("MULTIOMICS_VALIDATION_OUTPUT","")
if(nzchar(out))dir.create(out,recursive=TRUE,showWarnings=FALSE)
rows<-list()
scales<-list()
features<-list()
methods<-c("DESeq2","edgeR_QL","limma_voom")
for(i in seq_len(nrow(grid))) {
  scenario<-grid$scenario[[i]]
  for(rep in seq_len(reps)) {
    original<-rnaseq_stress_generate(grid[i,,drop=FALSE],rep)
    study<-rnaseq_dilute_counts(original,rep)
    strict<-rnaseq_independent_filter(study$x,"strict_prespecified")
    if(strict$retained_genes<60L)
      stop("Prespecified strict external filter leaves fewer than 60 genes.")
    for(method in methods) {
      prefix<-file.path(root,paste0(scenario,"_",rep,"_",method))
      fit_native<-function(x,tag) {
        target<-paste0(prefix,"_",tag)
        fitted<-switch(method,
          DESeq2=run_deseq2_counts(x,study$meta,target,
              design_formula=~batch+condition,
              contrast=c("condition","treated","control")),
          edgeR_QL=run_edger_ql_counts(x,study$meta,target,
              design_formula=~batch+condition,
              coefficient="conditiontreated"),
          limma_voom=run_voom_counts(x,study$meta,target,
              design_formula=~batch+condition,
              coefficient="conditiontreated"))
        tab<-switch(method,
          DESeq2=data.frame(feature=fitted$results$feature,
                q=fitted$results$padj,
                log2fc=fitted$results$log2FoldChange),
          edgeR_QL=data.frame(feature=fitted$results$feature,
                q=fitted$results$FDR,
                log2fc=fitted$results$logFC),
          limma_voom=data.frame(feature=fitted$results$feature,
                q=fitted$results[["adj.P.Val"]],
                log2fc=fitted$results$logFC))
        scale<-rnaseq_reference_sample_scales(method,fitted)
        unlink(target,recursive=TRUE)
        list(tab=tab,scale=scale)
      }
      base<-fit_native(study$x,"baseline_native")
      native<-fit_native(strict$x,"strict_native")
      fixed<-rnaseq_locked_scale_fit(method,strict$x,study$meta,
                                    base$scale)
      if(!isTRUE(fixed$normalization_locked_verified)||
         !identical(names(fixed$effective_scales),
                    names(base$scale)))
        stop("Failed to verify frozen normalizations.")
      for(arm in c("reference_native","strict_frozen","strict_native")) {
        tab<-switch(arm,
          reference_native=base$tab,
          strict_frozen=fixed$results,
          strict_native=native$tab)
        m<-mixed_truth_feature_metrics(tab,study$truth,
             method,scenario,rep,study$seed)
        m$arm<-arm
        m$external_retained_genes<-
          if(arm=="reference_native") ncol(study$x) else ncol(strict$x)
        m$excluded_planted_genes<-
          if(arm=="reference_native") 0L else
            sum(study$truth[strict$discarded]!=0)
        m$nominal_fdr_certified<-FALSE
        if(m$tested_features<60L||m$finite_q<40L)
          stop("Severe gene loss under normalization intervention.")
        rows[[length(rows)+1L]]<-m
      }
      delta<-log2(native$scale[rownames(study$x)])-
             log2(base$scale[rownames(study$x)])
      if(any(!is.finite(delta)))stop("Invalid differential normalization.")
      delta<-delta-mean(delta)
      scales[[length(scales)+1L]]<-data.frame(
          scenario=scenario,replicate=rep,method=method,
          seed=study$seed,
          mean_abs_centered_log2_scale_shift=mean(abs(delta)),
          max_abs_centered_log2_scale_shift=max(abs(delta)),
          frozen_scales_verified=TRUE)
      features[[length(features)+1L]]<-data.frame(
          scenario=scenario,replicate=rep,method=method,
          n_original=ncol(study$x),
          n_external_strict=ncol(strict$x),
          n_planted=sum(study$truth!=0),
          n_planted_removed=sum(study$truth[strict$discarded]!=0),
          stringsAsFactors=FALSE)
      rm(base,native,fixed)
    }
    cat(sprintf("NORM LOCK %s %d/%d completed\n",scenario,rep,reps))
  }
}
records<-do.call(rbind,rows)
summary<-rnaseq_scale_path_summary(records,minimum=8L)
paired_uncertainty<-rnaseq_scale_path_paired_uncertainty(records,minimum=8L,
                              bootstrap_repeats=2000L,seed=20261015L)
scale_data<-do.call(rbind,scales)
feature_data<-do.call(rbind,features)
if(any(!is.finite(summary$power_total_change))||
   any(!is.finite(summary$fdr_total_change))||
   nrow(paired_uncertainty)!=nrow(summary)*6L||
   any(!is.finite(paired_uncertainty$mean_paired_difference))||
   any(paired_uncertainty$percentile_95_lower>paired_uncertainty$percentile_95_upper)||
   any(paired_uncertainty$clinical_certification)||
   any(paired_uncertainty$nominal_fdr_certified)||
   any(!summary$certified==FALSE))
  stop("Missing path-decomposition results.")
if(nzchar(out)) {
  utils::write.csv(records,file.path(out,"rnaseq_norm_lock_per_study.csv"),
                   row.names=FALSE)
  utils::write.csv(summary,file.path(out,"rnaseq_norm_lock_decomposition.csv"),
                   row.names=FALSE)
  utils::write.csv(paired_uncertainty,
      file.path(out,"rnaseq_norm_lock_paired_uncertainty.csv"),
      row.names=FALSE)
  utils::write.csv(scale_data,file.path(out,"rnaseq_norm_lock_scale_shift.csv"),
                   row.names=FALSE)
  utils::write.csv(feature_data,file.path(out,"rnaseq_norm_lock_features.csv"),
                   row.names=FALSE)
  writeLines(c(
    paste("github_sha:",Sys.getenv("GITHUB_SHA","not_recorded")),
    paste("independent_studies_per_scenario:",reps),
    "Same simulated subjects, count matrices, batch and planted signals per triple.",
    "Native: baseline original internal filters and its own sample scales.",
    "Frozen: strict external filter and all internal native filters; baseline scales locked.",
    "Strict native: strict external filter and its own recalculated sample scales.",
    "TOTAL = STRICT NATIVE - REFERENCE NATIVE",
    "PATH FILTER/REFIT = STRICT FROZEN - REFERENCE NATIVE",
    "PATH RENORMALIZATION = STRICT NATIVE - STRICT FROZEN",
    "Both path terms include downstream model refitting; not a pure causal intervention.",
    "Do not choose filtering/scaling on observed q-values or benchmark power.",
    "Decomposition describes paired empirical FDP; zero FDP does not establish FDR control.",
    "Bootstrap: paired independent simulated STUDIES (not genes), 2000 draws, percentile 95% CI.",
    "Bootstrap CI is descriptive, synthetic-generator-specific, and not evidence of nominal FDR control.",
    "No external validation; no nominal FDR/clinical certification.",
    paste("utc_timestamp:",format(Sys.time(),tz="UTC",usetz=TRUE)),
    capture.output(sessionInfo())),
    file.path(out,"rnaseq_norm_lock_environment.txt"))
}
print(summary,row.names=FALSE)
print(paired_uncertainty,row.names=FALSE)
unlink(root,recursive=TRUE)
cat("REAL BIOCONDUCTOR FROZEN-SCALE PATH DECOMPOSITION: PASS (execution only; no FDR certification)\n")
