#!/usr/bin/env Rscript
# Real RNA-seq STRESS experiment, not a confirmatory/scientific certificate.
# PR: eight studies per prespecified stress scenario (software smoke only).
# Manual/monthly: 120+ studies per scenario to reduce Monte Carlo uncertainty.
# ALL reference methods fit the SAME synthetic independent cohorts.
source("multiomics-engine/advanced_methods.R")
source("multiomics-engine/rnaseq_stress_designs.R")
source("multiomics-engine/mixed_truth_reference_metrics.R")
for(pkg in c("DESeq2","edgeR","limma","statmod"))
  if(!requireNamespace(pkg,quietly=TRUE))
    stop("Missing actual Bioconductor R implementation: ",pkg)
raw <- Sys.getenv("MULTIOMICS_RNASEQ_STRESS_REPLICATES","8")
repeats <- suppressWarnings(as.integer(raw))
if(length(raw)!=1L || !grepl("^[0-9]+$",raw) ||
   !is.finite(repeats) || repeats<8L || repeats>1000L)
  stop("MULTIOMICS_RNASEQ_STRESS_REPLICATES must be an integer between 8 and 1000.")
grid <- rnaseq_stress_grid()
root <- tempfile("rnaseq_reference_stress_")
dir.create(root)
output <- Sys.getenv("MULTIOMICS_VALIDATION_OUTPUT","")
if(nzchar(output)) dir.create(output,recursive=TRUE,showWarnings=FALSE)
rows <- list()
coverage_rows <- list()
design_rows <- list()
for(i in seq_len(nrow(grid))) {
  scenario <- grid$scenario[[i]]
  for(rep in seq_len(repeats)) {
    study <- rnaseq_stress_generate(grid[i,,drop=FALSE],rep)
    prefix <- file.path(root,paste0(scenario,"_",sprintf("%04d",rep)))
    d_path <- paste0(prefix,"_DESeq2")
    e_path <- paste0(prefix,"_edgeR")
    v_path <- paste0(prefix,"_voom")
    d <- run_deseq2_counts(study$x,study$meta,d_path,
       design_formula=~batch+condition,
       contrast=c("condition","treated","control"))
    e <- run_edger_ql_counts(study$x,study$meta,e_path,
       design_formula=~batch+condition,
       coefficient="conditiontreated")
    v <- run_voom_counts(study$x,study$meta,v_path,
       design_formula=~batch+condition,
       coefficient="conditiontreated")
    tabs <- list(
      DESeq2=data.frame(feature=d$results$feature,
        q=d$results$padj,log2fc=d$results$log2FoldChange),
      edgeR_QL=data.frame(feature=e$results$feature,
        q=e$results$FDR,log2fc=e$results$logFC),
      limma_voom=data.frame(feature=v$results$feature,
        q=v$results[["adj.P.Val"]],log2fc=v$results$logFC))
    for(method in names(tabs)) {
      metric <- mixed_truth_feature_metrics(tabs[[method]],study$truth,
        method,scenario,rep,study$seed)
      if(metric$total_planted_features!=400L ||
         metric$n_true!=grid$n_alternatives[[i]] ||
         metric$tested_features<300L || metric$finite_q<250L)
        stop("Severe filtering, missing q-values or misaligned planted truth: ",
             scenario," / ",method)
      rows[[length(rows)+1L]] <- metric
    }
    # DESeq2's raw Wald SE provides only an unshrunk approximate interval.
    # Other methods' inferential uncertainty is deliberately not imputed.
    tab <- d$results
    truth <- unname(study$truth[tab$feature])
    eligible <- is.finite(tab$log2FoldChange) &
      is.finite(tab$lfcSE) & tab$lfcSE>0 & truth!=0
    cov <- if(any(eligible)) mean(
      abs(tab$log2FoldChange[eligible]-truth[eligible]) <=
        stats::qnorm(.975)*tab$lfcSE[eligible]) else NA_real_
    coverage_rows[[length(coverage_rows)+1L]] <- data.frame(
      scenario=scenario,replicate=rep,seed=study$seed,
      available_signal_se=sum(eligible),
      total_planted_signals=sum(study$truth!=0),
      raw_deseq2_approx_wald95_coverage=cov)
    design_rows[[length(design_rows)+1L]] <- data.frame(
      scenario=scenario,replicate=rep,seed=study$seed,
      n_control=study$n_control,n_treated=study$n_treated,
      observed_depth_ratio=study$observed_depth_ratio,
      design_condition_number=study$design_condition_number,
      planted_abs_log2fc=study$planted_abs_log2fc)
    # Prevent reference model RDS objects filling disk during long campaigns.
    unlink(c(d_path,e_path,v_path),recursive=TRUE)
    rm(d,e,v,study,tabs)
    if(rep%%10L==0L || rep==repeats)
      cat(sprintf("RNA STRESS %s replicate %d/%d complete\n",
                  scenario,rep,repeats))
  }
}
records <- do.call(rbind,rows)
summary <- mixed_truth_study_summary(records,min_replicates=8L,
  bootstrap_repeats=1000L,seed=20261011L)
coverage <- do.call(rbind,coverage_rows)
design <- do.call(rbind,design_rows)
# Flags are recorded for independent review; never reinterpret a successful
# execution as nominal FDR control, especially for 8-replicate PR smoke runs.
summary$mc_standard_error <- vapply(seq_len(nrow(summary)),function(i) {
  row <- summary[i,,drop=FALSE]
  values <- records$false_discovery_proportion[
    records$scenario==row$scenario & records$method==row$method]
  stats::sd(values)/sqrt(length(values))
},numeric(1))
# Bootstrap percentile intervals are descriptive only: all-zero FDP yields
# [0,0] even when just eight independent studies were simulated. Use a
# simultaneous distribution-free bound for nominal-q review decisions.
summary$nominal_fdr_compatible_with_conservative95 <-
  summary$conservative_fdr_95_low<=.05 &
  summary$conservative_fdr_95_high>=.05
summary$review_flag <- ifelse(summary$conservative_fdr_95_low>.05,
    "simulated_fdr_lower_bound_above_nominal_review",
    ifelse(summary$conservative_fdr_95_high<=.05,
      "scenario_specific_upper_bound_below_nominal_NOT_clinical",
      "fdr_calibration_unresolved_insufficient_independent_studies"))
summary$bootstrap_is_not_finite_sample_calibration <- TRUE
summary$scientific_certified <- FALSE
if(nzchar(output)) {
  utils::write.csv(records,file.path(output,"rnaseq_stress_per_study.csv"),
                   row.names=FALSE)
  utils::write.csv(summary,file.path(output,"rnaseq_stress_summary.csv"),
                   row.names=FALSE)
  utils::write.csv(coverage,file.path(output,"rnaseq_stress_deseq2_ci.csv"),
                   row.names=FALSE)
  utils::write.csv(design,file.path(output,"rnaseq_stress_design_diagnostics.csv"),
                   row.names=FALSE)
  utils::write.csv(grid,file.path(output,"rnaseq_stress_prespecified_grid.csv"),
                   row.names=FALSE)
  writeLines(c(
    paste("github_sha:",Sys.getenv("GITHUB_SHA","not_recorded")),
    paste("replicates_per_scenario:",repeats),
    "familywise_BH_each_method_each_synthetic_cohort: 0.05",
    "FDR: mean(V/max(R,1)) across independently generated cohorts",
    "Bootstrap 95% percentile intervals: DESCRIPTIVE ONLY; may be [0,0] after zero discoveries",
    "Conservative 95% simultaneous Hoeffding-Bonferroni bounds: independent iid synthetic cohorts",
    "The bound family corrects across every scenario-method cell in the same CSV",
    "Power denominator: all planted genes, including those removed by filtering",
    "Effect bias: only available retained planted features (selection caveat)",
    "DESeq2 coverage: unshrunk normal-approximation Wald only",
    "A bootstrap interval is NOT a confidence interval for external transport",
    "NO automatic nominal FDR validation or method ranking from these results",
    paste("utc_timestamp:",format(Sys.time(),tz="UTC",usetz=TRUE)),
    capture.output(sessionInfo())),
    file.path(output,"rnaseq_stress_environment.txt"))
}
print(summary,row.names=FALSE)
unlink(root,recursive=TRUE)
cat("REAL RNA-seq STRESS GRID: executed; inspect FDR/power flags and uncertainty; NOT a confirmatory certificate\n")
