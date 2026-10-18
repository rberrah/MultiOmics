#!/usr/bin/env Rscript
# Prespecified fixed-horizon precision decision. Never choose a new stopping
# point or a winning RNA-seq method from observed outcomes. Only independent
# synthetic DATASETS (not genes) count as Monte Carlo replications.
source("multiomics-engine/null_reference_metrics.R")
source("multiomics-engine/mixed_truth_reference_metrics.R")

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=1L || !nzchar(args[[1L]]))
  stop("Usage: Rscript multiomics-engine/validation/rnaseq_precision_decision.R OUTPUT_DIR")
out <- args[[1L]]
null_path <- file.path(out,"bioconductor_null_replicates.csv")
mixed_path <- file.path(out,"bioconductor_mixed_truth_per_study.csv")
if (!file.exists(null_path) || !file.exists(mixed_path))
  stop("Missing actual native reference result files; never substitute demo metrics.")

null <- utils::read.csv(null_path,stringsAsFactors=FALSE)
mixed <- utils::read.csv(mixed_path,stringsAsFactors=FALSE)
expected_null <- c("balanced_base","balanced_depth_dispersion",
                   "unequal_partial_batch")
expected_mixed <- c("one_pct","five_pct","twenty_pct")
expected_methods <- c("DESeq2","edgeR_QL","limma_voom")
check_design <- function(tab,scenarios,replicates,label) {
  if(!all(c("scenario","method","replicate","seed") %in% names(tab)) ||
     anyNA(tab[,c("scenario","method","replicate","seed")]) ||
     !setequal(unique(tab$scenario),scenarios) ||
     !setequal(unique(tab$method),expected_methods))
    stop("Incorrect original scenario/method manifest for ",label)
  cells <- split(tab,list(tab$scenario,tab$method),drop=TRUE)
  if(length(cells)!=9L ||
     any(vapply(cells,nrow,integer(1))!=replicates) ||
     any(vapply(cells,function(x) anyDuplicated(x$replicate)>0L ||
                anyDuplicated(x$seed)>0L,logical(1))) ||
     any(tab$replicate<=1000L) ||
     any(tab$replicate>1000L+replicates))
    stop("Missing, repeated, or reused simulation seeds for ",label)
  # Each method must have been applied to the SAME independent simulation
  # for this scenario/replicate. Cross-method rows are not independent.
  for(scenario in scenarios) {
    local <- tab[tab$scenario==scenario,,drop=FALSE]
    for(rep in seq.int(1001L,1000L+replicates)) {
      selected <- local[local$replicate==rep,,drop=FALSE]
      if(nrow(selected)!=3L || length(unique(selected$seed))!=1L)
        stop("R fits are not paired on exactly one same-source study.")
    }
  }
  invisible(TRUE)
}
check_design(null,expected_null,128L,"global null")
check_design(mixed,expected_mixed,64L,"mixed truth")

null_summary <- null_reference_summary(
  null,min_replicates=128L,confidence=0.95,
  simultaneous_family_alpha=0.05)
mixed_summary <- mixed_truth_study_summary(
  mixed,min_replicates=64L,bootstrap_repeats=1000L,seed=20261018L)
stopifnot(nrow(null_summary)==9L,nrow(mixed_summary)==9L,
  all(null_summary$replicates==128L),
  all(null_summary$simultaneous_family_cells==9L),
  all(null_summary$simultaneous_family_confidence==.95),
  all(!null_summary$scientific_certified),
  all(!mixed_summary$scientific_certified),
  all(!mixed_summary$scenario_calibration_upper_below_5pct |
      mixed_summary$conservative_fdr_95_high<=0.05))

utils::write.csv(null_summary,
  file.path(out,"rnaseq_precision_null_simultaneous_95.csv"),row.names=FALSE)
utils::write.csv(mixed_summary,
  file.path(out,"rnaseq_precision_mixed_95.csv"),row.names=FALSE)

display <- function(x) formatC(x,format="f",digits=4)
lines <- c(
  "# MultiOmics — RNA-seq frozen Monte Carlo precision audit",
  "",
  paste0("Source SHA: \`",Sys.getenv("GITHUB_SHA","local_unverified"),"\`"),
  "",
  "**This report is a computation audit, NOT proof of universal 5% FDR control.**",
  "",
  "- Global null: 128 new independent datasets per scenario, three methods;",
  "- Mixed truth: 64 new independent datasets per scenario, three methods;",
  "- Nominal BH threshold q=0.05; no post-hoc method selection;",
  "- Seeds offset by 1000 from short-screen fixtures, per-study outputs retained;",
  "- Confidence intervals under the global null are exact two-sided binomial",
  "  per cell and simultaneous 95% family-wise Clopper–Pearson–Bonferroni",
  "  across the **nine scenario/method comparisons**;",
  "- Mixed truth: replicate means of false-discovery proportion plus",
  "  95% simultaneous Hoeffding–Bonferroni bounds; very conservative;",
  "- Biological certification / real-world clinical performance: NOT assessed.",
  "",
  "## Global null: mean FDP = chance of ≥1 false call",
  "",
  "| Scenario | Method | False-call studies / 128 | FDR | Simultaneous 95% CI | Status |",
  "|---|---|---:|---:|---:|---|")
for(i in seq_len(nrow(null_summary))) {
  x <- null_summary[i,,drop=FALSE]
  status <- if(x$simultaneous_inflation_detected)
    "Inflation evidence in this simulated generator"
    else if(x$simultaneous_inconclusive)
      "Inconclusive vs nominal 0.05"
    else "Below nominal in this simulated generator"
  lines <- c(lines,sprintf("| %s | %s | %d / %d | %s | [%s, %s] | %s |",
    x$scenario,x$method,x$datasets_with_any_false_discovery,
    x$replicates,display(x$empirical_global_null_fdr),
    display(x$simultaneous_ci_low),display(x$simultaneous_ci_high),status))
}
lines <- c(lines,"","## Planted-truth FDR and recovery",
  "",
  "| Scenario | Method | Mean FDP | Simultaneous conservative 95% upper | Mean power |",
  "|---|---|---:|---:|---:|")
for(i in seq_len(nrow(mixed_summary))) {
  x <- mixed_summary[i,,drop=FALSE]
  lines<-c(lines,sprintf("| %s | %s | %s | %s | %s |",
    x$scenario,x$method,display(x$estimated_fdr),
    display(x$conservative_fdr_95_high),display(x$mean_power)))
}
lines <- c(lines,"",
  "## Scientific interpretation / remaining uncertainty",
  "",
  "The same datasets are reused across **methods**, but independent across",
  "replicate IDs within each prespecified scenario. Bonferroni coverage",
  "permits cross-method dependence; it does not validate the generator.",
  "",
  "A numerical BH replay passing only proves correct post-processing of",
  "the tested p-values. NB model specification, dispersion estimation,",
  "independent filtering, real biological batch effects, MNAR, paired",
  "studies and cross-platform transport still require separate qualification.",
  "",
  "Never change production significance thresholds or choose a model",
  "based on whichever of these nine cells has the lowest observed FDR.",
  "")
writeLines(lines,file.path(out,"RNASEQ_MONTE_CARLO_PRECISION_REPORT.md"))
writeLines(c(
  paste("source_sha:",Sys.getenv("GITHUB_SHA","local_unverified")),
  "null_replicates_per_scenario: 128",
  "mixed_replicates_per_scenario: 64",
  "replicate_seed_offset: 1000",
  "method_comparisons: DESeq2, edgeR_QL, limma_voom on the same generated datasets",
  "scenario_method_cells_per_family: 9",
  "confidence_family: 0.95",
  "stopping_rule: fixed sample, no data-adaptive extension",
  "nominal_bh_q: 0.05",
  "inferential_scope: fixed negative-binomial synthetic mechanisms only",
  "clinical_validity: false",
  capture.output(sessionInfo())
),file.path(out,"rnaseq_precision_run_manifest.txt"))
cat("MONTE CARLO PRECISION REPORT: generated 9 null + 9 mixed-truth cells; no FDR certification\\n")
