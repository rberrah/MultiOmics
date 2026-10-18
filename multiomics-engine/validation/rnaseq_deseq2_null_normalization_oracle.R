#!/usr/bin/env Rscript
# Investigate (not "correct") the inflation of DESeq2 global-null FDR from
# the existing MultiOmics negative-binomial simulator. This uses a new,
# predeclared seed block, EXACT SAME counts and study design in both arms:
#   (1) normal DESeq2 median-ratio size-factor estimation;
#   (2) synthetic-generator TRUE sequencing depths, known only in simulation.
# Both arms use untouched DESeq2 default Wald dispersion fitting, baseline
# low-count prefilter and independent filtering; NO production change.
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/reference_inference_audit.R")
if(!requireNamespace("DESeq2",quietly=TRUE))
  stop("Native pinned Bioconductor DESeq2 required; no model substitutes.")
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L || !nzchar(args[[1L]]))
  stop("Usage: Rscript multiomics-engine/validation/rnaseq_deseq2_null_normalization_oracle.R <empty output directory>")
out <- args[[1L]]
if(dir.exists(out) && length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Do not overwrite earlier scientific diagnostic evidence.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
scenario_names <- c("balanced_base","balanced_depth_dispersion",
                    "unequal_partial_batch")
repeats <- 32L
offset <- 2000L
records <- list()
for(scenario in scenario_names) {
  for(index in seq_len(repeats)) {
    rep <- offset+index
    study <- rnaseq_oracle_null_study(scenario,rep)
    stopifnot(all(study$counts>=0),
              all(study$counts==floor(study$counts)),
              identical(rownames(study$counts),rownames(study$meta)),
              study$planted_condition_log2fc==0)
    dds <- DESeq2::DESeqDataSetFromMatrix(
       countData=t(study$counts),colData=study$meta,
       design=~batch+condition)
    # Exactly the production adapter's prespecified simple count prefilter;
    # do not alter the gene family or condition labels between arms.
    keep <- rowSums(DESeq2::counts(dds)>=10L)>=max(2L,ceiling(ncol(dds)*.20))
    dds <- dds[keep,]
    if(nrow(dds)<300L)stop("Unexpected original 400-gene benchmark filtering.")
    sf_true <- study$depth/exp(mean(log(study$depth)))
    names(sf_true) <- rownames(study$meta)
    if(any(!is.finite(sf_true)) || any(sf_true<=0))
      stop("Simulation oracle sequencing-depth factors invalid.")
    for(arm in c("DESeq2_default","DESeq2_oracle_true_depth")) {
      local <- dds
      if(arm=="DESeq2_oracle_true_depth")
        DESeq2::sizeFactors(local) <- sf_true
      local <- DESeq2::DESeq(local,quiet=TRUE)
      factors <- as.numeric(DESeq2::sizeFactors(local))
      if(length(factors)!=length(sf_true) || any(!is.finite(factors))||
         any(factors<=0))
        stop("Native DESeq2 output lost or invalidated sequencing size factors.")
      if(arm=="DESeq2_oracle_true_depth" &&
         max(abs(log(factors/sf_true)))>1e-10)
        stop("DESeq2 oracle arm did not preserve the known simulation depths.")
      res <- DESeq2::results(local,contrast=c("condition","treated","control"),
                             alpha=.05)
      p <- as.numeric(res$pvalue)
      q <- as.numeric(res$padj)
      if(length(p)!=nrow(dds)||length(q)!=nrow(dds)||
         any(!is.na(p)&(p<0|p>1|!is.finite(p)))||
         any(!is.na(q)&(q<0|q>1|!is.finite(q)))||
         any(is.finite(q)&!is.finite(p))||
         sum(is.finite(q))<300L)
        stop("Native DESeq2 Wald null-result p-value structure invalid.")
      adj <- is.finite(q)
      deviation <- max(abs(independent_bh_adjust(p[adj])-q[adj]))
      if(!is.finite(deviation)||deviation>1e-8)
        stop("Independent BH q-replay fails under null-oracle normalization.")
      records[[length(records)+1L]] <- data.frame(
        scenario=scenario,replicate=as.integer(rep),seed=study$seed,
        method=arm,source_genes=length(study$dispersion),
        modeled_genes=nrow(dds),
        finite_p=sum(is.finite(p)),finite_q=sum(is.finite(q)),
        p_below_005=sum(is.finite(p)&p<=.05),
        discoveries_q05=sum(is.finite(q)&q<=.05),
        any_false_discovery=as.integer(any(is.finite(q)&q<=.05)),
        size_factor_log_rmse=sqrt(mean((log(factors)-log(sf_true))^2)),
        median_signed_size_factor_log_error=stats::median(log(factors/sf_true)),
        oracle_factors_from_simulation_only=arm=="DESeq2_oracle_true_depth",
        independent_native_BH_verified=TRUE,
        stringsAsFactors=FALSE)
    }
    cat("DESEQ2 ORACLE",scenario,index,"/",repeats,
        "paired original counts PASSED\n")
  }
}
rows <- do.call(rbind,records)
summary <- rnaseq_null_oracle_summary(rows,repeats)
if(nrow(rows)!=length(scenario_names)*repeats*2L || nrow(summary)!=6L)
  stop("Actual null-normalization diagnostic has missing cohorts.")
utils::write.csv(rows,file.path(out,"deseq2_normalization_oracle_per_study.csv"),
                 row.names=FALSE)
utils::write.csv(summary,
  file.path(out,"deseq2_normalization_oracle_global_null_summary.csv"),
  row.names=FALSE)
writeLines(c(
  "schema: multiomics.deseq2-null-size-factor-oracle-diagnostic.v1",
  paste("github_sha:",Sys.getenv("GITHUB_SHA","local_unverified")),
  "original_null_generator: same NB means, depths, gene-specific batch effects, and 3 scenario designs as PR43",
  "replicate_seed_offset: 2000",
  "independent_studies_per_scenario: 32",
  "same_original_count_matrix_in_both_arms: true",
  "default_reference: production DESeq2 median ratio / default dispersion Wald / BH",
  "oracle_only_diagnostic: known positive sample depth factors normalized to geometric mean one",
  "both_arms_preserve: prespecified count filter, gene families, DESeq2 Wald test and independent filtering",
  "direct_bh_validation: independently coded BH formula (not R p.adjust)",
  "paired_scientific_question: does source-depth normalization contribute to inflated null p-values?",
  "no_causal_attribution_if_still_inflated: correct_true_depth_is_only_one_diagnostic",
  "nominal_bh_q: .05",
  "sample_size_precision: 32 studies per scenario; underpowered to certify 5% FDR",
  "no_native_method_optimisation_or_production_defaults_changed: true",
  "clinical_and_universal_FDR_certified: false",
  capture.output(sessionInfo())
),file.path(out,"deseq2_normalization_oracle_method_manifest.txt"))
print(summary,row.names=FALSE)
cat("PASS: diagnostic produced authentic DESeq2 paired null-genotype controls; neither arm FDR certified.\n")
