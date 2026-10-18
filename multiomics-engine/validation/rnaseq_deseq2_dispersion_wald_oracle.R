#!/usr/bin/env Rscript
# Independent frozen-generator diagnostic: what happens to the actual DESeq2
# Wald global-null p-values when *true* NB dispersion is supplied?
# The third arm also supplies true sample-depth factors. These truths are
# known ONLY for synthetic counts. No production model or p-values are changed.
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/rnaseq_deseq2_dispersion_wald_metrics.R")
source("multiomics-engine/reference_inference_audit.R")
if(!requireNamespace("DESeq2",quietly=TRUE))
  stop("Real pinned DESeq2 required: cannot substitute a mock.")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L||!nzchar(args[[1L]]))
  stop("Usage: Rscript multiomics-engine/validation/rnaseq_deseq2_dispersion_wald_oracle.R EMPTY_OUTPUT_DIR")
out<-args[[1L]]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Refuse to overwrite a prior scientific benchmark run.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
scenarios<-c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
methods<-c("native_default","oracle_true_dispersion",
           "oracle_true_depth_and_dispersion")
repetitions<-32L
offset<-3000L
rows<-list()
for(scenario in scenarios) for(iteration in seq_len(repetitions)) {
  replicate<-offset+iteration
  study<-rnaseq_oracle_null_study(scenario,replicate)
  dds<-DESeq2::DESeqDataSetFromMatrix(
    countData=t(study$counts),colData=study$meta,
    design=~batch+condition)
  # Identical to production run_deseq2_counts, across ALL arms.
  keep<-rowSums(DESeq2::counts(dds)>=10L)>=max(2L,ceiling(ncol(dds)*.2))
  dds<-dds[keep,]
  if(nrow(dds)<300L||anyDuplicated(rownames(dds)))
    stop("Unexpected gene filtering or missing simulation features.")
  truth<-study$dispersion[match(rownames(dds),colnames(study$counts))]
  if(length(truth)!=nrow(dds)||any(!is.finite(truth))||any(truth<=0))
    stop("Native R feature order differs from the simulated true dispersion.")
  sf_true<-study$depth/exp(mean(log(study$depth)))
  names(sf_true)<-rownames(study$meta)
  for(method in methods) {
    fit<-dds
    if(method=="native_default") {
      fit<-DESeq2::DESeq(fit,quiet=TRUE)
    } else {
      # Supported native API: nbinomWaldTest uses previously supplied
      # size factors and per-gene negative-binomial dispersion estimates.
      # Only the nuisance inputs change; original counts and test stay fixed.
      if(method=="oracle_true_dispersion")
        fit<-DESeq2::estimateSizeFactors(fit)
      else DESeq2::sizeFactors(fit)<-sf_true
      DESeq2::dispersions(fit)<-truth
      if(max(abs(log(DESeq2::dispersions(fit)/truth)))>1e-12)
        stop("Truth-dispersion oracle setter changed supplied NB dispersion.")
      fit<-DESeq2::nbinomWaldTest(fit,betaPrior=FALSE,quiet=TRUE)
    }
    est<-as.numeric(DESeq2::dispersions(fit))
    factor<-as.numeric(DESeq2::sizeFactors(fit))
    if(length(est)!=length(truth)||any(!is.finite(est))||any(est<=0)||
       length(factor)!=length(sf_true)||any(!is.finite(factor))||
       any(factor<=0))
      stop("Native DESeq2 size factors or dispersion became nonfinite.")
    if(method!="native_default"&&
       max(abs(log(est/truth)))>1e-10)
      stop("True-dispersion arm no longer uses the original NB truth.")
    if(method=="oracle_true_depth_and_dispersion"&&
       max(abs(log(factor/sf_true)))>1e-10)
      stop("True-depth plus dispersion arm lost the known simulated depths.")
    res<-DESeq2::results(fit,contrast=c("condition","treated","control"),
                        alpha=.05)
    p<-as.numeric(res$pvalue)
    q<-as.numeric(res$padj)
    stat<-as.numeric(res$stat)
    if(length(p)!=nrow(dds)||length(q)!=nrow(dds)||
       length(stat)!=nrow(dds)||
       any(!is.na(p)&(!is.finite(p)|p<0|p>1))||
       any(!is.na(q)&(!is.finite(q)|q<0|q>1))||
       any(is.finite(q)&!is.finite(p))||
       sum(is.finite(q))<300L)
      stop("True-dispersion native Wald results have invalid p/q families.")
    active<-is.finite(q)
    replay<-independent_bh_adjust(p[active])
    discrepancy<-max(abs(replay-q[active]))
    if(!is.finite(discrepancy)||discrepancy>1e-8)
      stop("Independent BH integrity failed after oracle dispersion modification.")
    if(sum(is.finite(stat))<300L)
      stop("Too few finite true-null DESeq2 Wald coefficients.")
    rows[[length(rows)+1L]]<-data.frame(
      scenario=scenario,method=method,
      replicate=as.integer(replicate),seed=as.integer(study$seed),
      n_genes_initial=400L,n_genes_modelled=nrow(dds),
      finite_p=sum(is.finite(p)),finite_q=sum(active),
      raw_p_below_001=sum(is.finite(p)&p<=.001),
      raw_p_below_01=sum(is.finite(p)&p<=.01),
      raw_p_below_05=sum(is.finite(p)&p<=.05),
      bh_discoveries=sum(active&q<=.05),
      fdr_any_discovery=as.integer(any(active&q<=.05)),
      wald_stat_mean=mean(stat[is.finite(stat)]),
      wald_stat_squared_mean=mean(stat[is.finite(stat)]^2),
      dispersion_log_rmse=sqrt(mean(log(est/truth)^2)),
      dispersion_log_median_ratio=stats::median(log(est/truth)),
      depth_factor_log_rmse=sqrt(mean(log(factor/sf_true)^2)),
      native_bh_integrity_verified=TRUE,
      original_true_condition_log2fc=0,
      synthetic_truth_is_known_only_in_simulation=TRUE,
      scientific_calibration_certified=FALSE,
      stringsAsFactors=FALSE)
  }
  cat("NB DISPERSION/WALD ORACLE",scenario,iteration,"/",repetitions,
      "three paired native fits complete\n",flush=TRUE)
}
observations<-do.call(rbind,rows)
summary<-rnaseq_dispersion_wald_summary(observations,
                                      replicates=repetitions,seed_offset=offset)
if(nrow(observations)!=288L||nrow(summary)!=9L)
  stop("Not all predeclared 96×3 fits completed.")
utils::write.csv(observations,
  file.path(out,"deseq2_null_dispersion_oracle_per_study.csv"),row.names=FALSE)
utils::write.csv(summary,
  file.path(out,"deseq2_null_dispersion_oracle_summary.csv"),row.names=FALSE)
writeLines(c(
  "schema: multiomics.deseq2-null-dispersion-wald-diagnostic.v1",
  paste("github_sha:",Sys.getenv("GITHUB_SHA","local_unverified")),
  "original_frozen_negative_binomial_generator: exact reused parameters",
  "fresh_replicates: 3001 through 3032 (32 per scenario)",
  "identical_counts_batch_labels_gene_filters_within_each_triplet: true",
  "source_400_genes_and_40_subjects: true",
  "native_DESeq2_Wald_fits: 288",
  "native_default: DESeq2::DESeq with estimated depth/dispersion",
  "oracle_true_dispersion: DESeq2::estimateSizeFactors, exact supplied NB alpha, DESeq2::nbinomWaldTest",
  "oracle_true_depth_and_dispersion: known original depth plus exact supplied NB alpha, nbinomWaldTest",
  "independent_BH_integrity_checked_with_no_p_adjust_call: true",
  "per_dataset_global_null_FDR_equals_FWER: true",
  "exact_CIs: individual 95pct; no simultaneous calibration claim",
  "normal_null_Wald_stat_squared_expectation_reference: 1",
  "study-level_rare_event_precision: only 32 independent studies per scenario",
  "original_batch_model: ~ batch + condition",
  "production_DESeq2_settings_changed: false",
  "universal_FDR_or_molecular_PK_certification: false",
  capture.output(sessionInfo())
),file.path(out,"deseq2_null_dispersion_oracle_manifest.txt"))
print(summary,row.names=FALSE)
cat("PASS real DESeq2 paired nuisance-dispersion diagnostic; observed FDR estimates are NOT automatically calibrated.\n")
