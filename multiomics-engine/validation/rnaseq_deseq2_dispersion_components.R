#!/usr/bin/env Rscript
# Four-way native DESeq2 *dispersion-component* ablation on identical
# synthetic global-null counts: gene-wise Cox-Reid estimate, fitted
# dispersion-mean trend, final MAP shrinkage and simulator true alpha.
#
# We intentionally disable Cook's filter in ALL arms in results() so we can
# isolate p/Wald/dispersion effects. We do NOT alter production DESeq2.
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/rnaseq_deseq2_dispersion_component_metrics.R")
source("multiomics-engine/reference_inference_audit.R")
if(!requireNamespace("DESeq2",quietly=TRUE) ||
   !requireNamespace("S4Vectors",quietly=TRUE))
  stop("Exact native DESeq2/Bioc packages required (no mock fits).")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L||!nzchar(args[[1L]]))
  stop("Usage: Rscript rnaseq_deseq2_dispersion_components.R <empty evidence dir>")
out<-args[[1L]]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Do not overwrite prior scientific evidence.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
scenarios<-c("balanced_base","balanced_depth_dispersion",
             "unequal_partial_batch")
arms<-c("map_final","gene_wise_CR","fitted_trend","true_nb_alpha")
reps<-32L
first_rep<-4001L
records<-vector("list",length(scenarios)*reps*length(arms))
position<-0L
for(scenario in scenarios) for(idx in seq_len(reps)) {
  id<-first_rep+idx-1L
  source_study<-rnaseq_oracle_null_study(scenario,id)
  counts<-source_study$counts
  dds<-DESeq2::DESeqDataSetFromMatrix(
    countData=t(counts),colData=source_study$meta,
    design=~batch+condition)
  # Frozen original production-adapter low-count rule. NO source-gene
  # removal or reordering between the four dispersion arms.
  keep<-rowSums(DESeq2::counts(dds)>=10L)>=
    max(2L,ceiling(ncol(dds)*.2))
  dds<-dds[keep,,drop=FALSE]
  if(nrow(dds)<300L||anyDuplicated(rownames(dds))>0L)
    stop("Unexpected low-count filtering of the source null study.")
  truth<-source_study$dispersion[
    match(rownames(dds),colnames(counts))]
  if(any(!is.finite(truth))||any(truth<=0)||
     length(truth)!=nrow(dds))
    stop("Invalid or incorrectly aligned original gene-wise NB true alphas.")
  # Split the documented native DESeq2 fitting pipeline. Estimate once
  # per cohort, not independently in every arm.
  estimated<-DESeq2::estimateDispersions(
    DESeq2::estimateSizeFactors(dds),quiet=TRUE)
  meta<-S4Vectors::mcols(estimated)
  components<-list(
    map_final=as.numeric(DESeq2::dispersions(estimated)),
    gene_wise_CR=as.numeric(meta$dispGeneEst),
    fitted_trend=as.numeric(meta$dispFit),
    true_nb_alpha=as.numeric(truth))
  if(!all(vapply(components,length,integer(1))==nrow(dds))||
     any(vapply(components,function(z)any(!is.finite(z)|z<=0),
                logical(1))))
    stop("Native MAP/gene-wise/trend alphas missing or nonfinite.")
  if(anyDuplicated(rownames(estimated))>0L||
     !identical(rownames(estimated),rownames(dds)))
    stop("Native DESeq2 gene identity drift during dispersion estimation.")
  size_factors<-as.numeric(DESeq2::sizeFactors(estimated))
  if(length(size_factors)!=ncol(dds)||
     any(!is.finite(size_factors)|size_factors<=0))
    stop("Normalized native sample depths changed unexpectedly.")
  for(arm in arms){
    model<-estimated
    # Change ONLY per-gene dispersions on the already fitted same
    # size factors and design. Wald API is exactly the same for all arms.
    DESeq2::dispersions(model)<-components[[arm]]
    if(max(abs(log(DESeq2::dispersions(model)/components[[arm]])))>1e-10)
      stop("One DESeq2 component arm did not preserve the supplied dispersion.")
    model<-DESeq2::nbinomWaldTest(model,betaPrior=FALSE,quiet=TRUE)
    if(max(abs(log(DESeq2::sizeFactors(model)/size_factors)))>1e-10)
      stop("Wald run altered the shared size factors.")
    # No Cook cutoff in any arm to isolate dispersion effect.
    # DESeq2's independentFiltering remains enabled in every arm.
    result<-DESeq2::results(model,
      contrast=c("condition","treated","control"),alpha=.05,
      cooksCutoff=FALSE,independentFiltering=TRUE)
    p<-as.numeric(result$pvalue)
    q<-as.numeric(result$padj)
    w<-as.numeric(result$stat)
    if(length(p)!=nrow(dds)||length(q)!=nrow(dds)||
       length(w)!=nrow(dds)||
       any(!is.na(p)&(!is.finite(p)|p<0|p>1))||
       any(!is.na(q)&(!is.finite(q)|q<0|q>1))||
       any(is.finite(q)&!is.finite(p))||
       sum(is.finite(q))<300L||sum(is.finite(w))<300L)
      stop("Native Wald statistic or p/q family invalid.")
    eligible<-is.finite(q)
    replay_error<-max(abs(independent_bh_adjust(p[eligible])-q[eligible]))
    if(!is.finite(replay_error)||replay_error>1e-8)
      stop("Numerical BH mismatch between independent oracle and native DESeq2.")
    position<-position+1L
    records[[position]]<-data.frame(
      scenario=scenario,replicate=as.integer(id),
      seed=as.integer(source_study$seed),arm=arm,
      source_features=400L,modeled_features=nrow(dds),
      finite_p=sum(is.finite(p)),finite_q=sum(eligible),
      false_calls_bh05=sum(eligible&q<=.05),
      any_false_call=as.integer(any(eligible&q<=.05)),
      raw_p_under_005=sum(is.finite(p)&p<=.05),
      raw_p_under_001=sum(is.finite(p)&p<=.001),
      mean_wald_stat_squared=mean(w[is.finite(w)]^2),
      dispersion_log_rmse=sqrt(mean(log(components[[arm]]/truth)^2)),
      median_log_dispersion_error=stats::median(log(components[[arm]]/truth)),
      size_factor_log_rmse=sqrt(mean(log(
        size_factors/(source_study$depth/
          exp(mean(log(source_study$depth)))))^2)),
      native_BH_replay_ok=TRUE,
      native_bh_error_max=replay_error,
      statistical_model="DESeq2_native_nbinomWaldTest_no_Cook_filter",
      science_calibrated=FALSE,
      stringsAsFactors=FALSE)
  }
  cat("DISPERSION COMPONENT",scenario,idx,"/",reps,
      "four original-matrix native Wald tests completed\n",flush=TRUE)
}
dat<-do.call(rbind,records)
summary<-rnaseq_deseq2_dispersion_component_summary(dat,
          per_scenario=reps,replicate_start=first_rep)
if(nrow(dat)!=384L||nrow(summary)!=12L)
  stop("Not every planned original global-null gene-wise/trend/MAP/truth experiment ran.")
utils::write.csv(dat,file.path(out,
  "deseq2_dispersion_component_null_per_study.csv"),row.names=FALSE)
utils::write.csv(summary,file.path(out,
  "deseq2_dispersion_component_null_summary.csv"),row.names=FALSE)
writeLines(c(
  "schema: MultiOmics.DESeq2.dispersion_components.v1",
  paste("source_SHA:",Sys.getenv("GITHUB_SHA","local_unverified")),
  "real_native_R_DESeq2_Wald_fits: 384",
  "fixed_replicates: 4001 to 4032; 32 per each of three original null scenarios",
  "sample_and_feature_source: same 400-source-gene 40-subject NB generator for all four arms",
  "native_gene_wise_CR: DESeq2 mcols dispGeneEst, Cox-Reid adjusted",
  "native_fitted_trend: DESeq2 mcols dispFit",
  "native_final_MAP: DESeq2 dispersions after estimateDispersions",
  "simulation_true_alpha: original generator gene-specific exact negative-binomial dispersion",
  "size_factors: native DESeq2 estimates shared exactly in all four arms",
  "Wald: same exported nbinomWaldTest betaPrior FALSE for all arms",
  "Cook_filter: FALSE in all four arms to isolate dispersion component",
  "independent_filter: TRUE in all four arms; BH independently checked",
  "original_production_algorithm_modified: false",
  "substantive_FDR_validity_or_causal_mechanism_certified: false",
  capture.output(sessionInfo())),
  file.path(out,"deseq2_dispersion_component_method_manifest.txt"))
print(summary,row.names=FALSE)
cat("PASS real native DESeq2 four-component ablation; nominal FDR unproven.\n")
