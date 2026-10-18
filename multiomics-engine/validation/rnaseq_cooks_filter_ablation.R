#!/usr/bin/env Rscript
# DESeq2 results() nuisance-screening ablation. One original native DESeq2
# WALD model fit is reused for FOUR Cook/independent-filter result tables.
# Outlier count replacement is explicitly disabled for ALL four result arms.
# 192 whole null studies (96 fresh original NB, 96 independently generated
# Gamma-Poisson studies re-used transparently as paired diagnostic cohorts).
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/rnaseq_independent_gamma_poisson_null.R")
source("multiomics-engine/rnaseq_cooks_filter_null_metrics.R")
source("multiomics-engine/reference_inference_audit.R")
if(!requireNamespace("DESeq2",quietly=TRUE))
  stop("Native pinned Bioconductor DESeq2 is required.")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L||!nzchar(args[[1L]]))
  stop("Usage: Rscript rnaseq_cooks_filter_ablation.R <empty evidence dir>")
out<-args[[1L]]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Never overwrite prior native scientific calibration results.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
original<-c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
gp<-c("gp_balanced_heterogeneous","gp_depth_heterogeneous",
      "gp_unbalanced_batch")
arms<-c("cooks_on_if_on","cooks_on_if_off",
        "cooks_off_if_on","cooks_off_if_off")
records<-list()
for(sc in c(original,gp)) {
  from_original<-sc%in%original
  reps<-if(from_original)6001:6032 else 5001:5032
  for(id in reps) {
    study<-if(from_original)rnaseq_oracle_null_study(sc,id)
      else independent_gp_null_study(sc,id)
    if(!identical(rownames(study$counts),rownames(study$meta))||
       nrow(study$counts)!=40L||ncol(study$counts)!=400L||
       any(table(study$meta$condition,study$meta$batch)==0L))
      stop("Generator source counts, subject labels or design are invalid.")
    dds<-DESeq2::DESeqDataSetFromMatrix(
      countData=t(study$counts),colData=study$meta,
      design=~batch+condition)
    # Exact original model prefilter; never change gene family by arm.
    keep<-rowSums(DESeq2::counts(dds)>=10L)>=
      max(2L,ceiling(ncol(dds)*.2))
    dds<-dds[keep,,drop=FALSE]
    if(nrow(dds)<300L||anyDuplicated(rownames(dds)))
      stop("Source gene family was unexpectedly discarded.")
    # Disabling count replacement lets results() Cook's distance SCREEN
    # be compared alone; this is intentionally not the full default
    # production procedure when at least seven replicates are available.
    fitted<-DESeq2::DESeq(
      dds,test="Wald",minReplicatesForReplace=Inf,quiet=TRUE)
    if(!"cooks"%in%names(SummarizedExperiment::assays(fitted)))
      stop("No native Cook's distances were fitted.")
    model_hash<-paste(study$scenario,study$seed,nrow(fitted),sep=":")
    outputs<-list()
    for(arm in arms) {
      keep_cooks<-arm%in%c("cooks_on_if_on","cooks_on_if_off")
      keep_if<-arm%in%c("cooks_on_if_on","cooks_off_if_on")
      if(keep_cooks) {
        res<-DESeq2::results(fitted,contrast=c("condition","treated","control"),
          independentFiltering=keep_if,alpha=.05)
      } else {
        res<-DESeq2::results(fitted,contrast=c("condition","treated","control"),
          cooksCutoff=FALSE,independentFiltering=keep_if,alpha=.05)
      }
      p<-as.numeric(res$pvalue)
      q<-as.numeric(res$padj)
      if(length(p)!=nrow(dds)||length(q)!=nrow(dds)||
         any(!is.na(p)&(!is.finite(p)|p<0|p>1))||
         any(!is.na(q)&(!is.finite(q)|q<0|q>1))||
         any(is.finite(q)&!is.finite(p))||
         sum(is.finite(q))<100L)
        stop("Native DESeq2 failed p/q-family validity after Cook/filter screening.")
      used<-is.finite(q)
      if(!length(which(used)))stop("No eligible native BH p-values.")
      deviation<-max(abs(independent_bh_adjust(p[used])-q[used]))
      if(!is.finite(deviation)||deviation>1e-8)
        stop("Independent BH replay failed with Cook/filter screening.")
      outputs[[arm]]<-list(p=p,q=q,eligible=used,
        keep_cooks=keep_cooks,keep_if=keep_if,
        bh_error=deviation)
    }
    on_if<-outputs[["cooks_on_if_on"]]
    on_noif<-outputs[["cooks_on_if_off"]]
    off_if<-outputs[["cooks_off_if_on"]]
    off_noif<-outputs[["cooks_off_if_off"]]
    equal_p<-function(a,b) {
      same_na<-identical(is.na(a),is.na(b))
      both<-which(is.finite(a)&is.finite(b))
      same_na&&length(both)>0L&&max(abs(a[both]-b[both]))<1e-12
    }
    if(!equal_p(on_if$p,on_noif$p)||
       !equal_p(off_if$p,off_noif$p)||
       any(is.finite(on_if$p)&!is.finite(off_if$p))||
       any(is.finite(on_noif$p)&!is.finite(off_noif$p))||
       any(on_if$eligible&!on_noif$eligible)||
       any(off_if$eligible&!off_noif$eligible))
      stop("Cook's cutoff or independent filtering changed the wrong p/q family.")
    cook_exclusions<-sum(is.finite(off_noif$p)&!is.finite(on_noif$p))
    for(arm in arms) {
      x<-outputs[[arm]]
      source_name<-if(from_original)"original_smooth_NB" else
                      "independent_heterogeneous_GP"
      records[[length(records)+1L]]<-data.frame(
        scenario=sc,source_generator=source_name,
        replicate=as.integer(id),seed=as.integer(study$seed),arm=arm,
        source_genes=400L,modeled_genes=nrow(dds),
        total_native_fit_count=1L,
        finite_p=sum(is.finite(x$p)),finite_q=sum(x$eligible),
        discovery_count=sum(x$eligible & x$q<=.05),
        any_discovery=as.integer(any(x$eligible & x$q<=.05)),
        p_removed_by_cooks=if(x$keep_cooks)as.integer(cook_exclusions) else 0L,
        q_removed_by_filter=as.integer(sum(is.finite(x$p)&!is.finite(x$q))),
        raw_p_005=sum(is.finite(x$p)&x$p<=.05),
        native_bh_replay_ok=TRUE,
        native_bh_max_deviation=x$bh_error,
        same_model_pvalue_check=TRUE,
        common_native_model_identifier=model_hash,
        original_count_replacement_disabled=TRUE,
        FDR_certified=FALSE,
        stringsAsFactors=FALSE)
    }
    cat("COOK/IF",source_name,sc,id,
        "one true native fitted model / four BH tables PASS\n",flush=TRUE)
  }
}
rows<-do.call(rbind,records)
summary<-rnaseq_cook_if_summary(rows)
if(nrow(rows)!=768L||nrow(summary)!=24L)
  stop("Incomplete 192-original-fits × 4-results study.")
utils::write.csv(rows,file.path(out,
  "cooks_independent_filter_per_study.csv"),row.names=FALSE)
utils::write.csv(summary,file.path(out,
  "cooks_independent_filter_24_cell_summary.csv"),row.names=FALSE)
writeLines(c(
  "schema: MultiOmics.DESeq2.Cook_and_independent_filter_global_null.v1",
  paste("source_SHA:",Sys.getenv("GITHUB_SHA","local_unverified")),
  "native_DESeq2_refitted_models: 192; each has four results() calls",
  "original_NB_generator_replicates: 6001–6032 (fresh)",
  "independent_GammaPoisson_replicates: 5001–5032 (re-use of prior independent generator studies, newly fitted and paired across screens)",
  "condition_true_effect: ZERO for every tested gene",
  "Cook_count_replacement: disabled in every model; no post-hoc model refit",
  "Cook_screening: results() native default versus cooksCutoff=FALSE",
  "independent_filtering: native true versus false",
  "all_arms_use_identical_size_factor_dispersion_Wald_and_source_design: true",
  "original_gene_low_count_prefilter: 10 counts in at least 20% samples",
  "native_results_alpha: .05",
  "BenjaminiHochberg_integrity: independently reconstructed for all 768 q families",
  "global_null_FDR_equals_study_FWER: true",
  "individual_exact_binomial_uncertainty: 95pct; never a simultaneous validity guarantee",
  "genuine_clinical_FDR_calibration_certified: false",
  "production_method_or_results_changed: false",
  capture.output(sessionInfo())),
  file.path(out,"cooks_independent_filter_method_manifest.txt"))
print(summary,row.names=FALSE)
cat("PASS: 192 complete real R fits, 768 native screening outcomes, no inferential certification\n")
