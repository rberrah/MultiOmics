#!/usr/bin/env Rscript
source("multiomics-engine/rnaseq_deseq2_null_oracle.R")
source("multiomics-engine/rnaseq_deseq2_null_tail_metrics.R")
source("multiomics-engine/reference_inference_audit.R")
if(!requireNamespace("DESeq2",quietly=TRUE))stop("Pinned native DESeq2 required.")
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript rnaseq_deseq2_null_tail_diagnostic.R EMPTY_OUTPUT_DIR")
out<-args[[1L]]
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))
  stop("Refuse to overwrite prior scientific evidence.")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
scenarios<-c("balanced_base","balanced_depth_dispersion","unequal_partial_batch")
methods<-c("native_default","oracle_true_dispersion")
rows<-list()
for(scenario in scenarios)for(replicate in 6001:6032){
  study<-rnaseq_oracle_null_study(scenario,replicate)
  dds<-DESeq2::DESeqDataSetFromMatrix(countData=t(study$counts),
    colData=study$meta,design=~batch+condition)
  keep<-rowSums(DESeq2::counts(dds)>=10L)>=max(2L,ceiling(ncol(dds)*.2))
  dds<-dds[keep,]
  truth<-study$dispersion[match(rownames(dds),colnames(study$counts))]
  if(nrow(dds)<300L||any(!is.finite(truth))||any(truth<=0))
    stop("Generator/fitted-feature contract changed.")
  for(method in methods){
    fit<-dds
    if(method=="native_default"){
      fit<-DESeq2::DESeq(fit,quiet=TRUE)
    } else {
      fit<-DESeq2::estimateSizeFactors(fit)
      DESeq2::dispersions(fit)<-truth
      fit<-DESeq2::nbinomWaldTest(fit,betaPrior=FALSE,quiet=TRUE)
    }
    res<-DESeq2::results(fit,contrast=c("condition","treated","control"),
      alpha=.05,cooksCutoff=FALSE,independentFiltering=FALSE)
    p<-as.numeric(res$pvalue);q<-as.numeric(res$padj);z<-as.numeric(res$stat)
    active<-is.finite(p)
    replay<-rep(NA_real_,length(p))
    replay[active]<-independent_bh_adjust(p[active])
    if(any(is.finite(q)!=active)||max(abs(q[active]-replay[active]))>1e-8)
      stop("Native BH family differs from independent replay.")
    met<-rnaseq_null_tail_metrics(p,z)
    rows[[length(rows)+1L]]<-cbind(data.frame(
      scenario=scenario,method=method,replicate=replicate,seed=study$seed,
      n_genes_modelled=nrow(dds),bh_discoveries=sum(q<=.05,na.rm=TRUE),
      stringsAsFactors=FALSE),met)
  }
  cat("NULL TAIL",scenario,replicate,"complete\n",flush=TRUE)
}
all<-do.call(rbind,rows)
summary<-rnaseq_null_tail_summary(all)
utils::write.csv(all,file.path(out,"per_study_null_tail_metrics.csv"),row.names=FALSE)
utils::write.csv(summary$summary,file.path(out,"scenario_method_tail_summary.csv"),row.names=FALSE)
utils::write.csv(summary$paired,file.path(out,"paired_default_minus_true_dispersion.csv"),row.names=FALSE)
writeLines(c(
 "schema: multiomics.deseq2-null-tail-diagnostic.v1",
 paste("github_sha:",Sys.getenv("GITHUB_SHA","local_unverified")),
 "replicates: 6001-6032; 32 independent studies per scenario",
 "methods: native_default; oracle_true_dispersion",
 "cooksCutoff: false for raw-tail isolation",
 "independentFiltering: false for raw-tail isolation",
 "genes_are_not_independent_calibration_replicates: true",
 "production_inference_changed: false",
 "universal_FDR_certified: false",
 capture.output(sessionInfo())
),file.path(out,"manifest.txt"))
cat("PASS: 192 native DESeq2 whole-study fits; null p/Wald tails retained without pooled-gene calibration claims.\n")
