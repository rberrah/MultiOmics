#!/usr/bin/env Rscript
# Optional but real reference-package execution: DESeq2, edgeR QL and limma-voom.
# It must run with the full Bioconductor toolchain (see scheduled CI workflow).
# This controlled simulation is a smoke test for adapter direction & selective
# detection, not global calibration or proof of validity on patient cohorts.
source(file.path("multiomics-engine", "advanced_methods.R"))
for (pkg in c("DESeq2","edgeR","limma","statmod")) {
  if (!requireNamespace(pkg, quietly=TRUE)) stop("Missing required package: ",pkg)
}
set.seed(20261008)
genes <- paste0("G",sprintf("%03d",seq_len(300)))
ids <- paste0("S",sprintf("%02d",seq_len(40)))
group <- factor(rep(c("control","treated"),each=20),levels=c("control","treated"))
batch <- factor(rep(c("B1","B2"),times=20))
meta <- data.frame(subject_id=paste0("BIO",seq_along(ids)),condition=group,batch=batch,row.names=ids)
base <- exp(stats::runif(length(genes),log(60),log(300)))
dispersion <- stats::runif(length(genes),0.08,0.18)
is_signal <- seq_along(genes)<=25
library_factor <- stats::runif(length(ids),0.90,1.10)
counts <- vapply(seq_along(ids),function(i) {
  fold <- if(group[i]=="treated") ifelse(is_signal,4.0,1) else rep(1,length(genes))
  stats::rnbinom(length(genes),mu=base*fold*library_factor[i],size=1/dispersion)
},numeric(length(genes)))
rownames(counts) <- genes
colnames(counts) <- ids
x <- t(counts)
out_dir <- tempfile("bioconductor_reference_")
dir.create(out_dir)
d <- run_deseq2_counts(x,meta,file.path(out_dir,"deseq2"),
    design_formula=~ batch+condition,contrast=c("condition","treated","control"),
    confirmatory_contract=TRUE)
v <- run_voom_counts(x,meta,file.path(out_dir,"voom"),
    design_formula=~ batch+condition,coefficient="conditiontreated")
e <- run_edger_ql_counts(x,meta,file.path(out_dir,"edger"),
    design_formula=~ batch+condition,coefficient="conditiontreated")
stopifnot(nrow(d$results)>200L,nrow(v$results)>200L,nrow(e$results)>200L)
stopifnot(identical(d$confirmatory_input_contract$status,
                    "eligible_for_reference_fitting_only"),
          d$confirmatory_input_contract$independent_subjects==length(ids),
          d$confirmatory_input_contract$residual_df==37L,
          identical(d$confirmatory_input_contract$certified,FALSE))
cat("REAL DESeq2 UNDER R1 STRICT CONFIRMATORY INPUT CONTRACT: PASS (no biological certificate)\n")
reference <- compare_rnaseq_methods(d$results,v$results)
stopifnot(identical(reference$status,"ok"))
method_specs <- list(
  DESeq2=list(tab=d$results,effect="log2FoldChange",q="padj"),
  limma_voom=list(tab=v$results,effect="logFC",q="adj.P.Val"),
  edgeR_QL=list(tab=e$results,effect="logFC",q="FDR")
)
for (name in names(method_specs)) {
  spec <- method_specs[[name]]
  tab <- spec$tab
  signal <- tab$feature %in% genes[is_signal]
  null <- !signal
  effect <- as.numeric(tab[[spec$effect]])
  fdr <- as.numeric(tab[[spec$q]])
  detected <- is.finite(fdr) & fdr<=0.05
  signal_detection <- sum(detected & signal)/sum(signal)
  false_discoveries <- sum(detected & null)
  signal_direction <- mean(effect[signal]>0,na.rm=TRUE)
  stopifnot(signal_detection>=0.60, signal_direction>=0.85,
    false_discoveries<=20L)
  cat(sprintf("%s | signal recall %.3f | false discoveries %d | direction %.3f\n",
    name,signal_detection,false_discoveries,signal_direction))
}
stopifnot(is.finite(reference$spearman_statistics),
    reference$spearman_statistics >= 0.70)
cat(sprintf("DESeq2↔voom rank-statistic agreement %.3f | PASS\n",
  reference$spearman_statistics))

# A non-integer / transformed matrix must not be silently coerced for DESeq2.
bad <- x
bad[1,1] <- bad[1,1]+0.25
err <- try(run_deseq2_counts(bad,meta,file.path(out_dir,"reject"),
  design_formula=~batch+condition),silent=TRUE)
stopifnot(inherits(err,"try-error"))
cat("DESeq2 non-integer input rejection: PASS\n")

# Replicated negative-control screen: all genes have zero true treatment effect.
# This guards against gross false-positive inflation, not formal FDR validation.
null_replicates <- 4L
null_summary <- setNames(vector("list",3L),c("DESeq2","limma_voom","edgeR_QL"))
for (name in names(null_summary))
  null_summary[[name]] <- c(discoveries=0,raw_below_005=0,tests=0)
for (rep in seq_len(null_replicates)) {
  set.seed(20261009 + rep)
  null_x <- t(vapply(seq_along(ids),function(i)
    stats::rnbinom(length(genes),mu=base*library_factor[i],size=1/dispersion),
    numeric(length(genes))))
  rownames(null_x) <- ids
  colnames(null_x) <- genes
  null_dir <- file.path(out_dir,paste0("null_",rep))
  dn <- run_deseq2_counts(null_x,meta,file.path(null_dir,"deseq2"),
    design_formula=~batch+condition,contrast=c("condition","treated","control"))
  vn <- run_voom_counts(null_x,meta,file.path(null_dir,"voom"),
    design_formula=~batch+condition,coefficient="conditiontreated")
  en <- run_edger_ql_counts(null_x,meta,file.path(null_dir,"edger"),
    design_formula=~batch+condition,coefficient="conditiontreated")
  methods <- list(DESeq2=list(tab=dn$results,q="padj",p="pvalue"),
    limma_voom=list(tab=vn$results,q="adj.P.Val",p="P.Value"),
    edgeR_QL=list(tab=en$results,q="FDR",p="PValue"))
  for (name in names(methods)) {
    info <- methods[[name]]
    p <- as.numeric(info$tab[[info$p]])
    q <- as.numeric(info$tab[[info$q]])
    stopifnot(length(p)==length(genes),length(q)==length(genes),
      all(is.na(p)|(p>=0 & p<=1)),all(is.na(q)|(q>=0 & q<=1)))
    null_summary[[name]]["discoveries"] <- null_summary[[name]]["discoveries"] +
      sum(q<=0.05,na.rm=TRUE)
    null_summary[[name]]["raw_below_005"] <- null_summary[[name]]["raw_below_005"] +
      sum(p<0.05,na.rm=TRUE)
    null_summary[[name]]["tests"] <- null_summary[[name]]["tests"] + sum(is.finite(p))
  }
}
for (name in names(null_summary)) {
  s <- null_summary[[name]]
  rate <- s["raw_below_005"]/s["tests"]
  stopifnot(s["tests"]>=0.9*length(genes)*null_replicates,
    rate<=0.10,s["discoveries"]<=12L)
  cat(sprintf("%s NULL | raw p<0.05 %.3f | q<=0.05 %d / %d | PASS (gross-regression bounds)\n",
    name,rate,as.integer(s["discoveries"]),as.integer(s["tests"])))
}

unlink(out_dir,recursive=TRUE)
