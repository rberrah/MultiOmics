#!/usr/bin/env Rscript
# REAL R fits, no mocked reference-engine outputs.
source(file.path("multiomics-engine","advanced_methods.R"))
source(file.path("multiomics-engine","reference_full_results.R"))
for(p in c("DESeq2","edgeR","limma","jsonlite","statmod"))
  if(!requireNamespace(p,quietly=TRUE))stop("Missing genuine R package ",p)
set.seed(20261013)
n<-28L;genes<-168L
samples<-sprintf("S%02d",seq_len(n))
groups<-factor(rep(c("control","treated"),each=n/2),levels=c("control","treated"))
batch<-factor(rep(rep(c("batch1","batch2"),each=7),times=2))
meta<-data.frame(condition=groups,batch=batch,row.names=samples)
input<-vapply(seq_len(genes),function(j){
  mu<-60*exp(rnorm(1,0,.25))*
    ifelse(groups=="treated" & j<=20,2.1,1)
  stats::rnbinom(n,mu=mu,size=1/0.14)
},numeric(n))
dimnames(input)<-list(samples,sprintf("ENSG_SYN%04d",seq_len(genes)))
out<-tempfile("pmx_reference_real_");dir.create(out)

fits<-list(
 DESeq2=run_deseq2_counts(input,meta,file.path(out,"deseq"),
                         ~batch+condition,c("condition","treated","control")),
 voom=run_voom_counts(input,meta,file.path(out,"voom"),~batch+condition),
 limma=run_limma_matrix(log2(input+1),meta,file.path(out,"limma"),~batch+condition)
)
for(method in names(fits)){
  full<-fits[[method]]$results
  stopifnot(nrow(full)>100L)
  exported<-reference_complete_table(full,original_features=ncol(input))
  stopifnot(exported$referenceTableComplete,
            exported$referenceTableRows==nrow(full),
            exported$referenceTableExportedRows==nrow(full),
            identical(exported$referenceTableOriginalFeatures,genes))
  # Simulates actual Plumber json serializer; CSV MUST still parse losslessly.
  wire<-jsonlite::toJSON(exported,auto_unbox=TRUE,na="null")
  received<-jsonlite::fromJSON(wire,simplifyVector=FALSE)
  stopifnot(isTRUE(received$referenceTableComplete))
  table<-utils::read.csv(text=received$referenceTableCsv,check.names=FALSE,
                         stringsAsFactors=FALSE,na.strings="NA")
  stopifnot(nrow(table)==nrow(full),nrow(table)>100L,
            identical(table$feature,full$feature),
            identical(names(table),names(full)))
  for(name in intersect(c("padj","pvalue","adj.P.Val","P.Value",
                          "log2FoldChange","logFC","t"),names(full))){
    stopifnot(identical(which(is.na(table[[name]])),which(is.na(full[[name]]))),
              isTRUE(all.equal(table[[name]],full[[name]],tolerance=1e-9)))
  }
  cat("REAL COMPLETE REFERENCE FIT:",method,nrow(table),"of",genes,
      "input genes preserved in HTTP JSON CSV\n")
}
cat("REAL R DESeq2+VOOM+LIMMA FULL REFERENCE EXPORT: PASS\n")
