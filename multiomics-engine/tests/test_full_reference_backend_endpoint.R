#!/usr/bin/env Rscript
# Full application-level integration, not just a serialization unit test.
# Real RNA-seq counts => run_backend_analysis => Plumber-style JSON =>
# all adjusted feature rows preserved and associated to the right method.
source(file.path("multiomics-engine","server.R"))
for(p in c("DESeq2","jsonlite"))if(!requireNamespace(p,quietly=TRUE))
  stop("Real reference test requires ",p)
set.seed(30361013)
n<-28L;ng<-172L
sample_ids<-sprintf("ASSAY%03d",seq_len(n))
condition<-rep(c("control","treated"),each=n/2)
batch<-rep(rep(c("B1","B2"),each=7),times=2)
counts<-vapply(seq_len(ng),function(k){
  effect<-if(k<=24)ifelse(condition=="treated",1.8,1) else 1
  stats::rnbinom(n,mu=80*effect,size=1/.15)
},numeric(n))
rownames(counts)<-sample_ids
colnames(counts)<-sprintf("GENE_%04d",seq_len(ng))
metadata<-data.frame(
  subject_id=sample_ids,sample_id=sample_ids,assay_id=sample_ids,
  omic="transcriptomics",condition=condition,
  timepoint="T0",batch=batch,technical_replicate="1",
  sample_type="biological",stringsAsFactors=FALSE
)
measurements<-data.frame(feature_id=colnames(counts),
                          t(counts),check.names=FALSE)
csv_text<-function(x)paste0(paste(capture.output(
  utils::write.csv(x,row.names=FALSE,quote=TRUE)),collapse="\n"),"\n")
payload<-list(
 protocol=list(objective="groups",organism="human",
    designType="independent",analysisIntent="exploratory",
    longitudinal=FALSE,outcomeType="none",rnaCountMethod="deseq2",
    covariateColumns=list()),
 dataTypes=list(transcriptomics="raw_counts"),
 identifierTypes=list(transcriptomics="ensembl_gene"),
 columnMapping=list(subject_id="subject_id",sample_id="sample_id",
    assay_id="assay_id",omic="omic",condition="condition",
    timepoint="timepoint",batch="batch",
    technical_replicate="technical_replicate",sample_type="sample_type"),
 metadataCsv=csv_text(metadata),
 matrices=list(transcriptomics=csv_text(measurements))
)
response<-run_backend_analysis(payload)
stopifnot(identical(response$status,"ok"))
m<-response$methods$transcriptomics_differential_deseq2
if(is.null(m)||m$status!="ok")stop("Real DESeq2 run did not succeed")
stopifnot(isTRUE(m$referenceTableComplete),
          identical(m$referenceTableStatus,"complete"),
          m$referenceTableRows > 100L,
          m$referenceTableRows==m$referenceTableExportedRows,
          length(m$top)==100L,
          identical(m$referenceTableOriginalFeatures,ng),
          m$referenceTableExcludedBeforeFit>=0L)
# Explicitly exercise the existing serializer contract used by Plumber.
json<-jsonlite::toJSON(response,auto_unbox=TRUE,dataframe="rows",na="null")
again<-jsonlite::fromJSON(json,simplifyVector=FALSE)
mm<-again$methods$transcriptomics_differential_deseq2
stopifnot(isTRUE(mm$referenceTableComplete),
          mm$referenceTableRows==m$referenceTableRows)
rows<-utils::read.csv(text=mm$referenceTableCsv,check.names=FALSE,
                       stringsAsFactors=FALSE,na.strings="NA")
stopifnot(nrow(rows)==mm$referenceTableRows,
          nrow(rows)>100L,
          length(unique(rows$feature))==nrow(rows),
          all(rows$feature%in%colnames(counts)),
          all(c("log2FoldChange","pvalue","padj")%in%names(rows)),
          any(rows$feature=="GENE_0172"))
cat("PMX REAL BACKEND INTEGRATION: PASS, true DESeq2 fitted rows",
    nrow(rows),"top preview",length(m$top),"original",ng,"\n")
