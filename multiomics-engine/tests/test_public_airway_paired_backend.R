#!/usr/bin/env Rscript
# Public paired biological study through the REAL run_backend_analysis API
# (not merely a wrapper around DESeq2).
# GEO GSE52778 / Bioconductor airway, four paired cell lines, 8 RNA-seq assays.
suppressPackageStartupMessages({library(airway);library(DESeq2)})
source(file.path("multiomics-engine","server.R"))
stopifnot(requireNamespace("jsonlite",quietly=TRUE))
data("airway",package="airway")
all_counts<-SummarizedExperiment::assay(airway)
stopifnot(ncol(all_counts)==8L,nrow(all_counts)>60000)
mask<-rowSums(all_counts>=10)>=max(2L,ceiling(ncol(all_counts)*0.2))
counts<-all_counts[mask,,drop=FALSE]
ids<-colnames(counts)
subjects<-as.character(airway$cell)
cond<-as.character(airway$dex)
stopifnot(length(unique(subjects))==4L,identical(sort(unique(cond)),c("trt","untrt")),
          all(table(subjects)==2L),nrow(counts)>10000)
sheet<-data.frame(
 subject_id=subjects,sample_id=ids,assay_id=ids,
 omic="transcriptomics",condition=cond,timepoint="T0",batch="RUN1",
 technical_replicate="1",sample_type="biological",stringsAsFactors=FALSE)
expression<-data.frame(feature_id=rownames(counts),counts,check.names=FALSE)
csv<-function(x)paste0(paste(capture.output(utils::write.csv(
  x,row.names=FALSE,quote=TRUE)),collapse="\n"),"\n")
payload<-list(
 protocol=list(objective="groups",designType="paired",
   analysisIntent="confirmatory",longitudinal=FALSE,rnaCountMethod="deseq2",
   rnaPrimaryContrast=as.list(c("condition","trt","untrt")),
   covariateColumns=list()),
 dataTypes=list(transcriptomics="raw_counts"),
 identifierTypes=list(transcriptomics="ensembl_gene"),
 columnMapping=as.list(stats::setNames(names(sheet),names(sheet))),
 metadataCsv=csv(sheet),
 matrices=list(transcriptomics=csv(expression))
)
output<-run_backend_analysis(payload)
stopifnot(identical(output$status,"ok"))
r<-output$methods$transcriptomics_differential_deseq2
stopifnot(r$status=="ok",isTRUE(r$referenceTableComplete),
          r$pairedDesignContract$independent_biological_subjects==4L,
          r$pairedDesignContract$biological_pairs==4L,
          r$pairedDesignContract$model_samples==8L,
          r$pairedDesignContract$contrast_was_prespecified,
          r$referenceTableRows==sum(mask),
          r$referenceTableExportedRows==sum(mask))
observed<-utils::read.csv(text=r$referenceTableCsv,check.names=FALSE,
                          stringsAsFactors=FALSE,na.strings="NA")
stopifnot(nrow(observed)==sum(mask),any(observed$feature=="ENSG00000103196"))
# True independently assembled native DESeq2 model, blocked by the cell line.
m<-data.frame(subject_id=factor(subjects),condition=factor(cond,
              levels=sort(unique(cond))),row.names=ids)
dds<-DESeq2::DESeqDataSetFromMatrix(countData=counts,colData=m,
                                   design=~subject_id+condition)
dds<-DESeq2::DESeq(dds,quiet=TRUE)
native<-DESeq2::results(dds,contrast=c("condition","trt","untrt"),alpha=.05)
ref<-as.data.frame(native)[match(observed$feature,rownames(native)),,drop=FALSE]
stopifnot(identical(observed$feature,rownames(ref)))
for(k in c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")){
  same_na<-identical(which(is.na(observed[[k]])),which(is.na(ref[[k]])))
  numeric_equal<-isTRUE(all.equal(observed[[k]],ref[[k]],tolerance=1e-8))
  if(!same_na||!numeric_equal){
    kept<-which(is.finite(observed[[k]])&is.finite(ref[[k]]))
    maximum<-if(length(kept))max(abs(observed[[k]][kept]-ref[[k]][kept])) else NA_real_
    stop("Airway paired wrapper-native comparison failed on ",k,
      "; NA positions identical=",same_na,"; maximal absolute discrepancy=",maximum,
      ". Check contrast coding, design rank and input precision.")
  }
}
marker<-observed[observed$feature=="ENSG00000103196",,drop=FALSE]
stopifnot(nrow(marker)==1L,is.finite(marker$padj),
          marker$padj<.05,marker$log2FoldChange>0)
# Test the actual REST JSON representation and ensure full reference CSV
# remains byte-for-byte preserved.
json<-jsonlite::toJSON(output,auto_unbox=TRUE,dataframe="rows",na="null")
parsed<-jsonlite::fromJSON(json,simplifyVector=FALSE)
stopifnot(isTRUE(parsed$methods$transcriptomics_differential_deseq2$referenceTableComplete),
          identical(parsed$methods$transcriptomics_differential_deseq2$referenceTableCsv,
                    r$referenceTableCsv))
report<-list(accession="GSE52778",data="airway",n_samples=8L,
  independent_biological_pairs=4L,design="~ subject_id + condition",
  treatment_vs_reference=c("condition","trt","untrt"),
  prefiltered_genes=sum(mask),excluded_before_fitting=sum(!mask),
  complete_reference_rows=nrow(observed),
  CRISPLD2_log2_fold_change=unname(marker$log2FoldChange),
  CRISPLD2_q=unname(marker$padj),
  backend="run_backend_analysis",comparison="direct native DESeq2",
  software_integration_test=TRUE,independent_cohort_validation=FALSE,
  universal_FDR_validation=FALSE,clinical_validation=FALSE)
reportDir<-Sys.getenv("PMX_PAIRED_AIRWAY_REPORT_DIR","")
if(nzchar(reportDir)){
  dir.create(reportDir,recursive=TRUE,showWarnings=FALSE)
  jsonlite::write_json(report,file.path(reportDir,"paired_airway_backend.json"),
                       pretty=TRUE,auto_unbox=TRUE,null="null")
  # The browser R-primary numerical + figure integration test uses precisely
  # this PUBLIC backend request/response, rather than synthetic rows.
  writeLines(jsonlite::toJSON(payload,auto_unbox=TRUE,dataframe="rows",na="null",
            null="null",digits=NA),
             file.path(reportDir,"paired_airway_original_request.json"),useBytes=TRUE)
  writeLines(jsonlite::toJSON(output,auto_unbox=TRUE,dataframe="rows",na="null",
            null="null",digits=NA),
             file.path(reportDir,"paired_airway_original_response.json"),useBytes=TRUE)

}
cat("PUBLIC AIRWAY PAIRED BACKEND: PASS;",
    nrow(observed),"genes; four paired cell lines; CRISPLD2 positive;",
    "complete CSV and HTTP JSON; native DESeq2 matched\n")
