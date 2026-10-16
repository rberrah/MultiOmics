#!/usr/bin/env Rscript
# Base R only — does not silently mock an installed Bioconductor method.
source(file.path("multiomics-engine","reference_full_results.R"))
expect_error<-function(expr,pattern){
  e<-tryCatch({force(expr);NULL},error=function(e)conditionMessage(e))
  stopifnot(is.character(e),grepl(pattern,e))
}
n<-257L
f<-data.frame(
 feature=c("GENE,quoted","GENE\"quoted","GENE\nnewline",
           sprintf("GENE_%04d",seq_len(n-3L))),
 log2FoldChange=seq(-2,2,length.out=n),
 pvalue=seq(0,1,length.out=n),padj=p.adjust(seq(0,1,length.out=n),"BH"),
 stringsAsFactors=FALSE,check.names=FALSE
)
f$pvalue[c(10,20)]<-NA_real_
f$padj[c(11,21)]<-NA_real_
full<-reference_complete_table(f,original_features=300L)
stopifnot(full$referenceTableComplete,full$referenceTableStatus=="complete",
          full$referenceTableRows==257L,full$referenceTableExportedRows==257L,
          full$referenceTableOriginalFeatures==300L,
          full$referenceTableExcludedBeforeFit==43L,
          identical(unlist(full$referenceTableColumns),names(f)),
          length(strsplit(full$referenceTableCsv,"\n",fixed=TRUE)[[1]])>250L)
again<-utils::read.csv(text=full$referenceTableCsv,
                       check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(again)==nrow(f),identical(again$feature,f$feature),
          isTRUE(all.equal(again$log2FoldChange,f$log2FoldChange,tolerance=1e-9)),
          identical(which(is.na(again$padj)),which(is.na(f$padj))),
          identical(which(is.na(again$pvalue)),which(is.na(f$pvalue))),
          isTRUE(all.equal(again$padj[!is.na(again$padj)],
                           f$padj[!is.na(f$padj)],tolerance=1e-12)))
# Even if the native UI preview stops after 100, the reference export is full.
stopifnot(nrow(again)>100L && any(again$feature=="GENE_0254"))
small<-reference_complete_table(f,max_wire_bytes=300L)
stopifnot(identical(small$referenceTableStatus,"too_large"),
          !isTRUE(small$referenceTableComplete),
          is.null(small$referenceTableCsv),
          small$referenceTableRows==nrow(f))
d<-f;d$feature[2]<-""
expect_error(reference_complete_table(d),"missing")
d<-f;d$feature[2]<-d$feature[1]
expect_error(reference_complete_table(d),"duplicated")
d<-f;d$feature[2]<-NA_character_
expect_error(reference_complete_table(d),"missing")
expect_error(reference_complete_table(f,original_features=100L),
             "more feature rows")
if(requireNamespace("jsonlite",quietly=TRUE)){
  wire<-jsonlite::toJSON(full,auto_unbox=TRUE,na="null")
  parsed<-jsonlite::fromJSON(wire,simplifyVector=FALSE)
  stopifnot(isTRUE(parsed$referenceTableComplete),
            identical(parsed$referenceTableCsv,full$referenceTableCsv))
}
cat("REFERENCE FULL RESULTS: PASS, 257/257 rows, nonstandard identifiers, exact NA/q, payload failure visible\n")
