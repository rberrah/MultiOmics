#!/usr/bin/env Rscript
# TRUE server-level paired RNA-seq: 4 independent patients x 2 assays, DESeq2
# and limma-voom with a subject block. All failures must fail CLOSED before
# publishing any apparent "independent groups" scientific result.
source(file.path("multiomics-engine","server.R"))
for(pkg in c("DESeq2","edgeR","limma","jsonlite","statmod"))
  stopifnot(requireNamespace(pkg,quietly=TRUE))
set.seed(1301026)
pairs<-4L; n<-2L*pairs; genes<-156L
ids<-paste0("S",seq_len(n))
subj<-rep(paste0("SUBJECT",seq_len(pairs)),each=2)
group<-rep(c("control","treated"),times=pairs)
count_matrix<-vapply(seq_len(genes),function(j){
  effect<-if(j<=21)ifelse(group=="treated",2.0,1) else 1
  stats::rnbinom(n,mu=(63+j%%20)*effect,size=11)
},numeric(n))
dimnames(count_matrix)<-list(ids,paste0("GENE",seq_len(genes)))
csv<-function(x)paste0(paste(capture.output(
  utils::write.csv(x,row.names=FALSE,quote=TRUE)),collapse="\n"),"\n")
metadata<-data.frame(
  subject_id=subj,sample_id=ids,assay_id=ids,omic="transcriptomics",
  condition=group,timepoint="T0",batch="RUN1",
  technical_replicate="1",sample_type="biological",stringsAsFactors=FALSE
)
matrix_df<-data.frame(feature_id=colnames(count_matrix),t(count_matrix),
                      check.names=FALSE)
payload<-function(m=metadata,x=matrix_df,mode="exploratory"){
  list(protocol=list(objective="groups",designType="paired",
    analysisIntent=mode,longitudinal=FALSE,rnaCountMethod="both",
    rnaPrimaryContrast=as.list(c("condition","treated","control")),
    covariateColumns=list()),
    dataTypes=list(transcriptomics="raw_counts"),
    identifierTypes=list(transcriptomics="ensembl_gene"),
    columnMapping=as.list(stats::setNames(names(m),names(m))),
    metadataCsv=csv(m),matrices=list(transcriptomics=csv(x)))
}
expect_error<-function(expr,pattern){
  msg<-tryCatch({force(expr);NULL},error=function(e)conditionMessage(e))
  if(!is.character(msg)||!grepl(pattern,msg,ignore.case=TRUE))
    stop("Expected guarded error /",pattern,"/, received ",as.character(msg))
}
r<-run_backend_analysis(payload(mode="confirmatory"))
stopifnot(identical(r$status,"ok"))
d<-r$methods$transcriptomics_differential_deseq2
v<-r$methods$transcriptomics_differential_voom
stopifnot(d$status=="ok",v$status=="ok",
          isTRUE(d$referenceTableComplete),isTRUE(v$referenceTableComplete),
          d$referenceTableRows>100,v$referenceTableRows>100,
          d$referenceTableRows==d$referenceTableExportedRows,
          v$referenceTableRows==v$referenceTableExportedRows,
          length(d$top)==100L,length(v$top)==100L)
guard<-d$pairedDesignContract
stopifnot(guard$independent_biological_subjects==4L,
          guard$biological_pairs==4L,guard$model_samples==8L,
          guard$residual_degrees_of_freedom==3L,
          guard$contrast_was_prespecified,
          identical(unlist(guard$explicit_primary_contrast),
                    c("condition","treated","control")),
          !isTRUE(guard$certified))
stopifnot(isTRUE(d$assayIntegrityContract$assay_loss_detected==FALSE))
de <- utils::read.csv(text=d$referenceTableCsv,check.names=FALSE,stringsAsFactors=FALSE)
vo <- utils::read.csv(text=v$referenceTableCsv,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(nrow(de)==d$referenceTableRows,nrow(vo)==v$referenceTableRows,
          all(c("pvalue","padj","log2FoldChange")%in%names(de)),
          all(c("P.Value","adj.P.Val","logFC")%in%names(vo)),
          length(unique(de$feature))==nrow(de))
# Verify against a completely separate direct DESeq2 fit with the CORRECT
# paired design. This guards against accidental independent-condition fits.
col<-data.frame(subject_id=factor(subj),condition=factor(group,
  levels=c("control","treated")),row.names=ids)
dds<-DESeq2::DESeqDataSetFromMatrix(
  countData=t(count_matrix),colData=col,design=~subject_id+condition)
keep<-rowSums(DESeq2::counts(dds)>=10)>=2
dds<-DESeq2::DESeq(dds[keep,],quiet=TRUE)
direct<-DESeq2::results(dds,contrast=c("condition","treated","control"),alpha=.05)
compare<-as.data.frame(direct)[match(de$feature,rownames(direct)),,drop=FALSE]
stopifnot(identical(de$feature,rownames(compare)))
for(k in c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")){
  stopifnot(identical(which(is.na(de[[k]])),which(is.na(compare[[k]]))),
            isTRUE(all.equal(de[[k]],compare[[k]],tolerance=1e-7)))
}
wire<-jsonlite::toJSON(r,auto_unbox=TRUE,na="null",dataframe="rows")
received<-jsonlite::fromJSON(wire,simplifyVector=FALSE)
stopifnot(isTRUE(received$methods$transcriptomics_differential_deseq2$referenceTableComplete),
          received$methods$transcriptomics_differential_deseq2$pairedDesignContract$biological_pairs==4L)
cat("PAIRED RNA-SEQ REAL: PASS; 4 pairs, 8 assays; DESeq2/voom full tables, independent native ~subject+condition replay\n")

# Missing second condition for a subject, no automatic dropping.
bad<-metadata[-nrow(metadata),,drop=FALSE]
expect_error(run_backend_analysis(payload(m=bad,x=matrix_df)),
             "assay IDs|exactly one|pairs")
# Duplicated biological pairing cell.
bad<-metadata;bad$condition[2]<-"control"
expect_error(run_backend_analysis(payload(m=bad)), "exactly one|paired")
# One subject per assay is NOT a valid matched design.
bad<-metadata;bad$subject_id<-paste0("UNPAIRED",seq_len(n))
expect_error(run_backend_analysis(payload(m=bad)), "independent subjects|exactly one")
# Only three pairs — no silent crossing below the minimum guard.
bad<-metadata[1:6,,drop=FALSE]
x<-matrix_df[,c(1,match(bad$assay_id,names(matrix_df))),drop=FALSE]
expect_error(run_backend_analysis(payload(m=bad,x=x)),">=4")
# Batch cannot be fitted if perfectly aliased with subject blocks.
bad<-metadata;bad$batch<-rep(paste0("B",seq_len(pairs)),each=2)
expect_error(run_backend_analysis(payload(m=bad)), "confounded|design")
# Missing prespecified contrast blocks the confirmatory route.
unplanned<-payload(mode="confirmatory")
unplanned$protocol$rnaPrimaryContrast<-NULL
expect_error(run_backend_analysis(unplanned),"explicit|contrast")
# Technical/data integrity guards should never round or drop fractional counts.
badmatrix<-matrix_df;badmatrix[1,2]<-as.numeric(badmatrix[1,2])+.5
expect_error(run_backend_analysis(payload(x=badmatrix)),"fractional|integer")
# A paired request with normalized values cannot be silently fitted as counts.
wrong<-payload();wrong$dataTypes$transcriptomics<-"normalized"
expect_error(run_backend_analysis(wrong),"raw-count")
cat("PAIRED RNA-SEQ NEGATIVE DESIGNS: PASS; incomplete, duplicate, unpaired, low-N, confounded, no contrast, noninteger, wrong value type\n")
