#!/usr/bin/env Rscript
# Security and provenance regression for a RESEARCHER-OPTED-IN private
# directory on the same local machine running the R reference backend.
# Big scientific tables may exceed the HTTP JSON budget, but must never be
# silently reported as complete or exposed through a public web route.
source(file.path("multiomics-engine","reference_full_results.R"))
stopifnot(requireNamespace("digest",quietly=TRUE))
root<-tempfile("pmx_private_")
dir.create(root,mode="0700")
tryCatch({
  n<-420L
  features<-c("=HYPERLINK(\"https://evil.example\")",
              "ENSG,quoted","Gene\nmultiline",paste0("GENE_",seq_len(n-3L)))
  tab<-data.frame(feature=features,log2FoldChange=seq(-3,3,length.out=n),
                  pvalue=seq_len(n)/(n+1),
                  padj=stats::p.adjust(seq_len(n)/(n+1),"BH"),
                  check.names=FALSE,stringsAsFactors=FALSE)
  tab$padj[c(3,30)]<-NA_real_
  # No archive configured: MUST make it crystal clear there is no file.
  unretained<-reference_complete_table(tab,max_wire_bytes=300L,
                                        local_archive_dir="")
  stopifnot(identical(unretained$referenceTableStatus,"too_large"),
            !isTRUE(unretained$referenceTableComplete),
            is.null(unretained$referenceTableCsv),
            is.null(unretained$referenceTableLocalArchiveId),
            grepl("NOT retained",unretained$referenceTableMessage))
  # Explicit local opt-in: full result is preserved beyond the wire budget.
  saved<-reference_complete_table(tab,max_wire_bytes=300L,
                                  local_archive_dir=root)
  stopifnot(identical(saved$referenceTableStatus,"too_large_archived_locally"),
            !isTRUE(saved$referenceTableComplete),
            identical(saved$referenceTableRows,n),
            identical(saved$referenceTableExportedRows,0L),
            is.null(saved$referenceTableCsv),
            is.character(saved$referenceTableLocalArchiveId),
            !grepl("[/\\\\]",saved$referenceTableLocalArchiveId))
  retained<-file.path(root,saved$referenceTableLocalArchiveId,"fitted_reference.csv")
  stopifnot(file.exists(retained),identical(
    unname(file.info(retained)$size),saved$referenceTableLocalArchiveBytes),
    identical(digest::digest(file=retained,algo="sha256",serialize=FALSE),
              saved$referenceTableLocalArchiveSha256))
  reread<-utils::read.csv(retained,check.names=FALSE,stringsAsFactors=FALSE,
                          na.strings="NA")
  stopifnot(nrow(reread)==n,
            identical(reread$feature,features),
            identical(which(is.na(reread$padj)),which(is.na(tab$padj))),
            isTRUE(all.equal(reread$log2FoldChange,tab$log2FoldChange,
                             tolerance=1e-12)))
  # Don't retain excessive files if the researcher sets a small local quota.
  quota<-reference_complete_table(tab,max_wire_bytes=300L,
                  local_archive_dir=root,max_local_bytes=256L)
  stopifnot(identical(quota$referenceTableStatus,"too_large"),
            is.null(quota$referenceTableLocalArchiveId),
            grepl("NOT retained",quota$referenceTableMessage))
  invalid<-reference_complete_table(tab,max_wire_bytes=300L,
                  local_archive_dir=file.path(root,"not-present"))
  stopifnot(identical(invalid$referenceTableStatus,"too_large"),
            is.null(invalid$referenceTableLocalArchiveId))
  stopifnot(length(list.dirs(root,recursive=FALSE,full.names=TRUE))==1L)
  cat("PRIVATE OVERSIZED REFERENCE CSV: PASS, opt-in only, original IDs, NA/q, SHA256, quota and disabled cases\n")
},finally={unlink(root,recursive=TRUE,force=TRUE)})
