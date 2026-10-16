#!/usr/bin/env Rscript
# Genuine PUBLIC biological benchmark: airway / Himes et al., GEO GSE52778.
# Four primary airway smooth-muscle cell lines, each observed untreated and
# dexamethasone-treated (8 RNA-seq samples, paired within cell line).
#
# This benchmark tests the R method wrapper and COMPLETE export against an
# independently configured reference DESeq2 invocation, with a predeclared
# known glucocorticoid-responsive marker (CRISPLD2 / ENSG00000103196).
# It does NOT constitute independent statistical algorithm certification.
suppressPackageStartupMessages({
  library(airway)
  library(DESeq2)
})
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("jsonlite is required for the scientific reproducibility report")
source(file.path("multiomics-engine","advanced_methods.R"))
source(file.path("multiomics-engine","reference_full_results.R"))
set.seed(20261013)

data("airway",package="airway")
source_counts <- SummarizedExperiment::assay(airway)
stopifnot(ncol(source_counts)==8L, nrow(source_counts)>50000L)
subjects <- colnames(source_counts)
cell <- factor(as.character(airway$cell))
condition <- factor(as.character(airway$dex),levels=c("untrt","trt"))
stopifnot(identical(levels(condition),c("untrt","trt")),
          identical(as.integer(table(cell)),rep(2L,4L)),
          identical(as.integer(table(condition)),c(4L,4L)))
metadata <- data.frame(condition=condition,cell=cell,row.names=subjects)
# Prespecified expression threshold is the SAME documented wrapper rule:
# >=10 reads in >= max(2, ceiling(20% of samples)) (i.e., 2/8).
# Applying it here first keeps the replay reasonably sized; it is not a
# data-driven significance filter and should not be concealed in papers.
preselected <- rowSums(source_counts>=10)>=max(2L,ceiling(ncol(source_counts)*.2))
counts <- source_counts[preselected,,drop=FALSE]
stopifnot(nrow(counts)>10000L,nrow(counts)<nrow(source_counts),
          "ENSG00000103196"%in%rownames(counts))
matrix <- t(counts)
stopifnot(identical(rownames(matrix),rownames(metadata)))
formula <- stats::as.formula("~ cell + condition")
contrast <- c("condition","trt","untrt")
stopifnot(qr(stats::model.matrix(formula,data=metadata))$rank==
            ncol(stats::model.matrix(formula,data=metadata)))

scratch <- tempfile("pmx_public_airway_")
dir.create(scratch)
tryCatch({
  wrapper <- run_deseq2_counts(matrix,metadata,file.path(scratch,"wrapper"),
    design_formula=formula,contrast=contrast)
  stopifnot(nrow(wrapper$results)>10000L,
            all(c("feature","log2FoldChange","pvalue","padj")%in%names(wrapper$results)))

  # An INDEPENDENT call to Bioconductor's documented method. Use exactly the
  # same design, contrast and prefilter; do not compare to an unpaired ~dex
  # fit, nor to q-values adjusted over a different number of tested features.
  dds <- DESeq2::DESeqDataSetFromMatrix(countData=counts,colData=metadata,design=formula)
  dds <- DESeq2::DESeq(dds,quiet=TRUE)
  direct <- DESeq2::results(dds,contrast=contrast,alpha=.05)
  df <- data.frame(feature=rownames(direct),as.data.frame(direct),
                   row.names=NULL,check.names=FALSE)
  a <- wrapper$results
  b <- df[match(a$feature,df$feature),,drop=FALSE]
  stopifnot(identical(a$feature,b$feature),nrow(a)==nrow(b))
  fields <- c("baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
  numerical_max_abs <- list()
  for(field in fields){
    x <- a[[field]];y <- b[[field]]
    stopifnot(identical(which(is.na(x)),which(is.na(y))))
    valid <- which(is.finite(x)&is.finite(y))
    stopifnot(length(valid)>1000L,
              isTRUE(all.equal(x[valid],y[valid],tolerance=1e-8)))
    numerical_max_abs[[field]] <- if(length(valid))max(abs(x[valid]-y[valid])) else NA_real_
  }
  marker <- "ENSG00000103196"
  marker_row <- a[a$feature==marker,,drop=FALSE]
  stopifnot(nrow(marker_row)==1L,
            is.finite(marker_row$log2FoldChange),
            marker_row$log2FoldChange>0,
            is.finite(marker_row$padj),marker_row$padj<.05)

  export <- reference_complete_table(a,original_features=nrow(source_counts))
  stopifnot(isTRUE(export$referenceTableComplete),
            identical(export$referenceTableStatus,"complete"),
            export$referenceTableRows==nrow(a),
            export$referenceTableExportedRows==nrow(a),
            export$referenceTableOriginalFeatures==nrow(source_counts),
            export$referenceTableExcludedBeforeFit==nrow(source_counts)-nrow(a))
  reloaded <- utils::read.csv(text=export$referenceTableCsv,check.names=FALSE,
                        stringsAsFactors=FALSE,na.strings="NA")
  stopifnot(nrow(reloaded)==nrow(a),
            identical(reloaded$feature,a$feature),
            identical(which(is.na(reloaded$padj)),which(is.na(a$padj))),
            isTRUE(all.equal(reloaded$padj,a$padj,tolerance=1e-12)),
            isTRUE(all.equal(reloaded$log2FoldChange,a$log2FoldChange,tolerance=1e-12)))
  rwire <- jsonlite::toJSON(export,auto_unbox=TRUE,na="null")
  parsed <- jsonlite::fromJSON(rwire,simplifyVector=FALSE)
  stopifnot(isTRUE(parsed$referenceTableComplete),
            parsed$referenceTableRows==nrow(a),
            identical(parsed$referenceTableCsv,export$referenceTableCsv))

  n_significant <- sum(is.finite(a$padj)&a$padj<.05)
  report <- list(
    benchmark="Airway public 8-sample paired dexamethasone RNA-seq",
    accession="GSE52778",
    citation="Himes et al. 2014, PLOS ONE 9:e99625, PMID 24926665",
    source_package="Bioconductor airway",
    source_version=as.character(utils::packageVersion("airway")),
    DESeq2_version=as.character(utils::packageVersion("DESeq2")),
    R_version=R.version.string,
    statistical_design="~ cell + condition",
    contrast=as.list(contrast),
    independent_biological_units=4L,
    samples=8L,
    samples_per_condition=list(untrt=4L,trt=4L),
    raw_input_genes=nrow(source_counts),
    prespecified_prefilter="counts>=10 in >=2 of 8 samples",
    fitted_genes=nrow(a),
    excluded_genes=nrow(source_counts)-nrow(a),
    complete_reference_rows=nrow(reloaded),
    crosscheck="wrapper vs separate native DESeq2 invocation; same method and model",
    comparison_max_absolute_difference=numerical_max_abs,
    expected_marker=list(
      symbol="CRISPLD2",ensembl=marker,
      treatment_log2_fold_change=unname(marker_row$log2FoldChange),
      adjusted_p=unname(marker_row$padj),direction="higher with dexamethasone"
    ),
    adjusted_p_below_0_05=n_significant,
    evidence_scope="data/model integration, sign, numerical reproducibility, completeness",
    statistical_certification=FALSE,
    clinical_validation=FALSE,
    cross_method_validation=FALSE,
    independent_cohort_validation=FALSE
  )
  output <- Sys.getenv("PMX_AIRWAY_REPORT_DIR","")
  if(nzchar(output)){
    dir.create(output,recursive=TRUE,showWarnings=FALSE)
    jsonlite::write_json(report,file.path(output,"airway_public_reference_replay.json"),
                         pretty=TRUE,auto_unbox=TRUE,null="null")
  }
  cat("PUBLIC AIRWAY BIOLOGICAL REPLAY: PASS; paired ~cell+condition; ",
      nrow(a),"/",nrow(source_counts)," genes; ",
      n_significant," FDR < 0.05; CRISPLD2 log2FC=",
      format(marker_row$log2FoldChange,digits=4),"; q=",
      format(marker_row$padj,digits=4),"; CSV+JSON complete\n",sep="")
},finally={
  unlink(scratch,recursive=TRUE,force=TRUE)
})
