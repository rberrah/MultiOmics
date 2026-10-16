#!/usr/bin/env Rscript
# Frozen-environment replay for the RNA-seq reference methods.
# This verifies software-level numerical reproduction, NOT nominal FDR,
# biological discovery, model validity across populations, or clinical use.
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L || !args[[1L]] %in% c("capture","replay"))
  stop("Usage: Rscript rnaseq_replay.R capture|replay OUTPUT_DIRECTORY",call.=FALSE)
mode <- args[[1L]]
out <- normalizePath(args[[2L]],winslash="/",mustWork=FALSE)
dir.create(out,recursive=TRUE,showWarnings=FALSE)
required <- c("DESeq2","edgeR","limma","statmod","BiocManager","jsonlite","renv")
missing <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop(paste("Missing reference packages:",paste(missing,collapse=", ")))
stopifnot(startsWith(as.character(getRversion()),"4.5."),
          identical(as.character(BiocManager::version()),"3.22"))
versions <- setNames(vapply(required,function(p) as.character(utils::packageVersion(p)),character(1)),required)

source("multiomics-engine/advanced_methods.R")
study <- function(){
  set.seed(20261013L)
  n <- 24L; p <- 160L
  ids <- sprintf("S%03d",seq_len(n))
  features <- sprintf("GENE%04d",seq_len(p))
  condition <- factor(rep(c("control","treated"),each=n/2L),
                      levels=c("control","treated"))
  batch <- factor(rep(c("B1","B2"),length.out=n),levels=c("B1","B2"))
  meta <- data.frame(subject_id=sprintf("ANON%03d",seq_len(n)),
                     condition=condition,batch=batch,
                     row.names=ids,stringsAsFactors=FALSE)
  mu_base <- stats::runif(p,90,360)
  planted <- c(rep(log(1.6),12L),rep(0,p-12L))
  batch_effect <- c(rep(log(1.12),p/2L),rep(0,p/2L))
  libraries <- rep(c(0.85,1.05,1.2,0.92),length.out=n)
  counts <- matrix(0,n,p,dimnames=list(ids,features))
  for(i in seq_len(n)){
    mu <- mu_base*libraries[[i]]*exp(
      (as.character(condition[[i]])=="treated")*planted+
      (as.character(batch[[i]])=="B2")*batch_effect)
    counts[i,] <- stats::rnbinom(p,mu=mu,size=8)
  }
  list(counts=counts,metadata=meta,
       seed=20261013L,planted_genes=features[seq_len(12L)],
       note="Synthetic fixed truth, not a real clinical or biological benchmark")
}
input_file <- file.path(out,"synthetic_input.rds")
if(mode=="capture"){
  if(file.exists(input_file))stop("Refusing to overwrite prior reference input")
  input <- study()
  saveRDS(input,input_file,version=3)
}else{
  if(!file.exists(input_file))stop("No previously captured independent test input")
  input <- readRDS(input_file)
}
counts <- input$counts; metadata <- input$metadata
formula <- ~batch+condition
contrast <- c("condition","treated","control")
contract <- rnaseq_confirmatory_preflight(counts,metadata,formula,contrast)
stopifnot(identical(contract$status,"eligible_for_reference_fitting_only"),
          !isTRUE(contract$certified),contract$independent_subjects==24L)
bad <- counts;bad[1L,1L] <- bad[1L,1L]+0.5
stopifnot(inherits(try(rnaseq_confirmatory_preflight(
  bad,metadata,formula,contrast),silent=TRUE),"try-error"))

subdir <- tempfile(paste0("r1-",mode,"-"))
dir.create(subdir,recursive=TRUE)
deseq <- run_deseq2_counts(counts,metadata,
  file.path(subdir,"deseq2"),formula,contrast,confirmatory_contract=TRUE)
edger <- run_edger_ql_counts(counts,metadata,
  file.path(subdir,"edger"),formula,coefficient="conditiontreated")
voom <- run_voom_counts(counts,metadata,
  file.path(subdir,"voom"),formula,coefficient="conditiontreated")
extract <- function(df,columns){
  stopifnot(all(columns %in% colnames(df)))
  x <- df[,columns,drop=FALSE]
  x <- x[order(x$feature),,drop=FALSE]
  rownames(x) <- NULL
  if(!nrow(x))stop("No tested features")
  x
}
results <- list(
  DESeq2=extract(deseq$results,c("feature","log2FoldChange","pvalue","padj")),
  edgeR=extract(edger$results,c("feature","logFC","PValue","FDR")),
  limma_voom=extract(voom$results,c("feature","logFC","P.Value","adj.P.Val")))
reference_file <- file.path(out,"reference_results.rds")
if(mode=="capture"){
  if(file.exists(reference_file))stop("Reference results already exist")
  saveRDS(results,reference_file,version=3)
  renv::snapshot(project=getwd(),
    lockfile=file.path(out,"renv.lock"),
    packages=required,prompt=FALSE,force=FALSE,
    repos=BiocManager::repositories(version="3.22"))
  writeLines(capture.output(sessionInfo()),
             file.path(out,"sessionInfo.txt"))
  jsonlite::write_json(list(
    schema="pmx-rnaseq-reference-replay/v1",
    status="executed_in_pinned_bioconductor_container",
    r_version=as.character(getRversion()),
    bioconductor_version=as.character(BiocManager::version()),
    packages=as.list(versions),
    sample_count=nrow(counts),gene_count=ncol(counts),
    synthetic_truth=TRUE,
    confirmatory_input_eligible_not_certified=TRUE,
    clinical_validation=FALSE,universal_fdr_validation=FALSE,
    external_biological_validation=FALSE,
    limitations=c("One synthetic dataset; this is numerical replay only",
      "No FDR calibration or biological validation is inferred",
      "Lockfile was captured from the executed container; restore required to test portability")),
    file.path(out,"runtime_manifest.json"),auto_unbox=TRUE,pretty=TRUE)
  cat("REFERENCE CAPTURE OK: ",nrow(counts)," synthetic samples and ",
      ncol(counts)," genes; 3 real R reference methods\n",sep="")
}else{
  manifest <- jsonlite::read_json(file.path(out,"runtime_manifest.json"),simplifyVector=TRUE)
  stopifnot(identical(manifest$schema,"pmx-rnaseq-reference-replay/v1"),
    identical(as.character(getRversion()),manifest$r_version),
    identical(as.character(BiocManager::version()),manifest$bioconductor_version),
    !isTRUE(manifest$clinical_validation),
    !isTRUE(manifest$universal_fdr_validation))
  for(pkg in required)
    if(!identical(versions[[pkg]],manifest$packages[[pkg]]))
      stop(paste("PACKAGE DRIFT:",pkg,versions[[pkg]],manifest$packages[[pkg]]))
  isolated <- Sys.getenv("PMX_RESTORED_LIBRARY","")
  if(nzchar(isolated)) {
    if(!dir.exists(isolated))stop("Replayed R library does not exist")
    for(pkg in required) {
      # Requiring a package in the base image is NOT a successful restore.
      if(!file.exists(file.path(isolated,pkg,"DESCRIPTION")))
        stop(paste("PACKAGE NOT RESTORED IN PRIVATE LIBRARY:",pkg))
    }
  }
  previous <- readRDS(reference_file)
  for(k in names(results)){
    if(!identical(results[[k]]$feature,previous[[k]]$feature))
      stop(paste("FEATURE REPLAY DRIFT:",k))
    vars <- setdiff(colnames(results[[k]]),"feature")
    for(v in vars){
      a <- results[[k]][[v]];b <- previous[[k]][[v]]
      if(!isTRUE(all.equal(a,b,tolerance=1e-9,check.attributes=FALSE)))
        stop(paste("REFERENCE NUMERICAL DRIFT:",k,v))
    }
  }
  writeLines(paste("REPLAY PASS",Sys.time(),
    paste(names(results),collapse=",")),file.path(out,"replay-result.txt"))
  cat("INDEPENDENT REPLAY PASS: identical package versions and numerical results (tolerance 1e-9)\n")
}

unlink(subdir,recursive=TRUE)
