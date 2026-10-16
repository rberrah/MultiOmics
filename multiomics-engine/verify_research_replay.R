#!/usr/bin/env Rscript
# Independent, strict replay check for the reference R API (no Python needed).
# Usage: Rscript reproduction/verify_replay.R reproduction/
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L)stop("Usage: Rscript verify_replay.R REPRODUCTION_DIRECTORY")
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("jsonlite required to verify original native R results")
root <- normalizePath(args[[1]],mustWork=TRUE)
original <- jsonlite::fromJSON(file.path(root,"reference_expected.json"),simplifyVector=FALSE)
replayed <- jsonlite::fromJSON(file.path(root,"reference_replayed.json"),simplifyVector=FALSE)
errors <- character()
add_error <- function(message)errors <<- c(errors,message)
same <- function(a,b,path,rtol=1e-7,atol=1e-9) {
  if(is.null(a)||is.null(b)){
    if(!is.null(a)||!is.null(b))add_error(paste0(path,": missing/null mismatch"))
    return(invisible(NULL))
  }
  if(is.list(a)&&is.list(b)){
    # Named JSON objects require the exact same names and values. Arrays
    # preserve sequence and length. Names are not a feature identifier.
    named_a <- !is.null(names(a));named_b <- !is.null(names(b))
    if(named_a!=named_b){
      add_error(paste0(path,": object/array mismatch"))
    } else if(named_a){
      if(!setequal(names(a),names(b)))
        add_error(paste0(path,": missing/extra keys: ",paste(setdiff(union(names(a),names(b)),intersect(names(a),names(b))),collapse=", ")))
      for(name in intersect(names(a),names(b)))
        same(a[[name]],b[[name]],paste0(path,".",name),rtol,atol)
    }else{
      if(length(a)!=length(b))add_error(paste0(path,": array size ",length(a)," != ",length(b)))
      for(i in seq_len(min(length(a),length(b))))
        same(a[[i]],b[[i]],paste0(path,"[",i,"]"),rtol,atol)
    }
    return(invisible(NULL))
  }
  if(is.numeric(a)&&is.numeric(b) && length(a)==1L&&length(b)==1L){
    if(!is.finite(a)||!is.finite(b)||abs(a-b)>atol+rtol*max(abs(a),abs(b)))
      add_error(paste0(path,": numeric mismatch ",a," vs ",b))
    return(invisible(NULL))
  }
  if(!identical(a,b))add_error(paste0(path,": value/type mismatch"))
  invisible(NULL)
}
methods <- original$methods
if(is.null(methods)||!is.list(methods)||
   !any(vapply(methods,function(m)identical(m$status,"ok"),logical(1))))
  add_error("No successful original native R method")
if(is.null(replayed$methods)||!is.list(replayed$methods))
  add_error("No replayed R methods")
same(original$methods,replayed$methods,"methods")
same(original$preprocessing,replayed$preprocessing,"preprocessing")
for(field in c("analysisMode","analysedLayers"))
  same(original[[field]],replayed[[field]],field)
report <- list(
  scope="native_R_reference_methods_and_preprocessing_only",
  equivalence_if_pass="within_numeric_tolerance",
  relative_tolerance=1e-7,absolute_tolerance=1e-9,
  browser_exploratory_statistics_included=FALSE,
  browser_figures_included=FALSE,
  n_discrepancies=length(errors),
  discrepancies=as.list(head(errors,30)),
  result=if(length(errors))"FAIL" else "PASS"
)
writeLines(jsonlite::toJSON(report,auto_unbox=TRUE,pretty=TRUE,null="null"),
           file.path(root,"r_numeric_comparison_report.json"))
if(length(errors)){
  for(message in head(errors,30))cat("FAIL: ",message,"\n",sep="",file=stderr())
  stop("Native R numerical replay mismatched in ",length(errors)," places.")
}
cat("PASS: Native R methods, complete tables and preprocessing match original backend response.\n")
cat("This does not certify JS-only statistics or pixel-identical figures.\n")
