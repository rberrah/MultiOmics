#!/usr/bin/env Rscript
# Local, read-only onboarding diagnostic for MultiOmics reference R methods.
# This DOES NOT install packages, start a service, certify a study, or
# imply that observed package versions match the frozen research image.

profiles <- list(
  rnaseq = list(required=c("plumber","jsonlite"),
                alternatives=list(c("DESeq2"),c("edgeR","limma"))),
  groups = list(required=c("plumber","jsonlite","limma"),
                alternatives=list()),
  longitudinal = list(required=c("plumber","jsonlite","lmerTest"),
                      alternatives=list()),
  mofa2 = list(required=c("plumber","jsonlite","MOFA2"),
               alternatives=list()),
  diablo = list(required=c("plumber","jsonlite","mixOmics"),
                alternatives=list()),
  rawms = list(required=c("plumber","jsonlite","xcms",
                          "MsExperiment","Spectra","mzR"),
               alternatives=list())
)

readiness <- function(profile, installed) {
  spec <- profiles[[profile]]
  if(is.null(spec)) stop("Unknown reference profile: ",profile)
  required_missing <- spec$required[!vapply(spec$required,function(pkg)
    isTRUE(installed[[pkg]]),logical(1))]
  alternative_ok <- !length(spec$alternatives) ||
    any(vapply(spec$alternatives,function(group)
      all(vapply(group,function(pkg)isTRUE(installed[[pkg]]),logical(1))),
      logical(1)))
  optional_groups <- if(length(spec$alternatives))
    vapply(spec$alternatives,paste,collapse=" + ",character(1)) else character()
  list(ok=!length(required_missing)&&alternative_ok,
       missing=required_missing,
       requires_one_of=optional_groups,
       alternatives_ok=alternative_ok)
}

self_test <- function() {
  installed <- setNames(rep(TRUE,16),c(
    "plumber","jsonlite","DESeq2","edgeR","limma","lmerTest","MOFA2",
    "mixOmics","xcms","MsExperiment","Spectra","mzR","survival",
    "BiocParallel","fgsea","BiocManager"))
  stopifnot(all(vapply(names(profiles),function(n)
    readiness(n,installed)$ok,logical(1))))
  absent <- installed; absent[["DESeq2"]] <- FALSE
  stopifnot(readiness("rnaseq",absent)$ok)
  absent[["edgeR"]] <- FALSE
  stopifnot(!readiness("rnaseq",absent)$ok)
  absent <- installed; absent[["jsonlite"]] <- FALSE
  stopifnot(!readiness("mofa2",absent)$ok)
  absent <- installed; absent[["xcms"]] <- FALSE
  stopifnot(!readiness("rawms",absent)$ok)
  tryCatch({readiness("unknown",installed);stop("Accepted unknown profile")},
           error=function(e)stopifnot(grepl("Unknown reference profile",conditionMessage(e))))
  cat("PASS: doctor profile selection, fallback packages and fail-closed checks\n")
}

args <- commandArgs(trailingOnly=TRUE)
if(identical(args,"--self-test")) {
  self_test()
  quit(save="no",status=0L)
}
if(length(args)>1L || (length(args)==1L && !args[[1L]] %in% names(profiles))) {
  stop("Usage: Rscript multiomics-engine/reference_doctor.R [rnaseq|groups|longitudinal|mofa2|diablo|rawms|--self-test]",call.=FALSE)
}
profile <- if(length(args)) args[[1L]] else "rnaseq"
all_packages <- unique(c(unlist(lapply(profiles,function(x)x$required)),
                         unlist(lapply(profiles,function(x)x$alternatives))))
availability <- setNames(vapply(all_packages,function(pkg)
  requireNamespace(pkg,quietly=TRUE),logical(1)),all_packages)
report <- readiness(profile,availability)
cat("MultiOmics R reference diagnostic (local; no data uploaded)\n")
cat("Profile: ",profile,"\n",sep="")
cat("R version: ",as.character(getRversion())," (",R.version$platform,")\n",sep="")
for(pkg in profiles[[profile]]$required) {
  cat("  ",pkg,": ",if(availability[[pkg]]) "installed" else "MISSING","\n",sep="")
}
if(length(report$requires_one_of)) {
  cat("Alternative scientific methods (at least one group needed): ",
      paste(report$requires_one_of,collapse=" OR "),"\n",sep="")
  for(group in profiles[[profile]]$alternatives)
    for(pkg in group)
      cat("  ",pkg,": ",if(availability[[pkg]]) "installed" else "MISSING","\n",sep="")
}
if(!report$ok) {
  cat("NOT READY: missing packages or no applicable reference model.\n")
  cat("Install dependencies from the repository root with:\n")
  cat("  Rscript multiomics-engine/install_backend_dependencies.R\n")
  cat("Then rerun this diagnostic; packages are never installed automatically by the doctor.\n")
  quit(save="no",status=1L)
}
cat("READY TO START backend for this method family, NOT scientific validation.\n")
cat("  Rscript multiomics-engine/run_backend.R\n")
cat("Then in MultiOmics choose R required, leave http://127.0.0.1:8787 and click Check.\n")
cat("For frozen-environment work also run:\n")
cat("  Rscript multiomics-engine/install_backend_dependencies.R --verify-frozen\n")
quit(save="no",status=0L)
