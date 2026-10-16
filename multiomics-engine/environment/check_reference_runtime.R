#!/usr/bin/env Rscript
# Truthful fail-closed reference R environment check. No fictitious package versions.
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[[1L]] else "verify"
if(!mode%in%c("probe","verify"))stop("Use probe or verify")
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite required for locked scientific environment")
config <- jsonlite::fromJSON("multiomics-engine/environment/reference-runtime.json",
                            simplifyVector=TRUE)
stopifnot(identical(config$schema,"multiomics.r-reference-runtime.v1"))
if(!identical(as.character(getRversion()),config$rVersion))
  stop("R drift: required ",config$rVersion,", current ",getRversion())
if(!requireNamespace("BiocManager",quietly=TRUE) ||
   !identical(as.character(BiocManager::version()),config$bioconductorVersion))
  stop("Bioconductor release drift: required ",config$bioconductorVersion)
packages <- as.character(config$requiredPackages)
if(anyDuplicated(packages)||!length(packages))stop("Invalid package lock manifest")
missing<-packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop("Missing reference R packages: ",paste(missing,collapse=", "))
versions<-vapply(packages,function(p)as.character(utils::packageVersion(p)),character(1))
if(identical(mode,"verify")){
  if(!isTRUE(config$strictLockReady))
    stop("Scientific release blocked: reference package versions have not been pinned")
  locked<-unlist(config$versionLock,use.names=TRUE)
  if(!setequal(names(locked),packages)||anyNA(locked)||
     any(!nzchar(locked)))
    stop("Reference package version lock is missing entries")
  if(!identical(sort(names(locked)),sort(names(versions)))||
      any(versions[names(locked)]!=locked)){
    incorrect<-names(locked)[versions[names(locked)]!=locked]
    stop("Scientific R package version drift: ",paste(incorrect,collapse=", "))
  }
}
dir.create("tmp/reference-runtime",recursive=TRUE,showWarnings=FALSE)
installed<-utils::installed.packages()
all_packages<-setNames(as.list(installed[,"Version"]),installed[,"Package"])
report<-list(
  schema="multiomics.r-reference-report.v1",
  mode=mode,
  r_version=as.character(getRversion()),
  bioconductor_version=as.character(BiocManager::version()),
  platform=R.version$platform,
  image=config$container,
  required_package_versions=as.list(versions),
  all_installed_versions=all_packages,
  strict_lock_verified=identical(mode,"verify") && isTRUE(config$strictLockReady),
  github_sha=Sys.getenv("GITHUB_SHA","local"),
  limitations=c("R package versions pinned for declared reference packages only",
     "Full OS and transitive dependency set is anchored by immutable container digest",
     "User datasets and clinical conclusions are not certified")
)
jsonlite::write_json(report,"tmp/reference-runtime/environment.json",pretty=TRUE,auto_unbox=TRUE)
cat("MULTIOMICS_ACTUAL_REFERENCE_VERSIONS=",jsonlite::toJSON(as.list(versions),
 auto_unbox=TRUE),"\n",sep="")
cat("REFERENCE R ",mode,": R ",config$rVersion,", Bioconductor ",
    config$bioconductorVersion,", ",length(versions)," actual package versions\n",sep="")
