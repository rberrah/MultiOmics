#!/usr/bin/env Rscript
# Profile-specific reproducibility lock for authentic Nutrimouse methods.
args<-commandArgs(trailingOnly=TRUE)
if(length(args)<1L||!args[[1L]]%in%c("diablo","mofa"))
 stop("Usage: Rscript verify_nutrimouse_runtime.R diablo|mofa [probe|verify]")
profile<-args[[1L]]
mode<-if(length(args)>1L)args[[2L]]else "verify"
if(!mode%in%c("probe","verify"))stop("Mode must be probe or verify")
if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is required")
lock<-jsonlite::fromJSON("multiomics-engine/environment/nutrimouse-package-pins.json",
 simplifyVector=FALSE)
if(!identical(lock$schema,"multiomics.nutrimouse-r-lock.v1"))
 stop("Unknown Nutrimouse reference environment schema")
if(!identical(as.character(getRversion()),lock$rVersion))
 stop("R drift in Nutrimouse profile")
if(!requireNamespace("BiocManager",quietly=TRUE)||
   !identical(as.character(BiocManager::version()),lock$bioconductorVersion))
 stop("Bioconductor release drift in Nutrimouse profile")
p<-lock$profiles[[profile]]
pkgs<-unlist(p$requiredPackages,use.names=FALSE)
if(!length(pkgs)||anyDuplicated(pkgs))stop("Invalid Nutrimouse package list")
missing<-pkgs[!vapply(pkgs,requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop("Missing actual R methods: ",paste(missing,collapse=", "))
actual<-vapply(pkgs,function(x)as.character(utils::packageVersion(x)),character(1))
if(identical(mode,"verify")){
 if(!isTRUE(p$strictLockReady))stop("Method package version lock pending, cannot publish")
 target<-unlist(p$versionLock,use.names=TRUE)
 if(!setequal(names(target),pkgs)||
    any(actual[names(target)]!=target))
   stop("Nutrimouse ",profile," dependency drift: ",
        paste(names(target)[actual[names(target)]!=target],collapse=", "))
}
destination<-if(profile=="diablo")
 "tmp/multiomics/real-nutrimouse-diablo" else
 "tmp/multiomics/mofa-nutrimouse-real"
dir.create(destination,recursive=TRUE,showWarnings=FALSE)
installed<-utils::installed.packages()
report<-list(profile=profile,mode=mode,r=as.character(getRversion()),
 bioconductor=as.character(BiocManager::version()),
 image=lock$container,
 selected_package_versions=as.list(actual),
 all_installed_package_versions=setNames(as.list(installed[,"Version"]),installed[,"Package"]),
 github_sha=Sys.getenv("GITHUB_SHA","local"),
 exact_package_lock_verified=identical(mode,"verify")&&isTRUE(p$strictLockReady),
 caution="Software reproducibility only; neither external biomarker validation nor universal model superiority")
jsonlite::write_json(report,file.path(destination,"reference_environment.json"),
 pretty=TRUE,auto_unbox=TRUE)
cat("MULTIOMICS_ACTUAL_",toupper(profile),"_VERSIONS=",
    jsonlite::toJSON(as.list(actual),auto_unbox=TRUE),"\n",sep="")
cat("Nutrimouse ",profile," ",mode," environment PASS\n",sep="")
