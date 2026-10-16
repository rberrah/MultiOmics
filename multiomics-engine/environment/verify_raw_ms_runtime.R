#!/usr/bin/env Rscript
# Fail closed on native R/Bioconductor/MS package and transitive version drift.
# Every version and the full hash were measured from actual real faahKO CI.
args<-commandArgs(trailingOnly=TRUE)
mode<-if(length(args)) args[[1L]] else "verify"
if(!mode%in%c("verify","probe"))stop("Mode must be verify or probe")
if(!requireNamespace("jsonlite",quietly=TRUE)||!requireNamespace("digest",quietly=TRUE))
  stop("jsonlite and digest are required to verify genuine native MS runtime")
lock<-jsonlite::read_json("multiomics-engine/environment/raw-ms-reference-runtime.json",
                         simplifyVector=TRUE)
if(!identical(lock$schema,"multiomics.raw-ms-reference-runtime.v1"))
  stop("Incorrect native MS runtime schema")
if(!identical(as.character(getRversion()),lock$rVersion))
  stop("R version drift, wanted ",lock$rVersion)
if(!requireNamespace("BiocManager",quietly=TRUE)||
   !identical(as.character(BiocManager::version()),lock$bioconductorVersion))
  stop("Bioconductor version drift, wanted ",lock$bioconductorVersion)
required<-unlist(lock$requiredPackageVersions,use.names=TRUE)
if(!length(required)||anyDuplicated(names(required)))
  stop("Invalid native MS package requirements")
missing<-names(required)[!vapply(names(required),requireNamespace,logical(1),quietly=TRUE)]
if(length(missing))stop("Missing native MS packages: ",paste(missing,collapse=", "))
actual<-vapply(names(required),function(p)as.character(utils::packageVersion(p)),character(1))
installed<-utils::installed.packages()
names_all<-installed[,"Package"]
if(anyDuplicated(names_all)||!length(names_all))
  stop("Duplicate or missing installed R package names")
versions_all<-setNames(as.character(installed[,"Version"]),names_all)
ordered<-sort(names(versions_all),method="radix")
canonical<-paste0(paste0(ordered,"=",unname(versions_all[ordered]),collapse="\n"),"\n")
fingerprint<-digest::digest(canonical,algo="sha256",serialize=FALSE)
out<-Sys.getenv("MULTIOMICS_RAW_MS_OUTPUT","tmp/multiomics/faahko-xcms-evidence")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
jsonlite::write_json(list(
  schema="multiomics.raw-ms-actual-runtime.v1",mode=mode,
  r_version=as.character(getRversion()),
  bioconductor_version=as.character(BiocManager::version()),
  image=lock$container,
  direct_package_versions=as.list(actual),
  installed_package_count=length(versions_all),
  installed_package_versions=as.list(versions_all),
  full_installed_library_sha256=fingerprint,
  source_sha=Sys.getenv("GITHUB_SHA","local"),
  strict_check_executed=identical(mode,"verify"),
  disclaimer="Software and package version replay only; does not validate instrument parameters, metabolite identification or statistical FDR."
),file.path(out,"raw_ms_frozen_environment.json"),pretty=TRUE,auto_unbox=TRUE)
if(identical(mode,"verify")){
  if(!isTRUE(lock$strictLockReady))stop("MS package lock is not ready")
  changed<-names(required)[actual[names(required)]!=unname(required)]
  if(length(changed))stop("Native MS direct version drift: ",paste(changed,collapse=", "))
  if(!identical(as.integer(length(versions_all)),
                as.integer(lock$fullInstalledLibrary$packageCount))||
     !identical(fingerprint,lock$fullInstalledLibrary$sha256SortedPackageVersionLines))
    stop("Native MS transitive dependency drift: ",length(versions_all),
         " packages; SHA-256=",fingerprint)
}
cat("PASS genuine native MS ",mode," | R ",as.character(getRversion()),
    " | Bioconductor ",BiocManager::version(),
    " | packages ",length(versions_all),
    " | sha256 ",fingerprint,"\n",sep="")
