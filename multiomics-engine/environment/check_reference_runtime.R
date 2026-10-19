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
all_version_names<-installed[,"Package"]
if(anyDuplicated(all_version_names))stop("Ambiguous duplicate installed R package names")
all_versions<-setNames(as.character(installed[,"Version"]),all_version_names)
sorted_names<-sort(names(all_versions),method="radix")
canonical_inventory<-paste0(
  paste(paste0(sorted_names,"=",unname(all_versions[sorted_names])),collapse="\n"),
  "\n")
if(!requireNamespace("digest",quietly=TRUE))
  stop("digest required for SHA-256 of the actual R package closure")
inventory_hash<-digest::digest(canonical_inventory,algo="sha256",serialize=FALSE)
inventory_path <- "multiomics-engine/environment/reference-runtime-inventory.json"
if(identical(mode,"verify")){
  info<-config$fullInstalledLibrary
  if(!file.exists(inventory_path))
    stop("Missing transparent full R package inventory lock.")
  frozen_inventory<-jsonlite::fromJSON(inventory_path,simplifyVector=TRUE)
  if(!identical(frozen_inventory$schema,"multiomics.r-full-package-inventory.v1") ||
     !identical(as.integer(frozen_inventory$packageCount),
                as.integer(info$packageCount)) ||
     !identical(frozen_inventory$rVersion,config$rVersion) ||
     !identical(frozen_inventory$bioconductorVersion,config$bioconductorVersion) ||
     !identical(frozen_inventory$container,config$container))
    stop("Full R inventory metadata disagrees with reference runtime contract.")
  expected_all<-unlist(frozen_inventory$packages,use.names=TRUE)
  missing_all<-setdiff(names(expected_all),names(all_versions))
  extra_all<-setdiff(names(all_versions),names(expected_all))
  common_all<-intersect(names(expected_all),names(all_versions))
  changed_all<-common_all[expected_all[common_all]!=all_versions[common_all]]
  dir.create("tmp/reference-runtime",recursive=TRUE,showWarnings=FALSE)
  jsonlite::write_json(list(
    schema="multiomics.r-reference-drift-diagnostic.v1",
    expected_sha256=info$sha256SortedPackageVersionLines,
    observed_sha256=inventory_hash,
    expected_count=length(expected_all),
    observed_count=length(all_versions),
    missing=as.list(missing_all),
    extra=as.list(extra_all),
    changed=as.list(setNames(
      lapply(changed_all,function(p)list(expected=expected_all[[p]],observed=all_versions[[p]])),
      changed_all))
  ),"tmp/reference-runtime/drift.json",pretty=TRUE,auto_unbox=TRUE,null="null")
  if(length(missing_all)||length(extra_all)||length(changed_all) ||
     !identical(as.integer(info$packageCount),as.integer(length(all_versions))) ||
     !identical(inventory_hash,info$sha256SortedPackageVersionLines)) {
    changed_text<-if(length(changed_all))
      paste(paste0(changed_all," ",expected_all[changed_all],"->",all_versions[changed_all]),collapse=", ")
      else "none"
    stop("Frozen transitive R library drift: expected SHA-256 ",
         info$sha256SortedPackageVersionLines,
         ", observed ",inventory_hash,
         ", packages ",length(all_versions),
         "; changed=",changed_text,
         "; missing=",paste(missing_all,collapse=","),
         "; extra=",paste(extra_all,collapse=","))
  }
}
cat("MULTIOMICS_R_LIBRARY_SHA256=",inventory_hash,
    " COUNT=",length(all_versions),"\n",sep="")
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
  full_installed_library_sha256=inventory_hash,
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
