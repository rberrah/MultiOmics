#!/usr/bin/env Rscript
# Ordinary installation remains convenient but does not freeze package versions.
# For verifiable reference work use --verify-frozen, which NEVER updates packages.
if ("--verify-frozen" %in% commandArgs(trailingOnly=TRUE)) {
  rs <- file.path(R.home("bin"),"Rscript")
  status <- system2(rs,
    c("multiomics-engine/environment/check_reference_runtime.R","verify"))
  if(!identical(as.integer(status),0L))
    stop("Frozen R reference environment mismatch. Refusing to install/upgrade silently.")
  cat("PASS: installed scientific reference environment matches pinned versions.\n")
  quit(status=0L)
}
cran <- c("plumber","jsonlite","BiocManager","lmerTest")
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos="https://cloud.r-project.org")

bioc <- c(
  "DESeq2","edgeR","limma","fgsea","MOFA2","mixOmics",
  "xcms","MsExperiment","Spectra","mzR","BiocParallel"
)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly=TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask=FALSE, update=FALSE)

cat("Backend package availability:\n")
for (pkg in c("plumber","jsonlite","DESeq2","edgeR","limma","lmerTest","fgsea","MOFA2","mixOmics","xcms","MsExperiment","Spectra","mzR","BiocParallel")) {
  cat(sprintf("  %-14s %s\n", pkg, if (requireNamespace(pkg, quietly=TRUE)) "OK" else "MISSING"))
}
