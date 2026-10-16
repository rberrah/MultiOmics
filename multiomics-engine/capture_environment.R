#!/usr/bin/env Rscript
# Record the actual local R/Bioconductor environment. Never invent versions.
args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "multiomics-r-environment.json"
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("jsonlite is required to write this environment report")
}
packages <- c(
  "DESeq2", "edgeR", "limma", "fgsea", "MOFA2", "mixOmics",
  "BiocManager", "BiocParallel", "plumber", "jsonlite",
  "lmerTest", "xcms", "MsExperiment", "Spectra", "mzR"
)
pkg_versions <- setNames(lapply(packages, function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE))
    as.character(utils::packageVersion(pkg))
  else NULL
}), packages)
bioc_version <- if (requireNamespace("BiocManager", quietly = TRUE))
  as.character(BiocManager::version()) else NULL
info <- list(
  schema_version = "1.0.0",
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  r_version = as.character(getRversion()),
  platform = R.version$platform,
  bioconductor_version = bioc_version,
  package_versions = pkg_versions,
  installed_package_table = local({
    installed <- utils::installed.packages()
    lapply(seq_len(nrow(installed)), function(i)
      list(package = installed[i, "Package"], version = installed[i, "Version"],
           repository = installed[i, "Repository"], built = installed[i, "Built"]))
  }),
  renv_lockfile_present = file.exists("renv.lock") ||
    file.exists("multiomics-engine/renv.lock"),
  note = "Inventory only. This is not a reproducible renv.lock and does not guarantee reinstall."
)
writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE, null = "null"), output)
cat("Wrote actual R environment inventory:", output, "\n")
