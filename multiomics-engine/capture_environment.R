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
    fields <- colnames(installed)
    # The installed.packages() matrix does not guarantee a Repository field.
    # Never assign a fictitious source for locally compiled/installed packages.
    get_field <- function(i,name) {
      if(name %in% fields) as.character(installed[i,name]) else NA_character_
    }
    lapply(seq_len(nrow(installed)), function(i) list(
      package = get_field(i,"Package"),
      version = get_field(i,"Version"),
      repository = get_field(i,"Repository"),
      built = get_field(i,"Built")
    ))
  }),
  installed_inventory_available_columns = colnames(utils::installed.packages()),
  installed_repository_metadata_available =
    "Repository" %in% colnames(utils::installed.packages()),
  renv_lockfile_present = file.exists("renv.lock") ||
    file.exists("multiomics-engine/renv.lock"),
  note = "Inventory only. This is not a reproducible renv.lock and does not guarantee reinstall."
)
writeLines(jsonlite::toJSON(info, auto_unbox = TRUE, pretty = TRUE, null = "null"), output)
cat("Wrote actual R environment inventory:", output, "\n")
