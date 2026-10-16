#!/usr/bin/env Rscript
# R Markdown bootstrap / Préparer le carnet R Markdown.
# Run from the fully extracted MultiOmics research ZIP:
#   Rscript render_notebook.R
# Only missing REPORTING packages are installed. Original scientific package
# versions and R sources are not modified or silently upgraded.
if (!file.exists("ANALYSE_REPRODUCTIBLE.Rmd"))
  stop("Run from the extracted ZIP root (missing ANALYSE_REPRODUCTIBLE.Rmd).")
if (!requireNamespace("rmarkdown", quietly = TRUE)) {
  message("Installing the missing R Markdown report package into your R library...")
  install.packages("rmarkdown", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("rmarkdown", quietly = TRUE))
  stop("The rmarkdown package could not be installed; check internet and R library permissions.")
if (!rmarkdown::pandoc_available())
  stop("Pandoc is unavailable. Install RStudio or Pandoc to generate the HTML report.")
message("Replaying native R when available, generating R figures and knitting the methods notebook.")
rmarkdown::render("ANALYSE_REPRODUCTIBLE.Rmd",
                  output_file = "analyse_reproductible.html",
                  params = list(execute_replay = TRUE),
                  envir = new.env(parent = globalenv()),
                  quiet = FALSE)
stopifnot(file.exists("analyse_reproductible.html"))
message("Report generated locally: analyse_reproductible.html")
