#!/usr/bin/env Rscript
# Research ZIP entry point: exactly replays native R results where available,
# creates scientifically similar figure equivalents from the original
# plot-ready values, and reports these DIFFERENT levels of fidelity.
# From extracted ZIP root: Rscript run_all.R
args <- commandArgs(trailingOnly=TRUE)
if(length(args))stop("No arguments: run this script from the extracted ZIP root.")
if(!file.exists("provenance/analysis_complete.json"))
  stop("Run from the extracted ZIP root: missing provenance/analysis_complete.json")
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("jsonlite is required; restore the R environment before starting")
status <- list(
  schema_version="1.0.0",
  r_numerical_replay="NOT_AVAILABLE_BROWSER_RESULTS",
  r_figures="NOT_CHECKED",
  r_figures_equivalence="SCIENTIFICALLY_SIMILAR_NOT_PIXEL_IDENTICAL",
  messages=character()
)
rscript <- file.path(R.home("bin"),if(.Platform$OS.type=="windows")"Rscript.exe" else "Rscript")
call_r <- function(script,arguments) {
  result <- suppressWarnings(system2(rscript,c(shQuote(script),vapply(arguments,shQuote,character(1))),
                                     stdout=TRUE,stderr=TRUE))
  code <- attr(result,"status")
  if(is.null(code))code <- 0L
  if(code != 0L)stop("Failed: ",script,"\n",paste(result,collapse="\n"))
  invisible(result)
}
# A replay is possible only if the original request/response and the exact
# R engine source snapshot were stored when the browser analysis was run.
if(file.exists("reproduction/run_reference_backend.R")) {
  status$r_numerical_replay <- "NOT_CHECKED"
  tryCatch({
    call_r("reproduction/run_reference_backend.R",c("reproduction"))
    call_r("reproduction/verify_replay.R",c("reproduction"))
    status$r_numerical_replay <- "PASS_R_REFERENCE_ONLY"
  },error=function(e){
    status$r_numerical_replay <- "FAIL_OR_ENVIRONMENT_MISSING"
    status$messages <- c(status$messages,conditionMessage(e))
  })
}
tryCatch({
  call_r("figures_r/generate_figures.R",
         c("provenance/analysis_complete.json","figures_r/output"))
  status$r_figures <- "GENERATED_SIMILAR_PDF_SVG"
},error=function(e){
  status$r_figures <- "FAIL"
  status$messages <- c(status$messages,conditionMessage(e))
})
status$browser_numeric_r_replay <- "NOT_SUPPORTED_JS_IS_CANONICAL_FOR_BROWSER"
status$original_browser_svg <- "SNAPSHOT_ONLY"
writeLines(jsonlite::toJSON(status,pretty=TRUE,auto_unbox=TRUE),
           "research_reproduction_status.json")
cat("R numeric replay:",status$r_numerical_replay,"\n")
cat("R figure generation:",status$r_figures,"\n")
cat("Details: research_reproduction_status.json\n")
if(identical(status$r_figures,"FAIL") || identical(status$r_numerical_replay,"FAIL_OR_ENVIRONMENT_MISSING"))
  quit(status=1L)
