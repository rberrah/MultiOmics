#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
`%||%` <- function(x, y) if (is.null(x) || !length(x) || all(is.na(x))) y else x
host <- if (length(args) >= 1L) args[[1]] else Sys.getenv("PMX_MULTIOMICS_HOST", "127.0.0.1")
port <- if (length(args) >= 2L) as.integer(args[[2]]) else as.integer(Sys.getenv("PMX_MULTIOMICS_PORT", "8787"))

if (!requireNamespace("plumber", quietly=TRUE)) stop("Install plumber first: install.packages('plumber')")
if (!requireNamespace("jsonlite", quietly=TRUE)) stop("Install jsonlite first: install.packages('jsonlite')")

loopback_hosts <- c("127.0.0.1", "localhost", "::1")
if (!host %in% loopback_hosts && tolower(Sys.getenv("PMX_MULTIOMICS_ALLOW_REMOTE", "no")) != "yes") {
  stop(
    paste(
      "Refusing to expose the multi-omics R backend on a non-loopback host.",
      "For a deliberate remote deployment set PMX_MULTIOMICS_ALLOW_REMOTE=yes",
      "and place the service behind TLS, authentication, request-size limits and an origin allow-list."
    ),
    call. = FALSE
  )
}

split_origins <- function(value) {
  value <- trimws(unlist(strsplit(value, ",", fixed=TRUE), use.names=FALSE))
  unique(value[nzchar(value)])
}

default_origins <- c(
  "https://rberrah.github.io",
  "http://localhost:5173",
  "http://127.0.0.1:5173",
  "http://localhost:4173",
  "http://127.0.0.1:4173"
)
configured_origins <- split_origins(Sys.getenv("PMX_MULTIOMICS_ALLOWED_ORIGINS", ""))
allowed_origins <- unique(c(default_origins, configured_origins))

reference_packages <- c(
  "plumber","jsonlite","DESeq2","edgeR","limma","lmerTest","fgsea","MOFA2","mixOmics",
  "xcms","MsExperiment","Spectra","mzR","BiocParallel"
)
package_manifest <- function() {
  versions <- lapply(reference_packages, function(pkg) {
    if (!requireNamespace(pkg, quietly=TRUE)) return(NULL)
    as.character(utils::packageVersion(pkg))
  })
  names(versions) <- reference_packages
  list(
    r = R.version.string,
    platform = R.version$platform,
    packages = versions,
    library_paths = as.list(.libPaths())
  )
}

api <- plumber::plumb(file.path("multiomics-engine", "server.R"))

# The annotated server historically emitted Access-Control-Allow-Origin: *.
# This second filter is intentionally placed before endpoint execution so an
# unrelated web page cannot use a browser to submit research data to the local
# service. For accepted origins it also overwrites the permissive header with
# the exact requesting origin.
api <- plumber::pr_filter(
  api,
  "origin_guard",
  function(req, res) {
    origin <- req$HTTP_ORIGIN
    res$setHeader("Cache-Control", "no-store")
    if (is.null(origin) || !nzchar(origin)) return(plumber::forward())
    if (!origin %in% allowed_origins) {
      res$setHeader("Access-Control-Allow-Origin", "null")
      res$status <- 403
      return(list(
        status = "error",
        message = "Origin is not allowed by the PMx multi-omics backend. Configure PMX_MULTIOMICS_ALLOWED_ORIGINS explicitly if needed."
      ))
    }
    res$setHeader("Access-Control-Allow-Origin", origin)
    res$setHeader("Vary", "Origin")
    res$setHeader("Access-Control-Allow-Methods", "GET,POST,OPTIONS")
    res$setHeader("Access-Control-Allow-Headers", "Content-Type")
    if (identical(tolower(req$HTTP_ACCESS_CONTROL_REQUEST_PRIVATE_NETWORK %||% ""), "true")) {
      res$setHeader("Access-Control-Allow-Private-Network", "true")
    }
    # Browsers cannot POST JSON until their OPTIONS preflight succeeds with
    # the requesting origin's exact value. Never return "*" for local data.
    if (identical(req$REQUEST_METHOD, "OPTIONS")) {
      res$status <- 204
      return(list(status="ok"))
    }
    plumber::forward()
  },
  serializer = plumber::serializer_unboxed_json()
)

# Explicit environment endpoint for audit trails. It contains versions only;
# no study data, file paths from the submitted payload, or secrets are exposed.
api <- plumber::pr_get(
  api,
  "/environment",
  function() package_manifest(),
  serializer = plumber::serializer_unboxed_json()
)

# Attach the exact R/package environment directly to every /run response.
# This ensures the browser JSON/HTML export remains self-describing even if the
# standalone /environment endpoint was never called. Plumber postroute hooks
# may replace the handler value; all non-/run routes are returned unchanged.
api <- plumber::pr_hook(
  api,
  "postroute",
  function(req, value) {
    if (req$PATH_INFO %in% c("/run", "/external-validation") && is.list(value)) {
      value$environment <- package_manifest()
    }
    value
  }
)

manifest <- package_manifest()
cat(sprintf("PMx multi-omics reference R backend: http://%s:%d\n", host, port))
cat(sprintf("R runtime: %s\n", manifest$r))
cat(sprintf("Allowed browser origins: %s\n", paste(allowed_origins, collapse=", ")))
for (pkg in names(manifest$packages)) {
  if (!is.null(manifest$packages[[pkg]])) cat(sprintf("  %-14s %s\n", pkg, manifest$packages[[pkg]]))
}
api$run(host=host, port=port, swagger=TRUE)
