#!/usr/bin/env Rscript
# Normalize the small bootstrap layer that basilisk upgrades without exact
# versions before running the already-pinned MOFA2/mofapy2 scientific stack.
for(p in c("MOFA2","basilisk","reticulate","jsonlite"))
  if(!requireNamespace(p,quietly=TRUE)) stop("Missing authentic method dependency: ",p)

config <- jsonlite::read_json(
  "multiomics-engine/environment/mofa-python-pins.json",
  simplifyVector=FALSE
)
if(!identical(config$schema,"multiomics.mofa-basilisk-python.v1"))
  stop("Unknown Basilisk environment lock schema")

bootstrap <- c("pip","wheel","setuptools","packaging")
expected <- unlist(config$pinnedPackages[bootstrap],use.names=TRUE)
if(length(expected)!=length(bootstrap) || anyNA(expected))
  stop("MOFA Python lock must declare exact bootstrap package versions")

env <- tryCatch(
  get("mofa_env",asNamespace("MOFA2"),inherits=FALSE),
  error=function(e) stop("Internal MOFA2 Basilisk environment unavailable: ",conditionMessage(e))
)
proc <- basilisk::basiliskStart(env)
on.exit(basilisk::basiliskStop(proc),add=TRUE)

result <- basilisk::basiliskRun(proc,function(expected){
  expected_json <- jsonlite::toJSON(as.list(expected),auto_unbox=TRUE)
  Sys.setenv(MULTIOMICS_EXPECTED_BOOTSTRAP_JSON=expected_json)
  reticulate::py_run_string("
import json, os, re, subprocess, sys
from importlib import metadata
pins = json.loads(os.environ['MULTIOMICS_EXPECTED_BOOTSTRAP_JSON'])
req = [f\"{name}=={version}\" for name, version in pins.items()]
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', '--no-user', *req])
actual = {}
for package in metadata.distributions():
    name = package.metadata.get('Name')
    if name:
        actual[re.sub(r'[-_.]+','-',name).lower()] = package.version
multiomics_bootstrap_json = json.dumps({k: actual.get(k) for k in pins}, sort_keys=True)
")
  as.character(reticulate::py$multiomics_bootstrap_json)
}, expected=as.list(expected))

actual <- jsonlite::fromJSON(result,simplifyVector=TRUE)
for(name in names(expected)){
  if(!identical(unname(actual[[name]]),unname(expected[[name]])))
    stop("Failed to normalize ",name,": expected ",expected[[name]],
         ", got ",actual[[name]])
}
cat("PASS: normalized MOFA2 Basilisk bootstrap packages: ",
    paste(paste0(names(expected),"=",unname(expected)),collapse=", "), "\n",sep="")
