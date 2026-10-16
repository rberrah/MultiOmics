#!/usr/bin/env Rscript
# Capture the REAL isolated Python dependency closure used by MOFA2::run_mofa
# with use_basilisk=TRUE. Does not inspect ambient system Python.
args <- commandArgs(trailingOnly=TRUE)
mode <- if(length(args)) args[[1L]] else "probe"
if(!mode %in% c("probe","verify")) stop("Mode must be probe or verify")
for(p in c("MOFA2","basilisk","reticulate","jsonlite","digest"))
  if(!requireNamespace(p,quietly=TRUE))stop("Missing authentic method dependency: ",p)
config <- jsonlite::read_json("multiomics-engine/environment/mofa-python-pins.json",
                              simplifyVector=FALSE)
if(!identical(config$schema,"multiomics.mofa-basilisk-python.v1"))
  stop("Unknown Basilisk environment lock schema")
if(!identical(as.character(utils::packageVersion("MOFA2")),config$mofaRVersion))
  stop("Python lock belongs to another MOFA2 R release")
target <- tryCatch(as.character(get(".mofapy2_version",asNamespace("MOFA2"),inherits=FALSE)),
                   error=function(e)NA_character_)
if(!identical(target,config$mofapy2Target))
  stop("The actual MOFA2 R namespace targets a different mofapy2 release")
env <- tryCatch(get("mofa_env",asNamespace("MOFA2"),inherits=FALSE),
                error=function(e)stop("Internal MOFA2 Basilisk environment unavailable: ",conditionMessage(e)))
proc <- basilisk::basiliskStart(env)
on.exit(basilisk::basiliskStop(proc),add=TRUE)
actualJSON <- basilisk::basiliskRun(proc,function(){
  reticulate::py_run_string("
import sys, json, re
from importlib import metadata
records = []
for package in metadata.distributions():
    name = package.metadata.get('Name')
    if not name:
        raise RuntimeError('Basilisk distribution without a package name')
    canonical = re.sub(r'[-_.]+', '-', name).lower()
    records.append({'name':canonical, 'version':package.version})
records.sort(key=lambda r:r['name'])
if len({r['name'] for r in records}) != len(records):
    raise RuntimeError('Duplicate normalized Basilisk package names')
multiomics_basilisk_json = json.dumps({'python':sys.version.split()[0],
                                      'packages':records},sort_keys=True)
")
  as.character(reticulate::py$multiomics_basilisk_json)
})
actual <- jsonlite::fromJSON(actualJSON,simplifyVector=FALSE)
versions <- vapply(actual$packages,function(item)item$version,character(1))
names(versions) <- vapply(actual$packages,function(item)item$name,character(1))
if(!identical(unname(versions["mofapy2"]),config$mofapy2Target))
  stop("Basilisk installed an unexpected mofapy2 version: ",versions["mofapy2"])
if(!startsWith(actual$python,paste0(config$pythonMajorMinor,".")))
  stop("Basilisk Python interpreter version drift: ",actual$python)
canonical <- paste0(paste(paste0(names(versions),"=",unname(versions)),
                          collapse="\n"),"\n")
fingerprint <- digest::digest(canonical,algo="sha256",serialize=FALSE)
if(identical(mode,"verify")){
  if(!isTRUE(config$strictPythonLockReady) ||
     !identical(as.integer(length(versions)),
                as.integer(config$packageCount)) ||
     !identical(fingerprint,config$fullPythonLibrarySha256))
    stop("Basilisk transitive Python environment drift or pending lock: ",fingerprint)
}
dest <- "tmp/multiomics/mofa-nutrimouse-real"
dir.create(dest,recursive=TRUE,showWarnings=FALSE)
record <- list(
  schema="multiomics.mofa-basilisk-provenance.v1",
  source="actual MOFA2 basiliskStart(mofa_env) Python environment, not ambient Python",
  runMode=mode,
  rPackageVersion=as.character(utils::packageVersion("MOFA2")),
  pythonVersion=actual$python,
  mofapy2Version=unname(versions["mofapy2"]),
  installedPythonPackages=actual$packages,
  installedPackageCount=length(versions),
  installedPythonLibrarySha256=fingerprint,
  strictLockVerified=identical(mode,"verify")&&isTRUE(config$strictPythonLockReady),
  sourceCommit=Sys.getenv("GITHUB_SHA","local"),
  limitation="Dependency fingerprint proves observed package versions, not offline reinstallability or cross-hardware bitwise numerical identity."
)
jsonlite::write_json(record,file.path(dest,"basilisk_python_environment.json"),
                     pretty=TRUE,auto_unbox=TRUE)
cat("MOFA2_BASILISK_PYTHON_SHA256=",fingerprint,
    " PACKAGES=",length(versions)," PYTHON=",actual$python,
    " MOFAPY2=",versions["mofapy2"],"\n",sep="")
cat("PASS: captured isolated MOFA2 Python runtime, mode ",mode,"\n",sep="")
