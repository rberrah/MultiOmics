#!/usr/bin/env Rscript
# Normalize the small dependency layer that Basilisk or unconstrained
# transitive resolution can upgrade before running the already-pinned
# MOFA2/mofapy2 scientific stack.
for(p in c("MOFA2","basilisk","reticulate","jsonlite"))
  if(!requireNamespace(p,quietly=TRUE)) stop("Missing authentic method dependency: ",p)

config <- jsonlite::read_json(
  "multiomics-engine/environment/mofa-python-pins.json",
  simplifyVector=FALSE
)
if(!identical(config$schema,"multiomics.mofa-basilisk-python.v1"))
  stop("Unknown Basilisk environment lock schema")

bootstrap <- c("pip","wheel","setuptools","packaging","python-dateutil")
expected <- unlist(config$pinnedPackages[bootstrap],use.names=TRUE)
if(length(expected)!=length(bootstrap) || anyNA(expected))
  stop("MOFA Python lock must declare exact normalized package versions")

env <- tryCatch(
  get("mofa_env",asNamespace("MOFA2"),inherits=FALSE),
  error=function(e) stop("Internal MOFA2 Basilisk environment unavailable: ",conditionMessage(e))
)
proc <- basilisk::basiliskStart(env)
on.exit(basilisk::basiliskStop(proc),add=TRUE)

result <- basilisk::basiliskRun(proc,function(expected){
  if(!requireNamespace("reticulate",quietly=TRUE) ||
     !requireNamespace("jsonlite",quietly=TRUE))
    stop("Missing reticulate/jsonlite inside Basilisk worker")
  py_run <- getExportedValue("reticulate","py_run_string")
  py_eval <- getExportedValue("reticulate","py_eval")
  to_json <- getExportedValue("jsonlite","toJSON")
  expected_json <- as.character(to_json(as.list(expected),auto_unbox=TRUE))
  expected_literal <- as.character(to_json(expected_json,auto_unbox=TRUE))
  code <- paste0(
    "import json, re, subprocess, sys\n",
    "from importlib import metadata\n",
    "pins = json.loads(", expected_literal, ")\n",
    "req = [f'{name}=={version}' for name, version in pins.items()]\n",
    "subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', '--no-user', *req])\n",
    "actual = {}\n",
    "for package in metadata.distributions():\n",
    "    name = package.metadata.get('Name')\n",
    "    if name:\n",
    "        actual[re.sub(r'[-_.]+','-',name).lower()] = package.version\n",
    "multiomics_bootstrap_json = json.dumps({k: actual.get(k) for k in pins}, sort_keys=True)\n"
  )
  py_run(code)
  as.character(py_eval("multiomics_bootstrap_json",convert=TRUE))
}, expected=as.list(expected))

actual <- jsonlite::fromJSON(result,simplifyVector=TRUE)
for(name in names(expected)){
  if(!identical(unname(actual[[name]]),unname(expected[[name]])))
    stop("Failed to normalize ",name,": expected ",expected[[name]],
         ", got ",actual[[name]])
}
cat("PASS: normalized MOFA2 Basilisk drift-sensitive packages: ",
    paste(paste0(names(expected),"=",unname(expected)),collapse=", "), "\n",sep="")
