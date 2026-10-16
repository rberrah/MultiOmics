// @ts-nocheck
// A faithful R replay is possible ONLY for a successful reference R backend.
// Browser differential statistics and Svelte/SVG figures are NOT generated in R.
import {REFERENCE_R_SOURCE_FILES} from './reference-r-sources.js';
import {REFERENCE_R_REPLAY_VERIFY_PY} from './reference-replay-verifier-source.js';
import {RESEARCH_R_NUMERIC_VERIFIER_SOURCE} from './r-numeric-verifier-source.js';

const REPLAY_R = String.raw`#!/usr/bin/env Rscript
# Replay the exact JSON request with the same source functions that serviced /run.
# Unzip the research dossier, then execute:
# Rscript reproduction/run_reference_backend.R reproduction/
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1L) stop("Usage: Rscript run_reference_backend.R REPRODUCTION_DIRECTORY")
root <- normalizePath(args[[1]],mustWork=TRUE)
setwd(root)
if(!requireNamespace("jsonlite",quietly=TRUE))
  stop("Install jsonlite in the exact reference backend R environment.")
if(!file.exists("reference_request.json"))stop("Missing original R backend request")
if(!requireNamespace("digest",quietly=TRUE))stop("digest is needed to attest source SHA-256")
expected <- jsonlite::fromJSON("reference_expected.json",simplifyVector=FALSE)
original_runtime <- expected$runtime
if(is.null(original_runtime)||!identical(as.character(getRversion()),original_runtime$RVersion))
  stop("R version differs from the originally executed backend; exact replay blocked.")
if(!is.null(original_runtime$packageVersions)){
  for(pkg in names(original_runtime$packageVersions)){
    required_version <- as.character(original_runtime$packageVersions[[pkg]])
    if(!requireNamespace(pkg,quietly=TRUE) ||
       !identical(as.character(utils::packageVersion(pkg)),required_version))
      stop("R package missing or at another version: ",pkg," (required ",required_version,")")
  }
}
expected_sha <- original_runtime$referenceSourceSha256
if(is.null(expected_sha)||!length(expected_sha))
  stop("Original backend did not attest its source hashes. Exact R replay is NOT verifiable. Update and rerun the backend.")
for(name in names(expected_sha)){
  source_path <- file.path("multiomics-engine",name)
  if(!file.exists(source_path))stop("Missing original backend source: ",name)
  actual_hash <- digest::digest(file=source_path,algo="sha256")
  if(!identical(actual_hash,expected_sha[[name]]))
    stop("R source SHA-256 differs from original backend: ",name)
}
request <- jsonlite::fromJSON("reference_request.json",simplifyVector=FALSE)
source(file.path("multiomics-engine","server.R"),local=FALSE)
actual <- run_backend_analysis(request)
json <- jsonlite::toJSON(actual,auto_unbox=TRUE,dataframe="rows",
                         na="null",null="null",digits=NA)
writeLines(json,"reference_replayed.json",useBytes=TRUE)
capture.output(sessionInfo(),file="R_session_info_replayed.txt")
cat("R reference backend calculations replayed. Compare reference_replayed.json",
    "with reference_expected.json and retained R CSV tables.\n")
cat("Browser-generated statistics and Svelte SVG figures are NOT R-replayed.\n")
`;

const REPLAY_GUIDE = `# Reproduction du calcul R / Native R replay

1. Extraire l'ensemble du dossier ZIP sur votre ordinateur. **Les fichiers contiennent potentiellement des données confidentielles.**
2. Restaurer la version de R et les packages indiqués dans le résultat original. Le script exige également que les SHA-256 du code R correspondent aux empreintes transmises par le backend d'origine. Si elles manquent ou divergent, le rejeu échoue sans attribuer de fausse équivalence.
3. Depuis la racine du ZIP : \`Rscript reproduction/run_reference_backend.R reproduction/\`.
4. Exécuter \`Rscript reproduction/verify_replay.R reproduction/\` pour comparer les méthodes R, toutes les valeurs et le prétraitement sans Python. Le vérificateur Python indépendant reste disponible en option. Un écart interdit toute revendication de reproductibilité.
5. Vérifier les tableaux R complets dans \`reference_r/\`. Les résultats trop volumineux conservés uniquement sur le disque R local doivent être récupérés séparément.

**Périmètre :** les fonctions R incluses ci-dessous sont une copie des sources exactes de l'API R lors de la préparation du site. Le fichier \`reference_request.json\` est l'entrée exacte de l'appel R effectué pendant l'analyse. Il faut aussi les mêmes dépendances, graînes aléatoires et versions des bases externes. Les accès en direct à Reactome ne sont pas figés.

**Limite majeure :** les tableaux exploratoires du navigateur et les figures SVG de Svelte ne sont PAS générés en R. Le script R ne peut pas prétendre en recréer les mêmes valeurs et graphismes. Les SVG exportés sont des captures vectorielles de l'affichage initial et non de nouvelles figures générées par R. L'équivalence numérique et graphique totale exige de migrer ces calculs et rendus vers une seule implémentation canonique, puis une validation sur plusieurs jeux publics.

**Status :** source-exact replay candidate; statistical/numerical equivalence NOT certified until independently checked on the exact deployment and package lock.
`;

export function referenceReplayStatus(replay){
  if(!replay?.requestJson || !replay?.responseJson) return {
    status:'not_available',scope:'browser_only',
    reason:'No successful R request/response was captured. Browser results are not R-generated.'
  };
  return {
    status:'r_reference_source_replay_available',
    scope:'r_backend_only',
    numerical_equivalence_verified:false,
    figure_equivalence_verified:false,
    limitations:['browser_statistics_use_javascript','interactive_svg_figures_use_svelte',
      'r_package_versions_and_external_mappings_must_be_restored']
  };
}
export async function referenceReplayEntries(replay){
  const status=referenceReplayStatus(replay);
  if(status.status==='not_available')return {entries:[],status};
  if(typeof replay.requestJson!=='string'||typeof replay.responseJson!=='string')
    throw Error('Invalid R replay request or response');
  JSON.parse(replay.requestJson);
  const response=JSON.parse(replay.responseJson);
  const expectedSources=response?.runtime?.referenceSourceSha256;
  let sourceAttestation='not_available';
  if(expectedSources && typeof expectedSources==='object'){
    const localHashes={};
    for(const [name,source] of Object.entries(REFERENCE_R_SOURCE_FILES)){
      const bytes=new TextEncoder().encode(source);
      const raw=await globalThis.crypto.subtle.digest('SHA-256',bytes);
      localHashes[name]=[...new Uint8Array(raw)].map(x=>x.toString(16).padStart(2,'0')).join('');
    }
    sourceAttestation=Object.keys(localHashes).length===Object.keys(expectedSources).length &&
      Object.entries(localHashes).every(([name,hash])=>expectedSources[name]===hash)
      ? 'sha256_match' : 'sha256_mismatch';
  }
  status.source_attestation=sourceAttestation;
  status.exact_runtime_source_verified=sourceAttestation==='sha256_match';
  const entries=[
    {name:'reproduction/run_reference_backend.R',data:REPLAY_R},
    {name:'reproduction/reference_request.json',data:replay.requestJson+'\n'},
    {name:'reproduction/reference_expected.json',data:replay.responseJson+'\n'},
    {name:'reproduction/README.md',data:REPLAY_GUIDE},
    {name:'reproduction/verify_replay.py',data:REFERENCE_R_REPLAY_VERIFY_PY},
    {name:'reproduction/verify_replay.R',data:RESEARCH_R_NUMERIC_VERIFIER_SOURCE}
  ];
  for(const [name,source] of Object.entries(REFERENCE_R_SOURCE_FILES))
    entries.push({name:'reproduction/multiomics-engine/'+name,data:source});
  return {entries,status};
}
