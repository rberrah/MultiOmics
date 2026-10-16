/**
 * PMx Explain portable project, schema v1.
 * Captures the EXACT original input bytes (as UTF-8 text) and their SHA-256
 * hashes plus parameters and an optional analysis snapshot in ONE file.
 *
 * NOTE: CSV/TSV are text-only in this contract. Raw mzML/mzXML/binary MS data
 * must be kept outside this project bundle and are never silently included.
 * Imported snapshots are never executed or presented as newly validated results:
 * the user explicitly reruns the deterministic workflow.
 *
 * No fetch, browser storage, JavaScript evaluation or connected-app access.
 */
export const PORTABLE_PROJECT_FORMAT = 'pmx-multiomics-portable-project';
export const PORTABLE_PROJECT_SCHEMA = 1;
export const PORTABLE_PROJECT_MAX_BYTES = 80 * 1024 * 1024;
const ALL_INPUT_KEYS = ['metadata','transcriptomics','proteomics','metabolomics','msOriginal','msUserAnnotations'];
const ALLOWED_SETTINGS = [
  'studyName','objective','organism','unitType','studySetting','designType',
  'groupCount','timepointCount','sampleOverlap','technicalReplicatesExpected',
  'batchKnown','outcomeType','outcomeTimepoint','covariatesAvailable',
  'partialOmicsExpected','selectedCovariates','groupVariable','outcome',
  'subjectCount','paired','longitudinal','transcriptomicsPlatform',
  'transcriptomicsValues','transcriptomicsIdType','proteomicsPlatform',
  'proteomicsValues','proteomicsIdType','metabolomicsPlatform',
  'metabolomicsValues','metabolomicsIdType','msBlankFilter','msBlankFold',
  'msQcRsdFilter','msQcRsdThreshold','msDriftCorrection','msMnarStrategy',
  'analysisIntent','validateDiabloHoldout','confirmIndependentAssays','msAutoMetadata','demoLoaded',
  'columnMapping','metaboliteAnnotations','msAnnotationCsv',
  'msAnnotationFields','msImportFormat','annotationFileName'
];
const textEncoder = new TextEncoder();

/** @param {string} text */
export async function sha256Content(text) {
  const cryptoApi = globalThis.crypto;
  if (!cryptoApi?.subtle) throw new Error('Secure SHA-256 is unavailable. Use HTTPS or a secure local browser.');
  const hash = await cryptoApi.subtle.digest('SHA-256',textEncoder.encode(text));
  return Array.from(new Uint8Array(hash)).map(byte=>byte.toString(16).padStart(2,'0')).join('');
}
/** @param {Record<string,any>} settings */
function onlyPublicSettings(settings) {
  /** @type {Record<string,any>} */
  const safe = {};
  for (const name of ALLOWED_SETTINGS) {
    if (settings && Object.prototype.hasOwnProperty.call(settings,name))
      safe[name] = settings[name];
  }
  return JSON.parse(JSON.stringify(safe));
}
/** @param {string} name */
function safeName(name) {
  const cleaned=String(name||'input.csv').replace(/[\x00-\x1F\x7F]/g,'')
    .split(/[\\/]/).pop()?.slice(0,180) || '';
  return cleaned || 'input.csv';
}
/** @param {number} size */
function assertAllowedSize(size) {
  if(!Number.isFinite(size)||size>PORTABLE_PROJECT_MAX_BYTES)
    throw new Error('Project exceeds the 80 MiB portable text limit. Export the original large files separately and use the normal JSON results export.');
}
/** @param {{files:Record<string,File|null>, settings:Record<string,any>, result?:any, extraFiles?:Record<string,File|null>}} args */
export async function createPortableProject({files,settings,result=null,extraFiles={}}) {
  if(!files?.metadata)throw new Error('A sample sheet is required to export a reproducible project.');
  const full={...files,...extraFiles};
  /** @type {Record<string, {filename:string,mime:string,sizeBytes:number,sha256:string,content:string}>} */
  const out={};
  const size=ALL_INPUT_KEYS.reduce((n,key)=>n+(full[key]?.size||0),0);
  assertAllowedSize(size);
  for(const key of ALL_INPUT_KEYS) {
    const file=full[key];
    if(!file)continue;
    const content=await file.text();
    out[key]={
      filename:safeName(file.name),mime:String(file.type||'text/plain'),
      sizeBytes:textEncoder.encode(content).byteLength,
      sha256:await sha256Content(content),content
    };
  }
  return {
    format:PORTABLE_PROJECT_FORMAT,
    schemaVersion:PORTABLE_PROJECT_SCHEMA,
    exportedAt:new Date().toISOString(),
    privacy:'SENSITIVE: contains original study matrices, subjects and analysis results. Store and share only within an approved scientific data-management environment.',
    limitations:'Only uploaded text matrices/metadata and optional original text MS exports are portable here. No backend R binaries, mzML files, software package lockfiles, or remote database snapshots. Re-run methods on verified compatible installations.',
    settings:onlyPublicSettings(settings),
    inputs:out,
    archivedResult:result,
    externalSettingsExcluded:['backend URL','backend credentials','automatic API connections'],
    rerunPolicy:'Import validates SHA-256 and restores only input files and protocol settings; it never runs the project or sends data to R / APIs without user action.'
  };
}
/**
 * @param {string | any} payload
 * @returns {Promise<{files:Record<string,File>,settings:Record<string,any>,archivedResultPresent:boolean,archiveDate:string,originalDigests:Record<string,string>}>}
 */
export async function parsePortableProject(payload) {
  if(typeof payload==='string')assertAllowedSize(textEncoder.encode(payload).byteLength);
  const value=typeof payload==='string'?JSON.parse(payload):payload;
  if(!value||value.format!==PORTABLE_PROJECT_FORMAT||
     value.schemaVersion!==PORTABLE_PROJECT_SCHEMA)
    throw new Error('Unsupported portable project: expected PMx multiomics schema v1.');
  if(!value.inputs||typeof value.inputs!=='object'||!value.inputs.metadata)
    throw new Error('Portable project is missing its original sample sheet.');
  const inputKeys=Object.keys(value.inputs);
  if(inputKeys.some(key=>!ALL_INPUT_KEYS.includes(key)))
    throw new Error('Portable project contains an unknown input category.');
  const total=inputKeys.reduce((n,key)=>n+Number(value.inputs[key]?.sizeBytes||0),0);
  assertAllowedSize(total);
  /** @type {Record<string,File>} */
  const restored={};
  for(const key of inputKeys) {
    const input=value.inputs[key];
    if(!input||typeof input.content!=='string'||typeof input.sha256!=='string'||
       !/^[a-f0-9]{64}$/.test(input.sha256)||typeof input.filename!=='string')
      throw new Error('Malformed portable input: '+key);
    const actualBytes=textEncoder.encode(input.content).byteLength;
    if(actualBytes!==input.sizeBytes ||
       await sha256Content(input.content)!==input.sha256)
      throw new Error('Original input integrity check failed (SHA-256): '+key);
    restored[key]=new File([input.content],safeName(input.filename),{
      type:typeof input.mime==='string'?input.mime:'text/plain'
    });
  }
  return {
    files:restored,
    settings:onlyPublicSettings(value.settings||{}),
    archivedResultPresent:!!value.archivedResult,
    archiveDate:typeof value.exportedAt==='string'?value.exportedAt:'',
    originalDigests:Object.fromEntries(inputKeys.map(key=>[key,value.inputs[key].sha256]))
  };
}
