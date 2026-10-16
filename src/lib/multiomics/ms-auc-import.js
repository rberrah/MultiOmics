// @ts-nocheck
import { parseDelimited } from './deterministic.js';

// Exact sample and feature identifiers are preserved. No implicit averaging,
// sample-group guessing, MS peak annotation, concentration conversion or
// biological interpretation occurs in this importer.
const aliases = {
  assay: ['assay_id','injection_id','injection_name','file_name','filename','raw_file','raw_file_name','sample_name','sample_id','sample','injection','sampleid'],
  feature: ['feature_id','feature_name','metabolite_id','metabolite_name','compound_id','compound_name','analyte_name','metabolite','compound','feature','analyte','peak_id','peak_name','peak','name','identifier'],
  value: ['peak_area','peak_area_auc','peak_intensity','integrated_peak_area','auc','area','integrated_area','peakarea','area_under_curve','ms_auc','abundance','intensity'],
  subject: ['subject_id','subject','individual_id','patient_id','animal_id'],
  sample: ['sample_id','sample','specimen_id','biological_sample'],
  condition: ['condition','group','treatment','class','status'],
  batch: ['batch','batch_id','run_batch','plate'],
  timepoint: ['timepoint','time_point','time','visit'],
  type: ['sample_type','injection_type','sample_class','qc_type'],
  order: ['injection_order','run_order','injection_number','order'],
  replicate: ['technical_replicate','replicate','technical_rep']
};
const annotationKeys = new Set([
  'mz','m_z','m_over_z','mass_to_charge','retention_time','rt','rt_min',
  'retention_time_min','formula','molecular_formula','adduct','ion_mode',
  'polarity','annotation','confidence','score','compound_name','ion_formula',
  'charge','isotope','msms_score','library_match','identification_level',
  'hmdb_id','chebi_id','kegg_id','pubchem_id','inchikey'
]);
const norm = (x)=>String(x??'').trim().toLowerCase().replace(/[^a-z0-9]+/g,'_').replace(/^_+|_+$/g,'');
const find = (headers, keys) => headers.find((h)=>keys.includes(norm(h))) || null;
const escapeCsv=(s)=>'"'+String(s??'').replace(/"/g,'""')+'"';
const csv=(header,rows)=>[header.map(escapeCsv).join(','),...rows.map((row)=>row.map(escapeCsv).join(','))].join('\n')+'\n';
function parseArea(raw) {
  if(String(raw??'').trim()==='') return null;
  const v=Number(String(raw).trim().replace(',','.'));
  if(!Number.isFinite(v)||v<0) throw new Error('MS AUC: peak areas must be non-negative finite numbers; invalid value: '+raw);
  return v;
}
function metadataRowsFromRecords(records, cols) {
  if(!cols.condition) return null;
  const assays=new Map();
  for(const record of records){
    const assay=String(record[cols.assay]??'').trim();
    const v=(key, fallback='')=>cols[key] ? String(record[cols[key]]??'').trim() : fallback;
    const entry={
      subject_id:v('subject',assay),sample_id:v('sample',assay),assay_id:assay,omic:'metabolomics',
      condition:v('condition'),timepoint:v('timepoint'),batch:v('batch'),
      technical_replicate:v('replicate','1'),sample_type:v('type','biological'),injection_order:v('order')
    };
    if(!entry.subject_id||!entry.sample_id)throw new Error('MS AUC: missing subject or sample identifier');
    const old=assays.get(assay);
    if(old && JSON.stringify(old)!==JSON.stringify(entry))
      throw new Error('MS AUC: contradictory metadata for injection '+assay);
    assays.set(assay,entry);
  }
  if(!assays.size)return null;
  return [...assays.values()];
}
/**
 * Convert an MS vendor peak-area export to the existing feature × injection
 * contract. Long exports require a feature ID, injection/sample ID, and AUC.
 * Optional per-row condition/group columns can generate a traceable metadata
 * table. Unknown feature IDs are preserved and never guessed.
 */
export function convertMsAucExport(text) {
  const parsed=parseDelimited(text);
  const heads=parsed.headers;
  const cols=Object.fromEntries(Object.entries(aliases).map(([k,options])=>[k,find(heads,options)]));
  if (!cols.feature || !cols.assay || !cols.value || new Set([cols.feature,cols.assay,cols.value]).size!==3) {
    // Most MS exports are already feature × injection matrices. Drop only
    // standard peak-annotation columns (m/z, retention time, adduct...) so
    // they cannot accidentally masquerade as biological sample intensities.
    // Generic matrices without a recognized feature header are passed through.
    if (!cols.feature) return {format:'matrix',matrixCsv:text,metadataCsv:null,
      annotationsCsv:null,annotationColumns:[],sourceColumns:heads,features:null,assays:null};
    const annotationHeaders=heads.filter((header)=>header!==cols.feature&&annotationKeys.has(norm(header)));
    const sampleHeaders=heads.filter((header)=>header!==cols.feature&&!annotationKeys.has(norm(header)));
    if (!sampleHeaders.length)throw new Error('MS AUC: no injection intensity columns were found in the wide peak table');
    const seen=new Set();
    const records=parsed.rows.map((row)=>{
      const feature=String(row[cols.feature]||'').trim();
      if(!feature)throw new Error('MS AUC: a wide peak has no feature identifier');
      if(seen.has(feature))throw new Error('MS AUC: duplicate feature identifier in wide export: '+feature);
      seen.add(feature);
      return [feature,...sampleHeaders.map((header)=>parseArea(row[header])??'')];
    });
    const annotationsCsv=annotationHeaders.length
      ? csv(['feature_id',...annotationHeaders],
          parsed.rows.map((row)=>[String(row[cols.feature]||'').trim(),
            ...annotationHeaders.map((header)=>row[header]??'')]))
      : null;
    return {format:'wide_ms_auc',matrixCsv:csv(['feature_id',...sampleHeaders],records),
      metadataCsv:null,annotationsCsv,annotationColumns:annotationHeaders,sourceColumns:heads,
      features:records.length,assays:sampleHeaders.length,
      observed:records.reduce((n,row)=>n+row.slice(1).filter((v)=>v!=='').length,0),metadataGenerated:false};
  }
  if(!parsed.rows.length)throw new Error('MS AUC: export contains no observations');
  const features=[], assays=[], vals=new Map();
  const seenFeature=new Set(),seenAssay=new Set();
  let observed=0,missing=0;
  for(const row of parsed.rows){
    const feature=String(row[cols.feature]??'').trim();
    const assay=String(row[cols.assay]??'').trim();
    if(!feature||!assay)throw new Error('MS AUC: each row must contain a feature and a sample/injection ID');
    const value=parseArea(row[cols.value]);
    const key=JSON.stringify([feature,assay]);
    if(vals.has(key))throw new Error('MS AUC: duplicate feature × injection: '+feature+' / '+assay+'; resolve technical duplicates before import');
    vals.set(key,value);
    if(!seenFeature.has(feature)){features.push(feature);seenFeature.add(feature);}
    if(!seenAssay.has(assay)){assays.push(assay);seenAssay.add(assay);}
    if(value===null)missing++;else observed++;
  }
  const matrixCsv=csv(['feature_id',...assays],features.map((feature)=>[
    feature,...assays.map((assay)=>vals.get(JSON.stringify([feature,assay]))??'')
  ]));
  const metadata=metadataRowsFromRecords(parsed.rows,cols);
  // Keep analytical annotations independently of quantitative peak areas. For
  // long-form data they remain injection-specific, not falsely merged by feature.
  const annotationHeaders=heads.filter((header)=>annotationKeys.has(norm(header)));
  const annotationsCsv=annotationHeaders.length
    ? csv(['feature_id','assay_id',...annotationHeaders],
        parsed.rows.map((row)=>[String(row[cols.feature]||'').trim(),
          String(row[cols.assay]||'').trim(),...annotationHeaders.map((header)=>row[header]??'')]))
    : null;
  const metadataCsv=metadata?csv(
    ['subject_id','sample_id','assay_id','omic','condition','timepoint','batch','technical_replicate','sample_type','injection_order'],
    metadata.map((m)=>Object.values(m))
  ):null;
  return {format:'long_ms_auc',matrixCsv,metadataCsv,annotationsCsv,
    annotationColumns:annotationHeaders,sourceColumns:heads,features:features.length,assays:assays.length,
    observed,missing,metadataGenerated:Boolean(metadata)};
}
