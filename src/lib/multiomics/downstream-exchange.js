// @ts-nocheck
// Conservative, opt-in exchange from the EXACT original portable project inputs.
// Unsupported metadata/matrices skip the exchange; never relabel processed data as raw.
const TYPES=['transcriptomics','proteomics','metabolomics'];
const utf8=new TextEncoder();
function parseCsv(text){
  if(typeof text!=='string'||!text.trim())throw Error('Empty CSV');
  text=text.replace(/^\uFEFF/,'');
  const rows=[];let row=[],cell='',quoted=false;
  for(let i=0;i<text.length;i++){
    const ch=text[i];
    if(quoted){
      if(ch==='"'&&text[i+1]==='"'){cell+='"';i++;}
      else if(ch==='"')quoted=false;
      else cell+=ch;
    }else if(ch==='"'&&cell==='')quoted=true;
    else if(ch===','){row.push(cell);cell='';}
    else if(ch==='\n'){row.push(cell.replace(/\r$/,''));rows.push(row);row=[];cell='';}
    else cell+=ch;
  }
  if(quoted)throw Error('Unclosed CSV quote');
  if(cell||row.length){row.push(cell);rows.push(row);}
  if(!rows.length||rows[0].some(x=>!x)||new Set(rows[0]).size!==rows[0].length)throw Error('Invalid CSV header');
  if(rows.slice(1).some(x=>x.length!==rows[0].length))throw Error('Ragged CSV');
  return {header:rows[0],rows:rows.slice(1)};
}
function csvCell(v){
  const s=String(v??'');
  return /[",\r\n]/.test(s)?'"'+s.replaceAll('"','""')+'"':s;
}
function csv(rows){return rows.map(r=>r.map(csvCell).join(',')).join('\r\n')+'\r\n';}
async function digest(s){
  const b=await crypto.subtle.digest('SHA-256',utf8.encode(s));
  return [...new Uint8Array(b)].map(x=>x.toString(16).padStart(2,'0')).join('');
}
export async function buildDownstreamExchange(portableProject,{individualSignatures=null}={}){
  const input=portableProject?.inputs||{};
  const metadata=parseCsv(input.metadata?.content);
  const h=metadata.header;
  for(const key of ['assay_id','sample_id','subject_id','omic'])
    if(!h.includes(key))throw Error('Exchange requires explicit sample sheet column '+key);
  const sampleRows=metadata.rows.map(r=>{
    const sample=Object.fromEntries(h.map((key,i)=>[key,r[i]]));
    // Native demos call this technical_replicate; downstream contract uses _id.
    // Explicit downstream column wins; never synthesize a biological sample.
    sample.technical_replicate_id=sample.technical_replicate_id||sample.technical_replicate||'';
    return sample;
  });
  const specimenOwners=new Map();
  const specimenContext=new Map();
  const sampleIds=new Set(),sampleOutput=[[
    'assay_id','sample_id','subject_id','omic','tissue','collection_time',
    'collection_time_unit','condition','batch','technical_replicate_id','source_dataset_id',
    // Backwards-compatible appended columns: categorical time is NOT numeric PK time.
    'timepoint','timepoint_unit','time_origin'
  ]];
  for(const s of sampleRows){
    if(!TYPES.includes(s.omic)||!s.assay_id||!s.sample_id||!s.subject_id)
      throw Error('Invalid sample sheet for downstream exchange');
    const owner=specimenOwners.get(s.sample_id);
    if(owner!==undefined && owner!==s.subject_id)
      throw Error('One sample_id cannot belong to multiple biological subjects');
    specimenOwners.set(s.sample_id,s.subject_id);
    const previous=specimenContext.get(s.sample_id)||{};
    for(const field of ['tissue','collection_time','collection_time_unit','timepoint','timepoint_unit','time_origin']){
      if(previous[field] && s[field] && previous[field]!==s[field])
        throw Error('Conflicting '+field+' for shared specimen '+s.sample_id);
    }
    specimenContext.set(s.sample_id,{...previous,...Object.fromEntries(['tissue','collection_time','collection_time_unit','timepoint','timepoint_unit','time_origin'].map(field=>[field,s[field]||previous[field]||'']))});
    const key=s.omic+'\u0000'+s.assay_id;
    if(sampleIds.has(key))throw Error('Duplicate assay_id for '+s.omic);
    sampleIds.add(key);
    sampleOutput.push(sampleOutput[0].map(k=>s[k]||''));
  }
  const supplied=new Set(TYPES.filter(omic=>!!input[omic]));
  for(const sample of sampleRows)
    if(!supplied.has(sample.omic))throw Error('Sample sheet declares '+sample.omic+' without its original matrix');
  const entries=[],features=[['omic','feature_id','id_namespace','display_name','organism','taxon_id','annotation_database','annotation_version','tissue_scope','mapping_status']];
  const matrices=[];let total=0;
  async function add(path,content){
    const sha256=await digest(content);
    entries.push({name:'exchange/'+path,data:content});
    return {path,sha256};
  }
  const sampleSpec=await add('samples.csv',csv(sampleOutput));
  const featureKeys=new Set();
  const declared=portableProject.settings||{};
  for(const omic of TYPES){
    if(!input[omic])continue;
    const matrix=parseCsv(input[omic].content);
    if(matrix.header[0]!=='feature_id')throw Error('Original '+omic+' matrix must start with feature_id');
    const assayIds=matrix.header.slice(1);
    if(!assayIds.length||new Set(assayIds).size!==assayIds.length)
      throw Error('Duplicate or absent '+omic+' assay IDs');
    const expected=sampleRows.filter(s=>s.omic===omic).map(s=>s.assay_id);
    if(expected.length!==assayIds.length||assayIds.some(s=>!expected.includes(s)))
      throw Error('Sample sheet and '+omic+' matrix columns differ');
    const seen=new Set();
    const scale=declared[omic+'Values']||'source_values_unspecified';
    for(const r of matrix.rows){
      const id=r[0];
      if(!id||seen.has(id))throw Error('Missing/duplicate original feature in '+omic);
      seen.add(id);
      for(const v of r.slice(1)){
        if(v==='NA')continue;
        if(!v||!Number.isFinite(Number(v)))throw Error('Non-numeric/missing value in '+omic+'; use literal NA');
        if(omic==='transcriptomics'&&scale==='raw_counts'&&(!Number.isSafeInteger(Number(v))||Number(v)<0))
          throw Error('RNA raw_counts must be non-negative safe integers');
      }
      const k=omic+'\u0000'+id;
      if(featureKeys.has(k))throw Error('Feature collision');
      featureKeys.add(k);
      const idtype=declared[omic+'IdType']||'unresolved_source_namespace';
      features.push([omic,id,idtype,'',declared.organism||'','','','','','unresolved']);
    }
    const file=await add('matrices/'+omic+'.csv',csv([matrix.header,...matrix.rows]));
    matrices.push({omic,file,value_type:scale,orientation:'feature_rows_sample_columns',normalization:null,missing_value:'NA'});
    total++;
  }
  if(!total)throw Error('No original omic matrix available');
  const featuresSpec=await add('features.csv',csv(features));
  const manifest={
    schema_version:'1.0.0',
    analysis_id:'source-'+(input.metadata.sha256||'unknown').slice(0,16),
    study:{organism:declared.organism||'unknown',biological_unit:declared.unitType||'other',data_access:'private'},
    files:{samples:sampleSpec,features:featuresSpec},
    matrices,
    software:{
      producer:'MultiOmics browser source exchange',
      // Injected by Vite from the checked-out Git commit. Unknown is explicit
      // for unversioned browser/local builds and must not be forged.
      commit:(import.meta.env && import.meta.env.VITE_MULTIOMICS_SOURCE_SHA)||'unknown',
      r_version:null,python_version:null,dependency_lockfiles:[]
    },
    validation:{status:'exploratory',methods:[]}
  };
  // Optional *individual* biological scores from an explicit, traceable
  // upstream scoring result. Never synthesize such scores from differential
  // group statistics or invent fitted models from the original feature matrix.
  // Technical schema qualification != scientific/predictive validation.
  if(individualSignatures!==null){
    if(!Array.isArray(individualSignatures)||individualSignatures.length===0)
      throw Error('Individual signature export needs actual per-assay score rows');
    const columns=[
      'subject_id','sample_id','assay_id','omic','signature_id','score',
      'score_scale','fit_scope','reference_model_sha256',
      'feature_definition_sha256','source_exchange_analysis_id'
    ];
    const byAssay=new Map(sampleRows.map(s=>[s.omic+'\\u0000'+s.assay_id,s]));
    const seen=new Set();
    const output=[columns];
    for(const row of individualSignatures){
      if(!row||typeof row!=='object'||Array.isArray(row))
        throw Error('Invalid individual signature row');
      const assay=byAssay.get(row.omic+'\\u0000'+row.assay_id);
      if(!assay||row.subject_id!==assay.subject_id||row.sample_id!==assay.sample_id)
        throw Error('Individual signature subject/specimen/assay does not match original sample sheet');
      if(typeof row.signature_id!=='string'||!row.signature_id.trim()||
         !/^[A-Za-z0-9_.:/-]{1,160}$/.test(row.signature_id))
        throw Error('Invalid individual signature identifier');
      if(typeof row.score!=='number'||!Number.isFinite(row.score))
        throw Error('Individual score must be a measured finite number, not a missing or group contrast');
      if(typeof row.score_scale!=='string'||!row.score_scale.trim()||
         /^(raw_counts|group_logfc|log2foldchange|p_value|q_value|unknown)$/i.test(row.score_scale))
        throw Error('Individual signatures need an explicit per-assay normalized scoring scale');
      if(!['training_frozen_test_projection','out_of_fold_training'].includes(row.fit_scope))
        throw Error('Unfrozen/in-sample/transductive signature scopes are excluded');
      if(!['reference_model_sha256','feature_definition_sha256'].every(key=>
        typeof row[key]==='string'&&/^[a-f0-9]{64}$/.test(row[key])))
        throw Error('Frozen scoring reference and feature SHA-256 are mandatory');
      if(row.source_exchange_analysis_id!==manifest.analysis_id)
        throw Error('Individual signature exchange analysis ID mismatch');
      // The schema alone cannot prove that a model was genuinely frozen:
      // the caller must supply independently audited scoring provenance.
      if(!['native_r_individual_score','external_frozen_projection'].includes(row.score_provenance))
        throw Error('Individual score needs native assay score or frozen-model projection provenance');
      if(['pValue','qValue','effect','logFC','log2FC','padj'].some(key=>Object.hasOwn(row,key)))
        throw Error('Group differential statistics cannot masquerade as individual covariates');
      const unique=row.omic+'\\u0000'+row.assay_id+'\\u0000'+row.signature_id;
      if(seen.has(unique))throw Error('Duplicate individual signature per assay');
      seen.add(unique);
      output.push(columns.map(key=>row[key]));
    }
    manifest.files.individual_signatures=await add('individual_signatures.csv',csv(output));
    manifest.validation.individual_signature_export={
      status:'schema_checked_only_not_biologically_validated',
      fit_scopes:[...new Set(individualSignatures.map(r=>r.fit_scope))],
      n_assay_signature_rows:individualSignatures.length,
      source:'explicit_per_assay_scoring_result',
      no_pk_effect_estimated:true
    };
  }
  entries.push({name:'exchange/manifest.json',data:JSON.stringify(manifest,null,2)+'\n'});
  return {entries,manifest};
}
