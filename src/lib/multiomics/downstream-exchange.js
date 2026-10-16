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
export async function buildDownstreamExchange(portableProject){
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
      commit:import.meta.env?.VITE_MULTIOMICS_SOURCE_SHA||'unknown',
      r_version:null,python_version:null,dependency_lockfiles:[]
    },
    validation:{status:'exploratory',methods:[]}
  };
  entries.push({name:'exchange/manifest.json',data:JSON.stringify(manifest,null,2)+'\n'});
  return {entries,manifest};
}
