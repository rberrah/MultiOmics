// @ts-nocheck
// Friendly preflight only. The local R backend must independently revalidate
// original assay CSVs, exact pairs and statistical rank; this frontend preview
// never declares an experiment statistically certified.
const canon=(x)=>String(x??'').trim();
const isRNA=(x)=>['transcriptomics','transcriptome','rna','rnaseq','rna_seq','rna-seq','rna sequencing'].includes(
  canon(x).toLowerCase().replace(/\s+/g,' '));
const biological=(x)=>['biological','biological_sample','sample'].includes(canon(x).toLowerCase());
function result(code,extras={}){return {ok:code==='ok',code,...extras}}

/** Preview independent subject numbers and pairing using mapped source columns.
 * Repeated technical assays can share sample_id, but distinct paired conditions
 * MUST have distinct biological sample IDs. */
export function auditPairedRnaStudy({rows,mapping,matrixAssays,valueType}){
  if(valueType!=='raw_counts')return result('raw_counts_only');
  if(!Array.isArray(rows)||!rows.length)return result('missing_metadata');
  const fields=['subject_id','sample_id','assay_id','condition','omic','sample_type'];
  if(fields.some(k=>!canon(mapping?.[k])))return result('missing_columns');
  const assays=new Set(matrixAssays||[]);
  if(!assays.size||assays.size!==(matrixAssays||[]).length)return result('missing_or_duplicate_assays');
  const scoped=rows.filter(row=>isRNA(row[mapping.omic]));
  if(!scoped.length)return result('missing_transcriptomics_rows');
  const observed=new Set(),cells=new Map(),biologicalSamples=new Map(),groups=new Set();
  for(const row of scoped){
    if(!biological(row[mapping.sample_type]))return result('non_biological_sample');
    const subject=canon(row[mapping.subject_id]),sample=canon(row[mapping.sample_id]),
      assay=canon(row[mapping.assay_id]),condition=canon(row[mapping.condition]);
    if(!subject||!sample||!assay||!condition)return result('missing_pairing_value');
    if(observed.has(assay))return result('duplicate_assay');
    observed.add(assay);groups.add(condition);
    const key=JSON.stringify([subject,condition]);
    const old=cells.get(key);
    if(old!==undefined&&old!==sample)return result('duplicate_pair_cell');
    cells.set(key,sample);
    const previous=biologicalSamples.get(sample);
    if(previous!==undefined&&previous!==key)return result('shared_sample_between_conditions');
    biologicalSamples.set(sample,key);
  }
  if(observed.size!==assays.size||[...observed].some(id=>!assays.has(id)))
    return result('matrix_metadata_mismatch');
  const groupNames=[...groups].sort();
  if(groupNames.length!==2)return result('two_conditions_only');
  const subjects=[...new Set(scoped.map(row=>canon(row[mapping.subject_id])))];
  if(subjects.length<4)return result('min_four_pairs',{nPairs:subjects.length});
  for(const subject of subjects)
    for(const condition of groupNames)
      if(!cells.has(JSON.stringify([subject,condition])))
        return result('missing_pair',{nPairs:subjects.length});
  if(biologicalSamples.size!==subjects.length*2)
    return result('invalid_paired_samples');
  return result('ok',{
    nPairs:subjects.length,
    nBiologicalSamples:biologicalSamples.size,
    nAssays:assays.size,
    conditionNames:groupNames
  });
}

export function validPrimaryRnaContrast(groups,comparison,reference){
  return Array.isArray(groups)&&groups.length===2&&
    comparison!==reference&&Boolean(comparison)&&Boolean(reference)&&
    groups.includes(comparison)&&groups.includes(reference);
}

/** Fictional documentation scaffold, NEVER a real patient/sample source. */
export function pairedExampleSampleSheetCsv(){
  const header=['subject_id','sample_id','assay_id','omic','condition','timepoint','batch','technical_replicate','sample_type'];
  const records=[header];
  for(let i=1;i<=4;i++){
    for(const [suffix,group] of [['A','control'],['B','treated']]){
      const subject='DEMO_SUBJECT_'+String(i).padStart(2,'0');
      const id=subject+'_'+suffix;
      records.push([subject,id,id,'transcriptomics',group,'T0','RUN1','1','biological']);
    }
  }
  return records.map(row=>row.join(',')).join('\n')+'\n';
}
