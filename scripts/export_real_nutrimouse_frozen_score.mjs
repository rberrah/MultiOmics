#!/usr/bin/env node
// Independent Node verification of real R-fitted Nutrimouse heldout signatures.
// Recomputes train-only scaling from the ORIGINAL archived transcriptomic matrix,
// checks frozen model/definition hashes, then exercises the production exporter
// and its strict nonclinical exploratory manifest. No synthetic biological scores.
import assert from 'node:assert/strict';
import {readFile,mkdir,writeFile} from 'node:fs/promises';
import {createHash,webcrypto} from 'node:crypto';
import {join} from 'node:path';
import {buildDownstreamExchange,parseIndividualSignatureCsv} from '../src/lib/multiomics/downstream-exchange.js';
if(!globalThis.crypto) globalThis.crypto=webcrypto;

const source=process.argv[2];
if(!source)throw Error('Usage: node scripts/export_real_nutrimouse_frozen_score.mjs real_R_output_dir');
const read=async name=>readFile(join(source,name),'utf8');
const hash=value=>createHash('sha256').update(value).digest('hex');
const [meta,rna,lipid,scoreText,modelText,featureText,trainText]=await Promise.all([
  'metadata.csv','transcriptomics.csv','metabolomics.csv',
  'individual_signatures.csv','frozen_model.json','feature_definition.json',
  'training_marker_observations.csv'
].map(read));

// Source CSV is written by base-R's write.csv with proper RFC 4180 quoting.
// Re-parse independently rather than relying solely on the exporter contract.
function parse(text){
  const rows=[];let row=[],cell='',quoted=false;
  for(let i=0;i<text.length;i++){
    const ch=text[i];
    if(quoted){
      if(ch==='"' && text[i+1]==='"'){cell+='"';i++;}
      else if(ch==='"')quoted=false;
      else cell+=ch;
    }else if(ch==='"' && !cell)quoted=true;
    else if(ch===','){row.push(cell);cell='';}
    else if(ch==='\n'){row.push(cell.replace(/\r$/,''));rows.push(row);row=[];cell='';}
    else cell+=ch;
  }
  if(quoted)throw Error('Unclosed input CSV quote');
  if(cell||row.length){row.push(cell);rows.push(row);}
  assert(rows.length>1);
  const [header,...data]=rows;
  assert.equal(header.length,new Set(header).size);
  for(const r of data)assert.equal(r.length,header.length,'R CSV ragged input');
  return {header,rows:data};
}
const model=JSON.parse(modelText),definition=JSON.parse(featureText);
const signatures=parseIndividualSignatureCsv(scoreText);
const original=parse(rna),originalMetadata=parse(meta),training=parse(trainText);
const ids=original.header.slice(1);
const markerRow=original.rows.filter(row=>row[0]===definition.feature_id);
assert.equal(markerRow.length,1,'Exact prespecified real gene feature required');
assert.equal(model.feature_id,definition.feature_id);
assert.equal(model.method,'one_feature_train_only_standardization');
assert.equal(model.phenotype_labels_used_for_fitting,false);
assert.equal(model.heldout_values_used_for_fitting,false);
assert.equal(model.train_n,30);
assert.equal(model.heldout_n,10);
assert.equal(model.pinned_source_commit,'ef3e760526623d9e91e945e4b50be48d764efc40');
assert.equal(model.training_measurements_sha256,hash(trainText));
assert.equal(signatures.length,10,'Ten genuine heldout assays, never group-level scores');
assert.equal(model.training_subject_ids.length,30);
assert.equal(model.heldout_subject_ids.length,10);
assert.equal(new Set([...model.training_subject_ids,...model.heldout_subject_ids]).size,40,
             'Biological subjects must never cross training and heldout');
assert.equal(originalMetadata.rows.length,80,'Two assay types, 40 original animals');
assert.equal(original.rows.length,120);
assert.equal(parse(lipid).rows.length,21);
const toOriginal=Object.fromEntries(ids.map((id,i)=>[id,Number(markerRow[0][i+1])]));
const trainingVals=model.training_subject_ids.map((id)=>{
  const assay='RNA_'+id.replace('NUTRIMOUSE_','');
  const value=toOriginal[assay];
  assert(Number.isFinite(value),'Cannot reconstruct training measurement');
  return value;
});
const mean=trainingVals.reduce((a,b)=>a+b,0)/trainingVals.length;
const sd=Math.sqrt(trainingVals.reduce((a,v)=>a+(v-mean)**2,0)/(trainingVals.length-1));
assert(Math.abs(mean-model.training_mean)<1e-10,'Frozen R training mean mismatch');
assert(Math.abs(sd-model.training_sample_sd)<1e-10,'Frozen R training SD mismatch');
for(const item of signatures){
  assert(model.heldout_subject_ids.includes(item.subject_id));
  // Biological identity is the animal; do not claim shared RNA/lipid aliquot.
  assert.equal(item.sample_id,'SPEC_RNA_'+item.subject_id.replace('NUTRIMOUSE_',''));
  assert.equal(item.omic,'transcriptomics');
  assert.equal(item.assay_id,'RNA_'+item.subject_id.replace('NUTRIMOUSE_',''));
  assert.equal(item.score_scale,'training_centered_z_score');
  assert.equal(item.fit_scope,'training_frozen_test_projection');
  assert.equal(item.score_provenance,'native_r_individual_score');
  assert.equal(item.reference_model_sha256,hash(modelText));
  assert.equal(item.feature_definition_sha256,hash(featureText));
  const expected=(toOriginal[item.assay_id]-mean)/sd;
  assert(Math.abs(expected-item.score)<1e-9,'Heldout score did not use frozen train preprocessing');
}
const portable={
  format:'pmx-multiomics-portable-project',schemaVersion:1,
  inputs:{
    metadata:{content:meta,sha256:hash(meta)},
    transcriptomics:{content:rna},metabolomics:{content:lipid}
  },
  settings:{
    organism:'mouse',unitType:'animal',
    transcriptomicsValues:'normalized_expression_not_raw_RNAseq_counts',
    metabolomicsValues:'liver_fatty_acid_relative_percentages_not_LCMS_areas',
    transcriptomicsIdType:'original_nutrimouse_gene_symbol',
    metabolomicsIdType:'original_nutrimouse_fatty_acid_label'
  }
};
const produced=await buildDownstreamExchange(portable,{individualSignatures:signatures});
assert.equal(produced.manifest.analysis_id,'source-'+hash(meta).slice(0,16));
assert.equal(produced.manifest.validation.status,'exploratory');
assert.equal(produced.manifest.validation.individual_signature_export.status,
             'schema_checked_only_not_biologically_validated');
assert.equal(produced.manifest.validation.individual_signature_export.no_pk_effect_estimated,true);
assert.equal(produced.manifest.validation.individual_signature_export.n_assay_signature_rows,10);
const output=join(source,'exchange');
await mkdir(output,{recursive:true});
for(const entry of produced.entries){
  const rel=entry.name.replace(/^exchange\//,'');
  const target=join(output,rel);
  const slash=target.lastIndexOf('/');
  await mkdir(target.slice(0,slash),{recursive:true});
  await writeFile(target,entry.data);
}
console.log('REAL measured Nutrimouse: R 30-train/10-heldout frozen biological assay score independently rechecked, valid exploratory ZIP contract, no PK claim');
