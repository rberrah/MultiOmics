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
const [meta,rna,lipid,scoreText,modelText,featureText,trainText,
  panelModelText,panelFeatureText,panelTrainText]=await Promise.all([
  'metadata.csv','transcriptomics.csv','metabolomics.csv',
  'individual_signatures.csv','frozen_model.json','feature_definition.json',
  'training_marker_observations.csv',
  'frozen_panel_model.json','panel_feature_definition.json','training_panel_observations.csv'
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
const panelModel=JSON.parse(panelModelText),panelDefinition=JSON.parse(panelFeatureText);
const signatures=parseIndividualSignatureCsv(scoreText);
const singleScores=signatures.filter(x=>x.signature_id==='nutrimouse_liver_Cyp3a11_trainonly_z');
const panelScores=signatures.filter(x=>x.signature_id==='nutrimouse_four_gene_trainonly_composite');
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
assert.equal(signatures.length,20,'Two scores per each of ten real heldout animal assays');
assert.equal(singleScores.length,10);
assert.equal(panelScores.length,10);
assert.equal(new Set(signatures.map(x=>x.signature_id)).size,2);
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
for(const item of singleScores){
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
// Reproduce the SECOND R-fitted four-gene composite independently from
// original source features. Selection must use TRAINING variance only.
assert.equal(panelModel.method,'four_gene_train_only_variance_selection_and_scaling');
assert.equal(panelModel.pinned_source_commit,model.pinned_source_commit);
assert.equal(panelModel.phenotype_labels_used_for_selection_or_fitting,false);
assert.equal(panelModel.heldout_values_used_for_selection_or_fitting,false);
assert.equal(panelModel.clinical_or_PK_prediction_validated,false);
assert.deepEqual(panelModel.training_subject_ids,model.training_subject_ids);
assert.deepEqual(panelModel.heldout_subject_ids,model.heldout_subject_ids);
assert.equal(panelModel.training_measurements_sha256,hash(panelTrainText));
assert.deepEqual(panelDefinition.selected_gene_ids,panelModel.feature_ids);
assert.deepEqual(panelDefinition.feature_weights,[.25,.25,.25,.25]);
assert.equal(panelModel.feature_ids.length,4);
assert.equal(panelModel.feature_ids[0],definition.feature_id);
assert.equal(new Set(panelModel.feature_ids).size,4);
const trainingSubjects=model.training_subject_ids;
const heldoutSubjects=model.heldout_subject_ids;
const trainingAssays=trainingSubjects.map(id=>'RNA_'+id.replace('NUTRIMOUSE_',''));
const geneRows=original.rows.map(row=>{
  const values=Object.fromEntries(ids.map((assay,j)=>[assay,Number(row[j+1])]));
  const training=trainingAssays.map(assay=>values[assay]);
  assert(training.every(Number.isFinite),'All training gene measurements must be genuine and finite');
  const mean=training.reduce((a,b)=>a+b,0)/training.length;
  const sd=Math.sqrt(training.reduce((a,v)=>a+(v-mean)**2,0)/(training.length-1));
  return {name:row[0],values,mean,sd};
});
const eligible=geneRows.filter(g=>g.name!==definition.feature_id && Number.isFinite(g.sd) && g.sd>1e-12);
eligible.sort((a,b)=>b.sd-a.sd||(a.name<b.name?-1:a.name>b.name?1:0));
assert(eligible.length>=3,'Need at least 3 training-variable other genes');
assert.deepEqual(panelModel.feature_ids.slice(1),eligible.slice(0,3).map(g=>g.name),
  'Panel gene selection must match independent TRAIN-only variance ranking');
const panelRows=panelModel.feature_ids.map(geneName=>{
  const record=geneRows.find(g=>g.name===geneName);
  assert(record,'Frozen selected gene missing from original source');
  return record;
});
for(let i=0;i<4;i++){
  assert(Math.abs(panelRows[i].mean-panelModel.training_means[i])<1e-10,
    'Panel gene center was not computed exclusively from train animals');
  assert(Math.abs(panelRows[i].sd-panelModel.training_sample_sds[i])<1e-10,
    'Panel gene scale was not computed exclusively from train animals');
}
for(const item of panelScores){
  assert(heldoutSubjects.includes(item.subject_id));
  const suffix=item.subject_id.replace('NUTRIMOUSE_','');
  const assay='RNA_'+suffix;
  assert.equal(item.assay_id,assay);
  assert.equal(item.sample_id,'SPEC_RNA_'+suffix);
  assert.equal(item.omic,'transcriptomics');
  assert.equal(item.score_scale,'mean_of_four_training_centered_gene_z_scores');
  assert.equal(item.fit_scope,'training_frozen_test_projection');
  assert.equal(item.score_provenance,'native_r_individual_score');
  assert.equal(item.reference_model_sha256,hash(panelModelText));
  assert.equal(item.feature_definition_sha256,hash(panelFeatureText));
  const expected=panelRows.reduce((acc,g,i)=>
    acc+(g.values[assay]-panelModel.training_means[i])/panelModel.training_sample_sds[i],0)/4;
  assert(Math.abs(expected-item.score)<1e-9,
    'Genuine heldout multigene score disagrees with independent frozen Node calculation');
}
assert.equal(new Set([...singleScores,...panelScores].map(x=>x.subject_id)).size,10);
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
assert.equal(produced.manifest.validation.individual_signature_export.n_assay_signature_rows,20);
assert.deepEqual(produced.manifest.validation.individual_signature_export.fit_scopes,
  ['training_frozen_test_projection']);
const output=join(source,'exchange');
await mkdir(output,{recursive:true});
for(const entry of produced.entries){
  const rel=entry.name.replace(/^exchange\//,'');
  const target=join(output,rel);
  const slash=target.lastIndexOf('/');
  await mkdir(target.slice(0,slash),{recursive:true});
  await writeFile(target,entry.data);
}
console.log('REAL Nutrimouse: ten measured heldout animals, genuine frozen 1-gene + 4-gene train-only scores; independently verified all R fits and real ZIP; NOT biologically or PK validated');
