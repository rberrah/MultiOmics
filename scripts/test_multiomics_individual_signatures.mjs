import assert from 'node:assert/strict';
import {webcrypto} from 'node:crypto';
import {mkdirSync,writeFileSync} from 'node:fs';
import {resolve,join,dirname} from 'node:path';
import {buildDownstreamExchange} from '../src/lib/multiomics/downstream-exchange.js';
import {buildResearchBundle} from '../src/lib/multiomics/research-export-package.js';
if(!globalThis.crypto)globalThis.crypto=webcrypto;

// Explicitly synthetic, per-assay scoring outputs for testing the SCHEMA only.
// This does NOT attest that a fitted model exists, or claim biological validity.
const metadata={content:[
 'assay_id,sample_id,subject_id,omic,tissue,timepoint',
 'RNA1,S01,P01,transcriptomics,blood,T0',
 'RNA2,S02,P02,transcriptomics,blood,T12'
].join('\n')+'\n',sha256:'a'.repeat(64)};
const project={
 format:'pmx-multiomics-portable-project',schemaVersion:1,archivedResult:null,
 inputs:{metadata,transcriptomics:{content:'feature_id,RNA1,RNA2\nENSG1,10,20\n'}},
 settings:{organism:'human',unitType:'participant',transcriptomicsValues:'raw_counts',transcriptomicsIdType:'ensembl_gene'}
};
const id='source-'+metadata.sha256.slice(0,16);
const base={
 subject_id:'P01',sample_id:'S01',assay_id:'RNA1',omic:'transcriptomics',
 signature_id:'CYP3A4_frozen_training_model',score:0.75,score_scale:'z_score',
 fit_scope:'training_frozen_test_projection',
 reference_model_sha256:'b'.repeat(64),
 feature_definition_sha256:'c'.repeat(64),
 source_exchange_analysis_id:id,
 score_provenance:'external_frozen_projection'
};
const scores=[base,{...base,subject_id:'P02',sample_id:'S02',assay_id:'RNA2',score:0.25,fit_scope:'out_of_fold_training'}];
const valid=await buildDownstreamExchange(project,{individualSignatures:scores});
assert.equal(valid.manifest.files.individual_signatures.path,'individual_signatures.csv');
assert.equal(valid.manifest.validation.status,'exploratory');
assert.equal(valid.manifest.validation.individual_signature_export.status,'schema_checked_only_not_biologically_validated');
const csv=valid.entries.find(x=>x.name==='exchange/individual_signatures.csv').data;
assert(csv.includes('training_frozen_test_projection'));
assert(csv.includes('out_of_fold_training'));
assert(csv.includes('external_frozen_projection'));
assert(!csv.includes('log2FC'));
const noScore=await buildDownstreamExchange(project);
assert(!noScore.manifest.files.individual_signatures);
const result={
 engine:{name:'synthetic-contract-fixture-only',version:'0'},
 protocol:{objective:'explore'},
 layers:{transcriptomics:{rows:[{feature:'ENSG1',effect:1,pValue:null,qValue:null}]}},
 individualSignatures:scores
};
const archive=await buildResearchBundle({result,portableProject:project});
assert.equal(archive.manifest.downstream_exchange.individual_signatures,'explicit_per_assay_schema_checked');
assert(archive.manifest.files['exchange/individual_signatures.csv']);
const badCases=[
 [{...base,subject_id:'P99'},/does not match/],
 [{...base,score:NaN},/finite number/],
 [{...base,fit_scope:'in_sample_transductive'},/Unfrozen/],
 [{...base,reference_model_sha256:'unknown'},/SHA-256/],
 [{...base,source_exchange_analysis_id:'other'},/analysis ID mismatch/],
 [{...base,qValue:0.001},/Group differential/],
 [{...base,score_scale:'group_logfc'},/normalized scoring scale/],
 [{...base,score_provenance:'group_effect'},/provenance/]
];
for(const [bad,re] of badCases){
 await assert.rejects(buildDownstreamExchange(project,{individualSignatures:[bad]}),re);
}
await assert.rejects(buildDownstreamExchange(project,{individualSignatures:[base,base]}),/Duplicate/);
const badResult={...result,individualSignatures:[{...base,score:Infinity}]};
await assert.rejects(buildResearchBundle({result:badResult,portableProject:project}),/finite number/);
const output=process.argv[2];
if(output){
 const target=resolve(output);
 mkdirSync(target,{recursive:true});
 for(const entry of valid.entries){
   const name=entry.name.replace(/^exchange\//,'');
   const dest=join(target,name);
   mkdirSync(dirname(dest),{recursive:true});
   writeFileSync(dest,entry.data);
 }
}
console.log('PASS: optional individual per-assay schema; exact sample linkage, scopes, SHA syntax, invalid-group leakage guards, ZIP inclusion');
