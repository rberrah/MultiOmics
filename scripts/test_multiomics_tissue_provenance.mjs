import assert from 'node:assert/strict';
import {webcrypto} from 'node:crypto';
import {buildDownstreamExchange} from '../src/lib/multiomics/downstream-exchange.js';
import {validateSourceLocation,summarizeSourceLocation} from '../src/lib/multiomics/tissue-provenance.js';
if(!globalThis.crypto)globalThis.crypto=webcrypto;

const columns=[
 'assay_id','sample_id','subject_id','omic','tissue','specimen_type',
 'tissue_ontology_id','cell_type_or_bulk','cell_composition_status',
 'cell_composition_method','sample_matrix','assay_platform','source_preprocessing',
 'pharmacology_target_organ','relative_to_dose','source_protocol_id'
];
const rows=[
 ['A1','S1','P01','transcriptomics','blood','PBMC','UBERON:0013756','bulk',
  'estimated','CIBERSORTx','PBMC','RNA-seq','CPM','liver','pre_dose','protocol1'],
 ['A2','S2','P02','metabolomics','blood','plasma','UBERON:0013756','cell_free',
  'not_measured','','plasma','LC-MS','log intensity','liver','pre_dose','protocol2'],
 ['A3','S3','P03','transcriptomics','liver','liver_biopsy','UBERON:0002107','bulk',
  'unknown','','liver tissue','RNA-seq','TPM','liver','not_applicable','protocol3']
];
const csv=[columns.join(','),...rows.map(r=>r.join(','))].join('\n')+'\n';
const root={
 inputs:{
  metadata:{content:csv,sha256:'a'.repeat(64)},
  transcriptomics:{content:'feature_id,A1,A3\nCYP3A4,12,120\n'},
  metabolomics:{content:'feature_id,A2\nBILE_ACID_1,0.1\n'}
 },
 settings:{organism:'human',unitType:'participant',
   transcriptomicsValues:'normalized_log2',metabolomicsValues:'normalized_log2'}
};
const exchange=await buildDownstreamExchange(root);
assert.equal(exchange.manifest.schema_version,'1.0.0');
assert.equal(exchange.manifest.validation.source_location_provenance.n_assays,3);
assert.equal(exchange.manifest.validation.source_location_provenance.n_cross_organ_proxy_not_validated,2,
 'PBMC and plasma should NOT be presented as direct human liver CYP activity');
assert.equal(exchange.manifest.validation.source_location_provenance.biological_source_to_target_organ_translation_validated,false);
assert.equal(exchange.manifest.validation.source_location_provenance.n_explicit_predose,2);
assert.equal(exchange.manifest.validation.source_location_provenance.n_cell_composition_measured_or_estimated,1);
assert.equal(exchange.manifest.validation.source_location_provenance.status,
 'source_specimen_and_cellular_origin_declared_not_functionally_validated');
const file=exchange.entries.find(e=>e.name==='exchange/samples.csv').data;
const headers=file.split(/\r?\n/)[0].split(',');
for(const col of columns)assert(headers.includes(col),'the original source column '+col+' must survive downstream exchange');
assert(file.includes('UBERON:0002107'));
assert(file.includes('UBERON:0013756'));
assert(file.includes('pre_dose'));
const basic=structuredClone(root);
basic.inputs.metadata.content=[
 'assay_id,sample_id,subject_id,omic',
 'A1,S1,P01,transcriptomics',
 'A3,S3,P03,transcriptomics'].join('\n')+'\n';
basic.inputs.metabolomics=null;
basic.inputs.transcriptomics.content='feature_id,A1,A3\nCYP3A4,12,120\n';
const missing=await buildDownstreamExchange(basic);
assert.equal(missing.summary.source_location_provenance.status,
 'incomplete_source_provenance_no_organ_proxy_claim');
assert.equal(missing.summary.source_location_provenance.n_original_tissue_declared,0);
assert.equal(missing.summary.source_location_provenance.n_explicit_predose,0);
assert.equal(missing.summary.source_location_provenance.direct_enzyme_or_transporter_activity_inferred,false);
assert.throws(()=>validateSourceLocation({cell_type_or_bulk:'PBMC'}),/cell_type_or_bulk/);
assert.throws(()=>validateSourceLocation({tissue_ontology_id:'liver'}),/ontology/);
assert.throws(()=>validateSourceLocation({relative_to_dose:'T0'}),/relative_to_dose/);
assert.throws(()=>validateSourceLocation({cell_composition_status:'unknown',cell_composition_method:'deconvolution'}),/Cell composition/);
const clash=structuredClone(root);
const cc=[columns.join(','), ...rows.map(r=>r.join(',')),
 ['A4','S1','P01','proteomics','liver','liver_biopsy','UBERON:0002107','bulk',
 'unknown','','liver','mass_spectrometry','normalized','liver','pre_dose','protocol4'].join(',')].join('\n')+'\n';
clash.inputs.metadata.content=cc;
clash.inputs.proteomics={content:'feature_id,A4\nCYP3A4,1\n'};
await assert.rejects(buildDownstreamExchange(clash),/Conflicting tissue/);
console.log('PASS: sample origin tissue/ontology/cell composition/PK time preserved, source mismatch rejected, blood-to-liver CYP proxy unvalidated');
