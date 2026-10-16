import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import {
  createPortableProject,parsePortableProject,
  sha256Content,PORTABLE_PROJECT_FORMAT,PORTABLE_PROJECT_SCHEMA
} from '../src/lib/multiomics/portable-project.js';

if(!globalThis.File)globalThis.File=File;
const meta = new File([
  'subject_id;sample_id;assay_id;omic;condition\n' +
  'α01;S01;A01;proteomics;control\n' +
  'α02;S02;A02;proteomics;treated\n'
],'original_metadata.csv',{type:'text/csv'});
const prot = new File(['feature_id,A01,A02\r\nP0001,12.0,13.2\r\nP0002,5.0,7.1\r\n'],
  'original_proteins.csv',{type:'text/csv'});
const raw = new File(['row-id\tarea\tqc_label\nA01\t14532\tQC1\n'],
  'original-MS-vendor-export.tsv',{type:'text/tab-separated-values'});
const bundle = await createPortableProject({
  files:{metadata:meta,proteomics:prot,transcriptomics:null,metabolomics:null},
  extraFiles:{msOriginal:raw},
  settings:{
    studyName:'Confidential study',objective:'groups',organism:'human',
    selectedCovariates:['age'],columnMapping:{subject_id:'subject_id'},
    transcriptomicsValues:'raw_counts',proteomicsValues:'log_intensity',
    demoLoaded:true,referenceBackendUrl:'http://private-server:8787',
    useReactome:true,resolveIdentifiers:true
  },
  result:{
    layers:{proteomics:{rows:[{feature:'P0001',effect:0.25,pValue:0.04}]}},
    scientificAssurance:{confirmatoryReadiness:{certified:false}}
  }
});
assert.equal(bundle.format,PORTABLE_PROJECT_FORMAT);
assert.equal(bundle.schemaVersion,PORTABLE_PROJECT_SCHEMA);
assert.equal(bundle.settings.demoLoaded,true);
assert.ok(!JSON.stringify(bundle.settings).includes('private-server'));
assert.ok(!('useReactome' in bundle.settings));
assert.ok(!('resolveIdentifiers' in bundle.settings));
const serialized=JSON.stringify(bundle);
assert.ok(!serialized.includes('private-server'));
const restored=await parsePortableProject(serialized);
assert.equal(restored.files.metadata.name,'original_metadata.csv');
assert.equal(await restored.files.metadata.text(),await meta.text());
assert.equal(await restored.files.proteomics.text(),await prot.text());
assert.equal(await restored.files.msOriginal.text(),await raw.text());
assert.equal(restored.settings.demoLoaded,true);
assert.equal(restored.archivedResultPresent,true);
assert.equal(restored.originalDigests.metadata,await sha256Content(await meta.text()));
assert.equal(restored.settings.columnMapping.subject_id,'subject_id');
assert.ok(!('archivedResult' in restored),'past results must not be rehydrated as new evidence');

const tampered=JSON.parse(serialized);
tampered.inputs.proteomics.content=tampered.inputs.proteomics.content.replace('12.0','999.0');
await assert.rejects(parsePortableProject(tampered),/SHA-256/);
const missing=JSON.parse(serialized);
delete missing.inputs.metadata;
await assert.rejects(parsePortableProject(missing),/sample sheet/);
const fake=JSON.parse(serialized);
fake.inputs.malware={content:'!',sizeBytes:1,sha256:'0'.repeat(64)};
await assert.rejects(parsePortableProject(fake),/unknown input category/);
const version=JSON.parse(serialized);
version.schemaVersion=999;
await assert.rejects(parsePortableProject(version),/Unsupported portable project/);
await assert.rejects(createPortableProject({files:{metadata:null},settings:{}}),/sample sheet/);
const long=JSON.parse(serialized);
long.inputs.metadata.sizeBytes+=1;
await assert.rejects(parsePortableProject(long),/SHA-256/);
console.log('multiomics portable project SHA-256 integrity, original inputs, privacy and guarded restore: PASS');
