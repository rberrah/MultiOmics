import assert from 'node:assert/strict';
import {webcrypto} from 'node:crypto';
import {buildDownstreamExchange} from '../src/lib/multiomics/downstream-exchange.js';

if (!globalThis.crypto) globalThis.crypto=webcrypto;
const meta=[
  'assay_id,sample_id,subject_id,omic,tissue,timepoint,technical_replicate,time_origin',
  'R1,S1,P01,transcriptomics,blood,T0,run1,',
  'R1rep,S1,P01,transcriptomics,blood,T0,run2,',
  'R2,S2,P01,transcriptomics,blood,T12,run1,',
  'P1,S1,P01,proteomics,blood,T0,,',
  'P2,S2,P01,proteomics,blood,T12,,'
].join('\n')+'\n';
const root={
  inputs:{
    metadata:{content:meta,sha256:'a'.repeat(64)},
    transcriptomics:{content:'feature_id,R1,R1rep,R2\nENSG01,2,3,5\n'},
    proteomics:{content:'feature_id,P1,P2\nP123,0.3,0.5\n'}
  },
  settings:{
    organism:'human',unitType:'participant',
    transcriptomicsValues:'raw_counts',proteomicsValues:'log_intensity'
  }
};
const exchange=await buildDownstreamExchange(root);
assert.equal(exchange.manifest.schema_version,'1.0.0','V1 remains readable');
assert.equal(exchange.summary.status,'source_exchange_ready_not_pk_ready');
assert.equal(exchange.summary.structural_exchange_valid,true);
assert.equal(exchange.summary.pk_ready,false);
assert.match(exchange.summary.interpretation,/canonicalized, not biologically renormalized/i);
assert.equal(exchange.manifest.matrices.length,2);
const sampleFile=exchange.entries.find(e=>e.name==='exchange/samples.csv');
const [headerLine,...body]=sampleFile.data.trim().split(/\r?\n/);
const header=headerLine.split(',');
const records=body.map(line=>Object.fromEntries(header.map((col,i)=>[col,line.split(',')[i]])));
assert.deepEqual(records.map(x=>x.timepoint),['T0','T0','T12','T0','T12']);
assert(records.every(x=>x.collection_time===''),'categorical T12 is not fabricated numeric hours');
assert(records.every(x=>x.collection_time_unit===''),'no fabricated time units');
assert(records.every(x=>x.time_origin===''),'unknown dose/inclusion anchor is empty');
assert.deepEqual(records.slice(0,3).map(x=>x.technical_replicate_id),['run1','run2','run1']);
assert.equal(new Set(records.map(x=>x.subject_id)).size,1,'repeated assays are not new subjects');
assert.equal(new Set(records.map(x=>x.sample_id)).size,2,'shared cross-omic specimens');
assert.equal(records[0].timepoint,records[3].timepoint);
const digest=Buffer.from(await crypto.subtle.digest('SHA-256',new TextEncoder().encode(sampleFile.data))).toString('hex');
assert.equal(exchange.manifest.files.samples.sha256,digest);
assert.match(exchange.manifest.software.commit,/^(?:unknown|[a-f0-9]{40})$/);
const singleton={...root,inputs:{metadata:{content:'assay_id,sample_id,subject_id,omic,timepoint\nR1,S1,P01,transcriptomics,T12\n',sha256:'b'.repeat(64)},transcriptomics:{content:'feature_id,R1\nENSG01,2\n'}}};
const single=await buildDownstreamExchange(singleton);
assert.equal(single.manifest.matrices.length,1,'mono-omics must stay supported');
assert(single.entries.find(e=>e.name==='exchange/samples.csv').data.includes('T12'));
const withTime={...root,inputs:{...root.inputs,metadata:{...root.inputs.metadata,content:meta.replaceAll('blood,T12,','blood,T12,')}}};
withTime.inputs.metadata.content=withTime.inputs.metadata.content.replace('tissue,timepoint,technical_replicate,time_origin','tissue,timepoint,technical_replicate,time_origin,collection_time,collection_time_unit')
 .replaceAll('run1,\n','run1,,\n').replaceAll('run2,\n','run2,,\n');
const conflict={...root,inputs:{...root.inputs,metadata:{...root.inputs.metadata,content:meta.replace('P1,S1,P01,proteomics,blood,T0','P1,S1,P01,proteomics,blood,T12')}}};
await assert.rejects(buildDownstreamExchange(conflict),/Conflicting timepoint/);
const anchored={...root,inputs:{...root.inputs,metadata:{...root.inputs.metadata,
  content:meta.replace('time_origin\n','time_origin\n')
    .replaceAll('blood,T0,run1,','blood,T0,run1,after_inclusion')
    .replaceAll('blood,T0,run2,','blood,T0,run2,after_inclusion')
    .replaceAll('blood,T12,run1,','blood,T12,run1,after_inclusion')
    .replaceAll('blood,T0,,','blood,T0,,after_inclusion')
    .replaceAll('blood,T12,,','blood,T12,,after_inclusion')
}}};
assert((await buildDownstreamExchange(anchored)).entries.find(x=>x.name==='exchange/samples.csv').data.includes('after_inclusion'));
console.log('PASS: issue #20 categorical time, assay-vs-subject cardinality, mono/multi-omic, SHA-256 and negative guards');


// The website accepts mapped source columns, TSV/semicolon input and transposed
// matrices. The Omics2PK exchange must normalize the SAME accepted study
// rather than becoming a narrower hidden export contract.
const mappedMeta=[
  'patient_id\tspecimen\trun\tmodality\ttissue_name\tvisit\tvisit_unit\tanchor',
  'P10\tS10\tA10\trna\tblood\t0\th\tfirst_dose',
  'P11\tS11\tA11\trna\tblood\t12\th\tfirst_dose'
].join('\n')+'\n';
const mappedProject={
  inputs:{
    metadata:{content:mappedMeta,sha256:'d'.repeat(64)},
    transcriptomics:{content:'assay;ENSG_A;ENSG_B\nA10;10;\nA11;12;15\n'}
  },
  settings:{
    organism:'human',unitType:'participant',
    transcriptomicsValues:'normalized',transcriptomicsIdType:'ensembl_gene',
    columnMapping:{
      subject_id:'patient_id',sample_id:'specimen',assay_id:'run',omic:'modality',
      tissue:'tissue_name',collection_time:'visit',collection_time_unit:'visit_unit',
      time_origin:'anchor'
    }
  }
};
const normalized=await buildDownstreamExchange(mappedProject);
assert.equal(normalized.manifest.schema_version,'1.0.0');
assert.equal(normalized.manifest.matrices[0].orientation,'feature_rows_sample_columns');
assert.equal('source_orientation' in normalized.manifest.matrices[0],false,
  'source orientation must not break the strict V1 manifest schema');
assert.equal(normalized.summary.source_orientations.transcriptomics,'assay_rows_feature_columns');
assert.equal(normalized.summary.independent_subjects,2);
assert.equal(normalized.summary.unique_specimens,2);
assert.equal(normalized.summary.numeric_collection_time_complete,true);
assert.equal(normalized.summary.pk_temporal_numeric_compatible,true);
assert.equal(normalized.summary.explicit_time_origin_complete,true);
assert.equal(normalized.summary.tissue_present,true);
const normalizedSamples=normalized.entries.find(e=>e.name==='exchange/samples.csv').data;
assert.match(normalizedSamples,/A10,S10,P10,transcriptomics,blood,0,h/);
assert.match(normalizedSamples,/A11,S11,P11,transcriptomics,blood,12,h/);
const normalizedMatrix=normalized.entries.find(e=>e.name==='exchange/matrices/transcriptomics.csv').data;
assert.equal(normalizedMatrix,
  'feature_id,A10,A11\r\nENSG_A,10,12\r\nENSG_B,NA,15\r\n');
assert.equal(normalized.manifest.matrices[0].missing_value,'NA');
assert.equal(normalized.summary.pk_dose_concentration_data_included,false);
assert.match(normalized.summary.interpretation,/dose\/concentration/i);
assert.equal(normalized.manifest.validation.downstream_pk_readiness.status,'source_exchange_only_not_pk_ready');
assert.equal(normalized.manifest.validation.downstream_pk_readiness.structural_exchange_valid,true);
assert.equal(normalized.manifest.validation.downstream_pk_readiness.pk_temporal_numeric_compatible,true);
assert.equal(normalized.manifest.validation.downstream_pk_readiness.shared_molecular_pk_time_origin_verified,false);
assert.equal(normalized.manifest.validation.downstream_pk_readiness.clinical_or_pk_effect_validated,false);

const invalidUnitProject=structuredClone(mappedProject);
invalidUnitProject.inputs.metadata.content=mappedMeta.replaceAll('\th\tfirst_dose','\tfortnight\tfirst_dose');
const invalidUnit=await buildDownstreamExchange(invalidUnitProject);
assert.equal(invalidUnit.summary.numeric_collection_time_complete,true,
  'numeric-looking time plus arbitrary unit remains syntactically complete');
assert.equal(invalidUnit.summary.pk_temporal_numeric_compatible,false,
  'unsupported time units must never be advertised as Omics2PK-compatible');
assert.equal(invalidUnit.manifest.validation.downstream_pk_readiness.pk_temporal_numeric_compatible,false);

const negativeTimeProject=structuredClone(mappedProject);
negativeTimeProject.inputs.metadata.content=mappedMeta.replace('\t0\th\tfirst_dose','\t-1\th\tfirst_dose');
const negativeTime=await buildDownstreamExchange(negativeTimeProject);
assert.equal(negativeTime.summary.numeric_collection_time_complete,true);
assert.equal(negativeTime.summary.pk_temporal_numeric_compatible,false,
  'negative molecular times must not pass the PK consumer preflight');
assert.equal(negativeTime.manifest.validation.downstream_pk_readiness.pk_temporal_numeric_compatible,false);

const jsOnlyNumericProject=structuredClone(mappedProject);
jsOnlyNumericProject.inputs.metadata.content=mappedMeta.replace('\t0\th\tfirst_dose','\t0x10\th\tfirst_dose');
const jsOnlyNumeric=await buildDownstreamExchange(jsOnlyNumericProject);
assert.equal(jsOnlyNumeric.summary.numeric_collection_time_complete,true,
  'JS Number accepts hexadecimal syntax, exposing why the stricter consumer check is required');
assert.equal(jsOnlyNumeric.summary.pk_temporal_numeric_compatible,false,
  'hexadecimal JS-only number syntax must not be advertised as Python Omics2PK-compatible');


const ambiguousProject={
  inputs:{
    metadata:{content:'subject_id,sample_id,assay_id,omic\nP1,S1,A1,transcriptomics\nP2,S2,A2,transcriptomics\n',sha256:'e'.repeat(64)},
    transcriptomics:{content:'feature_id,A1,A2\nA1,1,2\nA2,3,4\n'}
  },
  settings:{organism:'human',unitType:'participant',
    transcriptomicsValues:'normalized',transcriptomicsIdType:'unknown'}
};
await assert.rejects(
  ()=>buildDownstreamExchange(ambiguousProject),
  /Ambiguous transcriptomics matrix orientation/,
  'an exchange must never guess when assay identifiers match both axes'
);
