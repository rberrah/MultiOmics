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
