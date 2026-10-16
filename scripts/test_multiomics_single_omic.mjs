import assert from 'node:assert/strict';
import { readFile } from 'node:fs/promises';
import { File } from 'node:buffer';
import { parseDelimited, runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';
import { convertMsAucExport } from '../src/lib/multiomics/ms-auc-import.js';

const read = async (name) => new File(
  [await readFile(new URL('../static/multiomics/' + name, import.meta.url))],
  name, { type: 'text/csv' }
);
const originalMetadata = await read('demo_metadata.csv');
const parsed = parseDelimited(await originalMetadata.text());
const mapping = Object.fromEntries([
  'subject_id','sample_id','assay_id','omic','condition','timepoint','batch',
  'technical_replicate','outcome','survival_time','survival_event'
].map((key) => [key,key]));
const inputs = {
  transcriptomics: { file: await read('demo_transcriptomics.csv'), valueType: 'raw_counts' },
  proteomics: { file: await read('demo_proteomics.csv'), valueType: 'log_intensity' },
  metabolomics: { file: await read('demo_metabolomics.csv'), valueType: 'peak_area' }
};
for (const [layer, entry] of Object.entries(inputs)) {
  const files = { metadata: originalMetadata, transcriptomics:null,proteomics:null,metabolomics:null, [layer]:entry.file };
  const base = {
    files, metadataRows:parsed.rows, columnMapping:mapping,
    dataTypes:{ [layer]:entry.valueType }, useReactome:false, resolveIdentifiers:false
  };
  const commonProtocol = {
    organism:'human',designType:'repeated',longitudinal:true,objective:'time',
    groupCount:'2',sampleOverlap:'same_specimen',batchKnown:'yes'
  };
  const differential = await runDeterministicAnalysis({ ...base, protocol:commonProtocol });
  assert.equal(differential.engine.analysisMode,'single_omic',layer);
  assert.deepEqual(differential.engine.analysedLayers,[layer],layer);
  assert.deepEqual(Object.keys(differential.layers),[layer],layer);
  assert.equal(differential.supervisedIntegration,null,layer);
  assert.equal(differential.crossOmics.testedPairs,0,layer);
  assert.equal(differential.crossOmics.pairs.length,0,layer);
  assert.equal(differential.metadataSummary.subjects,8,layer);
  assert.equal(differential.metadataSummary.overlap.allMatched,8,layer);
  assert.ok(differential.layers[layer].rows.length,layer);
  assert.ok(differential.layers[layer].qc?.featuresAfter > 0,layer);
  const exploration = await runDeterministicAnalysis({
    ...base,
    protocol:{...commonProtocol,objective:'explore',longitudinal:false,designType:'independent'}
  });
  assert.equal(exploration.engine.analysisMode,'single_omic',layer);
  assert.equal(exploration.exploration.analysisMode,'single_omic',layer);
  assert.equal(exploration.layers[layer].mode,'exploratory-single-omic-pca',layer);
  assert.ok(exploration.exploration.components.length > 0,layer);
  assert.equal(exploration.supervisedIntegration,null,layer);
}

const wide=convertMsAucExport([
  'feature_id,mz,rt,adduct,compound_name,A01,A02',
  'PK001,123.45,2.4,[M+H]+,Unidentified,1234,2345',
  'PK002,555.55,4.7,[M-H]-,Candidate,333,444'
].join('\n'));
assert.equal(wide.format,'wide_ms_auc');
assert.equal(wide.features,2);
assert.ok(wide.annotationsCsv.includes('123.45'));
assert.ok(wide.annotationsCsv.includes('Unidentified'));
assert.ok(wide.annotationsCsv.includes('PK002'));
assert.ok(!wide.matrixCsv.includes('123.45'));
assert.ok(!wide.matrixCsv.includes('compound_name'));
assert.deepEqual(wide.annotationColumns,['mz','rt','adduct','compound_name']);
const long=convertMsAucExport([
  'feature_id,assay_id,auc,mz,rt,adduct,condition,subject_id,sample_id',
  'PK001,A01,1234,123.45,2.4,[M+H]+,A,S1,S1',
  'PK001,A02,2345,123.46,2.5,[M+H]+,B,S2,S2'
].join('\n'));
assert.equal(long.format,'long_ms_auc');
assert.equal(long.features,1);
assert.equal(long.assays,2);
assert.equal(long.metadataGenerated,true);
assert.ok(long.annotationsCsv.includes('123.46'));
assert.ok(long.annotationsCsv.includes('A02'));
assert.ok(!long.matrixCsv.includes('123.45'));
assert.throws(()=>convertMsAucExport([
  'feature_id,assay_id,auc','PK001,A01,100','PK001,A01,200'
].join('\n')), /duplicate feature/);
console.log('Single-omics RNA/proteomics/metabolomics exploration + differential and MS provenance: PASS');
