import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { convertMsAucExport } from '../src/lib/multiomics/ms-auc-import.js';
import { parseDelimited, runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';
import { buildMultiomicsVisualizationData } from '../src/lib/multiomics/visualization-data.js';
import { focusMetabolicRegion } from '../src/lib/multiomics/metabolic-network.js';
import { parseMetaboliteAnnotations } from '../src/lib/multiomics/metabolite-annotation.js';

const row = (feature, effect, qValue, effectScale = 'log2') => ({ feature, effect, qValue, effectScale });
const result = {
  layers: {
    transcriptomics: {
      rows: [
        row('LDHA', -2.2, 0.001), row('LDHB', 1.4, 0.020),
        row('GAPDH', 0.7, 0.050), row('CS', 1.2, 0.030),
        row('HK1', 1.5, 0.020), row('HK2', -1.2, 0.030),
        row('ENSG_UNKNOWN', 3.5, 0.001),
        row('G6PD', 2, 0.005, 'logit')
      ]
    },
    metabolomics: {
      rows: [
        row('Lactate', -1.3, 0.03), row('Pyruvate', 0.5, 0.08),
        row('Citrate', 1.5, 0.01),
        row('Fructose-6-phosphate', 4, 0.001, 'untyped')
      ]
    }
  }
};

const visuals = await buildMultiomicsVisualizationData({
  files: {}, metadataRows: [], columnMapping: {}, dataTypes: {}, analysisResult: result
});
assert.equal(visuals.metabologramPathways.length, 4);
const glycolysis = visuals.metabologramPathways.find((x) => x.id === 'glycolysis');
const tca = visuals.metabologramPathways.find((x) => x.id === 'tca');
const ppp = visuals.metabologramPathways.find((x) => x.id === 'ppp');
assert.deepEqual(new Set(glycolysis.transcriptomics.map((x) => x.feature)),
  new Set(['LDHA', 'LDHB', 'GAPDH']));
assert.ok(glycolysis.metabolomics.some((x) => x.feature === 'Lactate'));
assert.ok(tca.metabolomics.some((x) => x.feature === 'Citrate'));
assert.ok(tca.transcriptomics.some((x) => x.feature === 'CS'));
assert.equal(ppp.transcriptomics.length, 0, 'A logit estimate must not be presented as a log2 fold change');
const hk = visuals.centralCarbon.enzymes.find((x) => x.id === 'hk');
assert.equal(hk.measurement.status, 'discordant_isoforms');
assert.equal(hk.measurement.effect, null, 'A gene-family summary cannot hide opposing effects');
assert.ok(!glycolysis.transcriptomics.some((x) => x.feature === 'HK1' || x.feature === 'HK2'));
assert.ok(!visuals.centralCarbon.metabolites.find((x) => x.id === '3pg').aliases.includes('2,3pg'), '3-PG is not 2,3-PG');
assert.ok(!glycolysis.metabolomics.some((x) => x.feature === 'Fructose-6-phosphate'));
assert.equal(glycolysis.coverageMetabolites, glycolysis.metabolomics.length);
assert.equal(glycolysis.meanTranscriptomicLog2Fc,
  glycolysis.transcriptomics.reduce((acc, x) => acc + x.effect, 0) / glycolysis.transcriptomics.length);
assert.ok(visuals.centralCarbon.edges.some(([a,b]) => a === 'pyruvate' && b === 'acetylcoa'));
assert.ok(visuals.centralCarbon.edges.some(([a,b]) => a === 'acetylcoa' && b === 'citrate'));
assert.ok(visuals.centralCarbon.edges.some(([a,b]) => a === 'glutamine' && b === 'glutamate'));
assert.ok(!visuals.centralCarbon.edges.some(([a,b]) => a === 'pyruvate' && b === 'citrate'));
assert.ok(!visuals.centralCarbon.edges.some(([a,b]) => a === 'g3p' && b === 'serine'));
assert.match(visuals.methodologicalBoundary.join(' '), /No metabolic flux/i);

// ChEBI identities are curated from EMBL-EBI rather than inferred from the ID
// text. Demo-like data must show lactate/succinate on the fixed map; tryptophan
// and kynurenine remain measured but outside this intentionally limited map.
const example = await buildMultiomicsVisualizationData({
  files: {}, metadataRows: [], columnMapping: {}, dataTypes: {},
  analysisResult: {
    layers: {
      transcriptomics: { rows: [
        row('GAPDH', -0.06, 0.4), row('STAT1', 1.01, 0.02),
        row('IDO1', 1.92, 0.01), row('KYNU', 1.09, 0.03),
        row('CYP3A5', -0.8, 0.01)
      ] },
      metabolomics: { rows: [
        row('CHEBI:16828', -1.2, 0.04), row('CHEBI:16946', 0.37, 0.04),
        row('CHEBI:24996', -0.37, 0.08), row('CHEBI:30031', -0.49, 0.03)
      ] }
    }
  }
});
assert.equal(example.centralCarbon.measuredMetabolites, 2,
  'ChEBI lactate and succinate must appear in central-carbon map');
assert.equal(example.centralCarbon.measuredTranscripts, 1,
  'Only GAPDH among these demo genes is part of fixed central-carbon schematic');
assert.equal(example.centralCarbon.metabolites.find((n) => n.id === 'lactate').measurement.feature, 'CHEBI:24996');
assert.equal(example.centralCarbon.metabolites.find((n) => n.id === 'succinate').measurement.feature, 'CHEBI:30031');
assert.deepEqual(new Set(example.centralCarbon.offMapMetabolites.map((n) => n.feature)),
  new Set(['CHEBI:16828', 'CHEBI:16946']));
assert.deepEqual(new Set(example.centralCarbon.offMapTranscripts.map((n) => n.feature)),
  new Set(['STAT1', 'IDO1', 'KYNU', 'CYP3A5']));
assert.equal(example.metabologram.metabolomics.find((n) => n.feature === 'CHEBI:16828').labelFr,
  'L-tryptophane');
assert.equal(example.metabologram.metabolomics.find((n) => n.feature === 'CHEBI:16946').labelFr,
  'L-kynurénine');
assert.deepEqual(example.metabologramPathways.find((n) => n.id === 'glycolysis').metabolomics.map((n) => n.feature),
  ['CHEBI:24996']);
assert.deepEqual(example.metabologramPathways.find((n) => n.id === 'tca').metabolomics.map((n) => n.feature),
  ['CHEBI:30031']);

const untouched = JSON.stringify(result);
const again = await buildMultiomicsVisualizationData({
  files: {}, metadataRows: [], columnMapping: {}, dataTypes: {}, analysisResult: result
});
assert.deepEqual(again, visuals, 'Visualization data must be deterministic');
assert.equal(JSON.stringify(result), untouched, 'Visualizations must not mutate analysis results');

// Expanded network: strict ID matching, region selection, informative non-links.
assert.ok(example.focusedMetabolicNetwork.coverage.metabolites >= 80);
assert.ok(example.focusedMetabolicNetwork.coverage.genes >= 90);
assert.ok(example.focusedMetabolicNetwork.coverage.regions >= 12);
assert.ok(example.focusedMetabolicNetwork.coverage.links >= 90);
const graph=example.focusedMetabolicNetwork;
const id=(i)=>graph.audit.find((r)=>r.input===i);
assert.equal(id('CHEBI:16828').status, 'verified_id');
assert.equal(id('CHEBI:16946').status, 'verified_id');
assert.equal(id('CHEBI:24996').status, 'verified_id');
assert.equal(id('CHEBI:30031').status, 'verified_id');
const auto=focusMetabolicRegion(graph,'auto',1);
assert.ok(auto.some((r)=>r.id==='kynurenine'));
assert.ok(auto.some((r)=>r.id==='glycolysis'));
assert.ok(auto.some((r)=>r.id==='tca'));
assert.ok(!auto.some((r)=>r.id==='serotonin'),
  'Do not show a second unrelated branch for the same measured precursor');
assert.ok(auto.find((r)=>r.id==='kynurenine').nodes.some((n)=>n.id==='formylkyn'));
assert.ok(!focusMetabolicRegion(graph,'kynurenine',0)[0].nodes.some((n)=>n.id==='formylkyn'));
assert.ok(focusMetabolicRegion(graph,'kynurenine',2)[0].nodes.length >
  focusMetabolicRegion(graph,'kynurenine',0)[0].nodes.length);

const strictResult=await buildMultiomicsVisualizationData({
  files:{},metadataRows:[],columnMapping:{},dataTypes:{},
  analysisResult:{layers:{metabolomics:{rows:[
    row('CHEBI:99999999',1,0.01), row('unrecognizable metabolite',2,0.01),
    row('CHEBI:16828',1,0.01,'logit'),row('L-kynurenine',1,0.01)
  ]},transcriptomics:{rows:[]}}}
});
const strict=strictResult.focusedMetabolicNetwork;
assert.equal(strict.audit.find((r)=>r.input==='CHEBI:99999999').status,'unmapped_identifier');
assert.equal(strict.audit.find((r)=>r.input==='unrecognizable metabolite').status,'unrecognized');
assert.equal(strict.audit.find((r)=>r.input==='L-kynurenine').status,'name_only');
assert.equal(strict.audit.find((r)=>r.input==='CHEBI:16828').usable,false,
  'Identity evidence must never replace a missing log2 effect scale');
assert.equal(strict.nodes.filter((n)=>n.id==='tryptophan')[0].measurement,null);
assert.deepEqual(focusMetabolicRegion(graph,'auto',1),focusMetabolicRegion(graph,'auto',1),
  'Same input must produce same focused subgraph');


// MS peak AUC is an integrated chromatographic area, not a PK exposure AUC.
// A single file with group annotations must import without inventing other
// omics layers or averaging duplicate observations.
const aucRows = ['feature_id,assay_id,auc,subject_id,sample_id,condition'];
for (let i=0; i<8; i++) {
  for (const [id,control,treatment] of [
    ['CHEBI:24996',[102,107,98,110],[350,375,390,370]],
    ['CHEBI:30031',[430,410,448,425],[220,233,215,228]]
  ]) aucRows.push([
    id,'INJ'+(i+1),i<4?control[i]:treatment[i-4],'SUBJ'+(i+1),'S'+(i+1),i<4?'control':'treatment'
  ].join(','));
}
const msInput = aucRows.join('\n');
const converted=convertMsAucExport(msInput);
assert.equal(converted.format,'long_ms_auc');
assert.equal(converted.features,2);
assert.equal(converted.assays,8);
assert.equal(converted.observed,16);
assert.ok(converted.metadataGenerated);
const matrixParsed=parseDelimited(converted.matrixCsv);
assert.equal(matrixParsed.rows[0].INJ1,'102');
assert.equal(matrixParsed.rows[0].INJ8,'370');
assert.throws(()=>convertMsAucExport(msInput+'\n'+aucRows[1]),/duplicate feature/);
assert.throws(()=>convertMsAucExport('feature_id,assay_id,auc\nX,A,-2'),/non-negative/);
assert.equal(convertMsAucExport('feature_id,assay_id,auc\nX,A,1').metadataCsv,null,
  'Without explicit group metadata, do not guess conditions.');
const wide=convertMsAucExport('feature_id,m/z,Retention Time,MS1,MS2\nCHEBI:24996,90.0,2.3,1000,2000\nCHEBI:30031,118.0,3.4,800,600');
assert.equal(wide.format,'wide_ms_auc');
assert.equal(wide.assays,2);
assert.equal(wide.features,2);
assert.equal(wide.observed,4);
assert.deepEqual(parseDelimited(wide.matrixCsv).headers,['feature_id','MS1','MS2'],
  'm/z and retention time must never become sample measurement columns');
assert.throws(()=>convertMsAucExport('feature_id,MS1\nA,-10'),/non-negative/);
assert.throws(()=>convertMsAucExport('feature_id,MS1\nA,10\nA,11'),/duplicate feature/);
const vendor=convertMsAucExport('Metabolite Name,Sample Name,Peak Area AUC,Group\nLactate,S1,15320,control\nLactate,S2,23500,treated');
assert.equal(vendor.format,'long_ms_auc');
assert.equal(vendor.features,1);
assert.equal(vendor.assays,2);
assert.ok(vendor.metadataCsv?.includes('treated'));
const msMetadata=new File([converted.metadataCsv],'metadata.csv',{type:'text/csv'});
const msMatrix=new File([converted.matrixCsv],'ms_auc.csv',{type:'text/csv'});
const msMapping={subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',omic:'omic',
  condition:'condition',timepoint:'timepoint',batch:'batch',technical_replicate:'technical_replicate',
  sample_type:'sample_type',injection_order:'injection_order'};
const msResult=await runDeterministicAnalysis({
  files:{metadata:msMetadata,metabolomics:msMatrix},
  metadataRows:parseDelimited(converted.metadataCsv).rows,
  columnMapping:msMapping,
  protocol:{organism:'human',objective:'groups',longitudinal:false,designType:'independent',
    studySetting:'animal',groupCount:'2',sampleOverlap:'same_specimen',
    msBlankFilter:'no',msQcRsdFilter:'no',msDriftCorrection:'no'},
  dataTypes:{metabolomics:'peak_area'},
  resolveIdentifiers:false,useReactome:false
});
assert.equal(msResult.engine.analysisMode,'single_omic');
assert.equal(msResult.crossOmics.testedPairs,0);
assert.ok(msResult.layers.metabolomics.rows.length===2);
assert.ok(msResult.layers.metabolomics.rows.every((r)=>r.effectScale==='log2'),
  'Raw MS peak area must be log transformed by the declared type, never interpreted as absolute concentration.');
const msVisuals=await buildMultiomicsVisualizationData({
  files:{metabolomics:msMatrix},metadataRows:parseDelimited(converted.metadataCsv).rows,
  columnMapping:msMapping,dataTypes:{metabolomics:'peak_area'},analysisResult:msResult
});
assert.equal(msVisuals.metabologram.metabolomics.length,2);
assert.equal(msVisuals.metabologram.transcriptomics.length,0);
assert.equal(msVisuals.centralCarbon.measuredMetabolites,2);

const native=await buildMultiomicsVisualizationData({
  files:{}, metadataRows:[], columnMapping:{}, dataTypes:{},
  analysisResult:{layers:{
    metabolomics:{rows:[row('CHEBI:24996',410,0.01,'as_supplied'),row('CHEBI:30031',-70,0.02,'as_supplied')]},
    transcriptomics:{rows:[row('GAPDH',0.15,0.04,'as_supplied')]}
  }}
});
assert.equal(native.metabologram.metabolomics.length,2,
  'Native-scale effects must not silently disappear from the metabologram.');
assert.equal(native.metabologram.metabolomics[0].effectScale,'as_supplied');
assert.equal(native.metabologram.metabolomics[0].effect,410,
  'A native 410-unit effect must never be relabelled or rescaled into log2FC.');
assert.equal(native.centralCarbon.metabolites.find((n)=>n.id==='lactate').measurement.effect,410);
assert.equal(native.focusedMetabolicNetwork.nodes.find((n)=>n.id==='lactate').measurement.effectScale,'as_supplied');


// Optional feature-to-ChEBI table is an explicit user claim, never an inferred
// identification. Statistical results remain exactly as originally fitted.
const sampleAnnotations=parseMetaboliteAnnotations(
  'feature_id,chebi_id\nMS_peak_A,CHEBI:24996\nMS_peak_B,CHEBI:30031\n');
assert.equal(sampleAnnotations.length,2);
assert.equal(sampleAnnotations[0].status,'user_declared');
assert.throws(()=>parseMetaboliteAnnotations('feature_id,chebi_id\nA,CHEBI:unknown'),/CHEBI:number/);
assert.throws(()=>parseMetaboliteAnnotations('feature_id,chebi_id\nA,CHEBI:24996\nA,CHEBI:30031'),/duplicate feature_id/);
assert.throws(()=>parseMetaboliteAnnotations('feature_id,unknown\nA,CHEBI:24996'),/feature_id and chebi_id/);
const anonymousInput={layers:{metabolomics:{rows:[
  row('MS_peak_A',420,0.03,'as_supplied'),row('MS_peak_B',-90,0.03,'as_supplied'),
  row('MS_unidentified',70,0.02,'as_supplied')
]},transcriptomics:{rows:[]}}};
const unannotated=await buildMultiomicsVisualizationData({
  files:{},metadataRows:[],columnMapping:{},dataTypes:{},analysisResult:anonymousInput
});
assert.equal(unannotated.centralCarbon.measuredMetabolites,0);
const annotated=await buildMultiomicsVisualizationData({
  files:{},metadataRows:[],columnMapping:{},dataTypes:{},
  analysisResult:anonymousInput,metaboliteAnnotations:sampleAnnotations
});
assert.equal(annotated.centralCarbon.measuredMetabolites,2);
const userMatch=annotated.focusedMetabolicNetwork.audit.find((row)=>row.input==='MS_peak_A');
assert.equal(userMatch.status,'user_declared');
assert.equal(userMatch.chebi,'CHEBI:24996');
assert.equal(annotated.focusedMetabolicNetwork.audit.find((row)=>row.input==='MS_unidentified').id,null);
assert.equal(annotated.centralCarbon.metabolites.find((node)=>node.id==='lactate').measurement.effect,420);
assert.deepEqual(anonymousInput.layers.metabolomics.rows.map((row)=>row.feature),
  ['MS_peak_A','MS_peak_B','MS_unidentified'],
  'Chemical annotation must never alter fitted feature IDs');

const duplicatedPeaks=await buildMultiomicsVisualizationData({
  files:{},metadataRows:[],columnMapping:{},dataTypes:{},
  analysisResult:{layers:{metabolomics:{rows:[
    row('MS_peak_A',12,0.01,'log2'),row('MS_peak_B',-15,0.02,'log2')
  ]}}},
  metaboliteAnnotations:parseMetaboliteAnnotations(
    'feature_id,chebi_id\nMS_peak_A,CHEBI:24996\nMS_peak_B,CHEBI:24996\n')
});
assert.equal(duplicatedPeaks.centralCarbon.metabolites.find((n)=>n.id==='lactate').measurement.status,
  'multiple_features','Do not pick an arbitrary MS adduct among two features for one chemical identity');
assert.equal(duplicatedPeaks.focusedMetabolicNetwork.nodes.find((n)=>n.id==='lactate').measurement.effect,null);
const unsupportedChemical=await buildMultiomicsVisualizationData({
  files:{},metadataRows:[],columnMapping:{},dataTypes:{},
  analysisResult:{layers:{metabolomics:{rows:[row('MYFEATURE',2,0.01)]}}},
  metaboliteAnnotations:parseMetaboliteAnnotations('feature_id,chebi_id\nMYFEATURE,CHEBI:99999999\n')
});
assert.equal(unsupportedChemical.focusedMetabolicNetwork.audit[0].status,'declared_not_in_network');
assert.equal(unsupportedChemical.focusedMetabolicNetwork.audit[0].id,null);

console.log('multiomics Figure 4-inspired visualizations PASS');
