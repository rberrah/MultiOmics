// Scientific regression: neither an independent held-out assay nor its label
// may influence training-sample LFQ/MS transformation before nested CV.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { prepareFoldIsolatedPredictionMatrix } from '../src/lib/multiomics/deterministic-core.js';
import { runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';

function simpleMatrix(otherAssayMultiplier = 1) {
  const assays = ['S01', 'S02', 'S03'];
  const features = Array.from({length: 16}, (_,i)=>'F' + i);
  const values = new Map(features.map((feature,i) => [
    feature,
    new Map(assays.map((assay,j)=>[
      assay, Math.round((1100+i*90+j*45+(i%3)*13) *
        (j===2?otherAssayMultiplier:1))
    ]))
  ]));
  return {assays,features,values,transposed:false};
}
function get(matrix, feature, assay) {
  return matrix.values.get(feature).get(assay);
}

for (const [layer, valueType] of [
  ['proteomics','lfq_intensity'],
  ['metabolomics','peak_area']
]) {
  const source = simpleMatrix();
  const altered = simpleMatrix(1000000);
  const one = prepareFoldIsolatedPredictionMatrix(source,layer,valueType);
  const two = prepareFoldIsolatedPredictionMatrix(altered,layer,valueType);
  assert.ok(one && two);
  assert.equal(one.features.length,source.features.length);
  for (const assay of ['S01','S02']) for (const feature of source.features) {
    assert.equal(get(one,feature,assay),get(two,feature,assay),
      'Changing only a held-out sample must not alter other samples');
  }
  for (const assay of source.assays) {
    const vals = source.features.map(feature=>get(one,feature,assay)).sort((a,b)=>a-b);
    assert.ok(Math.abs((vals[7]+vals[8])/2)<1e-12,
      'assay-local median should be zero');
  }
  assert.equal(get(source,'F0','S01'),1100,'source data are not modified');
  assert.match(one.predictionPreparation,/no cohort-wide feature filtering/);
  const negative = simpleMatrix();
  negative.values.get('F0').set('S01',-1);
  assert.equal(prepareFoldIsolatedPredictionMatrix(negative,layer,valueType),null,
    'invalid negative raw-intensity input must be refused by isolated route');
}
assert.equal(prepareFoldIsolatedPredictionMatrix(simpleMatrix(),'metabolomics','normalized'),null);
const asLog = prepareFoldIsolatedPredictionMatrix(simpleMatrix(),'proteomics','log_intensity');
assert.equal(get(asLog,'F0','S01'),1100,'declared log values are not secretly retransformed');

const subjects = Array.from({length:30},(_,i)=>'P'+String(i+1).padStart(2,'0'));
const features = Array.from({length:20},(_,i)=>'FEATURE_'+i);
const metadataRows = [];
function createOmic(layer,prefix,magnitude) {
  const assays = subjects.map((subject)=>prefix+'_'+subject);
  const lines = [['feature_id',...assays].join(',')];
  for (let j=0;j<features.length;j++) {
    const vals=subjects.map((subject,i)=>{
      const signal=j<3&&i>=15?1.7:1;
      return String(Math.round(
        magnitude*(1+0.05*j+0.012*(i%7)+0.008*((i*j)%11))*signal
      ));
    });
    lines.push([features[j],...vals].join(','));
  }
  subjects.forEach((subject,i)=>metadataRows.push({
    subject_id:subject,sample_id:subject,assay_id:assays[i],
    omic:layer,condition:'cohort',timepoint:'T0',batch:'',
    outcome:i>=15?'B':'A',sample_type:'biological'
  }));
  return new File([lines.join('\n')],layer+'.csv',{type:'text/csv'});
}
const protein=createOmic('proteomics','PROT',100000);
const metabolite=createOmic('metabolomics','MS',300000);
const result=await runDeterministicAnalysis({
  files:{
    metadata:new File(['fixture'],'metadata.csv',{type:'text/csv'}),
    transcriptomics:null,proteomics:protein,metabolomics:metabolite
  },
  metadataRows,
  columnMapping:{
    subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',
    omic:'omic',condition:'condition',timepoint:'timepoint',batch:'batch',
    outcome:'outcome',sample_type:'sample_type'
  },
  protocol:{
    organism:'human',objective:'outcome',outcomeType:'binary',
    designType:'independent',longitudinal:false,studySetting:'synthetic_test',
    outcomeTimepoint:'T0',batchKnown:'no',groupCount:'1',covariateColumns:[],
    msMnarStrategy:'none',msDriftCorrection:false,msQcRsdFilter:false,
    msBlankFilter:'off'
  },
  dataTypes:{proteomics:'lfq_intensity',metabolomics:'peak_area'},
  useReactome:false,resolveIdentifiers:false
});
assert.equal(result.predictiveOutcome.status,'ok');
assert.equal(result.predictiveOutcome.preprocessingLeakageRisk,false,
  'raw-intensity omics must not use cohort-fitted predictive preprocessing');
assert.equal(result.predictiveOutcome.rawIntensityScreening,true);
assert.ok(result.predictiveOutcome.technicalQcCaveat.includes('pooled-QC'));
assert.equal(result.predictiveOutcome.predictions.length,subjects.length);
assert.equal(new Set(result.predictiveOutcome.predictions.map(x=>x.subjectId)).size,subjects.length);
assert.ok(Number.isFinite(result.predictiveOutcome.metrics.auc));
const brier=result.predictiveOutcome.metrics.brier;
const p=result.predictiveOutcome.predictions;
const levels=[...new Set(p.map(item=>item.truth))].sort();
assert.equal(levels.length,2);
const independentlyComputedBrier=p.reduce((sum,item)=>
  sum+(item.prediction-(item.truth===levels[1]?1:0))**2,0)/p.length;
assert.ok(Number.isFinite(brier)&&brier>=0&&brier<=1,
  'Brier must be finite and in [0,1] for genuine held-out binary probability predictions');
assert.ok(Math.abs(brier-independentlyComputedBrier)<1e-12,
  'Brier is the exact mean of the independent outer-fold squared probability errors');

console.log('multiomics LFQ/MS assay-local prediction, held-out invariance and nested CV: PASS');
