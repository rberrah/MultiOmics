import assert from 'node:assert/strict';
import {projectNativeRNumbers} from '../src/lib/multiomics/r-primary-results.js';

const originalLayer=()=>({
  effectScale:'log2',
  qc:{preprocessingAudit:{declaredValueType:'log_intensity'}},
  rows:[{feature:'browser-only',effect:99,pValue:0,qValue:0}],
  steps:['browser exploratory fit']
});
const result=(protocol, methods)=>({
  protocol,
  layers:{proteomics:originalLayer()},
  referenceBackend:{status:'ok',methods}
});

const longCsv=[
  'feature,term,estimate,std_error,df,statistic,p_value,q_bh,n_observations,n_subjects,random_effect_structure,status',
  'P1,conditionB:time,0.25,0.08,31,3.125,0.004,0.016,48,24,random_intercept,ok',
  'P2,conditionB:time,-0.10,0.05,30,-2,0.05,0.10,48,24,random_intercept_and_slope,ok'
].join('\n');
const longMethod={status:'ok',method:'R lmerTest group-by-time slope difference',
  referenceTableComplete:true,referenceTableRows:2,referenceTableCsv:longCsv};
const projectedLong=projectNativeRNumbers(result(
  {analysisIntent:'confirmatory',objective:'time',outcomeType:'none'},
  {proteomics_longitudinal:longMethod}
));
assert.equal(projectedLong.layers.proteomics.mode,'reference-r');
assert.equal(projectedLong.layers.proteomics.primaryReferenceMethodKey,'proteomics_longitudinal');
assert.equal(projectedLong.layers.proteomics.effectScale,'group_by_time_slope');
assert.equal(projectedLong.layers.proteomics.rows[0].effect,0.25);
assert.equal(projectedLong.layers.proteomics.rows[0].qValue,0.016);
assert.equal(projectedLong.layers.proteomics.rows[0].term,'conditionB:time');
assert.equal(projectedLong.layers.proteomics.rows[0].nSubjects,24);
assert.equal(projectedLong.browserExploratoryLayers.proteomics.rows[0].feature,'browser-only');

const survCsv=[
  'feature,effect,standardError,pValue,qValue,hazardRatio,ciLow,ciHigh,proportionalHazardsP,globalPhP,nSubjects,nEvents,missingFraction,status',
  'M1,0.693147,0.2,0.0005,0.001,2.0,0.30,1.08,0.45,0.40,80,35,0,ph_test_not_rejected',
  'M2,-0.223144,0.1,0.025,0.05,0.8,-0.42,-0.02,0.20,0.15,80,35,0,ph_test_not_rejected'
].join('\n');
const survMethod={status:'ok',method:'R survival::coxph Efron + cox.zph',
  referenceTableComplete:true,referenceTableRows:2,referenceTableCsv:survCsv};
const projectedSurv=projectNativeRNumbers(result(
  {analysisIntent:'confirmatory',objective:'outcome',outcomeType:'survival'},
  {proteomics_survival:survMethod}
));
assert.equal(projectedSurv.layers.proteomics.primaryReferenceMethodKey,'proteomics_survival');
assert.equal(projectedSurv.layers.proteomics.effectScale,'log_hazard');
assert.equal(projectedSurv.layers.proteomics.rows[0].hazardRatio,2);
assert.equal(projectedSurv.layers.proteomics.rows[0].nEvents,35);
assert.equal(projectedSurv.layers.proteomics.rows[0].proportionalHazardsP,0.45);

const unsupported=result(
  {analysisIntent:'confirmatory',objective:'outcome',outcomeType:'binary'},
  {proteomics_survival:survMethod}
);
const unchanged=projectNativeRNumbers(unsupported);
assert.equal(unchanged.layers.proteomics.rows[0].feature,'browser-only');
assert.equal(unchanged.displayedReferenceMethodCount,undefined);

const incomplete=result(
  {analysisIntent:'confirmatory',objective:'time'},
  {proteomics_longitudinal:{...longMethod,referenceTableRows:3}}
);
assert.throws(()=>projectNativeRNumbers(incomplete),/Incomplete R reference CSV/);

console.log('PASS: native R longitudinal and survival tables replace browser numbers only for supported objectives');
