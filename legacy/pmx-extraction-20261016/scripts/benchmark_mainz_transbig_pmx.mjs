// Real PMx independent-study validation: MAINZ training only;
// TRANSBIG has NO outcomes in this Node process.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {fitFrozenBinaryPredictor,scoreFrozenBinaryPredictor}
  from '../src/lib/multiomics/deterministic.js';

const dest=process.argv[2]||'tmp/multiomics/mainz-transbig';
const [trainText,testText,manifestText]=await Promise.all([
  fs.readFile(path.join(dest,'train.json'),'utf8'),
  fs.readFile(path.join(dest,'test-unlabelled.json'),'utf8'),
  fs.readFile(path.join(dest,'source-contract.json'),'utf8')
]);
const training=JSON.parse(trainText);
const test=JSON.parse(testText);
const manifest=JSON.parse(manifestText);
assert.equal(training.features.length,manifest.all_shared_probes);
assert.equal(training.features.length,22283);
assert.deepEqual(training.features,test.features);
assert.equal(training.rows.length,manifest.train_n);
assert.equal(test.rows.length,manifest.external_n);
assert.ok(training.rows.every(x=>['ER_positive','ER_negative'].includes(x.outcome)));
assert.ok(test.rows.every(x=>!Object.hasOwn(x,'outcome')));
assert.equal(manifest.global_probes_removed,0);
assert.equal(manifest.expression_duplicate_audit.suspected_pairs_above_0995,0);

function train(rows){
  return fitFrozenBinaryPredictor({
    features:training.features,rows,
    positiveLabel:'ER_positive',negativeLabel:'ER_negative',maxFeatures:12
  });
}
const model=train(training.rows);
assert.equal(model.fullTrainingSampleCount,manifest.train_n);
assert.equal(model.requiredFeatures.length,12);
assert.ok(model.requiredFeatures.every(x=>training.features.includes(x.feature)));
const predictions=scoreFrozenBinaryPredictor(model,test);
assert.equal(predictions.length,test.rows.length);
assert.ok(predictions.every(x=>Number.isFinite(x.prediction) && x.prediction>=0 && x.prediction<=1));
const frozen={
  ...model,
  source:'MAINZ trained, TRANSBIG scored without outcome labels',
  sourceProbeCount:manifest.all_shared_probes,
  trainingInputSha256:createHash('sha256').update(trainText).digest('hex'),
  intendedUse:'Research benchmarking only; not for a medical decision',
  clinicalValidation:false,
  sourcePreprocessingCaveat:manifest.preprocessing
};
await fs.writeFile(path.join(dest,'frozen-pmx-model.json'),
  JSON.stringify(frozen,null,2)+'\n');
const reopened=JSON.parse(await fs.readFile(path.join(dest,'frozen-pmx-model.json'),'utf8'));
assert.deepEqual(scoreFrozenBinaryPredictor(reopened,test),predictions,
  'Replayed frozen PMx model must reproduce every TRANSBIG probability.');
assert.throws(()=>scoreFrozenBinaryPredictor(reopened,{
  features:test.features,rows:[{...test.rows[0],id:training.rows[0].id}]
}),/Training subject leaked/);
const predictionsCsv=['subject_id,prediction',...predictions.map(p=>p.subject_id+','+p.prediction)].join('\n')+'\n';
await fs.writeFile(path.join(dest,'transbig-frozen-predictions.csv'),predictionsCsv);

// External labels remain unread. Ten complete null refits invalidate the
// supervision while rerunning ranking, feature scaling, inner CV and ridge fit.
const labels=training.rows.map(r=>r.outcome);
function permute(seed){
  let state=seed>>>0;
  const rng=()=>((state=(Math.imul(state,1664525)+1013904223)>>>0)+.5)/4294967296;
  const shuffled=labels.slice();
  for(let i=shuffled.length-1;i>0;i--){
    const j=Math.floor(rng()*(i+1));
    [shuffled[i],shuffled[j]]=[shuffled[j],shuffled[i]];
  }
  return shuffled;
}
const nullMatrix=[];
for(let i=0;i<10;i++){
  const randomized=permute(20261009+65537*i);
  const nullRows=training.rows.map((r,j)=>({...r,outcome:randomized[j]}));
  const nullModel=train(nullRows);
  const scored=scoreFrozenBinaryPredictor(nullModel,test);
  assert.deepEqual(scored.map(x=>x.subject_id),predictions.map(x=>x.subject_id));
  nullMatrix.push(scored.map(x=>x.prediction));
  console.log('MAINZ->TRANSBIG PMx full training-label null refit '+(i+1)+'/10 OK');
}
await fs.writeFile(path.join(dest,'train-label-null-scores.csv'),
  ['subject_id,'+nullMatrix.map((_,i)=>'perm_'+i).join(','),
  ...predictions.map((p,i)=>[p.subject_id,...nullMatrix.map(column=>column[i])].join(','))].join('\n')+'\n');
console.log('CROSS STUDY PMx FROZEN PREDICTIONS PASS: MAINZ n='+training.rows.length+
  ' TRANSBIG n='+predictions.length+' shared probes='+training.features.length+
  ' selected='+model.requiredFeatures.length+' lambda='+model.lambda);
