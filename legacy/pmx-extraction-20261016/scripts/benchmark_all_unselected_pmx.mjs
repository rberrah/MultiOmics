// NATIVE PMx high-dimensional single-omic holdout, 12,625 unselected probes.
// The source's outcome CSV is deliberately NEVER opened in this process.
// Posthoc metric assessment is delegated to an independently invoked R script.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {fitFrozenBinaryPredictor,scoreFrozenBinaryPredictor}
  from '../src/lib/multiomics/deterministic.js';

const dest=process.argv[2]||'tmp/multiomics/all-native';
const [trainText,testText,sourceText]=await Promise.all([
  fs.readFile(path.join(dest,'train.json'),'utf8'),
  fs.readFile(path.join(dest,'test-unlabelled.json'),'utf8'),
  fs.readFile(path.join(dest,'source-contract.json'),'utf8')
]);
const training=JSON.parse(trainText),test=JSON.parse(testText),source=JSON.parse(sourceText);
assert.equal(training.features.length,12625);
assert.deepEqual(training.features,test.features);
assert.equal(new Set(training.features).size,12625);
assert.equal(training.rows.length,source.train_n);
assert.equal(test.rows.length,source.heldout_n);
assert.equal(training.rows.length+test.rows.length,79);
assert.ok(training.rows.every(r=>['BCR/ABL','NEG'].includes(r.outcome)));
assert.ok(test.rows.every(r=>!Object.hasOwn(r,'outcome')));
assert.ok(training.rows.every(r=>r.values.length===12625));
assert.ok(test.rows.every(r=>r.values.length===12625));

const fit=(rows)=>fitFrozenBinaryPredictor({
  features:training.features,rows,
  positiveLabel:'BCR/ABL',negativeLabel:'NEG',maxFeatures:12
});
const model=fit(training.rows);
assert.equal(model.requiredFeatures.length,12);
assert.ok(Number.isFinite(model.lambda));
assert.equal(model.fullTrainingSampleCount,source.train_n);
const probabilities=scoreFrozenBinaryPredictor(model,test);
assert.equal(probabilities.length,source.heldout_n);
assert.ok(probabilities.every(p=>Number.isFinite(p.prediction) &&
  p.prediction>=0 && p.prediction<=1));
const savedModel={
  ...model,trainingInputSha256:createHash('sha256').update(trainText).digest('hex'),
  sourceEvidence:'Bioconductor ALL 12,625 original probes with internal stratified split',
  outcomeLeakageCheck:'heldout labels physically unavailable to this Node stage',
  studyCertification:false
};
await fs.writeFile(path.join(dest,'frozen-pmx-model.json'),
  JSON.stringify(savedModel,null,2)+'\n');
const replay=JSON.parse(await fs.readFile(path.join(dest,'frozen-pmx-model.json'),'utf8'));
assert.deepEqual(scoreFrozenBinaryPredictor(replay,test),probabilities,
  'Model reload must replay every heldout probability exactly.');
assert.throws(()=>scoreFrozenBinaryPredictor(replay,{
  features:test.features,rows:[{...test.rows[0],id:training.rows[0].id}]
}),/Training subject leaked/);
assert.throws(()=>scoreFrozenBinaryPredictor(replay,{
  features:test.features.slice(1),rows:test.rows
}),/Missing frozen predictor feature|does not match/);
await fs.writeFile(path.join(dest,'frozen-pmx-predictions.csv'),
  ['subject_id,prediction',...probabilities.map(p=>p.subject_id+','+p.prediction)].join('\n')+'\n');

// Ten full train-label refits: selection, scaling, penalty tuning and fit
// are recomputed each time WITHOUT opening test outcomes.
const originalLabels=training.rows.map(r=>r.outcome);
function shuffled(seed){
  let x=seed>>>0;
  const next=()=>((x=(Math.imul(x,1664525)+1013904223)>>>0)+.5)/4294967296;
  const a=originalLabels.slice();
  for(let i=a.length-1;i>0;i--){
    const j=Math.floor(next()*(i+1));
    [a[i],a[j]]=[a[j],a[i]];
  }
  return a;
}
const nullOutputs=[];
for(let b=0;b<10;b++){
  const labels=shuffled(20261009+9973*b);
  const changed=training.rows.map((r,i)=>({...r,outcome:labels[i]}));
  const nullModel=fit(changed);
  const frozen=scoreFrozenBinaryPredictor(nullModel,test);
  assert.deepEqual(frozen.map(p=>p.subject_id),probabilities.map(p=>p.subject_id));
  nullOutputs.push(frozen.map(p=>p.prediction));
  console.log('ALL native training-label null refit '+(b+1)+'/10 PASS');
}
const cols=nullOutputs.map((_,j)=>'perm_'+j);
await fs.writeFile(path.join(dest,'pmx-train-label-null-predictions.csv'),
  ['subject_id,'+cols.join(','),
    ...probabilities.map((r,i)=>[r.subject_id,...nullOutputs.map(a=>a[i])].join(','))].join('\n')+'\n');
console.log('ALL PMx NATIVE FROZEN PREDICTIONS PASS: 12625 measured probes, train='+training.rows.length+
  ', heldout='+probabilities.length+', features frozen='+model.requiredFeatures.length+
  ', lambda='+model.lambda+', 10 negative refits.');
