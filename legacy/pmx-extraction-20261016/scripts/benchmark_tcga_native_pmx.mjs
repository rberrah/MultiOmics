// Real PMx frozen ridge classifier evaluated on publisher-heldout TCGA.
// No reference-R training model: calls the same core ridgeGlmFit used by
// PMx internal cross-validation. The official holdout label file is NOT read.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {fitFrozenBinaryPredictor,scoreFrozenBinaryPredictor} from
  '../src/lib/multiomics/deterministic.js';

const dir=process.argv[2]||'tmp/multiomics/tcga-native';
const [trainText,testText]=await Promise.all([
  fs.readFile(path.join(dir,'train.json'),'utf8'),
  fs.readFile(path.join(dir,'test-unlabelled.json'),'utf8')
]);
const train=JSON.parse(trainText),test=JSON.parse(testText);
assert.equal(train.rows.length,105);
assert.equal(test.rows.length,70);
assert.deepEqual(train.features,test.features);
assert.equal(new Set(train.features).size,train.features.length);
assert.ok(train.features.length>=50);
assert.ok(train.rows.every(x=>['Her2','LumA'].includes(x.outcome)));
assert.ok(test.rows.every(x=>!Object.hasOwn(x,'outcome')),'Holdout outcomes reached prediction code');
const model=fitFrozenBinaryPredictor({
  features:train.features,rows:train.rows,
  positiveLabel:'Her2',negativeLabel:'LumA',maxFeatures:12
});
assert.equal(model.schema,'pmx_frozen_binary_ridge_v1');
assert.equal(model.fullTrainingSampleCount,105);
assert.equal(model.requiredFeatures.length,12);
assert.ok(model.coefficients.every(Number.isFinite));
assert.ok(model.requiredFeatures.every(x=>train.features.includes(x.feature)));
const predictions=scoreFrozenBinaryPredictor(model,test);
assert.equal(predictions.length,70);
assert.ok(predictions.every(x=>Number.isFinite(x.prediction) &&
  x.prediction>=0 && x.prediction<=1));
const unseenIds=new Set(test.rows.map(r=>r.id));
assert.equal(unseenIds.size,70);
assert.ok(predictions.every(x=>unseenIds.has(x.subject_id)));
assert.ok(predictions.every(x=>!model.trainingSubjectIds.includes(x.subject_id)));
// A renamed test feature or a repeated training subject must NOT silently
// receive a score, even if it would be numerically possible.
const renamedFeatures=test.features.slice();
renamedFeatures[train.features.indexOf(model.requiredFeatures[0].feature)]='UNKNOWN_FEATURE';
assert.throws(()=>scoreFrozenBinaryPredictor(model,{
  features:renamedFeatures,rows:test.rows
}),/Missing frozen predictor feature/);
assert.throws(()=>scoreFrozenBinaryPredictor(model,{
  features:test.features,rows:[{...test.rows[0],id:train.rows[0].id}]
}),/Training subject leaked/);
const alteredTest=test.rows.map((row,i)=>i===0?{
  ...row,values:row.values.slice().map((x,j)=>
    j===train.features.indexOf(model.requiredFeatures[0].feature)&&x!==null ?x+8:x)
}:row);
const other=scoreFrozenBinaryPredictor(model,{features:test.features,rows:alteredTest});
assert.ok(other.slice(1).every((item,i)=>item.prediction===predictions[i+1].prediction),
  'Changing one heldout subject must not alter another subject prediction.');
assert.deepEqual(model,fitFrozenBinaryPredictor({
  features:train.features,rows:train.rows.slice().reverse(),
  positiveLabel:'Her2',negativeLabel:'LumA',maxFeatures:12
}), 'Frozen model must not depend on import row ordering.');
const digest=createHash('sha256').update(trainText).digest('hex');
const modelReport={
  ...model,trainingInputSha256:digest,
  intendedUse:'Binary research benchmarking, never clinical deployment',
  clinicalValidation:false,
  dataProvenance:'mixOmics source-native public filtered/normalized example',
  fittingProvenance:'100% training-only including feature selection, imputation, scaling, lambda CV',
  heldoutLabelsAccessed:false,
  caution:'Heldout test is from same TCGA source; preprocessing and preselection of source upstream are outside PMx control.'
};
await fs.writeFile(path.join(dir,'frozen-pmx-model.json'),
  JSON.stringify(modelReport,null,2)+'\n');
// Prove the persisted model (not only the in-memory fit) can be replayed
// verbatim after a fresh JSON parse, without labels or feature retraining.
const replayed=JSON.parse(await fs.readFile(path.join(dir,'frozen-pmx-model.json'),'utf8'));
assert.deepEqual(scoreFrozenBinaryPredictor(replayed,test),predictions,
  'Reloaded frozen PMx coefficients must exactly replay all 70 predictions.');
await fs.writeFile(path.join(dir,'frozen-pmx-predictions.csv'),
  ['subject_id,prediction',...predictions.map(x=>x.subject_id+','+x.prediction)].join('\n')+'\n');

// Test the NATIVE PMx training process under destroyed supervision.
// Each permutation reruns PMx's own feature selection, stratified penalty
// tuning, ridge fit and frozen scoring. Holdout labels are still unavailable.
const baseLabels=train.rows.map(x=>x.outcome);
const shuffled=(seed)=>{
  let state=seed>>>0;
  const random=()=>((state=(Math.imul(state,1664525)+1013904223)>>>0)+.5)/4294967296;
  const result=baseLabels.slice();
  for(let j=result.length-1;j>0;j--){
    const k=Math.floor(random()*(j+1));
    [result[k],result[j]]=[result[j],result[k]];
  }
  return result;
};
const negativePredictions=[];
for(let perm=0;perm<30;perm++){
  const labels=shuffled(20261011+perm*32771);
  const permTrain=train.rows.map((r,i)=>({...r,outcome:labels[i]}));
  const nullModel=fitFrozenBinaryPredictor({
    features:train.features,rows:permTrain,
    positiveLabel:'Her2',negativeLabel:'LumA',maxFeatures:12
  });
  const nullScores=scoreFrozenBinaryPredictor(nullModel,test);
  assert.deepEqual(nullScores.map(p=>p.subject_id),predictions.map(p=>p.subject_id));
  negativePredictions.push(nullScores.map(p=>p.prediction));
}
await fs.writeFile(path.join(dir,'pmx-train-label-null-predictions.csv'),
  ['subject_id,'+negativePredictions.map((_,i)=>'perm_'+i).join(','),
  ...predictions.map((p,i)=>[
    p.subject_id,...negativePredictions.map(scores=>scores[i])
  ].join(','))].join('\n')+'\n');

console.log('PMx frozen native benchmark PREDICTIONS PASS: train=105 scored_holdout=70 selected='+model.requiredFeatures.length+
  ' training-only lambda='+model.lambda);
