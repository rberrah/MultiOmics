// True cross-TECHNOLOGY PMx benchmark: Affymetrix MAINZ -> Agilent NKI.
// NKI ER labels are never opened here. Predeclared rank PRIMARY and raw
// SENSITIVITY are always scored and assessed; no test-label model choice.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {createHash} from 'node:crypto';
import {fitFrozenBinaryPredictor,scoreFrozenBinaryPredictor}
  from '../src/lib/multiomics/deterministic.js';

const out=process.argv[2]||'tmp/multiomics/mainz-nki';
const manifest=JSON.parse(await fs.readFile(path.join(out,'source-contract.json'),'utf8'));
assert.equal(manifest.raw_sample_id_overlap,0);
assert.equal(manifest.suspicious_profile_pairs,0);
assert.ok(manifest.common_gene_count>=5000);
assert.equal(manifest.primary,'primary_rank');
assert.equal(manifest.sensitivity,'sensitivity_raw');

function shuffle(rows,seed){
  let state=seed>>>0;
  const random=()=>((state=(Math.imul(state,1664525)+1013904223)>>>0)+.5)/4294967296;
  const labels=rows.map(r=>r.outcome);
  for(let j=labels.length-1;j>0;j--){
    const k=Math.floor(random()*(j+1));
    [labels[k],labels[j]]=[labels[j],labels[k]];
  }
  return rows.map((r,i)=>({...r,outcome:labels[i]}));
}

async function evaluateArm(arm,permCount){
  const dir=path.join(out,arm);
  const trainText=await fs.readFile(path.join(dir,'train.json'),'utf8');
  const train=JSON.parse(trainText);
  const evals=JSON.parse(await fs.readFile(path.join(dir,'test-unlabelled.json'),'utf8'));
  assert.equal(train.features.length,manifest.common_gene_count);
  assert.deepEqual(train.features,evals.features);
  assert.equal(new Set(train.features).size,train.features.length);
  assert.equal(train.rows.length,manifest.mainz_n);
  assert.equal(evals.rows.length,manifest.nki_n);
  assert.ok(train.rows.every(r=>['ER_positive','ER_negative'].includes(r.outcome)));
  assert.ok(evals.rows.every(r=>!Object.hasOwn(r,'outcome')),
    'NKI outcomes reached supervised native PMx training/scoring');
  assert.ok(train.rows.every(r=>r.values.length===train.features.length));
  assert.ok(evals.rows.every(r=>r.values.length===train.features.length));

  const fit=rows=>fitFrozenBinaryPredictor({
    features:train.features,rows,
    negativeLabel:'ER_negative',positiveLabel:'ER_positive',maxFeatures:12
  });
  const frozen=fit(train.rows);
  assert.equal(frozen.fullTrainingSampleCount,manifest.mainz_n);
  assert.equal(frozen.requiredFeatures.length,12);
  assert.ok(frozen.coefficients.every(Number.isFinite));
  const predictions=scoreFrozenBinaryPredictor(frozen,evals);
  assert.equal(predictions.length,manifest.nki_n);
  assert.ok(predictions.every(r=>Number.isFinite(r.prediction)&&
    r.prediction>=0&&r.prediction<=1));
  const enriched={
    ...frozen,
    scientificScope:'Cross-platform research-only independent study benchmark',
    trainingStudy:'MAINZ Affymetrix one-channel',
    evaluationStudy:'NKI Agilent/Rosetta two-channel, heldout',
    transformContract:arm==='primary_rank'?manifest.primary_transform:manifest.sensitivity_transform,
    fullGenePanel:manifest.common_gene_count,
    trainingJsonSha256:createHash('sha256').update(trainText).digest('hex'),
    externalLabelsAccessed:false,clinicalValidation:false
  };
  await fs.writeFile(path.join(dir,'frozen-pmx-model.json'),
    JSON.stringify(enriched,null,2)+'\n');
  const replay=JSON.parse(await fs.readFile(path.join(dir,'frozen-pmx-model.json'),'utf8'));
  assert.deepEqual(scoreFrozenBinaryPredictor(replay,evals),predictions,
    'Serialized native PMx model must replay NKI probabilities exactly');
  assert.throws(()=>scoreFrozenBinaryPredictor(replay,{
    features:evals.features,rows:[{...evals.rows[0],id:train.rows[0].id}]
  }),/Training subject leaked/);
  assert.throws(()=>scoreFrozenBinaryPredictor(replay,{
    features:evals.features.slice(1),rows:evals.rows
  }),/Missing frozen predictor feature|does not match/);
  await fs.writeFile(path.join(dir,'external-predictions.csv'),
    ['subject_id,prediction',...predictions.map(x=>x.subject_id+','+x.prediction)].join('\n')+'\n');

  const nullPredictions=[];
  for(let i=0;i<permCount;i++){
    const permTrain=shuffle(train.rows,20261009+i*31013);
    const permFit=fit(permTrain);
    const scored=scoreFrozenBinaryPredictor(permFit,evals);
    assert.deepEqual(scored.map(x=>x.subject_id),predictions.map(x=>x.subject_id));
    nullPredictions.push(scored.map(x=>x.prediction));
    console.log('MAINZ->NKI '+arm+' native full-refit label null '+(i+1)+'/'+permCount+' PASS');
  }
  await fs.writeFile(path.join(dir,'native-train-null-predictions.csv'),
    ['subject_id,'+nullPredictions.map((_,i)=>'perm_'+i).join(','),
     ...predictions.map((r,i)=>[r.subject_id,...nullPredictions.map(x=>x[i])].join(','))].join('\n')+'\n');
  console.log('MAINZ -> NKI PMx FROZEN '+arm+' PASS: MAINZ='+train.rows.length+
    ' NKI='+predictions.length+' genes='+train.features.length+
    ' selected='+frozen.requiredFeatures.length+' lambda='+frozen.lambda+
    ' nullRefits='+permCount);
  return {arm,trainingCount:train.rows.length,scoredCount:predictions.length,
    genes:train.features.length,lambda:frozen.lambda};
}

const rank=await evaluateArm('primary_rank',10);
const raw=await evaluateArm('sensitivity_raw',10);
await fs.writeFile(path.join(out,'frozen-models-summary.json'),
  JSON.stringify({rank,raw,comparisonProtocol:'Primary was frozen as within-patient rank BEFORE scoring, raw-scale control always reported; no NKI labels accessed.'},null,2)+'\n');
