// Score the independent TRANSBIG calibration cohort with PMx frozen MAINZ
// native ridge models. This process MUST NOT open TRANSBIG ER labels or NKI
// outcomes, and cannot refit any feature, penalty or coefficient.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import {scoreFrozenBinaryPredictor} from '../src/lib/multiomics/deterministic.js';

const out=process.argv[2]||'tmp/multiomics/mainz-nki';
const provenance=JSON.parse(await fs.readFile(path.join(out,'transbig-calibration-contract.json'),'utf8'));
assert.ok(provenance.mapped_gene_panel>=5000);
assert.equal(provenance.labels_blinded_to_js,true);
const results={};
for(const arm of ['primary_rank','sensitivity_raw']){
  const loc=path.join(out,arm);
  const model=JSON.parse(await fs.readFile(path.join(loc,'frozen-pmx-model.json'),'utf8'));
  const train=JSON.parse(await fs.readFile(path.join(loc,'train.json'),'utf8'));
  const cal=JSON.parse(await fs.readFile(path.join(loc,'calibration-unlabelled.json'),'utf8'));
  const heldout=await fs.readFile(path.join(loc,'external-predictions.csv'),'utf8');
  assert.equal(model.fullGenePanel,provenance.mapped_gene_panel);
  assert.deepEqual(cal.features,train.features);
  assert.equal(cal.rows.length,provenance.eligible_calibration_participants);
  assert.ok(cal.rows.every(r=>!Object.hasOwn(r,'outcome')));
  assert.ok(cal.rows.every(r=>r.values.length===provenance.mapped_gene_panel));
  const preds=scoreFrozenBinaryPredictor(model,cal);
  assert.equal(preds.length,cal.rows.length);
  assert.ok(preds.every(p=>Number.isFinite(p.prediction) && p.prediction>=0 && p.prediction<=1));
  assert.ok(preds.every(p=>!train.rows.some(t=>t.id===p.subject_id)));
  const csv=['subject_id,prediction',
    ...preds.map(p=>p.subject_id+','+p.prediction)].join('\n')+'\n';
  await fs.writeFile(path.join(loc,'transbig-calibration-frozen-scores.csv'),csv);
  assert.equal(await fs.readFile(path.join(loc,'external-predictions.csv'),'utf8'),heldout,
    'Scoring calibration subjects must NOT mutate NKI frozen probabilities.');
  results[arm]={n:preds.length,modelSha256:model.trainingJsonSha256};
  console.log('PMx MAINZ -> TRANSBIG BLIND CALIBRATION SCORES '+arm+' PASS: n='+preds.length);
}
await fs.writeFile(path.join(out,'transbig-scoring-contract.json'),
 JSON.stringify({arms:results,scoring:'Frozen native PMx models used unchanged; no NKI or TRANSBIG labels read by Node.'},null,2)+'\n');
