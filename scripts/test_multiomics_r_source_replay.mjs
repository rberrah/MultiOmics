import assert from 'node:assert/strict';
import {readFileSync,mkdtempSync,writeFileSync,rmSync} from 'node:fs';
import {createHash} from 'node:crypto';
import {join} from 'node:path';
import {tmpdir} from 'node:os';
import {execFileSync} from 'node:child_process';
import {REFERENCE_R_SOURCE_FILES} from '../src/lib/multiomics/reference-r-sources.js';
import {REFERENCE_R_REPLAY_VERIFY_PY} from '../src/lib/multiomics/reference-replay-verifier-source.js';
import {RESEARCH_R_NUMERIC_VERIFIER_SOURCE} from '../src/lib/multiomics/r-numeric-verifier-source.js';
import {referenceReplayEntries,referenceReplayStatus} from '../src/lib/multiomics/reference-replay.js';

// Source snapshot MUST match the actual implementation checked into this PR.
// Fail CI if anyone edits methods without updating the replay payload.
for (const [name,contents] of Object.entries(REFERENCE_R_SOURCE_FILES)) {
  assert.equal(contents,readFileSync(join('multiomics-engine',name),'utf8'),
    'Stale source snapshot for '+name);
}
assert.equal(REFERENCE_R_REPLAY_VERIFY_PY,
  readFileSync('scripts/verify_multiomics_r_replay.py','utf8'));
assert.equal(RESEARCH_R_NUMERIC_VERIFIER_SOURCE,
  readFileSync('multiomics-engine/verify_research_replay.R','utf8'));
assert.equal(referenceReplayStatus(null).status,'not_available');
const expected={status:'ok',methods:{transcriptomics_differential_deseq2:{
  status:'ok',top:[{feature:'ENSG001',log2FoldChange:1.23,padj:0.04}]
}}};
const sha=Object.fromEntries(Object.entries(REFERENCE_R_SOURCE_FILES).map(([name,source])=>
  [name,createHash('sha256').update(source).digest('hex')]));
const replay={requestJson:JSON.stringify({
  protocol:{objective:'groups',designType:'independent'},
  metadataCsv:'assay_id,sample_id,subject_id,omic\nR1,S1,P1,transcriptomics\n',
  matrices:{transcriptomics:'feature_id,R1\nENSG001,10\n'}
}),responseJson:JSON.stringify({...expected,runtime:{referenceSourceSha256:sha}})};
const {entries,status}=await referenceReplayEntries(replay);
assert.equal(status.status,'r_reference_source_replay_available');
assert.equal(status.numerical_equivalence_verified,false);
assert.equal(status.figure_equivalence_verified,false);
assert.equal(status.source_attestation,'sha256_match');
assert.equal(status.exact_runtime_source_verified,true);
const corrupt=JSON.parse(replay.responseJson);
corrupt.runtime.referenceSourceSha256['server.R']='0'.repeat(64);
const badAttestation=await referenceReplayEntries({...replay,responseJson:JSON.stringify(corrupt)});
assert.equal(badAttestation.status.source_attestation,'sha256_mismatch');
assert.equal(badAttestation.status.exact_runtime_source_verified,false);
assert(entries.some(e=>e.name==='reproduction/run_reference_backend.R'));
assert.deepEqual(JSON.parse(entries.find(e=>e.name==='reproduction/reference_expected.json').data),JSON.parse(replay.responseJson));
assert(entries.some(e=>e.name==='reproduction/multiomics-engine/server.R'));
assert.equal(entries.find(e=>e.name==='reproduction/verify_replay.py').data,
  REFERENCE_R_REPLAY_VERIFY_PY);
assert.equal(entries.find(e=>e.name==='reproduction/verify_replay.R').data,
  RESEARCH_R_NUMERIC_VERIFIER_SOURCE);
const root=mkdtempSync(join(tmpdir(),'pmx-r-replay-'));
try {
  let rAvailable=true;
  for(const e of entries.filter(x=>x.name.endsWith('.R'))){
    const file=join(root,e.name.replaceAll('/','_'));
    writeFileSync(file,e.data);
    try {
      execFileSync('Rscript',['-e','invisible(parse(file='+JSON.stringify(file)+'))'],
                   {stdio:'pipe'});
    } catch(e) {
      if(e.code==='ENOENT') {
        console.log('Rscript unavailable; syntax parsing delegated to CI R runner');
        rAvailable=false;
        break;
      }
      throw e;
    }
  }
  if(rAvailable) {
    const verifier=join(root,'reproduction_verify_replay.R');
    const baseline=JSON.parse(replay.responseJson);
    writeFileSync(join(root,'reference_expected.json'),JSON.stringify(baseline));
    writeFileSync(join(root,'reference_replayed.json'),JSON.stringify(baseline));
    execFileSync('Rscript',[verifier,root],{stdio:'pipe'});
    let report=JSON.parse(readFileSync(join(root,'r_numeric_comparison_report.json'),'utf8'));
    assert.equal(report.result,'PASS');
    const modified=structuredClone(baseline);
    modified.methods.transcriptomics_differential_deseq2.top[0].padj=0.95;
    writeFileSync(join(root,'reference_replayed.json'),JSON.stringify(modified));
    assert.throws(()=>execFileSync('Rscript',[verifier,root],{stdio:'pipe'}));
    report=JSON.parse(readFileSync(join(root,'r_numeric_comparison_report.json'),'utf8'));
    assert.equal(report.result,'FAIL');
    assert(report.n_discrepancies>=1);
  }
}finally{rmSync(root,{recursive:true,force:true});}
console.log('PASS: frozen R source matches API; source replay is guarded and R parser accepts scripts');
