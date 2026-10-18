import assert from 'node:assert/strict';
import { evaluateReferenceBackendHealth as evaluate } from '../src/lib/multiomics/reference-backend-preflight.js';

const base = {
  status:'ok',engine:'PMx Explain reference R backend',version:'1.4.0',
  capabilities:{reference_analysis:true,
    mofa2_orientation_contract:'mofa2-feature-rows-sample-columns-v2'},
  runtime:{RVersion:'4.5.2',RPlatform:'x86_64-pc-linux-gnu'},
  packages:{DESeq2:true,edgeR:true,limma:true,lmerTest:true,
    survival:true,MOFA2:true,mixOmics:true,fgsea:true}
};
const clone = () => structuredClone(base);
const confirm = {intent:'confirmatory',objective:'groups',designType:'independent',
  hasTranscriptomics:true,transcriptomicsValues:'raw_counts'};
const reject = (health,context,code) => {
  const result = evaluate(health,context);
  assert.equal(result.ok,false);
  assert.equal(result.code,code);
  assert.ok(result.message.length && result.messageEn.length);
};

assert.equal(evaluate(clone(),confirm).code,'verified_method_available');
assert.equal(evaluate(clone(),{}).code,'verified_exploratory');
reject(null,confirm,'invalid_health_contract');
reject({status:'ok'},confirm,'invalid_health_contract');
reject({...clone(),engine:'unrelated local service'},confirm,'invalid_health_contract');
reject({...clone(),version:'unknown'},confirm,'invalid_health_contract');
reject({...clone(),runtime:{}},confirm,'invalid_health_contract');
reject({...clone(),capabilities:{reference_analysis:'true'}},confirm,'invalid_health_contract');
const incomplete=clone(); delete incomplete.packages.limma;
reject(incomplete,confirm,'incomplete_package_manifest');
const none=clone();
for(const key of Object.keys(none.packages))none.packages[key]=false;
reject(none,confirm,'no_scientific_methods');
const unavailable=clone();
unavailable.packages.DESeq2=false;
unavailable.packages.edgeR=false;
reject(unavailable,confirm,'missing_rnaseq_methods');
const voom=clone();
voom.packages.DESeq2=false;
assert.equal(evaluate(voom,confirm).ok,true);
const deseq=clone();
deseq.packages.edgeR=false;
assert.equal(evaluate(deseq,confirm).ok,true);
const processed=clone();
processed.packages.limma=false;
reject(processed,{...confirm,transcriptomicsValues:'processed'},'missing_limma');
const longitudinal=clone();
longitudinal.packages.lmerTest=false;
reject(longitudinal,{intent:'confirmatory',objective:'time'},'missing_lmerTest');
const outcome=clone();
outcome.packages.survival=false;
outcome.packages.mixOmics=false;
reject(outcome,{intent:'confirmatory',objective:'outcome'},'missing_outcome_methods');
reject(clone(),{intent:'confirmatory',objective:'explore'},'no_confirmatory_model');

// A positive health check alone never reports validated scientific inference.
const success=evaluate(clone(),confirm);
assert.equal(success.ok,true);
assert.match(success.messageEn,/assumptions/);
assert.equal('scientificallyValidated' in success,false);
console.log('PASS: reference R preflight identity, package gates, study context and negative controls');
