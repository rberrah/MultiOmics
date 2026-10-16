import assert from 'node:assert/strict';
import {
  MOFA2_ORIENTATION_CONTRACT,
  mofa2BackendCompatibility,
  validateMofa2MethodResult
} from '../src/lib/multiomics/reference-backend-integrity.js';

const ctx={objective:'explore',omicCount:2};
const old={status:'ok',version:'1.3.0',capabilities:{reference_analysis:true}};
assert.equal(mofa2BackendCompatibility(old,ctx).compatible,false);
assert.equal(mofa2BackendCompatibility(null,ctx).compatible,false);
assert.equal(mofa2BackendCompatibility(old,{objective:'groups',omicCount:2}).compatible,true);
assert.equal(mofa2BackendCompatibility(old,{objective:'explore',omicCount:1}).compatible,true);
const newHealth={capabilities:{mofa2_orientation_contract:MOFA2_ORIENTATION_CONTRACT}};
assert.equal(mofa2BackendCompatibility(newHealth,ctx).compatible,true);

const oldResult={status:'ok',methods:{
 mofa2:{method:'MOFA2',status:'ok',summary:{inputAudit:{nSharedSubjects:30}}},
 differential:{status:'ok',method:'DESeq2',summary:{featureCount:100}}
}};
const blocked=validateMofa2MethodResult(oldResult);
assert.equal(blocked.methods.mofa2.status,'blocked');
assert.equal(blocked.methods.mofa2.summary,undefined);
assert.equal(blocked.methods.differential.status,'ok');
assert.equal(oldResult.methods.mofa2.status,'ok','validation must not mutate original data');
assert.equal(blocked.mofa2OrientationBlocked,true);
for(const broken of [
 {implementationContract:MOFA2_ORIENTATION_CONTRACT,matrixOrientation:'feature_rows_subject_columns_for_MOFA2',subjectIdentityPreserved:false},
 {implementationContract:MOFA2_ORIENTATION_CONTRACT,matrixOrientation:'subject_rows_feature_columns',subjectIdentityPreserved:true},
 {implementationContract:'outdated',matrixOrientation:'feature_rows_subject_columns_for_MOFA2',subjectIdentityPreserved:true}
]){
  assert.equal(validateMofa2MethodResult({methods:{mofa2:{status:'ok',summary:broken}}}).methods.mofa2.status,'blocked');
}
const valid={status:'ok',methods:{mofa2:{status:'ok',summary:{
 implementationContract:MOFA2_ORIENTATION_CONTRACT,
 matrixOrientation:'feature_rows_subject_columns_for_MOFA2',
 subjectIdentityPreserved:true,inputAudit:{nSharedSubjects:48}
}}}};
assert.equal(validateMofa2MethodResult(valid),valid);
assert.equal(validateMofa2MethodResult({status:'ok',methods:{mofa2:{status:'blocked'}}}).methods.mofa2.status,'blocked');
console.log('MOFA2 old backend prevention, orientation-contract validation, unchanged other methods: PASS');
