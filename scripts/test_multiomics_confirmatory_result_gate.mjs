import assert from 'node:assert/strict';
import {assessConfirmatoryResult as audit} from '../src/lib/multiomics/confirmatory-result-gate.js';

const method=(count=2)=>({status:'ok',referenceTableComplete:true,
  referenceTableRows:count,referenceTableCsv:'feature,logFC,P.Value,adj.P.Val\nX,0.3,0.1,0.2\nY,0.2,0.2,0.3\n'});
const validLayer=(layer)=>({
  mode:'reference-r',numericalProvenance:'exact_native_r_reference_table',
  primaryReferenceMethodKey:layer+'_differential',
  primaryReferenceTableRows:2,rows:[{feature:'X'},{feature:'Y'}]
});
const valid=()=>({
  protocol:{analysisIntent:'confirmatory',objective:'groups'},
  layers:{proteomics:validLayer('proteomics'),metabolomics:validLayer('metabolomics')},
  numericSourceByLayer:{
    proteomics:'native_r_complete_table',metabolomics:'native_r_complete_table'},
  displayedReferenceMethodCount:2,
  referenceBackend:{status:'ok',referenceMethodExecuted:true,
    methods:{
      proteomics_differential:method(),
      metabolomics_differential:method()
    }}
});
const reject=(v,code,layers=[])=>{
  const a=audit(v);
  assert.equal(a.ok,false);
  assert.equal(a.code,code);
  for(const layer of layers) assert.ok(a.layers.includes(layer),layer);
  assert.ok(a.message && a.messageEn);
};
assert.equal(audit(valid()).code,'native_tables_complete_not_certified');
assert.match(audit(valid()).messageEn,/independent review/);
const exploratory=valid();exploratory.protocol.analysisIntent='exploratory';
assert.equal(audit(exploratory).code,'exploratory_not_certified');
const missingResult=valid();missingResult.referenceBackend.status='unavailable';
reject(missingResult,'reference_r_not_executed');
const unrelatedFit=valid();unrelatedFit.referenceBackend.referenceMethodExecuted=false;
reject(unrelatedFit,'reference_r_not_executed');
const onlyFirst=valid();
onlyFirst.layers.metabolomics={rows:[{feature:'X'}],mode:'browser'};
onlyFirst.numericSourceByLayer.metabolomics='browser_exploratory';
onlyFirst.displayedReferenceMethodCount=1;
reject(onlyFirst,'incomplete_native_layer_projection',['metabolomics']);
const missingCSV=valid();delete missingCSV.referenceBackend.methods.metabolomics_differential.referenceTableCsv;
reject(missingCSV,'incomplete_native_layer_projection',['metabolomics']);
const incomplete=valid();incomplete.referenceBackend.methods.proteomics_differential.referenceTableComplete=false;
reject(incomplete,'incomplete_native_layer_projection',['proteomics']);
const forgedCount=valid();forgedCount.referenceBackend.methods.metabolomics_differential.referenceTableRows=17;
reject(forgedCount,'incomplete_native_layer_projection',['metabolomics']);
const detached=valid();detached.layers.proteomics.primaryReferenceMethodKey='metabolomics_differential';
reject(detached,'incomplete_native_layer_projection',['proteomics']);
const mixed=valid();mixed.layers.transcriptomics=validLayer('transcriptomics');
mixed.layers.transcriptomics.primaryReferenceMethodKey='transcriptomics_differential_deseq2';
mixed.referenceBackend.methods.transcriptomics_differential_deseq2=method();
mixed.numericSourceByLayer.transcriptomics='native_r_complete_table';
mixed.displayedReferenceMethodCount=3;
assert.equal(audit(mixed).ok,true);
for(const objective of ['time','outcome','explore','unknown']) {
  const unsupported=valid();unsupported.protocol.objective=objective;
  reject(unsupported,'native_projection_not_supported');
}
const none=valid();none.layers={};
reject(none,'no_omic_layers');
console.log('PASS: fail-closed native confirmatory result audit, multi-omic layer coverage and negative controls');
