import assert from 'node:assert/strict';
import {assessConfirmatoryResult as audit} from '../src/lib/multiomics/confirmatory-result-gate.js';

const method=(count=2)=>({status:'ok',referenceTableComplete:true,
  referenceTableRows:count,referenceTableCsv:'feature,effect,pValue,qValue\nX,0.3,0.1,0.2\nY,0.2,0.2,0.3\n'});
const validLayer=(layer,key)=>({
  mode:'reference-r',numericalProvenance:'exact_native_r_reference_table',
  primaryReferenceMethodKey:key,
  primaryReferenceTableRows:2,rows:[{feature:'X'},{feature:'Y'}]
});
const base=(protocol, keyForLayer)=>({
  protocol:{analysisIntent:'confirmatory',...protocol},
  layers:{
    proteomics:validLayer('proteomics',keyForLayer('proteomics')),
    metabolomics:validLayer('metabolomics',keyForLayer('metabolomics'))
  },
  numericSourceByLayer:{proteomics:'native_r_complete_table',metabolomics:'native_r_complete_table'},
  displayedReferenceMethodCount:2,
  referenceBackend:{status:'ok',referenceMethodExecuted:true,methods:{
    [keyForLayer('proteomics')]:method(),
    [keyForLayer('metabolomics')]:method()
  }}
});
const group=()=>base({objective:'groups'},layer=>layer+'_differential');
const longitudinal=()=>base({objective:'time'},layer=>layer+'_longitudinal');
const survival=()=>base({objective:'outcome',outcomeType:'survival'},layer=>layer+'_survival');

const reject=(v,code,layers=[])=>{
  const a=audit(v);
  assert.equal(a.ok,false);
  assert.equal(a.code,code);
  for(const layer of layers) assert.ok(a.layers.includes(layer),layer);
  assert.ok(a.message && a.messageEn);
};

for(const fixture of [group(),longitudinal(),survival()]){
  const a=audit(fixture);
  assert.equal(a.ok,true);
  assert.equal(a.code,'native_tables_complete_not_certified');
  assert.match(a.messageEn,/independent review/);
}

const exploratory=group();exploratory.protocol.analysisIntent='exploratory';
assert.equal(audit(exploratory).code,'exploratory_not_certified');

const missingResult=survival();missingResult.referenceBackend.status='unavailable';
reject(missingResult,'reference_r_not_executed');

const unrelatedFit=longitudinal();unrelatedFit.referenceBackend.referenceMethodExecuted=false;
reject(unrelatedFit,'reference_r_not_executed');

const onlyFirst=survival();
onlyFirst.layers.metabolomics={rows:[{feature:'X'}],mode:'browser'};
onlyFirst.numericSourceByLayer.metabolomics='browser_exploratory';
onlyFirst.displayedReferenceMethodCount=1;
reject(onlyFirst,'incomplete_native_layer_projection',['metabolomics']);

const missingCSV=longitudinal();
delete missingCSV.referenceBackend.methods.metabolomics_longitudinal.referenceTableCsv;
reject(missingCSV,'incomplete_native_layer_projection',['metabolomics']);

const incomplete=group();
incomplete.referenceBackend.methods.proteomics_differential.referenceTableComplete=false;
reject(incomplete,'incomplete_native_layer_projection',['proteomics']);

const forgedCount=survival();
forgedCount.referenceBackend.methods.metabolomics_survival.referenceTableRows=17;
reject(forgedCount,'incomplete_native_layer_projection',['metabolomics']);

const wrongFamily=survival();
wrongFamily.layers.proteomics.primaryReferenceMethodKey='proteomics_differential';
wrongFamily.referenceBackend.methods.proteomics_differential=method();
reject(wrongFamily,'incomplete_native_layer_projection',['proteomics']);

const mixed=group();
mixed.layers.transcriptomics=validLayer('transcriptomics','transcriptomics_differential_deseq2');
mixed.referenceBackend.methods.transcriptomics_differential_deseq2=method();
mixed.numericSourceByLayer.transcriptomics='native_r_complete_table';
mixed.displayedReferenceMethodCount=3;
assert.equal(audit(mixed).ok,true);

for(const protocol of [
  {objective:'outcome',outcomeType:'binary'},
  {objective:'outcome',outcomeType:'continuous'},
  {objective:'explore'},
  {objective:'unknown'}
]) {
  const unsupported=group();unsupported.protocol={analysisIntent:'confirmatory',...protocol};
  reject(unsupported,'native_projection_not_supported');
}

const none=group();none.layers={};
reject(none,'no_omic_layers');

console.log('PASS: fail-closed native confirmatory group, longitudinal and survival result audit with negative controls');
