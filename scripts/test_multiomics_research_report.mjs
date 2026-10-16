import assert from 'node:assert/strict';
import {buildResearchMethodsManifest} from '../src/lib/multiomics/research-report.js';
const example={
  engine:{version:'3.2'},
  generatedAt:'2026-10-12T15:00:00Z',
  protocol:{objective:'groups',designType:'independent',analysisIntent:'confirmatory',
    organism:'mouse',sampleOverlap:'same_specimen',covariateColumns:['batch']},
  metadataSummary:{subjects:12},
  inputManifest:{sourceSha256:'a'.repeat(64)},
  layers:{
    transcriptomics:{inferenceMethod:'Welch screening',
      qc:{featuresBefore:300,featuresAfter:220,medianMissingFraction:0,
      sampleMetrics:[
        {assayId:'PRIVATE_ASSAY_001',missingFraction:0},
        {assayId:'PRIVATE_ASSAY_002',missingFraction:.38},
        {assayId:'PRIVATE_ASSAY_003',missingFraction:0}],
      preprocessingAudit:{declaredValueType:'raw_counts'},warnings:['review']}},
    proteomics:{qc:{featuresBefore:100,featuresAfter:90,
      medianMissingFraction:0}},
    metabolomics:{qc:{medianMissingFraction:.25,msQc:{mnar:{imputedValues:10}}}}
  },
  referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{status:'ok'},
    proteomics_differential:{status:'blocked'}}},
  scientificAssurance:{status:'independent_review_required'},
  // This payload is private and must never appear in methods-only audit.
  privatePatient:{name:'SENSITIVE PATIENT IDENTIFIER',data:[9911,4423]},
  patients:[{id:'PERSON_001',expression:3e8}],
  fullRawMatrix:[1,3,4]
};
const audit=buildResearchMethodsManifest(example);
assert.equal(audit.format,'pmx-omics-methods-audit/v1');
assert.equal(audit.clinical_validation,false);
assert.equal(audit.scientific_certification,false);
assert.equal(audit.omics.transcriptomics.omitted_features,80);
assert.equal(audit.omics.transcriptomics.median_missing_fraction,0);
assert.equal(audit.omics.transcriptomics.maximum_sample_missing_fraction,.38);
assert.equal(audit.omics.transcriptomics.sample_profiles_with_missing_features,1);
assert.equal(audit.omics.transcriptomics.missingness_review_required,true);
assert.equal(audit.omics.transcriptomics.reference_executed,true);
assert.equal(audit.omics.proteomics.reference_executed,false);
assert.equal(audit.omics.metabolomics.omitted_features,null);
assert.equal(audit.omics.metabolomics.mnar_imputation_applied,true);
assert.equal(audit.independent_subjects,12);
assert.equal(audit.input_manifest_present,true);
const bytes=JSON.stringify(audit);
for(const value of ['SENSITIVE PATIENT IDENTIFIER','PERSON_001','PRIVATE_ASSAY_002','9911','4423','300000000'])
  assert.equal(bytes.includes(value),false,'Methods-only export leaked a raw sensitive field');
assert.throws(()=>buildResearchMethodsManifest({protocol:{}}),/completed multi-omics analysis/);
console.log('MULTIOMICS METHODS-ONLY EXPORT AND PRIVACY TEST: PASS');
