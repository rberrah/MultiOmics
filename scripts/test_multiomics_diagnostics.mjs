import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';
import { evaluateConfirmatoryReadiness } from '../src/lib/multiomics/scientific-assurance.js';

function matrixFile(name, assayIds, shift = 0) {
  const header = ['feature_id', ...assayIds].join(',');
  const rows = [];
  for (let g = 1; g <= 24; g += 1) {
    const values = assayIds.map((_, index) => {
      const groupEffect = index >= assayIds.length / 2 && g <= 6 ? 1.5 : 0;
      return (3 + g * 0.07 + index * 0.03 + groupEffect + shift).toFixed(3);
    });
    rows.push([`F${g}`, ...values].join(','));
  }
  return new File([[header, ...rows].join('\n')], name, { type: 'text/csv' });
}

const subjects = ['P1','P2','P3','P4','P5','P6'];
const conditions = ['A','A','A','B','B','B'];
const rnaAssays = subjects.map((id) => `R_${id}`);
const protAssays = subjects.map((id) => `P_${id}`);
const metadataRows = [];
for (let i = 0; i < subjects.length; i += 1) {
  metadataRows.push({
    subject_id: subjects[i], sample_id: `${subjects[i]}_T0`, assay_id: rnaAssays[i],
    omic: 'transcriptomics', condition: conditions[i], timepoint: 'T0', batch: ''
  });
  metadataRows.push({
    subject_id: subjects[i], sample_id: `${subjects[i]}_T0`, assay_id: protAssays[i],
    omic: 'proteomics', condition: conditions[i], timepoint: 'T0', batch: ''
  });
}

const files = {
  metadata: new File(['placeholder'], 'samples.csv', { type: 'text/csv' }),
  transcriptomics: matrixFile('rna.csv', rnaAssays),
  proteomics: matrixFile('protein.csv', protAssays, 0.4),
  metabolomics: null
};

const columnMapping = {
  subject_id: 'subject_id', sample_id: 'sample_id', assay_id: 'assay_id', omic: 'omic',
  condition: 'condition', timepoint: 'timepoint', batch: 'batch', technical_replicate: '',
  outcome: '', survival_time: '', survival_event: '', sample_type: '', injection_order: ''
};

const result = await runDeterministicAnalysis({
  files,
  metadataRows,
  columnMapping,
  protocol: {
    objective: 'groups',
    designType: 'independent',
    longitudinal: false,
    batchKnown: 'no',
    partialOmicsExpected: 'no',
    covariateColumns: [],
    organism: 'human'
  },
  dataTypes: { transcriptomics: 'log_expression', proteomics: 'log_intensity' },
  identifierTypes: { transcriptomics: 'unknown', proteomics: 'unknown' },
  resolveIdentifiers: false,
  useReactome: false
});

assert.ok(result.preAnalysisDiagnostics);
assert.equal(result.engine.version, '1.8.0');
assert.equal(result.preAnalysisDiagnostics.groupBalance.subjectsPerGroup.A, 3);
assert.equal(result.preAnalysisDiagnostics.groupBalance.subjectsPerGroup.B, 3);
assert.equal(result.preAnalysisDiagnostics.groupBalance.smallestToLargestRatio, 1);
assert.equal(result.preAnalysisDiagnostics.groupBalance.status, 'good');
assert.equal(result.preAnalysisDiagnostics.omicsOverlap.subjectsInEveryLoadedLayer, 6);
assert.equal(result.preAnalysisDiagnostics.omicsOverlap.fractionInEveryLoadedLayer, 1);
assert.equal(result.preAnalysisDiagnostics.power.status, 'not_calculated');
assert.match(result.preAnalysisDiagnostics.power.note, /effect size/i);

// Small independent groups remain analysable but cannot be presented as a
// high-confidence confirmatory experiment simply because a p/q value exists.
assert.equal(result.preAnalysisDiagnostics.blockers.length, 0);
assert.ok(result.preAnalysisDiagnostics.warnings.some((item) => /fewer than five independent subjects/i.test(item)));
assert.equal(result.preAnalysisDiagnostics.status, 'review_required');

// Advanced methods are eligibility-gated independently of whether their R
// package happens to be installed on the machine.
assert.equal(result.preAnalysisDiagnostics.methodEligibility.mofa2.status, 'not_recommended');
assert.match(result.preAnalysisDiagnostics.methodEligibility.mofa2.reasons.join(' '), /more than 15 shared samples/i);
assert.equal(result.preAnalysisDiagnostics.methodEligibility.diablo.status, 'eligible_with_internal_cv');
assert.equal(result.preAnalysisDiagnostics.methodEligibility.diablo.smallestClass, 3);
assert.match(result.preAnalysisDiagnostics.methodEligibility.diablo.interpretation, /external validation/i);
assert.equal(result.metadataSummary.preAnalysisDiagnostics.status, result.preAnalysisDiagnostics.status);


// Confirmatory eligibility is method- and layer-specific. No automatic label
// may imply a publication certificate or validated clinical biomarker.
const candidate = {
  ...result,
  protocol:{...result.protocol,analysisIntent:'confirmatory'},
  referenceBackend:{status:'ok',methods:{
    transcriptomics_differential:{status:'ok',method:'limma'},
    proteomics_differential:{status:'ok',method:'limma'}
  }}
};
const noR = evaluateConfirmatoryReadiness({...candidate,referenceBackend:null});
assert.equal(noR.status,'blocked');
assert.equal(noR.certified,false);
const allR = evaluateConfirmatoryReadiness(candidate);
assert.equal(allR.status,'independent_review_required');
assert.equal(allR.certified,false);
assert.ok(allR.checks.some(item=>item.code==='preregistration' && item.status==='review'));
const partialR = evaluateConfirmatoryReadiness({
  ...candidate,referenceBackend:{status:'ok',methods:{transcriptomics_differential:{status:'ok'}}}
});
assert.equal(partialR.status,'blocked','a reference for only one omic must never pass');
assert.ok(partialR.checks.some(item=>item.code==='reference_proteomics' && item.status==='blocked'));
assert.equal(evaluateConfirmatoryReadiness(candidate,{demo:true}).status,'blocked');
assert.equal(evaluateConfirmatoryReadiness({
  ...candidate,protocol:{...candidate.protocol,objective:'outcome',outcomeType:'binary'}
}).status,'blocked');

// Native longitudinal and survival contracts are eligible only for independent
// review, never automatic certification. Assumption diagnostics stay visible.
const longitudinalCandidate={
  ...candidate,
  protocol:{...candidate.protocol,objective:'time',designType:'repeated',longitudinal:true},
  referenceBackend:{status:'ok',methods:{
    transcriptomics_longitudinal:{status:'ok',summary:{randomSlopeFallbacks:1}},
    proteomics_longitudinal:{status:'ok',summary:{randomSlopeFallbacks:0}}
  }}
};
const longitudinalReady=evaluateConfirmatoryReadiness(longitudinalCandidate);
assert.equal(longitudinalReady.status,'independent_review_required');
assert.equal(longitudinalReady.certified,false);
assert.ok(longitudinalReady.checks.some(x=>
  x.code==='longitudinal_assumptions_transcriptomics'&&x.status==='review'));
assert.ok(longitudinalReady.checks.some(x=>
  x.code==='longitudinal_random_slope_fallback_transcriptomics'&&x.status==='review'));

const survivalCandidate={
  ...candidate,
  protocol:{...candidate.protocol,objective:'outcome',outcomeType:'survival',
    designType:'independent',longitudinal:false},
  referenceBackend:{status:'ok',methods:{
    transcriptomics_survival:{status:'ok',summary:{statusCounts:{ph_assumption_warning:1}}},
    proteomics_survival:{status:'ok',summary:{statusCounts:{ph_assumption_warning:0}}}
  }}
};
const survivalReady=evaluateConfirmatoryReadiness(survivalCandidate);
assert.equal(survivalReady.status,'independent_review_required');
assert.equal(survivalReady.certified,false);
assert.ok(survivalReady.checks.some(x=>
  x.code==='survival_assumptions_transcriptomics'&&x.status==='review'));
assert.ok(survivalReady.checks.some(x=>
  x.code==='survival_ph_warning_transcriptomics'&&x.status==='review'));

const rawLongitudinal={
  ...longitudinalCandidate,
  layers:{
    ...longitudinalCandidate.layers,
    transcriptomics:{
      ...longitudinalCandidate.layers.transcriptomics,
      qc:{
        ...longitudinalCandidate.layers.transcriptomics.qc,
        preprocessingAudit:{
          ...longitudinalCandidate.layers.transcriptomics.qc.preprocessingAudit,
          declaredValueType:'raw_counts'
        }
      }
    }
  }
};
assert.equal(evaluateConfirmatoryReadiness(rawLongitudinal).status,'blocked');
assert.ok(evaluateConfirmatoryReadiness(rawLongitudinal).checks.some(x=>
  x.code==='raw_count_reference_transcriptomics'&&x.status==='blocked'));

const missingnessBlocks = evaluateConfirmatoryReadiness({
  ...candidate,
  layers: {
    ...candidate.layers,
    proteomics: {
      ...candidate.layers.proteomics,
      qc: {
        ...candidate.layers.proteomics.qc,
        differentialMissingness:{flaggedFeatures:2}
      }
    }
  }
});
assert.equal(missingnessBlocks.status,'blocked');
assert.ok(missingnessBlocks.checks.some(x=>
  x.code==='missingness_imbalance_proteomics' && x.status==='blocked'),
  'Independent R success must not silently override differential missingness');
const balancedMnarRisk = evaluateConfirmatoryReadiness({
  ...candidate,
  layers:{
    ...candidate.layers,
    proteomics:{
      ...candidate.layers.proteomics,
      qc:{
        ...candidate.layers.proteomics.qc,
        medianMissingFraction:0.26,
        differentialMissingness:{flaggedFeatures:0}
      }
    }
  }
});
assert.equal(balancedMnarRisk.status,'blocked',
  'Equal group missingness rates can conceal opposite-tail MNAR');
assert.ok(balancedMnarRisk.checks.some(x=>x.code==='high_missingness_proteomics'&&x.status==='blocked'));
const mildMissingReview = evaluateConfirmatoryReadiness({
  ...candidate,layers:{
    ...candidate.layers,
    proteomics:{
      ...candidate.layers.proteomics,
      qc:{...candidate.layers.proteomics.qc,medianMissingFraction:0.08}
    }
  }
});
assert.ok(mildMissingReview.checks.some(x=>x.code==='missingness_mechanism_proteomics' && x.status==='review'),
  'Low missingness also needs assumption review even when group rates match');
const mnarImputed = evaluateConfirmatoryReadiness({
  ...candidate,
  layers:{
    ...candidate.layers,
    metabolomics:{
      qc:{
        preprocessingAudit:{declaredValueType:'log_abundance'},
        msQc:{mnar:{imputedValues:8}}
      }
    }
  },
  referenceBackend:{status:'ok',methods:{
    transcriptomics_differential:{status:'ok'},
    proteomics_differential:{status:'ok'},
    metabolomics_differential:{status:'ok'}
  }}
});
assert.ok(mnarImputed.checks.some(x=>x.code==='mnar_imputation_metabolomics' &&
  x.status==='blocked'),'Imputed assumed MNAR values must be explicitly blocked');

// A successful raw-count DESeq2 status CANNOT bypass the independent-subject
// R1 preflight; only an exact server-side audit permits the review pathway.
const rawCountCandidate = {
  ...candidate,
  protocol:{...candidate.protocol,rnaPrimaryContrast:['condition','treated','control']},
  layers:{transcriptomics:{
    qc:{preprocessingAudit:{declaredValueType:'raw_counts'},
      medianMissingFraction:0,differentialMissingness:{flaggedFeatures:0}}
  }},
  referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{status:'ok',method:'DESeq2'}
  }}
};
const missingR1Audit=evaluateConfirmatoryReadiness(rawCountCandidate);
assert.equal(missingR1Audit.status,'blocked');
assert.ok(missingR1Audit.checks.some(x=>
  x.code==='rnaseq_confirmatory_input_contract'&&x.status==='blocked'));
const withR1Audit={
  ...rawCountCandidate,
  referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      status:'ok',method:'DESeq2',
      confirmatoryInputContract:{
        status:'eligible_for_reference_fitting_only',certified:false,
        independent_subjects:40,rejected_silent_subject_loss:true,
        condition_contrast:['condition','treated','control']
      },
      assayIntegrityContract:{
        status:'assay_identifiers_exact',assay_loss_detected:false,
        measured_assays:40,certified:false
      }
    }
  }}
};
const passedR1Audit=evaluateConfirmatoryReadiness(withR1Audit);
assert.equal(passedR1Audit.status,'independent_review_required');
assert.equal(passedR1Audit.certified,false);
assert.ok(passedR1Audit.checks.some(x=>
  x.code==='rnaseq_confirmatory_input_contract'&&x.status==='pass'));
const reversedReferenceEffect={
  ...withR1Audit,referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      ...withR1Audit.referenceBackend.methods.transcriptomics_differential_deseq2,
      confirmatoryInputContract:{
        ...withR1Audit.referenceBackend.methods.transcriptomics_differential_deseq2.confirmatoryInputContract,
        condition_contrast:['condition','control','treated']
      }
    }
  }}
};
assert.equal(evaluateConfirmatoryReadiness(reversedReferenceEffect).status,'blocked',
  'A reversed R effect must never be reported as matching a prespecified RNA contrast');
const missingAssayEvidence={
  ...withR1Audit,referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      ...withR1Audit.referenceBackend.methods.transcriptomics_differential_deseq2,
      assayIntegrityContract:null
    }
  }}
};
assert.equal(evaluateConfirmatoryReadiness(missingAssayEvidence).status,'blocked',
  'A strict DESeq2 contract must not hide lost assay CSV rows');
const brokenAssayEvidence={
  ...withR1Audit,referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      ...withR1Audit.referenceBackend.methods.transcriptomics_differential_deseq2,
      assayIntegrityContract:{
        status:'assay_identifiers_exact',assay_loss_detected:true,
        measured_assays:39,certified:false
      }
    }
  }}
};
assert.equal(evaluateConfirmatoryReadiness(brokenAssayEvidence).status,'blocked');

for(const faulty of [
  {certified:true},{independent_subjects:5},
  {rejected_silent_subject_loss:false},{status:'unverified'}
]){
  const bad={...withR1Audit,referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      status:'ok',confirmatoryInputContract:{
        ...withR1Audit.referenceBackend.methods.transcriptomics_differential_deseq2.confirmatoryInputContract,
        ...faulty
      }
    }
  }}};
  assert.equal(evaluateConfirmatoryReadiness(bad).status,'blocked',
    'R1 eligibility cannot be accepted from unsupported certificate or identity');
}

// Paired confirmatory preparation is not automatically rejected once the REAL
// server-side subject-blocked contract is verified. It still requires human
// review and must never claim that DESeq2 q-values are universally calibrated.
const pairedRContract={
  status:'eligible_for_reference_fitting_only',
  pairing:'same_subject_two_conditions',
  biological_pairs:4,independent_biological_subjects:4,
  model_samples:8,residual_degrees_of_freedom:3,
  sample_identifiers_exact:true,
  all_subject_condition_cells_observed_once:true,
  contrast_was_prespecified:true,
  explicit_primary_contrast:['condition','treated','control'],
  certified:false
};
const pairedRReview={
  ...withR1Audit,
  protocol:{...withR1Audit.protocol,designType:'paired'},
  metadataSummary:{...withR1Audit.metadataSummary,subjects:4,conditions:['control','treated']},
  referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      status:'ok',method:'DESeq2',
      confirmatoryInputContract:pairedRContract,
      pairedDesignContract:pairedRContract,
      assayIntegrityContract:{
        status:'assay_identifiers_exact',assay_loss_detected:false,
        measured_assays:8,certified:false
      }
    }
  }}
};
const validPaired=evaluateConfirmatoryReadiness(pairedRReview);
assert.equal(validPaired.status,'independent_review_required');
assert.equal(validPaired.certified,false);
assert.ok(validPaired.checks.some(x=>x.code==='design'&&x.status==='pass'));
assert.ok(validPaired.checks.some(x=>
  x.code==='rnaseq_confirmatory_input_contract'&&x.status==='pass'));
assert.ok(validPaired.checks.some(x=>
  x.code==='rnaseq_null_fdr_calibration'&&x.status==='review'));
for(const altered of [
  {biological_pairs:3},
  {model_samples:7},
  {all_subject_condition_cells_observed_once:false},
  {sample_identifiers_exact:false},
  {contrast_was_prespecified:false},
  {explicit_primary_contrast:['condition','control','treated']},
  {certified:true},
  {residual_degrees_of_freedom:1}
]){
  const contract={...pairedRContract,...altered};
  const bad={...pairedRReview,referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      ...pairedRReview.referenceBackend.methods.transcriptomics_differential_deseq2,
      confirmatoryInputContract:contract,pairedDesignContract:contract
    }
  }}};
  assert.equal(evaluateConfirmatoryReadiness(bad).status,'blocked',
    'Missing/invalid subject-blocked evidence cannot reach independent review');
}
const missingPairedAssays={
  ...pairedRReview,referenceBackend:{status:'ok',methods:{
    transcriptomics_differential_deseq2:{
      ...pairedRReview.referenceBackend.methods.transcriptomics_differential_deseq2,
      assayIntegrityContract:{
        status:'assay_identifiers_exact',assay_loss_detected:false,
        measured_assays:7,certified:false
      }
    }
  }}
};
assert.equal(evaluateConfirmatoryReadiness(missingPairedAssays).status,'blocked');
const wrongPairedSubjectCount={
  ...pairedRReview,metadataSummary:{...pairedRReview.metadataSummary,subjects:8}
};
assert.equal(evaluateConfirmatoryReadiness(wrongPairedSubjectCount).status,'blocked');

assert.equal(evaluateConfirmatoryReadiness(result).status,'not_requested');

console.log('multiomics pre-analysis diagnostics: PASS');
