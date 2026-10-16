import fs from 'node:fs/promises';
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { parseDelimited, runDeterministicAnalysis, resolveMetaboliteIdentifier, resolveGeneIdentifier, resolveProteinIdentifier } from '../src/lib/multiomics/deterministic.js';

const root = new URL('../static/multiomics/', import.meta.url);

async function load(name) {
  const text = await fs.readFile(new URL(name, root), 'utf8');
  return new File([text], name, { type: 'text/csv' });
}

const metadata = await load('demo_metadata.csv');
const transcriptomics = await load('demo_transcriptomics.csv');
const proteomics = await load('demo_proteomics.csv');
const metabolomics = await load('demo_metabolomics.csv');

const parsed = parseDelimited(await metadata.text());
const mapping = {
  subject_id: 'subject_id',
  sample_id: 'sample_id',
  assay_id: 'assay_id',
  omic: 'omic',
  condition: 'condition',
  timepoint: 'timepoint',
  batch: 'batch',
  technical_replicate: 'technical_replicate',
  outcome: 'outcome',
  survival_time: 'survival_time',
  survival_event: 'survival_event'
};

const result = await runDeterministicAnalysis({
  files: { metadata, transcriptomics, proteomics, metabolomics },
  metadataRows: parsed.rows,
  columnMapping: mapping,
  protocol: {
    organism: 'human',
    objective: 'time',
    longitudinal: true,
    designType: 'repeated',
    studySetting: 'clinical_interventional',
    groupCount: '2',
    sampleOverlap: 'same_specimen'
  },
  dataTypes: {
    transcriptomics: 'raw_counts',
    proteomics: 'log_intensity',
    metabolomics: 'peak_area'
  },
  useReactome: false,
  resolveIdentifiers: false
});

assert.equal(result.metadataSummary.subjects, 8);
assert.equal(result.metadataSummary.samples, 16);
assert.equal(result.metadataSummary.assays, 49);
assert.equal(result.layers.transcriptomics.replicateGroups.length, 1);
assert.ok(result.layers.transcriptomics.rows.length >= 5);
assert.ok(result.layers.proteomics.rows.length >= 4);
assert.ok(result.layers.metabolomics.rows.length >= 4);
assert.equal(result.metadataSummary.batchAudit.transcriptomics.status, 'single_batch');
assert.equal(result.metadataSummary.batchAudit.proteomics.status, 'single_batch');
assert.equal(result.metadataSummary.batchAudit.metabolomics.status, 'single_batch');

 
// Identifier integrity: contradictory mappings and duplicate matrix variables
// must fail before statistics can silently overwrite rows.
const basicDemoOptions = {
  files:{metadata,transcriptomics,proteomics,metabolomics},
  metadataRows:parsed.rows,columnMapping:mapping,
  protocol:{organism:'human',objective:'time',longitudinal:true,designType:'repeated',
    studySetting:'clinical_interventional',groupCount:'2',sampleOverlap:'same_specimen'},
  dataTypes:{transcriptomics:'raw_counts',proteomics:'log_intensity',metabolomics:'peak_area'},
  useReactome:false,resolveIdentifiers:false
};
await assert.rejects(() => runDeterministicAnalysis({
  ...basicDemoOptions,
  metadataRows:[...parsed.rows, {...parsed.rows[0]}]
}), /duplicate assay_id/i);
const subjectAtFirstSample = parsed.rows[0].subject_id;
await assert.rejects(() => runDeterministicAnalysis({
  ...basicDemoOptions,
  metadataRows:parsed.rows.map((row,i)=>i===0
    ? {...row,subject_id:'INVALID_' + subjectAtFirstSample} : row)
}), /conflicting sample_id/i);
const originalRnaText = await transcriptomics.text();
const originalRnaLines = originalRnaText.trimEnd().split(/\r?\n/);
const duplicateRnaFeature = new File(
  [originalRnaText.trimEnd() + '\n' + originalRnaLines[1] + '\n'],
  'duplicate_feature.csv',{type:'text/csv'});
await assert.rejects(() => runDeterministicAnalysis({
  ...basicDemoOptions,files:{...basicDemoOptions.files,transcriptomics:duplicateRnaFeature}
}), /duplicate feature_id/i);


// Every layer must expose the declared input scale and the exact preprocessing decision.
assert.equal(result.layers.transcriptomics.qc.preprocessingAudit.declaredValueType, 'raw_counts');
assert.equal(result.layers.transcriptomics.qc.preprocessingAudit.outputScale, 'log2');
assert.equal(result.layers.transcriptomics.qc.preprocessingAudit.status, 'compatible');
assert.equal(result.layers.transcriptomics.qc.preprocessingAudit.nonIntegerValues, 0);
assert.equal(result.layers.metabolomics.qc.preprocessingAudit.declaredValueType, 'peak_area');
assert.equal(result.layers.metabolomics.qc.preprocessingAudit.medianCentered, true);
assert.ok(result.layers.metabolomics.qc.preprocessingAudit.pseudoCount > 0);
assert.match(result.layers.metabolomics.qc.preprocessingAudit.policy, /does not infer or silently apply/i);

for (const layer of ['transcriptomics','proteomics','metabolomics']) {
  assert.equal(result.layers[layer].mode, 'random-intercept-longitudinal-model');
  assert.deepEqual(result.layers[layer].groupSizes, [4,4]);
  assert.ok(result.layers[layer].selected.length > 0);
}

const ido = result.layers.transcriptomics.rows.find((row) => row.feature === 'IDO1');
assert.ok(ido);
assert.ok(Number.isFinite(ido.effect));
assert.ok(ido.foldRatio > 1);
assert.ok(result.crossOmics.testedPairs > 0);
assert.ok(result.crossOmics.pairs.some((pair) => ['sign_reversal','gained_in_comparison','strengthened_in_comparison'].includes(pair.pattern)));



// Outcome models with repeated omics must use an explicitly selected molecular time point.
await assert.rejects(
  () => runDeterministicAnalysis({
    files: { metadata, transcriptomics, proteomics, metabolomics },
    metadataRows: parsed.rows,
    columnMapping: mapping,
    protocol: {
      organism: 'human',
      objective: 'outcome',
      outcomeType: 'binary',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_interventional',
      groupCount: '1',
      sampleOverlap: 'same_specimen',
      batchKnown: 'yes',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'raw_counts',
      proteomics: 'log_intensity',
      metabolomics: 'peak_area'
    },
    useReactome: false,
    resolveIdentifiers: false
  }),
  /explicit outcomeTimepoint/
);

const demoOutcome = await runDeterministicAnalysis({
  files: { metadata, transcriptomics, proteomics, metabolomics },
  metadataRows: parsed.rows,
  columnMapping: mapping,
  protocol: {
    organism: 'human',
    objective: 'outcome',
    outcomeType: 'binary',
    outcomeTimepoint: 'T0',
    longitudinal: false,
    designType: 'independent',
    studySetting: 'clinical_interventional',
    groupCount: '1',
    sampleOverlap: 'same_specimen',
    batchKnown: 'yes',
    covariateColumns: []
  },
  dataTypes: {
    transcriptomics: 'raw_counts',
    proteomics: 'log_intensity',
    metabolomics: 'peak_area'
  },
  useReactome: false,
  resolveIdentifiers: false
});
assert.equal(demoOutcome.layers.transcriptomics.mode, 'outcome-binary');
assert.equal(demoOutcome.protocol.outcomeTimepoint, 'T0');

// Declared measurement types are audited instead of silently guessed.
{
  const originalRna = await transcriptomics.text();
  const fractionalRna = new File(
    [originalRna.replace('STAT1,110,', 'STAT1,110.5,')],
    'fractional_counts.csv',
    { type: 'text/csv' }
  );
  const fractional = await runDeterministicAnalysis({
    files: { metadata, transcriptomics: fractionalRna, proteomics, metabolomics },
    metadataRows: parsed.rows,
    columnMapping: mapping,
    protocol: {
      organism: 'human', objective: 'time', longitudinal: true, designType: 'repeated',
      studySetting: 'clinical_interventional', groupCount: '2', sampleOverlap: 'same_specimen'
    },
    dataTypes: {
      transcriptomics: 'raw_counts',
      proteomics: 'log_intensity',
      metabolomics: 'peak_area'
    },
    useReactome: false,
    resolveIdentifiers: false
  });
  const fractionalAudit = fractional.layers.transcriptomics.qc.preprocessingAudit;
  assert.equal(fractionalAudit.status, 'review_required');
  assert.equal(fractionalAudit.nonIntegerValues, 1);
  assert.ok(fractional.layers.transcriptomics.qc.warnings.some((warning) => /non-integer/i.test(warning)));

  const negativeRna = new File(
    [originalRna.replace('STAT1,110,', 'STAT1,-1,')],
    'negative_counts.csv',
    { type: 'text/csv' }
  );
  await assert.rejects(
    () => runDeterministicAnalysis({
      files: { metadata, transcriptomics: negativeRna, proteomics, metabolomics },
      metadataRows: parsed.rows,
      columnMapping: mapping,
      protocol: {
        organism: 'human', objective: 'time', longitudinal: true, designType: 'repeated',
        studySetting: 'clinical_interventional', groupCount: '2', sampleOverlap: 'same_specimen'
      },
      dataTypes: {
        transcriptomics: 'raw_counts',
        proteomics: 'log_intensity',
        metabolomics: 'peak_area'
      },
      useReactome: false,
      resolveIdentifiers: false
    }),
    /raw_counts.*negative value/i
  );
}

// Advanced metabolomics MS QC: blanks, pooled-QC drift, RSD and explicit MNAR handling.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,sample_type,injection_order'];
  const rnaHeader = ['feature_id'];
  const metHeader = ['feature_id'];
  const rnaSignal = ['RNA_MS_SIGNAL'];
  const rnaAnchor = ['RNA_MS_ANCHOR'];
  const metSignal = ['MS_TRUE_SIGNAL'];
  const metBlank = ['MS_BLANK_CONTAM'];
  const metUnstable = ['MS_UNSTABLE_QC'];
  const metMnar = ['MS_MNAR'];
  let injection = 1;

  // Biological samples: two groups, same biological subjects across RNA and metabolomics.
  for (const condition of ['control','treatment']) {
    for (let s = 1; s <= 6; s += 1) {
      const subject = (condition === 'control' ? 'MC' : 'MT') + s;
      const rnaId = 'MSR_' + subject;
      rows.push([subject,subject,rnaId,'transcriptomics',condition,'T0','B1','biological',''].join(','));
      rnaHeader.push(rnaId);
      rnaSignal.push(String(5 + (condition === 'treatment' ? 1.5 : 0) + s*0.03));
      rnaAnchor.push(String(3 + s*0.08));

      const metId = 'MSM_' + subject;
      const drift = 1 + injection * 0.025;
      rows.push([subject,subject,metId,'metabolomics',condition,'T0','B1','biological',injection].join(','));
      metHeader.push(metId);
      metSignal.push(String((condition === 'treatment' ? 180 : 100) * drift));
      metBlank.push(String(10 * drift));
      metUnstable.push(String(90 * drift));
      metMnar.push((s <= 2 && condition === 'treatment') ? '' : String(40 * drift));
      injection += 1;
    }
  }

  // Pooled QCs interspersed across the sequence.
  for (let q = 1; q <= 6; q += 1) {
    const id = 'MS_QC_' + q;
    const order = injection++;
    const drift = 1 + order * 0.025;
    rows.push(['QC' + q,'QC' + q,id,'metabolomics','','','B1','pooled_qc',order].join(','));
    metHeader.push(id);
    metSignal.push(String(100 * drift));
    metBlank.push(String(10 * drift));
    metUnstable.push(String((q % 2 ? 45 : 170) * drift));
    metMnar.push(String(40 * drift));
  }

  // Blanks: one feature is deliberately contaminated.
  for (let b = 1; b <= 3; b += 1) {
    const id = 'MS_BLANK_' + b;
    const order = injection++;
    rows.push(['BL' + b,'BL' + b,id,'metabolomics','','','B1','blank',order].join(','));
    metHeader.push(id);
    metSignal.push('1');
    metBlank.push('4');
    metUnstable.push('1');
    metMnar.push('0.5');
  }

  const meta = csvFile('ms_qc_metadata.csv', rows.join('\n'));
  const rna = csvFile('ms_qc_rna.csv', [rnaHeader.join(','),rnaSignal.join(','),rnaAnchor.join(',')].join('\n'));
  const metabolomics = csvFile('ms_qc_metabolomics.csv', [
    metHeader.join(','),metSignal.join(','),metBlank.join(','),metUnstable.join(','),metMnar.join(',')
  ].join('\n'));
  const parsedMeta = parseDelimited(await meta.text());

  const ms = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:null,metabolomics},
    metadataRows:parsedMeta.rows,
    columnMapping:{
      ...mapping,
      sample_type:'sample_type',
      injection_order:'injection_order'
    },
    protocol:{
      organism:'human',objective:'groups',longitudinal:false,designType:'independent',
      studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen',
      batchKnown:'yes',covariateColumns:[],
      msBlankFilter:'yes',msBlankFold:5,
      msQcRsdFilter:'yes',msQcRsdThreshold:0.30,
      msDriftCorrection:'yes',msMnarStrategy:'left_censored'
    },
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'peak_area'},
    useReactome:false,
    resolveIdentifiers:false
  });

  assert.equal(ms.metadataSummary.subjects, 12);
  assert.equal(ms.metadataSummary.technicalQcAssays, 9);
  const msQc = ms.layers.metabolomics.qc.msQc;
  assert.equal(msQc.blankAssays, 3);
  assert.equal(msQc.qcAssays, 6);
  assert.ok(msQc.blankFilteredFeatures >= 1);
  assert.ok(msQc.qcRsdFilteredFeatures >= 1);
  assert.ok(msQc.driftCorrection.applied);
  assert.ok(msQc.driftCorrection.correctedFeatures >= 2);
  assert.ok(msQc.mnar.imputedValues >= 2);
  assert.ok(ms.layers.metabolomics.rows.some((row) => row.feature === 'MS_TRUE_SIGNAL'));
  assert.ok(!ms.layers.metabolomics.rows.some((row) => row.feature === 'MS_BLANK_CONTAM'));
  assert.ok(!ms.layers.metabolomics.rows.some((row) => row.feature === 'MS_UNSTABLE_QC'));

  const msFlagOnly = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:null,metabolomics},
    metadataRows:parsedMeta.rows,
    columnMapping:{
      ...mapping,
      sample_type:'sample_type',
      injection_order:'injection_order'
    },
    protocol:{
      organism:'human',objective:'groups',longitudinal:false,designType:'independent',
      studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen',
      batchKnown:'yes',covariateColumns:[],
      msBlankFilter:'flag',msBlankFold:5,
      msQcRsdFilter:'no',msDriftCorrection:'no',msMnarStrategy:'none'
    },
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'peak_area'},
    useReactome:false,
    resolveIdentifiers:false
  });
  const flagQc = msFlagOnly.layers.metabolomics.qc.msQc;
  assert.equal(flagQc.blankMode, 'flag');
  assert.ok(flagQc.blankFlaggedFeatures >= 1);
  assert.equal(flagQc.blankRemovedFeatures, 0);
  assert.ok(msFlagOnly.layers.metabolomics.rows.some((row) => row.feature === 'MS_BLANK_CONTAM'));
  assert.ok(msFlagOnly.layers.metabolomics.qc.differentialMissingness.flaggedFeatures >= 1);
  assert.ok(msFlagOnly.layers.metabolomics.qc.differentialMissingness.topFeatures.some((item) => item.feature === 'MS_MNAR'));
  assert.ok(msFlagOnly.layers.metabolomics.qc.warnings.some((warning) => /missingness imbalance/i.test(warning)));
}

// Survival prediction: nested CV must produce out-of-sample Harrell C-index.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,survival_time,survival_event'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaSignal = ['SURV_GENE'];
  const rnaNoise = ['SURV_RNA_NOISE'];
  const protSignal = ['SURV_PROTEIN'];
  const protNoise = ['SURV_PROT_NOISE'];
  let assay = 1;

  for (let s = 1; s <= 30; s += 1) {
    const subject = 'SV' + String(s).padStart(2,'0');
    const risk = (s - 15.5) / 8;
    const event = s % 5 === 0 ? 0 : 1;
    const survivalTime = Math.max(1, 40 - 8 * risk + (s % 3));
    for (const layer of ['transcriptomics','proteomics']) {
      const id = (layer === 'transcriptomics' ? 'SR' : 'SP') + assay++;
      rows.push([subject,subject,id,layer,'cohort','T0',survivalTime.toFixed(2),event].join(','));
      if (layer === 'transcriptomics') {
        rnaHeader.push(id);
        rnaSignal.push(String(5 + risk));
        rnaNoise.push(String(3 + ((s * 7) % 11) / 10));
      } else {
        protHeader.push(id);
        protSignal.push(String(2 + 0.8 * risk));
        protNoise.push(String(4 + ((s * 13) % 17) / 20));
      }
    }
  }

  const meta = csvFile('survival_metadata.csv', rows.join('\n'));
  const rna = csvFile('survival_rna.csv', [rnaHeader.join(','),rnaSignal.join(','),rnaNoise.join(',')].join('\n'));
  const protein = csvFile('survival_protein.csv', [protHeader.join(','),protSignal.join(','),protNoise.join(',')].join('\n'));
  const parsedMeta = parseDelimited(await meta.text());

  const survival = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:parsedMeta.rows,
    columnMapping:{...mapping,survival_time:'survival_time',survival_event:'survival_event'},
    protocol:{
      organism:'human',objective:'outcome',outcomeType:'survival',outcomeTimepoint:'T0',
      longitudinal:false,designType:'independent',studySetting:'synthetic_test',
      groupCount:'1',sampleOverlap:'same_specimen',batchKnown:'no',covariateColumns:[]
    },
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });

  assert.equal(survival.layers.transcriptomics.inferenceStatus,
    'reference_R_survival_and_ph_diagnostics_required');
  assert.ok(survival.layers.transcriptomics.rows.every(x=>
    x.pValue===null && x.qValue===null && x.ciLow===null && x.ciHigh===null));
  assert.equal(survival.predictiveOutcome.status, 'ok');
  assert.equal(survival.predictiveOutcome.preprocessingLeakageRisk, false);
  assert.equal(survival.predictiveOutcome.predictions.length, survival.predictiveOutcome.subjects, 'CV cannot silently discard failed outer folds');
  assert.equal(new Set(survival.predictiveOutcome.predictions.map((x)=>x.subjectId)).size, survival.predictiveOutcome.subjects, 'every subject has one held-out prediction');
  assert.equal(survival.predictiveOutcome.outcomeType, 'survival');
  assert.ok(survival.predictiveOutcome.predictions.length >= 24);
  assert.ok(Number.isFinite(survival.predictiveOutcome.metrics.cIndex));
  assert.ok(survival.predictiveOutcome.metrics.cIndex > 0.75);
  assert.ok(survival.predictiveOutcome.foldSummaries.every(x => x.testSubjects > 0 && x.trainingSubjects > x.testSubjects));
  assert.equal(survival.predictiveOutcome.nuisanceAdjustment.policy, 'fold-local');

  // Two assays for the same subject must never overwrite conflicting endpoints.
  await assert.rejects(() => runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:parsedMeta.rows.map((row)=>
      row.subject_id === 'SV01' && row.omic === 'proteomics'
        ? {...row,survival_time:'999'} : row),
    columnMapping:{...mapping,survival_time:'survival_time',survival_event:'survival_event'},
    protocol:{organism:'human',objective:'outcome',outcomeType:'survival',outcomeTimepoint:'T0',
      longitudinal:false,designType:'independent',studySetting:'synthetic_test',
      groupCount:'1',sampleOverlap:'same_specimen',batchKnown:'no',covariateColumns:[]},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,resolveIdentifiers:false
  }), /conflicting outcome for subject/i);

}

// Predictive nuisance handling must be estimated inside each training fold.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,outcome,age'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaSignal = ['PRED_SIGNAL'];
  const rnaBatch = ['PRED_BATCH'];
  const protSignal = ['PRED_PROT'];
  const protAge = ['PRED_AGE'];
  let assay = 1;

  for (let s = 1; s <= 24; s += 1) {
    const subject = 'PV' + String(s).padStart(2,'0');
    const outcome = s <= 12 ? 'A' : 'B';
    const batch = s % 2 ? 'B1' : 'B2';
    const age = 40 + s;
    for (const layer of ['transcriptomics','proteomics']) {
      const id = (layer === 'transcriptomics' ? 'PVR' : 'PVP') + assay++;
      rows.push([subject,subject,id,layer,'cohort','T0',batch,outcome,age].join(','));
      if (layer === 'transcriptomics') {
        rnaHeader.push(id);
        rnaSignal.push(String(3 + (outcome === 'B' ? 1.8 : 0) + (s % 3) * 0.03));
        rnaBatch.push(String(batch === 'B2' ? 8 : 2));
      } else {
        protHeader.push(id);
        protSignal.push(String(4 + (outcome === 'B' ? 1.1 : 0) + (s % 4) * 0.02));
        protAge.push(String(age / 10));
      }
    }
  }

  const meta = csvFile('predictive_nuisance_metadata.csv', rows.join('\n'));
  const rna = csvFile('predictive_nuisance_rna.csv', [rnaHeader.join(','),rnaSignal.join(','),rnaBatch.join(',')].join('\n'));
  const protein = csvFile('predictive_nuisance_protein.csv', [protHeader.join(','),protSignal.join(','),protAge.join(',')].join('\n'));
  const parsedMeta = parseDelimited(await meta.text());
  const pred = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:parsedMeta.rows,
    columnMapping:{...mapping,outcome:'outcome'},
    protocol:{
      organism:'human',objective:'outcome',outcomeType:'binary',outcomeTimepoint:'T0',
      longitudinal:false,designType:'independent',studySetting:'synthetic_test',
      groupCount:'1',sampleOverlap:'same_specimen',batchKnown:'yes',covariateColumns:['age']
    },
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(pred.predictiveOutcome.status, 'ok');
  assert.equal(pred.predictiveOutcome.predictions.length, pred.predictiveOutcome.subjects, 'CV cannot silently discard failed outer folds');
  assert.equal(new Set(pred.predictiveOutcome.predictions.map((x)=>x.subjectId)).size, pred.predictiveOutcome.subjects, 'every subject has one held-out prediction');
  assert.equal(pred.predictiveOutcome.nuisanceAdjustment.policy, 'fold-local');
  assert.equal(pred.predictiveOutcome.preprocessingLeakageRisk, false, 'log-scale predictors should bypass whole-cohort feature filtering');
  assert.deepEqual(pred.predictiveOutcome.nuisanceAdjustment.covariates, ['age']);
  assert.ok(Number.isFinite(pred.predictiveOutcome.metrics.auc));
  assert.ok(pred.predictiveOutcome.metrics.auc > 0.80);

  const inconsistentRows = parsedMeta.rows.map((row)=>
    row.subject_id === 'PV01' && row.omic === 'proteomics'
      ? {...row,outcome:'B'} : row);
  await assert.rejects(() => runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:inconsistentRows,columnMapping:{...mapping,outcome:'outcome'},
    protocol:{organism:'human',objective:'outcome',outcomeType:'binary',outcomeTimepoint:'T0',
      longitudinal:false,designType:'independent',studySetting:'synthetic_test',
      groupCount:'1',sampleOverlap:'same_specimen',batchKnown:'yes',covariateColumns:['age']},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,resolveIdentifiers:false
  }), /conflicting outcome for subject/i);

}


// One-layer integer-count RNA must also avoid outcome-blind whole-cohort QC
// before nested CV; library-size normalization is per assay, never learned
// from validation subjects.
{
  const ids = Array.from({length:30},(_,i)=>'RC'+String(i+1).padStart(2,'0'));
  const rows = ids.map((id,i)=>({
    subject_id:id,sample_id:id,assay_id:id,omic:'transcriptomics',
    condition:'cohort',timepoint:'T0',batch:'',
    outcome:i<15?'control':'treated'
  }));
  const genes = Array.from({length:18},(_,j)=>'COUNT_GENE_'+String(j+1).padStart(2,'0'));
  const matrix = [
    ['feature_id',...ids].join(','),
    ...genes.map((gene,j)=>[
      gene,...ids.map((id,i)=>String(60+j*4+((i*7+j*13)%17) +
        (j<2 && i>=15 ? 135 : 0)))
    ].join(','))
  ].join('\n');
  const rna = new File([matrix],'rna_integer_counts.csv',{type:'text/csv'});
  const meta = new File(['fixture'],'count_metadata.csv',{type:'text/csv'});
  const pred = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:null,metabolomics:null},
    metadataRows:rows,
    columnMapping:{...mapping,outcome:'outcome'},
    protocol:{organism:'human',objective:'outcome',outcomeType:'binary',
      designType:'independent',longitudinal:false,studySetting:'synthetic_test',
      outcomeTimepoint:'T0',groupCount:'1',batchKnown:'no',covariateColumns:[]},
    dataTypes:{transcriptomics:'raw_counts'},
    useReactome:false,resolveIdentifiers:false
  });
  assert.equal(pred.predictiveOutcome.status,'ok');
  assert.equal(pred.predictiveOutcome.preprocessingLeakageRisk,false,
    'integer-count per-assay CPM must bypass cohort-wide feature filtering');
  assert.equal(pred.predictiveOutcome.predictions.length,30);
  assert.equal(pred.engine.analysisMode,'single_omic');
}

console.log('multiomics deterministic engine: PASS');
console.log(JSON.stringify({
  subjects: result.metadataSummary.subjects,
  assays: result.metadataSummary.assays,
  rnaTop: result.layers.transcriptomics.rows.slice(0,3).map((x) => [x.feature, x.foldRatio]),
  proteinTop: result.layers.proteomics.rows.slice(0,3).map((x) => [x.feature, x.foldRatio]),
  metaboliteTop: result.layers.metabolomics.rows.slice(0,3).map((x) => [x.feature, x.foldRatio])
}, null, 2));


function csvFile(name, text) {
  return new File([text], name, { type: 'text/csv' });
}

// Paired design: the same subjects are measured in both biological conditions.
{
  const subjects = ['S1','S2','S3','S4'];
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaValues = ['PAIR_GENE'];
  const protValues = ['PAIR_PROTEIN'];
  let assay = 1;
  for (const subject of subjects) {
    for (const condition of ['A','B']) {
      for (const layer of ['transcriptomics','proteomics']) {
        const sample = `${subject}_${condition}`;
        const id = `${layer === 'transcriptomics' ? 'R' : 'P'}${assay++}`;
        rows.push([subject,sample,id,layer,condition,''].join(','));
        if (layer === 'transcriptomics') {
          rnaHeader.push(id);
          rnaValues.push(String(condition === 'A' ? Number(subject.slice(1)) : Number(subject.slice(1)) + 4));
        } else {
          protHeader.push(id);
          protValues.push(String(condition === 'A' ? Number(subject.slice(1))*2 : Number(subject.slice(1))*2 + 3));
        }
      }
    }
  }
  const meta = csvFile('paired_metadata.csv', rows.join('\n'));
  const rna = csvFile('paired_rna.csv', [rnaHeader.join(','), rnaValues.join(',')].join('\n'));
  const protein = csvFile('paired_protein.csv', [protHeader.join(','), protValues.join(',')].join('\n'));
  const metaParsed = parseDelimited(await meta.text());
  const paired = await runDeterministicAnalysis({
    files: { metadata: meta, transcriptomics: rna, proteomics: protein, metabolomics: null },
    metadataRows: metaParsed.rows,
    columnMapping: mapping,
    protocol: { organism:'human', objective:'groups', longitudinal:false, designType:'paired', studySetting:'clinical_interventional', groupCount:'2', sampleOverlap:'same_specimen' },
    dataTypes: { transcriptomics:'log_expression', proteomics:'log_intensity', metabolomics:'concentration' },
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(paired.layers.transcriptomics.mode, 'paired-permutation');
  assert.deepEqual(paired.layers.transcriptomics.groupSizes, [4,4]);
  assert.ok(paired.layers.transcriptomics.rows[0].effect > 0);
  assert.equal(paired.crossOmics.testedPairs, 0);
}

// Longitudinal design with 3 time points: all times contribute through an individual slope.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint'];
  const headers = { transcriptomics:['feature_id'], proteomics:['feature_id'] };
  const values = { transcriptomics:['SLOPE_GENE'], proteomics:['SLOPE_PROTEIN'] };
  let assay = 1;
  for (const condition of ['control','treatment']) {
    for (let s = 1; s <= 4; s += 1) {
      const subject = `${condition[0].toUpperCase()}${s}`;
      for (const time of [0,6,12]) {
        for (const layer of ['transcriptomics','proteomics']) {
          const sample = `${subject}_T${time}`;
          const id = `${layer === 'transcriptomics' ? 'R' : 'P'}${assay++}`;
          rows.push([subject,sample,id,layer,condition,`T${time}`].join(','));
          headers[layer].push(id);
          const base = s + (layer === 'proteomics' ? 2 : 0);
          const slope = condition === 'treatment' ? 0.20 : 0.02;
          values[layer].push(String(base + slope*time));
        }
      }
    }
  }
  const meta = csvFile('slope_metadata.csv', rows.join('\n'));
  const rna = csvFile('slope_rna.csv', [headers.transcriptomics.join(','), values.transcriptomics.join(',')].join('\n'));
  const protein = csvFile('slope_protein.csv', [headers.proteomics.join(','), values.proteomics.join(',')].join('\n'));
  const metaParsed = parseDelimited(await meta.text());
  const longitudinalResult = await runDeterministicAnalysis({
    files: { metadata: meta, transcriptomics: rna, proteomics: protein, metabolomics: null },
    metadataRows: metaParsed.rows,
    columnMapping: mapping,
    protocol: { organism:'human', objective:'time', longitudinal:true, designType:'repeated', studySetting:'clinical_interventional', groupCount:'2', sampleOverlap:'same_specimen' },
    dataTypes: { transcriptomics:'log_expression', proteomics:'log_intensity', metabolomics:'concentration' },
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(longitudinalResult.layers.transcriptomics.mode, 'random-intercept-longitudinal-model');
  assert.deepEqual(longitudinalResult.layers.transcriptomics.groupSizes, [4,4]);
  assert.ok(longitudinalResult.layers.transcriptomics.rows[0].effect > 0.15);
}

// Multi-group longitudinal designs use an omnibus condition × time interaction.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint'];
  const headers = { transcriptomics:['feature_id'], proteomics:['feature_id'] };
  const values = { transcriptomics:['MULTI_SLOPE_GENE'], proteomics:['MULTI_SLOPE_PROTEIN'] };
  let assay = 1;
  const slopeByCondition = { control: 0.02, treatmentA: 0.18, treatmentB: -0.10 };

  for (const condition of Object.keys(slopeByCondition)) {
    for (let s = 1; s <= 5; s += 1) {
      const subject = (condition === 'control' ? 'C' : condition === 'treatmentA' ? 'A' : 'B') + s;
      for (const time of [0,6,12]) {
        for (const layer of ['transcriptomics','proteomics']) {
          const sample = `${subject}_T${time}`;
          const id = `${layer === 'transcriptomics' ? 'MR' : 'MP'}${assay++}`;
          rows.push([subject,sample,id,layer,condition,`T${time}`].join(','));
          headers[layer].push(id);
          const base = 5 + s * 0.05 + (layer === 'proteomics' ? 2 : 0);
          values[layer].push(String(base + slopeByCondition[condition] * time));
        }
      }
    }
  }

  const meta = csvFile('multigroup_longitudinal_metadata.csv', rows.join('\n'));
  const rna = csvFile('multigroup_longitudinal_rna.csv', [headers.transcriptomics.join(','), values.transcriptomics.join(',')].join('\n'));
  const protein = csvFile('multigroup_longitudinal_protein.csv', [headers.proteomics.join(','), values.proteomics.join(',')].join('\n'));
  const metaParsed = parseDelimited(await meta.text());

  const multiLong = await runDeterministicAnalysis({
    files: { metadata: meta, transcriptomics: rna, proteomics: protein, metabolomics: null },
    metadataRows: metaParsed.rows,
    columnMapping: mapping,
    protocol: {
      organism:'human', objective:'time', longitudinal:true, designType:'repeated',
      studySetting:'synthetic_test', groupCount:'3', sampleOverlap:'same_specimen'
    },
    dataTypes: { transcriptomics:'log_expression', proteomics:'log_intensity', metabolomics:'concentration' },
    useReactome:false,
    resolveIdentifiers:false
  });

  assert.equal(multiLong.layers.transcriptomics.mode, 'random-intercept-longitudinal-omnibus');
  assert.deepEqual(multiLong.layers.transcriptomics.groupSizes, [5,5,5]);
  assert.match(multiLong.layers.transcriptomics.contrast, /omnibus condition × time interaction across 3 groups/);
  const row = multiLong.layers.transcriptomics.rows.find((item) => item.feature === 'MULTI_SLOPE_GENE');
  assert.ok(row);
  assert.equal(row.statistic, null);
  assert.equal(row.pValue, null);
  assert.equal(row.qValue, null);
  assert.equal(row.inferentialStatus, 'uncalibrated_reference_R_required');
  assert.deepEqual(Object.keys(row.groupSlopes), ['control','treatmentA','treatmentB']);
  assert.ok(row.effect > 2);
}

// High-dimensional branch: >500 features must switch from permutations to analytic inference.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaRows = [];
  const protRows = [];
  const rnaAssays = [];
  const protAssays = [];
  let assay = 1;
  for (const condition of ['control','treatment']) {
    for (let s = 1; s <= 8; s += 1) {
      const subject = `${condition[0].toUpperCase()}${s}`;
      for (const layer of ['transcriptomics','proteomics']) {
        const id = `${layer === 'transcriptomics' ? 'HR' : 'HP'}${assay++}`;
        rows.push([subject,subject,id,layer,condition,''].join(','));
        if (layer === 'transcriptomics') rnaAssays.push({id,condition,s});
        else protAssays.push({id,condition,s});
      }
    }
  }
  rnaHeader.push(...rnaAssays.map(x=>x.id));
  protHeader.push(...protAssays.map(x=>x.id));
  for (let feature = 0; feature < 501; feature += 1) {
    const name = feature === 0 ? 'HIGH_DIM_SIGNAL' : `GENE_${feature}`;
    const values = rnaAssays.map(({condition,s}) => {
      const baseline = 5 + s*0.03 + (feature % 7)*0.01;
      return String(baseline + (feature === 0 && condition === 'treatment' ? 2.5 : 0));
    });
    rnaRows.push([name,...values].join(','));
  }
  for (let feature = 0; feature < 10; feature += 1) {
    const name = feature === 0 ? 'PROT_SIGNAL' : `PROT_${feature}`;
    const values = protAssays.map(({condition,s}) => {
      const baseline = 3 + s*0.02 + feature*0.01;
      return String(baseline + (feature === 0 && condition === 'treatment' ? 1.2 : 0));
    });
    protRows.push([name,...values].join(','));
  }

  const meta = csvFile('highdim_metadata.csv', rows.join('\n'));
  const rna = csvFile('highdim_rna.csv', [rnaHeader.join(','),...rnaRows].join('\n'));
  const protein = csvFile('highdim_protein.csv', [protHeader.join(','),...protRows].join('\n'));
  const metaParsed = parseDelimited(await meta.text());
  const highdim = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:metaParsed.rows,
    columnMapping:mapping,
    protocol:{organism:'human',objective:'groups',longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen'},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(highdim.layers.transcriptomics.inferenceMethod, 'OLS with HC3 robust standard errors');
  assert.equal(highdim.layers.proteomics.inferenceMethod, 'OLS with HC3 robust standard errors');
  const signal = highdim.layers.transcriptomics.rows.find(x=>x.feature === 'HIGH_DIM_SIGNAL');
  assert.ok(signal && signal.effect > 2 && signal.pValue < 1e-6);
}

// ChEBI mapping is conservative: explicit lipid alias is queried, but only exact returned names are accepted.
{
  const mockFetch = async () => new Response(JSON.stringify({
    results: [{ _source: { chebi_accession:'CHEBI:17351', name:'linoleic acid', ascii_name:'linoleic acid', stars:3 } }],
    total:1,
    number_pages:1
  }), { status:200, headers:{'Content-Type':'application/json'} });
  const mappingResult = await resolveMetaboliteIdentifier('C18.2n.6', { fetchFn: mockFetch });
  assert.equal(mappingResult.resolved, 'CHEBI:17351');
  assert.equal(mappingResult.method, 'chebi_exact_after_alias');

  const wrongFetch = async () => new Response(JSON.stringify({
    results: [{ _source: { chebi_accession:'CHEBI:99999', name:'malic enzyme substrate', ascii_name:'malic enzyme substrate' } }],
    total:1,
    number_pages:1
  }), { status:200, headers:{'Content-Type':'application/json'} });
  const unresolved = await resolveMetaboliteIdentifier('definitely-not-an-exact-metabolite', { fetchFn: wrongFetch });
  assert.equal(unresolved.resolved, null);
  assert.equal(unresolved.status, 'unresolved');

  const unichemFetch = async (url) => {
    assert.match(String(url), /unichem\/rest\/verbose_inchikey/);
    return new Response(JSON.stringify([
      { name:'chebi', src_id:7, src_compound_id:['CHEBI:34967'] },
      { name:'hmdb', src_id:18, src_compound_id:['HMDB0015623'] }
    ]), { status:200, headers:{'Content-Type':'application/json'} });
  };
  const unichem = await resolveMetaboliteIdentifier('JZCPYUJPEARBJL-UHFFFAOYSA-N', { fetchFn: unichemFetch });
  assert.equal(unichem.resolved, 'CHEBI:34967');
  assert.equal(unichem.method, 'unichem_inchikey_to_chebi');
}

// Gene and protein identifier mapping remains conservative and testable with mocked official APIs.
{
  const geneFetch = async (url) => {
    assert.match(String(url), /rest\.ensembl\.org\/lookup\/symbol\/homo_sapiens\/TP53/);
    return new Response(JSON.stringify({ id:'ENSG00000141510', display_name:'TP53', biotype:'protein_coding' }), {
      status:200,
      headers:{'Content-Type':'application/json'}
    });
  };
  const gene = await resolveGeneIdentifier('TP53', { identifierType:'gene_symbol', organism:'human', fetchFn:geneFetch });
  assert.equal(gene.resolved, 'ENSG00000141510');
  assert.equal(gene.method, 'ensembl_lookup_symbol');

  const proteinFetch = async (url) => {
    assert.match(String(url), /rest\.uniprot\.org\/uniprotkb\/search/);
    return new Response(JSON.stringify({
      results:[{ primaryAccession:'P04637', entryType:'UniProtKB reviewed (Swiss-Prot)' }]
    }), { status:200, headers:{'Content-Type':'application/json'} });
  };
  const protein = await resolveProteinIdentifier('TP53', { identifierType:'gene_symbol', organism:'human', fetchFn:proteinFetch });
  assert.equal(protein.resolved, 'P04637');
  assert.equal(protein.method, 'uniprot_gene_search');
}

// Batch safety: a biological condition completely nested inside technical batch must be refused.
{
  const metadataText = [
    'subject_id,sample_id,assay_id,omic,condition,timepoint,batch',
    'C1,C1,RC1,transcriptomics,control,,BATCH_CONTROL',
    'C1,C1,PC1,proteomics,control,,BATCH_CONTROL',
    'C2,C2,RC2,transcriptomics,control,,BATCH_CONTROL',
    'C2,C2,PC2,proteomics,control,,BATCH_CONTROL',
    'T1,T1,RT1,transcriptomics,treatment,,BATCH_TREATMENT',
    'T1,T1,PT1,proteomics,treatment,,BATCH_TREATMENT',
    'T2,T2,RT2,transcriptomics,treatment,,BATCH_TREATMENT',
    'T2,T2,PT2,proteomics,treatment,,BATCH_TREATMENT'
  ].join('\n');
  const meta = csvFile('batch_confounded_metadata.csv', metadataText);
  const rna = csvFile('batch_confounded_rna.csv', [
    'feature_id,RC1,RC2,RT1,RT2',
    'GENE_A,1,1.1,5,5.1'
  ].join('\n'));
  const protein = csvFile('batch_confounded_protein.csv', [
    'feature_id,PC1,PC2,PT1,PT2',
    'PROT_A,2,2.1,6,6.1'
  ].join('\n'));
  const metaParsed = parseDelimited(await meta.text());

  await assert.rejects(
    () => runDeterministicAnalysis({
      files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
      metadataRows:metaParsed.rows,
      columnMapping:mapping,
      protocol:{organism:'human',objective:'groups',longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen'},
      dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
      useReactome:false,
      resolveIdentifiers:false
    }),
    /Technical batch confounding prevents identifiable biological inference/
  );
}


// Balanced batch adjustment: batch effect must be removed while preserving the biological contrast.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,outcome,age'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaValues = ['BATCH_GENE'];
  const protValues = ['BATCH_PROTEIN'];
  let assay = 1;
  for (const condition of ['control','treatment']) {
    for (const batch of ['B1','B2']) {
      for (let s = 1; s <= 4; s += 1) {
        const subject = condition[0].toUpperCase() + batch + s;
        const age = 30 + s;
        for (const layer of ['transcriptomics','proteomics']) {
          const id = (layer === 'transcriptomics' ? 'BR' : 'BP') + assay++;
          rows.push([subject,subject,id,layer,condition,'',batch,'1','',age].join(','));
          const biological = condition === 'treatment' ? 2 : 0;
          const batchEffect = batch === 'B2' ? 8 : 0;
          const value = 5 + biological + batchEffect + s * 0.01;
          if (layer === 'transcriptomics') {
            rnaHeader.push(id);
            rnaValues.push(String(value));
          } else {
            protHeader.push(id);
            protValues.push(String(value + 1));
          }
        }
      }
    }
  }
  const meta = csvFile('balanced_batch_metadata.csv', rows.join('\n'));
  const rna = csvFile('balanced_batch_rna.csv', [rnaHeader.join(','),rnaValues.join(',')].join('\n'));
  const protein = csvFile('balanced_batch_protein.csv', [protHeader.join(','),protValues.join(',')].join('\n'));
  const metaParsed = parseDelimited(await meta.text());
  const adjusted = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:metaParsed.rows,
    columnMapping:mapping,
    protocol:{organism:'human',objective:'groups',longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen',batchKnown:'yes',covariateColumns:['age']},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(adjusted.metadataSummary.batchAudit.transcriptomics.status, 'multiple_batches_adjusted');
  assert.equal(adjusted.layers.transcriptomics.adjustment.applied, true);
  assert.ok(Math.abs(adjusted.layers.transcriptomics.rows[0].effect - 2) < 0.1);
}

// True unsupervised multi-block branch: shared latent structure across two omics.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaSignal = ['RNA_SHARED'];
  const rnaNoise = ['RNA_NOISE'];
  const protSignal = ['PROT_SHARED'];
  const protNoise = ['PROT_NOISE'];
  let assay = 1;
  for (let s = 1; s <= 8; s += 1) {
    const subject = 'E' + s;
    const latent = s - 4.5;
    for (const layer of ['transcriptomics','proteomics']) {
      const id = (layer === 'transcriptomics' ? 'ER' : 'EP') + assay++;
      rows.push([subject,subject,id,layer,'','',''].join(','));
      if (layer === 'transcriptomics') {
        rnaHeader.push(id);
        rnaSignal.push(String(latent * 2));
        rnaNoise.push(String((s % 3) * 0.1));
      } else {
        protHeader.push(id);
        protSignal.push(String(latent * 1.5));
        protNoise.push(String(((s + 1) % 4) * 0.1));
      }
    }
  }
  const meta = csvFile('explore_metadata.csv', rows.join('\n'));
  const rna = csvFile('explore_rna.csv', [rnaHeader.join(','),rnaSignal.join(','),rnaNoise.join(',')].join('\n'));
  const protein = csvFile('explore_protein.csv', [protHeader.join(','),protSignal.join(','),protNoise.join(',')].join('\n'));
  const metaParsed = parseDelimited(await meta.text());
  const explored = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:metaParsed.rows,
    columnMapping:mapping,
    protocol:{organism:'human',objective:'explore',longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'1',sampleOverlap:'same_specimen',batchKnown:'no',covariateColumns:[]},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(explored.exploration.subjects, 8);
  assert.ok(explored.exploration.components.length >= 1);
  const loadingNames = explored.exploration.components[0].topLoadings.map(x=>x.feature);
  assert.ok(loadingNames.includes('RNA_SHARED'));
  assert.ok(loadingNames.includes('PROT_SHARED'));
  assert.equal(explored.layers.transcriptomics.mode, 'exploratory-multiblock-pca');
}

// Outcome branch: continuous, binary, count and survival models with explicit covariate adjustment.
{
  const makeOutcomeDataset = (type) => {
    const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,outcome,age,survival_time,survival_event'];
    const rnaHeader = ['feature_id'];
    const protHeader = ['feature_id'];
    const rnaSignal = ['OUTCOME_GENE'];
    const rnaNoise = ['OUTCOME_NOISE'];
    const protSignal = ['OUTCOME_PROTEIN'];
    const protNoise = ['PROT_NOISE'];
    let assay = 1;
    for (let s = 1; s <= 24; s += 1) {
      const subject = 'O' + s;
      const age = 40 + (s % 8);
      const latent = (s - 12.5) / 4;
      const batch = s % 2 === 0 ? 'B2' : 'B1';
      const batchEffect = batch === 'B2' ? 1.5 : 0;
      let outcome = '';
      let survivalTime = '';
      let survivalEvent = '';
      if (type === 'continuous') outcome = String(2 * latent + 0.05 * age);
      if (type === 'binary') outcome = latent > 0 ? 'responder' : 'nonresponder';
      if (type === 'multiclass') outcome = s <= 8 ? 'class_A' : s <= 16 ? 'class_B' : 'class_C';
      if (type === 'count') outcome = String(Math.max(0, Math.round(4 + latent * 1.5)));
      if (type === 'survival') {
        survivalTime = String(Math.max(1, 30 - latent * 4 + (s % 3)));
        survivalEvent = s % 5 === 0 ? '0' : '1';
      }
      for (const layer of ['transcriptomics','proteomics']) {
        const id = (layer === 'transcriptomics' ? 'OR' : 'OP') + assay++;
        rows.push([subject,subject,id,layer,'','',batch,'1',outcome,age,survivalTime,survivalEvent].join(','));
        const signalValue = latent + 0.02 * age + batchEffect;
        if (layer === 'transcriptomics') {
          rnaHeader.push(id);
          rnaSignal.push(String(signalValue));
          rnaNoise.push(String((s % 5) * 0.03));
        } else {
          protHeader.push(id);
          protSignal.push(String(signalValue * 0.8));
          protNoise.push(String(((s + 2) % 6) * 0.02));
        }
      }
    }
    return {
      meta: csvFile('outcome_' + type + '_metadata.csv', rows.join('\n')),
      rna: csvFile('outcome_' + type + '_rna.csv', [rnaHeader.join(','),rnaSignal.join(','),rnaNoise.join(',')].join('\n')),
      protein: csvFile('outcome_' + type + '_protein.csv', [protHeader.join(','),protSignal.join(','),protNoise.join(',')].join('\n'))
    };
  };

  for (const type of ['continuous','binary','multiclass','count','survival']) {
    const ds = makeOutcomeDataset(type);
    const parsedMeta = parseDelimited(await ds.meta.text());
    const outcomeResult = await runDeterministicAnalysis({
      files:{metadata:ds.meta,transcriptomics:ds.rna,proteomics:ds.protein,metabolomics:null},
      metadataRows:parsedMeta.rows,
      columnMapping:mapping,
      protocol:{organism:'human',objective:'outcome',outcomeType:type,longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'1',sampleOverlap:'same_specimen',batchKnown:'yes',covariateColumns:['age']},
      dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
      useReactome:false,
      resolveIdentifiers:false
    });
    assert.equal(outcomeResult.layers.transcriptomics.mode, 'outcome-' + type);
    assert.match(outcomeResult.layers.transcriptomics.adjustment.method, /direct nuisance adjustment/);
    assert.ok(outcomeResult.layers.transcriptomics.adjustment.columns.some(x=>x.startsWith('batch=')));
    assert.ok(outcomeResult.layers.transcriptomics.adjustment.columns.includes('age'));
    assert.equal(outcomeResult.predictiveOutcome.status, 'ok');
    assert.equal(outcomeResult.predictiveOutcome.preprocessingLeakageRisk, false, 'as-supplied log-scale predictive route should be fold-isolated');
    assert.equal(outcomeResult.predictiveOutcome.predictions.length, outcomeResult.predictiveOutcome.subjects, 'CV cannot silently discard failed outer folds');
    assert.equal(new Set(outcomeResult.predictiveOutcome.predictions.map((x)=>x.subjectId)).size, outcomeResult.predictiveOutcome.subjects, 'every subject has one held-out prediction');
    assert.ok(outcomeResult.predictiveOutcome.predictions.length >= 12);
    assert.ok(outcomeResult.predictiveOutcome.foldSummaries.every(x=>x.testSubjects > 0 && x.trainingSubjects > x.testSubjects));
    if (type === 'survival') {
      assert.ok(Number.isFinite(outcomeResult.predictiveOutcome.metrics.cIndex));
      assert.ok(outcomeResult.predictiveOutcome.metrics.events >= 6);
    }
    const signal = outcomeResult.layers.transcriptomics.rows.find(x=>x.feature === 'OUTCOME_GENE');
    assert.ok(signal && Number.isFinite(signal.effect));
    if (type === 'survival') {
      assert.equal(signal.pValue,null,'exploratory Cox must not report unchecked p-values');
      assert.equal(signal.qValue,null);
      assert.equal(signal.ciLow,null);
      assert.equal(signal.ciHigh,null);
      assert.equal(signal.inferentialStatus,'exploratory_cox_reference_R_and_PH_diagnostics_required');
      assert.equal(signal.tieMethod,'breslow');
      assert.equal(outcomeResult.layers.transcriptomics.inferenceStatus,'reference_R_survival_and_ph_diagnostics_required');
    } else {
      assert.ok(Number.isFinite(signal.pValue));
    }
  }
}

// QC: low-information features are filtered and QC PCA is produced.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaSignal = ['RNA_SIGNAL'];
  const rnaSecond = ['RNA_SECOND'];
  const rnaLow = ['RNA_LOW'];
  const protSignal = ['PROT_SIGNAL'];
  const protMissing = ['PROT_MISSING'];
  let assay = 1;
  for (const condition of ['control','treatment']) {
    for (let s = 1; s <= 5; s += 1) {
      const subject = condition[0].toUpperCase() + 'Q' + s;
      for (const layer of ['transcriptomics','proteomics']) {
        const id = (layer === 'transcriptomics' ? 'QR' : 'QP') + assay++;
        rows.push([subject,subject,id,layer,condition,'','B1'].join(','));
        if (layer === 'transcriptomics') {
          rnaHeader.push(id);
          rnaSignal.push(String(100 + (condition === 'treatment' ? 100 : 0) + s));
          rnaSecond.push(String(70 + (condition === 'treatment' ? 15 : 0) + 4 * s));
          rnaLow.push(s === 1 && condition === 'control' ? '1' : '0');
        } else {
          protHeader.push(id);
          protSignal.push(String(5 + (condition === 'treatment' ? 2 : 0) + s * 0.1));
          protMissing.push(s <= 2 && condition === 'control' ? String(2+s) : '');
        }
      }
    }
  }
  const meta = csvFile('qc_metadata.csv', rows.join('\n'));
  const rna = csvFile('qc_rna.csv', [rnaHeader.join(','),rnaSignal.join(','),rnaSecond.join(','),rnaLow.join(',')].join('\n'));
  const protein = csvFile('qc_protein.csv', [protHeader.join(','),protSignal.join(','),protMissing.join(',')].join('\n'));
  const metaParsed = parseDelimited(await meta.text());
  const qcResult = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:metaParsed.rows,
    columnMapping:mapping,
    protocol:{organism:'human',objective:'groups',longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen',batchKnown:'yes',covariateColumns:[]},
    dataTypes:{transcriptomics:'raw_counts',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.ok(qcResult.layers.transcriptomics.qc.featuresAfter < qcResult.layers.transcriptomics.qc.featuresBefore);
  assert.ok(qcResult.layers.transcriptomics.qc.pca.scores.length >= 3);
  assert.equal(qcResult.layers.transcriptomics.mode, 'adjusted-two-group-model');
  assert.match(qcResult.layers.transcriptomics.inferenceMethod, /HC3/);
}

// Longitudinal repeated-measures branch: random intercept and condition × time interaction.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rnaSignal = ['TIME_GENE'];
  const protSignal = ['TIME_PROTEIN'];
  let assay = 1;
  for (const condition of ['control','treatment']) {
    for (let s = 1; s <= 6; s += 1) {
      const subject = condition[0].toUpperCase() + 'L' + s;
      const randomIntercept = s * 0.4;
      for (const time of [0,6,12]) {
        for (const layer of ['transcriptomics','proteomics']) {
          const id = (layer === 'transcriptomics' ? 'LR' : 'LP') + assay++;
          rows.push([subject,subject + '_T' + time,id,layer,condition,'T' + time,'B1'].join(','));
          const slope = condition === 'treatment' ? 0.20 : 0.02;
          const value = randomIntercept + slope * time + (layer === 'proteomics' ? 1 : 0);
          if (layer === 'transcriptomics') { rnaHeader.push(id); rnaSignal.push(String(value)); }
          else { protHeader.push(id); protSignal.push(String(value)); }
        }
      }
    }
  }
  const meta = csvFile('lmm_metadata.csv', rows.join('\n'));
  const rna = csvFile('lmm_rna.csv', [rnaHeader.join(','),rnaSignal.join(',')].join('\n'));
  const protein = csvFile('lmm_protein.csv', [protHeader.join(','),protSignal.join(',')].join('\n'));
  const parsedMeta = parseDelimited(await meta.text());
  const longitudinal = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rna,proteomics:protein,metabolomics:null},
    metadataRows:parsedMeta.rows,
    columnMapping:mapping,
    protocol:{organism:'human',objective:'time',longitudinal:true,designType:'repeated',studySetting:'synthetic_test',groupCount:'2',sampleOverlap:'same_specimen',batchKnown:'yes',covariateColumns:[]},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(longitudinal.layers.transcriptomics.mode, 'random-intercept-longitudinal-model');
  const timeGene = longitudinal.layers.transcriptomics.rows.find((row) => row.feature === 'TIME_GENE');
  assert.ok(timeGene && timeGene.effect > 1);
  assert.ok(Number.isFinite(timeGene.intraclassCorrelation));
}

// Partial-block exploratory integration retains subjects missing one entire layer.
{
  const rows = ['subject_id,sample_id,assay_id,omic,condition,timepoint,batch'];
  const rnaHeader = ['feature_id'];
  const protHeader = ['feature_id'];
  const rna = ['RNA_PARTIAL'];
  const protein = ['PROT_PARTIAL'];
  let assay = 1;
  for (let s = 1; s <= 8; s += 1) {
    const subject = 'M' + s;
    const rnaId = 'MR' + assay++;
    rows.push([subject,subject,rnaId,'transcriptomics','','',''].join(','));
    rnaHeader.push(rnaId); rna.push(String(s));
    if (s <= 6) {
      const protId = 'MP' + assay++;
      rows.push([subject,subject,protId,'proteomics','','',''].join(','));
      protHeader.push(protId); protein.push(String(s * 2));
    }
  }
  const meta = csvFile('partial_metadata.csv', rows.join('\n'));
  const rnaFile = csvFile('partial_rna.csv', [rnaHeader.join(','),rna.join(',')].join('\n'));
  const protFile = csvFile('partial_protein.csv', [protHeader.join(','),protein.join(',')].join('\n'));
  const parsedMeta = parseDelimited(await meta.text());
  const partial = await runDeterministicAnalysis({
    files:{metadata:meta,transcriptomics:rnaFile,proteomics:protFile,metabolomics:null},
    metadataRows:parsedMeta.rows,
    columnMapping:mapping,
    protocol:{organism:'human',objective:'explore',longitudinal:false,designType:'independent',studySetting:'synthetic_test',groupCount:'1',sampleOverlap:'partial',batchKnown:'no',covariateColumns:[],partialOmicsExpected:'yes'},
    dataTypes:{transcriptomics:'log_expression',proteomics:'log_intensity',metabolomics:'concentration'},
    useReactome:false,
    resolveIdentifiers:false
  });
  assert.equal(partial.exploration.subjects, 8);
  assert.equal(partial.exploration.completeSubjects, 6);
  assert.ok(partial.exploration.imputedCells > 0);
  assert.match(partial.exploration.missingDataPolicy, /partial blocks/);
}

console.log('multiomics paired / longitudinal / identifier-resolution branches: PASS');

// Standalone display contract: does not change any inference engine result.
await import('./test_multiomics_visualizations.mjs');
