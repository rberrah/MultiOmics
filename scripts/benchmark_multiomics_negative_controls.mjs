import fs from 'node:fs/promises';
import { File } from 'node:buffer';
import { parseDelimited, runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';

const root = process.argv[2] || 'tmp/public-benchmarks';
const reportPath = `${root}/benchmark-report.json`;

async function asFile(filePath) {
  const text = await fs.readFile(filePath, 'utf8');
  return new File([text], filePath.split('/').pop(), { type: 'text/csv' });
}

function mulberry32(seed) {
  let a = seed >>> 0;
  return () => {
    a |= 0;
    a = (a + 0x6D2B79F5) | 0;
    let t = Math.imul(a ^ (a >>> 15), 1 | a);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function shuffled(values, seed) {
  const random = mulberry32(seed);
  const out = [...values];
  for (let i = out.length - 1; i > 0; i -= 1) {
    const j = Math.floor(random() * (i + 1));
    [out[i], out[j]] = [out[j], out[i]];
  }
  return out;
}

function mean(values) {
  return values.reduce((sum, value) => sum + value, 0) / values.length;
}

function median(values) {
  const sorted = [...values].sort((a, b) => a - b);
  if (!sorted.length) return NaN;
  const mid = Math.floor(sorted.length / 2);
  return sorted.length % 2 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;
}

function requireControl(condition, label, detail = '') {
  if (!condition) throw new Error(`NEGATIVE CONTROL FAILED: ${label}${detail ? ` — ${detail}` : ''}`);
  console.log(`NEGATIVE CONTROL PASS: ${label}${detail ? ` — ${detail}` : ''}`);
}

const dir = `${root}/tcga-her2-luma`;
const metadataFile = await asFile(`${dir}/metadata.csv`);
const parsed = parseDelimited(await metadataFile.text());
const rna = await asFile(`${dir}/transcriptomics.csv`);
const protein = await asFile(`${dir}/proteomics.csv`);

const subjectRows = new Map();
for (const row of parsed.rows) {
  if (!subjectRows.has(row.subject_id)) subjectRows.set(row.subject_id, row);
}
const subjects = [...subjectRows.keys()].sort();
const originalLabels = subjects.map((subject) => subjectRows.get(subject)?.condition || '');
const classCounts = Object.fromEntries([...new Set(originalLabels)].map((label) => [label, originalLabels.filter((x) => x === label).length]));
requireControl(Object.keys(classCounts).length === 2, 'public TCGA negative control uses a binary target', JSON.stringify(classCounts));

const mapping = {
  subject_id: 'subject_id',
  sample_id: 'sample_id',
  assay_id: 'assay_id',
  omic: 'omic',
  condition: 'condition',
  timepoint: 'timepoint',
  batch: 'batch',
  technical_replicate: 'technical_replicate',
  outcome: 'outcome'
};

const seeds = [1103, 2213, 3301, 4421, 5503];
const permutations = [];
for (const seed of seeds) {
  let labels = shuffled(originalLabels, seed);
  // Avoid a chance near-identity or near-inversion permutation because both can
  // preserve genuine subtype signal even though labels were technically shuffled.
  let agreement = labels.filter((label, i) => label === originalLabels[i]).length / labels.length;
  if (agreement < 0.30 || agreement > 0.70) {
    labels = shuffled(originalLabels, seed + 7919);
    agreement = labels.filter((label, i) => label === originalLabels[i]).length / labels.length;
  }

  const bySubject = new Map(subjects.map((subject, i) => [subject, labels[i]]));
  const permutedRows = parsed.rows.map((row) => ({
    ...row,
    outcome: bySubject.get(row.subject_id) || ''
  }));

  const result = await runDeterministicAnalysis({
    files: {
      metadata: metadataFile,
      transcriptomics: rna,
      proteomics: protein,
      metabolomics: null
    },
    metadataRows: permutedRows,
    columnMapping: mapping,
    protocol: {
      organism: 'human',
      objective: 'outcome',
      outcomeType: 'binary',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      sampleOverlap: 'same_specimen',
      batchKnown: 'no',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'normalized'
    },
    identifierTypes: {
      transcriptomics: 'unknown',
      proteomics: 'unknown',
      metabolomics: 'unknown'
    },
    useReactome: false,
    resolveIdentifiers: false
  });

  requireControl(result.predictiveOutcome?.status === 'ok', `permutation ${seed} completes subject-level nested CV (upstream preprocessing remains global)`);
  const auc = Number(result.predictiveOutcome?.metrics?.auc);
  const qValues = Object.values(result.layers || {})
    .flatMap((layer) => layer?.rows || [])
    .map((row) => Number(row.qValue))
    .filter(Number.isFinite);
  const q05 = qValues.filter((q) => q <= 0.05).length;
  const q10 = qValues.filter((q) => q <= 0.10).length;
  permutations.push({ seed, labelAgreementWithTruth: agreement, auc, testedFeatures: qValues.length, q05, q10 });
}

const aucs = permutations.map((item) => item.auc).filter(Number.isFinite);
const q05Fractions = permutations.map((item) => item.testedFeatures ? item.q05 / item.testedFeatures : 0);
const q10Fractions = permutations.map((item) => item.testedFeatures ? item.q10 / item.testedFeatures : 0);
const aucMean = mean(aucs);
const aucMedian = median(aucs);
const q05MedianFraction = median(q05Fractions);
const q10MedianFraction = median(q10Fractions);

requireControl(
  aucs.length === seeds.length && aucMean >= 0.35 && aucMean <= 0.65,
  'permuted TCGA labels collapse mean out-of-sample discrimination toward chance',
  `mean AUC=${aucMean.toFixed(3)}, median AUC=${aucMedian.toFixed(3)}`
);
requireControl(
  q05MedianFraction <= 0.02,
  'permuted labels do not generate a large feature-wise FDR-positive fraction',
  `median q≤0.05 fraction=${(100 * q05MedianFraction).toFixed(2)}%`
);
requireControl(
  q10MedianFraction <= 0.05,
  'exploratory q≤0.10 findings remain sparse under permuted labels',
  `median q≤0.10 fraction=${(100 * q10MedianFraction).toFixed(2)}%`
);

const report = JSON.parse(await fs.readFile(reportPath, 'utf8'));
report.negativeControls = {
  tcgaHer2LumALabelPermutation: {
    purpose: 'Detect information leakage and gross anti-conservative behaviour by destroying the true phenotype-to-omics relationship while preserving matrices and class counts.',
    permutations,
    summary: {
      meanAuc: aucMean,
      medianAuc: aucMedian,
      medianQ05Fraction: q05MedianFraction,
      medianQ10Fraction: q10MedianFraction
    },
    acceptanceCriteria: {
      meanAuc: '0.35 to 0.65',
      medianQ05FeatureFraction: '<= 0.02',
      medianQ10FeatureFraction: '<= 0.05'
    },
    interpretation: 'Passing this control does not prove absence of every possible leakage mechanism, but failure is a hard stop for predictive claims.'
  }
};
await fs.writeFile(reportPath, JSON.stringify(report, null, 2) + '\n');
console.log('PUBLIC MULTI-OMICS NEGATIVE CONTROLS: PASS');
