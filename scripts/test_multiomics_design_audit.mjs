import assert from 'node:assert/strict';
import {
  buildContrastPlan, assessBaselineMismatch, assessReplication,
  variancePartition, assessPValueCalibration, auditFeatureTable,
  detectNonBiologicalSignals, reproducibilityAcrossBlocks, auditDesign
} from '../src/lib/multiomics/design-audit.js';

const codes = (findings) => findings.map((f) => f.code);
const has = (findings, code) => codes(findings).includes(code);
const bilingual = (findings) => findings.every((f) =>
  typeof f.message === 'string' && f.message.length > 0 &&
  typeof f.messageEn === 'string' && f.messageEn.length > 0 &&
  ['error', 'warning', 'info'].includes(f.severity));

/* ---- 1. contrast plan ------------------------------------------------- */

const PLAN_SPEC = {
  arms: ['Untreated', 'Vehicle', 'MicrosomalBlank', 'TAC_0.05', 'TAC_0.5', 'MET_0.05', 'MET_0.5'],
  controlOf: {
    'TAC_0.05': 'Vehicle', 'TAC_0.5': 'Vehicle',
    'MET_0.05': 'MicrosomalBlank', 'MET_0.5': 'MicrosomalBlank'
  },
  untreated: 'Untreated'
};
const plan = buildContrastPlan(PLAN_SPEC);
assert.deepEqual(plan.controls.sort(), ['MicrosomalBlank', 'Vehicle']);
assert.equal(plan.treatments.length, 4);

// every treatment is compared with its OWN control, never with the other one
for (const c of plan.contrasts.filter((x) => x.kind === 'treatment_vs_own_control')) {
  assert.equal(c.reference, PLAN_SPEC.controlOf[c.test], c.id);
}
// the two controls are compared with each other, and each with untreated
assert.ok(plan.contrasts.some((c) => c.kind === 'control_vs_control'));
assert.equal(plan.contrasts.filter((c) => c.kind === 'control_vs_untreated').length, 2);
// arms with different controls get a matrix-corrected head-to-head, not a raw one
const cross = plan.contrasts.find((c) => c.test === 'TAC_0.5' && c.reference === 'MET_0.5')
  || plan.contrasts.find((c) => c.test === 'MET_0.5' && c.reference === 'TAC_0.5');
assert.equal(cross.kind, 'head_to_head_corrected');
assert.match(cross.expression, /\(.+-.+\) - \(.+-.+\)/);
// arms that share a control get a plain difference
const shared = plan.contrasts.find((c) => c.kind === 'head_to_head_shared_control');
assert.ok(shared && !/\(/.test(shared.expression));
assert.ok(bilingual(plan.findings));

// an arm with no declared control is an error, not a silent default
const loose = buildContrastPlan({ arms: ['A', 'B', 'Ctrl'], controlOf: { A: 'Ctrl' } });
assert.ok(has(loose.findings, 'arm_without_control'));
assert.throws(() => buildContrastPlan({ arms: ['only'] }));

/* ---- 2. baseline mismatch --------------------------------------------- */

// the real study: the control-vs-control contrast dwarfs every treatment
const mismatch = assessBaselineMismatch(plan, {
  'MicrosomalBlank__vs__Vehicle': 1381,
  'TAC_0.5__vs__Vehicle': 1,
  'MET_0.5__vs__MicrosomalBlank': 943
});
assert.equal(mismatch.ok, false);
assert.equal(mismatch.code, 'baseline_mismatch');
assert.match(mismatch.findings[0].messageEn, /do not start from the same biological baseline/);
assert.ok(bilingual(mismatch.findings));

// negative control: comparable baselines must NOT raise the finding
const okBaseline = assessBaselineMismatch(plan, {
  'MicrosomalBlank__vs__Vehicle': 5,
  'TAC_0.5__vs__Vehicle': 300,
  'MET_0.5__vs__MicrosomalBlank': 400
});
assert.equal(okBaseline.ok, true);
assert.equal(okBaseline.findings.length, 0);

// missing counts must not fabricate a verdict
assert.equal(assessBaselineMismatch(plan, {}).code, 'not_assessable');

/* ---- 3. replication ---------------------------------------------------- */

const replicated = [];
for (const s of ['E1', 'E2', 'E3']) {
  for (const c of ['ctrl', 'treat']) {
    for (const t of ['D14', 'D28']) replicated.push({ subject_id: s, condition: c, timepoint: t });
  }
}
const rep = assessReplication(replicated);
assert.equal(rep.ok, true);
assert.equal(rep.code, 'replicated');
assert.equal(rep.minUnits, 3);
assert.equal(rep.findings.length, 0);

// one library per cell -> inference must be refused
const unreplicated = [
  { subject_id: 'L1', condition: 'ctrl', timepoint: 'D28' },
  { subject_id: 'L2', condition: 'treat', timepoint: 'D28' }
];
const norep = assessReplication(unreplicated);
assert.equal(norep.ok, false);
assert.equal(norep.code, 'no_replication');
assert.match(norep.findings[0].messageEn, /assumed rather than an observed dispersion/);

// technical replicates of ONE subject must not count as replication
const technical = [
  { subject_id: 'S1', condition: 'ctrl', timepoint: 'D28' },
  { subject_id: 'S1', condition: 'ctrl', timepoint: 'D28' },
  { subject_id: 'S1', condition: 'ctrl', timepoint: 'D28' },
  { subject_id: 'S2', condition: 'treat', timepoint: 'D28' }
];
assert.equal(assessReplication(technical).code, 'no_replication');
assert.equal(assessReplication([{ subject_id: 'a', condition: 'x', timepoint: 'T' },
  { subject_id: 'b', condition: 'x', timepoint: 'T' }]).code, 'minimal_replication');
assert.equal(assessReplication([]).code, 'no_metadata');

/* ---- 4. variance partition --------------------------------------------- */

// batch effect deliberately as large as the treatment effect
const nS = 12;
const batch = Array.from({ length: nS }, (_, i) => 'B' + (i % 3));
const cond = Array.from({ length: nS }, (_, i) => (i < 6 ? 'ctrl' : 'treat'));
const matrix = [];
for (let f = 0; f < 60; f += 1) {
  matrix.push(Array.from({ length: nS }, (_, i) =>
    (cond[i] === 'treat' ? 2 : 0) + (i % 3) * 2 + ((f * 7 + i * 13) % 5) * 0.01));
}
const vp = variancePartition(matrix, { condition: cond, batch });
assert.equal(vp.nFeatures, 60);
const byName = Object.fromEntries(vp.factors.map((f) => [f.factor, f.varianceFractionMarginal]));
assert.ok(byName.batch > 0.2, 'batch should carry real variance');
assert.ok(byName.condition > 0.2, 'condition should carry real variance');
assert.ok(has(vp.findings, 'nuisance_factor_dominant'));
assert.ok(bilingual(vp.findings));

// negative control: a nuisance factor carrying nothing must not be flagged
const flat = Array.from({ length: 40 }, (_, f) =>
  Array.from({ length: nS }, (_, i) => (cond[i] === 'treat' ? 5 : 0) + ((f + i) % 3) * 0.001));
const vp2 = variancePartition(flat, { condition: cond, batch });
assert.ok(!has(vp2.findings, 'nuisance_factor_dominant'));

/* ---- 5. p-value calibration -------------------------------------------- */

const uniform = Array.from({ length: 5000 }, (_, i) => (i + 0.5) / 5000);
assert.equal(assessPValueCalibration(uniform).code, 'calibrated');

// conservative: the fixed-dispersion signature seen in the real transcriptomics
const conservative = Array.from({ length: 5000 }, (_, i) => 0.3 + 0.7 * ((i + 0.5) / 5000));
const cons = assessPValueCalibration(conservative);
assert.equal(cons.code, 'conservative');
assert.ok(cons.fractionBelow05 < 0.025);
assert.match(cons.findings[0].messageEn, /not calibrated/);

// no flat tail at all: every test "significant-ish", the unmodelled-nuisance signature
const inflated = Array.from({ length: 5000 }, (_, i) => ((i + 0.5) / 5000) * 0.15);
assert.equal(assessPValueCalibration(inflated).code, 'inflated');
assert.equal(assessPValueCalibration([0.1, 0.2]).code, 'too_few_tests');

/* ---- 6. feature-table hygiene ------------------------------------------ */

const ids = ['P55036', 'Reverse_sp|Q8WZ42|TITIN_HUMAN', 'CON__P02769', 'Glutathione', 'Lysine', 'Alanine'];
const tbl = [
  [10, 11, 12], [9, 8, 10], [5, 5, 5], [0, 0, 0], [0, 0, 0], [3, 4, 5]
];
const fa = auditFeatureTable(ids, tbl);
assert.deepEqual(fa.decoys, ['Reverse_sp|Q8WZ42|TITIN_HUMAN']);
assert.deepEqual(fa.contaminants, ['CON__P02769']);
assert.deepEqual(fa.ubiquitousAbsent.sort(), ['Glutathione', 'Lysine']);
assert.ok(has(fa.findings, 'decoy_features_present'));
assert.ok(has(fa.findings, 'zeros_are_missing'));
assert.ok(bilingual(fa.findings));

// negative control: a clean table raises nothing
const clean = auditFeatureTable(['A', 'B'], [[1, 2, 3], [4, 5, 6]]);
assert.equal(clean.findings.length, 0);
assert.equal(clean.decoys.length, 0);

/* ---- 7. non-biological signal detector --------------------------------- */

// Reproduces the NADPH confound: the pool contains a cofactor the blank lacks.
const nbIds = ['NAD', 'Niacinamide', 'Hypoxanthine', 'Inosine', 'Adenosine', 'RealProtein'];
const groupOf = ['blank', 'blank', 'blank', 'pool', 'pool', 'pool'];
const nbMatrix = [
  [0, 0, 0, 120, 140, 130],        // NAD: present in the pool arm only
  [0.01, 0.01, 0.01, 9.4, 10.1, 9.8], // nicotinamide: ~1000x
  [0.002, 0.003, 0.002, 1.1, 1.3, 1.2],
  [0.01, 0.012, 0.011, 0.2, 0.21, 0.19],
  [0.5, 0.52, 0.49, 2.9, 3.1, 3.0],
  [10, 11, 10.5, 12, 11.5, 12.5]   // a genuine, modest biological change
];
const nb = detectNonBiologicalSignals({
  featureIds: nbIds, matrix: nbMatrix, groupOf, testArm: 'pool', referenceArm: 'blank'
});
const flagged = nb.candidates.map((c) => c.feature);
assert.ok(flagged.includes('NAD'));
assert.ok(flagged.includes('Niacinamide'));
assert.ok(!flagged.includes('RealProtein'), 'a modest genuine change must not be flagged');
assert.ok(nb.families.some((f) => /NADPH/.test(f.family)), 'the cofactor family must be recognised');
assert.ok(has(nb.findings, 'formulation_family_detected'));
assert.ok(bilingual(nb.findings));

// negative control: ordinary biology raises nothing
const calm = detectNonBiologicalSignals({
  featureIds: ['X', 'Y'], matrix: [[10, 11, 12, 14, 15, 13], [5, 5, 6, 6, 7, 6]],
  groupOf, testArm: 'pool', referenceArm: 'blank'
});
assert.equal(calm.candidates.length, 0);
assert.equal(calm.families.length, 0);
assert.ok(has(detectNonBiologicalSignals({
  featureIds: ['X'], matrix: [[1, 2, 3, 4, 5, 6]], groupOf, testArm: 'nope', referenceArm: 'blank'
}).findings, 'arms_not_found'));

/* ---- 8. reproducibility across blocks ---------------------------------- */

const n = 400;
const seedy = (k) => ((Math.sin(k * 12.9898) * 43758.5453) % 1 + 1) % 1 - 0.5;
const signal = Array.from({ length: n }, (_, i) => (i < 40 ? 2.5 : 0));
const repro = {
  featureIds: Array.from({ length: n }, (_, i) => 'F' + i),
  effectByBlock: {
    D14: signal.map((v, i) => v + seedy(i + 1) * 0.4),
    D21: signal.map((v, i) => v + seedy(i + 1000) * 0.4),
    D28: signal.map((v, i) => v + seedy(i + 2000) * 0.4)
  }
};
const rr = reproducibilityAcrossBlocks(repro);
assert.equal(rr.correlations.length, 3);
assert.ok(rr.medianPearson > 0.5, 'a real effect must correlate across blocks');
assert.equal(rr.consistent.length, 40);
assert.ok(has(rr.findings, 'effects_reproducible'));

// pure noise: the key negative control -- an unreplicated layer whose
// apparent "signature" does not survive comparison across experiments
const noiseOnly = {
  featureIds: repro.featureIds,
  effectByBlock: {
    D14: Array.from({ length: n }, (_, i) => seedy(i + 7) * 3),
    D21: Array.from({ length: n }, (_, i) => seedy(i + 777) * 3),
    D28: Array.from({ length: n }, (_, i) => seedy(i + 77777) * 3)
  }
};
const nr = reproducibilityAcrossBlocks(noiseOnly);
assert.ok(Math.abs(nr.medianPearson) < 0.2);
assert.ok(has(nr.findings, 'effects_not_reproducible'));
assert.match(nr.findings[0].messageEn, /must not be presented as a signature/);
assert.ok(has(reproducibilityAcrossBlocks({ featureIds: ['A'], effectByBlock: { only: [1] } })
  .findings, 'too_few_blocks'));

/* ---- 9. aggregate ------------------------------------------------------ */

const audit = auditDesign({
  armPlan: PLAN_SPEC,
  significantByContrast: {
    'MicrosomalBlank__vs__Vehicle': 1381,
    'TAC_0.5__vs__Vehicle': 1,
    'MET_0.5__vs__MicrosomalBlank': 943
  },
  metadata: unreplicated,
  featureIds: ids,
  matrix: tbl,
  pValues: conservative,
  nonBiological: {
    featureIds: nbIds, matrix: nbMatrix, groupOf, testArm: 'pool', referenceArm: 'blank'
  },
  reproducibility: noiseOnly
});
assert.equal(audit.ok, false);
for (const code of ['baseline_mismatch', 'no_replication', 'decoy_features_present',
  'zeros_are_missing', 'formulation_family_detected', 'effects_not_reproducible']) {
  assert.ok(audit.blocking.includes(code), 'missing blocking finding: ' + code);
}
assert.ok(has(audit.findings, 'pvalues_conservative'));
assert.ok(bilingual(audit.findings));
// errors are ordered before warnings and info
const ranks = audit.findings.map((f) => ({ error: 3, warning: 2, info: 1 }[f.severity]));
assert.deepEqual(ranks, ranks.slice().sort((a, b) => b - a));

// a clean study must pass with no blocking finding
const cleanAudit = auditDesign({
  armPlan: { arms: ['ctrl', 'treat'], controlOf: { treat: 'ctrl' } },
  significantByContrast: { 'treat__vs__ctrl': 120 },
  metadata: replicated,
  featureIds: ['A', 'B'],
  matrix: [[1, 2, 3], [4, 5, 6]],
  pValues: uniform
});
assert.equal(cleanAudit.ok, true);
assert.equal(cleanAudit.blocking.length, 0);

// an empty call must not throw and must not invent findings
const empty = auditDesign({});
assert.equal(empty.ok, true);
assert.equal(empty.findings.length, 0);

console.log('PASS: design audit — contrast plan, baseline mismatch, replication, variance partition, p-value calibration, feature-table hygiene, non-biological signals and between-block reproducibility, each with negative controls');
