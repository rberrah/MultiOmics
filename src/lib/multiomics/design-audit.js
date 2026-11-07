// @ts-nocheck
// Design audit: checks that run BEFORE a model is chosen, and that describe
// why a result may not mean what it appears to mean.
//
// Every check here exists because it caught a real error in a real study:
// an exposure experiment in which (a) each treatment arm had its own
// preparation control, (b) one of those controls was itself strongly
// bioactive, (c) the two preparations differed by a 2 mM cofactor whose
// degradation products dominated the metabolome, (d) the supplied protein
// table still contained reverse-database decoys, (e) zeros in a targeted
// panel meant "not quantified", and (f) one layer had no replication at all
// yet still carried p-values.
//
// Pure functions, no dependencies. Messages are bilingual (fr, en) to match
// the rest of src/lib/multiomics.

const FINDING = (code, severity, message, messageEn, detail = {}) =>
  ({ code, severity, message, messageEn, detail });

const num = (v) => (Number.isFinite(Number(v)) ? Number(v) : null);
const mean = (a) => (a.length ? a.reduce((s, x) => s + x, 0) / a.length : NaN);
const uniq = (a) => Array.from(new Set(a));

/* ------------------------------------------------------------------ *
 * 1. Contrast plan: every arm against ITS OWN control
 * ------------------------------------------------------------------ */

/**
 * Build the full contrast set implied by an arm -> control mapping.
 *
 * Treatment arms are compared with their own control, the controls are
 * compared with each other and with an untreated arm when one exists, and
 * arms that do not share a control get a matrix-corrected head-to-head
 * ((A - controlA) - (B - controlB)) instead of a raw difference.
 *
 * @param {{arms: string[], controlOf: Record<string,string>, untreated?: string|null}} spec
 */
export function buildContrastPlan(spec) {
  const arms = Array.isArray(spec?.arms) ? spec.arms.filter(Boolean) : [];
  const controlOf = spec?.controlOf && typeof spec.controlOf === 'object' ? spec.controlOf : {};
  const untreated = spec?.untreated ?? null;
  if (arms.length < 2) throw new Error('buildContrastPlan needs at least 2 arms');

  const controls = uniq(arms.map((a) => controlOf[a]).filter(Boolean));
  const unmapped = arms.filter((a) => !controlOf[a] && !controls.includes(a) && a !== untreated);
  const treatments = arms.filter((a) => Boolean(controlOf[a]));

  const contrasts = [];
  for (const arm of treatments) {
    contrasts.push({
      id: `${arm}__vs__${controlOf[arm]}`, kind: 'treatment_vs_own_control',
      test: arm, reference: controlOf[arm], expression: `${arm} - ${controlOf[arm]}`
    });
  }
  for (let i = 0; i < controls.length; i += 1) {
    for (let j = i + 1; j < controls.length; j += 1) {
      contrasts.push({
        id: `${controls[i]}__vs__${controls[j]}`, kind: 'control_vs_control',
        test: controls[i], reference: controls[j],
        expression: `${controls[i]} - ${controls[j]}`
      });
    }
  }
  if (untreated) {
    for (const c of controls) {
      contrasts.push({
        id: `${c}__vs__${untreated}`, kind: 'control_vs_untreated',
        test: c, reference: untreated, expression: `${c} - ${untreated}`
      });
    }
  }
  for (let i = 0; i < treatments.length; i += 1) {
    for (let j = i + 1; j < treatments.length; j += 1) {
      const a = treatments[i]; const b = treatments[j];
      const sameControl = controlOf[a] === controlOf[b];
      contrasts.push(sameControl
        ? { id: `${a}__vs__${b}`, kind: 'head_to_head_shared_control',
            test: a, reference: b, expression: `${a} - ${b}` }
        : { id: `${a}__vs__${b}__corrected`, kind: 'head_to_head_corrected',
            test: a, reference: b,
            expression: `(${a} - ${controlOf[a]}) - (${b} - ${controlOf[b]})`,
            note: 'arms do not share a control; raw difference is not interpretable' });
    }
  }

  const findings = [];
  if (unmapped.length) {
    findings.push(FINDING('arm_without_control', 'error',
      'Ces bras n’ont pas de contrôle déclaré : ' + unmapped.join(', ') +
      '. Déclarez le contrôle apparié de chaque bras avant toute comparaison.',
      'These arms have no declared control: ' + unmapped.join(', ') +
      '. Declare the matched control of every arm before comparing anything.',
      { arms: unmapped }));
  }
  if (controls.length > 1) {
    findings.push(FINDING('multiple_controls', 'info',
      'Le plan contient ' + controls.length + ' contrôles distincts (' + controls.join(', ') +
      '). Chaque bras est comparé au sien, et les contrôles sont comparés entre eux.',
      'The design contains ' + controls.length + ' distinct controls (' + controls.join(', ') +
      '). Each arm is compared with its own, and the controls are compared with each other.',
      { controls }));
  }
  return { contrasts, controls, treatments, untreated, findings };
}

/**
 * Raise a baseline-mismatch finding when a control-vs-control contrast is
 * itself larger than the treatment effects it is supposed to anchor.
 *
 * @param {ReturnType<typeof buildContrastPlan>} plan
 * @param {Record<string, number>} significantByContrast  contrast id -> count
 */
export function assessBaselineMismatch(plan, significantByContrast) {
  const counts = significantByContrast && typeof significantByContrast === 'object'
    ? significantByContrast : {};
  // Counts may be keyed in either direction (A__vs__B or B__vs__A); a contrast
  // and its reverse have the same number of significant features.
  const get = (id) => {
    const direct = num(counts[id]);
    if (direct !== null) return direct;
    const parts = String(id).split('__vs__');
    if (parts.length !== 2) return null;
    return num(counts[`${parts[1]}__vs__${parts[0]}`]);
  };
  const controlContrasts = plan.contrasts
    .filter((c) => c.kind === 'control_vs_control' || c.kind === 'control_vs_untreated')
    .map((c) => ({ ...c, n: get(c.id) }))
    .filter((c) => c.n !== null);
  const treatmentContrasts = plan.contrasts
    .filter((c) => c.kind === 'treatment_vs_own_control')
    .map((c) => ({ ...c, n: get(c.id) }))
    .filter((c) => c.n !== null);
  if (!controlContrasts.length || !treatmentContrasts.length) {
    return { ok: true, code: 'not_assessable', findings: [] };
  }
  const worst = controlContrasts.reduce((m, c) => (c.n > m.n ? c : m), controlContrasts[0]);
  const biggestTreatment = treatmentContrasts.reduce((m, c) => (c.n > m.n ? c : m), treatmentContrasts[0]);
  if (worst.n > biggestTreatment.n) {
    return {
      ok: false, code: 'baseline_mismatch',
      findings: [FINDING('baseline_mismatch', 'error',
        'Le contraste entre contrôles (' + worst.id + ', ' + worst.n + ' variables) dépasse le plus grand effet de traitement (' +
        biggestTreatment.id + ', ' + biggestTreatment.n + '). Les bras ne partent donc pas du même état biologique et une comparaison directe entre eux n’est pas interprétable.',
        'The control-versus-control contrast (' + worst.id + ', ' + worst.n + ' features) exceeds the largest treatment effect (' +
        biggestTreatment.id + ', ' + biggestTreatment.n + '). The arms therefore do not start from the same biological baseline, and a direct comparison between them is not interpretable.',
        { controlContrast: worst.id, controlCount: worst.n,
          treatmentContrast: biggestTreatment.id, treatmentCount: biggestTreatment.n })]
    };
  }
  return { ok: true, code: 'baselines_comparable', findings: [] };
}

/* ------------------------------------------------------------------ *
 * 2. Replication
 * ------------------------------------------------------------------ */

/**
 * Count independent units per analysis cell. Technical replicates and repeated
 * assays of the same subject do not count.
 *
 * @param {Array<{subject_id:string,condition:string,timepoint?:string}>} rows
 */
export function assessReplication(rows) {
  const list = Array.isArray(rows) ? rows : [];
  const cells = new Map();
  for (const r of list) {
    const key = String(r?.condition ?? '') + '@@' + String(r?.timepoint ?? '');
    if (!cells.has(key)) cells.set(key, new Set());
    cells.get(key).add(String(r?.subject_id ?? ''));
  }
  const perCell = Array.from(cells, ([key, set]) => {
    const [condition, timepoint] = key.split('@@');
    return { condition, timepoint, independentUnits: set.size };
  });
  if (!perCell.length) {
    return { ok: false, code: 'no_metadata', perCell, findings: [FINDING('no_metadata', 'error',
      'Aucune métadonnée exploitable.', 'No usable metadata.')] };
  }
  const minUnits = Math.min(...perCell.map((c) => c.independentUnits));
  if (minUnits < 2) {
    const bad = perCell.filter((c) => c.independentUnits < 2);
    return { ok: false, code: 'no_replication', perCell, minUnits,
      findings: [FINDING('no_replication', 'error',
        'Certaines cellules du plan ne contiennent qu’une seule unité indépendante (' +
        bad.map((c) => c.condition + (c.timepoint ? '/' + c.timepoint : '')).join(', ') +
        '). Aucune variance intra-groupe ne peut être estimée : les p-values et FDR affichés dépendraient d’une dispersion supposée, pas observée. Utilisez l’analyse de reproductibilité entre blocs.',
        'Some design cells contain a single independent unit (' +
        bad.map((c) => c.condition + (c.timepoint ? '/' + c.timepoint : '')).join(', ') +
        '). No within-group variance can be estimated, so any p-value or FDR would reflect an assumed rather than an observed dispersion. Use the between-block reproducibility analysis instead.',
        { cells: bad })] };
  }
  if (minUnits < 3) {
    return { ok: true, code: 'minimal_replication', perCell, minUnits,
      findings: [FINDING('minimal_replication', 'warning',
        'Deux unités indépendantes seulement dans au moins une cellule : puissance très faible et variance mal estimée.',
        'Only two independent units in at least one cell: very low power and a poorly estimated variance.',
        { minUnits })] };
  }
  return { ok: true, code: 'replicated', perCell, minUnits, findings: [] };
}

/* ------------------------------------------------------------------ *
 * 3. Variance partition
 * ------------------------------------------------------------------ */

function betweenGroupSS(values, groups) {
  const grand = mean(values);
  const byGroup = new Map();
  for (let i = 0; i < values.length; i += 1) {
    const g = groups[i];
    if (!byGroup.has(g)) byGroup.set(g, []);
    byGroup.get(g).push(values[i]);
  }
  let ss = 0;
  for (const arr of byGroup.values()) {
    const m = mean(arr);
    ss += arr.length * (m - grand) * (m - grand);
  }
  return ss;
}

/**
 * Mean share of per-feature variance explained by each declared factor.
 *
 * This is the MARGINAL (one-way) decomposition: each factor is evaluated on
 * its own. When factors are confounded the fractions can sum above 1, and
 * that is itself reported, because a confounded nuisance factor is exactly
 * what makes an unblocked test misleading.
 *
 * @param {number[][]} matrix features x samples
 * @param {Record<string,string[]>} factors factor name -> value per sample
 */
export function variancePartition(matrix, factors) {
  const rows = Array.isArray(matrix) ? matrix : [];
  const names = Object.keys(factors || {});
  if (!rows.length || !names.length) return { factors: [], confounded: false, nFeatures: 0 };
  const acc = Object.fromEntries(names.map((n) => [n, []]));
  let used = 0;
  for (const row of rows) {
    const vals = []; const idx = [];
    row.forEach((v, i) => { const n = num(v); if (n !== null) { vals.push(n); idx.push(i); } });
    if (vals.length < 4) continue;
    const grand = mean(vals);
    const total = vals.reduce((s, v) => s + (v - grand) * (v - grand), 0);
    if (!(total > 0)) continue;
    used += 1;
    for (const name of names) {
      const g = idx.map((i) => String(factors[name][i]));
      acc[name].push(betweenGroupSS(vals, g) / total);
    }
  }
  const out = names
    .map((name) => ({ factor: name, varianceFractionMarginal: used ? mean(acc[name]) : NaN }))
    .sort((a, b) => b.varianceFractionMarginal - a.varianceFractionMarginal);
  const sum = out.reduce((s, f) => s + (Number.isFinite(f.varianceFractionMarginal) ? f.varianceFractionMarginal : 0), 0);
  const findings = [];
  const biggest = out[0];
  const nuisance = out.filter((f) => /batch|replicate|serie|series|run|plate|experiment/i.test(f.factor));
  for (const f of nuisance) {
    if (f.varianceFractionMarginal >= 0.5 * (biggest?.varianceFractionMarginal ?? 0) &&
        f.varianceFractionMarginal > 0.05) {
      findings.push(FINDING('nuisance_factor_dominant', 'warning',
        'Le facteur « ' + f.factor + ' » explique ' + (100 * f.varianceFractionMarginal).toFixed(1) +
        ' % de la variance par variable, soit autant que le traitement. Il doit figurer dans le modèle ; un test non ajusté sur ce facteur serait biaisé et sous-puissant.',
        'Factor "' + f.factor + '" explains ' + (100 * f.varianceFractionMarginal).toFixed(1) +
        '% of per-feature variance, comparable to the treatment itself. It must be in the model; a test that does not block on it is both biased and underpowered.',
        { factor: f.factor, fraction: f.varianceFractionMarginal }));
    }
  }
  if (sum > 1.1) {
    findings.push(FINDING('factors_confounded', 'warning',
      'Les parts de variance marginales totalisent ' + sum.toFixed(2) +
      ' (> 1) : les facteurs déclarés sont partiellement confondus et ne peuvent pas être séparés par ce plan.',
      'Marginal variance fractions sum to ' + sum.toFixed(2) +
      ' (> 1): the declared factors are partly confounded and cannot be separated by this design.',
      { sum }));
  }
  return { factors: out, confounded: sum > 1.1, nFeatures: used, findings };
}

/* ------------------------------------------------------------------ *
 * 4. p-value calibration
 * ------------------------------------------------------------------ */

/**
 * A calibrated null gives a uniform p-value distribution: about 5% below 0.05.
 * Far fewer means the test is conservative (typically a fixed dispersion) and
 * its FDR is not interpretable; far more with a flat tail means an unmodelled
 * nuisance factor.
 */
export function assessPValueCalibration(pValues) {
  const p = (Array.isArray(pValues) ? pValues : []).map(num)
    .filter((v) => v !== null && v >= 0 && v <= 1);
  if (p.length < 200) {
    return { ok: true, code: 'too_few_tests', n: p.length, findings: [] };
  }
  const frac = p.filter((v) => v < 0.05).length / p.length;
  const upper = p.filter((v) => v > 0.5).length / p.length;
  if (frac < 0.025) {
    return { ok: false, code: 'conservative', n: p.length, fractionBelow05: frac, findings: [
      FINDING('pvalues_conservative', 'warning',
        'Seulement ' + (100 * frac).toFixed(2) + ' % des tests ont p < 0,05, contre 5 % attendus sous l’hypothèse nulle. Le test est conservateur — typiquement une dispersion fixée à la main — et le FDR affiché n’est pas calibré.',
        'Only ' + (100 * frac).toFixed(2) + '% of tests have p < 0.05 versus the 5% expected under the null. The test is conservative — typically a hand-set dispersion — and the displayed FDR is not calibrated.',
        { fractionBelow05: frac })] };
  }
  if (frac > 0.25 && upper < 0.35) {
    return { ok: false, code: 'inflated', n: p.length, fractionBelow05: frac, findings: [
      FINDING('pvalues_inflated', 'warning',
        (100 * frac).toFixed(1) + ' % des tests ont p < 0,05 et la distribution n’a pas de queue plate. Cela traduit souvent un facteur de nuisance non modélisé ou une non-indépendance des échantillons.',
        (100 * frac).toFixed(1) + '% of tests have p < 0.05 and the distribution has no flat tail. This usually indicates an unmodelled nuisance factor or non-independent samples.',
        { fractionBelow05: frac })] };
  }
  return { ok: true, code: 'calibrated', n: p.length, fractionBelow05: frac, findings: [] };
}

/* ------------------------------------------------------------------ *
 * 5. Feature-table hygiene
 * ------------------------------------------------------------------ */

const DECOY_RE = /(^|[|;_])(reverse|decoy|rev)[_|]/i;
const CONTAMINANT_RE = /(^|[|;])CON__|(^|[|;])contam/i;

// Analytes that are present in every living cell. If a targeted panel reports
// them as absent in every sample, the zeros encode "not quantified".
const UBIQUITOUS = ['atp', 'adp', 'amp', 'nad', 'nadh', 'nadp', 'glutathione',
  'glutamate', 'glutamic acid', 'glutamine', 'lysine', 'aspartate',
  'aspartic acid', 'citrate', 'citric acid', 'cysteine', 'taurine', 'alanine',
  'serine', 'glycine', 'lactate', 'lactic acid', 'pyruvate', 'pyruvic acid'];

const normName = (s) => String(s ?? '').toLowerCase().replace(/[_\-\s]+/g, ' ').trim();

/**
 * @param {string[]} featureIds
 * @param {Array<Array<number|null>>} matrix features x samples
 * @param {{zerosAreMissing?: boolean}} [options]
 */
export function auditFeatureTable(featureIds, matrix, options = {}) {
  const ids = Array.isArray(featureIds) ? featureIds.map(String) : [];
  const rows = Array.isArray(matrix) ? matrix : [];
  const decoys = ids.filter((id) => DECOY_RE.test(id));
  const contaminants = ids.filter((id) => CONTAMINANT_RE.test(id));

  const alwaysZero = []; const ubiquitousAbsent = []; let zeroCells = 0; let cells = 0;
  rows.forEach((row, i) => {
    const vals = (row || []).map(num);
    cells += vals.length;
    const zeros = vals.filter((v) => v === 0).length;
    zeroCells += zeros;
    const present = vals.filter((v) => v !== null && v !== 0).length;
    if (vals.length && present === 0) {
      alwaysZero.push(ids[i] ?? String(i));
      const n = normName(ids[i]);
      if (UBIQUITOUS.some((u) => n === u || n.includes(u))) ubiquitousAbsent.push(ids[i]);
    }
  });

  const findings = [];
  if (decoys.length) {
    findings.push(FINDING('decoy_features_present', 'error',
      decoys.length + ' entrées leurres (reverse/decoy) sont encore dans la table et seraient comptées comme des résultats. Elles doivent être retirées avant toute analyse.',
      decoys.length + ' reverse-database decoy entries are still in the table and would be counted as results. Remove them before any analysis.',
      { features: decoys.slice(0, 50), count: decoys.length }));
  }
  if (contaminants.length) {
    findings.push(FINDING('contaminant_features_present', 'warning',
      contaminants.length + ' entrées contaminants (CON__) sont présentes.',
      contaminants.length + ' contaminant entries (CON__) are present.',
      { features: contaminants.slice(0, 50), count: contaminants.length }));
  }
  if (ubiquitousAbsent.length) {
    findings.push(FINDING('zeros_are_missing', 'error',
      'Des analytes présents dans toute cellule vivante sont rapportés à zéro dans tous les échantillons (' +
      ubiquitousAbsent.join(', ') + '). Ces zéros signifient « non quantifié » et doivent être traités comme des valeurs manquantes, pas comme des zéros.',
      'Analytes present in every living cell are reported as zero in every sample (' +
      ubiquitousAbsent.join(', ') + '). These zeros mean "not quantified" and must be treated as missing, not as true zeros.',
      { features: ubiquitousAbsent }));
  } else if (alwaysZero.length && options.zerosAreMissing !== true) {
    findings.push(FINDING('all_zero_features', 'warning',
      alwaysZero.length + ' variables valent zéro dans tous les échantillons. Un zéro partout signifie en général « jamais mesuré » : confirmez l’encodage des valeurs manquantes.',
      alwaysZero.length + ' features are zero in every sample. A zero everywhere usually means "never measured": confirm how missing values are encoded.',
      { count: alwaysZero.length, features: alwaysZero.slice(0, 50) }));
  }
  return {
    decoys, contaminants, alwaysZero, ubiquitousAbsent,
    zeroFraction: cells ? zeroCells / cells : 0, findings
  };
}

/* ------------------------------------------------------------------ *
 * 6. Non-biological signal detector
 * ------------------------------------------------------------------ */

// Families of analytes that enter a culture from the formulation rather than
// from the cells. If several members of one family move together and hugely,
// the difference is chemistry, not biology.
const FORMULATION_FAMILIES = [
  { family: 'NADPH / NADH cofactor and its degradation products',
    members: ['nad', 'nadh', 'nadp', 'nadph', 'nicotinamide', 'niacinamide',
      'nicotinic acid', 'adenosine', 'inosine', 'hypoxanthine', 'xanthine',
      'uric acid', 'urate', 'adp ribose', 'ribose', 'd ribose',
      'gluconic acid', '6 phosphogluconate'] },
  { family: 'microsome / membrane preparation buffer',
    members: ['sucrose', 'glycerol', 'mannitol', 'tris', 'hexose', 'glucose',
      'pyridoxal', 'pyridoxine'] },
  { family: 'solvent and vehicle',
    members: ['dmso', 'ethanol', 'methanol', 'acetonitrile', 'formic acid'] },
  { family: 'medium supplement',
    members: ['hydrocortisone', 'insulin', 'transferrin', 'selenite',
      'epidermal growth factor', 'triiodothyronine', 'penicillin',
      'streptomycin', 'g 418', 'geneticin'] }
];

/**
 * Flag features whose behaviour between two arms is chemically, not
 * biologically, plausible.
 *
 * @param {{featureIds:string[], matrix:Array<Array<number|null>>,
 *          groupOf:string[], testArm:string, referenceArm:string,
 *          options?:{maxAbsLog2FC?:number, presenceThreshold?:number}}} spec
 */
export function detectNonBiologicalSignals(spec) {
  const ids = (spec?.featureIds || []).map(String);
  const rows = spec?.matrix || [];
  const groupOf = spec?.groupOf || [];
  const test = String(spec?.testArm ?? '');
  const ref = String(spec?.referenceArm ?? '');
  const maxAbs = Number.isFinite(Number(spec?.options?.maxAbsLog2FC))
    ? Number(spec.options.maxAbsLog2FC) : 5;
  const presenceCut = Number.isFinite(Number(spec?.options?.presenceThreshold))
    ? Number(spec.options.presenceThreshold) : 0.8;

  const iTest = groupOf.map((g, i) => (String(g) === test ? i : -1)).filter((i) => i >= 0);
  const iRef = groupOf.map((g, i) => (String(g) === ref ? i : -1)).filter((i) => i >= 0);
  if (!iTest.length || !iRef.length) {
    return { candidates: [], families: [], findings: [FINDING('arms_not_found', 'error',
      'Bras introuvables dans les métadonnées.', 'Arms not found in the metadata.')] };
  }

  const candidates = [];
  rows.forEach((row, k) => {
    const pick = (idx) => idx.map((i) => num(row?.[i])).filter((v) => v !== null);
    const a = pick(iTest); const b = pick(iRef);
    if (!a.length || !b.length) return;
    const detFrac = (arr, idx) => (idx.length ? arr.filter((v) => v > 0).length / idx.length : 0);
    const dTest = detFrac(a, iTest); const dRef = detFrac(b, iRef);
    const mTest = mean(a.filter((v) => v > 0)); const mRef = mean(b.filter((v) => v > 0));
    const reasons = [];
    let lfc = null;
    if (Number.isFinite(mTest) && Number.isFinite(mRef) && mTest > 0 && mRef > 0) {
      lfc = Math.log2(mTest / mRef);
      if (Math.abs(lfc) > maxAbs) reasons.push('implausible_fold_change');
    }
    if (dTest >= presenceCut && dRef <= 1 - presenceCut) reasons.push('present_in_one_arm_only');
    if (dRef >= presenceCut && dTest <= 1 - presenceCut) reasons.push('absent_in_one_arm_only');
    if (reasons.length) {
      candidates.push({ feature: ids[k] ?? String(k), log2FoldChange: lfc,
        detectionTest: dTest, detectionReference: dRef, reasons });
    }
  });

  const families = [];
  for (const fam of FORMULATION_FAMILIES) {
    const hit = candidates.filter((c) => {
      const n = normName(c.feature);
      return fam.members.some((m) => n === m || n.includes(m));
    });
    if (hit.length >= 2) families.push({ family: fam.family, features: hit.map((h) => h.feature) });
  }

  const findings = [];
  if (candidates.length) {
    findings.push(FINDING('non_biological_candidates', 'warning',
      candidates.length + ' variable(s) se comportent de façon chimiquement, non biologiquement, plausible entre « ' + test + ' » et « ' + ref +
      ' » (variation extrême ou présence dans un seul bras). Vérifiez s’il s’agit de composants du véhicule, du tampon, d’un cofacteur ou du milieu.',
      candidates.length + ' feature(s) behave in a chemically rather than biologically plausible way between "' + test + '" and "' + ref +
      '" (extreme change, or present in one arm only). Check whether they are components of the vehicle, buffer, cofactor or medium.',
      { count: candidates.length, features: candidates.slice(0, 50) }));
  }
  for (const f of families) {
    findings.push(FINDING('formulation_family_detected', 'error',
      'Plusieurs membres d’une même famille de formulation varient ensemble (' + f.family + ' : ' + f.features.join(', ') +
      '). La différence entre ces deux bras inclut une différence de composition, pas seulement un effet biologique.',
      'Several members of one formulation family move together (' + f.family + ': ' + f.features.join(', ') +
      '). The difference between these two arms includes a difference in composition, not only a biological effect.',
      f));
  }
  return { candidates, families, findings };
}

/* ------------------------------------------------------------------ *
 * 7. Reproducibility instead of inference
 * ------------------------------------------------------------------ */

function pearson(x, y) {
  const n = Math.min(x.length, y.length);
  const xs = []; const ys = [];
  for (let i = 0; i < n; i += 1) {
    const a = num(x[i]); const b = num(y[i]);
    if (a !== null && b !== null) { xs.push(a); ys.push(b); }
  }
  if (xs.length < 3) return NaN;
  const mx = mean(xs); const my = mean(ys);
  let sxy = 0; let sxx = 0; let syy = 0;
  for (let i = 0; i < xs.length; i += 1) {
    const dx = xs[i] - mx; const dy = ys[i] - my;
    sxy += dx * dy; sxx += dx * dx; syy += dy * dy;
  }
  return sxx > 0 && syy > 0 ? sxy / Math.sqrt(sxx * syy) : NaN;
}

/**
 * When a design has no within-cell replication but repeats across blocks
 * (time points, sites, runs), those blocks are de facto independent
 * experiments. A real effect reproduces across them; noise does not.
 *
 * @param {{featureIds:string[], effectByBlock:Record<string,number[]>,
 *          options?:{minAbsEffect?:number}}} spec
 */
export function reproducibilityAcrossBlocks(spec) {
  const ids = (spec?.featureIds || []).map(String);
  const byBlock = spec?.effectByBlock || {};
  const blocks = Object.keys(byBlock);
  const minAbs = Number.isFinite(Number(spec?.options?.minAbsEffect))
    ? Number(spec.options.minAbsEffect) : 1;
  if (blocks.length < 2) {
    return { blocks, correlations: [], consistent: [], findings: [FINDING('too_few_blocks', 'error',
      'Au moins deux blocs sont nécessaires.', 'At least two blocks are required.')] };
  }
  const correlations = [];
  for (let i = 0; i < blocks.length; i += 1) {
    for (let j = i + 1; j < blocks.length; j += 1) {
      correlations.push({ blockA: blocks[i], blockB: blocks[j],
        pearson: pearson(byBlock[blocks[i]], byBlock[blocks[j]]) });
    }
  }
  const consistent = [];
  for (let k = 0; k < ids.length; k += 1) {
    const vals = blocks.map((b) => num(byBlock[b][k]));
    if (vals.some((v) => v === null)) continue;
    const big = vals.every((v) => Math.abs(v) > minAbs);
    const sameSign = vals.every((v) => v > 0) || vals.every((v) => v < 0);
    if (big && sameSign) {
      consistent.push({ feature: ids[k], effects: vals, meanEffect: mean(vals) });
    }
  }
  consistent.sort((a, b) => Math.abs(b.meanEffect) - Math.abs(a.meanEffect));

  const rs = correlations.map((c) => c.pearson).filter(Number.isFinite);
  const medianR = rs.length ? rs.slice().sort((a, b) => a - b)[Math.floor(rs.length / 2)] : NaN;
  const findings = [];
  if (Number.isFinite(medianR) && medianR < 0.1) {
    findings.push(FINDING('effects_not_reproducible', 'error',
      'La corrélation médiane des tailles d’effet entre blocs indépendants est de ' + medianR.toFixed(3) +
      ' : l’effet observé ne se reproduit pas d’une expérience à l’autre et ne doit pas être présenté comme une signature.',
      'The median between-block correlation of effect sizes is ' + medianR.toFixed(3) +
      ': the observed effect does not reproduce across independent experiments and must not be presented as a signature.',
      { medianPearson: medianR, correlations }));
  } else if (Number.isFinite(medianR)) {
    findings.push(FINDING('effects_reproducible', 'info',
      'Corrélation médiane entre blocs : ' + medianR.toFixed(3) + ' ; ' + consistent.length +
      ' variable(s) dépassent le seuil d’effet dans le même sens dans tous les blocs.',
      'Median between-block correlation: ' + medianR.toFixed(3) + '; ' + consistent.length +
      ' feature(s) exceed the effect threshold in the same direction in every block.',
      { medianPearson: medianR, consistentCount: consistent.length }));
  }
  return { blocks, correlations, medianPearson: medianR, consistent, findings };
}

/* ------------------------------------------------------------------ *
 * 8. Aggregate
 * ------------------------------------------------------------------ */

const SEVERITY_RANK = { error: 3, warning: 2, info: 1 };

/**
 * Run every applicable check and return one ordered list of findings.
 * Checks are skipped silently when their inputs are absent.
 */
export function auditDesign(input = {}) {
  const findings = [];
  let plan = null;
  if (input.armPlan) {
    plan = buildContrastPlan(input.armPlan);
    findings.push(...plan.findings);
    if (input.significantByContrast) {
      findings.push(...assessBaselineMismatch(plan, input.significantByContrast).findings);
    }
  }
  if (Array.isArray(input.metadata)) {
    findings.push(...assessReplication(input.metadata).findings);
  }
  if (input.matrix && input.factors) {
    findings.push(...(variancePartition(input.matrix, input.factors).findings || []));
  }
  if (Array.isArray(input.pValues)) {
    findings.push(...assessPValueCalibration(input.pValues).findings);
  }
  if (input.featureIds && input.matrix) {
    findings.push(...auditFeatureTable(input.featureIds, input.matrix,
      { zerosAreMissing: input.zerosAreMissing }).findings);
  }
  if (input.nonBiological) {
    findings.push(...detectNonBiologicalSignals(input.nonBiological).findings);
  }
  if (input.reproducibility) {
    findings.push(...reproducibilityAcrossBlocks(input.reproducibility).findings);
  }
  findings.sort((a, b) => (SEVERITY_RANK[b.severity] ?? 0) - (SEVERITY_RANK[a.severity] ?? 0));
  const blocking = findings.filter((f) => f.severity === 'error');
  return {
    ok: blocking.length === 0,
    blocking: blocking.map((f) => f.code),
    findings,
    plan
  };
}
