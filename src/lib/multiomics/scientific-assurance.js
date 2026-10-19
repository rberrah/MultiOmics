/**
 * Conservative scientific evidence gate for the local multi-omics prototype.
 * This is NOT a validation certificate: a successful run, a small p-value,
 * or an installed R package cannot prove a study design is publishable.
 * Never grant a "publication-ready" status automatically.
 */
/**
 * "Confirmatory" is a scientific-workflow request, not a certificate.
 * Rules deliberately never auto-authorise a publication or clinical claim.
 * R methods must match EVERY actual input layer and the declared design.
 */
/** @param {any} result @param {{demo?: boolean}} [options] */
export function evaluateConfirmatoryReadiness(result, { demo = false } = {}) {
  const protocol = result?.protocol || {};
  const requested = protocol.analysisIntent === 'confirmatory';
  const layers = Object.entries(result?.layers || {});
  const methods = result?.referenceBackend?.status === 'ok'
    ? result.referenceBackend.methods || {} : {};
  /** @type {Array<{code:string,status:'pass'|'review'|'blocked',fr:string,en:string}>} */
  const checks = [];
  /** @param {string} code @param {'pass'|'review'|'blocked'} status @param {string} fr @param {string} en */
  const add = (code, status, fr, en) => checks.push({ code, status, fr, en });
  if (!requested) return {
    requested: false, status: 'not_requested', checks: [],
    certified: false, scope: 'No automatic scientific certification.'
  };

  if (demo) add('demo', 'blocked',
    'Démo synthétique : aucune conclusion biologique confirmatoire.',
    'Synthetic demonstration: no confirmatory biological inference.');
  if (!layers.length) add('empty', 'blocked',
    'Aucune omique analysée.', 'No analysed omics layer.');
  const conditions = result?.metadataSummary?.conditions || [];
  const groupObjective = protocol.objective === 'groups';
  const longitudinalObjective = protocol.objective === 'time';
  const survivalObjective = protocol.objective === 'outcome' &&
    protocol.outcomeType === 'survival';
  const independentDesign = groupObjective &&
    protocol.designType === 'independent' && !protocol.longitudinal &&
    conditions.length === 2;
  const pairedMethod = methods.transcriptomics_differential_deseq2;
  const pairedContract = pairedMethod?.pairedDesignContract;
  // A paired study is eligible for *review* only with the genuine server-side
  // subject-blocked contract, not just a designType flag or an HTTP 200.
  const pairedDesign = groupObjective &&
    protocol.designType === 'paired' && !protocol.longitudinal &&
    conditions.length === 2 && layers.length === 1 &&
    layers[0][0] === 'transcriptomics' &&
    layers[0][1]?.qc?.preprocessingAudit?.declaredValueType === 'raw_counts' &&
    pairedMethod?.status === 'ok' &&
    pairedContract?.status === 'eligible_for_reference_fitting_only' &&
    pairedContract?.pairing === 'same_subject_two_conditions' &&
    pairedContract?.certified === false &&
    pairedContract?.sample_identifiers_exact === true &&
    pairedContract?.all_subject_condition_cells_observed_once === true &&
    Number.isFinite(pairedContract?.biological_pairs) &&
    pairedContract.biological_pairs >= 4 &&
    pairedContract.independent_biological_subjects === pairedContract.biological_pairs &&
    pairedContract.model_samples === 2 * pairedContract.biological_pairs &&
    Number.isFinite(pairedContract.residual_degrees_of_freedom) &&
    pairedContract.residual_degrees_of_freedom >= 2 &&
    Number.isFinite(result?.metadataSummary?.subjects) &&
    result.metadataSummary.subjects === pairedContract.biological_pairs;
  const longitudinalDesign = longitudinalObjective &&
    (protocol.longitudinal === true || protocol.designType === 'repeated');
  const survivalDesign = survivalObjective &&
    protocol.designType === 'independent' && !protocol.longitudinal;
  const supportedDesign = independentDesign || pairedDesign ||
    longitudinalDesign || survivalDesign;

  if (!supportedDesign) add('design', 'blocked',
    'Le plan déclaré n’a pas de contrat R natif confirmatoire qualifié. Les voies prises en charge sont : deux groupes indépendants, RNA-seq apparié strict, groupe × temps longitudinal répété et survie indépendante.',
    'The declared design has no qualified native-R confirmatory contract. Supported paths are two independent groups, strictly paired RNA-seq, repeated group-by-time longitudinal analysis, and independent survival.');
  else if (pairedDesign) add('design', 'pass',
    'RNA-seq apparié avec blocage sujet R vérifié (éligibilité technique seulement).',
    'R subject-blocked paired RNA-seq verified (technical eligibility only).');
  else if (longitudinalDesign) add('design', 'pass',
    'Plan longitudinal répété reconnu ; l’éligibilité exacte des visites et du contraste groupe × temps est vérifiée par lmerTest avant l’inférence.',
    'Repeated longitudinal design recognised; exact visit and group-by-time contrast eligibility is checked by lmerTest before inference.');
  else if (survivalDesign) add('design', 'pass',
    'Plan de survie indépendant reconnu ; l’unicité des sujets, les événements et les diagnostics PH sont vérifiés dans le backend R.',
    'Independent survival design recognised; unique subjects, event counts and PH diagnostics are checked in the R backend.');
  else add('design', 'pass',
    'Plan à deux groupes indépendants reconnu.',
    'Two-independent-group design recognised.');

  if (!Object.keys(methods).length) add('backend', 'blocked',
    'Aucun moteur R de référence exécuté. Lancez le backend R pour toutes les omiques.',
    'No executed reference R backend. Run the R backend for every omics layer.');
  for (const [layer, data] of layers) {
    const qc = data?.qc || {};
    const declared = qc.preprocessingAudit?.declaredValueType;
    const key = longitudinalObjective ? layer + '_longitudinal'
      : survivalObjective ? layer + '_survival'
      : layer === 'transcriptomics' && declared === 'raw_counts'
        ? layer + '_differential_deseq2'
        : layer + '_differential';
    const method = methods[key];
    const status = method?.status;
    if (status !== 'ok') add('reference_'+layer, 'blocked',
      layer + ' : modèle R primaire compatible non exécuté (' + key + ').',
      layer + ': matching primary R model not executed (' + key + ').');
    else add('reference_'+layer,'pass',
      layer + ' : modèle R de référence exécuté (' + key + ').',
      layer + ': reference R model executed (' + key + ').');

    if ((longitudinalObjective || survivalObjective) &&
        layer === 'transcriptomics' && declared === 'raw_counts') {
      add('raw_count_reference_'+layer,'blocked',
        layer + ' : les comptages RNA-seq bruts ne sont pas acceptés dans cette voie confirmatoire ; la transformation CPM-log2 automatique reste exploratoire.',
        layer + ': raw RNA-seq counts are not accepted in this confirmatory path; automatic CPM-log2 transformation remains exploratory.');
    }
    if (longitudinalObjective && status === 'ok') {
      const fallbacks = Number(method?.summary?.randomSlopeFallbacks || 0);
      add('longitudinal_assumptions_'+layer,'review',
        layer + ' : vérifier linéarité de l’effet du temps, structure aléatoire, résidus, données manquantes et indépendance des sujets.',
        layer + ': review linear time effect, random-effects structure, residuals, missing visits and subject independence.');
      if (fallbacks > 0) add('longitudinal_random_slope_fallback_'+layer,'review',
        layer + ' : ' + fallbacks + ' variable(s) ont nécessité un fallback de pente aléatoire vers intercept seul ; cette adaptation doit être examinée avant interprétation confirmatoire.',
        layer + ': ' + fallbacks + ' feature(s) required random-slope fallback to random intercept only; this adaptation must be reviewed before confirmatory interpretation.');
    }
    if (survivalObjective && status === 'ok') {
      const phWarnings = Number(method?.summary?.statusCounts?.ph_assumption_warning || 0);
      add('survival_assumptions_'+layer,'review',
        layer + ' : vérifier risques proportionnels, mécanisme de censure, covariables, événements suffisants, competing risks et dépendance temporelle éventuelle.',
        layer + ': review proportional hazards, censoring mechanism, covariates, event information, competing risks and possible time-dependent effects.');
      if (phWarnings > 0) add('survival_ph_warning_'+layer,'review',
        layer + ' : ' + phWarnings + ' variable(s) ont un diagnostic cox.zph < 0,05 ; ne pas interpréter leur hazard ratio comme effet proportionnel constant sans modèle de sensibilité.',
        layer + ': ' + phWarnings + ' feature(s) have cox.zph < 0.05; do not interpret their hazard ratio as a constant proportional effect without sensitivity modelling.');
    }
    // Raw-count confirmatory intent is NOT eligible just because a server
    // advertises a successful DESeq2 fit. It must include the exact strict
    // sample/subject, contrast, additive design and raw-count preflight.
    // This remains eligibility only, never a certificate.
    if (groupObjective && layer === 'transcriptomics' && declared === 'raw_counts') {
      const input = methods[key]?.confirmatoryInputContract;
      const assay = methods[key]?.assayIntegrityContract;
      const requestedContrast = protocol.rnaPrimaryContrast;
      const fittedContrast = protocol.designType === 'paired'
        ? pairedContract?.explicit_primary_contrast : input?.condition_contrast;
      const exactDirectionalContrast =
        Array.isArray(requestedContrast) &&
        requestedContrast.length === 3 &&
        requestedContrast[0] === 'condition' &&
        requestedContrast[1] !== requestedContrast[2] &&
        Array.isArray(fittedContrast) &&
        fittedContrast.length === 3 &&
        fittedContrast.every((v, i) => v === requestedContrast[i]);
      const assayExact = assay?.status === 'assay_identifiers_exact' &&
        assay?.assay_loss_detected === false && assay?.certified === false &&
        Number.isFinite(assay?.measured_assays);
      const pairedProof = pairedDesign &&
        input?.status === 'eligible_for_reference_fitting_only' &&
        input?.pairing === 'same_subject_two_conditions' &&
        pairedContract?.contrast_was_prespecified === true &&
        assayExact &&
        assay.measured_assays >= pairedContract.model_samples;
      const independentProof = independentDesign &&
        input?.status === 'eligible_for_reference_fitting_only' &&
        input?.certified === false &&
        input?.rejected_silent_subject_loss === true &&
        Number.isFinite(input?.independent_subjects) &&
        input.independent_subjects >= 6 &&
        assayExact &&
        assay.measured_assays >= input.independent_subjects;
      if (status !== 'ok' || !exactDirectionalContrast ||
          !(pairedProof || independentProof)) {
        add('rnaseq_confirmatory_input_contract', 'blocked',
          'RNA-seq : le contrat strict R (comptages, unités biologiques, appariement éventuel, identifiants et contraste) est absent ou invalide. Un DESeq2 réussi ne suffit pas.',
          'RNA-seq: missing/invalid R evidence for raw counts, biological units, paired design where applicable, assay IDs, or exact contrast. DESeq2 success alone is insufficient.');
      } else {
        add('rnaseq_confirmatory_input_contract', 'pass',
          'RNA-seq : audit strict des entrées R réussi (éligibilité technique seulement).',
          'RNA-seq: strict R input audit passed (technical eligibility only).');
      }
    }
    for (const warning of qc.warnings || []) {
      add('qc_'+layer, 'review',
        layer + ' : contrôle qualité à examiner — ' + warning,
        layer + ': quality control requires review — ' + warning);
    }
    if (layer === 'metabolomics') add('metabolomics_validation','review',
      'Métabolomique : vérifier blancs, pooled-QC, identification et analyses de sensibilité MS.',
      'Metabolomics: review blanks, pooled QC, compound identification and MS sensitivity analyses.');
    // Group-dependent missingness can distort the estimand even if a reference
    // package reproduces the same observed-data p-value. This is a risk flag,
    // not an algorithmic way to distinguish MAR, MCAR and MNAR mechanisms.
    // Differential missingness by group is NOT a sufficient MNAR screen:
    // both groups may have the same missingness rate but censor opposite
    // tails, causing severe null false-positive inflation. A conservative
    // abundance-level threshold is an operational review trigger, NOT proof
    // of MCAR when data fall below it.
    const missingFraction = Number(qc.medianMissingFraction || 0);
    if (missingFraction >= 0.20) add('high_missingness_'+layer,'blocked',
      layer + ' : le taux médian de valeurs manquantes par profil mesuré atteint au moins 20 %. Sans modèle d’observation et analyses de sensibilité préspécifiés, le parcours confirmatoire reste bloqué (seuil de précaution, non diagnostic MNAR).',
      layer + ': median missingness across measured sample profiles is at least 20%. Confirmatory preparation requires a prespecified observation model and sensitivity analyses (caution threshold, not a diagnosis of MNAR).');
    else if (missingFraction > 0) add('missingness_mechanism_'+layer,'review',
      layer + ' : certaines valeurs sont manquantes. Des taux identiques entre groupes n’excluent pas une dépendance aux valeurs non observées ; vérifier la sensibilité à MCAR, MAR et MNAR.',
      layer + ': some values are missing. Equal missingness rates between groups do not exclude dependence on unobserved values; review sensitivity to MCAR, MAR and MNAR.');
    const differentialMissingness = Number(qc.differentialMissingness?.flaggedFeatures || 0);
    if (differentialMissingness > 0) add('missingness_imbalance_'+layer,'blocked',
      layer + ' : ' + differentialMissingness + ' variable(s) présentent un déséquilibre de valeurs manquantes entre les groupes (≥30 points). Examiner le mécanisme, les exclusions et les analyses de sensibilité avant une conclusion confirmatoire.',
      layer + ': ' + differentialMissingness + ' feature(s) have differential missingness between groups (≥30 percentage points). Review missingness mechanisms, exclusions and sensitivity analyses before confirmatory claims.');
    const mnarImputationCount = Number(qc.msQc?.mnar?.imputedValues || 0);
    if (mnarImputationCount > 0) add('mnar_imputation_'+layer,'blocked',
      layer + ' : imputation de valeurs manquantes supposées censurées appliquée (' + mnarImputationCount + ' valeurs). Une analyse de sensibilité à des hypothèses MNAR alternatives est requise.',
      layer + ': left-censored missing-value imputation was applied (' + mnarImputationCount + ' values). Sensitivity to alternative MNAR assumptions is required.');

    if (groupObjective && layer === 'transcriptomics' && declared === 'raw_counts') add('rnaseq_design','review',
      'RNA-seq : vérifier indépendamment la dispersion, les facteurs d’ajustement et les contrastes DESeq2.',
      'RNA-seq: independently check dispersion, adjustment variables and DESeq2 contrasts.');
    if (groupObjective && layer === 'transcriptomics' && declared === 'raw_counts') add('rnaseq_null_fdr_calibration','review',
      'RNA-seq : des simulations nulles indépendantes ont montré un excès de faux appels dans certains scénarios DESeq2. Les q-values demandent une vérification de calibration adaptée au plan ; revue méthodologique indispensable.',
      'RNA-seq: independent null simulations showed excess false calls for DESeq2 in some scenarios. q-values need design-specific calibration assessment and independent methodological review.');
  }
  if ((result?.preAnalysisDiagnostics?.blockers || []).length) add('study_design', 'blocked',
    'Le diagnostic préalable contient des blocages.', 'Pre-analysis diagnostics report blockers.');
  if (result?.predictiveOutcome?.status === 'ok') add('prediction', 'review',
    'La prédiction interne ne constitue pas une validation externe du biomarqueur.',
    'Internal prediction is not external biomarker validation.');
  if (layers.length > 1) add('multiplicity','review',
    'Préspécifier les familles de tests et le critère principal ; des q-values par omique ne contrôlent pas automatiquement l’erreur globale.',
    'Prespecify testing families and the primary endpoint; per-layer BH q-values do not automatically control study-wide error.');
  add('preregistration','review',
    'Vérifier le protocole préspécifié, l’indépendance des réplicats, les exclusions et les facteurs de confusion.',
    'Review the prespecified protocol, replicate independence, exclusions and confounders.');
  add('external_review','review',
    'Une vérification humaine indépendante reste nécessaire avant une conclusion scientifique confirmatoire.',
    'Independent human review is still required before confirmatory scientific claims.');

  return {
    requested: true,
    status: checks.some(x => x.status === 'blocked') ? 'blocked' : 'independent_review_required',
    checks,
    certified: false,
    scope: 'Software eligibility checks; not confirmation of hypotheses, clinical validation or publication approval.'
  };
}

/** @param {any} profile */
function profileHasMissingMeasurement(profile) {
  return Number(profile?.missingFraction) > 0;
}

/** @param {any} result @param {{demo?: boolean}} [options] */
export function assessScientificAssurance(result, { demo = false } = {}) {
  const layers = Object.entries(result?.layers || {});
  const protocol = result?.protocol || {};
  const exploratory = protocol.objective === 'explore';
  /** @type {Array<{code:string,fr:string,en:string,level:string}>} */
  const notes = [];
  /** @param {string} code @param {string} fr @param {string} en @param {string} [level] */
  const add = (code, fr, en, level = 'caution') => notes.push({ code, fr, en, level });
  const n = Number(result?.metadataSummary?.subjects || 0);
  if (demo) add('synthetic_demo',
    'Données synthétiques : aucune découverte biologique n’est validée.',
    'Synthetic data: no biological finding is validated.');
  if (exploratory) add('exploratory_only',
    'ACP et visualisations : structures descriptives, sans test d’hypothèse biologique.',
    'PCA and plots describe patterns, not tests of a biological hypothesis.');
  else {
    add('independent_confirmation',
      'Les associations et les q-values du moteur navigateur doivent être confirmées par la méthode statistique de référence adaptée au plan.',
      'Browser-model associations and q-values need confirmation with the design-appropriate reference statistical method.',
      'requires_confirmation');
    if (n < 20) add('sample_size', 
      'Effectif inférieur à 20 sujets : estimation de variance, puissance et modèles multivariés potentiellement instables.',
      'Fewer than 20 subjects: variance, power and multivariable models may be unstable.',
      'requires_confirmation');
  }
  const matched = result?.referenceBackend?.status === 'ok'
    ? Object.entries(result.referenceBackend.methods || {}).filter(([, method]) => method?.status === 'ok').map(([name])=>name)
    : [];
  if (!exploratory && matched.length === 0) add('no_reference_engine',
    'Aucune méthode R de référence n’a été exécutée : les résultats principaux restent exploratoires.',
    'No R reference method was executed: primary statistical results remain exploratory.',
    'requires_confirmation');
  if (matched.length) add('reference_not_certificate',
    'Méthodes R exécutées : ' + matched.join(', ') + '. Leur exécution ne prouve pas que le design et les hypothèses sont corrects.',
    'Executed R methods: ' + matched.join(', ') + '. Execution does not prove the study design or assumptions are valid.');
  for (const [layer, data] of layers) {
    const qc = data?.qc || {};
    // A near-equal percentage of missing values in the biological groups
    // is NOT evidence against MNAR: opposite censored tails can fabricate
    // apparent treatment effects without a group-wise missingness imbalance.
    // Inspect ANY measured profile, not only the median across profiles:
    // the median can be zero despite real feature-wise missingness.
    const hasObservedMissingness =
      (Array.isArray(qc.sampleMetrics) &&
        qc.sampleMetrics.some(profileHasMissingMeasurement))
      || Number(qc.medianMissingFraction) > 0;
    if (hasObservedMissingness) add(layer + '_missingness_not_ignorable',
      'Certaines mesures sont manquantes. Même avec des taux similaires entre groupes, cela peut créer de fausses différences biologiques : vérifier la cause des absences et refaire une analyse de sensibilité. Les p-values/q-values seules ne suffisent pas.',
      'Some measurements are missing. Similar missing-data rates between groups can still create false biological differences: investigate why values are absent and run a sensitivity analysis. P-values/q-values alone are insufficient.',
      'requires_confirmation');
    if (qc?.inferenceTier?.level === 'screening') add(layer + '_count_screening',
      'RNA-seq en comptages : le calcul navigateur sur log2-CPM ne remplace pas DESeq2/limma-voom sur comptages bruts.',
      'RNA-seq counts: browser log2-CPM screening is not a substitute for count-aware DESeq2/limma-voom.',
      'requires_confirmation');
    for (const warning of qc.warnings || []) {
      add(layer + '_qc_warning_' + notes.length,
        layer + ' : ' + warning, layer + ': ' + warning, 'requires_confirmation');
    }
    if (layer === 'metabolomics') {
      const ms = qc.msQc || {};
      if (ms.qcAssays === 0 && ms.biologicalAssays > 0) add('ms_no_pooled_qc',
        'Aucun pooled-QC MS renseigné : stabilité analytique et dérive instrumentale non vérifiables.',
        'No pooled MS QC annotated: instrumental stability and drift cannot be assessed.',
        'requires_confirmation');
      if (ms.driftCorrection?.applied) add('ms_drift_sensitivity',
        'Correction de dérive MS appliquée : comparer aussi aux résultats sans correction.',
        'MS drift correction applied: include a no-drift-correction sensitivity analysis.');
      if (ms.mnar?.imputedValues > 0) add('ms_mnar_sensitivity',
        'Valeurs manquantes MS imputées : vérifier la stabilité des conclusions sans imputation.',
        'MS missing values imputed: evaluate findings without imputation.',
        'requires_confirmation');
    }
  }
  if (protocol.objective === 'outcome' && protocol.outcomeType === 'survival') add('survival_browser_inference_blocked',
    'Survie : le moteur navigateur affiche uniquement les coefficients Cox (gestion Breslow des ex æquo), sans p-value, q-value ni intervalle de confiance. La validation confirmatoire exige un modèle R indépendant, la vérification de l’hypothèse de risques proportionnels et l’examen du mécanisme de censure.',
    'Survival: browser Cox reports exploratory coefficients only (Breslow ties), with no p/q-values or confidence intervals. Confirmatory inference requires independent R survival fitting, proportional-hazards diagnostics and review of censoring assumptions.',
    'requires_confirmation');
    if (protocol.objective === 'time' && protocol.longitudinal) add('longitudinal_browser_uncalibrated',
    'Le modèle longitudinal navigateur ne fournit plus de p-value, q-value ou intervalle de confiance : une simulation nulle a révélé une inflation possible des faux positifs. Les effets affichés restent descriptifs. Pour conclure, utilisez lmerTest sur le plan complet.',
    'The browser longitudinal model no longer reports p-values, q-values or confidence intervals: null simulation suggested possible type-I inflation. Effects are descriptive only. Use lmerTest with the full study design for inference.',
    'requires_confirmation');
  if (result?.predictiveOutcome?.status === 'ok') {
    if (result.predictiveOutcome.preprocessingLeakageRisk !== false) add('cv_preprocessing_limit',
      'La validation croisée reste exploratoire : au moins une omique peut subir un prétraitement global avant la séparation des sujets.',
      'Cross-validation remains exploratory: at least one modality may undergo global preprocessing before splitting subjects.',
      'requires_confirmation');
    else add('cv_fold_local_preparation',
      'Pour les échelles prises en charge, le moteur ne filtre plus les variables sur la cohorte entière avant la CV. Les prétraitements réalisés avant import, la conception de l’étude et la validation externe restent à vérifier.',
      'For supported scales, the engine no longer filters features on the full cohort before CV. Pre-import processing, study design and external validation still need review.',
      'requires_confirmation');
    if (result.predictiveOutcome.rawIntensityScreening) add('raw_intensity_cv_qc_limit',
      'Prédiction LFQ/MS : transformation fixée et normalisation individuelle sans fuite inter-sujets, mais correction de dérive pooled-QC, blancs et filtres RSD non transférés à ce modèle. Le même panel d’analytes doit être mesuré chez les futurs sujets. Validation technique et externe indispensable.',
      'LFQ/MS prediction: fixed transformation and assay-local normalization without cross-subject leakage, but pooled-QC drift, blank and RSD corrections were not transferred to this model. Future samples require the same analyte panel. Technical and external validation are essential.',
      'requires_confirmation');
  }
  if (layers.length > 1 && !exploratory) add('separate_fdr_families',
    'Les q-values sont corrigées par couche omique, pas globalement sur toutes les couches et analyses de voies : préspécifier les familles de tests.',
    'BH q-values are adjusted within each omics layer, not globally across omics and pathway tests: prespecify the tested families.',
    'requires_confirmation');
  if (result?.reactome) {
    add('pathway_hypotheses',
      'Les enrichissements Reactome génèrent des hypothèses, sans preuve causale. Une sélection préalable des variables sur les mêmes données empêche de les interpréter comme validation indépendante.',
      'Reactome enrichments generate hypotheses, not causal evidence. Feature selection on the same data prevents independent confirmatory interpretation.');
    const responses = [result.reactome.combined, ...Object.values(result.reactome.perLayer || {})];
    const incomplete = responses.filter((response) => response?.assayUniverse?.multiplicityStatus !== 'full_assay_universe_family');
    if (incomplete.length) add('pathway_fdr_not_estimable',
      'Au moins un enrichissement Reactome n’a pas d’univers de voies complet ou cohérent : sa FDR spécifique à l’étude n’est pas calculable, et la FDR Reactome d’origine ne peut pas la remplacer.',
      'At least one Reactome enrichment has an incomplete/inconsistent pathway universe: study-specific FDR is unavailable and cannot be replaced by Reactome default FDR.',
      'requires_confirmation');
    else add('pathway_family_validated',
      'FDR Reactome calculée sur toutes les voies associées à l’univers mesuré ; les dépendances entre voies et le biais de sélection des variables restent des limites.',
      'Reactome FDR is calculated over all pathways in the measured assay universe; pathway dependence and prior feature selection remain limitations.');
  }
  if (!layers.length) add('no_results', 'Aucun résultat statistique exploitable.', 'No usable statistical result.', 'requires_confirmation');
  const status = exploratory ? 'descriptive' : notes.some((item) => item.level === 'requires_confirmation')
    ? 'needs_reference_confirmation' : 'exploratory_inference';
  const confirmatoryReadiness = evaluateConfirmatoryReadiness(result, { demo });
  return {
    status,
    confirmatoryReadiness,
    publicationReady: false,
    referenceMethodsExecuted: matched,
    notes,
    version: 'scientific-assurance-v1',
    limitation: 'Automatic statistical quality checks are not external validation or a publication approval.'
  };
}
