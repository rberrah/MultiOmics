(() => {
  if (!window.location.pathname.includes('/multiomics/tool')) return;

  let latestDiagnostics = window.__PMX_MULTIOMICS_DIAGNOSTICS__ || null;
  const lang = () => document.documentElement.lang === 'en' ? 'en' : 'fr';
  const clean = (value) => String(value ?? '').replace(/\s+/g, ' ').trim();
  // Diagnostics and table labels may contain strings imported from user files.
  // Only markup authored by PMx may be inserted with innerHTML.
  const escapeHtml = (value) => String(value ?? '').replace(/[&<>"']/g, (character) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
  })[character]);

  function selectContaining(value) {
    const option = [...document.querySelectorAll('select option')].find((item) => item.value === value);
    return option?.parentElement instanceof HTMLSelectElement ? option.parentElement : null;
  }

  const objective = () => selectContaining('explore')?.value || 'explore';
  const design = () => selectContaining('crossover')?.value || 'independent';
  const endpoint = () => selectContaining('survival')?.value || 'none';

  // Read ONLY the actual feature-wise inference grid. Reactome pathway FDR
  // belongs to a different hypothesis family and must never be pooled with
  // gene/protein/metabolite BH q-values in a novice-facing conclusion.
  function collectStatisticalEvidence() {
    const results = document.querySelector('[data-testid="multiomics-results"]');
    if (!results) return [];
    const evidence = [];
    for (const [index, grid] of [...results.querySelectorAll('[data-testid="multiomics-feature-inference"]')].entries()) {
      const header = grid.querySelector(':scope > .feature-head');
      if (!header || !/q\s*BH/i.test(header.textContent || '')) continue;
      let tested = 0;
      let significant = 0;
      for (const row of [...grid.children]) {
        if (row === header) continue;
        // The displayed three-digit value can round across a 0.05 boundary.
        // Only use the original full-precision numeric q-value from Svelte.
        const qText = row.lastElementChild?.getAttribute('data-scientific-q');
        if (qText == null || !qText.trim()) continue;
        const q = Number(qText);
        if (!Number.isFinite(q) || q < 0 || q > 1) continue;
        tested++;
        if (q <= 0.05) significant++;
      }
      if (!tested) continue;
      evidence.push({
        label: clean(grid.closest('article')?.querySelector('.layer-title .method-tag')?.textContent) || `Omic ${index + 1}`,
        qHeader: 'q BH',
        displayedRowsWithQ: tested,
        significantDisplayedRows: significant,
        // A fold ratio above 1 is not the same scale as a regression beta.
        // Never infer up/down direction from heterogeneous rendered effects.
        positiveSignificantDisplayedRows: 0,
        negativeSignificantDisplayedRows: 0
      });
    }
    return evidence;
  }

  function collectVisibleWarnings() {
    const results = document.querySelector('[data-testid="multiomics-results"]');
    if (!results) return [];
    const values = [];
    for (const node of results.querySelectorAll('.warning,.error,.alert,.caution,[role="alert"],.status.warn,.status.error')) {
      const text = clean(node.textContent);
      if (text && text.length <= 500) values.push(text);
    }
    return [...new Set(values)].slice(0, 20);
  }

  function translateDiagnostic(message, language) {
    const text = clean(message);
    if (!text || language === 'en') return text;
    const rules = [
      [/Group comparison was requested but fewer than two biological conditions are represented\.?/i, 'Une comparaison de groupes a été demandée, mais moins de deux conditions biologiques sont représentées.'],
      [/At least one biological group has fewer than two independent subjects[^.]*\.?/i, 'Au moins un groupe contient moins de deux sujets biologiquement indépendants : une inférence comparative défendable n’est pas estimable.'],
      [/At least one biological group has fewer than five independent subjects[^.]*\.?/i, 'Au moins un groupe contient moins de cinq sujets biologiquement indépendants : les résultats inférentiels doivent rester exploratoires et être lus avec les tailles d’effet et leur incertitude.'],
      [/One biological group contains less than half as many unique subjects as the largest group\.?/i, 'Un groupe contient moins de la moitié du nombre de sujets du groupe le plus grand.'],
      [/Group sizes are moderately unbalanced[^.]*\.?/i, 'Les tailles de groupes sont modérément déséquilibrées ; le modèle reste possible mais la précision peut différer entre groupes.'],
      [/Fewer than half of subjects are represented in every loaded omics layer\.?/i, 'Moins de la moitié des sujets sont présents dans toutes les couches omiques chargées.'],
      [/Not all subjects are represented in every omics layer[^.]*\.?/i, 'Tous les sujets ne sont pas présents dans chaque omique ; les analyses appariées utilisent uniquement le chevauchement réel.'],
      [/median missingness exceeds 30%\.?/i, 'la proportion médiane de données manquantes dépasse 30 %.'],
      [/median missingness is between 10% and 30%\.?/i, 'la proportion médiane de données manquantes est comprise entre 10 et 30 %.'],
      [/fewer than two subjects contain repeated time points[^.]*\.?/i, 'moins de deux sujets possèdent de véritables mesures répétées ; un effet populationnel longitudinal n’est pas estimable de façon défendable.'],
      [/At least one technical series is completely confounded with the biological design[^.]*\.?/i, 'Au moins une série technique est totalement confondue avec le plan biologique : l’effet biologique correspondant ne peut pas être séparé statistiquement.'],
      [/Technical-series annotation is incomplete in at least one omics layer\.?/i, 'L’annotation des séries techniques est incomplète dans au moins une couche omique.'],
      [/Survival analysis requires more usable subjects and at least two observed events\.?/i, 'L’analyse de survie nécessite davantage de sujets exploitables et au moins deux événements observés.'],
      [/Fewer than ten observed survival events are available[^.]*\.?/i, 'Moins de dix événements de survie sont observés : les associations de Cox et les performances prédictives internes peuvent être très instables.'],
      [/MOFA2 factor analysis requires more than 15 shared samples[^.]*\.?/i, 'MOFA2 nécessite plus de 15 sujets réellement communs aux couches pour être utilisé comme analyse principale défendable.'],
      [/Known technical series are confounded with the biological design\.?/i, 'Des séries techniques connues sont confondues avec le plan biologique.'],
      [/Technical-series annotation is incomplete\.?/i, 'L’annotation des séries techniques est incomplète.'],
      [/DIABLO is only applicable to a supervised categorical target in this tool\.?/i, 'DIABLO n’est applicable ici qu’à un critère catégoriel supervisé.'],
      [/Too few subjects are shared across all loaded omics blocks\.?/i, 'Trop peu de sujets sont communs à toutes les couches omiques chargées.'],
      [/At least three subjects per class are required for stratified cross-validation\.?/i, 'Au moins trois sujets par classe sont nécessaires pour une validation croisée stratifiée.'],
      [/Technical confounding must be resolved before supervised multi-omics discrimination\.?/i, 'La confusion technique doit être résolue avant une discrimination multi-omique supervisée.']
    ];
    for (const [pattern, replacement] of rules) {
      if (pattern.test(text)) return text.replace(pattern, replacement);
    }
    return text;
  }

  function confidenceAssessment(evidence, diagnostics, visibleWarnings) {
    const status = diagnostics?.status;
    if (status === 'blocked_or_redesign_required') return 'blocked';
    if (status === 'review_required') return 'caution';
    if (status === 'usable_with_cautions') return 'caution';
    if (status === 'ready') return evidence.length ? 'interpretable' : 'descriptive';
    return visibleWarnings.length ? 'caution' : evidence.length ? 'interpretable' : 'descriptive';
  }

  function confidenceCopy(language, status) {
    const copy = {
      fr: {
        interpretable: ['Interprétable', 'Aucun problème déterministe bloquant n’a été détecté. Cela ne démontre ni une puissance suffisante, ni une absence de biais, ni une validité externe.'],
        caution: ['Interpréter avec prudence', 'Le moteur a détecté au moins une limite qui doit rester visible dans la conclusion et le rapport.'],
        blocked: ['Interprétation forte bloquée', 'Le moteur a détecté un problème de design non compatible avec une conclusion biologique forte. Corrigez le plan ou l’annotation, ou limitez l’analyse à une description adaptée.'],
        descriptive: ['Exploratoire', 'Cette route ne fournit pas nécessairement de tests feature-by-feature avec FDR. Interprétez d’abord la structure, le contrôle qualité et les associations.']
      },
      en: {
        interpretable: ['Interpretable', 'No deterministic blocking issue was detected. This does not demonstrate adequate power, absence of bias or external validity.'],
        caution: ['Interpret with caution', 'The engine detected at least one limitation that must remain visible in the conclusion and report.'],
        blocked: ['Strong interpretation blocked', 'The engine detected a design problem incompatible with a strong biological conclusion. Correct the design/annotation or restrict the analysis to an appropriate descriptive route.'],
        descriptive: ['Exploratory', 'This route does not necessarily produce feature-wise FDR tests. Read structure, quality control and associations first.']
      }
    };
    return copy[language][status] || copy[language].descriptive;
  }

  function conclusionFor(language, currentObjective, evidence, diagnostics) {
    if (diagnostics?.status === 'blocked_or_redesign_required') {
      return language === 'en'
        ? 'The numerical outputs may still be useful for data checking, but the detected design problem prevents a defensible strong biological interpretation.'
        : 'Les sorties numériques peuvent rester utiles pour contrôler les données, mais le problème de design détecté empêche une interprétation biologique forte défendable.';
    }
    if (!evidence.length) {
      return language === 'en'
        ? 'No feature-wise adjusted-q/FDR rows are visible for this route. Read the multivariate structure, quality checks and pathway evidence instead of looking for a single “significant” variable.'
        : 'Aucune ligne de résultat avec q BH ajusté n’est visible dans ce parcours. Examinez la structure multivariée, les contrôles qualité et les voies biologiques, plutôt que de rechercher une seule variable « significative ».';
    }
    // Each omic has its own BH test family. Never pool apparent discovery
    // proportions across omics or imply the ten-row preview covers all tests.
    const layerCounts = evidence.map((item) =>
      `${item.label} : ${item.significantDisplayedRows}/${item.displayedRowsWithQ}`).join(' ; ');
    if (language === 'en') {
      const base = `Among the first feature rows displayed in each omic, q BH ≤ 0.05: ${layerCounts}. These are visible-row counts, not the discovery rate for all tested features. The main tables may use the looser exploratory q ≤ 0.10 threshold; BH corrections are performed per omic, not pooled here.`;
      if (currentObjective === 'groups') return `${base} Group associations are not proof of causality.`;
      if (currentObjective === 'outcome') return `${base} Prediction quality must be assessed on held-out subjects.`;
      if (currentObjective === 'time') return `${base} Temporal association is not proof of causality.`;
      return `${base} Explore covariation and data quality separately.`;
    }
    const base = `Parmi les premières variables affichées par omique, q BH ≤ 0,05 : ${layerCounts}. Il s’agit seulement des lignes visibles, pas du taux de découvertes pour toutes les variables testées. Les tableaux peuvent utiliser le seuil exploratoire q ≤ 0,10 ; chaque omique a sa propre correction BH, sans regroupement ici.`;
    if (currentObjective === 'groups') return `${base} Une association aux groupes n’est pas une preuve de causalité.`;
    if (currentObjective === 'outcome') return `${base} La qualité prédictive doit être évaluée hors échantillon.`;
    if (currentObjective === 'time') return `${base} Une association temporelle ne démontre pas la causalité.`;
    return `${base} Examinez séparément les relations entre omiques et la qualité des données.`;
  }

  function eligibilityLabel(status, language) {
    const map = language === 'en'
      ? { eligible: 'eligible', eligible_with_internal_cv: 'eligible + internal CV', not_recommended: 'not recommended', not_applicable: 'not applicable' }
      : { eligible: 'éligible', eligible_with_internal_cv: 'éligible + validation interne', not_recommended: 'non recommandé', not_applicable: 'non applicable' };
    return map[status] || status || '—';
  }

  function diagnosticsHtml(language, diagnostics) {
    if (!diagnostics) return '';
    const blockers = diagnostics.blockers || [];
    const warnings = diagnostics.warnings || [];
    const notes = diagnostics.notes || [];
    const statusTitle = language === 'en' ? 'Deterministic study-design decision' : 'Décision déterministe sur le plan d’étude';
    const statusMap = language === 'en'
      ? { ready: 'Ready', usable_with_cautions: 'Usable with cautions', review_required: 'Review required', blocked_or_redesign_required: 'Blocked / redesign required' }
      : { ready: 'Prêt', usable_with_cautions: 'Utilisable avec précautions', review_required: 'Revue requise', blocked_or_redesign_required: 'Bloqué / plan à corriger' };
    const itemList = (title, values, kind) => values.length
      ? `<div class="diagnostic-list ${kind}"><strong>${title}</strong><ul>${values.map((value) => `<li>${escapeHtml(translateDiagnostic(value, language))}</li>`).join('')}</ul></div>`
      : '';
    return `
      <div class="engine-diagnostics" data-testid="multiomics-engine-diagnostics">
        <div class="engine-diagnostics-head"><strong>${statusTitle}</strong><b>${escapeHtml(statusMap[diagnostics.status] || diagnostics.status)}</b></div>
        ${itemList(language === 'en' ? 'Blocking issues' : 'Problèmes bloquants', blockers, 'blocking')}
        ${itemList(language === 'en' ? 'Warnings' : 'Avertissements', warnings, 'warning')}
        ${itemList(language === 'en' ? 'Cautions / notes' : 'Précautions / notes', notes, 'note')}
      </div>`;
  }

  function methodEligibilityHtml(language, diagnostics) {
    const eligibility = diagnostics?.methodEligibility;
    if (!eligibility) return '';
    const cards = [
      ['MOFA2', eligibility.mofa2],
      ['DIABLO', eligibility.diablo]
    ].map(([name, method]) => {
      if (!method) return '';
      const reasons = (method.reasons || []).map((reason) => `<li>${escapeHtml(translateDiagnostic(reason, language))}</li>`).join('');
      const sampleInfo = name === 'DIABLO' && method.smallestClass != null
        ? (language === 'en' ? `Smallest class: ${method.smallestClass}` : `Plus petite classe : ${method.smallestClass}`)
        : (language === 'en' ? `Shared subjects: ${method.sharedSubjects ?? '—'}` : `Sujets communs : ${method.sharedSubjects ?? '—'}`);
      const interpretation = name === 'MOFA2'
        ? (language === 'en' ? 'Latent covariance factors are exploratory; they are not causal mechanisms or validated biomarkers.' : 'Les facteurs latents de covariance sont exploratoires ; ce ne sont ni des mécanismes causaux ni des biomarqueurs validés.')
        : (language === 'en' ? 'Internal repeated CV supports model tuning; independent external validation is still required before biomarker claims.' : 'La validation croisée interne répétée sert au tuning ; une validation externe indépendante reste nécessaire avant toute revendication de biomarqueur.');
      return `<article><div><strong>${name}</strong><b>${escapeHtml(eligibilityLabel(method.status, language))}</b></div><span>${escapeHtml(sampleInfo)}</span>${reasons ? `<ul>${reasons}</ul>` : ''}<p>${interpretation}</p></article>`;
    }).join('');
    return `<div class="method-eligibility" data-testid="multiomics-method-eligibility"><strong>${language === 'en' ? 'Advanced-method eligibility' : 'Éligibilité des méthodes avancées'}</strong><div>${cards}</div></div>`;
  }

  function methodologicalRules(language) {
    const rules = language === 'en'
      ? [
          ['Biological independence', 'Technical repeats are not counted as new biological subjects.'],
          ['Multiple testing', 'Feature-wise inference is read with adjusted q/FDR rather than raw p-values alone.'],
          ['Technical structure', 'A technical series completely confounded with biology cannot be statistically separated from biology.'],
          ['Cross-omics matching', 'Integration uses the actual subject overlap; missing layers are not silently invented.']
        ]
      : [
          ['Indépendance biologique', 'Les répétitions techniques ne sont pas comptées comme de nouveaux sujets biologiques.'],
          ['Tests multiples', 'L’inférence feature-by-feature se lit avec les q/FDR ajustés et non avec les seules p-values brutes.'],
          ['Structure technique', 'Une série technique totalement confondue avec la biologie ne peut pas être séparée statistiquement de la biologie.'],
          ['Appariement inter-omique', 'L’intégration utilise le chevauchement réel des sujets ; une couche absente n’est pas inventée.']
        ];
    if (objective() === 'outcome') rules.push(language === 'en'
      ? ['Prediction', 'Selection/tuning and performance assessment remain separated; performance is read on held-out or nested-CV data.']
      : ['Prédiction', 'Sélection/tuning et évaluation restent séparés ; la performance est lue sur les sujets hors-échantillon ou en validation croisée imbriquée.']);
    return rules;
  }

  function controlSnapshot() {
    const output = [];
    [...document.querySelectorAll('main select, main input, main textarea')].forEach((control, index) => {
      if (!(control instanceof HTMLInputElement || control instanceof HTMLSelectElement || control instanceof HTMLTextAreaElement)) return;
      if (control.type === 'file' || ['button','submit','reset'].includes(control.type)) return;
      const id = control.id || '';
      const label = id ? clean(document.querySelector(`label[for="${CSS.escape(id)}"]`)?.textContent) : '';
      let value = control.value;
      if (control instanceof HTMLInputElement && ['checkbox','radio'].includes(control.type)) value = control.checked;
      if (value === '' || value == null) return;
      output.push({ key: control.name || id || label || `control_${index + 1}`, label: label || null, value });
    });
    return output;
  }

  function fileSnapshot() {
    return [...document.querySelectorAll('main input[type="file"]')].flatMap((input, inputIndex) => {
      const id = input.id || '';
      const label = id ? clean(document.querySelector(`label[for="${CSS.escape(id)}"]`)?.textContent) : '';
      return [...(input.files || [])].map((file) => ({
        input: input.name || id || label || `file_input_${inputIndex + 1}`,
        name: file.name,
        size_bytes: file.size,
        last_modified_iso: Number.isFinite(file.lastModified) ? new Date(file.lastModified).toISOString() : null,
        type: file.type || null
      }));
    });
  }

  function buildManifest() {
    const evidence = collectStatisticalEvidence();
    return {
      schema: 'pmx-explain-multiomics-reproducibility-manifest/v2',
      generated_at: new Date().toISOString(),
      application: { name: 'PMx Explain Multi-omics', page: window.location.pathname, origin: window.location.origin, language: lang() },
      scientific_question: { objective: objective(), design: design(), endpoint_type: endpoint() },
      controls: controlSnapshot(),
      input_files: fileSnapshot(),
      pre_analysis_diagnostics: latestDiagnostics,
      displayed_statistical_evidence: evidence,
      visible_warnings: collectVisibleWarnings(),
      result_sections: [...document.querySelectorAll('[data-testid="multiomics-results"] h2,[data-testid="multiomics-results"] h3,[data-testid="multiomics-results"] h4')].map((node) => clean(node.textContent)).filter(Boolean),
      reproducibility_notes: [
        'This manifest records interface choices and local file metadata, not input-file contents.',
        'Retain immutable input files and cryptographic checksums in the study archive.',
        'A corrected q/FDR value is not a causal claim.',
        'Cross-validated performance is internal validation unless an independent external cohort was used.',
        'Advanced integration remains subject to method-specific eligibility, tuning and external validation.'
      ]
    };
  }

  function manifestMarkdown(manifest) {
    const lines = ['# PMx Explain multi-omics — reproducibility report','',`Generated: ${manifest.generated_at}`,`Objective: ${manifest.scientific_question.objective}`,`Design: ${manifest.scientific_question.design}`,`Endpoint type: ${manifest.scientific_question.endpoint_type}`,''];
    lines.push('## Pre-analysis decision','',`- Status: ${manifest.pre_analysis_diagnostics?.status || 'not available'}`);
    for (const value of manifest.pre_analysis_diagnostics?.blockers || []) lines.push(`- BLOCKER: ${value}`);
    for (const value of manifest.pre_analysis_diagnostics?.warnings || []) lines.push(`- WARNING: ${value}`);
    for (const value of manifest.pre_analysis_diagnostics?.notes || []) lines.push(`- NOTE: ${value}`);
    lines.push('','## Advanced-method eligibility','');
    for (const [name, method] of Object.entries(manifest.pre_analysis_diagnostics?.methodEligibility || {})) {
      if (name === 'global') continue;
      lines.push(`- ${name}: ${method.status}`);
      for (const reason of method.reasons || []) lines.push(`  - ${reason}`);
    }
    lines.push('','## Interface parameters','');
    for (const item of manifest.controls) lines.push(`- ${item.label || item.key}: ${String(item.value)}`);
    lines.push('','## Input files','');
    if (!manifest.input_files.length) lines.push('- No browser file metadata available (for example, the built-in demo may have been used).');
    for (const file of manifest.input_files) lines.push(`- ${file.input}: ${file.name} (${file.size_bytes} bytes; modified ${file.last_modified_iso || 'unknown'})`);
    lines.push('','## Displayed adjusted evidence','');
    if (!manifest.displayed_statistical_evidence.length) lines.push('- No visible feature-wise q/FDR table was detected for this route.');
    for (const table of manifest.displayed_statistical_evidence) lines.push(`- ${table.label}: ${table.significantDisplayedRows}/${table.displayedRowsWithQ} displayed rows with ${table.qHeader} ≤ 0.05.`);
    lines.push('','## Reproducibility requirements','');
    for (const note of manifest.reproducibility_notes) lines.push(`- ${note}`);
    lines.push('','## Interpretation boundary','','Association, prediction and causality are distinct claims. Pathway enrichment does not by itself prove pathway activation. Internal cross-validation does not establish external validity.');
    return lines.join('\n');
  }

  function download(name, text, type) {
    const href = URL.createObjectURL(new Blob([text], { type }));
    const anchor = document.createElement('a');
    anchor.href = href;
    anchor.download = name;
    document.body.appendChild(anchor);
    anchor.click();
    anchor.remove();
    setTimeout(() => URL.revokeObjectURL(href), 1000);
  }

  function render() {
    const results = document.querySelector('[data-testid="multiomics-results"]');
    if (!results) return;
    latestDiagnostics = window.__PMX_MULTIOMICS_DIAGNOSTICS__ || latestDiagnostics;
    const anchor = document.querySelector('[data-testid="multiomics-result-reading-guide"]') || document.querySelector('[data-testid="multiomics-plain-readiness"]') || results.firstElementChild || results;
    let panel = document.querySelector('[data-testid="multiomics-scientific-summary"]');
    if (!panel) {
      panel = document.createElement('section');
      panel.className = 'scientific-summary';
      panel.dataset.testid = 'multiomics-scientific-summary';
      anchor.insertAdjacentElement('afterend', panel);
    }

    const language = lang();
    const evidence = collectStatisticalEvidence();
    const visibleWarnings = collectVisibleWarnings();
    const confidence = confidenceAssessment(evidence, latestDiagnostics, visibleWarnings);
    const [confidenceTitle, confidenceText] = confidenceCopy(language, confidence);
    const conclusion = conclusionFor(language, objective(), evidence, latestDiagnostics);
    const diagnostics = diagnosticsHtml(language, latestDiagnostics);
    const eligibility = methodEligibilityHtml(language, latestDiagnostics);
    const rules = methodologicalRules(language).map(([title, text]) => `<article><strong>${title}</strong><span>${text}</span></article>`).join('');
    const evidenceHtml = evidence.length
      ? evidence.map((item) => `<li><strong>${escapeHtml(item.label)}</strong><span>${item.significantDisplayedRows}/${item.displayedRowsWithQ} ${language === 'en' ? 'displayed rows with' : 'lignes affichées avec'} ${escapeHtml(item.qHeader)} ≤ 0.05</span></li>`).join('')
      : `<li><strong>${language === 'en' ? 'No visible feature-wise FDR table' : 'Aucun tableau FDR feature-by-feature visible'}</strong><span>${language === 'en' ? 'This can be normal for exploratory multivariate routes.' : 'Cela peut être normal pour une route multivariée exploratoire.'}</span></li>`;

    const state = JSON.stringify({ language, objective: objective(), confidence, diagnostics: latestDiagnostics, evidence });
    if (panel.dataset.state === state) return;
    const wasOpen = panel.querySelector('details.scientific-extra')?.open === true;
    panel.innerHTML = `
      <div class="scientific-summary-head"><div><span>${language === 'en' ? 'Scientific synthesis' : 'Synthèse scientifique'}</span><h3>${language === 'en' ? 'What do these results actually support?' : 'Que permettent réellement de conclure ces résultats ?'}</h3></div><b class="confidence ${confidence}">${escapeHtml(confidenceTitle)}</b></div>
      <p class="scientific-conclusion">${escapeHtml(conclusion)}</p>
      <p class="scientific-confidence">${escapeHtml(confidenceText)}</p>
      <details class="scientific-extra"><summary>${language === 'en' ? 'Scientific checks and reproducible exports' : 'Vérifications et exports scientifiques'}</summary>
      ${diagnostics}
      ${eligibility}
      <div class="scientific-summary-grid">
        <article><span>${language === 'en' ? 'Adjusted evidence visible in the results' : 'Preuves corrigées visibles dans les résultats'}</span><ul>${evidenceHtml}</ul></article>
        <article><span>${language === 'en' ? 'Interpretation boundary' : 'Frontière d’interprétation'}</span><p>${language === 'en' ? 'Association ≠ prediction ≠ causality. Pathway enrichment organises evidence; it does not prove pathway activation. Internal cross-validation does not replace validation in an independent cohort.' : 'Association ≠ prédiction ≠ causalité. Un enrichissement de voie organise les preuves ; il ne prouve pas l’activation de la voie. Une validation croisée interne ne remplace pas une validation dans une cohorte indépendante.'}</p></article>
      </div>
      <details class="methodology-checklist"><summary>${language === 'en' ? 'Why the statistical route is defensible' : 'Pourquoi la route statistique est défendable'}</summary><div>${rules}</div></details>
      <div class="reproducibility-export" data-testid="multiomics-methods-manifest"><div><strong>${language === 'en' ? 'Reproduce and report the analysis' : 'Reproduire et rapporter l’analyse'}</strong><span>${language === 'en' ? 'Export scientific choices, deterministic design diagnostics, local file metadata, visible adjusted evidence and interpretation limits. Input data are not embedded.' : 'Exportez les choix scientifiques, les diagnostics déterministes du plan, les métadonnées locales des fichiers, les preuves corrigées visibles et les limites d’interprétation. Les données ne sont pas intégrées.'}</span></div><div class="reproducibility-actions"><button type="button" data-export-method-report>${language === 'en' ? 'Methods report (.md)' : 'Rapport méthodes (.md)'}</button><button type="button" data-export-manifest-json>${language === 'en' ? 'Manifest (.json)' : 'Manifeste (.json)'}</button></div></div>
      </details>`;
    if (wasOpen) panel.querySelector('details.scientific-extra').open = true;
    panel.dataset.state = state;

    panel.querySelector('[data-export-manifest-json]')?.addEventListener('click', () => {
      const manifest = buildManifest();
      download('multiomics-reproducibility-manifest.json', JSON.stringify(manifest, null, 2), 'application/json;charset=utf-8');
    });
    panel.querySelector('[data-export-method-report]')?.addEventListener('click', () => {
      const manifest = buildManifest();
      download('multiomics-methods-report.md', manifestMarkdown(manifest), 'text/markdown;charset=utf-8');
    });
  }

  function installStyle() {
    if (document.getElementById('multiomics-scientific-summary-style')) return;
    const style = document.createElement('style');
    style.id = 'multiomics-scientific-summary-style';
    style.textContent = `
      .scientific-extra{margin-top:10px}.scientific-extra>summary{cursor:pointer;font-size:.88rem;font-weight:700}
      .scientific-summary{margin:18px 0;padding:18px;border:1px solid var(--border-strong,#c9c9c9);border-radius:12px;background:var(--bg-primary,#fff)}.scientific-summary-head{display:flex;justify-content:space-between;align-items:flex-start;gap:16px}.scientific-summary-head>div>span{display:block;font-size:.75rem;text-transform:uppercase;letter-spacing:.08em;color:var(--text-secondary,#666)}.scientific-summary-head h3{margin:4px 0 0;font-size:1.1rem}.scientific-summary .confidence{border:1px solid var(--border-subtle,#ccc);border-radius:999px;padding:5px 9px;font-size:.8rem;white-space:nowrap}.scientific-summary .confidence.blocked{font-weight:800;border-width:2px}.scientific-conclusion{font-size:1rem;line-height:1.55;margin:14px 0 6px}.scientific-confidence{font-size:.9rem;color:var(--text-secondary,#555);line-height:1.5;margin:0 0 14px}
      .engine-diagnostics,.method-eligibility{margin:12px 0;padding:12px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-secondary,#f7f7f7)}.engine-diagnostics-head{display:flex;justify-content:space-between;gap:10px}.diagnostic-list{margin-top:9px}.diagnostic-list>strong{font-size:.86rem}.diagnostic-list ul{margin:5px 0 0;padding-left:20px}.diagnostic-list li{margin:4px 0;font-size:.88rem;line-height:1.42}.diagnostic-list.blocking{border-left:3px solid currentColor;padding-left:10px}.method-eligibility>strong{display:block;margin-bottom:8px}.method-eligibility>div{display:grid;grid-template-columns:1fr 1fr;gap:8px}.method-eligibility article{padding:10px;border:1px solid var(--border-subtle,#ddd);border-radius:8px;background:var(--bg-primary,#fff)}.method-eligibility article>div{display:flex;justify-content:space-between;gap:8px}.method-eligibility article>span{display:block;margin-top:5px;font-size:.82rem;color:var(--text-secondary,#555)}.method-eligibility article ul{padding-left:18px;margin:7px 0}.method-eligibility article p{margin:7px 0 0;font-size:.84rem;line-height:1.4;color:var(--text-secondary,#555)}
      .scientific-summary-grid{display:grid;grid-template-columns:1fr 1fr;gap:10px}.scientific-summary-grid>article{padding:12px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-secondary,#f7f7f7)}.scientific-summary-grid>article>span{font-weight:700}.scientific-summary-grid p{margin:7px 0 0;font-size:.89rem;line-height:1.48}.scientific-summary-grid ul{margin:8px 0 0;padding-left:18px}.scientific-summary-grid li{margin:0 0 6px}.scientific-summary-grid li strong,.scientific-summary-grid li span{display:block}.scientific-summary-grid li span{font-size:.84rem;color:var(--text-secondary,#555)}.scientific-summary details{margin-top:11px;border-top:1px solid var(--border-subtle,#ddd);padding-top:10px}.scientific-summary summary{cursor:pointer;font-weight:700}.methodology-checklist>div{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:8px;margin-top:9px}.methodology-checklist article{padding:10px;border:1px solid var(--border-subtle,#ddd);border-radius:8px}.methodology-checklist article strong,.methodology-checklist article span{display:block}.methodology-checklist article span{margin-top:4px;font-size:.85rem;line-height:1.42;color:var(--text-secondary,#555)}
      .reproducibility-export{display:flex;justify-content:space-between;align-items:center;gap:16px;margin-top:14px;padding:12px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-secondary,#f7f7f7)}.reproducibility-export>div:first-child strong,.reproducibility-export>div:first-child span{display:block}.reproducibility-export>div:first-child span{margin-top:4px;max-width:70ch;font-size:.85rem;line-height:1.4;color:var(--text-secondary,#555)}.reproducibility-actions{display:flex;gap:8px;flex-wrap:wrap}.reproducibility-actions button{white-space:nowrap}
      @media(max-width:820px){.scientific-summary-head,.reproducibility-export{flex-direction:column;align-items:stretch}.scientific-summary-grid,.methodology-checklist>div,.method-eligibility>div{grid-template-columns:1fr}.scientific-summary .confidence{white-space:normal;align-self:flex-start}}
    `;
    document.head.appendChild(style);
  }

  let queued = false;
  function refresh() {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      installStyle();
      render();
    });
  }

  window.addEventListener('pmx-multiomics-diagnostics', (event) => {
    latestDiagnostics = event.detail || latestDiagnostics;
    refresh();
  });
  const observer = new MutationObserver(refresh);
  observer.observe(document.documentElement, { childList: true, subtree: true, attributes: true, attributeFilter: ['lang'] });
  document.addEventListener('change', refresh);
  document.addEventListener('DOMContentLoaded', refresh, { once: true });
  refresh();
})();
