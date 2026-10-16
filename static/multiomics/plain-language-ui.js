(() => {
  const TOOL_PATH = '/multiomics/tool';
  const isTool = () => window.location.pathname.includes(TOOL_PATH);
  if (!isTool()) return;

  const translations = {
    fr: new Map([
      ['Audit des batches techniques', 'Vérification des séries techniques'],
      ['Structure des batches vérifiée avant l’inférence biologique', 'Vérification que les séries techniques ne se confondent pas avec les groupes ou les temps'],
      ['Batch technique', 'Série technique'],
      ['Batch et covariables sélectionnées entrent directement dans les modèles d’outcome pour :', 'Les séries techniques et facteurs d’ajustement choisis sont pris en compte directement dans les modèles pour :'],
      ['Une résidualisation OLS a été appliquée avant l’inférence pour :', 'Les différences techniques déclarées ont été corrigées avant la comparaison pour :'],
      ['Ajustement technique', 'Correction des différences techniques'],
      ['missing médian', 'données manquantes (médiane)'],
      ['assays suspects', 'mesures à vérifier'],
      ['réplicats r<0,80', 'répétitions techniques peu concordantes'],
      ['QC par couche', 'contrôle par type de données'],
      ['QC MS avancé', 'Contrôle qualité MS avancé'],
      ['Missing MNAR', 'Valeurs probablement sous la limite de détection (MNAR)'],
      ['Outcome / critère', 'Critère étudié'],
      ['Batch technique', 'Série technique'],
      ['Covariable explicite', 'Facteur d’ajustement'],
      ['Métadonnées échantillons', 'Tableau des échantillons'],
      ['Template générique de métadonnées', 'Modèle générique du tableau des échantillons'],
      ['Dictionnaire détaillé des métadonnées', 'Détail des colonnes du tableau des échantillons'],
      ['Faites correspondre métadonnées et matrices avant l’analyse', 'Faites correspondre le tableau des échantillons et les matrices avant l’analyse']
    ]),
    en: new Map([
      ['Technical batch audit', 'Technical-series check'],
      ['Batch structure checked before biological inference', 'Check that technical series are not confounded with groups or time points'],
      ['Technical adjustment', 'Correction for technical differences'],
      ['median missing', 'missing data (median)'],
      ['flagged assays', 'measurements to review'],
      ['replicates r<0.80', 'poorly concordant technical repeats'],
      ['per-layer QC', 'quality control by data type'],
      ['Sample metadata', 'Sample sheet'],
      ['Detailed metadata dictionary', 'Sample-sheet column guide']
    ])
  };

  function lang() {
    return document.documentElement.lang === 'en' ? 'en' : 'fr';
  }

  function replaceExactText(root = document) {
    const map = translations[lang()];
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT);
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    for (const node of nodes) {
      const raw = node.nodeValue || '';
      const trimmed = raw.trim();
      if (!trimmed) continue;
      const replacement = map.get(trimmed);
      if (!replacement) continue;
      const leading = raw.match(/^\s*/)?.[0] || '';
      const trailing = raw.match(/\s*$/)?.[0] || '';
      node.nodeValue = leading + replacement + trailing;
    }
  }

  function glossaryHtml(language) {
    return language === 'en'
      ? `<summary>Quick glossary · six terms used throughout the tool</summary>
         <div class="plain-glossary-grid">
           <article><strong>Sample sheet</strong><span>Describes who, which specimen, group, time point and technical measurement. Also called metadata.</span></article>
           <article><strong>Matrix</strong><span>The table of measured RNA, protein or metabolite values.</span></article>
           <article><strong>Technical series</strong><span>Plate, run, day or analytical batch that may create non-biological differences.</span></article>
           <article><strong>Adjustment factor</strong><span>A variable such as age, sex or centre that you explicitly ask the model to account for.</span></article>
           <article><strong>Adjusted result (q/FDR)</strong><span>A statistical result corrected because many genes, proteins or metabolites were tested at once.</span></article>
           <article><strong>Cross-validation</strong><span>Tests prediction on subjects that were not used to build the model.</span></article>
         </div>`
      : `<summary>Lexique express · six mots utilisés dans tout l’outil</summary>
         <div class="plain-glossary-grid">
           <article><strong>Tableau des échantillons</strong><span>Décrit qui a été étudié, quel prélèvement, quel groupe, quel temps et quelle mesure. C’est ce que certains logiciels appellent « metadata ».</span></article>
           <article><strong>Matrice</strong><span>Tableau contenant les valeurs RNA, protéines ou métabolites réellement mesurées.</span></article>
           <article><strong>Série technique</strong><span>Plaque, run, jour ou lot d’analyse pouvant créer des différences qui ne sont pas biologiques.</span></article>
           <article><strong>Facteur d’ajustement</strong><span>Variable comme l’âge, le sexe ou le centre que vous demandez explicitement au modèle de prendre en compte.</span></article>
           <article><strong>Résultat corrigé (q/FDR)</strong><span>Résultat statistique corrigé parce que beaucoup de gènes, protéines ou métabolites sont testés en même temps.</span></article>
           <article><strong>Validation croisée</strong><span>Évalue une prédiction sur des sujets qui n’ont pas servi à construire le modèle.</span></article>
         </div>`;
  }

  function renderQuickGlossary() {
    const guide = document.querySelector('[data-testid="multiomics-sample-sheet-guide"]');
    if (!guide) return;
    const language = lang();
    let glossary = document.querySelector('[data-testid="multiomics-quick-glossary"]');
    if (!glossary) {
      glossary = document.createElement('details');
      glossary.className = 'plain-glossary';
      glossary.dataset.testid = 'multiomics-quick-glossary';
      guide.insertAdjacentElement('afterend', glossary);
    }
    if (glossary.dataset.language !== language) {
      const wasOpen = glossary.hasAttribute('open');
      glossary.innerHTML = glossaryHtml(language);
      glossary.dataset.language = language;
      if (wasOpen) glossary.setAttribute('open', '');
    }
  }

  function statusLabel(raw, language) {
    const value = String(raw || '').trim().toLowerCase().replaceAll('_', ' ');
    const labels = language === 'en'
      ? {
          'single batch': 'one technical series',
          'multiple batches adjusted': 'multiple technical series accounted for',
          'not provided': 'not provided',
          'incomplete': 'incomplete information',
          'confounded': 'confounded with the biological design'
        }
      : {
          'single batch': 'une seule série technique',
          'multiple batches adjusted': 'plusieurs séries techniques prises en compte',
          'not provided': 'information non fournie',
          'incomplete': 'information incomplète',
          'confounded': 'confondue avec le plan biologique'
        };
    return labels[value] || value;
  }

  function extractSummary() {
    const cards = [...document.querySelectorAll('.computed-summary article')];
    const values = {};
    for (const card of cards) {
      const key = card.querySelector('span')?.textContent?.trim() || '';
      const value = card.querySelector('strong')?.textContent?.trim() || '';
      if (key) values[key] = value;
    }
    return values;
  }

  function extractQc() {
    return [...document.querySelectorAll('[data-testid="multiomics-qc"] .qc-card')].map((card) => {
      const layer = card.querySelector('.qc-card-head strong')?.textContent?.trim() || '';
      const metrics = [...card.querySelectorAll('.qc-metrics > div')].map((item) => ({
        value: item.querySelector('b')?.textContent?.trim() || '',
        label: item.querySelector('span')?.textContent?.trim() || ''
      }));
      return { layer, metrics };
    });
  }

  function renderReadiness() {
    const results = document.querySelector('[data-testid="multiomics-results"]');
    const summary = document.querySelector('.computed-summary');
    if (!results || !summary) return;

    let panel = document.querySelector('[data-testid="multiomics-plain-readiness"]');
    if (!panel) {
      panel = document.createElement('section');
      panel.dataset.testid = 'multiomics-plain-readiness';
      panel.className = 'plain-readiness';
      summary.insertAdjacentElement('beforebegin', panel);
    }

    const language = lang();
    const summaryValues = extractSummary();
    const overlapStrong = document.querySelector('.overlap-box strong')?.textContent?.trim() || '';
    const batchBoxes = [...document.querySelectorAll('.overlap-box')];
    const batchBox = batchBoxes.find((box) => /batch|série technique|technical series/i.test(box.textContent || ''));
    const batchRows = batchBox
      ? [...batchBox.querySelectorAll('.overlap-pairs span')].map((row) => {
          const text = row.textContent || '';
          const [layer, rest = ''] = text.split(':');
          const rawStatus = row.querySelector('b')?.textContent || rest;
          return `${layer.trim()}: ${statusLabel(rawStatus, language)}`;
        })
      : [];
    const qc = extractQc();

    const subjects = summaryValues[language === 'en' ? 'Subjects' : 'Sujets'] || '—';
    const conditions = summaryValues.Conditions || '—';
    const timepoints = summaryValues[language === 'en' ? 'Time points' : 'Temps'] || '—';

    const qcText = qc.length
      ? qc.map(({ layer, metrics }) => {
          const missing = metrics.find((m) => /missing|manquant/i.test(m.label));
          const flagged = metrics.find((m) => /review|vérifier|suspect/i.test(m.label));
          return `${layer}: ${missing?.value || '—'} ${language === 'en' ? 'missing (median)' : 'de données manquantes (médiane)'}, ${flagged?.value || '0'} ${language === 'en' ? 'measurement(s) to review' : 'mesure(s) à vérifier'}`;
        }).join('<br>')
      : (language === 'en' ? 'Quality-control details appear below.' : 'Le détail du contrôle qualité apparaît ci-dessous.');

    const title = language === 'en' ? 'Check before interpreting the results' : 'Vérification avant d’interpréter les résultats';
    const intro = language === 'en'
      ? 'These checks describe whether the uploaded study structure is coherent enough to interpret. They are not a substitute for a formal power calculation.'
      : 'Ces vérifications indiquent si la structure des données chargées est suffisamment cohérente pour être interprétée. Elles ne remplacent pas un calcul formel de puissance.';
    const power = language === 'en'
      ? 'Statistical power is not estimated from sample size alone. A defensible calculation requires an expected effect, variability and the intended statistical model; simulation is preferable for complex multi-omics designs.'
      : 'La puissance statistique n’est pas déduite du seul nombre de sujets. Un calcul défendable nécessite un effet attendu, une variabilité et le modèle statistique prévu ; pour un plan multi-omique complexe, une simulation est préférable.';

    const html = `
      <div class="plain-readiness-head">
        <div><span>${language === 'en' ? 'Study structure' : 'Structure de l’étude'}</span><h3>${title}</h3></div>
        <b>${subjects} ${language === 'en' ? 'subjects' : 'sujets'}</b>
      </div>
      <p>${intro}</p>
      <div class="plain-readiness-grid">
        <article><strong>${language === 'en' ? 'Groups and time' : 'Groupes et temps'}</strong><span>${language === 'en' ? 'Groups' : 'Groupes'}: ${conditions}<br>${language === 'en' ? 'Time points' : 'Temps'}: ${timepoints}</span></article>
        <article><strong>${language === 'en' ? 'Overlap between omics' : 'Chevauchement entre omiques'}</strong><span>${overlapStrong || (language === 'en' ? 'See overlap details below.' : 'Voir le détail du chevauchement ci-dessous.')}</span></article>
        <article><strong>${language === 'en' ? 'Technical series' : 'Séries techniques'}</strong><span>${batchRows.length ? batchRows.join('<br>') : (language === 'en' ? 'No blocking technical-series issue is displayed.' : 'Aucun problème bloquant de série technique n’est affiché.')}</span></article>
        <article><strong>${language === 'en' ? 'Missing data and measurements to review' : 'Données manquantes et mesures à vérifier'}</strong><span>${qcText}</span></article>
      </div>
      <details><summary>${language === 'en' ? 'About statistical power' : 'À propos de la puissance statistique'}</summary><p>${power}</p></details>
    `;
    if (panel.innerHTML !== html) panel.innerHTML = html;
  }

  function installStyle() {
    if (document.getElementById('multiomics-plain-language-style')) return;
    const style = document.createElement('style');
    style.id = 'multiomics-plain-language-style';
    style.textContent = `
      .plain-glossary{max-width:1180px;margin:10px auto 20px;padding:0 18px}.plain-glossary>summary{cursor:pointer;font-weight:700;padding:11px 13px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-secondary,#f7f7f7)}
      .plain-glossary-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:9px;margin-top:9px}.plain-glossary-grid article{padding:11px;border:1px solid var(--border-subtle,#ddd);border-radius:9px;background:var(--bg-primary,#fff)}.plain-glossary-grid strong,.plain-glossary-grid span{display:block}.plain-glossary-grid span{margin-top:5px;font-size:.88rem;line-height:1.4;color:var(--text-secondary,#555)}
      .plain-readiness{margin:18px 0;padding:18px;border:1px solid var(--border-strong,#c9c9c9);background:var(--bg-secondary,#f7f7f7);border-radius:12px}
      .plain-readiness-head{display:flex;align-items:flex-start;justify-content:space-between;gap:16px;margin-bottom:8px}
      .plain-readiness-head span{font-size:.76rem;text-transform:uppercase;letter-spacing:.08em;color:var(--text-secondary,#666)}
      .plain-readiness-head h3{margin:4px 0 0;font-size:1.08rem}
      .plain-readiness-head b{white-space:nowrap;border:1px solid var(--border-subtle,#ddd);padding:5px 8px;border-radius:999px;font-size:.82rem}
      .plain-readiness>p{margin:8px 0 14px;color:var(--text-secondary,#555);max-width:80ch}
      .plain-readiness-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px}
      .plain-readiness-grid article{padding:12px;border:1px solid var(--border-subtle,#ddd);background:var(--bg-primary,#fff);border-radius:9px}
      .plain-readiness-grid strong,.plain-readiness-grid span{display:block}.plain-readiness-grid span{margin-top:6px;font-size:.9rem;line-height:1.45;color:var(--text-secondary,#555)}
      .plain-readiness details{margin-top:12px}.plain-readiness summary{cursor:pointer;font-weight:650}.plain-readiness details p{margin:8px 0 0;color:var(--text-secondary,#555)}
      @media(max-width:760px){.plain-glossary-grid,.plain-readiness-grid{grid-template-columns:1fr}.plain-readiness-head{flex-direction:column}.plain-readiness-head b{white-space:normal}}
    `;
    document.head.appendChild(style);
  }

  let scheduled = false;
  function refresh() {
    if (scheduled) return;
    scheduled = true;
    requestAnimationFrame(() => {
      scheduled = false;
      if (!isTool()) return;
      installStyle();
      replaceExactText(document.body);
      renderQuickGlossary();
      renderReadiness();
    });
  }

  const observer = new MutationObserver(refresh);
  observer.observe(document.documentElement, {
    childList: true,
    subtree: true,
    characterData: true,
    attributes: true,
    attributeFilter: ['lang']
  });
  window.addEventListener('popstate', refresh);
  document.addEventListener('DOMContentLoaded', refresh, { once: true });
  refresh();
})();
