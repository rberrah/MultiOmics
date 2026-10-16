(() => {
  if (!window.location.pathname.includes('/multiomics/tool')) return;

  const lang = () => document.documentElement.lang === 'en' ? 'en' : 'fr';
  const LEVEL_KEY = 'pmx-multiomics-metabolite-confidence';

  const definitions = {
    mixed_unknown: {
      fr: ['Mixte / non documenté', 'Le niveau varie entre variables ou n’est pas documenté. La conclusion doit suivre le niveau de confiance de chaque variable ; à défaut, les noms de métabolites et les voies restent exploratoires.'],
      en: ['Mixed / undocumented', 'Confidence varies across features or is undocumented. Conclusions should follow feature-level evidence; otherwise named-metabolite and pathway claims remain exploratory.']
    },
    msi1: {
      fr: ['MSI niveau 1 · composé identifié', 'Identité confirmée avec un standard authentique analysé dans des conditions comparables et au moins deux propriétés orthogonales. C’est le niveau adapté aux affirmations portant sur un métabolite nommé.'],
      en: ['MSI level 1 · identified compound', 'Identity confirmed with an authentic standard under comparable analytical conditions and at least two orthogonal properties. This is the appropriate level for named-metabolite claims.']
    },
    msi2: {
      fr: ['MSI niveau 2 · composé annoté putativement', 'Annotation sans standard authentique mesuré dans le même laboratoire, typiquement soutenue par une similarité spectrale/bibliothèque. Le nom du composé reste putatif et l’interprétation de voie est génératrice d’hypothèses.'],
      en: ['MSI level 2 · putatively annotated compound', 'Annotation without an authentic standard measured in the same laboratory, typically supported by spectral/library similarity. The compound name remains putative and pathway interpretation is hypothesis-generating.']
    },
    msi3: {
      fr: ['MSI niveau 3 · classe de composé putative', 'Les données soutiennent surtout une classe chimique. Une structure moléculaire unique ne doit pas être présentée comme identifiée ; les voies dépendant d’un métabolite précis ne doivent pas être considérées comme confirmées.'],
      en: ['MSI level 3 · putatively characterised compound class', 'Evidence mainly supports a chemical class. A unique molecular structure should not be reported as identified; pathways depending on one exact metabolite should not be treated as confirmed.']
    },
    msi4: {
      fr: ['MSI niveau 4 · variable inconnue', 'Le signal peut être quantifié et comparé, mais il ne doit pas être transformé en métabolite nommé ni en preuve de voie biologique tant qu’une annotation structurale suffisante n’est pas obtenue.'],
      en: ['MSI level 4 · unknown feature', 'The signal can be quantified and compared, but it should not be turned into a named metabolite or pathway-level claim until sufficient structural annotation is obtained.']
    }
  };

  function currentLevel() {
    const select = document.getElementById('pmx-metabolomics-msi-level');
    return select instanceof HTMLSelectElement && definitions[select.value] ? select.value : 'mixed_unknown';
  }

  function publishLevel(value) {
    window.__PMX_METABOLOMICS_IDENTIFICATION_CONFIDENCE__ = value;
  }

  function persist(value) {
    publishLevel(value);
    try { localStorage.setItem(LEVEL_KEY, value); } catch {}
  }

  function restore() {
    try {
      const value = localStorage.getItem(LEVEL_KEY);
      return value && definitions[value] ? value : 'mixed_unknown';
    } catch {
      return 'mixed_unknown';
    }
  }

  function interpretationPolicy(level, language) {
    const label = definitions[level]?.[language]?.[0] || definitions.mixed_unknown[language][0];
    const text = definitions[level]?.[language]?.[1] || definitions.mixed_unknown[language][1];
    const pathway = level === 'msi1'
      ? (language === 'en' ? 'Named-metabolite pathway interpretation is permitted, but still depends on statistical evidence and pathway coverage.' : 'L’interprétation des voies à partir de métabolites nommés est permise, mais dépend toujours de la preuve statistique et de la couverture des voies.')
      : level === 'msi2'
        ? (language === 'en' ? 'Pathway interpretation is hypothesis-generating because the metabolite identity is putative.' : 'L’interprétation de voie reste génératrice d’hypothèses car l’identité du métabolite est putative.')
        : (language === 'en' ? 'Exact named-metabolite contribution to integrated pathway claims is suppressed at this confidence level.' : 'La contribution d’un métabolite nommé exact aux conclusions de voies intégrées est supprimée à ce niveau de confiance.');
    return { label, text, pathway };
  }

  function render() {
    const platformOption = [...document.querySelectorAll('select option')].find((option) => option.value === 'untargeted_lcms');
    const platformSelect = platformOption?.parentElement;
    const card = platformSelect?.closest('article');
    const host = document.querySelector('[data-testid="multiomics-simple-study"]') || card;
    if (!host) return;

    let panel = host.querySelector('[data-testid="metabolomics-identification-confidence"]');
    if (!panel) {
      panel = document.createElement('div');
      panel.dataset.testid = 'metabolomics-identification-confidence';
      panel.className = 'metabolomics-identification-confidence';
      host.appendChild(panel);
    }

    const language = lang();
    const existingSelect = document.getElementById('pmx-metabolomics-msi-level');
    const value = existingSelect ? currentLevel() : restore();
    publishLevel(value);
    const policy = interpretationPolicy(value, language);
    const state = `${language}|${value}`;
    if (panel.dataset.state === state && panel.querySelector('select')) return;

    const wasOpen = panel.querySelector('details.msi-expert')?.open === true;
    panel.innerHTML = language === 'en'
      ? `
        <details class="msi-expert">
          <summary>Metabolite identification confidence</summary>
        <strong>Metabolite identification confidence <span class="msi-badge">MSI</span></strong>
        <label for="pmx-metabolomics-msi-level">Confidence supporting the metabolite names in this matrix</label>
        <select id="pmx-metabolomics-msi-level" name="metabolomics_identification_confidence">
          <option value="mixed_unknown">Mixed / undocumented</option>
          <option value="msi1">MSI 1 · identified with authentic standard</option>
          <option value="msi2">MSI 2 · putatively annotated compound</option>
          <option value="msi3">MSI 3 · putatively characterised compound class</option>
          <option value="msi4">MSI 4 · unknown feature</option>
        </select>
        <small><b>${policy.label}.</b> ${policy.text}</small>
        <small class="msi-pathway"><b>Pathway rule:</b> ${policy.pathway}</small>
        <small data-testid="multiomics-ms-qc-evidence-boundary"><b>Instrument QC:</b> Metabolite annotation is separate from pooled-QC precision, solvent blanks, injection-order drift and missing-value mechanisms. If QC/blank injections are absent, those properties are not verified. No automatic filtering proves a biomarker. Even when QC and blanks exist, a given molecular feature can be absent from those controls. Its precision or background cannot then be declared satisfactory.</small>
        <details><summary>Why this is separate from ChEBI/HMDB/KEGG mapping</summary><p>A database identifier says which database concept a label points to; it does not prove that the measured LC-MS/GC-MS feature was experimentally identified as that molecule. Identifier resolution and analytical identification confidence are therefore recorded separately.</p></details>
        </details>
      `
      : `
        <details class="msi-expert">
          <summary>Confiance d’identification des métabolites</summary>
        <strong>Confiance d’identification des métabolites <span class="msi-badge">MSI</span></strong>
        <label for="pmx-metabolomics-msi-level">Niveau de preuve qui soutient les noms de métabolites de cette matrice</label>
        <select id="pmx-metabolomics-msi-level" name="metabolomics_identification_confidence">
          <option value="mixed_unknown">Mixte / non documenté</option>
          <option value="msi1">MSI 1 · identifié avec standard authentique</option>
          <option value="msi2">MSI 2 · composé annoté putativement</option>
          <option value="msi3">MSI 3 · classe de composé putative</option>
          <option value="msi4">MSI 4 · variable inconnue</option>
        </select>
        <small><b>${policy.label}.</b> ${policy.text}</small>
        <small class="msi-pathway"><b>Règle pour les voies :</b> ${policy.pathway}</small>
        <small data-testid="multiomics-ms-qc-evidence-boundary"><b>Contrôle qualité instrumental :</b> L’identification des métabolites est distincte de la précision pooled-QC, des blancs, de la dérive selon l’ordre d’injection et des mécanismes de valeurs manquantes. Sans injections QC ou blancs, ces propriétés ne sont pas vérifiées. Un filtrage automatique ne valide pas un biomarqueur. Même si des QC et des blancs existent, une variable peut être absente de ces injections : sa précision ou son bruit de fond ne peuvent alors pas être déclarés satisfaisants.</small>
        <details><summary>Pourquoi ceci est séparé du mapping ChEBI/HMDB/KEGG</summary><p>Un identifiant de base indique vers quel concept pointe un libellé ; il ne prouve pas que le signal LC-MS/GC-MS mesuré a été expérimentalement identifié comme cette molécule. La résolution d’identifiant et la confiance d’identification analytique sont donc enregistrées séparément.</p></details>
        </details>
      `;
    if (wasOpen) panel.querySelector('details.msi-expert').open = true;

    const select = panel.querySelector('select');
    if (select instanceof HTMLSelectElement) {
      select.value = value;
      select.addEventListener('change', () => {
        persist(select.value);
        panel.dataset.state = '';
        render();
      }, { once: true });
    }
    panel.dataset.state = state;
  }

  // The engine dispatches this synchronously before returning the result. Mutating
  // detail here therefore adds the user-declared analytical identification confidence
  // to the same result object that is subsequently downloaded as JSON/HTML.
  window.addEventListener('pmx-multiomics-analysis', (event) => {
    const result = event?.detail;
    if (!result || typeof result !== 'object') return;
    const level = currentLevel();
    const policyEn = interpretationPolicy(level, 'en');
    const suppressExactMetabolitePathways = level === 'msi3' || level === 'msi4';

    result.protocol = {
      ...(result.protocol || {}),
      metabolomicsIdentificationConfidence: level
    };
    result.metabolomicsAnnotationConfidence = {
      scheme: 'Metabolomics Standards Initiative (MSI) four-level reporting framework',
      declaredLevel: level,
      label: policyEn.label,
      pathwayInterpretation: policyEn.pathway,
      exactPathwayMappingSuppressed: suppressExactMetabolitePathways,
      note: 'Database identifier resolution is not evidence of analytical compound identification. Apply per-feature confidence when mixed levels are present.'
    };

    // MSI 3 describes a putative class and MSI 4 an unknown feature. In either
    // case an exact ChEBI/HMDB/KEGG label must not contribute as if a unique
    // metabolite structure had been identified. Statistical feature results are
    // retained, but integrated exact-metabolite Reactome claims are suppressed.
    if (suppressExactMetabolitePathways && result.reactome) {
      const perLayer = { ...(result.reactome.perLayer || {}) };
      if ('metabolomics' in perLayer) perLayer.metabolomics = null;
      const previousCombined = result.reactome.combined || {};
      result.reactome = {
        ...result.reactome,
        combined: {
          ...previousCombined,
          token: null,
          pathwaysFound: 0,
          pathways: [],
          suppressed: true,
          suppressionReason: `Exact metabolite pathway integration suppressed because declared confidence is ${level}.`
        },
        consensus: [],
        perLayer,
        metabolomicsConfidenceSuppressed: true,
        metabolomicsConfidenceReason: `Exact metabolite pathway integration suppressed because declared confidence is ${level}.`,
        backgroundCaveat: [
          result.reactome.backgroundCaveat,
          'Integrated Reactome consensus was suppressed because metabolite confidence did not support exact compound-level mapping.'
        ].filter(Boolean).join(' ')
      };
    }
  });

  function installStyle() {
    if (document.getElementById('metabolomics-confidence-style')) return;
    const style = document.createElement('style');
    style.id = 'metabolomics-confidence-style';
    style.textContent = `
      .metabolomics-identification-confidence{display:grid;gap:7px;margin-top:12px;padding:12px;border:1px solid var(--border-subtle,#d4d4d4);border-radius:9px;background:var(--bg-secondary,#f7f7f7)}
      .metabolomics-identification-confidence>strong{display:flex;align-items:center;gap:7px}.metabolomics-identification-confidence>label{font-size:.86rem;font-weight:650}.metabolomics-identification-confidence>small{display:block;line-height:1.42;color:var(--text-secondary,#555)}
      .metabolomics-identification-confidence select{width:100%}.msi-badge{font-size:.7rem;border:1px solid var(--border-subtle,#ccc);border-radius:999px;padding:2px 6px}.msi-pathway{padding-top:6px;border-top:1px solid var(--border-subtle,#ddd)}
      .metabolomics-identification-confidence details{font-size:.83rem}.metabolomics-identification-confidence summary{cursor:pointer;font-weight:650}.metabolomics-identification-confidence details p{margin:6px 0 0;line-height:1.4;color:var(--text-secondary,#555)}
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

  const observer = new MutationObserver(refresh);
  observer.observe(document.documentElement, { childList: true, subtree: true, attributes: true, attributeFilter: ['lang'] });
  document.addEventListener('DOMContentLoaded', refresh, { once: true });
  refresh();
})();
