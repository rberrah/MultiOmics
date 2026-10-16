(() => {
  if (!window.location.pathname.includes('/multiomics/tool')) return;

  const language = () => document.documentElement.lang === 'en' ? 'en' : 'fr';

  function findSelectWithOption(value) {
    const option = [...document.querySelectorAll('select option')].find((item) => item.value === value);
    return option?.parentElement instanceof HTMLSelectElement ? option.parentElement : null;
  }

  function currentObjective() {
    return findSelectWithOption('explore')?.value || 'explore';
  }

  function currentOutcomeType() {
    return findSelectWithOption('survival')?.value || 'none';
  }

  function currentDesign() {
    return findSelectWithOption('crossover')?.value || 'independent';
  }

  function objectiveCopy(lang, objective) {
    const copy = {
      fr: {
        explore: {
          question: 'Qu’est-ce qui varie ensemble entre mes omiques ?',
          short: 'Explorer sans imposer de groupe ou de critère cible.',
          action: 'L’outil recherche des structures communes entre RNA, protéines et métabolites, après avoir équilibré les couches pour qu’une omique ne domine pas uniquement parce qu’elle contient plus de variables.',
          science: 'ACP multi-blocs équilibrée sur variables standardisées, puis relations inter-omiques et annotation biologique.'
        },
        groups: {
          question: 'Qu’est-ce qui diffère entre mes groupes ?',
          short: 'Comparer contrôle/traitement, répondeurs/non-répondeurs ou plusieurs conditions.',
          action: 'L’outil analyse d’abord chaque omique avec le plan expérimental déclaré, corrige les tests multiples, puis cherche les signaux cohérents entre omiques et les voies biologiques communes.',
          science: 'Contrastes ajustés au design, correction BH-FDR, associations inter-omiques et convergence Reactome.'
        },
        outcome: {
          question: 'Qu’est-ce qui est associé à mon critère clinique ou expérimental ?',
          short: 'Relier les omiques à une réponse, un score, un événement ou une survie.',
          action: 'L’outil choisit le modèle adapté au type de critère, ajuste les facteurs déclarés et, lorsque les données le permettent, vérifie la capacité prédictive sur des sujets laissés de côté.',
          science: 'Régression adaptée au critère (linéaire, logistique, Poisson, ANCOVA ou Cox) avec BH-FDR et validation croisée imbriquée lorsqu’elle est faisable.'
        },
        time: {
          question: 'Qu’est-ce qui change au cours du temps ?',
          short: 'Étudier des visites répétées ou des trajectoires.',
          action: 'L’outil tient compte du fait que plusieurs mesures appartiennent au même sujet et teste si l’évolution diffère selon le groupe ou la condition.',
          science: 'Modèle longitudinal avec intercept aléatoire sujet et interaction condition × temps, avec séries techniques et facteurs d’ajustement explicites.'
        }
      },
      en: {
        explore: {
          question: 'What varies together across my omics layers?',
          short: 'Explore without imposing a group or target endpoint.',
          action: 'The tool searches for shared structure across RNA, proteins and metabolites after balancing the layers so that one omic does not dominate simply because it contains more variables.',
          science: 'Balanced multi-block PCA on standardized variables, followed by cross-omics relationships and biological annotation.'
        },
        groups: {
          question: 'What differs between my groups?',
          short: 'Compare control/treatment, responders/non-responders or several conditions.',
          action: 'The tool first analyses each omic according to the declared design, corrects for multiple testing, then searches for coherent signals across omics and shared biological pathways.',
          science: 'Design-adjusted contrasts, BH-FDR correction, cross-omics associations and Reactome convergence.'
        },
        outcome: {
          question: 'What is associated with my clinical or experimental endpoint?',
          short: 'Link omics to a response, score, event or survival endpoint.',
          action: 'The tool chooses a model that matches the endpoint type, adjusts declared factors and, when the data support it, checks prediction on held-out subjects.',
          science: 'Endpoint-specific regression (linear, logistic, Poisson, ANCOVA or Cox) with BH-FDR and nested cross-validation when feasible.'
        },
        time: {
          question: 'What changes over time?',
          short: 'Study repeated visits or trajectories.',
          action: 'The tool accounts for repeated measurements from the same subject and tests whether trajectories differ between groups or conditions.',
          science: 'Longitudinal model with a subject random intercept and condition × time interaction, plus explicit technical-series and covariate adjustment.'
        }
      }
    };
    return copy[lang]?.[objective] || copy[lang].explore;
  }

  function endpointLabel(lang, type) {
    const labels = {
      fr: { none:'aucun critère cible', binary:'critère oui/non', multiclass:'critère à plusieurs catégories', continuous:'critère numérique continu', survival:'temps jusqu’à un événement', count:'nombre d’événements' },
      en: { none:'no target endpoint', binary:'yes/no endpoint', multiclass:'multi-category endpoint', continuous:'continuous numeric endpoint', survival:'time-to-event endpoint', count:'event-count endpoint' }
    };
    return labels[lang]?.[type] || type;
  }

  function designLabel(lang, type) {
    const labels = {
      fr: { independent:'unités indépendantes', paired:'échantillons appariés', crossover:'comparaison intra-sujet', repeated:'mesures répétées' },
      en: { independent:'independent units', paired:'paired samples', crossover:'within-subject comparison', repeated:'repeated measurements' }
    };
    return labels[lang]?.[type] || type;
  }

  function renderGuide() {
    const workflow = document.querySelector('section.simple-steps, section.workflow');
    if (!workflow) return;
    let panel = document.querySelector('[data-testid="multiomics-beginner-guide"]');
    if (!panel) {
      panel = document.createElement('section');
      panel.className = 'licence-guide';
      panel.dataset.testid = 'multiomics-beginner-guide';
      workflow.insertAdjacentElement('afterend', panel);
    }

    const lang = language();
    const objective = currentObjective();
    const outcome = currentOutcomeType();
    const design = currentDesign();
    const state = `${lang}|${objective}|${outcome}|${design}`;
    if (panel.dataset.state === state) return;

    const cards = ['explore','groups','outcome','time'].map((key) => {
      const item = objectiveCopy(lang, key);
      const selected = key === objective;
      return `<button type="button" class="licence-question${selected ? ' active' : ''}" data-licence-objective="${key}" aria-pressed="${selected}">
        <span>${selected ? (lang === 'en' ? 'Selected' : 'Sélectionné') : (lang === 'en' ? 'Choose' : 'Choisir')}</span>
        <strong>${item.question}</strong>
        <small>${item.short}</small>
      </button>`;
    }).join('');

    const item = objectiveCopy(lang, objective);
    const context = objective === 'outcome'
      ? (lang === 'en' ? `Current endpoint: <b>${endpointLabel(lang, outcome)}</b>.` : `Critère actuellement déclaré : <b>${endpointLabel(lang, outcome)}</b>.`)
      : (lang === 'en' ? `Current design: <b>${designLabel(lang, design)}</b>.` : `Plan actuellement déclaré : <b>${designLabel(lang, design)}</b>.`);

    const wasOpen = panel.querySelector('details.licence-guide-body')?.open === true;
    panel.innerHTML = `
      <details class="licence-guide-body">
        <summary>${lang === 'en' ? 'Help me choose a question' : 'M’aider à choisir ma question'}</summary>
      <div class="licence-guide-head">
        <div><span>${lang === 'en' ? 'Guided mode · undergraduate level' : 'Mode guidé · niveau licence'}</span><h2>${lang === 'en' ? 'Start with the biological question, not the method name' : 'Commencez par la question biologique, pas par le nom d’une méthode'}</h2></div>
        <button type="button" class="licence-demo" data-licence-demo>${lang === 'en' ? 'Try a complete example' : 'Voir un exemple complet'}</button>
      </div>
      <p>${lang === 'en' ? 'Choose the sentence closest to your question. The tool selects the compatible analysis family and keeps the exact statistical method in the reproducible report.' : 'Choisissez la phrase la plus proche de votre question. L’outil sélectionne la famille d’analyse compatible et conserve la méthode statistique exacte dans le rapport reproductible.'}</p>
      <div class="licence-question-grid">${cards}</div>
      <div class="licence-plan" data-testid="multiomics-guided-plan">
        <div><span>${lang === 'en' ? 'What the tool will do' : 'Ce que l’outil va faire'}</span><strong>${item.question}</strong><p>${item.action}</p><small>${context}</small></div>
        <ol>
          <li><b>1</b>${lang === 'en' ? 'Check samples, missing data and technical structure.' : 'Vérifier les échantillons, les données manquantes et la structure technique.'}</li>
          <li><b>2</b>${lang === 'en' ? 'Analyse each omic with a method compatible with its data type and study design.' : 'Analyser chaque omique avec une méthode compatible avec son type de données et le plan d’étude.'}</li>
          <li><b>3</b>${lang === 'en' ? 'Combine evidence across omics instead of simply concatenating every variable.' : 'Combiner les preuves entre omiques plutôt que simplement concaténer toutes les variables.'}</li>
          <li><b>4</b>${lang === 'en' ? 'Report interpretable signals, uncertainty and biological pathways.' : 'Présenter les signaux interprétables, leur incertitude et les voies biologiques.'}</li>
        </ol>
      </div>
      <details><summary>${lang === 'en' ? 'Scientific details of the selected route' : 'Détails scientifiques de la route choisie'}</summary><p>${item.science}</p><p>${lang === 'en' ? 'The interface does not replace experimental-design judgement: biological independence, confounding, missing layers and sample size still determine what can be interpreted.' : 'L’interface ne remplace pas le jugement sur le plan expérimental : indépendance biologique, confusion, couches manquantes et effectif déterminent toujours ce qui peut être interprété.'}</p></details>
      </details>
    `;
    if (wasOpen) panel.querySelector('details.licence-guide-body').open = true;
    panel.dataset.state = state;
  }

  function conclusionCopy(lang, objective) {
    const copy = {
      fr: {
        explore: ['Vous pouvez identifier des structures communes et des variables qui évoluent ensemble.', 'Vous ne pouvez pas conclure qu’une variable en cause une autre.'],
        groups: ['Vous pouvez décrire des différences associées aux groupes après les ajustements déclarés.', 'Une différence entre groupes ne prouve pas à elle seule un mécanisme causal.'],
        outcome: ['Vous pouvez identifier des associations avec le critère et, si la validation croisée est disponible, évaluer une performance prédictive interne.', 'Une bonne prédiction ne démontre ni causalité ni validité dans une nouvelle population.'],
        time: ['Vous pouvez décrire quelles variables évoluent avec le temps et si les trajectoires diffèrent entre conditions.', 'Une association temporelle ne démontre pas que la variation moléculaire cause l’évolution clinique.']
      },
      en: {
        explore: ['You can identify shared structure and variables that move together.', 'You cannot conclude that one variable causes another.'],
        groups: ['You can describe group-associated differences after the declared adjustments.', 'A group difference alone does not prove a causal mechanism.'],
        outcome: ['You can identify endpoint associations and, when cross-validation is available, assess internal predictive performance.', 'Good prediction does not establish causality or validity in a new population.'],
        time: ['You can describe variables that change over time and whether trajectories differ between conditions.', 'A temporal association does not prove that a molecular change causes the clinical change.']
      }
    };
    return copy[lang]?.[objective] || copy[lang].explore;
  }

  function renderResultGuide() {
    const readiness = document.querySelector('[data-testid="multiomics-plain-readiness"]');
    if (!readiness) return;
    let panel = document.querySelector('[data-testid="multiomics-result-reading-guide"]');
    if (!panel) {
      panel = document.createElement('section');
      panel.className = 'licence-result-guide';
      panel.dataset.testid = 'multiomics-result-reading-guide';
      readiness.insertAdjacentElement('afterend', panel);
    }
    const lang = language();
    const objective = currentObjective();
    const state = `${lang}|${objective}`;
    if (panel.dataset.state === state) return;
    const [canSay, cannotSay] = conclusionCopy(lang, objective);
    panel.innerHTML = `
      <details class="licence-result-body">
        <summary>${lang === 'en' ? 'How to read this analysis' : 'Comment lire cette analyse'}</summary>
      <span class="licence-kicker">${lang === 'en' ? 'Recommended reading order' : 'Ordre de lecture conseillé'}</span>
      <h3>${lang === 'en' ? 'Understand the result in four passes' : 'Comprendre le résultat en quatre passages'}</h3>
      <div class="licence-result-grid">
        <article><b>1</b><strong>${lang === 'en' ? 'Quality first' : 'Qualité d’abord'}</strong><span>${lang === 'en' ? 'Check missing data, technical series and suspicious measurements before biological interpretation.' : 'Vérifiez les données manquantes, les séries techniques et les mesures suspectes avant toute interprétation biologique.'}</span></article>
        <article><b>2</b><strong>${lang === 'en' ? 'Within each omic' : 'Dans chaque omique'}</strong><span>${lang === 'en' ? 'Read direction/magnitude and adjusted q/FDR together. A small p-value alone is not enough.' : 'Lisez ensemble la direction/l’amplitude et le q/FDR corrigé. Une petite p-value seule ne suffit pas.'}</span></article>
        <article><b>3</b><strong>${lang === 'en' ? 'Across omics' : 'Entre les omiques'}</strong><span>${lang === 'en' ? 'Give more weight to signals that converge across RNA, proteins and metabolites than to an isolated feature.' : 'Donnez plus de poids aux signaux qui convergent entre RNA, protéines et métabolites qu’à une variable isolée.'}</span></article>
        <article><b>4</b><strong>${lang === 'en' ? 'Biological interpretation' : 'Interprétation biologique'}</strong><span>${lang === 'en' ? 'Pathways organize signals; enrichment is not proof that a pathway is activated or causal.' : 'Les voies biologiques organisent les signaux ; un enrichissement ne prouve pas qu’une voie est activée ni causale.'}</span></article>
      </div>
      <div class="licence-conclusion-grid"><article><strong>${lang === 'en' ? 'What you may conclude' : 'Ce que vous pouvez conclure'}</strong><span>${canSay}</span></article><article><strong>${lang === 'en' ? 'What you should not conclude' : 'Ce que vous ne devez pas conclure'}</strong><span>${cannotSay}</span></article></div>
      </details>
    `;
    panel.dataset.state = state;
  }

  function installStyle() {
    if (document.getElementById('multiomics-licence-guide-style')) return;
    const style = document.createElement('style');
    style.id = 'multiomics-licence-guide-style';
    style.textContent = `
      .licence-guide>details>summary,.licence-result-guide>details>summary{cursor:pointer;font-weight:700}
      .licence-guide{max-width:1180px;margin:18px auto 24px;padding:20px;border:1px solid var(--border-strong,#c9c9c9);border-radius:14px;background:var(--bg-secondary,#f7f7f7)}
      .licence-guide-head{display:flex;justify-content:space-between;align-items:flex-start;gap:18px}.licence-guide-head span,.licence-plan>div>span,.licence-kicker{display:block;font-size:.75rem;text-transform:uppercase;letter-spacing:.08em;color:var(--text-secondary,#666)}.licence-guide-head h2{margin:4px 0 0;font-size:1.28rem}.licence-guide>p{max-width:88ch;color:var(--text-secondary,#555);line-height:1.5}.licence-demo{border:1px solid var(--border-strong,#bbb);background:var(--bg-primary,#fff);color:inherit;border-radius:9px;padding:9px 12px;font-weight:700;cursor:pointer}
      .licence-question-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:10px;margin:14px 0}.licence-question{text-align:left;border:1px solid var(--border-subtle,#d4d4d4);background:var(--bg-primary,#fff);color:inherit;padding:13px;border-radius:10px;cursor:pointer}.licence-question:hover,.licence-question:focus-visible{border-color:currentColor}.licence-question.active{outline:2px solid currentColor;outline-offset:1px}.licence-question span{font-size:.72rem;text-transform:uppercase;letter-spacing:.06em;opacity:.7}.licence-question strong,.licence-question small{display:block}.licence-question strong{margin-top:5px}.licence-question small{margin-top:5px;line-height:1.4;color:var(--text-secondary,#555)}
      .licence-plan{display:grid;grid-template-columns:minmax(0,1fr) minmax(0,1.1fr);gap:18px;margin-top:16px;padding:16px;border:1px solid var(--border-subtle,#ddd);border-radius:11px;background:var(--bg-primary,#fff)}.licence-plan strong{display:block;margin-top:5px}.licence-plan p{margin:6px 0;line-height:1.5;color:var(--text-secondary,#555)}.licence-plan ol{list-style:none;margin:0;padding:0;display:grid;gap:7px}.licence-plan li{display:grid;grid-template-columns:26px 1fr;gap:8px;align-items:start;font-size:.9rem;line-height:1.4}.licence-plan li b,.licence-result-grid article>b{display:grid;place-items:center;width:24px;height:24px;border:1px solid var(--border-subtle,#ccc);border-radius:999px}.licence-guide details{margin-top:12px}.licence-guide summary{cursor:pointer;font-weight:700}.licence-guide details p{max-width:90ch;color:var(--text-secondary,#555);line-height:1.5}
      .licence-result-guide{margin:18px 0;padding:18px;border:1px solid var(--border-strong,#c9c9c9);border-radius:12px;background:var(--bg-primary,#fff)}.licence-result-guide h3{margin:4px 0 12px;font-size:1.08rem}.licence-result-grid{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:9px}.licence-result-grid article,.licence-conclusion-grid article{padding:11px;border:1px solid var(--border-subtle,#ddd);border-radius:9px}.licence-result-grid article>b{margin-bottom:8px}.licence-result-grid strong,.licence-result-grid span,.licence-conclusion-grid strong,.licence-conclusion-grid span{display:block}.licence-result-grid span,.licence-conclusion-grid span{margin-top:5px;font-size:.87rem;line-height:1.42;color:var(--text-secondary,#555)}.licence-conclusion-grid{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:9px;margin-top:10px}
      @media(max-width:900px){.licence-plan{grid-template-columns:1fr}.licence-result-grid{grid-template-columns:repeat(2,minmax(0,1fr))}}
      @media(max-width:760px){.licence-guide{margin-left:12px;margin-right:12px;padding:15px}.licence-guide-head{flex-direction:column}.licence-question-grid,.licence-result-grid,.licence-conclusion-grid{grid-template-columns:1fr}}
    `;
    document.head.appendChild(style);
  }

  function installInteractions() {
    if (document.documentElement.dataset.multiomicsLicenceGuide === '1') return;
    document.documentElement.dataset.multiomicsLicenceGuide = '1';
    document.addEventListener('click', (event) => {
      const target = event.target instanceof Element ? event.target : null;
      if (!target) return;
      const question = target.closest('[data-licence-objective]');
      if (question) {
        const value = question.getAttribute('data-licence-objective');
        const select = findSelectWithOption('explore');
        if (select && value && [...select.options].some((option) => option.value === value)) {
          select.value = value;
          select.dispatchEvent(new Event('change', { bubbles: true }));
          renderGuide();
          renderResultGuide();
        }
        return;
      }
      if (target.closest('[data-licence-demo]')) {
        const demo = document.querySelector('[data-testid="multiomics-load-demo"]');
        if (demo instanceof HTMLElement) demo.click();
      }
    });
    document.addEventListener('change', () => { renderGuide(); renderResultGuide(); });
  }

  let queued = false;
  function refresh() {
    if (queued) return;
    queued = true;
    requestAnimationFrame(() => {
      queued = false;
      installStyle();
      installInteractions();
      renderGuide();
      renderResultGuide();
    });
  }

  const observer = new MutationObserver(refresh);
  observer.observe(document.documentElement, { childList:true, subtree:true, attributes:true, attributeFilter:['lang'] });
  document.addEventListener('DOMContentLoaded', refresh, { once:true });
  refresh();
})();
