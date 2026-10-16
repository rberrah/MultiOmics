(() => {
  if (!window.location.pathname.includes('/multiomics/tool')) return;

  const language = () => document.documentElement.lang === 'en' ? 'en' : 'fr';

  // Values originating in uploaded metadata must never be re-interpreted as markup.
  const escapeHtml = (value) => String(value).replace(/[&<>"']/g, (char) => ({
    '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;'
  })[char]);

  function renderRunReadiness() {
    const runButton = document.querySelector('[data-testid="multiomics-run"]');
    if (!(runButton instanceof HTMLButtonElement)) return;
    const runBox = runButton.closest('.run-box');
    if (!runBox) return;
    const section = runBox.closest('section');
    const status = section?.querySelector('.status');
    if (!status) return;

    let panel = document.querySelector('[data-testid="multiomics-run-readiness"]');
    if (!panel) {
      panel = document.createElement('section');
      panel.className = 'licence-run-readiness';
      panel.dataset.testid = 'multiomics-run-readiness';
      runBox.insertAdjacentElement('beforebegin', panel);
    }

    const lang = language();
    const running = /analyse en cours|analyzing/i.test(runButton.textContent || '');
    const ready = !runButton.disabled && !running;
    const intent = document.querySelector('[data-testid="multiomics-analysis-intent"] select')?.value === 'confirmatory'
      ? 'confirmatory' : 'exploratory';
    const statusTitle = status.querySelector('strong')?.textContent?.trim() || '';
    const statusMessage = status.querySelector('span')?.textContent?.trim() || '';
    const state = `${lang}|${ready}|${running}|${intent}|${statusTitle}|${statusMessage}`;
    if (panel.dataset.state === state) return;

    const wasOpen = panel.querySelector('details.readiness-extra')?.open === true;
    panel.innerHTML = `
      <div class="licence-run-head">
        <div><span>${lang === 'en' ? 'Before running' : 'Avant de lancer'}</span><h3>${lang === 'en' ? 'Is the analysis ready?' : 'L’analyse est-elle prête ?'}</h3></div>
        <b class="${ready ? 'ready' : 'todo'}">${running ? (lang === 'en' ? 'Analysis running' : 'Analyse en cours') : ready ? (intent === 'confirmatory' ? (lang === 'en' ? 'Ready to run R checks' : 'Prêt pour les contrôles R') : (lang === 'en' ? 'Ready for exploration' : 'Prêt pour explorer')) : (lang === 'en' ? 'To complete' : 'À compléter')}</b>
      </div>
      <p><strong>${escapeHtml(statusTitle)}</strong> — ${escapeHtml(statusMessage)}</p>
      <details class="readiness-extra"><summary>${lang === 'en' ? 'See all checks' : 'Voir les conditions à vérifier'}</summary>
      <div class="licence-run-rules">
        <span>${lang === 'en' ? 'The tool unlocks the analysis only when:' : 'L’outil débloque l’analyse seulement lorsque :'}</span>
        <ul>
          <li>${lang === 'en' ? 'at least one omics matrix is loaded (one layer for single-omics, two or three for multi-omics);' : 'au moins une matrice omique est chargée (une seule pour l’analyse mono-omique, deux ou trois pour la multi-omique) ;'}</li>
          <li>${lang === 'en' ? 'subject, specimen and assay identifiers connect the sample sheet to the matrices;' : 'les identifiants sujet, prélèvement et mesure relient correctement le tableau des échantillons aux matrices ;'}</li>
          <li>${lang === 'en' ? 'the variables required by the chosen biological question are declared;' : 'les variables nécessaires à la question biologique choisie sont déclarées ;'}</li>
          <li>${lang === 'en' ? 'the declared study design is supported and has no deterministic blocking confounding.' : 'le plan d’étude déclaré est pris en charge et ne présente pas de confusion déterministe bloquante.'}</li>
        </ul>
      </div>
      </details>
      <small>${running
        ? (lang === 'en' ? 'The analysis is running. Wait for the results and review all warnings before interpreting them.' : 'L’analyse est en cours. Examinez les résultats et toutes les alertes avant de les interpréter.')
        : ready
        ? (lang === 'en' ? 'The input checks allow the calculation to start, not a scientific conclusion. Quality warnings and the study design still need review.' : 'Les contrôles d’entrée permettent de lancer le calcul, pas de valider un résultat scientifique. Le plan d’étude et les alertes de qualité restent à examiner.')
        : (lang === 'en' ? 'The message above shows what to resolve before starting. No reference method or scientific result has been validated yet.' : 'Le message ci-dessus indique ce qui bloque le lancement. Aucune méthode ni conclusion scientifique n’est encore validée.')}</small>
    `;
    if (wasOpen) panel.querySelector('details.readiness-extra').open = true;
    panel.dataset.state = state;
  }

  function installStyle() {
    if (document.getElementById('multiomics-run-readiness-style')) return;
    const style = document.createElement('style');
    style.id = 'multiomics-run-readiness-style';
    style.textContent = `
      .readiness-extra>summary{cursor:pointer;font-weight:700}.readiness-extra{margin-top:10px}
      .licence-run-readiness{margin:14px 0 18px;padding:15px;border:1px solid var(--border-subtle,#d4d4d4);border-radius:11px;background:var(--bg-secondary,#f7f7f7)}
      .licence-run-head{display:flex;justify-content:space-between;align-items:flex-start;gap:14px}.licence-run-head span{display:block;font-size:.74rem;text-transform:uppercase;letter-spacing:.08em;color:var(--text-secondary,#666)}.licence-run-head h3{margin:4px 0 0;font-size:1.02rem}.licence-run-head>b{border:1px solid var(--border-subtle,#ccc);border-radius:999px;padding:5px 9px;font-size:.8rem;white-space:nowrap}.licence-run-head>b.ready{font-weight:800}.licence-run-readiness>p{margin:10px 0;line-height:1.45}
      .licence-run-rules{display:grid;grid-template-columns:minmax(160px,.5fr) 1.5fr;gap:12px;align-items:start}.licence-run-rules>span{font-weight:700;font-size:.9rem}.licence-run-rules ul{margin:0;padding-left:18px}.licence-run-rules li{margin:0 0 5px;font-size:.88rem;line-height:1.4}.licence-run-readiness>small{display:block;margin-top:8px;color:var(--text-secondary,#555);line-height:1.4}
      @media(max-width:760px){.licence-run-head{flex-direction:column}.licence-run-rules{grid-template-columns:1fr}.licence-run-head>b{white-space:normal}}
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
      renderRunReadiness();
    });
  }

  const observer = new MutationObserver(refresh);
  observer.observe(document.documentElement, { childList:true, subtree:true, attributes:true, attributeFilter:['lang','disabled'] });
  document.addEventListener('change', refresh);
  document.addEventListener('DOMContentLoaded', refresh, { once:true });
  refresh();
})();
