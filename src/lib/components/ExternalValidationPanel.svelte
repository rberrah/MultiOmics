<script>
  /** @type {string} */
  export let backendUrl = 'http://127.0.0.1:8787';
  /** @type {string} */
  export let language = 'fr';

  /** @type {File | null} */
  let predictionsFile = null;
  /** @type {File | null} */
  let configFile = null;
  let status = 'idle';
  let error = '';
  /** @type {any} */
  let result = null;
  /** @type {any} */
  let parsedConfig = null;

  /** @param {string} fr @param {string} en */
  const t = (fr, en) => language === 'en' ? en : fr;

  /** @param {Event} event */
  function pickPredictions(event) {
    const input = event.currentTarget;
    predictionsFile = input instanceof HTMLInputElement ? input.files?.[0] || null : null;
    result = null;
    error = '';
  }

  /** @param {Event} event */
  async function pickConfig(event) {
    const input = event.currentTarget;
    configFile = input instanceof HTMLInputElement ? input.files?.[0] || null : null;
    result = null;
    error = '';
    parsedConfig = null;
    if (!configFile) return;
    try {
      parsedConfig = JSON.parse(await configFile.text());
    } catch {
      error = t('Le fichier de configuration JSON est invalide.', 'The JSON configuration file is invalid.');
    }
  }

  function baseUrl() {
    return String(backendUrl || '').trim().replace(/\/$/, '');
  }

  async function runValidation() {
    error = '';
    result = null;
    if (!predictionsFile || !parsedConfig) {
      error = t('Chargez le CSV de prédictions figées et sa configuration JSON.', 'Load the frozen-prediction CSV and its JSON configuration.');
      return;
    }
    if (parsedConfig.independent_cohort !== true) {
      error = t(
        'La configuration ne déclare pas cette cohorte comme indépendante. Corrigez la provenance avant de demander une validation externe.',
        'The configuration does not declare this cohort independent. Correct the provenance before requesting external validation.'
      );
      return;
    }
    const url = baseUrl();
    if (!url) {
      error = t('URL du backend R manquante.', 'R backend URL is missing.');
      return;
    }

    status = 'running';
    try {
      const response = await fetch(url + '/external-validation', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
        body: JSON.stringify({
          predictionsCsv: await predictionsFile.text(),
          config: parsedConfig
        })
      });
      const body = await response.json().catch(() => ({}));
      if (!response.ok || body.status === 'error') {
        throw new Error(body.message || ('HTTP ' + response.status));
      }
      result = body;
      status = 'done';
    } catch (cause) {
      status = 'error';
      error = cause instanceof Error ? cause.message : t('Validation externe impossible.', 'External validation failed.');
    }
  }

  function downloadResult() {
    if (!result) return;
    const blob = new Blob([JSON.stringify(result, null, 2)], { type: 'application/json' });
    const href = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = href;
    anchor.download = 'multiomics_external_validation.json';
    anchor.click();
    URL.revokeObjectURL(href);
  }
</script>

<div class="external-validation" data-testid="external-validation-panel">
  <div>
    <strong>{t('Validation externe depuis le navigateur', 'External validation from the browser')}</strong>
    <p>{t(
      'Chargez uniquement des prédictions produites par un modèle déjà figé. Cette étape ne réentraîne rien et appelle le backend R local.',
      'Upload only predictions produced by an already frozen model. This step does not refit anything and calls the local R backend.'
    )}</p>
  </div>

  <div class="files">
    <label>
      <span>{t('Prédictions figées (.csv)', 'Frozen predictions (.csv)')}</span>
      <input data-testid="external-validation-predictions" type="file" accept=".csv,text/csv" onchange={pickPredictions} />
    </label>
    <label>
      <span>{t('Configuration (.json)', 'Configuration (.json)')}</span>
      <input data-testid="external-validation-config" type="file" accept=".json,application/json" onchange={pickConfig} />
    </label>
  </div>

  {#if parsedConfig}
    <div class="summary" data-testid="external-validation-config-summary">
      <span><b>{t('Critère', 'Outcome')}</b> {parsedConfig.outcome_type || '—'}</span>
      <span><b>{t('Cohorte', 'Cohort')}</b> {parsedConfig.cohort_label || '—'}</span>
      <span><b>{t('Indépendante', 'Independent')}</b> {parsedConfig.independent_cohort === true ? t('oui', 'yes') : t('non', 'no')}</span>
      <span><b>Bootstrap</b> {parsedConfig.bootstrap_repetitions ?? '—'}</span>
    </div>
  {/if}

  <div class="actions">
    <button
      class="btn btn-primary"
      type="button"
      data-testid="external-validation-run"
      disabled={!predictionsFile || !parsedConfig || status === 'running'}
      onclick={runValidation}
    >
      {status === 'running' ? t('Validation…', 'Validating…') : t('Évaluer la cohorte indépendante', 'Evaluate independent cohort')}
    </button>
    {#if result}
      <button class="btn btn-outline" type="button" onclick={downloadResult}>{t('Télécharger le résultat', 'Download result')}</button>
    {/if}
  </div>

  {#if error}<p class="error" data-testid="external-validation-error">{error}</p>{/if}

  {#if result}
    <div class="result" data-testid="external-validation-result">
      <strong>{result.validation_status || result.status || t('Validation calculée', 'Validation computed')}</strong>
      <p>{result.cohort_label || parsedConfig?.cohort_label || ''}</p>
      {#if result.cohort_coverage?.complete === false || (result.validation_status || result.status) === 'external_validation_incomplete'}
        <p class="error" role="alert" data-testid="external-validation-incomplete">
          {t('Validation externe incomplète : certaines prédictions ou observations n’ont pas été évaluées. Le score ci-dessous porte seulement sur les données disponibles et ne permet pas de conclure sur toute la cohorte.',
            'Incomplete external validation: some predictions or outcomes were excluded. The score below covers evaluable records only and cannot describe the full cohort.')}
          {#if result.cohort_coverage}
            {result.cohort_coverage.evaluated_rows} / {result.cohort_coverage.submitted_rows} {t('lignes évaluées.', 'rows evaluated.')}
          {/if}
        </p>
      {/if}
      {#if result.subject_identity_audit?.status === 'not_verifiable_without_subject_id'}
        <p role="status" data-testid="external-validation-identity-warning">
          {t('Identifiants des sujets absents : l’absence de doublons ne peut pas être vérifiée. Ajoutez une colonne subject_id et documentez l’indépendance de la cohorte.',
            'No subject identifiers: repeated patients cannot be ruled out. Add a subject_id column and document cohort independence.')}
        </p>
      {/if}
      <pre>{JSON.stringify(result.metrics || result, null, 2)}</pre>
    </div>
  {/if}
</div>

<style>
  .external-validation { display:grid; gap:12px; margin-top:14px; padding:14px; border:1px solid var(--border-subtle); border-radius:var(--radius); background:var(--bg-secondary); }
  .external-validation p { margin:5px 0 0; color:var(--text-secondary); font-size:var(--text-sm); }
  .files { display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:10px; }
  .files label { display:grid; gap:5px; font-size:var(--text-sm); }
  .summary { display:flex; flex-wrap:wrap; gap:8px; }
  .summary span { padding:5px 7px; border:1px solid var(--border-subtle); border-radius:6px; font-size:var(--text-xs); }
  .actions { display:flex; flex-wrap:wrap; gap:8px; }
  .error { color:var(--danger,#b42318); }
  .result pre { max-height:320px; overflow:auto; padding:10px; background:var(--bg-primary); border:1px solid var(--border-subtle); font-size:var(--text-xs); }
  @media (max-width:700px) { .files { grid-template-columns:1fr; } }
</style>
