// @ts-nocheck
// Browser-side gate for the actual /health contract of the local R engine.
// Connectivity alone (HTTP 200, unrelated JSON, or a package list) MUST NOT
// be displayed as reference-computation readiness or scientific validation.

const backendNames = new Set([
  'PMx Explain reference R backend', // historically emitted by server.R
  'MultiOmics reference R backend'
]);
const fail = (code, message, messageEn) => ({ok:false,code,message,messageEn});
const object = (x) => x !== null && typeof x === 'object' &&
  !Array.isArray(x);
const installed = (packages, name) => packages[name] === true;

/**
 * Verifies that the response belongs to the expected R adapter, then
 * inspects minimum available methods for the declared *analysis context*.
 * This is NOT a package version lock, statistical study review, secure
 * authentication of the service, or a claim that a fitted method succeeded.
 *
 * @param {unknown} health JSON returned by GET /health
 * @param {{intent?:string,objective?:string,outcomeType?:string,designType?:string,
 *          hasTranscriptomics?:boolean,transcriptomicsValues?:string,
 *          hasProteomics?:boolean,hasMetabolomics?:boolean}} context
 */
export function evaluateReferenceBackendHealth(health, context = {}) {
  if (!object(health) || health.status !== 'ok' ||
      !backendNames.has(health.engine) ||
      !/^\d+\.\d+\.\d+(?:[-+][A-Za-z0-9.-]+)?$/.test(String(health.version||'')) ||
      !object(health.capabilities) || health.capabilities.reference_analysis !== true ||
      !object(health.packages) || !object(health.runtime) ||
      typeof health.runtime.RVersion !== 'string' ||
      !/^\d+\.\d+(?:\.\d+)?$/.test(health.runtime.RVersion)) {
    return fail(
      'invalid_health_contract',
      'Le serveur répond, mais ne présente pas le contrat attendu du moteur MultiOmics R. Vérifiez sa version avant d’envoyer des données.',
      'A server responded, but it does not expose the expected MultiOmics R contract. Verify the backend version before uploading data.'
    );
  }
  const packages = health.packages;
  const scientific = ['DESeq2','edgeR','limma','lmerTest','survival','MOFA2','mixOmics','fgsea'];
  if (scientific.some(name=>typeof packages[name] !== 'boolean')) {
    return fail(
      'incomplete_package_manifest',
      'La liste des packages R est incomplète : impossible de confirmer la disponibilité des méthodes.',
      'The R package inventory is incomplete: method availability cannot be assessed.'
    );
  }
  if (!scientific.some(name=>installed(packages,name))) {
    return fail('no_scientific_methods',
      'Aucune méthode statistique de référence n’est installée sur le backend R.',
      'No reference statistical method is installed on the R backend.');
  }
  if (context.intent !== 'confirmatory') {
    return {ok:true,code:'verified_exploratory',
      message:'Moteur MultiOmics R reconnu ; les méthodes réellement exécutées seront contrôlées dans les résultats.',
      messageEn:'MultiOmics R backend recognized; completed methods must still be checked in the results.'};
  }

  const objective = context.objective;
  if (objective === 'groups') {
    if (context.hasTranscriptomics && context.transcriptomicsValues === 'raw_counts') {
      if (!installed(packages,'DESeq2') &&
          !(installed(packages,'edgeR') && installed(packages,'limma'))) {
        return fail('missing_rnaseq_methods',
          'RNA-seq de comptages : installer DESeq2 ou edgeR + limma ; aucun modèle de référence applicable n’est disponible.',
          'Raw-count RNA-seq requires DESeq2 or edgeR + limma; no applicable reference model is available.');
      }
    }
    // All loaded layers must have a model: DESeq2 alone is not sufficient
    // when an additional processed proteomics/metabolomics layer is loaded.
    const needsLimma=Boolean(context.hasProteomics || context.hasMetabolomics ||
      (context.hasTranscriptomics && context.transcriptomicsValues !== 'raw_counts'));
    if (needsLimma && !installed(packages,'limma')) {
      return fail('missing_limma',
        'Comparaison de groupes : le package R limma requis pour les couches non-RNA-seq en comptages n’est pas disponible.',
        'Group comparisons require R limma for non-raw-count layers, but it is unavailable.');
    }
  } else if (objective === 'time') {
    if (!installed(packages,'lmerTest')) {
      return fail('missing_lmerTest',
        'Analyse longitudinale : le package lmerTest requis n’est pas disponible.',
        'Longitudinal analysis requires the R package lmerTest, which is unavailable.');
    }
  } else if (objective === 'outcome') {
    if (context.outcomeType !== 'survival') {
      return fail('unsupported_confirmatory_outcome',
        'Ce type d’outcome n’a pas encore de tableau R natif projeté et qualifié dans le mode confirmatoire. Utilisez le mode exploratoire jusqu’à l’ajout d’un contrat R spécifique.',
        'This outcome type does not yet have a qualified native R table projection for confirmatory mode. Use exploratory mode until a dedicated R contract is implemented.');
    }
    if (!installed(packages,'survival')) {
      return fail('missing_survival',
        'Survie : le package R survival requis pour le modèle de Cox de référence n’est pas disponible.',
        'Survival analysis requires the R survival package for the reference Cox model, but it is unavailable.');
    }
  } else {
    return fail('no_confirmatory_model',
      'L’exploration ne constitue pas un modèle confirmatoire : définissez d’abord une question statistique.',
      'Exploration is not a confirmatory model: first choose a statistical research question.');
  }
  return {ok:true,code:'verified_method_available',
    message:'Connexion R vérifiée et méthode potentiellement disponible ; le design, les hypothèses et l’exécution restent à contrôler.',
    messageEn:'R connection verified and a method may be available; study design, assumptions, and successful execution still need review.'};
}
