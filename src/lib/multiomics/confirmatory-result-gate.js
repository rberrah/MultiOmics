// @ts-nocheck
// Publication-facing guard. A successful HTTP response, a valid /health
// contract, or one R fit must NOT label a multi-layer browser result as
// confirmatory when some layers still contain exploratory JS values.
const KNOWN = ['transcriptomics','proteomics','metabolomics'];
const failure = (code, layers, message, messageEn) =>
  ({ok:false,code,layers,message,messageEn});

/**
 * @param {any} result Result after projectNativeRNumbers(), not the raw /run body
 * @returns {{ok:boolean,code:string,layers:string[],message:string,messageEn:string}}
 */
export function assessConfirmatoryResult(result) {
  const purpose = result?.protocol?.analysisIntent;
  if (purpose !== 'confirmatory') {
    return {ok:true,code:'exploratory_not_certified',layers:[],
      message:'Exploration ; aucune certification scientifique.',
      messageEn:'Exploratory analysis; no scientific certification.'};
  }
  const objective = result?.protocol?.objective;
  // Current page projects *group differential* complete R tables to the
  // displayed statistics. Time and outcome still show JS approximations even
  // when an R fit is stored separately. Reject their confirmatory presentation
  // until the canonical native projection for each is implemented and tested.
  if (objective !== 'groups') {
    return failure('native_projection_not_supported',[],
      'Le mode confirmatoire de cette question n’est pas encore pris en charge dans les tableaux affichés. Choisissez exploration ou une comparaison de groupes compatible avec les tableaux R complets.',
      'Confirmatory results for this objective are not yet supported in the displayed tables. Use exploration or a group comparison backed by complete native R tables.');
  }
  const expected=KNOWN.filter(layer=>Boolean(result?.layers?.[layer]));
  if (!expected.length) return failure('no_omic_layers',[],
    'Aucune couche omique n’a été analysée.',
    'No omics layer was analysed.');
  if (result.referenceBackend?.status !== 'ok' ||
      result.referenceBackend?.referenceMethodExecuted !== true) {
    return failure('reference_r_not_executed',expected,
      'Aucune méthode R de référence applicable n’a été exécutée avec succès.',
      'No applicable native R reference method completed successfully.');
  }
  const incomplete=expected.filter(layer=>{
    const item=result.layers[layer];
    const key=item?.primaryReferenceMethodKey;
    const method=typeof key==='string' ? result.referenceBackend.methods?.[key] : null;
    const count=Number(item?.primaryReferenceTableRows);
    return item?.numericalProvenance !== 'exact_native_r_reference_table' ||
      item?.mode !== 'reference-r' ||
      result.numericSourceByLayer?.[layer] !== 'native_r_complete_table' ||
      !key?.startsWith(layer+'_differential') ||
      method?.status !== 'ok' || method?.referenceTableComplete !== true ||
      typeof method?.referenceTableCsv !== 'string' ||
      method.referenceTableCsv.trim().length===0 ||
      !Number.isSafeInteger(count) || count<=0 ||
      !Number.isSafeInteger(Number(method.referenceTableRows)) ||
      Number(method.referenceTableRows)!==count ||
      !Array.isArray(item.rows) || item.rows.length!==count;
  });
  if (incomplete.length || result.displayedReferenceMethodCount !== expected.length) {
    const blocked=incomplete.length ? incomplete : expected;
    return failure('incomplete_native_layer_projection',blocked,
      'Analyse confirmatoire bloquée : toutes les omiques chargées doivent afficher leurs tableaux R complets. Couche(s) sans preuve R exhaustive : '+blocked.join(', ')+'. Aucune valeur exploratoire du navigateur ne remplace un calcul R.',
      'Confirmatory analysis blocked: every loaded layer needs its complete native R table. Missing reference evidence for: '+blocked.join(', ')+'. Browser approximations cannot replace R results.');
  }
  return {ok:true,code:'native_tables_complete_not_certified',layers:expected,
    message:'Tableaux R complets affichés pour chaque couche ; validité du plan, hypothèses statistiques et résultats biologiques à vérifier indépendamment.',
    messageEn:'Complete native R tables displayed for every layer; study design, statistical assumptions and biological findings still require independent review.'};
}
