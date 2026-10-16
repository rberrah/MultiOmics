// Compatibility gate for local R backends. Older MOFA2 adapters accepted
// sample × feature matrices while create_mofa(list) needs feature × sample.
// Reject those results instead of displaying biologically incorrect factors.
export const MOFA2_ORIENTATION_CONTRACT = 'mofa2-feature-rows-sample-columns-v2';

/** @param {{objective:string,omicCount:number}} context */
export function needsMofa2OrientationContract({objective,omicCount}) {
  return objective === 'explore' && Number(omicCount) >= 2;
}

/** @param {any} health @param {{objective:string,omicCount:number}} context */
export function mofa2BackendCompatibility(health,context) {
  if(!needsMofa2OrientationContract(context))
    return {compatible:true,reason:'not_applicable'};
  if(health?.capabilities?.mofa2_orientation_contract !== MOFA2_ORIENTATION_CONTRACT) {
    return {
      compatible:false,reason:'mofa2_backend_outdated',
      message:'Ce moteur R utilise une ancienne interface MOFA2. Mettez à jour le backend R avant toute intégration multi-omique. Les figures du navigateur restent disponibles.',
      messageEn:'This R backend uses an outdated MOFA2 adapter. Update the local R backend before multi-omics integration. Browser figures remain available.'
    };
  }
  return {compatible:true,reason:'verified_health_contract'};
}

/** @param {any} body */
export function validateMofa2MethodResult(body) {
  if(!body || typeof body!=='object')return body;
  const methods=body.methods;
  if(!methods || typeof methods!=='object')return body;
  const item=methods.mofa2;
  if(item?.status!=='ok')return body;
  const safe=item.summary?.implementationContract===MOFA2_ORIENTATION_CONTRACT &&
    item.summary?.matrixOrientation==='feature_rows_subject_columns_for_MOFA2' &&
    item.summary?.subjectIdentityPreserved===true;
  if(safe)return body;
  return {
    ...body,
    methods:{
      ...methods,
      mofa2:{
        method:'MOFA2',
        status:'blocked',
        message:'MOFA2 non affiché : ce résultat R ne prouve pas que les variables et les sujets ont été correctement orientés. Mettez à jour le backend et relancez.',
        messageEn:'MOFA2 result hidden: the backend did not prove correct feature and subject orientation. Update the R backend and rerun.',
        compatibilityStatus:'orientation_contract_missing_or_invalid'
      }
    },
    mofa2OrientationBlocked:true
  };
}
