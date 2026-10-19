// Methods-only, privacy-minimized reproducibility manifest: no patient-level
// rows, raw expression values, predictions, or exported PHI.
/** @param {unknown} x */
const asNumber=(x)=>typeof x==='number'&&Number.isFinite(x)?x:null;
/** @param {unknown} x */
const safe=(x)=>typeof x==='string'?x.slice(0,160):null;
/** @param {any} result */
export function buildResearchMethodsManifest(result){
  if(!result||typeof result!=='object'||!result.protocol||
    !result.engine||!result.layers || typeof result.layers!=='object')
    throw new Error('Methods export requires a completed multi-omics analysis.');
  const p=result.protocol;
  /** @type {Record<string, any>} */
  const layers={};
  for(const [layer,entry] of Object.entries(result.layers)){
    if(!['transcriptomics','proteomics','metabolomics'].includes(layer))
      continue;
    const qc=entry?.qc||{};
    // Median sample missingness can be zero even when one assay is missing
    // biologically relevant signals. Preserve the observed maximum and
    // the number of affected assays without exposing assay/subject IDs.
    const profiles=Array.isArray(qc?.sampleMetrics)?qc.sampleMetrics:[];
    /** @type {number[]} */
    const missingProfiles=profiles.flatMap((/** @type {any} */ item)=>{
      const fraction=asNumber(item?.missingFraction);
      return fraction!==null && fraction>=0 && fraction<=1 ? [fraction] : [];
    });
    const anyMissing=missingProfiles.some((fraction)=>fraction>0);
    const incompleteProfiles=missingProfiles.filter((fraction)=>fraction>0);
    const maxSampleMissingFraction=missingProfiles.length
      ? Math.max(...missingProfiles):null;

    const backendKey=p.objective==='time'
      ? layer+'_longitudinal'
      : p.objective==='outcome' && p.outcomeType==='survival'
        ? layer+'_survival'
        : layer==='transcriptomics' &&
          qc?.preprocessingAudit?.declaredValueType==='raw_counts'
          ? layer+'_differential_deseq2'
          : layer+'_differential';
    const backend=result.referenceBackend?.methods?.[backendKey];
    layers[layer]={
      reported_model:safe(entry?.inferenceMethod||entry?.mode)||'not_reported',
      declared_value_type:safe(qc?.preprocessingAudit?.declaredValueType)||'not_reported',
      features_input:asNumber(qc?.featuresBefore),
      features_after_filter:asNumber(qc?.featuresAfter),
      omitted_features:Number.isFinite(qc?.featuresBefore)&&Number.isFinite(qc?.featuresAfter) ? qc.featuresBefore-qc.featuresAfter : null,
      median_missing_fraction:asNumber(qc?.medianMissingFraction),
      maximum_sample_missing_fraction:maxSampleMissingFraction,
      sample_missingness_profiles_available:missingProfiles.length,
      sample_profiles_with_missing_features:incompleteProfiles.length,
      feature_level_missingness_evaluated:profiles.length>0,
      fdr_family:'BH adjustment within each omics-layer feature family; any UI q<=0.10 threshold is exploratory ranking and not universal FDR certification',
      statistical_calibration_certified:false,
      primary_numeric_source:entry?.numericalProvenance==='exact_native_r_reference_table'
        ? 'complete_native_R_reference_table'
        : 'browser_exploratory_model',
      null_rnaseq_calibration_review_required:
        layer==='transcriptomics' &&
        qc?.preprocessingAudit?.declaredValueType==='raw_counts',
      reference_method_key:backendKey,
      reference_executed:backend?.status==='ok',
      backend_confirmation_status:backend?.status==='ok'?
        'executed_not_statistically_certified':'not_executed',
      missingness_review_required:anyMissing ||
        (qc?.medianMissingFraction||0)>0 ||
        (qc?.differentialMissingness?.flaggedFeatures||0)>0,
      mnar_imputation_applied:(qc?.msQc?.mnar?.imputedValues||0)>0
    };
  }
  const sourceManifest=result.inputManifest||{};
  const runtime=result.referenceBackend?.environment||result.referenceBackend?.runtime||null;
  const status=result.scientificAssurance||{};
  const manifest={
    format:'pmx-omics-methods-audit/v1',
    scientific_certification:false,
    report_generated_at:safe(result.generatedAt)||null,
    engine_version:safe(result.engine.version)||'unknown',
    objective:safe(p.objective)||'unknown',
    intended_inference:safe(p.analysisIntent)||'exploratory',
    experimental_design:safe(p.designType)||'unknown',
    sample_overlap:safe(p.sampleOverlap)||'unknown',
    organism:safe(p.organism)||'unknown',
    outcome_type:safe(p.outcomeType)||'none',
    independent_subjects:asNumber(result.metadataSummary?.subjects),
    covariate_column_names:Array.isArray(p.covariateColumns)?
       p.covariateColumns.map(safe).filter(Boolean):[],
    input_manifest_present:Object.keys(sourceManifest).length>0,
    input_manifest_provenance:'Full file SHA-256 integrity is supplied by the portable project export, not guaranteed by this methods-only file.',
    omics:layers,
    scientific_assurance_status:safe(status?.status)||'not_available',
    clinical_validation:false,
    external_validation_claim:false,
    environment_capture_status:runtime?'partial_runtime_metadata_present':'runtime_lock_not_available',
    reference_runtime_summary:runtime?{
      r_version:safe(runtime?.rVersion||runtime?.r||runtime?.R),
      package_versions:Array.isArray(runtime?.packages)?
        runtime.packages.map((/** @type {any} */ x)=>({name:safe(x?.name),version:safe(x?.version)})):
        null
    }:null,
    reproducibility_limits:[
      'Methods-only audit is not a reproducible raw-data archive.',
      'Use the complete portable project export to preserve original inputs and their SHA-256 hashes; it can contain confidential data.',
      'An R environment lockfile/container and raw-MS vendor binaries are not part of this methods-only manifest.',
      'Sample sizes and QC exclusions must be independently checked by the investigator.',
      'A median of zero missingness can conceal missing observations in individual assays; the report records their count and the worst measured fraction.',
      'DESeq2 raw-count simulations showed excess null discoveries in some tested designs; assess statistical calibration before confirmatory claims.',
      'Reference-model execution, descriptive BH FDR, and public-benchmark recovery never certify a new biological cohort.'
    ]
  };
  for(const record of Object.values(manifest.omics)){
    if(record.omitted_features===null || !Number.isFinite(record.omitted_features))
      record.omitted_features=null;
  }
  return manifest;
}
