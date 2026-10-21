export const MULTIOMICS_DEMO_FILES = Object.freeze({
  metadata: 'demo_metadata.csv',
  transcriptomics: 'demo_transcriptomics.csv',
  proteomics: 'demo_proteomics.csv',
  metabolomics: 'demo_metabolomics.csv'
});

export const MULTIOMICS_DEMO_COLUMN_MAPPING = Object.freeze({
  subject_id: 'subject_id',
  sample_id: 'sample_id',
  assay_id: 'assay_id',
  omic: 'omic',
  condition: 'condition',
  timepoint: 'timepoint',
  batch: 'batch',
  technical_replicate: 'technical_replicate',
  outcome: 'outcome',
  sample_type: '',
  survival_time: '',
  survival_event: ''
});

export const MULTIOMICS_DEMO_UI_STATE = Object.freeze({
  studyName: 'Demo — treatment × time',
  subjectCount: 8,
  groupVariable: 'condition',
  outcome: 'outcome',
  paired: 'yes',
  longitudinal: 'yes',
  objective: 'time',
  organism: 'human',
  unitType: 'participant',
  studySetting: 'clinical_interventional',
  designType: 'repeated',
  groupCount: '2',
  timepointCount: '2',
  sampleOverlap: 'same_specimen',
  technicalReplicatesExpected: 'yes',
  batchKnown: 'yes',
  outcomeType: 'binary',
  outcomeTimepoint: '',
  covariatesAvailable: 'yes',
  partialOmicsExpected: 'no',
  transcriptomicsPlatform: 'bulk_rnaseq',
  transcriptomicsValues: 'raw_counts',
  transcriptomicsIdType: 'gene_symbol',
  proteomicsPlatform: 'label_free',
  proteomicsValues: 'log_intensity',
  proteomicsIdType: 'uniprot',
  metabolomicsPlatform: 'untargeted_lcms',
  metabolomicsValues: 'peak_area',
  metabolomicsIdType: 'chebi',
  metabolomicsIdentificationConfidence: 'mixed_unknown',
  msBlankFilter: 'flag',
  msBlankFold: 5,
  msQcRsdFilter: 'yes',
  msQcRsdThreshold: 0.30,
  msDriftCorrection: 'yes',
  msMnarStrategy: 'none',
  resolveIdentifiers: false,
  useReactome: false,
  referenceBackendMode: 'browser',
  analysisIntent: 'exploratory',
  rnaseqReference: '',
  rnaseqComparison: '',
  validateDiabloHoldout: false
});

export const MULTIOMICS_DEMO_PROTOCOL = Object.freeze({
  organism: MULTIOMICS_DEMO_UI_STATE.organism,
  objective: MULTIOMICS_DEMO_UI_STATE.objective,
  longitudinal: true,
  designType: MULTIOMICS_DEMO_UI_STATE.designType,
  studySetting: MULTIOMICS_DEMO_UI_STATE.studySetting,
  analysisIntent: MULTIOMICS_DEMO_UI_STATE.analysisIntent,
  rnaPrimaryContrast: null,
  unitType: MULTIOMICS_DEMO_UI_STATE.unitType,
  groupCount: MULTIOMICS_DEMO_UI_STATE.groupCount,
  timepointCount: MULTIOMICS_DEMO_UI_STATE.timepointCount,
  sampleOverlap: MULTIOMICS_DEMO_UI_STATE.sampleOverlap,
  technicalReplicatesExpected: MULTIOMICS_DEMO_UI_STATE.technicalReplicatesExpected,
  batchKnown: MULTIOMICS_DEMO_UI_STATE.batchKnown,
  outcomeType: MULTIOMICS_DEMO_UI_STATE.outcomeType,
  outcomeTimepoint: MULTIOMICS_DEMO_UI_STATE.outcomeTimepoint,
  covariatesAvailable: MULTIOMICS_DEMO_UI_STATE.covariatesAvailable,
  covariateColumns: [],
  partialOmicsExpected: MULTIOMICS_DEMO_UI_STATE.partialOmicsExpected,
  metabolomicsIdentificationConfidence: MULTIOMICS_DEMO_UI_STATE.metabolomicsIdentificationConfidence,
  msBlankFilter: MULTIOMICS_DEMO_UI_STATE.msBlankFilter,
  msBlankFold: MULTIOMICS_DEMO_UI_STATE.msBlankFold,
  msQcRsdFilter: MULTIOMICS_DEMO_UI_STATE.msQcRsdFilter,
  msQcRsdThreshold: MULTIOMICS_DEMO_UI_STATE.msQcRsdThreshold,
  msDriftCorrection: MULTIOMICS_DEMO_UI_STATE.msDriftCorrection,
  msMnarStrategy: MULTIOMICS_DEMO_UI_STATE.msMnarStrategy,
  validateDiabloHoldout: MULTIOMICS_DEMO_UI_STATE.validateDiabloHoldout
});

export const MULTIOMICS_DEMO_DATA_TYPES = Object.freeze({
  transcriptomics: MULTIOMICS_DEMO_UI_STATE.transcriptomicsValues,
  proteomics: MULTIOMICS_DEMO_UI_STATE.proteomicsValues,
  metabolomics: MULTIOMICS_DEMO_UI_STATE.metabolomicsValues
});

export const MULTIOMICS_DEMO_IDENTIFIER_TYPES = Object.freeze({
  transcriptomics: MULTIOMICS_DEMO_UI_STATE.transcriptomicsIdType,
  proteomics: MULTIOMICS_DEMO_UI_STATE.proteomicsIdType,
  metabolomics: MULTIOMICS_DEMO_UI_STATE.metabolomicsIdType
});

export const MULTIOMICS_DEMO_EXPECTED = Object.freeze({
  subjects: 8,
  biologicalSamples: 16,
  assays: 49,
  layers: 3,
  conditions: ['control','treatment'],
  timepoints: ['T0','T12']
});
