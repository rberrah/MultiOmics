<script>
  import { onMount } from 'svelte';
  import { resolve, asset } from '$app/paths';
  import { language } from '#lib/stores/language.js';
  import { runDeterministicAnalysis } from '#lib/multiomics/deterministic.js';
  import { buildResearchBundle, fullLayerCsv, referencePreviewAudit } from '#lib/multiomics/research-export-package.js';
  import { parseIndividualSignatureCsv } from '#lib/multiomics/downstream-exchange.js';
  import { auditPairedRnaStudy, pairedExampleSampleSheetCsv } from '#lib/multiomics/paired-study-client.js';
  import MultiomicsVisualizations from '#lib/components/MultiomicsVisualizations.svelte';
  import ExternalValidationPanel from '#lib/components/ExternalValidationPanel.svelte';
  import { buildMultiomicsVisualizationData } from '#lib/multiomics/visualization-data.js';
  import { projectNativeRNumbers } from '#lib/multiomics/r-primary-results.js';
  import { assessConfirmatoryResult } from '#lib/multiomics/confirmatory-result-gate.js';
  import { convertMsAucExport } from '#lib/multiomics/ms-auc-import.js';
  import { parseMetaboliteAnnotations } from '#lib/multiomics/metabolite-annotation.js';
  import { assessScientificAssurance } from '#lib/multiomics/scientific-assurance.js';
  import { mofa2BackendCompatibility,validateMofa2MethodResult } from '#lib/multiomics/reference-backend-integrity.js';
  import { evaluateReferenceBackendHealth } from '#lib/multiomics/reference-backend-preflight.js';
  import { createPortableProject, parsePortableProject, PORTABLE_PROJECT_MAX_BYTES } from '#lib/multiomics/portable-project.js';
  import { buildResearchMethodsManifest } from '#lib/multiomics/research-report.js';

  // A server-rendered button can be visible before its Svelte click handler
  // is attached. Public-site smoke tests must wait for actual hydration,
  // rather than treating visible SSR HTML as an operational application.
  onMount(() => {
    document.documentElement.dataset.pmxMultiomicsHydrated = 'true';
    return () => { delete document.documentElement.dataset.pmxMultiomicsHydrated; };
  });

  /** @param {string} fr @param {string} en */
  const t = (fr, en) => $language === 'en' ? en : fr;

  /** @type {'explore' | 'groups' | 'outcome' | 'time'} */
  let objective = 'explore';
  let organism = 'human';
  let studyName = '';
  let groupVariable = '';
  let outcome = '';
  let unitType = 'participant';
  let studySetting = 'clinical_observational';
  let designType = 'independent';
  let groupCount = '2';
  let timepointCount = '1';
  let sampleOverlap = 'same_specimen';
  let technicalReplicatesExpected = 'unknown';
  let batchKnown = 'unknown';
  let outcomeType = 'none';
  let outcomeTimepoint = '';
  let covariatesAvailable = 'yes';
  let partialOmicsExpected = 'no';
  /** @type {string[]} */
  let selectedCovariates = [];
  let demoLoaded = false;
  let demoLoading = false;
  let demoError = '';

  /** The simple controls and the expert form share the exact same analysis state. */
  /** @param {Event} event */
  function chooseSimpleStudySetting(event) {
    const selected = event.currentTarget;
    if (!(selected instanceof HTMLSelectElement)) return;
    studySetting = selected.value;
    unitType = studySetting === 'animal' ? 'animal'
      : ['cell', 'organoid'].includes(studySetting) ? 'culture'
      : studySetting === 'other' ? 'other' : 'participant';
  }

  /** @param {Event} event */
  function chooseSimpleDesign(event) {
    const value = event.currentTarget;
    if (!(value instanceof HTMLSelectElement)) return;
    designType = value.value;
    longitudinal = designType === 'repeated' ? 'yes' : 'no';
    paired = designType === 'paired' ? 'yes' : 'no';
    if (designType === 'repeated' && timepointCount === '1') timepointCount = '2';
  }

  /** @param {Event} event */
  function chooseSimpleObjective(event) {
    const value = event.currentTarget;
    if (!(value instanceof HTMLSelectElement)) return;
    objective = /** @type {'explore' | 'groups' | 'outcome' | 'time'} */ (value.value);
    if (objective === 'time' && designType === 'independent') {
      designType = 'repeated';
      longitudinal = 'yes';
      if (timepointCount === '1') timepointCount = '2';
    }
  }
  let analysisStatus = 'idle';
  let analysisError = '';
  let portableProjectBusy = false;
  let portableProjectMessage = '';
  let portableProjectError = '';
  let researchBundleMessage = '';
  let researchBundleError = '';
  // Optional user-supplied scores, computed *outside* the default browser
  // workflow using frozen/out-of-fold models. No automatic signature fitting.
  /** @type {Array<Record<string,any>>} */
  let individualSignatureRows = [];
  let individualSignatureFilename = '';
  let individualSignatureError = '';
  /** @type {any} */
  let analysisResult = null;
  /** @param {any} note */
  function isMissingnessRiskNote(note) {
    return String(note?.code ?? '').endsWith('_missingness_not_ignorable');
  }
  let helpTooltip = { visible: false, text: '', left: 0, top: 0, placement: 'above' };
  // External API requests and reference R data forwarding are opt-in.
  let useReactome = false;
  let resolveIdentifiers = false;
  let referenceBackendMode = 'browser';
  let analysisIntent = 'exploratory';
  // Explicit numerator and denominator for a prespecified RNA-seq contrast.
  // Never infer its scientific direction from alphabetical condition labels.
  let rnaseqReference = '';
  let rnaseqComparison = '';
  let validateDiabloHoldout = false;
  let referenceBackendUrl = 'http://127.0.0.1:8787';
  let referenceBackendStatus = 'unchecked';
  let referenceBackendMessage = '';
  /** @type {any} */
  let referenceBackendHealth = null;
  // Exact R request/response captured only when the local R backend actually fits a method.
  // Kept in memory, exported exclusively in the researcher-initiated confidential ZIP.
  /** @type {{requestJson:string,responseJson:string}|null} */
  let referenceReplay = null;

  let transcriptomicsPlatform = 'bulk_rnaseq';
  let transcriptomicsValues = 'raw_counts';
  let transcriptomicsIdType = 'ensembl_gene';
  let proteomicsPlatform = 'label_free';
  let proteomicsValues = 'lfq_intensity';
  let proteomicsIdType = 'uniprot';
  let metabolomicsPlatform = 'untargeted_lcms';
  let metabolomicsValues = 'peak_area';
  let msImportMessage = '';
  /** State is observable independently of the translated message. */
  let msImportState = 'idle'; // idle | loading | parsed | ready | error
  let msImportSequence = 0;
  /** Original vendor CSV/TSV remains unchanged for provenance and re-export. */
  /** @type {File|null} */
  let msOriginalFile = null;
  let msAnnotationCsv = '';
  /** @type {string[]} */
  let msAnnotationFields = [];
  let msImportFormat = '';
  /** @type {File|null} */
  let msUserAnnotationFile = null;
  let confirmIndependentAssays = false;
  let exploratorySheetGenerated = false;
  /** @type {{feature:string,chebi:string,status:string,provenance:string}[]} */
  let metaboliteAnnotations = [];
  let annotationFileName = '';
  let annotationError = '';
  let msAutoMetadata = false;
  let metabolomicsIdType = 'chebi';
  let msBlankFilter = 'flag';
  let msBlankFold = 5;
  let msQcRsdFilter = 'yes';
  let msQcRsdThreshold = 0.30;
  let msDriftCorrection = 'yes';
  let msMnarStrategy = 'none';

  /** @type {number | undefined} */
  let subjectCount;
  let paired = 'no';
  let longitudinal = 'no';

  /** @type {Record<'transcriptomics' | 'proteomics' | 'metabolomics' | 'metadata', File | null>} */
  let files = {
    transcriptomics: null,
    proteomics: null,
    metabolomics: null,
    metadata: null
  };

  /** @type {string[]} */
  let metadataHeaders = [];
  /** @type {Record<string, string>[]} */
  let metadataRows = [];
  let metadataDelimiter = '';
  let metadataError = '';

  /** @type {Record<'transcriptomics' | 'proteomics' | 'metabolomics', {headers:string[], sampleIds:string[], rowIds:string[], featureColumn:string, error:string}>} */
  let matrixInfo = {
    transcriptomics: { headers: [], sampleIds: [], rowIds: [], featureColumn: '', error: '' },
    proteomics: { headers: [], sampleIds: [], rowIds: [], featureColumn: '', error: '' },
    metabolomics: { headers: [], sampleIds: [], rowIds: [], featureColumn: '', error: '' }
  };

  /** @param {Event} event */
  function showGlobalHelp(event) {
    const raw = event.target;
    if (!(raw instanceof Element)) return;
    const target = raw.closest('.help-tip');
    if (!target) return;
    const text = target.getAttribute('data-tooltip') || '';
    if (!text) return;
    const rect = target.getBoundingClientRect();
    const halfWidth = Math.min(170, Math.max(120, window.innerWidth * 0.38));
    const center = rect.left + rect.width / 2;
    const left = Math.min(window.innerWidth - halfWidth - 8, Math.max(halfWidth + 8, center));
    const above = rect.top > 150;
    helpTooltip = {
      visible: true,
      text,
      left,
      top: above ? rect.top - 9 : rect.bottom + 9,
      placement: above ? 'above' : 'below'
    };
  }

  /** @param {Event} event */
  function hideGlobalHelp(event) {
    const raw = event.target;
    if (!(raw instanceof Element) || !raw.closest('.help-tip')) return;
    helpTooltip = { ...helpTooltip, visible: false };
  }

  const omicLayers = /** @type {const} */ (['transcriptomics', 'proteomics', 'metabolomics']);
  /** @param {string} layer */
  function omicLabel(layer) {
    if (layer === 'transcriptomics') return t('Transcriptomique', 'Transcriptomics');
    if (layer === 'proteomics') return t('Protéomique', 'Proteomics');
    if (layer === 'metabolomics') return t('Métabolomique', 'Metabolomics');
    return layer;
  }

  /** @type {Record<string, [string,string]>} */
  const fieldLabels = {
    subject_id: ['Sujet, animal ou culture', 'Subject, animal or culture'],
    sample_id: ['Prélèvement', 'Sample'],
    assay_id: ['Nom de la mesure dans le fichier', 'Measurement name in the file'],
    omic: ['Type de mesure', 'Measurement type'],
    condition: ['Groupe (ex. traité / témoin)', 'Group (e.g. treated / control)'],
    timepoint: ['Moment du prélèvement', 'Sampling time'],
    batch: ['Lot de mesures', 'Measurement batch'],
    technical_replicate: ['Mesure répétée au laboratoire', 'Repeated laboratory measurement'],
    outcome: ['Résultat étudié', 'Outcome of interest'],
    survival_time: ['Durée du suivi', 'Follow-up duration'],
    survival_event: ['Événement observé (0/1)', 'Observed event (0/1)'],
    sample_type: ['Type d’échantillon analysé', 'Sample type'],
    injection_order: ['Ordre de passage dans l’appareil', 'Instrument run order']
  };
  /** @param {string} key */
  const fieldLabel = (key) => {
    const labels = fieldLabels[key] || [key, key];
    return t(labels[0], labels[1]);
  };

  const fieldDefinitions = [
    {
      key: 'subject_id',
      label: 'Subject / experimental unit',
      required: true,
      aliases: ['subject_id', 'subject', 'patient_id', 'patient', 'participant_id', 'participant', 'individual_id', 'individual', 'donor_id', 'animal_id', 'unit_id', 'biological_unit', 'id_patient', 'id_sujet', 'sujet', 'participant_id', 'individu', 'id_individu', 'id_animal']
    },
    {
      key: 'sample_id',
      label: 'Biological sample',
      required: true,
      aliases: ['sample_id', 'sample', 'specimen_id', 'specimen', 'biospecimen_id', 'biosample_id', 'sample_name', 'id_echantillon', 'echantillon', 'id_prelevement', 'prelevement', 'nom_echantillon']
    },
    {
      key: 'assay_id',
      label: 'Assay / run identifier',
      required: true,
      aliases: ['assay_id', 'assay', 'run_id', 'run', 'measurement_id', 'library_id', 'file_id', 'assay_name', 'id_assay', 'id_mesure', 'mesure_id', 'id_run', 'id_fichier']
    },
    {
      key: 'omic',
      label: 'Omics layer',
      required: true,
      aliases: ['omic', 'omics', 'layer', 'modality', 'data_type', 'datatype', 'assay_type', 'omics_type', 'omique', 'type_omique', 'couche', 'modalite']
    },
    {
      key: 'condition',
      label: 'Condition / group',
      required: false,
      aliases: ['condition', 'group', 'treatment', 'arm', 'cohort', 'status', 'class', 'phenotype_group', 'condition_name', 'groupe', 'traitement', 'bras', 'cohorte', 'statut', 'condition_experimentale']
    },
    {
      key: 'timepoint',
      label: 'Time point',
      required: false,
      aliases: ['timepoint', 'time_point', 'time', 'visit', 'visit_id', 'visit_name', 'day', 'week', 'temps', 'temps_j', 'jour', 'semaine', 'visite', 'id_visite']
    },
    {
      key: 'batch',
      label: 'Technical batch',
      required: false,
      aliases: ['batch', 'batch_id', 'plate', 'run_batch', 'technical_batch', 'batch_name', 'lot', 'id_lot', 'plaque', 'batch_technique']
    },
    {
      key: 'technical_replicate',
      label: 'Technical replicate',
      required: false,
      aliases: ['technical_replicate', 'technical_rep', 'tech_rep', 'replicate', 'replicate_id', 'rep', 'replicat_technique', 'replicat', 'rep_technique', 'numero_replicat']
    },
    {
      key: 'outcome',
      label: 'Outcome / endpoint',
      required: false,
      aliases: ['outcome', 'endpoint', 'response', 'label', 'target', 'phenotype', 'clinical_outcome', 'reponse', 'critere', 'critere_jugement', 'phenotype_clinique', 'evenement']
    },
    {
      key: 'survival_time',
      label: 'Survival / follow-up time',
      required: false,
      aliases: ['survival_time', 'followup_time', 'follow_up_time', 'time_to_event', 'tte', 'duree_suivi', 'temps_survie', 'temps_evenement']
    },
    {
      key: 'survival_event',
      label: 'Survival event',
      required: false,
      aliases: ['survival_event', 'event', 'event_status', 'death', 'censor', 'status_event', 'evenement_survie', 'deces', 'censure']
    },
    {
      key: 'sample_type',
      label: 'Injection / sample type',
      required: false,
      aliases: ['sample_type','sample_class','sample_role','injection_type','qc_type','type_echantillon','type_injection','nature_echantillon']
    },
    {
      key: 'injection_order',
      label: 'MS injection order',
      required: false,
      aliases: ['injection_order','run_order','sequence_order','injection_number','order','ordre_injection','ordre_passage','numero_injection']
    }
  ];

  /** @type {Record<string,string>} */
  let columnMapping = Object.fromEntries(fieldDefinitions.map((field) => [field.key, '']));

  const objectives = {
    explore: {
      titleFr: 'Explorer la structure multi-omique partagée',
      titleEn: 'Explore the shared multi-omics structure',
      methodFr: 'ACP multi-blocs équilibrée',
      methodEn: 'Balanced multi-block PCA',
      detailFr: 'Les variables les plus informatives de chaque couche sont standardisées et pondérées par bloc afin d’extraire des axes latents communs sans variable cible.',
      detailEn: 'Top-variable features from each layer are standardized and block-balanced to extract shared latent axes without a target variable.'
    },
    groups: {
      titleFr: 'Comparer des groupes biologiques',
      titleEn: 'Compare biological groups',
      methodFr: 'Inférence déterministe ajustée',
      methodEn: 'Adjusted deterministic inference',
      detailFr: 'Ajustement explicite des batches/covariables, contrastes adaptés au design, BH-FDR, relations inter-omiques et convergence Reactome.',
      detailEn: 'Explicit batch/covariate adjustment, design-aware contrasts, BH-FDR, cross-omics relationships and Reactome convergence.'
    },
    outcome: {
      titleFr: 'Expliquer un outcome clinique ou expérimental',
      titleEn: 'Explain a clinical or experimental outcome',
      methodFr: 'Régression déterministe selon le type d’outcome',
      methodEn: 'Outcome-specific deterministic regression',
      detailFr: 'Régression linéaire, logistique, Poisson, ANCOVA multiclasse ou Cox selon l’outcome, plus validation prédictive nested-CV lorsque possible, avec covariables explicites et BH-FDR.',
      detailEn: 'Linear, logistic, Poisson, multiclass ANCOVA or Cox models are selected from the declared outcome, with nested-CV prediction when applicable, explicit covariates and BH-FDR.'
    },
    time: {
      titleFr: 'Étudier l’évolution au cours du temps',
      titleEn: 'Describe change over time',
      methodFr: 'Workflow longitudinal ajusté au design',
      methodEn: 'Design-aware longitudinal workflow',
      detailFr: 'Les mesures répétées utilisent un modèle à intercept aléatoire sujet avec interaction condition × temps, batch et covariables explicites.',
      detailEn: 'Repeated measurements use a subject random-intercept model with condition × time interaction and explicit batch/covariates.'
    }
  };

  const omicAliases = {
    transcriptomics: ['transcriptomics', 'transcriptome', 'transcriptomique', 'rna', 'rnaseq', 'rna_seq', 'gene_expression', 'expression_genique', 'mrna', 'arn'],
    proteomics: ['proteomics', 'proteome', 'proteomique', 'protein', 'proteins', 'proteine', 'proteines', 'lfq'],
    metabolomics: ['metabolomics', 'metabolome', 'metabolomique', 'metabolite', 'metabolites', 'met']
  };

  const databaseRegistry = [
    {
      name: 'Ensembl',
      scope: 'genes',
      scopeFr: 'gènes',
      role: 'Resolve Ensembl IDs and gene symbols, confirm species and canonical gene annotations.',
      roleFr: 'Résoudre les identifiants Ensembl et symboles de gènes, confirmer l’espèce et les annotations canoniques.',
      status: 'live when identifier resolution is enabled',
      statusFr: 'actif quand la résolution d’identifiants est activée',
      url: 'https://rest.ensembl.org/documentation/'
    },
    {
      name: 'UniProt',
      scope: 'proteins',
      scopeFr: 'protéines',
      role: 'Resolve protein accessions and cross-reference proteins to genes, Ensembl, Reactome and other resources.',
      roleFr: 'Résoudre les accessions protéiques et les relier aux gènes, à Ensembl, Reactome et d’autres ressources.',
      status: 'live when identifier resolution is enabled',
      statusFr: 'actif quand la résolution d’identifiants est activée',
      url: 'https://www.uniprot.org/help/id_mapping'
    },
    {
      name: 'ChEBI',
      scope: 'metabolites',
      scopeFr: 'métabolites',
      role: 'Resolve curated chemical entities, synonyms, structures and ontology relationships.',
      roleFr: 'Résoudre les entités chimiques curées, synonymes, structures et relations ontologiques.',
      status: 'live in the current engine',
      statusFr: 'actif dans le moteur actuel',
      url: 'https://www.ebi.ac.uk/chebi/tools'
    },
    {
      name: 'UniChem',
      scope: 'metabolites',
      scopeFr: 'métabolites',
      role: 'Cross-reference chemical identifiers between databases when the uploaded metabolite identifier is not ChEBI.',
      roleFr: 'Croiser les identifiants chimiques entre bases lorsque l’identifiant métabolite importé n’est pas ChEBI.',
      status: 'live for deterministic InChIKey → ChEBI resolution',
      statusFr: 'actif pour la résolution déterministe InChIKey → ChEBI',
      url: 'https://www.ebi.ac.uk/unichem/'
    },
    {
      name: 'KEGG',
      scope: 'pathways · optional',
      scopeFr: 'voies · optionnel',
      role: 'Optional academic-use cross-reference and pathway/reaction annotation. KEGG REST access is rate-limited and is not used as the core engine of the public MVP.',
      roleFr: 'Référence croisée et annotation de voies/réactions optionnelles pour usage académique. KEGG REST n’est pas le moteur central du prototype public.',
      status: 'optional registry target · not called automatically yet',
      statusFr: 'connecteur optionnel prévu · non appelé automatiquement',
      url: 'https://www.kegg.jp/kegg/rest/'
    },
    {
      name: 'Reactome',
      scope: 'pathways',
      scopeFr: 'voies',
      role: 'Map genes, proteins and ChEBI entities onto common pathways and reactions for integrated pathway interpretation.',
      roleFr: 'Mapper gènes, protéines et entités ChEBI sur des voies et réactions communes pour l’interprétation intégrée.',
      status: 'live in the current engine',
      statusFr: 'actif dans le moteur actuel',
      url: 'https://reactome.org/dev/analysis'
    },
    {
      name: 'STRING',
      scope: 'network',
      scopeFr: 'réseau',
      role: 'Build compact protein interaction modules after identifier resolution; not used as a substitute for the measured data.',
      roleFr: 'Construire des modules compacts d’interactions protéiques après résolution des identifiants, sans remplacer les données mesurées.',
      status: 'registry target · not called automatically yet',
      statusFr: 'connecteur prévu · non appelé automatiquement',
      url: 'https://string-db.org/help/api/'
    }
  ];

  /** @param {unknown} value */
  function normalise(value) {
    return String(value ?? '')
      .trim()
      .toLowerCase()
      .normalize('NFD')
      .replace(/[\u0300-\u036f]/g, '')
      .replace(/[^a-z0-9]+/g, '_')
      .replace(/^_|_$/g, '');
  }

  /**
   * @param {string} line
   * @param {string} delimiter
   */
  function splitDelimitedLine(line, delimiter) {
    const out = [];
    let current = '';
    let quoted = false;
    for (let i = 0; i < line.length; i += 1) {
      const char = line[i];
      if (char === '"') {
        if (quoted && line[i + 1] === '"') {
          current += '"';
          i += 1;
        } else {
          quoted = !quoted;
        }
      } else if (char === delimiter && !quoted) {
        out.push(current.trim());
        current = '';
      } else {
        current += char;
      }
    }
    out.push(current.trim());
    return out;
  }

  /** @param {string} text */
  function detectDelimiter(text) {
    const line = text.split(/\r?\n/).find((row) => row.trim()) ?? '';
    const candidates = [',', '\t', ';'];
    return candidates
      .map((delimiter) => ({ delimiter, count: splitDelimitedLine(line, delimiter).length }))
      .sort((a, b) => b.count - a.count)[0]?.delimiter ?? ',';
  }

  /** @param {string} text */
  function parseTable(text) {
    const delimiter = detectDelimiter(text);
    const lines = text.split(/\r?\n/).filter((line) => line.trim());
    if (!lines.length) return { delimiter, headers: [], rows: [] };
    const headers = splitDelimitedLine(lines[0], delimiter);
    const rows = lines.slice(1).map((line) => {
      const values = splitDelimitedLine(line, delimiter);
      return Object.fromEntries(headers.map((header, index) => [header, values[index] ?? '']));
    });
    return { delimiter, headers, rows };
  }

  /** @param {string[]} headers */
  function autoMapColumns(headers) {
    /** @type {Record<string,string>} */
    const next = {};
    const normalisedHeaders = headers.map((header) => ({ raw: header, normalised: normalise(header) }));
    for (const field of fieldDefinitions) {
      const aliases = new Set(field.aliases.map(normalise));
      const matches = normalisedHeaders.filter((header) => aliases.has(header.normalised));
      next[field.key] = matches.length === 1 ? matches[0].raw : '';
    }
    columnMapping = next;
  }

  /**
   * @param {File | null} file
   * @param {() => boolean} [isCurrent]
   */
  async function inspectMetadata(file, isCurrent = () => true) {
    if (!isCurrent()) return;
    metadataError = '';
    metadataHeaders = [];
    metadataRows = [];
    if (!file) return;
    try {
      const text = await file.text();
      if (!isCurrent()) return;
      const parsed = parseTable(text);
      metadataDelimiter = parsed.delimiter === '\t' ? 'tab' : parsed.delimiter === ';' ? 'semicolon' : 'comma';
      metadataHeaders = parsed.headers;
      metadataRows = parsed.rows;
      autoMapColumns(parsed.headers);
      if (!parsed.headers.length || !parsed.rows.length) metadataError = 'The metadata file appears empty.';
    } catch (error) {
      if (!isCurrent()) return;
      metadataError = error instanceof Error ? error.message : 'Unable to read metadata.';
    }
  }

  /**
   * @param {'transcriptomics' | 'proteomics' | 'metabolomics'} layer
   * @param {File | null} file
   * @param {() => boolean} [isCurrent]
   */
  async function inspectMatrix(layer, file, isCurrent = () => true) {
    if (!isCurrent()) return;
    if (!file) {
      matrixInfo = { ...matrixInfo, [layer]: { headers: [], sampleIds: [], rowIds: [], featureColumn: '', error: '' } };
      return;
    }
    try {
      const text = await file.text();
      if (!isCurrent()) return;
      const parsed = parseTable(text);
      const featureColumn = parsed.headers[0] ?? '';
      const sampleIds = parsed.headers.slice(1);
      const rowIds = parsed.rows.map((row) => row[featureColumn] ?? '').filter(Boolean);
      matrixInfo = {
        ...matrixInfo,
        [layer]: {
          headers: parsed.headers,
          sampleIds,
          rowIds,
          featureColumn,
          error: sampleIds.length ? '' : 'No assay columns detected.'
        }
      };
    } catch (error) {
      if (!isCurrent()) return;
      matrixInfo = {
        ...matrixInfo,
        [layer]: {
          headers: [],
          sampleIds: [],
          rowIds: [],
          featureColumn: '',
          error: error instanceof Error ? error.message : 'Unable to read matrix.'
        }
      };
    }
  }

  /**
   * @param {'transcriptomics' | 'proteomics' | 'metabolomics' | 'metadata'} layer
   * @param {Event} event
   */
  async function selectFile(layer, event) {
    const input = /** @type {HTMLInputElement} */ (event.currentTarget);
    const file = input.files?.[0] ?? null;
    demoLoaded = false;
    analysisResult = null;
    analysisStatus = 'idle';
    analysisError = '';
    if (layer !== 'metadata') {
      confirmIndependentAssays = false;
      exploratorySheetGenerated = false;
    }
    const currentMsImport = layer === 'metabolomics' ? ++msImportSequence : msImportSequence;
    const isCurrentMsImport = () => currentMsImport === msImportSequence;
    if (layer === 'metabolomics') {
      // Do not allow a prior matrix to be analyzed during an in-flight replacement.
      files = { ...files, metabolomics: null };
      msImportState = file ? 'loading' : 'idle';
      msOriginalFile = file;
      msAnnotationCsv = '';
      msAnnotationFields = [];
      msImportFormat = '';
      msImportMessage = '';
    }
    if (layer === 'metadata') {
      msAutoMetadata = false;
      exploratorySheetGenerated = false;
      files = { ...files, metadata: file };
      await inspectMetadata(file);
    } else if (layer === 'metabolomics' && file) {
      try {
        const sourceText = await file.text();
        if (currentMsImport !== msImportSequence) return;
        const converted = convertMsAucExport(sourceText);
        msImportState = 'parsed';
        msAnnotationCsv = converted.annotationsCsv || '';
        msAnnotationFields = converted.annotationColumns || [];
        msImportFormat = converted.format;
        const matrix = ['long_ms_auc','wide_ms_auc'].includes(converted.format)
          ? new File([converted.matrixCsv], file.name.replace(/\.[^.]+$/, '') + '_matrix.csv', { type:'text/csv' })
          : file;
        const shouldUpdateMetadata = converted.metadataCsv && (!files.metadata || msAutoMetadata);
        const autoMeta = shouldUpdateMetadata
          ? new File([converted.metadataCsv || ''], 'metadata_from_ms_auc.csv', { type:'text/csv' }) : null;
        const clearStaleAutoMetadata = msAutoMetadata && !converted.metadataCsv;

        files = {
          ...files,
          metabolomics: matrix,
          ...autoMeta
            ? { metadata: autoMeta }
            : clearStaleAutoMetadata ? { metadata: null } : {}
        };

        if (clearStaleAutoMetadata) {
          msAutoMetadata = false;
          await inspectMetadata(null, isCurrentMsImport);
          if (!isCurrentMsImport()) return;
        }
        if (autoMeta) {
          msAutoMetadata = true;
          await inspectMetadata(autoMeta, isCurrentMsImport);
          if (!isCurrentMsImport()) return;
          if (metadataError) throw new Error(metadataError);
          if (!files.transcriptomics && !files.proteomics && objective === 'explore') {
            // Only one measured modality: select an identifiable group contrast
            // rather than presenting single-layer MS as multi-block exploration.
            objective = 'groups';
            designType = 'independent';
            longitudinal = 'no';
            groupCount = String(new Set(metadataRows
              .filter((row)=>!['qc','blank'].includes(String(row.sample_type||'').toLowerCase()))
              .map((row)=>row.condition).filter(Boolean)).size || 2);
          }
        }
        if (converted.format !== 'matrix') metabolomicsValues = 'peak_area';
        msImportMessage = converted.format === 'long_ms_auc'
          ? t('Aires de pics MS importées : ', 'MS peak areas imported: ')
            + converted.features + t(' signaux MS × ', ' MS features × ') + converted.assays
            + t(' injections. ', ' injections. ')
            + (autoMeta ? t('Groupes lus dans le fichier.', 'Groups read from the file.')
              : t('Conservez ou fournissez votre tableau des échantillons.', 'Keep or upload sample metadata.'))
          : converted.format === 'wide_ms_auc'
            ? t('Aires de pics MS importées : les colonnes m/z, RT et annotations sont conservées séparément et non utilisées comme intensités. Fournissez les métadonnées des échantillons.', 'MS peak areas imported: m/z, RT and annotations are retained separately and not treated as intensities. Upload sample metadata.')
            : t('Matrice chargée directement. Vérifiez si ses valeurs représentent réellement des aires de pics ou des intensités déjà transformées.', 'Matrix loaded directly. Verify whether values represent raw peak areas or already transformed intensities.');
        await inspectMatrix(layer, matrix, isCurrentMsImport);
        if (!isCurrentMsImport()) return;
        if (matrixInfo.metabolomics?.error) throw new Error(matrixInfo.metabolomics.error);
        msImportState = 'ready';
      } catch (error) {
        if (currentMsImport !== msImportSequence) return;
        msImportState = 'error';
        const clearAutoMetadata = msAutoMetadata && files.metadata?.name === 'metadata_from_ms_auc.csv';

        files = {
          ...files,
          metabolomics: null,
          ...clearAutoMetadata ? { metadata: null } : {}
        };

        if (clearAutoMetadata) {
          msAutoMetadata = false;
          await inspectMetadata(null, isCurrentMsImport);
          if (!isCurrentMsImport()) return;
        }
        msOriginalFile = null;
        msAnnotationCsv = '';
        msAnnotationFields = [];
        msImportFormat = '';
        msImportMessage = '';
        analysisError = error instanceof Error ? error.message : 'Invalid MS AUC export';
        await inspectMatrix(layer, null, isCurrentMsImport);
      }
    } else if (layer === 'metabolomics' && !file) {
      const clearAutoMetadata = msAutoMetadata && files.metadata?.name === 'metadata_from_ms_auc.csv';
      if (clearAutoMetadata) {
        files = { ...files, metadata: null };
        msAutoMetadata = false;
        await inspectMetadata(null, isCurrentMsImport);
      }
      await inspectMatrix(layer, null, isCurrentMsImport);
    } else {
      files = { ...files, [layer]: file };
      await inspectMatrix(layer, file);
    }
  }

  /** @param {Event} event */
  async function selectMetaboliteAnnotations(event) {
    const input = /** @type {HTMLInputElement} */ (event.currentTarget);
    const file = input.files?.[0] || null;
    analysisResult = null;
    analysisStatus = 'idle';
    annotationFileName = '';
    annotationError = '';
    metaboliteAnnotations = [];
    msUserAnnotationFile = file;
    if (!file) return;
    try {
      metaboliteAnnotations = parseMetaboliteAnnotations(await file.text());
      annotationFileName = file.name;
    } catch (error) {
      annotationError = error instanceof Error ? error.message : String(error);
    }
  }

  /** @param {'transcriptomics' | 'proteomics' | 'metabolomics' | null} onlyLayer */
  async function loadDemo(onlyLayer = null) {
    /** @type {Array<['metadata'|'transcriptomics'|'proteomics'|'metabolomics', string]>} */
    const demoFiles = [
      ['metadata', 'demo_metadata.csv'],
      ['transcriptomics', 'demo_transcriptomics.csv'],
      ['proteomics', 'demo_proteomics.csv'],
      ['metabolomics', 'demo_metabolomics.csv']
    ];
    /** @type {Record<'metadata'|'transcriptomics'|'proteomics'|'metabolomics', File | null>} */
    const loaded = {
      metadata: null,
      transcriptomics: null,
      proteomics: null,
      metabolomics: null
    };
    // Static file paths are explicitly declared to Kit 3's AssetPath types.
    // No ad-hoc string cast; each path is validated against the static tree.
    /** @type {Record<string, string>} */
    const demoAssets = {
      'demo_metadata.csv': asset('multiomics/demo_metadata.csv'),
      'demo_transcriptomics.csv': asset('multiomics/demo_transcriptomics.csv'),
      'demo_proteomics.csv': asset('multiomics/demo_proteomics.csv'),
      'demo_metabolomics.csv': asset('multiomics/demo_metabolomics.csv')
    };
    for (const [layer, filename] of demoFiles) {
      if (onlyLayer && layer !== 'metadata' && layer !== onlyLayer) continue;
      const response = await fetch(demoAssets[filename]);
      if (!response.ok) throw new Error('Demo file unavailable: ' + filename);
      let text = await response.text();
      if (layer === 'metadata' && onlyLayer) {
        const parsed = parseTable(text);
        const rows = parsed.rows.filter((row) => row.omic === onlyLayer);
        /** @param {unknown} value */
        const quote = (value) => '"' + String(value ?? '').replaceAll('"','""') + '"';
        text = [parsed.headers.map(quote).join(','), ...rows.map((row) =>
          parsed.headers.map((key) => quote(row[key])).join(','))].join('\n') + '\n';
      }
      loaded[layer] = new File([text], filename, { type: 'text/csv' });
    }
    files = {
      metadata: loaded.metadata,
      transcriptomics: loaded.transcriptomics,
      proteomics: loaded.proteomics,
      metabolomics: loaded.metabolomics
    };
    studyName = onlyLayer ? 'Demo — single ' + onlyLayer : 'Demo — treatment × time';
    confirmIndependentAssays = false;
    exploratorySheetGenerated = false;
    msOriginalFile = null;
    msAnnotationCsv = '';
    msAnnotationFields = [];
    msImportFormat = '';
    msUserAnnotationFile = null;
    metaboliteAnnotations = [];
    annotationFileName = '';
    annotationError = '';
    subjectCount = 8;
    groupVariable = 'condition';
    outcome = 'outcome';
    paired = 'yes';
    longitudinal = 'yes';
    objective = 'time';
    studySetting = 'clinical_interventional';
    designType = 'repeated';
    groupCount = '2';
    timepointCount = '2';
    sampleOverlap = 'same_specimen';
    technicalReplicatesExpected = 'yes';
    batchKnown = 'yes';
    outcomeType = 'binary';
    partialOmicsExpected = 'no';
    transcriptomicsValues = 'raw_counts';
    transcriptomicsIdType = 'gene_symbol';
    proteomicsValues = 'log_intensity';
    proteomicsIdType = 'uniprot';
    metabolomicsValues = 'peak_area';
    metabolomicsIdType = 'chebi';
    resolveIdentifiers = false;
    // One-click demos must remain usable without network access or third-party
    // API availability. Reactome enrichment remains a separate explicit opt-in.
    useReactome = false;
    referenceBackendMode = 'browser';
    analysisIntent = 'exploratory';
    rnaseqReference = '';
    rnaseqComparison = '';
    demoLoaded = true;
    await inspectMetadata(files.metadata);
    await inspectMatrix('transcriptomics', files.transcriptomics);
    await inspectMatrix('proteomics', files.proteomics);
    await inspectMatrix('metabolomics', files.metabolomics);
    await runAnalysis();
  }

  /** @param {'transcriptomics' | 'proteomics' | 'metabolomics' | null} onlyLayer */
  async function loadDemoAndShowResults(onlyLayer = null) {
    if (demoLoading) return;
    demoLoading = true;
    demoError = '';
    analysisError = '';
    analysisResult = null;
    try {
      await loadDemo(onlyLayer);
      // runAnalysis reports some failures in analysisError without throwing.
      if (!analysisResult || analysisStatus !== 'done') {
        throw new Error(analysisError || 'No complete analysis results were produced.');
      }
      // The demonstration is an actual analysis, not a mockup.
      document.getElementById('analysis-results')?.scrollIntoView({ behavior: 'smooth', block: 'start' });
    } catch (error) {
      demoLoaded = false;
      analysisResult = null;
      analysisStatus = 'error';
      const detail = error instanceof Error ? error.message : String(error);
      demoError = t(
        'La démonstration n’a pas pu se terminer. Réessayez ; si le problème persiste, vérifiez la connexion et les fichiers de démonstration. Détail : ',
        'The demo could not finish. Try again; if the problem continues, check your connection and demo files. Detail: '
      ) + detail;
      analysisError = demoError;
    } finally {
      demoLoading = false;
    }
  }

  /**
   * @param {Record<string,string>} row
   * @param {string} key
   */
  function mappedValue(row, key) {
    const column = columnMapping[key];
    return column ? row[column] ?? '' : '';
  }

  /** @param {unknown} value */
  function canonicalOmic(value) {
    const key = normalise(value);
    for (const [omic, aliases] of Object.entries(omicAliases)) {
      if (aliases.map(normalise).includes(key)) return omic;
    }
    return '';
  }

  /**
   * @param {Record<string,string>[]} rows
   * @param {Record<string,string>} mapping
   * @param {string} key
   */
  function uniqueMappedFrom(rows, mapping, key) {
    const column = mapping[key];
    return new Set(rows.map((row) => column ? row[column] ?? '' : '').filter(Boolean));
  }

  /**
   * @param {Record<string,string>[]} rows
   * @param {Record<string,string>} mapping
   */
  function technicalReplicateGroups(rows, mapping) {
    if (!mapping.sample_id || !mapping.omic) return [];
    /** @type {Record<string,number>} */
    const counts = {};
    for (const row of rows) {
      const sample = row[mapping.sample_id] ?? '';
      const rawOmic = row[mapping.omic] ?? '';
      const omic = canonicalOmic(rawOmic) || normalise(rawOmic);
      if (!sample || !omic) continue;
      const key = `${sample}::${omic}`;
      counts[key] = (counts[key] ?? 0) + 1;
    }
    return Object.entries(counts).filter(([, count]) => count > 1).map(([key, count]) => ({ key, count }));
  }

  /** @param {'transcriptomics' | 'proteomics' | 'metabolomics'} layer */
  function expectedAssays(layer) {
    if (!columnMapping.assay_id || !columnMapping.omic) return new Set();
    return new Set(
      metadataRows
        .filter((row) => canonicalOmic(mappedValue(row, 'omic')) === layer)
        .map((row) => mappedValue(row, 'assay_id'))
        .filter(Boolean)
    );
  }

  /**
   * @param {'transcriptomics' | 'proteomics' | 'metabolomics'} layer
   */
  function matrixMatch(layer) {
    const expected = expectedAssays(layer);
    const columns = new Set(matrixInfo[layer].sampleIds);
    const rows = new Set(matrixInfo[layer].rowIds);

    if (!expected.size || !columns.size && !rows.size) {
      return {
        matched: 0,
        expected: expected.size,
        observed: columns.size,
        missing: [],
        extra: [],
        orientation: 'unknown'
      };
    }
    const columnMatched = [...expected].filter((id) => columns.has(id));
    const rowMatched = [...expected].filter((id) => rows.has(id));
    const transposed = rowMatched.length > columnMatched.length;
    const observed = transposed ? rows : columns;
    const matched = transposed ? rowMatched : columnMatched;
    return {
      matched: matched.length,
      expected: expected.size,
      observed: observed.size,
      missing: [...expected].filter((id) => !observed.has(id)),
      extra: [...observed].filter((id) => !expected.has(id)),
      orientation: transposed ? 'transposed' : 'features-by-assays'
    };
  }

  /**
   * @param {string} key
   * @param {string} value
   */
  function setMapping(key, value) {
    columnMapping = { ...columnMapping, [key]: value };
  }

  /** @param {string} column */
  function toggleCovariate(column) {
    selectedCovariates = selectedCovariates.includes(column)
      ? selectedCovariates.filter((item) => item !== column)
      : [...selectedCovariates, column];
  }

  /**
   * @param {'transcriptomics' | 'proteomics' | 'metabolomics'} layer
   * @param {{rowIds:string[]}} info
   */
  function inferFeatureIdType(layer, info) {
    const ids = info.rowIds.slice(0, 50).map((id) => id.trim()).filter(Boolean);
    if (!ids.length) return { type: 'unknown', confidence: 0 };
    /** @type {Record<string, (id:string) => boolean>} */
    const tests = {
      ensembl_gene: (id) => (/^ENSG\d+(?:\.\d+)?$/i).test(id),
      ensembl_protein: (id) => (/^ENSP\d+(?:\.\d+)?$/i).test(id),
      uniprot: (id) => (/^(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])(?:-\d+)?$/).test(id),
      hmdb: (id) => (/^HMDB\d+$/i).test(id),
      chebi: (id) => (/^CHEBI:\d+$/i).test(id),
      inchikey: (id) => (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/i).test(id),
      pubchem: (id) => (/^CID:?\d+$/i).test(id),
      kegg_compound: (id) => (/^C\d{5}$/i).test(id),
      entrez: (id) => (/^\d+$/).test(id)
    };
    let best = { type: 'unknown', confidence: 0 };
    for (const [type, test] of Object.entries(tests)) {
      const score = ids.filter(test).length / ids.length;
      if (score > best.confidence) best = { type, confidence: score };
    }
    return best.confidence >= 0.6 ? best : { type: 'unknown', confidence: best.confidence };
  }

  function downloadGeneratedTemplate() {
    const headers = ['subject_id', 'sample_id', 'assay_id', 'omic', 'condition'];
    if (longitudinal === 'yes' || Number(timepointCount) > 1) headers.push('timepoint');
    if (batchKnown !== 'no') headers.push('batch');
    if (technicalReplicatesExpected !== 'no') headers.push('technical_replicate');
    if (outcomeType !== 'none' && outcomeType !== 'survival') headers.push('outcome');
    if (outcomeType === 'survival') headers.push('survival_time', 'survival_event');
    if (['untargeted_lcms','targeted_lcms','gcms'].includes(metabolomicsPlatform)) headers.push('sample_type', 'injection_order');
    if (covariatesAvailable === 'yes') headers.push('covariate_1');

    const omics = ['transcriptomics', 'proteomics', 'metabolomics'];
    const nTime = longitudinal === 'yes' ? Math.max(2, Number(timepointCount) || 2) : 1;
    const exampleRows = [];
    for (let subject = 1; subject <= 2; subject += 1) {
      for (let time = 0; time < nTime; time += 1) {
        for (const omic of omics) {
          const prefix = omic === 'transcriptomics' ? 'RNA' : omic === 'proteomics' ? 'PROT' : 'MET';
          /** @type {Record<string,string>} */
          const values = {
            subject_id: `SUBJ${String(subject).padStart(3, '0')}`,
            sample_id: `SUBJ${String(subject).padStart(3, '0')}_T${time}`,
            assay_id: `${prefix}${String((subject - 1) * nTime + time + 1).padStart(3, '0')}`,
            omic,
            condition: subject === 1 ? 'control' : 'treatment',
            timepoint: `T${time}`,
            batch: `${prefix}_B1`,
            technical_replicate: '1',
            outcome: '',
            survival_time: '',
            survival_event: '',
            sample_type: 'biological',
            injection_order: String((subject - 1) * nTime * 3 + time * 3 + omics.indexOf(omic) + 1),
            covariate_1: ''
          };
          exampleRows.push(headers.map((header) => values[header] ?? '').join(','));
        }
      }
    }
    const csv = [headers.join(','), ...exampleRows].join('\n');
    const blob = new Blob([csv], { type: 'text/csv;charset=utf-8' });
    const href = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = href;
    anchor.download = 'multiomics_metadata_from_protocol.csv';
    anchor.click();
    URL.revokeObjectURL(href);
  }

  function backendBaseUrl() {
    return referenceBackendUrl.trim().replace(/\/$/, '');
  }

  /**
   * @param {string} url
   * @param {RequestInit} options
   * @param {number} timeoutMs
   */
  async function fetchWithTimeout(url, options = {}, timeoutMs = 1800) {
    const controller = new AbortController();
    const timer = setTimeout(() => controller.abort(), timeoutMs);
    try {
      return await fetch(url, { ...options, signal: controller.signal });
    } finally {
      clearTimeout(timer);
    }
  }

  async function checkReferenceBackend(showMessage = true) {
    if (referenceBackendMode === 'browser') {
      referenceBackendStatus = 'disabled';
      referenceBackendHealth = null;
      if (showMessage) referenceBackendMessage = t('Backend R désactivé.', 'R backend disabled.');
      return null;
    }
    const url = backendBaseUrl();
    if (!url) {
      referenceBackendStatus = 'unavailable';
      referenceBackendHealth = null;
      if (showMessage) referenceBackendMessage = t('URL du backend R manquante.', 'R backend URL is missing.');
      return null;
    }
    referenceBackendStatus = 'checking';
    try {
      const response = await fetchWithTimeout(url + '/health', { headers: { Accept: 'application/json' } }, 1800);
      if (!response.ok) throw new Error('HTTP ' + response.status);
      const health = await response.json();
      const preflight = evaluateReferenceBackendHealth(health, {
        intent: analysisIntent,
        objective,
        outcomeType,
        designType,
        hasTranscriptomics: Boolean(files.transcriptomics),
        hasProteomics: Boolean(files.proteomics),
        hasMetabolomics: Boolean(files.metabolomics),
        transcriptomicsValues
      });
      if (!preflight.ok) {
        referenceBackendStatus = 'error';
        referenceBackendHealth = null;
        referenceBackendMessage = t(preflight.message, preflight.messageEn);
        return null;
      }
      referenceBackendStatus = 'available';
      referenceBackendHealth = health;
      referenceBackendMessage = t(preflight.message, preflight.messageEn);
      return health;
    } catch (error) {
      referenceBackendStatus = 'unavailable';
      referenceBackendHealth = null;
      if (showMessage) {
        referenceBackendMessage = t(
          'Backend R indisponible ; le moteur navigateur reste utilisable.',
          'R backend unavailable; the browser engine remains usable.'
        );
      }
      return null;
    }
  }

  async function referenceBackendPayload() {
    if (!files.metadata) throw new Error('Metadata file is required.');
    /** @type {Record<string,string>} */
    const matrices = {};
    for (const layer of omicLayers) {
      if (files[layer]) matrices[layer] = await files[layer].text();
    }
    return {
      protocol: {
        organism,
        objective,
        longitudinal: longitudinal === 'yes',
        designType,
        studySetting,
        analysisIntent,
        rnaPrimaryContrast: analysisIntent === 'confirmatory' &&
          objective === 'groups' && files.transcriptomics &&
          transcriptomicsValues === 'raw_counts'
          ? ['condition', rnaseqComparison, rnaseqReference] : null,
        unitType,
        groupCount,
        timepointCount,
        sampleOverlap,
        technicalReplicatesExpected,
        batchKnown,
        outcomeType,
        outcomeTimepoint,
        covariatesAvailable,
        covariateColumns: selectedCovariates,
        partialOmicsExpected,
        msBlankFilter,
        msBlankFold: Number(msBlankFold),
        msQcRsdFilter,
        msQcRsdThreshold: Number(msQcRsdThreshold),
        msDriftCorrection,
        msMnarStrategy,
        validateDiabloHoldout
      },
      dataTypes: {
        transcriptomics: transcriptomicsValues,
        proteomics: proteomicsValues,
        metabolomics: metabolomicsValues
      },
      identifierTypes: {
        transcriptomics: transcriptomicsIdType,
        proteomics: proteomicsIdType,
        metabolomics: metabolomicsIdType
      },
      columnMapping,
      metadataCsv: await files.metadata.text(),
      matrices
    };
  }

  async function generateSingleOmicExploratorySheet() {
    if (!confirmIndependentAssays || objective !== 'explore' || omicsCount !== 1 || files.metadata)
      return;
    const layer = omicLayers.find((key) => Boolean(files[key]));
    if (!layer) return;
    // File selection is asynchronous: a fast researcher can click here while
    // file.text() is still being parsed. Never silently lose that click or
    // construct metadata from a stale/empty sample-column list.
    const selectedSource = files[layer];
    if (!matrixInfo[layer].sampleIds.length || matrixInfo[layer].error)
      await inspectMatrix(layer, selectedSource);
    if (files[layer] !== selectedSource || files.metadata) return;
    const assays = matrixInfo[layer].sampleIds;
    if (matrixInfo[layer].error || !assays.length ||
        assays.some((id) => !id) || new Set(assays).size !== assays.length) {
      analysisError = t('Colonnes de mesures absentes ou dupliquées. Corrigez le fichier avant la création du tableau.', 'Missing or duplicate assay IDs. Correct the matrix before creating a sample sheet.');
      return;
    }

    const header = [
      'subject_id',
      'sample_id',
      'assay_id',
      'omic',
      'condition',
      'timepoint',
      'batch',
      'technical_replicate'
    ];

    /** @param {unknown} value */
    const quote = (value) => '"' + String(value ?? '').replaceAll('"', '""') + '"';

    const text = [
      header.map(quote).join(','),
      ...assays.map((assay) => [assay, assay, assay, layer, '', 'T0', '', '1'].map(quote).join(','))
    ].join('\n') + '\n';

    const file = new File([text], 'exploratory_sample_sheet_generated.csv', { type: 'text/csv' });

    files = { ...files, metadata: file };
    msAutoMetadata = false;
    exploratorySheetGenerated = true;
    designType = 'independent';
    longitudinal = 'no';
    paired = 'no';
    timepointCount = '1';
    groupCount = '1';
    outcomeType = 'none';
    await inspectMetadata(file);
  }

  /** @param {File | null} file */
  async function sha256OfFile(file) {
    if (!file || !globalThis.crypto?.subtle) return null;
    try {
      const digest = await globalThis.crypto.subtle.digest('SHA-256', await file.arrayBuffer());
      return Array.from(new Uint8Array(digest), (byte) => byte.toString(16).padStart(2, '0')).join('');
    } catch {
      return null; // Integrity checks must not turn an otherwise valid analysis into a failed run.
    }
  }

  /** @param {any} result */
  async function attachMsProvenance(result) {
    if (!msOriginalFile && !msUserAnnotationFile) return;
    const source = msOriginalFile ? {
      name: msOriginalFile.name, bytes: msOriginalFile.size, sha256: await sha256OfFile(msOriginalFile),
      importFormat: msImportFormat || 'matrix', annotationColumns: msAnnotationFields,
      extractedAnnotationRows: msAnnotationCsv ? Math.max(0, msAnnotationCsv.trim().split(/\r?\n/).length - 1) : 0,
      quantitativeMatrixSha256: result?.inputIntegrity?.files?.metabolomics?.sha256 || null,
      note: 'Original vendor file must be archived alongside the derived matrix. A peak is not a verified chemical identity.'
    } : null;
    const crosswalk = msUserAnnotationFile ? {
      name: msUserAnnotationFile.name, bytes: msUserAnnotationFile.size,
      sha256: await sha256OfFile(msUserAnnotationFile),
      status: 'user_declared_unverified', affects: 'visualizations_only'
    } : null;
    result.msImportProvenance = { originalSource: source, userIdentifierCrosswalk: crosswalk };

    result.reproducibility = {
      ...result.reproducibility || {},
      msSource: result.msImportProvenance
    };
  }

  function downloadMsAnnotationTable() {
    if (!msAnnotationCsv) return;
    const link = document.createElement('a');
    const url = URL.createObjectURL(new Blob([msAnnotationCsv], { type:'text/csv;charset=utf-8' }));
    link.href = url;
    link.download = 'ms_technical_annotations_preserved.csv';
    link.click();
    URL.revokeObjectURL(url);
  }

  async function runReferenceBackend() {
    const health = await checkReferenceBackend(false);
    if (!health) {
      if (referenceBackendMode === 'required') {
        throw new Error(referenceBackendMessage || 'Reference R backend is required but unavailable at ' + backendBaseUrl() + '.');
      }
      return {
        status: 'unavailable',
        url: backendBaseUrl(),
        message: 'Reference R backend was not reachable; browser results were retained.'
      };
    }
    const activeOmicCount = Number(Boolean(files.transcriptomics)) + Number(Boolean(files.proteomics)) + Number(Boolean(files.metabolomics));
    const compatibility = mofa2BackendCompatibility(health,{objective,omicCount:activeOmicCount});
    if(!compatibility.compatible) {
      referenceBackendStatus = 'error';
      referenceBackendMessage = t(compatibility.message || 'Moteur R MOFA2 incompatible.',compatibility.messageEn || 'Incompatible MOFA2 R backend.');
      if(referenceBackendMode === 'required')throw new Error(referenceBackendMessage);
      return {
        status:'blocked',message:referenceBackendMessage,url:backendBaseUrl(),
        methods:{mofa2:{status:'blocked',message:referenceBackendMessage}},
        referenceMethodExecuted:false,
        compatibilityStatus:compatibility.reason
      };
    }
    referenceBackendStatus = 'running';
    referenceBackendMessage = t('Méthodes R de référence en cours…', 'Reference R methods running…');
    const payload = await referenceBackendPayload();
    const payloadJson = JSON.stringify(payload);
    const response = await fetch(backendBaseUrl() + '/run', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
      body: payloadJson
    });
    const body = await response.json().catch(() => ({}));
    if (!response.ok || body.status === 'error') {
      referenceBackendStatus = 'error';
      referenceBackendMessage = body.message || 'HTTP ' + response.status;
      if (referenceBackendMode === 'required') throw new Error(referenceBackendMessage);
      return {
        status: 'error',
        url: backendBaseUrl(),
        message: referenceBackendMessage
      };
    }
    const checkedBody = validateMofa2MethodResult(body);
    referenceBackendStatus = checkedBody.mofa2OrientationBlocked ? 'error' : 'done';
    const completed = Object.values(checkedBody.methods || {}).some((method) => method?.status === 'ok');
    // "Require R" must mean that a REAL subject-blocked reference model
    // completed, not merely that /health was reachable or some unrelated
    // method returned HTTP 200. Fail closed if pairing was not applied.
    if(referenceBackendMode === 'required' && objective === 'groups' &&
       designType === 'paired' && files.transcriptomics){
      const countMethods = [
        checkedBody.methods?.transcriptomics_differential_deseq2,
        checkedBody.methods?.transcriptomics_differential_voom
      ];

      const validPairedFit = countMethods.some((method) => method?.status === 'ok' && method.pairedDesignContract?.pairing === 'same_subject_two_conditions' && method.pairedDesignContract?.biological_pairs === pairedAudit?.nPairs && method.pairedDesignContract?.certified === false);

      if (!validPairedFit) {
        const diagnostics = countMethods.filter(Boolean).map((method) => method.method + ': ' + method.status + (method.message ? ' (' + method.message + ')' : '')).join('; ');

        throw new Error(t('Aucun modèle RNA-seq apparié valide n’a été exécuté sous R. Le résultat exploratoire du navigateur ne remplace pas une inférence appariée. Vérifiez les paires, le contraste et les messages R. ', 'No valid paired RNA-seq model was fitted in R. Browser exploration cannot substitute for paired inference. Check your pairs, contrast and R diagnostics. ') + diagnostics);
      }
    }
    referenceBackendMessage = checkedBody.mofa2OrientationBlocked
      ? t('Résultat MOFA2 masqué : ce backend n’a pas prouvé l’orientation correcte des variables et des sujets. Mettez à jour R.', 'MOFA2 result hidden: this backend did not prove the correct feature and subject orientation. Update R.')
      : completed
        ? t('Certaines méthodes R ont été exécutées : voir leur statut ci-dessous.', 'Some R methods ran: see their status below.')
        : t('Moteur R connecté, mais aucune méthode statistique de référence exécutée pour ce plan.', 'R backend connected, but no reference statistical method ran for this design.');
    if (completed && !checkedBody.mofa2OrientationBlocked) {
      referenceReplay = {requestJson:payloadJson,responseJson:JSON.stringify(body)};
    }
    return { ...checkedBody,url:backendBaseUrl(),referenceMethodExecuted:completed };
  }

  async function runAnalysis() {
    analysisError = '';
    analysisResult = null;
    referenceReplay = null;
    if (msImportState === 'loading' || msImportState === 'parsed') {
      analysisError = t('Terminez l’import MS avant de lancer l’analyse.', 'Complete the MS import before running an analysis.');
      return;
    }
    if(objective === 'groups' && designType === 'paired' && files.transcriptomics){
      if(!pairedAudit?.ok){
        analysisError = pairedProblemMessage(pairedAudit?.code);
        return;
      }
      if(referenceBackendMode !== 'required'){
        analysisError = t('L’analyse RNA-seq appariée nécessite « R requis » : le navigateur ne doit pas traiter deux prélèvements d’un même sujet comme indépendants.', 'Paired RNA-seq requires Require R: the browser must not treat two samples from one subject as independent.');
        return;
      }
    }

    if (analysisIntent === 'confirmatory' && referenceBackendMode !== 'required') {
      analysisError = t('Pour une analyse confirmatoire, choisissez « R requis » et vérifiez sa connexion. Aucune conclusion confirmatoire automatique ne sera produite.', 'For a confirmatory workflow, choose Require R and verify its connection. No automatic confirmatory approval will be given.');
      return;
    }
    // Confirmatory display is allowed only for objectives with a complete
    // native-R projection contract. Group, longitudinal group-by-time and
    // survival are supported; other outcome families remain exploratory.
    const confirmatoryObjectiveSupported = objective === 'groups' ||
      objective === 'time' ||
      (objective === 'outcome' && outcomeType === 'survival');
    if (analysisIntent === 'confirmatory' && !confirmatoryObjectiveSupported) {
      const gate = assessConfirmatoryResult({
        protocol:{analysisIntent,objective,outcomeType}
      });
      analysisError = t(gate.message, gate.messageEn);
      return;
    }
    if (!ready) {
      analysisError = t('Fournissez un tableau des échantillons et au moins une matrice RNA, protéines ou métabolites, puis vérifiez les colonnes et l’objectif.', 'Upload a sample sheet and at least one RNA, protein or metabolite matrix, then verify the columns and objective.');
      return;
    }
    if (analysisIntent === 'confirmatory' && objective === 'groups' &&
        files.transcriptomics && transcriptomicsValues === 'raw_counts' &&
        (rnaseqAvailableConditions.length !== 2 ||
         !rnaseqAvailableConditions.includes(rnaseqReference) ||
         !rnaseqAvailableConditions.includes(rnaseqComparison) ||
         rnaseqReference === rnaseqComparison)) {
      analysisError = t(
        'RNA-seq : choisissez explicitement le groupe de référence et le groupe comparé parmi les deux groupes présents dans les métadonnées. Leur ordre change le signe du log2FC.',
        'RNA-seq: explicitly select different reference and comparison groups from the two groups in the sample sheet. Reversing them changes the sign of log2FC.');
      return;
    }
    analysisStatus = 'running';
    try {
      analysisResult = await runDeterministicAnalysis({
        files,
        metadataRows,
        columnMapping,
        protocol: {
          organism,
          objective,
          longitudinal: longitudinal === 'yes',
          designType,
          studySetting,
          analysisIntent,
          rnaPrimaryContrast: analysisIntent === 'confirmatory' &&
            objective === 'groups' && files.transcriptomics &&
            transcriptomicsValues === 'raw_counts'
            ? ['condition', rnaseqComparison, rnaseqReference] : null,
          unitType,
          groupCount,
          timepointCount,
          sampleOverlap,
          technicalReplicatesExpected,
          batchKnown,
          outcomeType,
          outcomeTimepoint,
          covariatesAvailable,
          covariateColumns: selectedCovariates,
          partialOmicsExpected,
          msBlankFilter,
          msBlankFold: Number(msBlankFold),
          msQcRsdFilter,
          msQcRsdThreshold: Number(msQcRsdThreshold),
          msDriftCorrection,
          msMnarStrategy,
          validateDiabloHoldout
        },
        dataTypes: {
          transcriptomics: transcriptomicsValues,
          proteomics: proteomicsValues,
          metabolomics: metabolomicsValues
        },
        identifierTypes: {
          transcriptomics: transcriptomicsIdType,
          proteomics: proteomicsIdType,
          metabolomics: metabolomicsIdType
        },
        resolveIdentifiers,
        useReactome
      });
      if (referenceBackendMode !== 'browser') {
        const referenceBackend = await runReferenceBackend();
        analysisResult = projectNativeRNumbers({ ...analysisResult, referenceBackend });
        if (analysisIntent === 'confirmatory') {
          const gate = assessConfirmatoryResult(analysisResult);
          if (!gate.ok) throw new Error(t(gate.message, gate.messageEn));
        }
      }

      try {
        analysisResult.visualizations = await buildMultiomicsVisualizationData({
          files,
          metadataRows,
          columnMapping,
          dataTypes: {
            transcriptomics: transcriptomicsValues,
            proteomics: proteomicsValues,
            metabolomics: metabolomicsValues
          },
          analysisResult,
          metaboliteAnnotations
        });
      } catch (visualizationError) {
        analysisResult.visualizationError = visualizationError instanceof Error
          ? visualizationError.message
          : 'Visualization preparation failed.';
      }
      await attachMsProvenance(analysisResult);
      analysisResult.scientificAssurance = assessScientificAssurance(analysisResult, { demo: demoLoaded });
      analysisStatus = 'done';
    } catch (error) {
      // Never display an incomplete result when a later essential step fails.
      analysisResult = null;
      analysisStatus = 'error';
      analysisError = error instanceof Error ? error.message : 'Analysis failed.';
    }
  }


  function currentPortableSettings() {
    return {
            studyName, objective, organism, unitType, studySetting, designType,
            groupCount, timepointCount, sampleOverlap, technicalReplicatesExpected,
            batchKnown, outcomeType, outcomeTimepoint, covariatesAvailable,
            partialOmicsExpected, selectedCovariates, groupVariable, outcome,
            subjectCount, paired, longitudinal, transcriptomicsPlatform,
            transcriptomicsValues, transcriptomicsIdType, proteomicsPlatform,
            proteomicsValues, proteomicsIdType, metabolomicsPlatform,
            metabolomicsValues, metabolomicsIdType, msBlankFilter, msBlankFold,
            msQcRsdFilter, msQcRsdThreshold, msDriftCorrection, msMnarStrategy,
            analysisIntent, rnaseqReference, rnaseqComparison, validateDiabloHoldout, confirmIndependentAssays, msAutoMetadata,
            columnMapping, metaboliteAnnotations, msAnnotationCsv,
            msAnnotationFields, msImportFormat, annotationFileName,
            demoLoaded
          };
  }


  // One self-contained researcher dossier: full spreadsheets + methods +
  // provenance/hash manifest + source study for manual reproducible rerun.
  /** @param {Event} event */
  async function importIndividualSignatures(event) {
    individualSignatureRows = [];
    individualSignatureFilename = '';
    individualSignatureError = '';
    const input = event.currentTarget;
    const file = input instanceof HTMLInputElement ? input.files?.[0] : null;
    if (!file) return;
    try {
      individualSignatureRows = parseIndividualSignatureCsv(await file.text());
      individualSignatureFilename = file.name;
    } catch (error) {
      individualSignatureError = error instanceof Error ? error.message : 'Invalid score CSV.';
    }
  }

  async function downloadResearchBundle() {
    if (!analysisResult || portableProjectBusy) return;
    portableProjectBusy = true;
    researchBundleError = '';
    researchBundleMessage = '';
    try {
      if(individualSignatureError)
        throw Error('Individual signature import is invalid: '+individualSignatureError);
      const sourceProject = await createPortableProject({
        files,
        extraFiles:{ msOriginal:msOriginalFile, msUserAnnotations:msUserAnnotationFile },
        settings:currentPortableSettings(),
        result:null // historic outputs are never imported as new evidence
      });
      const dossier = await buildResearchBundle({
        result:individualSignatureRows.length
          ? {...analysisResult,individualSignatures:individualSignatureRows}
          : analysisResult,
        portableProject:sourceProject,
        referenceReplay,
        // Exact vector snapshots of CURRENTLY RENDERED browser figures only.
        figureSnapshots:[...document.querySelectorAll('#analysis-results svg[role="img"]')]
          .slice(0,12).map((svg,index)=>({name:'figure_'+String(index+1).padStart(2,'0')+'.svg',
            data:new XMLSerializer().serializeToString(svg)})),
        locale:$language === 'en' ? 'en' : 'fr'
      });
      const url=URL.createObjectURL(dossier.blob);
      const a=document.createElement('a');
      a.href=url;
      a.download=dossier.filename;
      a.click();
      URL.revokeObjectURL(url);
      researchBundleMessage = t(
        'Dossier ZIP préparé : tableaux complets, méthodes, sources et SHA-256. Il contient potentiellement des données confidentielles. Vérifiez les limites et sauvegardez-le de façon sécurisée.',
        'Research ZIP prepared: full tables, methods, source files and SHA-256. This may contain sensitive research data. Review the limitations and store securely.'
      );
    } catch (error) {
      researchBundleError = error instanceof Error ? error.message : 'Could not create the research dossier.';
    } finally {
      portableProjectBusy = false;
    }
  }

  // Export original study text and protocol alongside the current results.
  // This is a sensitive local file, never uploaded to a third-party service.
  async function exportPortableProject() {
    if (!analysisResult || portableProjectBusy) return;
    portableProjectBusy = true;
    portableProjectError = '';
    portableProjectMessage = '';
    try {
      const settings = currentPortableSettings();
      const packageData = await createPortableProject({
        files,
        extraFiles:{ msOriginal:msOriginalFile, msUserAnnotations:msUserAnnotationFile },
        settings,
        result:analysisResult
      });
      const blob = new Blob([JSON.stringify(packageData)],{type:'application/json'});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href=url;
      a.download='multiomics_project_reproducible_v1.json';
      a.click();
      URL.revokeObjectURL(url);
      portableProjectMessage = t(
        'Projet téléchargé avec les fichiers source, empreintes SHA-256, protocole et résultats. Contient potentiellement des données confidentielles.',
        'Project downloaded with original inputs, SHA-256 digests, protocol and results. May contain confidential research data.'
      );
    } catch (error) {
      portableProjectError = error instanceof Error ? error.message : 'Unable to export project.';
    } finally {
      portableProjectBusy = false;
    }
  }

  /** @param {Event} event */
  async function importPortableProject(event) {
    const input = /** @type {HTMLInputElement} */ (event.currentTarget);
    const source = input.files?.[0];
    if (!source || portableProjectBusy) return;
    portableProjectBusy = true;
    portableProjectMessage = '';
    portableProjectError = '';
    try {
      if (source.size > PORTABLE_PROJECT_MAX_BYTES * 2)
        throw new Error('Portable project is too large for the browser importer.');
      const imported = await parsePortableProject(await source.text());
      const next = imported.settings;
      // This explicitly restores only data and settings, never historical
      // results as fresh evidence. No external services run during import.
      files = {
        metadata:imported.files.metadata,
        transcriptomics:imported.files.transcriptomics||null,
        proteomics:imported.files.proteomics||null,
        metabolomics:imported.files.metabolomics||null
      };
      await inspectMetadata(files.metadata);
      for (const layer of omicLayers) await inspectMatrix(layer,files[layer]);
      const savedMapping = next.columnMapping || {};
      if (typeof savedMapping === 'object' && !Array.isArray(savedMapping)) {
        columnMapping = Object.fromEntries(Object.entries(savedMapping)
          .filter(([key,value])=>typeof value==='string' &&
            (!value || metadataHeaders.includes(value))));
      }
      studyName = next.studyName ?? '';
      objective = next.objective ?? 'explore';
      organism = next.organism ?? 'human';
      unitType = next.unitType ?? 'participant';
      studySetting = next.studySetting ?? 'clinical_observational';
      designType = next.designType ?? 'independent';
      groupCount = next.groupCount ?? '2';
      timepointCount = next.timepointCount ?? '1';
      sampleOverlap = next.sampleOverlap ?? 'same_specimen';
      technicalReplicatesExpected = next.technicalReplicatesExpected ?? 'unknown';
      batchKnown = next.batchKnown ?? 'unknown';
      outcomeType = next.outcomeType ?? 'none';
      outcomeTimepoint = next.outcomeTimepoint ?? '';
      covariatesAvailable = next.covariatesAvailable ?? 'yes';
      partialOmicsExpected = next.partialOmicsExpected ?? 'no';
      selectedCovariates = Array.isArray(next.selectedCovariates)
        ? /** @type {string[]} */ (next.selectedCovariates.filter((c) => typeof c === 'string' && metadataHeaders.includes(c)))
        : [];

      groupVariable = next.groupVariable ?? '';
      outcome = next.outcome ?? '';
      subjectCount = next.subjectCount;
      paired = next.paired ?? 'no';
      longitudinal = next.longitudinal ?? 'no';
      transcriptomicsPlatform = next.transcriptomicsPlatform ?? 'bulk_rnaseq';
      transcriptomicsValues = next.transcriptomicsValues ?? 'raw_counts';
      transcriptomicsIdType = next.transcriptomicsIdType ?? 'ensembl_gene';
      proteomicsPlatform = next.proteomicsPlatform ?? 'label_free';
      proteomicsValues = next.proteomicsValues ?? 'lfq_intensity';
      proteomicsIdType = next.proteomicsIdType ?? 'uniprot';
      metabolomicsPlatform = next.metabolomicsPlatform ?? 'untargeted_lcms';
      metabolomicsValues = next.metabolomicsValues ?? 'peak_area';
      metabolomicsIdType = next.metabolomicsIdType ?? 'chebi';
      msBlankFilter = next.msBlankFilter ?? 'flag';
      msBlankFold = next.msBlankFold ?? 5;
      msQcRsdFilter = next.msQcRsdFilter ?? 'yes';
      msQcRsdThreshold = next.msQcRsdThreshold ?? 0.30;
      msDriftCorrection = next.msDriftCorrection ?? 'yes';
      msMnarStrategy = next.msMnarStrategy ?? 'none';
      analysisIntent = next.analysisIntent ?? 'exploratory';
      rnaseqReference = next.rnaseqReference ?? '';
      rnaseqComparison = next.rnaseqComparison ?? '';
      validateDiabloHoldout = Boolean(next.validateDiabloHoldout);
      confirmIndependentAssays = Boolean(next.confirmIndependentAssays);
      msAutoMetadata = Boolean(next.msAutoMetadata);
      metaboliteAnnotations = Array.isArray(next.metaboliteAnnotations)?next.metaboliteAnnotations:[];
      msAnnotationCsv = next.msAnnotationCsv ?? '';
      msAnnotationFields = Array.isArray(next.msAnnotationFields)?next.msAnnotationFields:[];
      msImportFormat = next.msImportFormat ?? '';
      annotationFileName = next.annotationFileName ?? '';
      msOriginalFile = imported.files.msOriginal||null;
      msUserAnnotationFile = imported.files.msUserAnnotations||null;
      demoLoaded = Boolean(next.demoLoaded);
      // Import NEVER silently submits cohort data to R, Reactome or ChEBI.
      referenceBackendMode = 'browser';
      referenceBackendStatus = 'unchecked';
      referenceBackendMessage = '';
      referenceBackendHealth = null;
      useReactome = false;
      resolveIdentifiers = false;
      analysisStatus = 'idle';
      analysisError = '';
      analysisResult = null;
      portableProjectMessage = t(
        'Projet importé et intégrité SHA-256 vérifiée. Vérifiez le protocole puis relancez l’analyse ; aucune donnée n’a été envoyée à un serveur. Les résultats archivés ne sont pas présentés comme une nouvelle analyse.',
        'Project imported; SHA-256 verified. Review the protocol and rerun; no data was sent to a server. Archived results are not shown as freshly computed.'
      );
    } catch (error) {
      portableProjectError = error instanceof Error ? error.message : 'Invalid portable project.';
    } finally {
      portableProjectBusy = false;
      input.value = '';
    }
  }

  function downloadResearchMethodsAudit() {
    if(!analysisResult)return;
    const audit=buildResearchMethodsManifest(analysisResult);
    const blob=new Blob([JSON.stringify(audit,null,2)],{type:'application/json;charset=utf-8'});
    const href=URL.createObjectURL(blob);
    const link=document.createElement('a');
    link.href=href;
    link.download='multiomics_methods_scientific_audit.json';
    link.click();
    URL.revokeObjectURL(href);
  }

  function downloadAnalysisJson() {
    if (!analysisResult) return;
    const blob = new Blob([JSON.stringify(analysisResult, null, 2)], { type: 'application/json' });
    const href = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = href;
    anchor.download = 'multiomics_deterministic_results.json';
    anchor.click();
    URL.revokeObjectURL(href);
  }

  /** @param {any} qc */
  function poorReplicateCount(qc) {
    const correlations = /** @type {any[]} */ (qc?.replicateCorrelations || []);
    return correlations.filter((item) => Boolean(item?.warning)).length;
  }

  /** @param {unknown} value */
  function escapeHtml(value) {
    return String(value ?? '')
      .replaceAll('&', '&amp;')
      .replaceAll('<', '&lt;')
      .replaceAll('>', '&gt;')
      .replaceAll('"', '&quot;')
      .replaceAll("'", '&#039;');
  }

  function downloadReproducibleReport() {
    if (!analysisResult) return;
    const result = analysisResult;
    const layerSections = Object.entries(result.layers || {}).map(([layer, layerResult]) => {
      const layerRows = /** @type {any[]} */ (layerResult.rows || []);
      const rows = layerRows.slice(0, 50).map((row) =>
        '<tr><td>' + escapeHtml(row.feature) + '</td><td>' +
        escapeHtml(Number.isFinite(row.effect) ? Number(row.effect).toPrecision(4) : '') + '</td><td>' +
        escapeHtml(row.pValue == null ? '' : Number(row.pValue).toPrecision(4)) + '</td><td>' +
        escapeHtml(row.qValue == null ? '' : Number(row.qValue).toPrecision(4)) + '</td></tr>'
      ).join('');
      const qc = layerResult.qc || {};
      return '<section><h2>' + escapeHtml(omicLabel(layer)) + '</h2>' +
        '<p>Preview only: first ' + Math.min(layerRows.length, 50) + ' of ' + layerRows.length + ' features. Complete CSV tables are available in the research ZIP dossier.</p>' +
        '<p><strong>Model:</strong> ' + escapeHtml(layerResult.inferenceMethod || layerResult.mode || '') + '</p>' +
        '<p><strong>QC:</strong> ' + escapeHtml(qc.featuresAfter ?? '—') + '/' + escapeHtml(qc.featuresBefore ?? '—') +
        ' features retained; median missing ' + escapeHtml(Number.isFinite(qc.medianMissingFraction) ? (100 * qc.medianMissingFraction).toFixed(1) + '%' : '—') + '.</p>' +
        (qc.inferenceTier?.note ? '<p><strong>Inference tier:</strong> ' + escapeHtml(qc.inferenceTier.label || qc.inferenceTier.level) + ' — ' + escapeHtml(qc.inferenceTier.note) + '</p>' : '') +
        '<table><thead><tr><th>Feature</th><th>Effect</th><th>p</th><th>q BH</th></tr></thead><tbody>' + rows + '</tbody></table></section>';
    }).join('');

    const pathwayRows = /** @type {any[]} */ (result.reactome?.consensus || []);
    const pathways = pathwayRows.slice(0, 30).map((pathway) =>
      '<tr><td>' + escapeHtml(pathway.name) + '</td><td>' +
      escapeHtml(Number.isFinite(pathway.assayUniverseFdr) ? pathway.assayUniverseFdr.toPrecision(4) : 'not estimable') +
      '</td><td>' + escapeHtml(pathway.supportingLayers) + '</td></tr>'
    ).join('');

    const embeddedJson = JSON.stringify(result).replaceAll('<', '\\u003c');
    const html = '<!doctype html><html lang="' + ($language === 'en' ? 'en' : 'fr') + '"><head><meta charset="utf-8">' +
      '<meta name="viewport" content="width=device-width,initial-scale=1"><title>Multi-omics reproducible report</title>' +
      '<style>body{font:14px/1.55 system-ui,sans-serif;max-width:1100px;margin:40px auto;padding:0 24px;color:#161616}h1,h2{line-height:1.15}table{border-collapse:collapse;width:100%;margin:12px 0 28px}th,td{border:1px solid #ddd;padding:6px 8px;text-align:left}th{background:#f4f4f4}code,pre{font-family:ui-monospace,monospace}pre{white-space:pre-wrap;background:#f6f6f6;padding:12px}section{margin:32px 0}.muted{color:#666}</style></head><body>' +
      '<h1>PMx Explain — multi-omics reproducible report</h1>' +
      '<p class="muted">Generated ' + escapeHtml(result.generatedAt) + ' · engine ' + escapeHtml(result.engine?.version || 'unknown') + '</p>' +
      '<section><h2>Provenance</h2><pre>' + escapeHtml(JSON.stringify(result.inputManifest, null, 2)) + '</pre></section>' +
      '<section><h2>Protocol</h2><pre>' + escapeHtml(JSON.stringify(result.protocol, null, 2)) + '</pre></section>' +
      '<section><h2>Scientific evidence status (not an approval)</h2><pre>' + escapeHtml(JSON.stringify(result.scientificAssurance || {}, null, 2)) + '</pre></section>' +
      '<section><h2>Reproducibility methods audit (not a certificate)</h2><pre>' + escapeHtml(JSON.stringify(buildResearchMethodsManifest(result), null, 2)) + '</pre></section>' +
      '<section><h2>Sample structure</h2><pre>' + escapeHtml(JSON.stringify(result.metadataSummary, null, 2)) + '</pre></section>' +
      layerSections +
      (result.supervisedIntegration ? '<section><h2>Supervised multiblock integration</h2><pre>' + escapeHtml(JSON.stringify(result.supervisedIntegration, null, 2)) + '</pre></section>' : '') +
      (result.predictiveOutcome ? '<section><h2>Predictive validation</h2><pre>' + escapeHtml(JSON.stringify(result.predictiveOutcome, null, 2)) + '</pre></section>' : '') +
      (result.referenceBackend ? '<section><h2>Reference R backend</h2><pre>' + escapeHtml(JSON.stringify(result.referenceBackend, null, 2)) + '</pre></section>' : '') +
      '<section><h2>Reactome pathways</h2><table><thead><tr><th>Pathway</th><th>Assay-universe FDR</th><th>Supporting layers</th></tr></thead><tbody>' + pathways + '</tbody></table></section>' +
      '<section><h2>Interpretation limits</h2><p>This report separates observed data, statistical inference and external pathway knowledge. Associations are not causal claims. External validation is required for predictive use.</p></section>' +
      '<script type="application/json" id="multiomics-analysis-json">' + embeddedJson + '</scr' + 'ipt>' +
      '</body></html>';

    const blob = new Blob([html], { type: 'text/html;charset=utf-8' });
    const href = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href = href;
    link.download = 'multiomics_reproducible_report.html';
    link.click();
    URL.revokeObjectURL(href);
  }

  /** @param {'transcriptomics'|'proteomics'|'metabolomics'} layer */
  function downloadLayerCsv(layer) {
    const rows = analysisResult?.layers?.[layer]?.rows;
    if (!rows?.length) return;
    const methodKey = analysisResult.layers[layer]?.primaryReferenceMethodKey;
    const method = methodKey ? analysisResult.referenceBackend?.methods?.[methodKey] : null;
    const nativeCsv = method?.status === 'ok' && method?.referenceTableComplete === true &&
      typeof method?.referenceTableCsv === 'string' ? method.referenceTableCsv : null;
    const blob = new Blob([nativeCsv ?? fullLayerCsv(rows)], { type: 'text/csv;charset=utf-8' });
    const href = URL.createObjectURL(blob);
    const anchor = document.createElement('a');
    anchor.href = href;
    anchor.download = nativeCsv ? `${layer}_reference_R_complete.csv` : `${layer}_browser_exploratory.csv`;
    anchor.click();
    URL.revokeObjectURL(href);
  }

  /**
   * @param {number} value
   * @param {any[]} points
   * @param {'pc1'|'pc2'} key
   */
  function qcAxisPercent(value, points, key) {
    const values = points.map((point) => Number(point[key])).filter(Number.isFinite);
    if (!values.length) return 50;
    const min = Math.min(...values);
    const max = Math.max(...values);
    if (max === min) return 50;
    return 6 + 88 * (value - min) / (max - min);
  }

  /** @param {any} result @param {string} lang */
  function buildInterpretation(result, lang) {
    /** @param {string} fr @param {string} en */
    const pick = (fr, en) => lang === 'en' ? en : fr;
    const items = [];
    if (result.protocol?.objective === 'explore' && result.exploration?.components?.length) {
      const pc1 = result.exploration.components[0];
      const pct = Number.isFinite(pc1.explainedFraction) ? (100 * pc1.explainedFraction).toFixed(1) : '—';
      items.push({
        title: pick('Axe latent principal', 'Main latent axis'),
        text: pick(
          result.engine?.analysisMode === 'single_omic'
            ? 'PC1 résume ' + pct + ' % de la variance de cette omique chez ' + result.exploration.subjects + ' sujets. Les loadings décrivent les variables dominantes, sans prouver un effet de groupe.'
            : 'PC1 résume ' + pct + ' % de la variance multi-blocs pondérée chez ' + result.exploration.subjects + ' sujets communs.',
          result.engine?.analysisMode === 'single_omic'
            ? 'PC1 summarizes ' + pct + '% of this omics layer variance across ' + result.exploration.subjects + ' subjects. Loadings are descriptive and do not prove a group effect.'
            : 'PC1 summarizes ' + pct + '% of the balanced multi-block variance across ' + result.exploration.subjects + ' shared subjects.'
        )
      });
      items.push({
        title: pick('Attention au signe', 'Sign is arbitrary'),
        text: pick('Le signe d’un axe de PCA peut être inversé sans changer le résultat. Interprétez les loadings et les scores comme une structure descriptive, et non comme une preuve causale ou statistique de groupe.', 'A PCA axis can be sign-flipped without changing the result. Treat loadings and scores as descriptive structure, not causal evidence or a group-comparison test.')
      });
    } else if (result.protocol?.objective === 'outcome') {
      const mode = result.protocol?.outcomeType || 'continuous';
      /** @type {Record<string,string>} */
      const effectTexts = {
        binary: pick('exp(effect) est un odds ratio.', 'exp(effect) is an odds ratio.'),
        count: pick('exp(effect) est un rate ratio.', 'exp(effect) is a rate ratio.'),
        survival: pick('exp(effect) est un hazard ratio.', 'exp(effect) is a hazard ratio.'),
        continuous: pick('effect est la variation attendue de l’outcome pour une unité de la variable omique.', 'effect is the expected outcome change per one feature unit.'),
        multiclass: pick('effect est l’écart maximal entre moyennes ajustées des classes.', 'effect is the maximum difference between adjusted class means.')
      };
      const effectText = effectTexts[mode] || '';
      items.push({
        title: pick('Association avec l’outcome', 'Outcome association'),
        text: pick('Chaque variable est testée avec le modèle déterminé par le type d’outcome. ', 'Each feature is tested with the model selected from the outcome type. ') + effectText
      });
      items.push({
        title: pick('Covariables', 'Covariates'),
        text: result.protocol?.covariateColumns?.length
          ? pick('Le modèle ajuste explicitement : ', 'The model explicitly adjusts for: ') + result.protocol.covariateColumns.join(', ') + '.'
          : pick('Aucune covariable supplémentaire n’a été sélectionnée.', 'No additional covariates were selected.')
      });
      if (result.predictiveOutcome?.status === 'ok') {
        const metrics = result.predictiveOutcome.metrics || {};
        const predictionText = mode === 'survival' && Number.isFinite(metrics.cIndex)
          ? pick(
              'La performance pronostique est mesurée hors échantillon par un C-index de Harrell de ' + metrics.cIndex.toFixed(3) + '. Cette valeur provient uniquement des folds externes de la nested cross-validation.',
              'Prognostic performance is assessed out of sample with a Harrell C-index of ' + metrics.cIndex.toFixed(3) + '. This value uses only outer-fold predictions from nested cross-validation.'
            )
          : mode === 'binary' && Number.isFinite(metrics.auc)
            ? pick(
                'La discrimination hors échantillon est résumée par une AUC de ' + metrics.auc.toFixed(3) + ', calculée sur les folds externes.',
                'Out-of-sample discrimination is summarized by an AUC of ' + metrics.auc.toFixed(3) + ', computed on outer-fold predictions.'
              )
            : pick(
                'La performance prédictive affichée provient des folds externes de la nested cross-validation et doit rester distincte des p-values d’association.',
                'Displayed predictive performance comes from outer nested-CV folds and should be kept distinct from association p-values.'
              );
        items.push({
          title: pick('Prédiction hors échantillon', 'Out-of-sample prediction'),
          text: predictionText
        });
      }
    } else if (result.protocol?.objective === 'time') {
      items.push({
        title: pick('Interaction groupe × temps', 'Group × time interaction'),
        text: pick(
          'L’effet natif R est la différence de pente temporelle entre les deux groupes, par unité de temps déclarée. Ce coefficient n’est ni un fold ratio ni un changement moyen entre deux visites.',
          'The native R effect is the difference in time slopes between the two groups per declared time unit. This coefficient is neither a fold ratio nor a mean change between two visits.'
        )
      });
      items.push({
        title: pick('Hypothèses longitudinales', 'Longitudinal assumptions'),
        text: pick(
          'Le modèle lmerTest suppose ici une évolution temporelle linéaire, une appartenance de groupe stable et le plan répété déclaré. Les fits singuliers, visites incomplètes et replis de pente aléatoire restent à examiner.',
          'The lmerTest model assumes a linear time trend, stable group membership and the declared repeated-measures design. Singular fits, incomplete visits and random-slope fallbacks still require review.'
        )
      });
    } else {
      items.push({
        title: pick('Effet biologique', 'Biological effect'),
        text: pick('L’effet indique la direction et l’amplitude du contraste après prétraitement et ajustement des facteurs techniques/covariables déclarés. Pour une échelle log2, le fold ratio est 2^effect.', 'The effect gives the direction and magnitude of the contrast after preprocessing and declared technical/covariate adjustment. On a log2 scale, fold ratio equals 2^effect.')
      });
      items.push({
        title: pick('Significativité', 'Statistical evidence'),
        text: pick('La colonne q BH corrige les tests multiples. q ≤ 0,10 est utilisé ici comme seuil exploratoire ; il ne remplace ni la taille d’effet ni la plausibilité biologique.', 'BH q-values adjust for multiple testing. q ≤ 0.10 is used here as an exploratory threshold; it does not replace effect size or biological plausibility.')
      });
    }

    const adjustedEntries = Object.entries(result.metadataSummary?.adjustment || {})
      .filter(([, value]) => value?.applied);
    if (adjustedEntries.length) {
      const adjusted = adjustedEntries.map(([layer]) => omicLabel(layer));
      const directOutcome = adjustedEntries.some(([, value]) => String(value?.method || '').includes('direct nuisance'));
      items.push({
        title: pick('Ajustement technique', 'Technical adjustment'),
        text: directOutcome
          ? pick('Batch et covariables sélectionnées entrent directement dans les modèles d’outcome pour : ', 'Batch and selected covariates enter the outcome models directly for: ') + adjusted.join(', ') + '.'
          : pick('Une résidualisation OLS a été appliquée avant l’inférence pour : ', 'OLS residualisation was applied before inference for: ') + adjusted.join(', ') + '.'
      });
    }

    if (result.crossOmics?.significantPairs > 0) {
      items.push({
        title: pick('Relations inter-omiques', 'Cross-omics relationships'),
        text: pick(
          result.crossOmics.significantPairs + ' paire(s) présentent un changement de corrélation avec q ≤ 0,10. Cela signale une relation condition-dépendante, pas une causalité.',
          result.crossOmics.significantPairs + ' pair(s) show a correlation change at q ≤ 0.10. This indicates condition-dependent coupling, not causality.'
        )
      });
    }

    const topPathway = result.reactome?.consensus?.[0];
    if (topPathway) {
      items.push({
        title: pick('Voies biologiques', 'Biological pathways'),
        text: pick(
          'La voie la mieux classée est « ' + topPathway.name + ' », soutenue par ' + topPathway.supportingLayers + ' couche(s) à FDR ≤ 0,10. Le classement Reactome reste exploratoire et dépend de l’univers d’annotation.',
          'The top-ranked pathway is “' + topPathway.name + '”, supported by ' + topPathway.supportingLayers + ' layer(s) at FDR ≤ 0.10. Reactome ranking remains exploratory and depends on the annotation universe.'
        )
      });
    }

    items.push({
      title: pick('Ce que le résultat ne prouve pas', 'What the result does not prove'),
      text: pick('Aucune causalité, mécanisme non mesuré ou vérité clinique n’est déduite automatiquement. Les résultats doivent être confrontés au design, à la qualité des données, à la taille d’échantillon et à la littérature.', 'No causality, unmeasured mechanism or clinical truth is inferred automatically. Results must be interpreted against study design, data quality, sample size and external evidence.')
    });
    return items;
  }

  $: inferredTranscriptomicsId = inferFeatureIdType('transcriptomics', matrixInfo.transcriptomics);
  $: inferredProteomicsId = inferFeatureIdType('proteomics', matrixInfo.proteomics);
  $: inferredMetabolomicsId = inferFeatureIdType('metabolomics', matrixInfo.metabolomics);
  $: selectedObjective = objectives[objective];
  $: omicsCount = omicLayers.filter((key) => Boolean(files[key])).length;
  $: requiredMappingsComplete = fieldDefinitions.filter((field) => field.required).every((field) => Boolean(columnMapping[field.key]));
  $: mappedSubjects = uniqueMappedFrom(metadataRows, columnMapping, 'subject_id').size;
  $: mappedSamples = uniqueMappedFrom(metadataRows, columnMapping, 'sample_id').size;
  $: mappedAssays = uniqueMappedFrom(metadataRows, columnMapping, 'assay_id').size;
  $: availableOutcomeTimepoints = [...uniqueMappedFrom(metadataRows, columnMapping, 'timepoint')];
  $: rnaseqAvailableConditions = columnMapping.condition
    ? [...new Set(metadataRows.map((row) =>
      String(row[columnMapping.condition] ?? '').trim()).filter(Boolean))].sort()
    : [];
  $: if (availableOutcomeTimepoints.length <= 1 && outcomeTimepoint) outcomeTimepoint = '';
  $: replicateGroups = technicalReplicateGroups(metadataRows, columnMapping);
  $: canonicalMappedColumns = new Set(Object.values(columnMapping).filter(Boolean));
  $: availableCovariateColumns = metadataHeaders.filter((header) => !canonicalMappedColumns.has(header));
  $: if (selectedCovariates.some((header) => !availableCovariateColumns.includes(header))) {
    selectedCovariates = selectedCovariates.filter((header) => availableCovariateColumns.includes(header));
  }

  $: outcomeMappingComplete = objective !== 'outcome' || outcomeType !== 'none' && (availableOutcomeTimepoints.length <= 1 || Boolean(outcomeTimepoint)) && (outcomeType === 'survival'
    ? Boolean(columnMapping.survival_time && columnMapping.survival_event)
    : Boolean(columnMapping.outcome));
  /** @type {any} */
  let pairedAudit = null;
  $: pairedAudit = objective === 'groups' && designType === 'paired' && files.transcriptomics && files.metadata
    ? auditPairedRnaStudy({rows:metadataRows,mapping:columnMapping,
        matrixAssays:matrixInfo.transcriptomics.sampleIds,
        valueType:transcriptomicsValues}) : null;
  /** @param {string} code */
  function pairedProblemMessage(code) {
    const problems = /** @type {Record<string,string[]>} */ ({
      raw_counts_only:['Le modèle RNA-seq exige de vrais comptages bruts, pas des TPM/FPKM.', 'The RNA-seq model requires genuine raw counts, not TPM/FPKM.'],
      missing_columns:['Associez sujet, prélèvement, mesure, omique, groupe et type de prélèvement aux colonnes.', 'Map subject, sample, assay, omic, condition and sample type to the correct columns.'],
      missing_transcriptomics_rows:['Aucune ligne RNA reconnue dans les métadonnées.', 'No RNA sample rows recognized in metadata.'],
      non_biological_sample:['Les blancs et QC ne sont pas des prélèvements biologiques appariés.', 'Blanks and QC are not paired biological samples.'],
      matrix_metadata_mismatch:['Les identifiants RNA de la matrice et des métadonnées diffèrent.', 'RNA matrix assay IDs differ from metadata.'],
      two_conditions_only:['Deux conditions exactement sont nécessaires.', 'Exactly two conditions are required.'],
      min_four_pairs:['Au moins quatre paires biologiques sont nécessaires pour ce contrôle logiciel.', 'At least four biological pairs are required for this software guard.'],
      missing_pair:['Un sujet ne possède pas les deux conditions.', 'A subject is missing one of the paired conditions.']
    });
    const text = problems[code] || ['Paires manquantes, identifiants doublonnés ou incohérents : vérifiez chaque sujet et chaque mesure.', 'Missing pairs or duplicate/inconsistent identifiers: check every subject and assay.'];
    return t(text[0],text[1]);
  }
  function downloadPairedExampleSheet(){
    const blob = new Blob(['\uFEFF'+pairedExampleSampleSheetCsv()], {type:'text/csv;charset=utf-8'});
    const url = URL.createObjectURL(blob);
    const link = document.createElement('a');
    link.href=url;
    link.download='EXEMPLE_FICTIF_4_PAIRES_metadata.csv';
    link.click();
    URL.revokeObjectURL(url);
  }
  $: objectiveOperational = designType !== 'crossover';
  $: singleOmic = omicsCount === 1;
  $: ready = omicsCount >= 1 && Boolean(files.metadata) && requiredMappingsComplete && outcomeMappingComplete && objectiveOperational;
  $: interpretationItems = analysisResult ? buildInterpretation(analysisResult, String($language || 'fr')) : [];
  $: analysisPlan = singleOmic && objective === 'explore'
    ? t('contrôle qualité → prétraitement documenté → ACP et variables principales de cette omique (exploratoire, sans intégration multi-omique)',
      'quality control → documented preprocessing → PCA and leading features of this omic (exploratory, no multi-omics integration)')
    : objective === 'explore'
    ? t('prétraitement → agrégation des réplicats → ajustement batch/covariables → standardisation par couche → ACP multi-blocs équilibrée → loadings → Reactome',
        'preprocessing → replicate aggregation → batch/covariate adjustment → within-layer scaling → balanced multi-block PCA → loadings → Reactome')
    : objective === 'outcome'
      ? t('prétraitement → agrégation → régression selon l’outcome + batch/covariables → BH-FDR → nested CV prédictive → Reactome',
          'preprocessing → aggregation → outcome-specific regression + batch/covariates → BH-FDR → nested predictive CV → Reactome')
      : objective === 'time'
        ? t('agrégation des réplicats → prétraitement → modèle condition × temps à intercept aléatoire sujet + batch/covariables → BH-FDR → intégration inter-omique → Reactome',
            'replicate aggregation → preprocessing → condition × time random-intercept model + batch/covariates → BH-FDR → cross-omics integration → Reactome')
        : t('agrégation des réplicats → prétraitement → ajustement batch/covariables → contraste de groupes → inférence → BH-FDR → intégration inter-omique → Reactome',
            'replicate aggregation → preprocessing → batch/covariate adjustment → group contrast → inference → BH-FDR → cross-omics integration → Reactome');
</script>

<svelte:window
  onmouseover={showGlobalHelp}
  onfocusin={showGlobalHelp}
  onmouseout={hideGlobalHelp}
  onfocusout={hideGlobalHelp}
></svelte:window>

<svelte:head>
  <title>{t('Outil multi-omique — PMx Explain', 'Multi-omics tool — PMx Explain')}</title>
  <meta name="robots" content="noindex,nofollow,noarchive" />
  <meta name="description" content="Guided analysis of one or more transcriptomics, proteomics and metabolomics layers." />
</svelte:head>

<section class="hero">
  <a class="tool-back" href={resolve('multiomics/')}>← {t('Présentation de l’outil', 'Tool overview')}</a>
  <p class="eyebrow">{t('Prototype expérimental · omique unique ou multi-omique', 'Experimental prototype · single-omics or multi-omics')}</p>
  <h1>{t('Analyser une ou plusieurs omiques et comprendre les résultats.', 'Analyze one or several omics layers and understand the results.')}</h1>
  <p class="lede">{t('Choisissez ce que vous voulez comprendre, ajoutez vos données et consultez les graphiques. Les vérifications scientifiques restent accessibles.', 'Choose your question, add your data and explore the figures. All scientific checks remain available.')}</p>
  <div class="quick-start" data-testid="multiomics-quick-start">
    <div class="quick-start-copy">
      <strong>{t('Découvrir en 1 clic', 'Explore in one click')}</strong>
      <span>{t('Choisissez une omique ou les trois. Chaque bouton charge des données synthétiques et calcule leurs résultats. Aucun fichier à préparer.', 'Choose one omics layer or all three. Each button loads synthetic example data and runs a real analysis. No files to prepare.')}</span>
      <small>{t('Exemple pédagogique (8 sujets) : il illustre la méthode, sans valider une découverte biologique.', 'Educational example (8 subjects): it illustrates the workflow, not a validated biological discovery.')}</small>
    </div>
    <div class="quick-demo-actions" aria-label={t('Choix de démonstration', 'Demonstration choices')}>
      <button class="btn btn-outline quick-demo-single" type="button" data-testid="multiomics-quick-demo-rna" disabled={demoLoading} onclick={() => loadDemoAndShowResults('transcriptomics')}>{t('Transcriptomique', 'Transcriptomics')}</button>
      <button class="btn btn-outline quick-demo-single" type="button" data-testid="multiomics-quick-demo-protein" disabled={demoLoading} onclick={() => loadDemoAndShowResults('proteomics')}>{t('Protéomique', 'Proteomics')}</button>
      <button class="btn btn-outline quick-demo-single" type="button" data-testid="multiomics-quick-demo-metabolite" disabled={demoLoading} onclick={() => loadDemoAndShowResults('metabolomics')}>{t('Métabolomique', 'Metabolomics')}</button>
      <button class="btn btn-primary quick-demo-multi" type="button" data-testid="multiomics-quick-demo" disabled={demoLoading} onclick={() => loadDemoAndShowResults()}>
        {t('Multi-omiques · les trois ensemble', 'Multi-omics · all three together')}
      </button>
    </div>
    {#if demoLoading}<p role="status">{t('Analyse de la démonstration en cours…', 'Analyzing the demo…')}</p>{/if}
    {#if demoError}<p role="alert" data-testid="multiomics-demo-error" class="error">{demoError}</p>{/if}
  </div>
  <details class="simple-disclosure">
    <summary>{t('Confidentialité des données', 'Data privacy')}</summary>
    <p>{t('Les calculs du navigateur restent sur votre appareil. Si les bases publiques sont activées, des identifiants moléculaires sont transmis à leurs services. Un moteur R connecté reçoit les fichiers à son adresse affichée dans les paramètres avancés (par défaut : votre propre ordinateur).', 'Browser calculations stay on your device. If public databases are enabled, molecular identifiers are sent to their services. A connected R engine receives the files at its configured address (by default your own computer).')}</p>
  </details>
</section>

<section class="simple-steps" aria-label={t('Étapes principales', 'Main steps')}>
  <span><b>1</b> {t('Question', 'Question')}</span>
  <span><b>2</b> {t('Étude', 'Study')}</span>
  <span><b>3</b> {t('Fichiers (1 à 3 omiques)', 'Files (1 to 3 omics)')}</span>
  <span><b>4</b> {t('Résultats', 'Results')}</span>
</section>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('1 · Votre question', '1 · Your question')}</p>
      <h2>{t('À quelle question l’expérience doit-elle répondre ?', 'What should the experiment answer?')}</h2>
    </div>

  </div>

  <div class="form-grid">
    <label class="wide">
      <span>{t('Objectif principal', 'Main objective')}</span>
      <select value={objective} onchange={chooseSimpleObjective}>
        <option value="explore">{t('Découvrir ce qui varie dans mes données', 'Explore what varies in my data')}</option>
        <option value="groups">{t('Comparer des groupes (ex. traité / témoin)', 'Compare groups (e.g. treatment / control)')}</option>
        <option value="outcome">{t('Relier les mesures à un résultat (ex. réponse)', 'Relate measurements to an outcome (e.g. response)')}</option>
        <option value="time">{t('Observer les changements dans le temps', 'Track changes over time')}</option>
      </select>
    </label>

    <details class="inline-optional">
      <summary>{t('Nommer cette analyse (facultatif)', 'Name this analysis (optional)')}</summary>
      <label>
        <span>{t('Nom de l’étude', 'Study name')}</span>
        <input bind:value={studyName} type="text" placeholder={t('Ex. Réponse au traitement', 'e.g. Treatment response')} />
      </label>
    </details>

    <label>
      <span>{t('Organisme', 'Organism')}</span>
      <select bind:value={organism}>
        <option value="human">{t('Humain', 'Human')}</option>
        <option value="mouse">{t('Souris', 'Mouse')}</option>
        <option value="rat">Rat</option>
        <option value="other">{t('Autre / à définir', 'Other / to define')}</option>
      </select>
    </label>
  </div>

  <details class="simple-disclosure">
    <summary>{t('Quelle méthode sera utilisée ?', 'Which method will be used?')}</summary>
    <strong>{$language === 'en' ? selectedObjective.methodEn : selectedObjective.methodFr}</strong>
    <p>{$language === 'en' ? selectedObjective.detailEn : selectedObjective.detailFr}</p>
  </details>
</section>


<section class="panel simple-study" data-testid="multiomics-simple-study">
  <div class="section-head"><div>
      <p class="eyebrow">{t('2 · Votre expérience', '2 · Your experiment')}</p>
      <h2>{t('Comment les données ont-elles été recueillies ?', 'How were the data collected?')}</h2>
  </div></div>
  <p class="simple-hint">{t('Ces réponses permettent de choisir une analyse adaptée. Les réglages techniques sont juste en dessous.', 'These answers help select a suitable analysis. Technical settings are just below.')}</p>
  <div class="form-grid">
    <label>
      <span>{t('Quel type d’étude ?', 'Type of study?')}</span>
      <select value={studySetting} onchange={chooseSimpleStudySetting} aria-label={t('Type d’étude', 'Study type')}>
        <option value="clinical_observational">{t('Participants humains', 'Human participants')}</option>
        <option value="clinical_interventional">{t('Essai ou intervention chez l’humain', 'Human trial or intervention')}</option>
        <option value="animal">{t('Étude chez l’animal', 'Animal study')}</option>
        <option value="cell">{t('Cellules en laboratoire', 'Cell cultures')}</option>
        <option value="organoid">{t('Organoïdes / tissus', 'Organoids / tissues')}</option>
        <option value="other">{t('Autre', 'Other')}</option>
      </select>
    </label>
    <label>
      <span>{t('Comment comparer les échantillons ?', 'How are samples compared?')}</span>
      <select value={designType} onchange={chooseSimpleDesign} data-testid="multiomics-simple-design" aria-label={t('Organisation des échantillons', 'Sample organization')}>
        <option value="independent">{t('Groupes de sujets différents', 'Different subjects in each group')}</option>
        <option value="paired">{t('Deux mesures liées chez les mêmes sujets', 'Paired measurements from the same subjects')}</option>
        <option value="repeated">{t('Mêmes sujets suivis à plusieurs moments', 'Same subjects measured over time')}</option>
        <option value="crossover">{t('Traitements successifs chez les mêmes sujets', 'Crossover treatment study')}</option>
      </select>
    </label>
    {#if objective === 'groups' || objective === 'time'}
      <label>
        <span>{t('Combien de groupes ?', 'How many groups?')}</span>
        <select bind:value={groupCount}>
          <option value="1">{t('Un seul', 'One')}</option>
          <option value="2">{t('Deux', 'Two')}</option>
          <option value="3plus">{t('Trois ou plus', 'Three or more')}</option>
        </select>
      </label>
    {/if}
    {#if designType === 'repeated' || objective === 'time'}
      <label>
        <span>{t('Combien de moments de mesure ?', 'How many measurement times?')}</span>
        <select bind:value={timepointCount}>
          <option value="2">2</option>
          <option value="3">3</option>
          <option value="4">{t('4 ou plus', '4 or more')}</option>
        </select>
      </label>
    {/if}
    {#if objective === 'outcome'}
      <label>
        <span>{t('Quel résultat voulez-vous étudier ?', 'What outcome are you studying?')}</span>
        <select bind:value={outcomeType} data-testid="multiomics-simple-outcome">
          <option value="none">{t('Choisir…', 'Select…')}</option>
          <option value="binary">{t('Oui / non (ex. réponse)', 'Yes / no (e.g. response)')}</option>
          <option value="continuous">{t('Un nombre (ex. concentration)', 'A number (e.g. concentration)')}</option>
          <option value="multiclass">{t('Plusieurs catégories', 'Several categories')}</option>
          <option value="count">{t('Nombre d’événements', 'Number of events')}</option>
          <option value="survival">{t('Temps avant un événement', 'Time until an event')}</option>
        </select>
      </label>
      {#if availableOutcomeTimepoints.length > 1}
        <label>
          <span>{t('À quel moment les échantillons ont-ils été prélevés ?', 'At which time were samples collected?')}</span>
          <select bind:value={outcomeTimepoint} aria-label={t('Temps des mesures moléculaires', 'Omics time point used')}>
            <option value="">— {t('choix nécessaire', 'selection required')} —</option>
            {#each availableOutcomeTimepoints as timepoint}<option value={timepoint}>{timepoint}</option>{/each}
          </select>
        </label>
      {/if}
    {/if}
    <label>
      <span>{t('Les analyses portent-elles sur les mêmes prélèvements ?', 'Were the same samples used for each omics measurement?')}</span>
      <select bind:value={sampleOverlap}>
        <option value="same_specimen">{t('Oui, même prélèvement', 'Yes, same sample')}</option>
        <option value="same_subject">{t('Mêmes sujets, prélèvements différents', 'Same subjects, different samples')}</option>
        <option value="partial">{t('Seulement certains prélèvements', 'Only some matched samples')}</option>
        <option value="unpaired">{t('Non, prélèvements différents', 'No, different samples')}</option>
        <option value="unknown">{t('Je ne sais pas', 'I am not sure')}</option>
      </select>
    </label>
  </div>
  {#if studySetting === 'animal' && organism === 'human'}
    <p class="simple-warning" role="alert">{t('Vous avez choisi une étude animale mais l’organisme est encore « Humain ». Modifiez l’organisme dans la première étape avant d’analyser.', 'You selected an animal study but the organism is still set to Human. Change the organism in the first step before analysis.')}</p>
  {/if}
  {#if objective === 'groups' && designType === 'paired'}
    <div class="simple-hint" data-testid="multiomics-paired-help">
      <strong>{t('Un sujet = une paire', 'One subject = one pair')}</strong>
      <p>{t('Dans votre tableau, même subject_id pour les deux conditions, mais sample_id et assay_id distincts pour chaque prélèvement. Quatre paires signifient quatre unités biologiques indépendantes, pas huit.', 'Use the same subject_id for both conditions but unique sample_id and assay_id for each specimen. Four pairs mean four independent biological units, not eight.')}</p>
      <p>{t('Pour des comptages RNA-seq, choisissez « R requis » : modèle bloqué par sujet, sans supposer les observations indépendantes.', 'For RNA-seq counts choose Require R: the model blocks on subject, rather than assuming independent observations.')}</p>
      <button type="button" class="btn btn-outline btn-small" data-testid="multiomics-paired-template-download" onclick={downloadPairedExampleSheet}>{t('Télécharger un tableau fictif de 4 paires (CSV)', 'Download a fictional four-pair sample sheet (CSV)')}</button>
      <p>{t('Exemple de structure uniquement : remplacez tous les identifiants et ne confondez pas ce tableau avec des mesures réelles. Les noms assay_id doivent correspondre aux colonnes de votre matrice RNA.', 'File structure example only: replace every identifier. This contains no actual measurements. assay_id values must match the RNA matrix columns.')}</p>
      {#if pairedAudit?.ok}
        <p role="status" data-testid="multiomics-paired-count">{pairedAudit.nPairs} {t('sujets appariés', 'paired subjects')} · {pairedAudit.nBiologicalSamples} {t('prélèvements', 'samples')} · {pairedAudit.nAssays} {t('mesures', 'assays')}</p>
      {:else if pairedAudit}
        <p role="alert" data-testid="multiomics-paired-problem">{pairedProblemMessage(pairedAudit.code)}</p>
      {/if}
    </div>
  {/if}

  {#if designType === 'crossover'}
    <p class="simple-warning" role="alert">{t('Cette étude nécessite une analyse particulière qui n’est pas encore disponible. Le lancement sera bloqué.', 'This study needs a model that is not yet available; running is blocked.')}</p>
  {/if}
  {#if objective === 'time' && designType !== 'repeated'}
    <p class="simple-warning" role="status">{t('Pour un suivi dans le temps, vérifiez le type de comparaison et les mesures répétées dans les paramètres avancés.', 'For a time-course study, check the design and repeated measurements in advanced settings.')}</p>
  {/if}
  <p class="simple-hint">{t('Le modèle, les corrections statistiques et les hypothèses seront détaillés dans le rapport.', 'The model, statistical corrections and assumptions will be documented in the report.')}</p>
</section>

<details class="advanced-parameters" data-testid="multiomics-advanced-settings">
  <summary><strong>{t('Paramètres avancés', 'Advanced settings')}</strong> <span>{t('Plan expérimental détaillé, types de fichiers et contrôles techniques', 'Detailed study design, file types and technical checks')}</span></summary>

  <section class="panel">
    <div class="section-head">
      <div>
        <p class="eyebrow">{t('Étape 2 · Protocole', 'Step 2 · Protocol')}</p>
        <h2>{t('Décrivez le design avec des questions fermées', 'Describe the design with closed questions')}</h2>
      </div>
      <p>{t('Les réponses définissent la structure statistique avant toute sélection de méthode omique.', 'The answers define the statistical structure before any omics method is selected.')}</p>
    </div>

    <div class="form-grid">
      <label>
        <span>{t('Cadre de l’étude', 'Study setting')}</span>
        <select bind:value={studySetting}>
          <option value="clinical_observational">{t('Humain · cohorte observationnelle', 'Human · observational cohort')}</option>
          <option value="clinical_interventional">{t('Humain · intervention / essai', 'Human · intervention / trial')}</option>
          <option value="animal">{t('Expérience animale', 'Animal experiment')}</option>
          <option value="cell">{t('Culture cellulaire / in vitro', 'Cell culture / in vitro')}</option>
          <option value="organoid">{t('Organoïde / modèle ex vivo', 'Organoid / ex vivo model')}</option>
          <option value="other">{t('Autre', 'Other')}</option>
        </select>
      </label>

      <label>
        <span>{t('Unité biologique indépendante', 'Independent biological unit')}</span>
        <select bind:value={unitType}>
          <option value="participant">{t('Participant humain', 'Human participant')}</option>
          <option value="animal">Animal</option>
          <option value="culture">{t('Culture indépendante / réplicat biologique', 'Independent culture / biological replicate')}</option>
          <option value="other">{t('Autre unité expérimentale', 'Other experimental unit')}</option>
        </select>
      </label>

      <label>
        <span>{t('Structure du design', 'Design structure')}</span>
        <select bind:value={designType}>
          <option value="independent">{t('Groupes / unités indépendants', 'Independent groups / units')}</option>
          <option value="paired">{t('Échantillons appariés', 'Paired samples')}</option>
          <option value="crossover">{t('Crossover / comparaison intra-sujet', 'Crossover / within-subject comparison')}</option>
          <option value="repeated">{t('Mesures répétées', 'Repeated measurements')}</option>
        </select>
      </label>

      <label>
        <span>{t('Nombre de groupes comparés', 'Number of comparison groups')}</span>
        <select bind:value={groupCount}>
          <option value="1">{t('Un groupe / exploratoire', 'One group / exploratory')}</option>
          <option value="2">{t('Deux groupes', 'Two groups')}</option>
          <option value="3plus">{t('Trois groupes ou plus', 'Three or more groups')}</option>
        </select>
      </label>

      <label>
        <span>{t('Design longitudinal ?', 'Longitudinal design?')}</span>
        <select bind:value={longitudinal}>
          <option value="no">{t('Non', 'No')}</option>
          <option value="yes">{t('Oui', 'Yes')}</option>
        </select>
      </label>

      <label>
        <span>{t('Nombre de temps', 'Number of time points')}</span>
        <select bind:value={timepointCount}>
          <option value="1">1</option>
          <option value="2">2</option>
          <option value="3">3</option>
          <option value="4">{t('4 ou plus', '4 or more')}</option>
        </select>
      </label>

      <label>
        <span>{t('Les omiques proviennent-elles du même prélèvement biologique ?', 'Do the omics come from the same biological specimen?')}</span>
        <select bind:value={sampleOverlap}>
          <option value="same_specimen">{t('Oui, même prélèvement pour toutes les omiques', 'Yes, same specimen for all omics')}</option>
          <option value="same_subject">{t('Même sujet, prélèvements différents', 'Same subject, different specimens')}</option>
          <option value="partial">{t('Partiellement apparié', 'Partially matched')}</option>
          <option value="unpaired">{t('Échantillons différents / non appariés', 'Different / unpaired samples')}</option>
          <option value="unknown">{t('Inconnu', 'Unknown')}</option>
        </select>
      </label>

      <label>
        <span>{t('Réplicats techniques attendus ?', 'Technical replicates expected?')}</span>
        <select bind:value={technicalReplicatesExpected}>
          <option value="no">{t('Non', 'No')}</option>
          <option value="yes">{t('Oui', 'Yes')}</option>
          <option value="unknown">{t('Inconnu', 'Unknown')}</option>
        </select>
      </label>

      <label>
        <span>{t('Batches techniques connus ?', 'Known technical batches?')}</span>
        <select bind:value={batchKnown}>
          <option value="yes">{t('Oui', 'Yes')}</option>
          <option value="no">{t('Non', 'No')}</option>
          <option value="unknown">{t('Inconnu', 'Unknown')}</option>
        </select>
      </label>

      <label>
        <span>{t('Type d’outcome principal', 'Primary outcome type')}</span>
        <select bind:value={outcomeType}>
          <option value="none">{t('Pas d’outcome / exploratoire', 'No outcome / exploratory')}</option>
          <option value="binary">{t('Binaire', 'Binary')}</option>
          <option value="multiclass">{t('Multiclasse', 'Multiclass')}</option>
          <option value="continuous">{t('Continu', 'Continuous')}</option>
          <option value="survival">{t('Temps jusqu’à événement / survie', 'Time-to-event / survival')}</option>
          <option value="count">{t('Comptage', 'Count')}</option>
        </select>
      </label>

      {#if objective === 'outcome' && availableOutcomeTimepoints.length > 1}
        <label>
          <span>{t('Temps omique utilisé dans le modèle d’outcome', 'Omics time point used in the outcome model')}</span>
          <select bind:value={outcomeTimepoint}>
            <option value="">— {t('sélection requise', 'selection required')} —</option>
            {#each availableOutcomeTimepoints as timepoint}
              <option value={timepoint}>{timepoint}</option>
            {/each}
          </select>
          <small>{t('Choisissez explicitement le prélèvement moléculaire qui entre dans le modèle. Pour une prédiction pronostique, il s’agit souvent du baseline.', 'Explicitly choose the molecular sampling time used by the model. For prognostic prediction, this is often baseline.')}</small>
        </label>
      {/if}

      <label>
        <span>{t('Covariables importantes disponibles ?', 'Important covariates available?')}</span>
        <select bind:value={covariatesAvailable}>
          <option value="yes">{t('Oui', 'Yes')}</option>
          <option value="no">{t('Non', 'No')}</option>
        </select>
      </label>

      <label>
        <span>{t('Certains sujets peuvent-ils manquer une couche omique entière ?', 'Some subjects may miss an entire omics layer?')}</span>
        <select bind:value={partialOmicsExpected}>
          <option value="no">{t('Non / non attendu', 'No / not expected')}</option>
          <option value="yes">{t('Oui', 'Yes')}</option>
          <option value="unknown">{t('Inconnu', 'Unknown')}</option>
        </select>
      </label>

      <label>
      <span>{t('Nombre d’unités biologiques', 'Number of biological units')} <small>{t('contrôle optionnel', 'optional consistency check')}</small></span>
      <input bind:value={subjectCount} type="number" min="1" placeholder="e.g. 48" />
      </label>

      <label>
      <span>{t('Nom de la variable de groupe', 'Name of the group variable')} <small>{t('optionnel', 'optional')}</small></span>
      <input bind:value={groupVariable} type="text" placeholder="e.g. treatment" />
      </label>

      <label>
      <span>{t('Nom de l’outcome', 'Name of the outcome')} <small>{t('optionnel', 'optional')}</small></span>
      <input bind:value={outcome} type="text" placeholder="e.g. response" />
      </label>
    </div>

    <div class="identity-grid">
      <article>
        <code>subject_id</code>
        <strong>{t('Unité biologique', 'Biological unit')}</strong>
        <p>{t('Même participant / animal / culture indépendante entre les visites et prélèvements.', 'Same participant / animal / independent culture across visits and samples.')}</p>
      </article>
      <article>
        <code>sample_id</code>
        <strong>{t('Prélèvement biologique', 'Biological specimen')}</strong>
        <p>{t('Le prélèvement physique. Les dosages RNA, protéines et métabolites du même prélèvement partagent cet ID.', 'The physical specimen. RNA, protein and metabolite assays from the same specimen share this ID.')}</p>
      </article>
      <article>
        <code>assay_id</code>
        <strong>{t('Mesure technique', 'Technical measurement')}</strong>
        <p>{t('Le nom de colonne utilisé dans une matrice omique. Les réplicats techniques ont des assay_id distincts mais le même sample_id.', 'The column name used in an omics matrix. Technical replicates have distinct assay IDs but the same sample ID.')}</p>
      </article>
    </div>
  </section>

  <section class="panel">
    <div class="section-head">
      <div>
        <p class="eyebrow">{t('Étape 3 · Type de données', 'Step 3 · Data type')}</p>
        <h2>{t('Indiquez ce que représente chaque matrice', 'Tell the app what each matrix represents')}</h2>
      </div>
      <p>{t('Cela évite d’interpréter à tort les mêmes nombres comme counts bruts, intensités normalisées ou concentrations absolues.', 'This prevents the same numbers from being interpreted incorrectly as raw counts, normalised intensities or absolute concentrations.')}</p>
    </div>

  <div class="api-summary" data-testid="input-support-boundary">
    <span><strong>RNA</strong> {t('matrice de counts/TPM/expression : oui · FASTQ/BAM : prétraitement amont requis', 'counts/TPM/expression matrix: yes · FASTQ/BAM: upstream preprocessing required')}</span>
    <span><strong>{t('Protéines', 'Proteins')}</strong> {t('matrice LFQ/TMT/DIA/spectral counts : oui · fichiers MS vendor : quantification amont requise', 'LFQ/TMT/DIA/spectral-count matrix: yes · vendor MS files: upstream quantification required')}</span>
    <span><strong>{t('Métabolites', 'Metabolites')}</strong> {t('matrice : oui · mzML/mzXML : pipeline xcms local disponible · fichiers vendor : conversion mzML requise', 'matrix: yes · mzML/mzXML: local xcms pipeline available · vendor files: mzML conversion required')}</span>
    </div>
    <p class="note">{t('L’outil n’étiquette jamais FASTQ/BAM ou fichiers protéomiques bruts comme « analysés » : il attend les sorties quantitatives des pipelines spécialisés. La seule branche brute intégrée actuellement est la LC/GC-MS ouverte via xcms.', 'The tool never labels FASTQ/BAM or raw proteomics files as “analysed”: it expects quantitative outputs from specialised upstream pipelines. The only integrated raw-data branch currently is open-format LC/GC-MS through xcms.')}</p>

    <div class="omics-question-grid">
      <article>
        <h3>{t('Transcriptomique', 'Transcriptomics')}</h3>
        <label>
          <span>{t('Plateforme', 'Platform')}</span>
          <select bind:value={transcriptomicsPlatform}>
            <option value="bulk_rnaseq">Bulk RNA-seq</option>
            <option value="microarray">Microarray</option>
            <option value="targeted">{t('Panel d’expression ciblé', 'Targeted expression panel')}</option>
            <option value="processed">{t('Matrice déjà prétraitée', 'Already processed matrix')}</option>
          </select>
        </label>
        <label>
          <span>{t('Valeurs', 'Values')}</span>
          <select bind:value={transcriptomicsValues}>
            <option value="raw_counts">{t('Counts entiers bruts', 'Raw integer counts')}</option>
            <option value="tpm">{t('Abondance TPM / type FPKM', 'TPM / FPKM-like abundance')}</option>
            <option value="normalized">{t('Expression normalisée', 'Normalised expression')}</option>
            <option value="log_expression">{t('Expression log-transformée', 'Log-transformed expression')}</option>
            <option value="unknown">{t('Inconnu', 'Unknown')}</option>
          </select>
        </label>
        <label>
          <span>{t('Identifiant des variables', 'Feature identifier')}</span>
          <select bind:value={transcriptomicsIdType}>
            <option value="ensembl_gene">Ensembl gene ID</option>
            <option value="gene_symbol">{t('Symbole de gène', 'Gene symbol')}</option>
            <option value="entrez">Entrez Gene ID</option>
            <option value="unknown">{t('Inconnu / détection automatique', 'Unknown / detect automatically')}</option>
          </select>
        </label>
        {#if matrixInfo.transcriptomics.rowIds.length}
        <small>Detected from uploaded features: <strong>{inferredTranscriptomicsId.type}</strong> ({Math.round(inferredTranscriptomicsId.confidence * 100)}%).</small>
        {/if}
      </article>

      <article>
        <h3>{t('Protéomique', 'Proteomics')}</h3>
        <label>
          <span>{t('Plateforme', 'Platform')}</span>
          <select bind:value={proteomicsPlatform}>
            <option value="label_free">Label-free LC-MS/MS</option>
            <option value="tmt">TMT / multiplexed</option>
            <option value="dia">DIA</option>
            <option value="targeted">{t('Protéomique ciblée', 'Targeted proteomics')}</option>
            <option value="processed">{t('Matrice déjà prétraitée', 'Already processed matrix')}</option>
          </select>
        </label>
        <label>
          <span>{t('Valeurs', 'Values')}</span>
          <select bind:value={proteomicsValues}>
            <option value="lfq_intensity">{t('LFQ / intensité', 'LFQ / intensity')}</option>
            <option value="log_intensity">{t('Intensité logarithmique', 'Log intensity')}</option>
            <option value="spectral_count">{t('Comptages spectraux', 'Spectral counts')}</option>
            <option value="normalized">{t('Abondance normalisée', 'Normalised abundance')}</option>
            <option value="unknown">{t('Inconnu', 'Unknown')}</option>
          </select>
        </label>
        <label>
          <span>{t('Identifiant des variables', 'Feature identifier')}</span>
          <select bind:value={proteomicsIdType}>
            <option value="uniprot">UniProt accession</option>
            <option value="gene_symbol">{t('Symbole de gène', 'Gene symbol')}</option>
            <option value="ensembl_protein">Ensembl protein ID</option>
            <option value="unknown">{t('Inconnu / détection automatique', 'Unknown / detect automatically')}</option>
          </select>
        </label>
        {#if matrixInfo.proteomics.rowIds.length}
        <small>Detected from uploaded features: <strong>{inferredProteomicsId.type}</strong> ({Math.round(inferredProteomicsId.confidence * 100)}%).</small>
        {/if}
      </article>

      <article>
        <h3>{t('Métabolomique', 'Metabolomics')}</h3>
        <label>
          <span>{t('Acquisition', 'Acquisition')}</span>
          <select bind:value={metabolomicsPlatform}>
            <option value="untargeted_lcms">{t('LC-MS non ciblée', 'Untargeted LC-MS')}</option>
            <option value="targeted_lcms">{t('LC-MS ciblée', 'Targeted LC-MS')}</option>
            <option value="gcms">GC-MS</option>
            <option value="nmr">NMR</option>
            <option value="processed">{t('Matrice déjà prétraitée', 'Already processed matrix')}</option>
          </select>
        </label>
        <label>
          <span>{t('Valeurs', 'Values')}</span>
          <select bind:value={metabolomicsValues}>
            <option value="peak_area">{t('AUC MS (aire sous le pic chromatographique)', 'MS AUC (chromatographic peak area)')}</option>
            <option value="normalized">{t('Abondance normalisée', 'Normalised abundance')}</option>
            <option value="concentration">{t('Concentration absolue', 'Absolute concentration')}</option>
            <option value="log_abundance">{t('Abondance logarithmique', 'Log abundance')}</option>
            <option value="unknown">{t('Inconnu', 'Unknown')}</option>
          </select>
        </label>
      {#if ['untargeted_lcms','targeted_lcms','gcms'].includes(metabolomicsPlatform)}
          <div class="ms-qc-controls">
          <strong>{t('QC MS avancé', 'Advanced MS QC')}
            <span class="help-tip" tabindex="0" data-tooltip={t('Utilise sample_type pour distinguer biological / qc / blank et injection_order pour corriger la dérive. Les blanks et pooled-QC sont exclus de l’inférence biologique.', 'Uses sample_type to distinguish biological / qc / blank and injection_order for drift correction. Blank and pooled-QC injections are excluded from biological inference.')}>?</span>
            </strong>
            <label>
              <span>{t('Gestion des signaux présents dans les blanks', 'Blank-associated signal handling')}</span>
              <select bind:value={msBlankFilter}>
                <option value="flag">{t('Signaler seulement · recommandé par défaut', 'Flag only · recommended default')}</option>
                <option value="remove">{t('Exclure selon le ratio déclaré', 'Exclude using the declared ratio')}</option>
                <option value="no">{t('Désactiver', 'Disable')}</option>
              </select>
              <small>{t('Le ratio sert à repérer les variables potentiellement contaminées. Leur exclusion automatique est une décision plus forte : utilisez-la seulement si elle correspond à votre SOP et conservez une analyse de sensibilité.', 'The ratio flags potentially contaminated features. Automatic exclusion is a stronger decision: use it only when consistent with your SOP and retain a sensitivity analysis.')}</small>
            </label>
            <label>
              <span>{t('Ratio biologique / blank minimal', 'Minimum biological / blank ratio')}</span>
            <input type="number" min="1" step="0.5" bind:value={msBlankFold} />
            </label>
            <label>
              <span>{t('Filtre RSD pooled-QC', 'Pooled-QC RSD filter')}</span>
              <select bind:value={msQcRsdFilter}>
                <option value="yes">{t('Oui', 'Yes')}</option>
                <option value="no">{t('Non', 'No')}</option>
              </select>
            </label>
            <label>
              <span>{t('RSD QC maximal', 'Maximum QC RSD')}</span>
            <input type="number" min="0.05" max="1" step="0.05" bind:value={msQcRsdThreshold} />
            </label>
            <label>
              <span>{t('Correction de dérive', 'Drift correction')}</span>
              <select bind:value={msDriftCorrection}>
                <option value="yes">{t('Oui · ≥5 QC (linéaire JS), ≥8 QC (LOESS R)', 'Yes · ≥5 QC (linear JS), ≥8 QC (R LOESS)')}</option>
                <option value="no">{t('Non', 'No')}</option>
              </select>
            </label>
            <label>
              <span>{t('Missing MNAR', 'MNAR missingness')}</span>
              <select bind:value={msMnarStrategy}>
                <option value="none">{t('Aucune imputation', 'No imputation')}</option>
                <option value="left_censored">{t('Bas de distribution déterministe', 'Deterministic low-tail imputation')}</option>
              </select>
            </label>
          </div>
        {/if}

        <label>
          <span>{t('Identifiant des variables', 'Feature identifier')}</span>
          <select bind:value={metabolomicsIdType}>
            <option value="chebi">ChEBI</option>
            <option value="hmdb">HMDB</option>
            <option value="kegg_compound">KEGG compound</option>
            <option value="pubchem">PubChem CID</option>
            <option value="inchikey">InChIKey (UniChem → ChEBI)</option>
            <option value="name">{t('Nom du métabolite', 'Metabolite name')}</option>
            <option value="mz_rt">{t('m/z + temps de rétention uniquement', 'm/z + retention time only')}</option>
            <option value="unknown">{t('Inconnu / détection automatique', 'Unknown / detect automatically')}</option>
          </select>
        </label>
        {#if matrixInfo.metabolomics.rowIds.length}
        <small>Detected from uploaded features: <strong>{inferredMetabolomicsId.type}</strong> ({Math.round(inferredMetabolomicsId.confidence * 100)}%).</small>
        {/if}
      </article>
    </div>
  </section>
</details>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('Pour démarrer', 'Getting started')}</p>
      <h2>{t('Un exemple pour comprendre, un modèle si nécessaire', 'Try an example or use a template')}</h2>
    </div>
    <p>{t('Vous pouvez charger vos propres tableaux, sans les renommer si les colonnes sont reconnues.', 'You can upload your own tables without renaming columns if they are recognized.')}</p>
  </div>

  <div class="contract">
    <details class="templates-details">
      <summary>{t('Besoin d’un modèle de fichier ?', 'Need a file template?')}</summary>
      <div>
        <h3>{t('Modèles à télécharger', 'Downloadable templates')}</h3>
        <p>{t('Une ligne = une mesure. Ce format gère les couches omiques absentes, les temps répétés et les réplicats techniques sans changer de schéma.', 'One row = one assay. This handles missing omics layers, repeated time points and technical replicates without changing the schema.')}</p>
        <pre>subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,outcome,survival_time,survival_event,sample_type,injection_order,covariate_1</pre>
        <div class="actions">
          <button
            class="btn btn-primary"
            type="button"
            onclick={downloadGeneratedTemplate}
          >{t('Générer un template depuis mon protocole', 'Generate template from my protocol')}</button>

          <a
            class="btn btn-outline"
            href={asset('multiomics/metadata_template.csv')}
            download
          >{t('Template générique de métadonnées', 'Generic metadata template')}</a>

          <a
            class="btn btn-outline"
            href={asset('multiomics/transcriptomics_template.csv')}
            download
          >{t('Template matrice RNA', 'RNA matrix template')}</a>

          <a
            class="btn btn-outline"
            href={asset('multiomics/proteomics_template.csv')}
            download
          >{t('Template matrice protéines', 'Protein matrix template')}</a>

          <a
            class="btn btn-outline"
            href={asset('multiomics/metabolomics_template.csv')}
            download
          >{t('Template matrice métabolites', 'Metabolite matrix template')}</a>

          <a
            class="btn btn-outline"
            data-testid="raw-ms-manifest-template"
            href={asset('multiomics/templates/raw_ms_manifest.csv')}
            download
          >{t('Template MS brute', 'Raw MS manifest')}</a>

          <a
            class="btn btn-outline"
            data-testid="external-validation-template"
            href={asset('multiomics/templates/external_validation_predictions_binary.csv')}
            download
          >{t('Template validation externe', 'External validation template')}</a>
        </div>
      </div>
    </details>

    <div class="demo-card">
      <p class="eyebrow">{t('Exemple prêt à explorer', 'Ready-to-run example')}</p>
      <h3>{t('Traitement × temps, trois omiques', 'Treatment × time, three omics')}</h3>
      <p>{t('8 sujets, 2 visites, 3 types de mesures. La démo exécute les mêmes calculs que vos fichiers.', '8 subjects, 2 visits, 3 measurement types. The demo runs the same calculations as your own files.')}</p>
      <button class="btn btn-primary" type="button" data-testid="multiomics-load-demo" disabled={demoLoading} onclick={() => loadDemoAndShowResults()}>{t('Charger les trois omiques', 'Load all three omics')}</button>
      <div class="actions" data-testid="multiomics-single-demo-buttons">
        <button class="btn btn-outline" type="button" data-testid="multiomics-demo-rna" disabled={demoLoading} onclick={() => loadDemoAndShowResults('transcriptomics')}>{t('Démo RNA seul', 'RNA-only demo')}</button>
        <button class="btn btn-outline" type="button" data-testid="multiomics-demo-protein" disabled={demoLoading} onclick={() => loadDemoAndShowResults('proteomics')}>{t('Démo protéines seules', 'Proteomics-only demo')}</button>
        <button class="btn btn-outline" type="button" data-testid="multiomics-demo-metabolite" disabled={demoLoading} onclick={() => loadDemoAndShowResults('metabolomics')}>{t('Démo métabolites seuls', 'Metabolomics-only demo')}</button>
      </div>
      {#if demoError}<p class="error">{demoError}</p>{/if}
      <p class="simple-hint">{t('Démo synthétique : 8 sujets, groupes traité/témoin et deux temps. Les résultats montrent le fonctionnement, pas une découverte biologique validée. Choisissez RNA, protéines, métabolites ou les trois.', 'Synthetic demo: 8 subjects, treated/control groups and two visits. The outputs demonstrate the workflow, not a validated biological finding. Choose RNA, proteins, metabolites or all three.')}</p>
      <details class="simple-disclosure">
        <summary>{t('Télécharger les données de démonstration', 'Download demo data')}</summary>
        <div class="demo-links">
          <a
            href={asset('multiomics/demo_metadata.csv')}
            download
          >{t('métadonnées', 'metadata')}</a>

          <a
            href={asset('multiomics/demo_transcriptomics.csv')}
            download
          >RNA</a>

          <a
            href={asset('multiomics/demo_proteomics.csv')}
            download
          >{t('protéines', 'protein')}</a>

          <a
            href={asset('multiomics/demo_metabolomics.csv')}
            download
          >{t('métabolites', 'metabolites')}</a>
        </div>
      </details>
    </div>
  </div>

  <details class="optional-workflows" data-testid="advanced-workflows">
    <summary>{t('Options expertes : données MS brutes et validation externe', 'Expert options: raw MS data and external validation')}</summary>
    <div class="contract">
      <div class="demo-card">
        <p class="eyebrow">{t('Workflow avancé · données MS brutes', 'Advanced workflow · raw MS data')}</p>
        <h3>{t('mzML / mzXML → matrice de features via xcms', 'mzML / mzXML → feature matrix with xcms')}</h3>
        <p>{t('Le navigateur n’essaie pas de traiter directement les fichiers vendor. Convertissez Thermo/Waters/Agilent/Bruker en mzML avec ProteoWizard, puis exécutez le pipeline local R : détection centWave, alignement obiwarp, groupement, gap filling et matrice d’aires de pics. Cette étape crée des features m/z–RT ; elle n’identifie pas les métabolites.', 'The browser does not pretend to process vendor files directly. Convert Thermo/Waters/Agilent/Bruker files to mzML with ProteoWizard, then run the local R pipeline: centWave peak picking, obiwarp retention-time alignment, grouping, gap filling and peak-area matrix. This creates m/z–RT features; it does not identify metabolites.')}</p>
        <div class="actions">
          <a
            class="btn btn-outline"
            href={asset('multiomics/templates/raw_ms_manifest.csv')}
            download
          >{t('Manifeste exemple', 'Example manifest')}</a>

          <a
            class="btn btn-outline"
            href={asset('multiomics/templates/raw_ms_parameters.json')}
            download
          >{t('Paramètres xcms', 'xcms parameters')}</a>
        </div>
        <p class="note"><code>Rscript multiomics-engine/run_raw_ms.R raw_ms_manifest.csv raw_ms_output raw_ms_parameters.json</code></p>
        <p class="note">{t('La matrice metabolomics_peak_area.csv obtenue revient ensuite dans l’étape d’import ci-dessous, avec sample_type et injection_order pour le QC blank / pooled-QC / dérive.', 'The resulting metabolomics_peak_area.csv then enters the upload step below, together with sample_type and injection_order for blank / pooled-QC / drift QC.')}</p>
      </div>

      <div class="demo-card">
        <p class="eyebrow">{t('Workflow avancé · validation externe', 'Advanced workflow · external validation')}</p>
        <h3>{t('Évaluer un modèle figé sur une cohorte indépendante', 'Evaluate a frozen model on an independent cohort')}</h3>
        <p>{t('La validation externe est volontairement séparée de l’entraînement : aucune sélection de variable, aucun tuning et aucun ré-entraînement ne sont autorisés sur la cohorte de validation. L’outil calcule les métriques adaptées au type de critère et un bootstrap déterministe lorsque c’est estimable.', 'External validation is deliberately separated from training: no feature selection, tuning or refitting is allowed on the validation cohort. The validator computes outcome-specific metrics and a deterministic bootstrap when estimable.')}</p>
        <div class="actions">
          <a
            class="btn btn-outline"
            href={asset('multiomics/templates/external_validation_predictions_binary.csv')}
            download
          >{t('Prédictions exemple', 'Example predictions')}</a>

          <a
            class="btn btn-outline"
            href={asset('multiomics/templates/external_validation_binary.json')}
            download
          >{t('Configuration', 'Configuration')}</a>
        </div>
        <p class="note"><code>Rscript multiomics-engine/run_external_validation.R external_validation_predictions_binary.csv external_validation_binary.json external_validation_output</code></p>
        <p class="note">{t('Le logiciel n’emploie le terme « validation externe » que si l’indépendance de la cohorte est explicitement déclarée ; cette indépendance doit rester documentée par la provenance de l’étude.', 'The software uses the term “external validation” only when cohort independence is explicitly asserted; that independence must still be documented by study provenance.')}</p>
      <ExternalValidationPanel backendUrl={referenceBackendUrl} language={$language ?? 'fr'} />
      </div>
    </div>
  </details>

  <details class="dictionary">
    <summary>{t('Dictionnaire détaillé des métadonnées', 'Detailed metadata dictionary')}</summary>
    <div class="dictionary-table">
      <div><b>Column</b><b>Meaning</b><b>Example</b><b>Rule</b></div>
      <div><code>subject_id</code><span>{t('Unité biologique indépendante', 'Independent biological unit')}</span><span>P001</span><span>Same value across all visits/omics for one participant, animal or culture.</span></div>
      <div><code>sample_id</code><span>Physical biological specimen</span><span>P001_T0</span><span>Same value only when assays come from the same specimen.</span></div>
      <div><code>assay_id</code><span>Technical measurement / run</span><span>RNA001</span><span>Must match the corresponding matrix row/column identifier exactly.</span></div>
      <div><code>omic</code><span>Measured layer</span><span>transcriptomics</span><span>Canonical values: transcriptomics, proteomics, metabolomics.</span></div>
      <div><code>condition</code><span>Experimental group</span><span>treatment</span><span>Use one consistent vocabulary across subjects.</span></div>
      <div><code>timepoint</code><span>Visit / experimental time</span><span>T12</span><span>Required for longitudinal designs.</span></div>
      <div><code>batch</code><span>Technical batch</span><span>RNA_B1</span><span>Keep assay-specific batches even when different omics use different batches.</span></div>
      <div><code>technical_replicate</code><span>Repeated technical assay</span><span>1</span><span>Distinct assay_id, same sample_id + omic.</span></div>
      <div><code>outcome</code><span>{t('Phénotype / critère principal', 'Primary phenotype / endpoint')}</span><span>responder</span><span>{t('Utilisé pour les outcomes binaires, continus, multiclasse ou de comptage.', 'Used for binary, continuous, multiclass or count outcomes.')}</span></div>
      <div><code>survival_time</code><span>{t('Temps de suivi / survie', 'Survival / follow-up time')}</span><span>365</span><span>{t('Requis uniquement pour la branche survie.', 'Required only for the survival branch.')}</span></div>
      <div><code>survival_event</code><span>{t('Événement 0/1', 'Event indicator 0/1')}</span><span>1</span><span>{t('1 = événement observé, 0 = censuré.', '1 = observed event, 0 = censored.')}</span></div>
      <div><code>covariate_*</code><span>{t('Covariable explicite', 'Explicit covariate')}</span><span>age</span><span>{t('Toute colonne supplémentaire peut être sélectionnée manuellement pour ajustement.', 'Any extra column can be selected manually for adjustment.')}</span></div>
    </div>
  </details>

  <details class="aliases">
    <summary>{t('Comment fonctionne la reconnaissance automatique des colonnes', 'How automatic column recognition works')}</summary>
    <p>{t('Les noms de colonnes sont normalisés (casse, espaces, tirets et accents ignorés), puis comparés à une liste contrôlée d’alias. La reconnaissance automatique n’est acceptée que si une seule colonne correspond au champ.', 'Column names are normalised (case, spaces, hyphens and accents ignored), then compared with a controlled alias list. Automatic recognition is only accepted when exactly one column matches a field.')}</p>
    <div class="alias-grid">
      {#each fieldDefinitions as field}
        <article>
          <strong>{field.key}</strong>
          <span>{field.required ? 'required' : 'optional'}</span>
          <p>{field.aliases.join(', ')}</p>
        </article>
      {/each}
    </div>
    <p class="note">{t('La liste d’alias propose un mapping ; elle ne redéfinit jamais silencieusement l’étude. Les champs ambigus ou absents doivent être mappés manuellement.', 'The alias list proposes a mapping; it does not silently redefine the study. Ambiguous or missing fields must be mapped manually.')}</p>
  </details>
</section>

<section class="panel">
  <div class="section-head">
    <div>
      <p class="eyebrow">{t('3 · Vos fichiers', '3 · Your files')}</p>
      <h2>{t('Ajoutez vos tableaux de données', 'Add your data files')}</h2>
    </div>

  </div>

  <div class="uploads">
    <label class:loaded={files.metadata}>
      <strong>{t('Tableau des échantillons', 'Sample information table')}</strong>
      <span>{t('Qui a été prélevé, quand et dans quel groupe ?', 'Who was sampled, when, and in which group?')}</span>
      <input type="file" accept=".csv,.tsv,.txt" onchange={(event) => selectFile('metadata', event)} />
      <small>{files.metadata ? files.metadata.name : 'No file selected'}</small>
    </label>

    <label class:loaded={files.transcriptomics}>
      <strong>{t('Transcriptomique', 'Transcriptomics')}</strong>
      <span>{t('Quantités mesurées pour les gènes', 'Measured amounts of gene expression')}</span>
      <input type="file" accept=".csv,.tsv,.txt" onchange={(event) => selectFile('transcriptomics', event)} />
      <small>{files.transcriptomics ? files.transcriptomics.name : 'No file selected'}</small>
    </label>

    <label class:loaded={files.proteomics}>
      <strong>{t('Protéomique', 'Proteomics')}</strong>
      <span>{t('Quantités mesurées pour les protéines', 'Measured amounts of proteins')}</span>
      <input type="file" accept=".csv,.tsv,.txt" onchange={(event) => selectFile('proteomics', event)} />
      <small>{files.proteomics ? files.proteomics.name : 'No file selected'}</small>
    </label>

    <label class:loaded={files.metabolomics}>
      <strong>{t('Métabolomique', 'Metabolomics')}</strong>
      <span>{t('Aires de pics (AUC) LC-MS / GC-MS : matrice ou export long avec signal, injection et aire', 'LC-MS / GC-MS peak areas (AUC): matrix or long export with feature, injection and area')}</span>
      <input type="file" data-testid="multiomics-ms-auc-upload" accept=".csv,.tsv,.txt" onchange={(event) => selectFile('metabolomics', event)} />
      <small>{files.metabolomics ? files.metabolomics.name : 'No file selected'}</small>

      <a
        href={asset('multiomics/ms_peak_areas_example.csv')}
        download
      >{t('Exemple CSV AUC (données fictives)', 'Example AUC CSV (illustrative synthetic data)')}</a>

      <small>{t('AUC = aire intégrée du pic chromatographique, non AUC pharmacocinétique. Exportez un CSV/TSV avec une ligne par signal et injection (feature_id, assay_id, auc), ou une matrice signaux × injections. Une aire de pic n’identifie pas à elle seule une molécule.', 'AUC = chromatographic peak area, not pharmacokinetic exposure AUC. Export a CSV/TSV with one feature and injection per row (feature_id, assay_id, auc), or a feature × injection matrix. A peak area alone does not identify a molecule.')}</small>
    </label>
    <label class:loaded={annotationFileName}>
      <strong>{t('Annotations des pics MS (facultatif)', 'MS peak annotations (optional)')}</strong>
      <span>{t('Si vos pics ont des identifiants internes, associez-les explicitement à ChEBI.', 'If peaks have internal IDs, map them explicitly to ChEBI.')}</span>
      <input type="file" data-testid="multiomics-ms-annotation-upload" accept=".csv,.tsv,.txt" onchange={selectMetaboliteAnnotations} />
      <small>{annotationFileName || t('CSV : feature_id,chebi_id', 'CSV: feature_id,chebi_id')}</small>

      <a
        href={asset('multiomics/ms_peak_annotations_example.csv')}
        download
      >{t('Exemple d’annotations', 'Annotation template')}</a>

      <small>{t('Correspondances déclarées par l’utilisateur, non validées chimiquement. Aucun pic inconnu n’est identifié automatiquement ; les statistiques restent inchangées.', 'User-declared mappings, not independently verified chemical identities. Unknown peaks are never identified automatically; statistical fits remain unchanged.')}</small>
      {#if annotationError}<small class="annotation-error" data-testid="multiomics-ms-annotation-error">{annotationError}</small>{/if}
      {#if metaboliteAnnotations.length}<small data-testid="multiomics-ms-annotation-count">{metaboliteAnnotations.length} {t('correspondances importées', 'annotations imported')}</small>{/if}
    </label>
  </div>

  {#if singleOmic && !files.metadata && objective === 'explore'}
    <div class="ms-import-message" data-testid="multiomics-single-sheet-helper">
      <strong>{t('Exploration sans tableau d’échantillons ?', 'Exploring without a sample sheet?')}</strong>
      <p>{t('Possible seulement si chaque colonne correspond à un sujet biologique différent, mesuré une seule fois. Aucun groupe ni réplicat ne sera deviné. Ce raccourci ne permet pas une comparaison de groupes ou une analyse temporelle.', 'Only if every matrix column is a different biological subject measured once. No group or replicate information is guessed. This shortcut cannot run group or longitudinal inference.')}</p>
      <label>
        <input type="checkbox" bind:checked={confirmIndependentAssays} />
        {t('Je confirme : une colonne = un sujet biologique indépendant.', 'I confirm: one column = one independent biological subject.')}
      </label>
      <button type="button" class="btn btn-outline btn-small" disabled={!confirmIndependentAssays} data-testid="multiomics-generate-explore-metadata" onclick={generateSingleOmicExploratorySheet}>
        {t('Créer le tableau minimal pour explorer', 'Create minimal exploratory sample sheet')}
      </button>
    </div>
  {/if}
  {#if exploratorySheetGenerated}
    <p class="ms-import-message" data-testid="multiomics-generated-sheet-notice">{t('Tableau créé à partir de votre déclaration : identifiants conservés, aucun groupe inféré. Les résultats restent descriptifs.', 'Sample sheet created from your declaration: IDs preserved, no groups inferred. Results remain descriptive.')}</p>
  {/if}
  <p class="ms-import-message" role="status" aria-live="polite" data-testid="multiomics-ms-auc-import-status" data-import-state={msImportState}>
    {#if msImportState === 'loading'}{t('Lecture du fichier MS…', 'Reading MS file…')}
    {:else if msImportState === 'parsed'}{t('Fichier MS analysé : vérification des échantillons…', 'MS file parsed: validating samples…')}
    {:else if msImportState === 'ready'}{msImportMessage}
    {:else if msImportState === 'error'}{t('Échec de l’import MS : ', 'MS import failed: ')}{analysisError}
    {/if}
  </p>
  {#if singleOmic}<p class="ms-import-message" data-testid="multiomics-single-omic-notice">{t('Mode omique unique : exploration, groupes, temps ou outcome selon votre plan. Aucune intégration inter-omique ne sera calculée.', 'Single-omics mode: exploration, groups, time or outcome, depending on the study. No cross-omics integration will be computed.')}</p>{/if}
  {#if msOriginalFile && msAnnotationCsv}
    <div class="ms-import-message" data-testid="multiomics-ms-annotations-preserved">
      <span>{t('Annotations MS conservées séparément : ', 'MS annotations preserved separately: ')}{msAnnotationFields.join(', ')}</span>
      <button class="btn btn-outline btn-small" type="button" onclick={downloadMsAnnotationTable}>{t('Exporter les annotations techniques', 'Export technical annotations')}</button>
    </div>
  {/if}
  {#if files.transcriptomics || files.proteomics || files.metabolomics}
    <div class="simple-file-values" data-testid="multiomics-simple-values">
      <strong>{t('Que représentent les nombres de vos fichiers ?', 'What do the numbers in your files mean?')}</strong>
      <p>{t('Indiquez le format donné par votre laboratoire : il influence les calculs et ne peut pas être deviné.', 'Choose the format supplied by your laboratory: it affects calculations and cannot be guessed.')}</p>
      <div class="form-grid">
        {#if files.transcriptomics}
          <label><span>{t('Gènes : type de valeurs', 'Genes: value type')}</span>
            <select bind:value={transcriptomicsValues} aria-label={t('Type de valeurs des gènes', 'Gene value type')}>
              <option value="raw_counts">{t('Comptages bruts', 'Raw read counts')}</option>
              <option value="tpm">TPM / FPKM</option>
              <option value="normalized">{t('Valeurs normalisées', 'Normalized values')}</option>
              <option value="log_expression">{t('Valeurs logarithmiques', 'Log-transformed values')}</option>
              <option value="unknown">{t('Je ne sais pas (à vérifier)', 'Unsure (review required)')}</option>
            </select>
          </label>
        {/if}
        {#if files.proteomics}
          <label><span>{t('Protéines : type de valeurs', 'Proteins: value type')}</span>
            <select bind:value={proteomicsValues} aria-label={t('Type de valeurs des protéines', 'Protein value type')}>
              <option value="lfq_intensity">{t('Intensités (LFQ)', 'Intensities (LFQ)')}</option>
              <option value="log_intensity">{t('Intensités logarithmiques', 'Log-transformed intensities')}</option>
              <option value="spectral_count">{t('Comptages de spectres', 'Spectral counts')}</option>
              <option value="normalized">{t('Valeurs normalisées', 'Normalized values')}</option>
              <option value="unknown">{t('Je ne sais pas (à vérifier)', 'Unsure (review required)')}</option>
            </select>
          </label>
        {/if}
        {#if files.metabolomics}
          <label><span>{t('Métabolites : type de valeurs', 'Metabolites: value type')}</span>
            <select bind:value={metabolomicsValues} aria-label={t('Type de valeurs des métabolites', 'Metabolite value type')}>
              <option value="peak_area">{t('AUC MS : aire de pic / intensité relative', 'MS AUC: peak area / relative intensity')}</option>
              <option value="normalized">{t('Valeurs normalisées', 'Normalized values')}</option>
              <option value="concentration">{t('Concentrations', 'Concentrations')}</option>
              <option value="log_abundance">{t('Valeurs logarithmiques', 'Log-transformed values')}</option>
              <option value="unknown">{t('Je ne sais pas (à vérifier)', 'Unsure (review required)')}</option>
            </select>
          </label>
        {/if}
      </div>
      <small>{t('Pour les identifiants et filtres spécifiques, ouvrez les paramètres avancés au-dessus.', 'For identifier types and specific filters, open Advanced settings above.')}</small>
    </div>
  {/if}

  {#if metadataError}
    <p class="error">{metadataError}</p>
  {/if}

  {#if metadataHeaders.length}
    <div class="mapping">
      <div class="mapping-head">
        <div>
          <h3>{t('Vérifiez les colonnes reconnues', 'Check recognized columns')}</h3>
          <p>{t('Si une correspondance est incorrecte, sélectionnez la bonne colonne.', 'If a match is wrong, select the correct column.')}</p>
        </div>
        <span class:ok={requiredMappingsComplete}>{requiredMappingsComplete ? t('Colonnes principales reconnues', 'Main columns recognized') : t('Colonnes à compléter', 'Columns need attention')}</span>
      </div>

      <div class="mapping-grid">
        {#each fieldDefinitions as field}
          {#if field.required || field.key === 'condition' && ['groups', 'time'].includes(objective) || field.key === 'timepoint' && (longitudinal === 'yes' || objective === 'time') || field.key === 'outcome' && objective === 'outcome' || objective === 'outcome' && outcomeType === 'survival' && ['survival_time', 'survival_event'].includes(field.key)}
            <label>
            <span><strong>{fieldLabel(field.key)}</strong> <small>{field.required ? t('requis', 'required') : t('optionnel', 'optional')}</small></span>
            <select value={columnMapping[field.key]} onchange={(event) => setMapping(field.key, event.currentTarget.value)}>
                <option value="">— {t('non mappé', 'not mapped')} —</option>
                {#each metadataHeaders as header}
                  <option value={header}>{header}</option>
                {/each}
              </select>
            </label>
          {/if}
        {/each}
      </div>
      <details class="simple-disclosure" data-testid="multiomics-more-columns">
        <summary>{t('Autres colonnes et corrections (facultatif)', 'Other columns and adjustments (optional)')}</summary>
        <p>{t('Si votre étude dépend de l’âge, du sexe ou d’une autre variable, sélectionnez-la explicitement ci-dessous.', 'If your study depends on age, sex or another variable, explicitly select it below.')}</p>
        <div class="mapping-grid">
          {#each fieldDefinitions as field}
            {#if !(field.required || field.key === 'condition' && ['groups', 'time'].includes(objective) || field.key === 'timepoint' && (longitudinal === 'yes' || objective === 'time') || field.key === 'outcome' && objective === 'outcome' || objective === 'outcome' && outcomeType === 'survival' && ['survival_time', 'survival_event'].includes(field.key))}
              <label>
            <span><strong>{fieldLabel(field.key)}</strong> <small>{field.required ? t('requis', 'required') : t('optionnel', 'optional')}</small></span>
            <select value={columnMapping[field.key]} onchange={(event) => setMapping(field.key, event.currentTarget.value)}>
                  <option value="">— {t('non mappé', 'not mapped')} —</option>
                  {#each metadataHeaders as header}
                    <option value={header}>{header}</option>
                  {/each}
                </select>
              </label>
            {/if}
          {/each}
        </div>

        {#if covariatesAvailable === 'yes'}
          <div class="covariate-picker">
            <div>
              <h3>{t('Covariables à ajuster', 'Covariates to adjust')}</h3>
              <p>{t('Sélectionnez uniquement les variables de confusion ou facteurs de design à ajuster. Elles entrent explicitement dans le modèle ; elles ne sont pas devinées par l’outil.', 'Select only confounders or design variables that should be adjusted. They enter the model explicitly; the tool does not guess them.')}</p>
            </div>
            {#if availableCovariateColumns.length}
              <div class="covariate-options">
                {#each availableCovariateColumns as header}
                  <label>
                  <input type="checkbox" checked={selectedCovariates.includes(header)} onchange={() => toggleCovariate(header)} />
                    <span>{header}</span>
                  </label>
                {/each}
              </div>
            {:else}
              <p class="muted">{t('Aucune colonne supplémentaire disponible après le mapping canonique.', 'No extra columns remain after canonical mapping.')}</p>
            {/if}
            <p class="note">
            {t('Batches techniques :', 'Technical batches:')}
            {columnMapping.batch
                ? t(' la colonne batch mappée est ajustée automatiquement si plusieurs batches sont présents et non totalement confondus.', ' the mapped batch column is adjusted automatically when multiple non-confounded batches are present.')
                : t(' aucune colonne batch n’est mappée.', ' no batch column is mapped.')}
            </p>
          </div>
        {/if}
      </details>
    </div>

    <details class="advanced-results" data-testid="multiomics-sample-checks">
      <summary>{t('Détail des échantillons et correspondances', 'Sample counts and matching details')}</summary>
      <div class="validation">
        <article>
          <span>{t('Unités biologiques', 'Biological units')}</span>
          <strong>{mappedSubjects || '—'}</strong>
          {#if subjectCount && mappedSubjects && subjectCount !== mappedSubjects}
            <small class="warning">Declared {subjectCount}; metadata contains {mappedSubjects}.</small>
          {:else}
            <small>{t('subject_id uniques mappés.', 'Unique mapped subject IDs.')}</small>
          {/if}
        </article>
        <article>
          <span>{t('Prélèvements biologiques', 'Biological samples')}</span>
          <strong>{mappedSamples || '—'}</strong>
          <small>{t('Prélèvements uniques entre visites / conditions.', 'Unique specimens across visits / conditions.')}</small>
        </article>
        <article>
          <span>{t('Mesures techniques', 'Assays')}</span>
          <strong>{mappedAssays || '—'}</strong>
          <small>{t('Identifiants uniques de mesures techniques.', 'Unique technical measurement IDs.')}</small>
        </article>
        <article>
          <span>{t('Groupes de réplicats techniques', 'Technical replicate groups')}</span>
          <strong>{replicateGroups.length}</strong>
        <small>{replicateGroups.length ? replicateGroups.map((item) => item.key.replace('::', ' / ')).slice(0, 3).join(', ') : 'None detected from repeated sample + omic pairs.'}</small>
        </article>
      </div>

      <div class="matrix-checks">
        {#each omicLayers as layer}
          {@const match = matrixMatch(layer)}
          <article>
            <div>
              <strong>{omicLabel(layer)}</strong>
            <span>{matrixInfo[layer].sampleIds.length ? `${matrixInfo[layer].sampleIds.length} data columns` : 'not loaded'}</span>
            </div>
            {#if matrixInfo[layer].sampleIds.length}
            <p><b>{match.matched}/{match.expected}</b> expected assay IDs matched.</p>
              {#if match.orientation === 'transposed'}
                <small class="warning">{t('Les assay_id correspondent mieux aux lignes qu’aux colonnes : la matrice semble transposée.', 'Assay IDs match rows better than columns: the matrix appears transposed.')}</small>
              {/if}
            {#if match.missing.length}<small class="warning">Missing from matrix: {match.missing.slice(0, 5).join(', ')}</small>{/if}
            {#if match.extra.length && match.orientation !== 'transposed'}<small class="warning">Not declared in metadata: {match.extra.slice(0, 5).join(', ')}</small>{/if}
            {:else}
              <p class="muted">{t('Chargez une matrice pour valider les assay_id.', 'Load a matrix to validate assay IDs.')}</p>
            {/if}
          </article>
        {/each}
      </div>
    </details>
  {/if}

  <div class="status" class:ready>
    <strong>{omicsCount}/3 {t('types de mesures chargés', 'measurement types loaded')}</strong>
    <span>{ready
      ? t('Le contrat de données et l’objectif scientifique permettent de lancer l’analyse.', 'The data contract and scientific objective are sufficient to run the analysis.')
      : !objectiveOperational
        ? t('Les designs crossover restent bloqués tant que les effets de période et de séquence ne sont pas modélisés.', 'Crossover designs remain blocked until period and sequence effects are modelled.')
        : !outcomeMappingComplete
          ? t('Complétez le type et le mapping de l’outcome requis.', 'Complete the required outcome type and mapping.')
          : t('Chargez au moins une couche omique et mappez subject_id, sample_id, assay_id et omic.', 'Load at least one omics layer and map subject_id, sample_id, assay_id and omic.')}</span>
  </div>

  <details class="simple-disclosure portable-project-panel" data-testid="multiomics-project-import-panel">
    <summary>{t('Reprendre un projet enregistré', 'Reopen a saved project')}</summary>
    <p>{t('Réimportez votre projet PMx au format JSON. Les empreintes SHA-256 protègent les fichiers contre les modifications. Rien n’est envoyé à un serveur et vous choisissez quand relancer l’analyse.', 'Reimport your PMx JSON project. SHA-256 verifies the original files. Nothing is sent to a server; you choose when to rerun.')}</p>
    <label for="multiomics-project-import-file">{t('Fichier de projet (.json)', 'Project file (.json)')}</label>
    <input id="multiomics-project-import-file" data-testid="multiomics-project-import-file" type="file" accept=".json,application/json" disabled={portableProjectBusy} onchange={importPortableProject} />
    {#if portableProjectMessage}<p role="status" class="ms-import-message">{portableProjectMessage}</p>{/if}
    {#if portableProjectError}<p role="alert" class="annotation-error">{portableProjectError}</p>{/if}
  </details>
    <div class="run-box">
    <div>
      <p class="eyebrow">{t('4 · Analyser', '4 · Analyze')}</p>
      <h3>{t('Obtenir mes résultats', 'Get my results')}</h3>
      <p>{t('Par défaut, tout reste dans votre navigateur. Activer explicitement ChEBI/Reactome transmet certains identifiants moléculaires ; activer le moteur R transmet vos matrices et métadonnées à son adresse configurée.', 'By default, computation stays in your browser. Opting in to ChEBI/Reactome transmits some molecular identifiers; connecting the R engine sends matrices and metadata to its configured address.')}</p>
      <label class="reference-backend-box" data-testid="multiomics-analysis-intent">
        <strong>{t('Quel usage pour ces résultats ?', 'How will you use these results?')}</strong>
        <select bind:value={analysisIntent} aria-label={t('Objectif de validité', 'Scientific evidence objective')}
          onchange={() => { if (analysisIntent === 'confirmatory') referenceBackendMode = 'required'; }}>
          <option value="exploratory">{t('Explorer mes données (recommandé pour débuter)', 'Explore my data (recommended for beginners)')}</option>
          <option value="confirmatory">{t('Préparer une analyse confirmatoire (R requis + revue indépendante)', 'Prepare a confirmatory analysis (R required + independent review)')}</option>
        </select>
        <small>{analysisIntent === 'confirmatory'
          ? t('Le moteur R local devient obligatoire. Les tests contrôlent l’éligibilité, mais ne certifient jamais automatiquement vos résultats. Les fichiers seront envoyés au backend R configuré après lancement de l’analyse.', 'The configured local R engine is required. Checks assess eligibility but never automatically certify findings. Files are sent to the configured R backend when you run the analysis.')
          : t('Le navigateur applique des méthodes exploratoires et parfois des approximations ou simplifications statistiques : ses p-values, q-values et effets peuvent différer des méthodes R de référence. Les analyses R natives (DESeq2, edgeR/limma selon les données) sont préférables pour une inférence appropriée au modèle, mais leur validité dépend du plan, des hypothèses et du contrôle qualité. Aucune installation R n’est nécessaire pour explorer.', 'The browser uses exploratory models and sometimes statistical approximations or simplifications: its p-values, q-values and effects can differ from native reference R methods. Native R methods (DESeq2, edgeR/limma where applicable) are preferable for model-appropriate inference, but validity still depends on design, assumptions and quality checks. No R installation is needed to explore.')}</small>
      </label>
      {#if analysisIntent === 'confirmatory' && !(
        objective === 'groups' ||
        objective === 'time' ||
        (objective === 'outcome' && outcomeType === 'survival')
      )}
        <p role="alert" class="simple-hint" data-testid="multiomics-confirmatory-projection-unsupported">
          {t('Le mode confirmatoire affiche uniquement les questions disposant d’un tableau R natif complet : comparaison de groupes, interaction groupe × temps et survie. Les autres modèles d’outcome restent exploratoires tant que leur projection R complète n’est pas qualifiée.',
             'Confirmatory display is limited to objectives with a complete native R table: group comparison, group-by-time interaction and survival. Other outcome models remain exploratory until their full R projection is qualified.')}
        </p>
      {/if}
      {#if analysisIntent !== 'confirmatory'}
        <p role="note" class="simple-hint" data-testid="multiomics-browser-approximation-warning">
          {t('Important : les résultats du navigateur servent à repérer des tendances. Selon la méthode et les données, certains calculs simplifient les hypothèses du modèle ou utilisent des approximations ; les chiffres ne sont donc pas nécessairement ceux de DESeq2, edgeR ou limma en R. Pour une analyse de référence, exécutez le backend R, vérifiez le plan expérimental et téléchargez le code effectivement utilisé.',
             'Important: browser results are intended to explore trends. Depending on the data and method, calculations may simplify model assumptions or use approximations; the numbers therefore need not equal native DESeq2, edgeR or limma results. For reference analysis, run the R backend, review study design and download the code actually used.')}
        </p>
      {/if}
      {#if analysisIntent === 'confirmatory' && files.transcriptomics && transcriptomicsValues === 'raw_counts'}
        <p class="simple-hint" role="note" data-testid="multiomics-rnaseq-fdr-calibration-warning">
          {t('Point scientifique : certains scénarios RNA-seq simulés sans effet réel ont montré plus de faux positifs que les 5 % attendus sous DESeq2. La correction BH et le calcul R ne suffisent pas à certifier la calibration de votre étude : confrontez le modèle, le plan et les hypothèses avec un statisticien avant toute conclusion confirmatoire.',
            'Scientific limitation: some null RNA-seq simulations showed more false positives than the nominal 5% with DESeq2. BH adjustment and execution in R do not certify calibration for your study. Review the model, study design and assumptions with a statistician before confirmatory claims.')}
        </p>
      {/if}
      {#if analysisIntent === 'confirmatory' && objective === 'groups' && files.transcriptomics && transcriptomicsValues === 'raw_counts'}
        <fieldset class="reference-backend-box" data-testid="multiomics-rnaseq-contrast">
          <legend><strong>{t('Quel effet RNA-seq souhaitez-vous tester ?', 'Which RNA-seq group effect do you want to test?')}</strong></legend>
          <label>
            {t('Groupe de référence (dénominateur)', 'Reference group (denominator)')}
            <select bind:value={rnaseqReference} aria-label={t('Groupe RNA-seq de référence', 'RNA-seq reference group')}>
              <option value="">{t('Sélectionner un groupe', 'Select a group')}</option>
              {#each rnaseqAvailableConditions as group}
                <option value={group}>{group}</option>
              {/each}
            </select>
          </label>
          <label>
            {t('Groupe comparé (numérateur)', 'Comparison group (numerator)')}
            <select bind:value={rnaseqComparison} aria-label={t('Groupe RNA-seq comparé', 'RNA-seq comparison group')}>
              <option value="">{t('Sélectionner un groupe', 'Select a group')}</option>
              {#each rnaseqAvailableConditions as group}
                <option value={group}>{group}</option>
              {/each}
            </select>
          </label>
          <small>{t('Le log2FC sera calculé comme « groupe comparé / groupe de référence ». Les deux valeurs proviennent des métadonnées et doivent être différentes. Ce choix est conservé dans le projet reproductible.', 'The log2FC is comparison / reference. Both group names come from the sample sheet and must differ. Your choice is stored in the reproducible project.')}</small>
        </fieldset>
      {/if}
        <details class="advanced-parameters" data-testid="multiomics-local-r-setup">
          <summary>{t('Comment activer le moteur R ? (guide débutant)', 'How to enable R? (beginner guide)')}</summary>
          <p>{t('Pour explorer les données et suivre la démonstration, aucune installation n’est nécessaire. Pour les analyses de référence, GitHub Pages ne peut pas exécuter R : celui-ci doit être lancé sur votre propre ordinateur.', 'No installation is needed to explore data or run the demo. For reference analyses, GitHub Pages cannot execute R: you must start it on your own computer.')}</p>
          <ol>
            <li>{t('Installer R et Git sur votre ordinateur.', 'Install R and Git on your computer.')}</li>
            <li>{t('Ouvrir un terminal et copier ces commandes une par une :', 'Open a terminal and run these commands one by one:')}</li>
          </ol>
          <pre style="max-width:100%;overflow-x:auto;white-space:pre;" tabindex="0"><code>git clone https://github.com/rberrah/MultiOmics.git
cd MultiOmics
Rscript multiomics-engine/install_backend_dependencies.R
Rscript multiomics-engine/reference_doctor.R rnaseq
Rscript multiomics-engine/run_backend.R</code></pre>
        <p>{t('La commande de diagnostic vérifie les vrais packages sans installer ni envoyer de données. Remplacez « rnaseq » par « groups », « longitudinal », « mofa2 », « diablo » ou « rawms » selon le type d’analyse. Un contrôle vert n’est pas une certification statistique.', 'The diagnostic checks installed packages without installing software or uploading data. Replace “rnaseq” with “groups”, “longitudinal”, “mofa2”, “diablo” or “rawms” for your method family. A successful check is not statistical certification.')}</p>
        <p>{t('Laisser le terminal R ouvert, revenir ici et ouvrir « Options de calcul et connexions externes ». Choisir « Auto · utiliser si disponible » (le mode initial « Navigateur uniquement » désactive les tests R), conserver l’adresse http://127.0.0.1:8787 et cliquer sur « Tester ». Le statut doit indiquer que R est disponible. Ne sélectionner « R requis » que si le backend fonctionne et qu’une analyse de référence est indispensable.', 'Keep the R terminal open, return here and open Analysis options and external services. Choose Auto · use when available (the initial Browser only mode disables R checks), keep http://127.0.0.1:8787 and click Check. The status should report that R is available. Select Require R only when the backend is operational and reference analysis is essential.')}</p>
        <p>{t('Confidentialité : vos matrices ne sont envoyées au serveur R que lorsque vous lancez l’analyse avec ce moteur activé. Gardez le serveur sur votre ordinateur (127.0.0.1). Ne connectez pas de serveur distant non sécurisé à des données patients.', 'Privacy: your matrices are sent to the configured R server only when you launch analysis with that engine enabled. Keep the server on your computer (127.0.0.1). Never send patient data to an unsecured remote server.')}</p>
          <p><a href="https://github.com/rberrah/MultiOmics/blob/main/multiomics-engine/README.md" target="_blank" rel="noopener noreferrer">{t('Lire les instructions complètes et les limitations du moteur R', 'Read full R setup instructions and limitations')}</a></p>
      </details>
      <details class="simple-disclosure" data-testid="multiomics-run-options">
        <summary>{t('Options de calcul et connexions externes', 'Analysis options and external services')}</summary>
        <label class="inline-check">
          <input type="checkbox" bind:checked={resolveIdentifiers} />
          <span>{t('Résoudre les métabolites non canoniques sélectionnés avec ChEBI avant l’analyse de voies', 'Resolve selected non-canonical metabolite labels with ChEBI before pathway analysis')}</span>
        </label>
        <label class="inline-check">
          <input type="checkbox" bind:checked={useReactome} />
          <span>{t('Interroger Reactome uniquement avec les identifiants moléculaires sélectionnés', 'Query Reactome with selected molecular identifiers only')}</span>
        </label>

        <div class="reference-backend-box">
          <div>
            <strong>{t('Backend R de référence', 'Reference R backend')}</strong>
          <span class="help-tip" tabindex="0" data-tooltip={t('En mode Auto, l’outil teste le backend local. S’il répond, il exécute les implémentations de référence applicables : DESeq2, limma, lmerTest, fgsea, MOFA2 et/ou DIABLO. GitHub Pages ne peut pas démarrer R lui-même.', 'In Auto mode, the tool probes the local backend. If available, it runs applicable reference implementations: DESeq2, limma, lmerTest, fgsea, MOFA2 and/or DIABLO. GitHub Pages cannot start R itself.')}>?</span>
          </div>
          <div class="backend-controls">
          <select bind:value={referenceBackendMode} aria-label={t('Mode backend R', 'R backend mode')}>
              <option value="auto">{t('Auto · utiliser si disponible', 'Auto · use when available')}</option>
              <option value="browser">{t('Navigateur uniquement', 'Browser only')}</option>
              <option value="required">{t('R requis · bloquer si absent', 'Require R · fail if unavailable')}</option>
            </select>
          <input bind:value={referenceBackendUrl} aria-label={t('URL backend R', 'R backend URL')} />
          <button class="btn btn-outline btn-small" type="button" onclick={() => checkReferenceBackend(true)}>{t('Tester', 'Check')}</button>
          </div>
        <small class:backend-ok={referenceBackendStatus === 'available' || referenceBackendStatus === 'done'} class:backend-bad={referenceBackendStatus === 'error'}>
            {referenceBackendMessage || t('Par défaut : bridge local http://127.0.0.1:8787.', 'Default: local bridge http://127.0.0.1:8787.')}
          </small>


        <label class="inline-check" data-testid="diablo-holdout-option">
          <input type="checkbox" bind:checked={validateDiabloHoldout} />
            <span>{t('Vérifier DIABLO sur des sujets tenus à l’écart (analyse R supplémentaire)', 'Check DIABLO on held-out subjects (additional R analysis)')}</span>
          </label>
          <small>{t('Option avancée : entraînement et sélection de variables refaits dans chaque séparation. Requiert au moins 32 sujets indépendants, deux groupes, au moins 12 sujets/groupe et des matrices logarithmiques déjà préparées sans correction de batch. La métabolomique est actuellement exclue car son QC est fait avant les plis, même si les valeurs sont déjà logarithmiques. Les covariables et batches variables sont aussi refusés. Ce contrôle interne n’est pas une validation externe.', 'Advanced option: refit training and feature selection separately for every split. Requires ≥32 independent subjects, two classes, at least 12 per class and as-supplied log-scale matrices without batch correction. Metabolomics is currently excluded because MS QC precedes the splits, including for log-scale inputs. Varying batch and globally adjusted covariates are refused. Internal holdouts are not external validation.')}</small>
          {#if referenceBackendHealth?.packages}
            <div class="backend-packages">
              {#each Object.entries(referenceBackendHealth.packages) as [pkg, installed]}
                <span class:installed>{pkg}: {installed ? 'OK' : '—'}</span>
              {/each}
            </div>
          {/if}
        </div>
      </details>
    </div>
    <button class="btn btn-primary" type="button" data-testid="multiomics-run" disabled={!ready || analysisStatus === 'running' || msImportState === 'loading' || msImportState === 'parsed'} onclick={runAnalysis}>
      {analysisStatus === 'running' ? t('Analyse en cours…', 'Analyzing…') : t('Analyser mes données', 'Analyze my data')}
    </button>
  </div>
  {#if analysisError}<p class="error">{analysisError}</p>{/if}
</section>

{#if analysisResult}
<section class="panel demo-results" id="analysis-results" data-testid="multiomics-results">
    <div class="section-head">
      <div>
        <p class="eyebrow">
          {demoLoaded
            ? t('Démo synthétique · résultats réellement calculés', 'Synthetic demo · actually computed results')
            : t('Résultats · moteur déterministe', 'Analysis results · deterministic engine')}
        </p>

        <h2>
          {demoLoaded
            ? t('Résultats de la démonstration', 'Demonstration results')
            : t('Vos résultats', 'Your results')}
        </h2>

        <a
          href={resolve('multiomics/learn/')}
          class="simple-hint"
        >{t('Comprendre ces méthodes et figures', 'Learn how to interpret methods and figures')}</a>
      </div>
      <div class="result-actions">
      <button class="btn btn-primary" type="button" data-testid="multiomics-research-bundle-download" disabled={portableProjectBusy} onclick={downloadResearchBundle}>{portableProjectBusy ? t('Préparation du dossier…', 'Preparing dossier…') : t('Télécharger résultats + figures + code R (ZIP)', 'Download results + figures + R code (ZIP)')}</button>
      <button class="btn btn-outline" type="button" onclick={downloadReproducibleReport}>{t('Télécharger le rapport complet', 'Download full report')}</button>
        <details class="simple-disclosure export-detail">
          <summary>{t('Autres exports', 'Other downloads')}</summary>
        <div class="signature-optional-export" data-testid="multiomics-signature-optional">
            <label for="individual-signature-upload">{t('Option expert : joindre des scores biologiques individuels calculés et validés hors de MultiOmics (CSV)', 'Expert option: include independently calculated individual biological scores (CSV)')}</label>
          <input id="individual-signature-upload" data-testid="multiomics-signature-upload"
            type="file" accept=".csv,text/csv" onchange={importIndividualSignatures}/>
            <small>{t('Aucun score n’est inventé par MultiOmics. Assay, sujet, fenêtre de validation, modèle figé et SHA-256 doivent correspondre exactement. Une statistique de groupe ou un facteur MOFA ajusté sur les sujets test est refusé. Le ZIP ne prouve pas la validité du biomarqueur.', 'MultiOmics never invents scores. Assay, subject, fold, frozen model and SHA-256 must match exactly. Group statistics or transductive MOFA factors are rejected. ZIP export does not prove biomarker validity.')}</small>
          {#if individualSignatureFilename}<p role="status" data-testid="multiomics-signature-ready">{individualSignatureRows.length} {t('lignes lues','rows parsed')} — {individualSignatureFilename}. {t('Validation de provenance au téléchargement ZIP.', 'Provenance validation on ZIP download.')}</p>{/if}
          {#if individualSignatureError}<p role="alert" class="error" data-testid="multiomics-signature-error">{individualSignatureError}</p>{/if}
          </div>
        <button class="btn btn-outline" type="button" onclick={downloadAnalysisJson}>{t('Données détaillées (JSON)', 'Detailed data (JSON)')}</button>
        <button class="btn btn-outline" type="button" data-testid="multiomics-methods-audit-export" onclick={downloadResearchMethodsAudit}>{t('Méthodes et contrôles (JSON, sans données brutes)', 'Methods and checks (JSON, no raw data)')}</button>
        <button class="btn btn-outline" type="button" data-testid="multiomics-project-export" disabled={portableProjectBusy} onclick={exportPortableProject}>{portableProjectBusy ? t('Préparation…', 'Preparing…') : t('Enregistrer le projet complet', 'Save complete project')}</button>
          <small>{t('Contient les fichiers source, les réglages et les résultats. Peut inclure des données personnelles ou de recherche sensibles. Limite : 80 Mio de données source textuelles.', 'Includes original files, settings and results. May contain sensitive research or personal data. Limit: 80 MiB of source text.')}</small>
        {#if portableProjectMessage}<p role="status" class="ms-import-message">{portableProjectMessage}</p>{/if}
        {#if portableProjectError}<p role="alert" class="annotation-error">{portableProjectError}</p>{/if}
        </details>
        {#if analysisResult.reactome?.combined?.token}
        <a class="btn btn-outline" href={`https://reactome.org/PathwayBrowser/#DTAB=AN&ANALYSIS=${analysisResult.reactome.combined.token}`} target="_blank" rel="noreferrer">{t('Ouvrir dans Reactome ↗', 'Open in Reactome ↗')}</a>
        {/if}
      </div>
      <p class="simple-hint">{t('Le ZIP contient tous les tableaux, les méthodes, les sources et les paramètres. Données potentiellement sensibles : stockage sécurisé recommandé. Les résultats du navigateur restent exploratoires.', 'ZIP includes complete tables, methods, original data and settings. Potentially sensitive: store securely. Browser results remain exploratory.')}</p>
    <p class="simple-hint" data-testid="multiomics-replay-fidelity-note">
        {#if referenceReplay?.requestJson}
        {t('Le ZIP inclut la requête R originale et une copie des sources R, avec contrôle des empreintes SHA-256 si le backend les fournit. La concordance numérique exige un rejeu contrôlé. Le ZIP contient aussi un script R produisant des figures scientifiques similaires (heatmaps, métabologramme, carte métabolique) à partir des valeurs affichées, sans recalculer les statistiques. Les copies SVG du navigateur restent distinctes.',
           'The ZIP includes the original R request and R source snapshot, checked against the runtime SHA-256 fingerprints if provided. Numerical agreement still requires a verified replay. The ZIP also includes an R script generating scientifically similar heatmaps, circular effect diagrams and metabolic maps from the displayed data without refitting statistics. Browser SVG snapshots remain separate.')}
        {:else}
        {t('Cette analyse n’a pas de calcul R de référence réussi : le ZIP ne peut pas fournir un script R prétendument identique aux calculs JavaScript. Les figures sont des instantanés SVG.',
           'This analysis has no successful R reference computation: the ZIP cannot provide a purportedly identical R script for JavaScript calculations. Figures are SVG snapshots.')}
        {/if}
      </p>
      {#if analysisResult.referenceBackend?.methods}
      <div class="ms-import-message" role="note" data-testid="multiomics-reference-preview-warning">
          {t('Résultats R : le tableau affiché est un aperçu. Le ZIP ne contient les résultats R exhaustifs que pour les méthodes marquées « complet ». Toute méthode incomplète doit être exportée séparément depuis le moteur R local.', 'R results: on-screen tables are previews. The ZIP contains all R rows only for methods marked complete. Incomplete methods require a separate full export from the local R engine.')}
          <ul>
            {#each referencePreviewAudit(analysisResult) as ref (ref.method_key)}
            <li>{ref.reported_name} — {ref.complete_reference_table_exported
                  ? t('complet', 'complete')
                  : t('aperçu uniquement', 'preview only')}
              {#if ref.reference_table_rows !== null} ({ref.reference_table_rows} {t('lignes', 'rows')}){/if}
                {#if ref.reference_export_status === 'too_large_archived_locally'}
                <span> — {t('sauvegardé sur le disque du serveur R local, pas dans le ZIP', 'saved on the local R machine, not in this ZIP')}</span>
                  {#if ref.local_archive_id}
                    <code>{ref.local_archive_id}</code>

                    {#if ref.local_archive_sha256}
                      <small>SHA-256: {ref.local_archive_sha256}</small>
                    {/if}
                  {/if}
                {:else if ref.reference_export_status === 'too_large'}
                <span> — {t('trop volumineux : aucune archive conservée par défaut, reconfigurez R puis relancez', 'too large: no archive retained by default, configure R and rerun')}</span>
                {/if}
              </li>
            {/each}
          </ul>
        </div>
      {/if}
    {#if researchBundleMessage}<p role="status" data-testid="multiomics-research-bundle-status">{researchBundleMessage}</p>{/if}
    {#if researchBundleError}<p role="alert" data-testid="multiomics-research-bundle-error" class="error">{researchBundleError}</p>{/if}
    </div>

    {#if analysisResult.referenceBackend?.compatibilityStatus === 'mofa2_backend_outdated' || analysisResult.referenceBackend?.mofa2OrientationBlocked}
    <p role="alert" data-testid="multiomics-mofa-version-warning" class="annotation-error">
        {t('MOFA2 : résultat non valide avec cet ancien moteur R. Mettez à jour le backend R et relancez l’analyse. Les résultats exploratoires du navigateur restent disponibles.', 'MOFA2: results from this outdated R backend cannot be trusted. Update the local R backend and rerun. Browser exploration remains available.')}
      </p>
    {/if}
    {#if analysisResult.scientificAssurance}
    <div class="scientific-assurance" data-testid="multiomics-scientific-assurance">
        {#if analysisResult.scientificAssurance.confirmatoryReadiness?.requested}
          <div data-testid="multiomics-confirmatory-readiness">
          <strong>{analysisResult.scientificAssurance.confirmatoryReadiness.status === 'blocked'
                ? t('Analyse confirmatoire bloquée : conditions insuffisantes', 'Confirmatory analysis blocked: requirements not met')
            : t('Méthodes exécutées : revue scientifique indépendante encore requise', 'Methods executed: independent scientific review still required')}</strong>
            <p>{t('Cette vérification n’est pas une autorisation de publication ni une validation clinique. Corrigez les points bloquants ci-dessous avant d’interpréter les résultats comme confirmatoires.', 'This is not publication approval or clinical validation. Resolve the blockers below before interpreting the findings as confirmatory.')}</p>
            <details open>
              <summary>{t('Étapes de validation à vérifier', 'Validation requirements to review')}</summary>
            <ul>{#each analysisResult.scientificAssurance.confirmatoryReadiness.checks as check}
              <li>{check.status === 'blocked' ? t('À corriger : ', 'Fix: ') : check.status === 'review' ? t('À vérifier : ', 'Review: ') : t('Vérifié automatiquement : ', 'Automatically checked: ')}{$language === 'en' ? check.en : check.fr}</li>
            {/each}</ul>
            </details>
          </div>
        {/if}
        {#if analysisResult.protocol?.objective === 'time' && analysisResult.protocol?.longitudinal}
        <p data-testid="multiomics-longitudinal-inference-warning"><strong>{t('Longitudinal : pas de significativité calculée dans le navigateur.', 'Longitudinal: no browser-side statistical significance is reported.')}</strong>
            {t('Les changements et différences de pente sont descriptifs. Les p-values, q-values et intervalles de confiance exigent l’analyse de référence R avec lmerTest et un plan longitudinal vérifié.', 'Estimated changes and slopes are descriptive. P-values, q-values and confidence intervals require the lmerTest R reference analysis and a verified longitudinal design.')}
          </p>
        {/if}
      <strong>{t('Statut de validité scientifique', 'Scientific evidence status')} · {analysisResult.scientificAssurance.status === 'descriptive'
            ? t('Exploration descriptive', 'Descriptive exploration')
        : t('Inférence exploratoire : confirmation requise', 'Exploratory inference: confirmation required')}</strong>
      <p>{t('Exécution réussie ≠ validation scientifique. Une q-value faible ou un contrôle QC sans alerte ne certifient pas une découverte. Une conclusion publiable exige un design vérifié, une méthode de référence et, selon l’objectif, une validation indépendante.',
          'Successful execution is not scientific validation. Low q-values or apparently clean QC do not certify a finding. Publication claims require a checked design, appropriate reference methods and sometimes external validation.')}</p>
        {#if analysisResult.scientificAssurance.notes?.some(isMissingnessRiskNote)}
        <p role="status" class="annotation-error" data-testid="multiomics-missingness-caution">
          {t('Données manquantes détectées : elles peuvent créer de fausses différences, même si leur proportion est identique entre groupes. Les p-values et q-values ne suffisent pas : consultez les limites et prévoyez une analyse de sensibilité.',
             'Missing measurements detected: they can create false differences even with identical missing-data rates between groups. P-values and q-values alone are insufficient: review the limitations and plan a sensitivity analysis.')}
          </p>
        {/if}
        <details>
          <summary>{t('Voir les limites spécifiques à cette analyse', 'See limits specific to this analysis')}</summary>
          <ul>
            {#each analysisResult.scientificAssurance.notes as note}
              <li>{$language === 'en' ? note.en : note.fr}</li>
            {/each}
          </ul>
        </details>
      </div>
    {/if}

    <div class="computed-summary">
    <article data-testid="multiomics-analysis-mode"><span>{t('Type d’analyse', 'Analysis type')}</span><strong>{analysisResult.engine?.analysisMode === 'single_omic' ? t('Omique unique', 'Single omic') : t('Multi-omique', 'Multi-omics')}</strong></article>
    <article><span>{t('Sujets', 'Subjects')}</span><strong>{analysisResult.metadataSummary.subjects}</strong></article>
    <article><span>{t('Prélèvements biologiques', 'Biological samples')}</span><strong>{analysisResult.metadataSummary.samples}</strong></article>
    <article><span>{t('Mesures', 'Assays')}</span><strong>{analysisResult.metadataSummary.assays}</strong></article>
    <article><span>{t('Conditions', 'Conditions')}</span><strong>{analysisResult.metadataSummary.conditions.join(' / ') || '—'}</strong></article>
    <article><span>{t('Temps', 'Time points')}</span><strong>{analysisResult.metadataSummary.timepoints.join(' / ') || '—'}</strong></article>
    </div>

    {#if analysisResult.visualizations}
    <MultiomicsVisualizations result={analysisResult} language={$language ?? 'fr'} />
    {:else if analysisResult.visualizationError}
      <p class="note">{t('Les résultats statistiques sont disponibles, mais la préparation des visualisations a échoué : ', 'Statistical results are available, but visualization preparation failed: ')}{analysisResult.visualizationError}</p>
    {/if}

  <div class="simple-qc-status" data-testid="multiomics-quality-summary">
      <strong>{t('Méthodes réellement exécutées', 'Methods actually run')}</strong>
      <span>{t('Navigateur : ', 'Browser: ')}{Object.entries(analysisResult.layers || {}).map(([key, value]) => omicLabel(key) + ' — ' + (value.inferenceMethod || value.mode || 'QC')).join(' ; ')}</span>
    <small data-testid="multiomics-numeric-source">{analysisResult.displayedReferenceMethodCount > 0
      ? t('Les tableaux des couches portant « Native R » proviennent des résultats R exhaustifs ; les figures JS utilisent ces mêmes effets, sans refaire les tests. Les contrôles QC et les autres méthodes exploratoires restent distincts.',
          'Layers labeled Native R show full R reference values; JavaScript figures use these same effects without refitting. QC and other exploratory calculations remain distinct.')
      : t('Chiffres exploratoires JS : ils peuvent reposer sur un modèle simplifié ou des approximations et différer des p-values, q-values et effets obtenus en R. Ne pas les utiliser comme inférence confirmatoire. Le moteur R de référence applique les méthodes statistiques natives appropriées, sans garantir à lui seul une conclusion valide.',
          'Exploratory JS numbers may use simplified models or approximations and differ in p-values, adjusted p-values and effects from native R methods. They are not confirmatory inference. The reference R engine applies the appropriate native statistical methods, which alone do not guarantee valid conclusions.')}</small>
      {#if analysisResult.referenceBackend?.status === 'ok'}
        {#each Object.entries(analysisResult.referenceBackend.methods || {}) as [key, value]}
          <span>{key} : {value.method || key} — {value.status}</span>
        {/each}
      {:else if analysisResult.referenceBackend}
        <span>{t('R non exécuté : ', 'R not executed: ')}{analysisResult.referenceBackend.message || analysisResult.referenceBackend.status}</span>
      {:else}
        <span>{t('R non demandé pour cette analyse (ex. démonstration navigateur).', 'R not requested for this analysis (e.g. browser demo).')}</span>
      {/if}
      <strong>{t('Contrôle des données', 'Data quality check')}</strong>
      {#each Object.entries(analysisResult.layers || {}) as [layer, layerResult]}
      <span class:qc-alert={layerResult.qc?.warnings?.length}>{omicLabel(layer)} : {layerResult.qc?.warnings?.length ? t('points à vérifier', 'needs review') : layerResult.qc?.inferenceTier?.level === 'screening' ? t('exploration uniquement', 'exploratory only') : t('aucune alerte détectée', 'no warning detected')}</span>
      {/each}
    </div>
  <details class="advanced-results" data-testid="multiomics-quality-details">
      <summary>{t('Voir les vérifications détaillées', 'View detailed quality checks')}</summary>

      {#if analysisResult.metadataSummary.overlap?.pairwise?.length}
        <div class="overlap-box">
          <div>
            <p class="eyebrow">{t('Chevauchement réellement utilisé', 'Sample overlap actually used')}</p>
            <strong>{analysisResult.metadataSummary.overlap.allMatched} {t('sujet(s) présents dans toutes les couches chargées', 'subject(s) present in every loaded omics layer')}</strong>
          </div>
          <div class="overlap-pairs">
            {#each analysisResult.metadataSummary.overlap.pairwise as pair}
          <span>{omicLabel(pair.layerA)} ↔ {omicLabel(pair.layerB)}: <b>{pair.matchedSubjects}</b> matched</span>
            {/each}
          </div>
        </div>
      {/if}

      {#if analysisResult.metadataSummary.batchAudit}
        <div class="overlap-box">
          <div>
            <p class="eyebrow">{t('Audit des batches techniques', 'Technical batch audit')}</p>
            <strong>{t('Structure des batches vérifiée avant l’inférence biologique', 'Batch structure checked before biological inference')}</strong>
          </div>
          <div class="overlap-pairs">
            {#each omicLayers as layer}
              {#if analysisResult.metadataSummary.batchAudit[layer]}
                {@const audit = analysisResult.metadataSummary.batchAudit[layer]}
            <span>{omicLabel(layer)}: <b>{audit.status.replaceAll('_', ' ')}</b>{audit.batches?.length ? ' · ' + audit.batches.length + ' batch(es)' : ''}</span>
              {/if}
            {/each}
          </div>
        </div>
      {/if}

      {#if analysisResult.referenceBackend}
    <div class="integration-result reference-backend-result" data-testid="reference-backend-results">
          <div class="integration-head">
            <div>
              <p class="eyebrow">{t('Moteur R de référence', 'Reference R engine')}</p>
          <h3>{analysisResult.referenceBackend.status === 'ok'
                  ? analysisResult.referenceBackend.referenceMethodExecuted
                    ? t('Méthodes R exécutées (résultats distincts des figures du navigateur)', 'R methods executed (results separate from browser figures)')
                    : t('Moteur R connecté — aucune méthode statistique de référence exécutée', 'R backend connected — no reference statistical method executed')
            : t('Backend R non utilisé', 'R backend not used')}</h3>
            </div>
            <span>{analysisResult.referenceBackend.status}</span>
          </div>
          {#if analysisResult.referenceBackend.status === 'ok'}
            <p class="note">{analysisResult.referenceBackend.engine?.name} · v{analysisResult.referenceBackend.engine?.version}</p>
            {#if analysisResult.referenceBackend.runtime}
          <details class="simple-disclosure" data-testid="multiomics-r-runtime-provenance">
                <summary>{t('Versions R et packages utilisés (reproductibilité)', 'R and package versions used (reproducibility)')}</summary>
            <p>{t('Version R : ', 'R version: ')}{analysisResult.referenceBackend.runtime.RVersion} ·
              {analysisResult.referenceBackend.runtime.RPlatform}
                  {#if analysisResult.referenceBackend.runtime.BioconductorRelease}
                    · Bioconductor {analysisResult.referenceBackend.runtime.BioconductorRelease}
                  {/if}
                </p>
                <div class="backend-packages">
                  {#each Object.entries(analysisResult.referenceBackend.runtime.packageVersions || {}) as [pkg, version]}
                    <span>{pkg} {version}</span>
                  {/each}
                </div>
                <small>{t('Ces versions sont conservées dans le rapport et le projet téléchargé, mais ne remplacent pas un environnement R figé (renv/conteneur).', 'Versions are included in the report and exported project, but do not replace an R environment lock (renv/container).')}</small>
              </details>
            {/if}
            <div class="backend-method-grid">
              {#each Object.entries(analysisResult.referenceBackend.methods || {}) as [key, method]}
                <article>
                  <strong>{method.method || key}</strong>
                  <span class:installed={method.status === 'ok'}>{method.status}</span>
              {#if method.message}<small>{method.message}</small>{/if}
                  {#if key.endsWith('_longitudinal') && method.summary}
                    <small data-testid="multiomics-r-longitudinal-summary">
                   {t('Contraste :', 'Contrast:')} {method.summary.contrast} ·
                   {t('Sujets :', 'Subjects:')} {method.summary.nSubjects} ·
                   {t('Variables estimées :', 'Fitted features:')} {method.summary.estimableFeatures}/{method.summary.attemptedFeatures} ·
                   {t('Variables non estimables :', 'Non-estimable features:')} {method.summary.blockedFeatures} ·
                   {t('Réductions de pente aléatoire :', 'Random-slope fallbacks:')} {method.summary.randomSlopeFallbacks}
                    </small>
                    <small>{t('Contraste linéaire à deux groupes fixes uniquement. Autres trajectoires ou visites manquantes : méthode dédiée nécessaire.', 'Two fixed groups and linear slope contrast only. Other trajectories or missing visits need dedicated methods.')}</small>
                  {/if}
                  {#if key === 'diablo_heldout' && method.summary}
                    <small data-testid="diablo-holdout-summary">
                   {t('Erreur équilibrée sur les sujets tenus à l’écart :', 'Balanced error on held-out subjects:')}
                   {(method.summary.meanBER * 100).toFixed(1)} %
                   · {method.summary.folds?.length ?? 0} {t('séparations', 'splits')}
                    </small>
                    {#if method.summary.signatureStability}
                   <details class="simple-disclosure" data-testid="diablo-signature-stability">
                        <summary>{t('Les mêmes variables reviennent-elles selon les séparations ?', 'Do the same features recur across splits?')}</summary>
                        {#each Object.entries(method.summary.signatureStability) as [view, stability]}
                       <small>{view} · Jaccard {stability.meanPairwiseJaccard.toFixed(2)}
                         · {stability.repeatedlySelected?.length ?? 0} {t('variables retenues au moins deux fois', 'features selected at least twice')}</small>
                        {/each}
                        <small>{t('Jaccard : de 0 (aucune variable commune) à 1 (exactement les mêmes variables). Cette répétabilité interne ne valide pas un biomarqueur.', 'Jaccard: 0 means no overlapping features, 1 means identical selections. This internal repeatability does not validate biomarkers.')}</small>
                      </details>
                    {/if}
                    <small>{t('Évaluation interne : ces sujets restent issus de la même cohorte. Une validation externe indépendante est toujours nécessaire.', 'Internal evaluation: subjects come from the same original cohort. External validation on an independent cohort is still required.')}</small>
                  {/if}
                  {#if key === 'mofa2' && method.summary?.inputAudit}
                  <small data-testid="mofa2-input-audit">{t('MOFA2 :', 'MOFA2:')} {method.summary.inputAudit.nSharedSubjects}
                    {t('sujets partagés. Un facteur latent résume des variations entre molécules. Il peut refléter la biologie ou un effet de lot (batch), sans prouver une causalité ni un biomarqueur.', 'shared subjects. A latent factor summarizes patterns among measured molecules. It may reflect biology or a technical batch effect, not causal proof or a validated biomarker.')}</small>
                  {/if}
                  {#if key === 'metabolomics_ms_qc' && method.summary}
                    <small class="backend-ms-detail">
                  {method.summary.blank_injections ?? 0} blanks ·
                  {method.summary.qc_injections ?? 0} pooled-QC ·
                  {method.summary.blank_flagged_features ?? method.summary.blank_filtered_features ?? 0} {t('signalées', 'flagged')} ·
                  {method.summary.blank_removed_features ?? method.summary.blank_filtered_features ?? 0} {t('exclues', 'removed')} ·
                  {method.summary.drift_corrected_features ?? 0} {t('corrigées dérive', 'drift-corrected')} ·
                  {method.summary.qc_rsd_filtered_features ?? 0} {t('exclues RSD', 'RSD-removed')} ·
                  {method.summary.mnar_imputed_values ?? 0} {t('imputées MNAR', 'MNAR-imputed')}
                    </small>
                  {/if}
                </article>
              {/each}
            </div>
          {:else}
            <p class="muted">{analysisResult.referenceBackend.message}</p>
          {/if}
        </div>
      {/if}

  <div class="integration-result qc-result" data-testid="multiomics-qc">
        <div class="integration-head">
          <div>
            <p class="eyebrow">{t('Contrôle qualité', 'Quality control')}</p>
            <h3>{t('Qualité des données avant interprétation biologique', 'Data quality before biological interpretation')}</h3>
          </div>
          <span>{t('QC par couche', 'per-layer QC')}</span>
        </div>

        <div class="qc-grid">
          {#each Object.entries(analysisResult.layers) as [layer, result]}
            {#if result.qc}
              <article class="qc-card">
                <div class="qc-card-head">
                  <div>
                    <strong>{omicLabel(layer)}</strong>
                    <small>{result.qc.featuresAfter}/{result.qc.featuresBefore} {t('variables conservées', 'features retained')}</small>
                  </div>
                  <div class="qc-badges">
                    <span class:qc-warning={result.qc.warnings?.length}>{result.qc.warnings?.length ? t('à vérifier', 'review') : 'OK'}</span>
                    {#if result.qc.inferenceTier?.level === 'screening'}
                  <span class="qc-screening">{t('screening', 'screening')} <span class="help-tip" tabindex="0" data-tooltip={t('Pour des counts RNA-seq bruts, la branche navigateur sert au screening reproductible. Pour une inférence différentielle publication-grade, utilisez l’adaptateur R DESeq2/limma avec le même design.', 'For raw RNA-seq counts, the browser branch is intended for reproducible screening. For publication-grade differential inference, use the DESeq2/limma R adapter with the same design.')}>?</span></span>
                    {/if}
                  </div>
                </div>

                <div class="qc-metrics">
                  <div>
                    <b>
                      {Number.isFinite(result.qc.medianMissingFraction)
                        ? (100 * result.qc.medianMissingFraction).toFixed(1) + '%'
                        : '—'}
                    </b>

                    <span>{t('missing médian', 'median missing')}</span>
                  </div>

                  <div>
                    <b>{result.qc.outlierSamples?.length || 0}</b>
                    <span>{t('assays suspects', 'flagged assays')}</span>
                  </div>

                  <div>
                    <b>{poorReplicateCount(result.qc)}</b>
                    <span>{t('réplicats r<0,80', 'replicates r<0.80')}</span>
                  </div>

                  <div>
                    <b>{result.qc.differentialMissingness?.flaggedFeatures || 0}</b>
                    <span>{t('variables à missing différentiel', 'features with differential missingness')}</span>
                  </div>
                </div>

                <div
                  class="qc-bars"
                  aria-label={t('Missingness par assay', 'Missingness by assay')}
                >
                  {#each result.qc.sampleMetrics.slice(0, 16) as metric}
                <div title={metric.assayId + ' · missing=' + (metric.missingFraction == null ? 'NA' : (100 * metric.missingFraction).toFixed(1) + '%') + ' · detected=' + metric.detectedFeatures}>
                      <span>{metric.assayId}</span>
                  <i><em style={'width:' + Math.max(1, Math.min(100, 100 * (metric.missingFraction || 0))) + '%'}></em></i>
                    </div>
                  {/each}
                </div>

                {#if result.qc.pca?.scores?.length}
                  <div class="qc-pca">
                    <div class="qc-pca-head">
                  <b>PCA QC <span class="help-tip" tabindex="0" data-tooltip={t('Projection non supervisée utilisée uniquement pour contrôler la structure globale des échantillons après prétraitement : outliers, batch éventuel et séparation dominante. Elle ne constitue pas un test statistique.', 'Unsupervised projection used only to inspect the global sample structure after preprocessing: outliers, possible batch structure and dominant separation. It is not a statistical test.')}>?</span></b>
                      <small>
                    PC1 {Number.isFinite(result.qc.pca.explained?.[0]) ? (100 * result.qc.pca.explained[0]).toFixed(1) + '%' : '—'}
                    · PC2 {Number.isFinite(result.qc.pca.explained?.[1]) ? (100 * result.qc.pca.explained[1]).toFixed(1) + '%' : '—'}
                      </small>
                    </div>
                <div class="qc-scatter" role="img" aria-label={t('PCA des échantillons après prétraitement', 'PCA of samples after preprocessing')}>
                  <span class="axis-x"></span><span class="axis-y"></span>
                      {#each result.qc.pca.scores as point}
                        <button
                          type="button"
                          class="qc-dot"
                          style={'left:' + qcAxisPercent(point.pc1, result.qc.pca.scores, 'pc1') + '%;bottom:' + qcAxisPercent(point.pc2, result.qc.pca.scores, 'pc2') + '%'}
                          title={point.sampleId + ' · ' + (point.condition || '—') + ' · batch ' + (point.batch || '—')}
                          aria-label={point.sampleId + ', ' + (point.condition || '') + ', ' + (point.batch || '')}
                        ></button>
                      {/each}
                    </div>
                    <p class="muted">{result.qc.pca.method}</p>
                  </div>
                {/if}

                {#if result.qc.msQc?.applied}
                  <div class="ms-qc-summary" data-testid="ms-qc-summary">
                    <strong>{t('QC MS', 'MS QC')}</strong>
                    <span>{result.qc.msQc.blankAssays} blanks</span>
                    <span>{result.qc.msQc.qcAssays} pooled-QC</span>
                    <span>{result.qc.msQc.blankFlaggedFeatures ?? result.qc.msQc.blankFilteredFeatures ?? 0} {t('signalées par blank', 'blank-flagged')}</span>
                    <span>{result.qc.msQc.blankRemovedFeatures ?? result.qc.msQc.blankFilteredFeatures ?? 0} {t('exclues par blank', 'blank-removed')}</span>
                    <span>{result.qc.msQc.qcRsdFilteredFeatures} {t('retirées par RSD', 'RSD-filtered')}</span>
                    <span>{result.qc.msQc.driftCorrection?.correctedFeatures || 0} {t('corrigées pour dérive', 'drift-corrected')}</span>
                    <span>{result.qc.msQc.mnar?.imputedValues || 0} {t('MNAR imputées', 'MNAR-imputed')}</span>
                  </div>
                {/if}

                {#if result.qc.preprocessingAudit}
              <div class="api-summary" data-testid="preprocessing-audit">
                <span><strong>{t('Type déclaré', 'Declared type')}</strong> {result.qc.preprocessingAudit.declaredValueType}</span>
                <span><strong>{t('Échelle analysée', 'Analysis scale')}</strong> {result.qc.preprocessingAudit.outputScale}</span>
                <span><strong>{t('Compatibilité', 'Compatibility')}</strong> {result.qc.preprocessingAudit.status}</span>
                <span><strong>{t('Plage entrée', 'Input range')}</strong> {result.qc.preprocessingAudit.inputRange?.min ?? '—'} → {result.qc.preprocessingAudit.inputRange?.max ?? '—'}</span>
                  </div>
                {/if}

                <details>
                  <summary>{t('Règles QC et prétraitement', 'QC and preprocessing rules')}</summary>
                  <p>{result.qc.filterPolicy}</p>
              {#if result.qc.preprocessingAudit}<p class="note">{result.qc.preprocessingAudit.policy}</p>{/if}
              {#if result.qc.inferenceTier?.note}<p class="note">{result.qc.inferenceTier.note}</p>{/if}
                  <ul>
                {#each result.qc.preprocessingSteps || [] as step}<li>{step}</li>{/each}
                {#each result.qc.warnings || [] as warning}<li class="warning">{warning}</li>{/each}
                {#each result.qc.msQc?.warnings || [] as warning}<li class="warning">{warning}</li>{/each}
                  </ul>
                </details>
              </article>
            {/if}
          {/each}
        </div>
      </div>
    </details>

  <div class="interpretation-box" data-testid="multiomics-interpretation">
      <div class="integration-head">
        <div>
          <p class="eyebrow">{t('Aide à l’interprétation', 'Interpretation guide')}</p>
          <h3>{t('Comment interpréter ces résultats ?', 'How should these results be interpreted?')}</h3>
        </div>
        <span>{t('règles déterministes', 'deterministic rules')}</span>
      </div>
      <div class="interpretation-grid">
        {#each interpretationItems as item}
        <article>
          <strong>{item.title}</strong>
          <p>{item.text}</p>
        </article>
        {/each}
      </div>
    </div>

    <details class="dictionary interpretation-glossary">
      <summary>{t('Glossaire pour lire les sorties statistiques', 'Glossary for reading statistical outputs')}</summary>
      <div class="dictionary-table">
        <div>
          <b>{t('Terme', 'Term')}</b>
          <b>{t('Lecture', 'Meaning')}</b>
          <b>{t('À retenir', 'Key point')}</b>
          <b>{t('Piège à éviter', 'Common pitfall')}</b>
        </div>

        <div>
          <code>effect / β</code>
          <span>{t('Direction et amplitude du contraste ou coefficient de régression.', 'Direction and magnitude of the contrast or regression coefficient.')}</span>
          <span>{t('Toujours lire son échelle et le contraste exact.', 'Always read its scale and exact contrast.')}</span>
          <span>{t('Un grand effet n’implique pas automatiquement une faible q-value.', 'A large effect does not automatically imply a low q-value.')}</span>
        </div>

        <div>
          <code>fold ratio</code>
          <span>{t('Pour une échelle log2 : 2^effect.', 'On a log2 scale: 2^effect.')}</span>
          <span>{t('>1 : hausse ; <1 : baisse selon le contraste affiché.', '>1: increase; <1: decrease according to the displayed contrast.')}</span>
          <span>{t('Ne pas l’interpréter si les valeurs ne sont pas sur une échelle log2 déclarée.', 'Do not interpret it if values are not on a declared log2 scale.')}</span>
        </div>

        <div>
          <code>p</code>
          <span>{t('Évidence du test pour une variable avant correction multiple.', 'Test evidence for one feature before multiple-testing correction.')}</span>
          <span>{t('Utile pour le calcul de q, pas pour sélectionner seule des milliers de variables.', 'Useful for computing q, not for selecting among thousands of features by itself.')}</span>
          <span>{t('p < 0,05 n’est pas une probabilité que l’hypothèse soit vraie.', 'p < 0.05 is not the probability that a hypothesis is true.')}</span>
        </div>

        <div>
          <code>q BH</code>
          <span>{t('p-value ajustée par Benjamini–Hochberg pour contrôler le FDR.', 'Benjamini–Hochberg adjusted p-value for FDR control.')}</span>
          <span>{t('Le prototype utilise q ≤ 0,10 comme seuil exploratoire.', 'The prototype uses q ≤ 0.10 as an exploratory threshold.')}</span>
          <span>{t('Ce seuil ne remplace pas la taille d’effet ni la validation externe.', 'This threshold does not replace effect size or external validation.')}</span>
        </div>

        <div>
          <code>OR / RR / HR</code>
          <span>{t('exp(β) en logistique, Poisson ou Cox.', 'exp(β) in logistic, Poisson or Cox models.')}</span>
          <span>{t('1 = absence d’association ; >1 et <1 donnent les deux directions.', '1 = no association; >1 and <1 indicate opposite directions.')}</span>
          <span>{t('Une association n’est pas une causalité.', 'Association is not causation.')}</span>
        </div>

        <div>
          <code>95% CI</code>
          <span>{t('Intervalle de confiance Wald autour de l’effet estimé.', 'Wald confidence interval around the estimated effect.')}</span>
          <span>{t('La largeur renseigne sur la précision.', 'Width reflects estimation precision.')}</span>
          <span>{t('Avec petits effectifs ou séparation logistique, l’approximation peut être fragile.', 'With small samples or logistic separation, the approximation can be fragile.')}</span>
        </div>

        <div>
          <code>loading</code>
          <span>{t('Contribution d’une variable à un axe latent de l’ACP multi-blocs.', 'Feature contribution to a latent balanced multi-block PCA axis.')}</span>
          <span>{t('Interpréter surtout |loading| et les variables qui covarient sur un axe.', 'Focus on |loading| and features co-varying along an axis.')}</span>
          <span>{t('Le signe global d’un axe PCA est arbitraire.', 'The global sign of a PCA axis is arbitrary.')}</span>
        </div>

        <div>
          <code>Δr</code>
          <span>{t('Différence de corrélation inter-omique entre deux conditions.', 'Difference in cross-omics correlation between two conditions.')}</span>
          <span>{t('Peut révéler un couplage gagné, perdu ou inversé.', 'Can reveal gained, lost or reversed coupling.')}</span>
          <span>{t('Une corrélation condition-dépendante ne prouve pas une interaction causale.', 'Condition-dependent correlation does not prove a causal interaction.')}</span>
        </div>
      </div>
    </details>

  <details class="advanced-results" data-testid="multiomics-detailed-results">
      <summary>{t('Tableaux et statistiques détaillés', 'Detailed tables and statistics')}</summary>
      {#if analysisResult.exploration?.components?.length}
        <div class="integration-result">
          <div class="integration-head">
            <div>
          <p class="eyebrow">{analysisResult.engine?.analysisMode === 'single_omic' ? t('Exploration d’une omique', 'Single-omics exploration') : t('Exploration multi-omique', 'Multi-omics exploration')}</p>
          <h3>{analysisResult.engine?.analysisMode === 'single_omic' ? t('ACP : variation au sein de cette omique', 'PCA: variation within this omics layer') : t('Axes latents partagés entre les couches', 'Shared latent axes across omics layers')}</h3>
            </div>
        <span>{analysisResult.exploration.subjects} {analysisResult.engine?.analysisMode === 'single_omic' ? t('sujets', 'subjects') : t('sujets communs', 'shared subjects')}</span>
          </div>
          <div class="api-summary">
            {#each analysisResult.exploration.components as component}
          <span><strong>PC{component.component}</strong> {Number.isFinite(component.explainedFraction) ? (100 * component.explainedFraction).toFixed(1) + '%' : '—'} {t('de variance pondérée', 'balanced variance')}</span>
            {/each}
          </div>
          <p class="note">{analysisResult.exploration.method}</p>
          <div class="cross-table">
        <div class="cross-head"><b>{t('Variable', 'Feature')}</b><b>{t('Omique', 'Omics')}</b><b>PC1 loading</b><b>{t('Contribution', 'Contribution')}</b><b></b><b></b></div>
            {#each analysisResult.exploration.components[0].topLoadings.slice(0, 15) as loading}
              <div>
                <code>{loading.feature}</code>
                <span>{omicLabel(loading.layer)}</span>
                <strong>{loading.loading.toPrecision(3)}</strong>
                <span>{Math.abs(loading.loading).toPrecision(3)}</span>
            <span></span><span></span>
              </div>
            {/each}
          </div>
        </div>
      {/if}

      {#if analysisResult.supervisedIntegration}
        <div class="integration-result">
          <div class="integration-head">
            <div>
              <p class="eyebrow">{t('Intégration supervisée', 'Supervised integration')}</p>
          <h3>{t('Composante latente multiblocs liée à la cible', 'Target-linked multiblock latent component')}
  <span class="help-tip" tabindex="0" data-tooltip={t('Résumé supervisé descriptif : chaque couche est standardisée et équilibrée, puis une composante maximise la covariance avec la cible. Les poids indiquent les variables qui structurent cet axe. Ce n’est pas une validation prédictive.', 'Descriptive supervised summary: each block is standardized and balanced, then one component maximizes covariance with the target. Weights indicate features structuring this axis. This is not predictive validation.')}>?</span>
              </h3>
            </div>
            <span>{analysisResult.supervisedIntegration.subjects} {t('sujets', 'subjects')}</span>
          </div>
          <div class="api-summary">
        <span><strong>{Number.isFinite(analysisResult.supervisedIntegration.scoreTargetCorrelation) ? analysisResult.supervisedIntegration.scoreTargetCorrelation.toFixed(2) : '—'}</strong> r score-cible</span>
        <span><strong>{Number.isFinite(analysisResult.supervisedIntegration.scoreTargetR2) ? (100 * analysisResult.supervisedIntegration.scoreTargetR2).toFixed(1) + '%' : '—'}</strong> R² descriptif</span>
          </div>
          <p class="note">{analysisResult.supervisedIntegration.method}</p>
          <div class="cross-table">
        <div class="cross-head"><b>{t('Variable', 'Feature')}</b><b>{t('Omique', 'Omics')}</b><b>{t('Poids', 'Weight')}</b><b>|weight|</b><b></b><b></b></div>
            {#each analysisResult.supervisedIntegration.topWeights.slice(0, 15) as item}
              <div>
                <code>{item.feature}</code>
                <span>{omicLabel(item.layer)}</span>
                <strong>{item.weight.toPrecision(3)}</strong>
                <span>{item.absoluteWeight.toPrecision(3)}</span>
            <span></span><span></span>
              </div>
            {/each}
          </div>
          <p class="note">{analysisResult.supervisedIntegration.caveat}</p>
        </div>
      {/if}

      {#if analysisResult.predictiveOutcome}
        <div class="integration-result predictive-result">
          <div class="integration-head">
            <div>
              <p class="eyebrow">{t('Validation prédictive', 'Predictive validation')}</p>
          <h3>{t('Performance hors échantillon', 'Out-of-sample performance')}
  <span class="help-tip" tabindex="0" data-tooltip={t('La nested cross-validation sépare des folds externes pour mesurer la performance et des folds internes pour choisir la pénalisation. Ajustement batch/covariables, sélection de variables, centrage et scaling sont estimés uniquement sur les sujets d’entraînement puis appliqués au fold testé.', 'Nested cross-validation uses outer folds for performance assessment and inner folds for penalty tuning. Batch/covariate adjustment, feature selection, centering and scaling are estimated from training subjects only and then applied to the held-out fold.')}>?</span>
              </h3>
            </div>
        <span>{analysisResult.predictiveOutcome.status === 'ok' ? t('nested CV', 'nested CV') : t('non estimable', 'not estimable')}</span>
          </div>
          {#if analysisResult.predictiveOutcome.status === 'ok'}
            <div class="api-summary">
          {#if analysisResult.predictiveOutcome.metrics.auc != null}<span><strong>{analysisResult.predictiveOutcome.metrics.auc.toFixed(3)}</strong> AUC</span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.accuracy != null}<span><strong>{(100 * analysisResult.predictiveOutcome.metrics.accuracy).toFixed(1)}%</strong> accuracy</span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.rmse != null}<span><strong>{analysisResult.predictiveOutcome.metrics.rmse.toPrecision(3)}</strong> RMSE</span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.r2 != null}<span><strong>{analysisResult.predictiveOutcome.metrics.r2.toFixed(3)}</strong> R² CV</span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.brier != null}<span data-testid="multiomics-prediction-brier"><strong>{analysisResult.predictiveOutcome.metrics.brier.toFixed(3)}</strong> Brier <span class="help-tip" tabindex="0" data-tooltip={t('Moyenne des erreurs quadratiques entre probabilités prédites et issue binaire observée, calculée uniquement sur les prédictions des folds externes. Plus petit est meilleur, mais comparez à une référence et contrôlez la calibration sur une cohorte externe.', 'Mean squared error between held-out predicted probabilities and observed binary outcomes; lower is better, but it must be compared against an appropriate reference and externally checked for calibration.')}>?</span></span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.logLoss != null}<span><strong>{analysisResult.predictiveOutcome.metrics.logLoss.toFixed(3)}</strong> log-loss</span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.cIndex != null}<span><strong>{analysisResult.predictiveOutcome.metrics.cIndex.toFixed(3)}</strong> C-index <span class="help-tip" tabindex="0" data-tooltip={t('Indice de concordance de Harrell calculé uniquement sur les prédictions des folds externes. 0,5 correspond approximativement au hasard ; plus il est proche de 1, mieux le score de risque ordonne les temps d’événement observables.', 'Harrell concordance index computed only from outer-fold predictions. About 0.5 corresponds to chance; values closer to 1 mean the risk score better orders observable event times.')}>?</span></span>{/if}
          {#if analysisResult.predictiveOutcome.metrics.events != null}<span><strong>{analysisResult.predictiveOutcome.metrics.events}</strong> {t('événements', 'events')}</span>{/if}
            </div>
            {#if analysisResult.predictiveOutcome.metrics.auc != null}
          <p class="note" role="note" data-testid="multiomics-binary-calibration-caution">
                <strong>{t('Attention : une AUC élevée ne prouve pas la calibration des probabilités.', 'Caution: a high AUC does not establish well-calibrated probabilities.')}</strong>
                {t('Le score de Brier et la log-loss sont calculés sur les prédictions de validation croisée, mais leur interprétation exige une comparaison à un modèle de référence, une vérification des courbes de calibration et une cohorte externe indépendante. Une recalibration peut dégrader les performances sur une autre cohorte.', 'Brier and log-loss are computed on held-out cross-validation predictions. Interpretation requires a reference predictor, calibration plots and an independent external cohort; recalibration can degrade performance on another cohort.')}
              </p>
            {/if}
            <p class="note">{analysisResult.predictiveOutcome.method}. {analysisResult.predictiveOutcome.caveat}</p>
            {#if analysisResult.predictiveOutcome.rawIntensityScreening}
          <p class="note" data-testid="multiomics-ms-prediction-caveat">
                <strong>{t('Prédiction LFQ/MS : contrôle technique encore nécessaire.', 'LFQ/MS prediction: technical QC review still required.')}</strong>
                {t('Pour éviter la fuite entre sujets, le modèle utilise une transformation fixée et un centrage par échantillon. Il ne reprend pas automatiquement les corrections pooled-QC, blancs et dérive MS des analyses descriptives. Le même panel de variables doit être mesuré lors de la validation externe.', 'To avoid cross-subject leakage, prediction uses a fixed transform and assay-local centering. It does not reuse pooled-QC, blank and MS drift corrections from descriptive analyses. External validation must measure the same feature panel.')}
              </p>
            {/if}
            <details>
              <summary>{t('Détails des folds', 'Fold details')}</summary>
              <div class="fold-grid">
                {#each analysisResult.predictiveOutcome.foldSummaries as fold}
                  <span>Fold {fold.fold}: n train={fold.trainingSubjects}, n test={fold.testSubjects}, λ={fold.lambda}, p={fold.selectedFeatures}</span>
                {/each}
              </div>
            </details>
          {:else}
            <p class="muted">{analysisResult.predictiveOutcome.reason}</p>
          {/if}
        </div>
      {/if}

      {#if analysisResult.identifierResolution?.metabolomics}
        <div class="identifier-resolution">
          <div class="integration-head">
            <div>
              <p class="eyebrow">{t('Résolution des identifiants métabolites', 'Metabolite identifier resolution')}</p>
              <h3>{t('Libellés originaux → identifiants canoniques envoyés à l’analyse de voies', 'Original labels → canonical identifiers sent to pathway analysis')}</h3>
            </div>
            <span>{analysisResult.identifierResolution.metabolomics.resolvedCount} resolved · {analysisResult.identifierResolution.metabolomics.unresolvedCount} unresolved</span>
          </div>
          <div class="mapping-result-table">
        <div class="mapping-result-head"><b>Original</b><b>Query</b><b>Resolved ID</b><b>Status</b></div>
            {#each analysisResult.identifierResolution.metabolomics.mappings.slice(0, 20) as mapping}
              <div>
                <code>{mapping.original}</code>
                <span>{mapping.query || '—'}</span>
                <code>{mapping.resolved || '—'}</code>
                <span>{mapping.status}</span>
              </div>
            {/each}
          </div>
          <p class="note">{t('Seules les correspondances exactes ChEBI ou les alias lipidiques déterministes explicites sont acceptés automatiquement. Une correspondance incertaine reste non résolue au lieu d’être devinée.', 'Only exact ChEBI label matches or explicit deterministic lipid aliases are accepted automatically. Uncertain matches remain unresolved rather than being guessed.')}</p>
        </div>
      {/if}

      <div class="actual-layer-grid">
        {#each omicLayers as layer}
          {#if analysisResult.layers[layer]}
            {@const result = analysisResult.layers[layer]}
            <article>
              <div class="layer-title">
                <div>
                  <span class="method-tag">{omicLabel(layer)}</span>
                  <h3>{result.contrast || 'No valid contrast'}</h3>
              <small data-testid={'multiomics-number-source-'+layer}>{analysisResult.numericSourceByLayer?.[layer] === 'native_r_complete_table'
                      ? t('Chiffres : R natif (table complète)', 'Numbers: native R (complete table)')
                : t('Chiffres : exploration JS', 'Numbers: JavaScript exploration')}</small>
                </div>
            <button class="text-button" type="button" onclick={() => downloadLayerCsv(layer)}>CSV ↓</button>
              </div>

              {#if result.error}
                <p class="error">{result.error}</p>
              {:else}
                <p class="muted">
                  {result.mode === 'exploratory-multiblock-pca' || result.mode === 'exploratory-single-omic-pca'
                    ? result.inferenceMethod
                    : result.mode?.startsWith('outcome-')
                      ? `${result.mode} · ${result.inferenceMethod || 'model'}`
                      : result.mode?.includes('multi-group')
                        ? `${t('tailles des groupes', 'group sizes')}: ${result.groupSizes.join(' / ')} · ${result.inferenceMethod || result.mode}`
                        : result.mode === 'random-intercept-longitudinal-model'
                          ? `${result.mode} · ${result.inferenceMethod || 'mixed model'}`
                          : `n=${result.groupSizes[0]} vs ${result.groupSizes[1]} · ${result.mode} · ${result.inferenceMethod || 'inference'}`}
                </p>
                <div class="preprocess-list">
              {#each result.steps as step}<span>{step}</span>{/each}
              {#if result.replicateGroups?.length}<span>{result.replicateGroups.length} technical-replicate group(s) averaged on the transformed scale</span>{/if}
                </div>
                <p class="selection-rule">{result.selectionRule}</p>
            {#if result.inferencePolicy}<p class="selection-rule"><strong>Inference:</strong> {result.inferencePolicy}</p>{/if}
            <div class="feature-table" data-testid="multiomics-feature-inference">
                  <div class="feature-head">
                    <b>{t('Variable', 'Feature')}</b>
  <b>{result.effectScale === 'log2' ? 'Fold ratio' : t('Effet', 'Effect')}
    <span class="help-tip" tabindex="0" data-tooltip={result.effectScale === 'log2'
                          ? t('Rapport comparaison/référence calculé comme 2^effect sur une échelle log2. 2 = deux fois plus, 0,5 = deux fois moins.', 'Comparison/reference ratio computed as 2^effect on a log2 scale. 2 = twice as high; 0.5 = half as high.')
      : t('Direction et amplitude du contraste ou du coefficient estimé. Toujours lire le contraste affiché et son échelle.', 'Direction and magnitude of the estimated contrast or coefficient. Always read the displayed contrast and its scale.')}>?</span>
                    </b>
  <b>p
    <span class="help-tip" tabindex="0" data-tooltip={t('p-value brute du test pour cette variable, avant correction des tests multiples. Une petite p-value seule ne suffit pas lorsqu’on teste beaucoup de variables.', 'Raw p-value for this feature before multiple-testing correction. A small p-value alone is not sufficient when many features are tested.')}>?</span>
                    </b>
  <b>q BH
    <span class="help-tip" tabindex="0" data-tooltip={t('p-value ajustée par Benjamini–Hochberg pour contrôler le taux de faux positifs parmi les résultats retenus. Ici q ≤ 0,10 est un seuil exploratoire.', 'Benjamini–Hochberg adjusted p-value controlling the false discovery rate among selected results. Here q ≤ 0.10 is an exploratory threshold.')}>?</span>
                    </b>
                  </div>
                  {#each (result.numericalProvenance === 'exact_native_r_reference_table' ? result.selected : result.rows).slice(0, 10) as row}
                    <div>
                      <code>{row.feature}</code>
                      {#if row.foldRatio != null}
                        <span
                          class:negative={row.foldRatio != null && row.foldRatio < 1}
                        >
                          {row.foldRatio == null
                            ? row.effect == null ? '—' : row.effect.toFixed(3)
                            : `${row.foldRatio.toFixed(2)}×`}
                        </span>
                      {:else if row.exponentiatedEffect != null}
                        <span class:negative={row.exponentiatedEffect < 1}>
                          {row.exponentiatedEffect.toPrecision(3)}×
                          {#if row.exponentiatedCiLow != null && row.exponentiatedCiHigh != null}
                            <small>[{row.exponentiatedCiLow.toPrecision(3)}–{row.exponentiatedCiHigh.toPrecision(3)}]</small>
                          {/if}
                          <small>(β={row.effect == null ? '—' : row.effect.toPrecision(3)})</small>
                        </span>
                      {:else}
                        <span class:negative={row.effect < 0}>
                          {row.effect == null ? '—' : row.effect.toPrecision(3)}
                          {#if row.ciLow != null && row.ciHigh != null}
                            <small>[{row.ciLow.toPrecision(3)}–{row.ciHigh.toPrecision(3)}]</small>
                          {/if}
                        </span>
                      {/if}
                      <span>{row.pValue == null ? '—' : row.pValue.toPrecision(3)}</span>
                  <span data-scientific-q={row.qValue == null ? '' : String(row.qValue)}>{row.qValue == null ? '—' : row.qValue.toPrecision(3)}</span>
                    </div>
                  {/each}
                </div>
              {/if}
            </article>
          {/if}
        {/each}
      </div>

  {#if analysisResult.protocol?.analysisIntent !== 'confirmatory' && !['explore','outcome'].includes(analysisResult.protocol?.objective)}
        <div class="integration-result cross-result">
          <div class="integration-head">
            <div>
              <p class="eyebrow">{t('Intégration inter-omique directe', 'Direct cross-omics integration')}</p>
              <h3>{t('Quelles relations moléculaires changent entre les conditions biologiques ?', 'Which molecular relationships change between the biological conditions?')}</h3>
            </div>
            <span>Spearman + Fisher z + BH-FDR</span>
          </div>
          <div class="api-summary">
      <span><strong>{analysisResult.crossOmics?.testedPairs ?? 0}</strong> matched cross-omic pairs tested</span>
      <span><strong>{analysisResult.crossOmics?.significantPairs ?? 0}</strong> pairs with q ≤ 0.10</span>
          </div>
          {#if analysisResult.crossOmics?.pairs?.length}
            <div class="cross-table">
              <div class="cross-head">
                <b>{t('Paire', 'Pair')}</b>
  <b>Pattern <span class="help-tip" tabindex="0" data-tooltip={t('Résumé déterministe du changement de corrélation : gain, perte, renforcement, affaiblissement ou inversion de signe.', 'Deterministic summary of the correlation change: gained, lost, strengthened, weakened or sign reversal.')}>?</span></b>
  <b>r reference <span class="help-tip" tabindex="0" data-tooltip={t('Corrélation de Spearman entre les deux variables dans le groupe de référence.', 'Spearman correlation between the two features in the reference group.')}>?</span></b>
  <b>r comparison <span class="help-tip" tabindex="0" data-tooltip={t('Même corrélation de Spearman dans le groupe comparé.', 'The same Spearman correlation in the comparison group.')}>?</span></b>
  <b>Δr <span class="help-tip" tabindex="0" data-tooltip={t('Différence r_comparison − r_reference. Une grande valeur absolue indique un changement important de couplage entre les deux conditions.', 'Difference r_comparison − r_reference. A large absolute value indicates a strong change in coupling between conditions.')}>?</span></b>
  <b>q BH <span class="help-tip" tabindex="0" data-tooltip={t('FDR Benjamini–Hochberg appliqué à toutes les paires inter-omiques testées. Il indique si le changement de corrélation reste crédible après correction multiple.', 'Benjamini–Hochberg FDR across all tested cross-omics pairs. It indicates whether the correlation change remains credible after multiple-testing correction.')}>?</span></b>
              </div>
              {#each analysisResult.crossOmics.pairs.slice(0, 15) as pair}
                <div>
            <span><code>{pair.featureA}</code> <small>{pair.layerA}</small> ↔ <code>{pair.featureB}</code> <small>{pair.layerB}</small></span>
                  <span class="pattern-tag">{pair.pattern.replaceAll('_', ' ')}</span>
                  <span>{Number.isFinite(pair.rReference) ? pair.rReference.toFixed(2) : '—'}</span>
                  <span>{Number.isFinite(pair.rComparison) ? pair.rComparison.toFixed(2) : '—'}</span>
                  <strong>{Number.isFinite(pair.deltaR) ? pair.deltaR.toFixed(2) : '—'}</strong>
                  <span>{pair.qValue == null ? '—' : pair.qValue.toPrecision(3)}</span>
                </div>
              {/each}
            </div>
            <p class="note">{analysisResult.crossOmics.method}. Only subjects represented in both layers contribute to a pair; no imputation is used here.</p>
          {:else}
            <p class="muted">{t('Aucun test de changement de corrélation inter-omique n’était estimable avec les tailles de groupes et sujets appariés actuels.', 'No cross-omic correlation test was estimable with the current group sizes and matched subjects.')}</p>
          {/if}
        </div>
      {/if}

      <div class="integration-result">
        <div class="integration-head">
          <div>
            <p class="eyebrow">{t('Résultat intégré par voies', 'Integrated pathway result')}</p>
            <h3>{t('Reactome reçoit ensemble les identifiants de gènes, protéines et métabolites sélectionnés', 'Reactome receives the selected gene/protein/metabolite identifiers together')}</h3>
          </div>
          <span>{t('ORA déterministe', 'deterministic ORA')}</span>
        </div>

        {#if analysisResult.reactome}
          {#if analysisResult.reactome.combined?.assayUniverse?.multiplicityStatus !== 'full_assay_universe_family'}
        <p class="api-error" data-testid="multiomics-reactome-fdr-warning">
          {t('FDR sur l’univers mesuré non estimable : couverture de l’API ou univers incohérent. Aucun résultat de voie ne peut être déclaré significatif sur cette base.', 'Measured-universe FDR not estimable: incomplete API coverage or inconsistent background. No pathway may be called significant on this basis.')}
          {analysisResult.reactome.combined?.assayUniverse?.warning || ''}
            </p>
          {:else}
        <p class="note" data-testid="multiomics-reactome-fdr-family">
              {t('Correction BH sur', 'BH correction across')}
              <strong>{analysisResult.reactome.combined.assayUniverse.hypothesesTested}</strong>
              {t('voies testables dans l’univers mesuré (et non uniquement les voies déjà enrichies).', 'testable pathways within the measured universe, not only previously enriched pathways.')}
            </p>
          {/if}
          <div class="api-summary">
        <span><strong>{analysisResult.reactome.combined.pathwaysFound}</strong> pathways found</span>
        <span><strong>{analysisResult.reactome.combined.identifiersNotFound}</strong> identifiers not found</span>
          </div>
      <div class="pathway-table" data-testid="multiomics-pathways">
            <div class="pathway-head">
              <b>{t('Voie', 'Pathway')}</b>
  <b>{t('FDR univers assay', 'Assay-universe FDR')}
    <span class="help-tip" tabindex="0" data-tooltip={t('FDR calculée par test hypergéométrique puis corrigée sur toutes les voies présentes dans l’univers réellement mesuré. Aucune FDR locale n’est affichée si la réponse Reactome est incomplète ; la FDR Reactome d’origine ne sert pas de remplacement.', 'FDR from a local hypergeometric test adjusted over all pathways supported by the actual measured assay universe. No local FDR is displayed if the Reactome response is incomplete; default Reactome FDR is never substituted.')}>?</span>
              </b>

              <b>
                RNA

                <span
                  class="help-tip"
                  tabindex="0"
                  data-tooltip={t('FDR Reactome calculée uniquement avec les variables transcriptomiques sélectionnées.', 'Reactome FDR using only selected transcriptomic features.')}
                >?</span>
              </b>

              <b>
                {t('Protéine', 'Protein')}

                <span
                  class="help-tip"
                  tabindex="0"
                  data-tooltip={t('FDR Reactome calculée uniquement avec les variables protéomiques sélectionnées.', 'Reactome FDR using only selected proteomic features.')}
                >?</span>
              </b>

              <b>
                {t('Métabolite', 'Metabolite')}

                <span
                  class="help-tip"
                  tabindex="0"
                  data-tooltip={t('FDR Reactome calculée uniquement avec les métabolites sélectionnés.', 'Reactome FDR using only selected metabolites.')}
                >?</span>
              </b>

              <b>
                {t('Couches ≤0,10', 'Layers ≤0.10')}

                <span
                  class="help-tip"
                  tabindex="0"
                  data-tooltip={t('Nombre de couches omiques qui soutiennent séparément cette voie avec une FDR Reactome ≤0,10. Une valeur élevée indique une convergence inter-omique.', 'Number of omics layers independently supporting this pathway at Reactome FDR ≤0.10. A higher value indicates cross-omics convergence.')}
                >?</span>
              </b>
            </div>
            {#each analysisResult.reactome.consensus.slice(0, 15) as pathway}
              <div>
            <a href={`https://reactome.org/content/detail/${pathway.id}`} target="_blank" rel="noreferrer">{pathway.name}</a>
                <span>{Number.isFinite(pathway.assayUniverseFdr) ? pathway.assayUniverseFdr.toPrecision(3) : '—'}</span>
            <span>{pathway.layerEvidence.transcriptomics?.fdr != null ? pathway.layerEvidence.transcriptomics.fdr.toPrecision(2) : '—'}</span>
            <span>{pathway.layerEvidence.proteomics?.fdr != null ? pathway.layerEvidence.proteomics.fdr.toPrecision(2) : '—'}</span>
            <span>{pathway.layerEvidence.metabolomics?.fdr != null ? pathway.layerEvidence.metabolomics.fdr.toPrecision(2) : '—'}</span>
                <strong>{pathway.supportingLayers}</strong>
              </div>
            {/each}
          </div>
          <p class="note">{t('Le classement est déterministe : nombre de couches avec FDR de voie ≤ 0,10, puis FDR Reactome combinée, puis couverture de la voie. Il s’agit d’un classement exploratoire, pas d’une probabilité postérieure ni d’un score causal.', 'Ranking is deterministic: number of omics layers with pathway FDR ≤0.10, then combined Reactome FDR, then pathway coverage. This is an exploratory ranking, not a posterior probability or causal score.')}</p>
      <p class="note"><strong>{t('Univers d’enrichissement :', 'Enrichment background:')} <span class="help-tip" tabindex="0" data-tooltip={t('L’univers est l’ensemble des variables qui auraient pu être sélectionnées après QC. Utiliser cet univers évite de comparer un panel ciblé à tous les gènes ou métabolites connus de la base.', 'The background universe is the set of features that could have been selected after QC. Using this universe avoids comparing a targeted panel against every gene or metabolite known to the database.')}>?</span></strong> {analysisResult.reactome.backgroundPolicy}. {analysisResult.reactome.backgroundCaveat} {t('Une case vide signifie que cette FDR ne peut pas être estimée. La FDR Reactome d’origine reste disponible séparément dans le JSON, mais ne remplace jamais cette estimation.', 'A dash means this FDR cannot be estimated. The original Reactome FDR remains available separately in the JSON, but never replaces this estimate.')}</p>
        {:else if analysisResult.reactomeError}
          <div class="api-error">
            <strong>{t('Les statistiques locales sont terminées ; Reactome n’a pas pu être joint.', 'Local statistics completed; Reactome could not be reached.')}</strong>
            <p>{analysisResult.reactomeError}</p>
            <p>{t('Les matrices ne quittent toujours pas le navigateur. Relancez plus tard ou désactivez Reactome pour utiliser uniquement les résultats statistiques locaux.', 'The matrices have still not left the browser. Re-run later or disable Reactome to use the local statistical output only.')}</p>
          </div>
        {:else}
          <p class="muted">{t('L’interrogation de Reactome était désactivée pour cette analyse.', 'Reactome querying was disabled for this run.')}</p>
        {/if}
      </div>

    </details>
    <div class="evidence-layers">
      <article>
        <strong>{t('Observé', 'Observed')}</strong>
        <p>{t('Valeurs lues directement dans les matrices importées.', 'Values parsed from the uploaded matrices.')}</p>
      </article>
      <article>
        <strong>{t('Inférence statistique', 'Statistical inference')}</strong>
      <p>{analysisResult.protocol?.objective === 'outcome'
            ? t('Régression adaptée au type d’outcome avec batch/covariables comme termes de nuisance et correction Benjamini–Hochberg.', 'Outcome-specific regression with batch/covariates as nuisance terms and Benjamini–Hochberg correction.')
            : analysisResult.protocol?.objective === 'explore'
              ? t('ACP multi-blocs équilibrée sur les sujets communs, sans variable cible.', 'Balanced multi-block PCA across shared subjects without a target variable.')
    : t('Contrastes déterministes adaptés au design, différences de corrélations inter-omiques de Spearman et correction Benjamini–Hochberg.', 'Design-aware deterministic contrasts, cross-omics Spearman correlation differences and Benjamini–Hochberg correction.')}</p>
      </article>
      <article>
        <strong>{t('Connaissance externe', 'External knowledge')}</strong>
        <p>{t('Mapping et sur-représentation Reactome à partir des identifiants moléculaires sélectionnés.', 'Reactome pathway mapping and over-representation from the selected molecular identifiers.')}</p>
      </article>
      <article>
        <strong>{t('Non inféré', 'Not inferred')}</strong>
        <p>{t('Aucun récit LLM, aucune causalité ni aucun mécanisme biologique inventé n’est généré.', 'No LLM narrative, causal claim or invented biological mechanism is generated.')}</p>
      </article>
    </div>
  </section>
{/if}

<details class="advanced-parameters methodology-details" data-testid="multiomics-methodology">
  <summary><strong>{t('Méthodes, hypothèses et bases utilisées', 'Methods, assumptions and databases')}</strong> <span>{t('Consulter les calculs en détail', 'Inspect how results were calculated')}</span></summary>
  <section class="panel">
    <div class="section-head">
      <div>
        <p class="eyebrow">{t('Méthode de calcul', 'Analysis method')}</p>
        <h2>{t('Le moteur de décision est déterministe et inspectable', 'The decision engine is deterministic and inspectable')}</h2>
      </div>
      <p>{t('La branche statistique est sélectionnée à partir de règles explicites liées au design et reste inspectable.', 'The statistical branch is selected from explicit study-design rules and remains inspectable.')}</p>
    </div>

    <div class="plan">
      <div>
        <span class="method-tag">{t('Workflow sélectionné', 'Selected workflow')}</span>
        <h3>{analysisPlan}</h3>
        <p class="muted">Protocol: {studySetting.replaceAll('_', ' ')} · {designType} · {groupCount === '3plus' ? '≥3' : groupCount} group(s) · {sampleOverlap.replaceAll('_', ' ')}.</p>
      </div>
      <ul>
      <li><strong>Same specimen:</strong> assays are linked by <code>sample_id</code>, never by fuzzy name matching.</li>
      <li><strong>{t('Sujet répété :', 'Repeated subject:')}</strong> {t('deux temps utilisent un changement intra-sujet ; ≥3 temps numériques utilisent une pente individuelle.', 'two time points use within-subject change; ≥3 numeric-labelled times use an individual slope.')}</li>
      <li><strong>{t('Design apparié :', 'Paired design:')}</strong> {t('les différences intra-sujet utilisent un test de permutation par inversion de signe. Les crossovers restent refusés tant que période et séquence ne sont pas modélisées.', 'within-subject differences use a sign-flip permutation test. Crossover designs are refused until period/sequence effects are modelled.')}</li>
      <li><strong>Technical replicate:</strong> multiple assays with the same <code>sample_id + omic</code> are flagged before analysis.</li>
      <li><strong>{t('Batch technique :', 'Technical batch:')}</strong> {t('une confusion complète batch–condition/temps/outcome bloque l’inférence. Les batches multiples non confondus sont ajustés explicitement par résidualisation OLS variable par variable.', 'complete batch–condition/time/outcome confounding blocks inference. Multiple non-confounded batches are explicitly adjusted by feature-wise OLS residualisation.')}</li>
      <li><strong>{t('Covariables :', 'Covariates:')}</strong> {t('les colonnes sélectionnées sont intégrées explicitement à l’ajustement ou au modèle d’outcome ; aucune covariable n’est choisie automatiquement.', 'selected columns enter adjustment or outcome models explicitly; no covariate is chosen automatically.')}</li>
      <li><strong>{t('Omiques partielles :', 'Partial omics:')}</strong> {t('une couche absente est distinguée d’une valeur manquante dans une matrice observée.', 'absent layers are distinguished from missing values inside an observed matrix.')}</li>
      <li><strong>{t('Outcome :', 'Outcome:')}</strong> {t('le type déclaré sélectionne régression linéaire, logistique, Poisson, ANOVA multiclasse ou Cox.', 'the declared type selects linear, logistic, Poisson, multiclass ANOVA or Cox regression.')}</li>
      </ul>
    </div>

    <div class="database-bridge">
      <div class="section-head compact">
        <div>
          <p class="eyebrow">{t('Identifiant → biologie', 'Identifier → biology bridge')}</p>
          <h3>{t('Les bases publiques apportent le contexte biologique après validation des mesures', 'Public databases provide meaning after the measurements are validated')}</h3>
        </div>
        <p>{t('Seuls les identifiants moléculaires nécessaires aux requêtes d’annotation, de voies ou de réseaux sont envoyés aux API publiques. Les identifiants sujets, métadonnées et matrices d’abondance restent hors de ces appels.', 'Only molecular identifiers needed for annotation/pathway/network queries should be sent to public APIs. Subject IDs, metadata and abundance matrices remain outside these calls.')}</p>
      </div>
      <div class="database-grid">
        {#each databaseRegistry as db}
          <article>
            <span>{$language === 'en' ? db.scope : db.scopeFr} · {$language === 'en' ? db.status : db.statusFr}</span>
            <h3>{db.name}</h3>
            <p>{$language === 'en' ? db.role : db.roleFr}</p>
            <a href={db.url} target="_blank" rel="noreferrer">{t('API officielle / documentation ↗', 'Official API / documentation ↗')}</a>
          </article>
        {/each}
      </div>
      <div class="resolution-flow">
        <code>{transcriptomicsIdType}</code>
        <b>+</b>
        <code>{proteomicsIdType}</code>
        <b>+</b>
        <code>{metabolomicsIdType}</code>
        <span>{t('→ moteur actuel : résolution ChEBI conservatrice des métabolites (optionnelle) → intégration de voies Reactome. Ensembl, UniProt, UniChem, KEGG et STRING sont des connecteurs explicitement prévus, pas des appels automatiques cachés.', '→ current engine: conservative ChEBI resolution for metabolites (optional) → Reactome pathway integration. Ensembl, UniProt, UniChem, KEGG and STRING are explicit registry targets, not hidden automatic calls.')}</span>
      </div>

      <details class="aliases">
        <summary>{t('Suite de validation publique du moteur déterministe', 'Public validation suite used to test the deterministic engine')}</summary>
        <div class="alias-grid">
        <article><strong>Nutrimouse</strong><span>engine truth test</span><p>PPARα-dependent transcript/metabolite signal, CYP3A11 and lipid-pathway recovery.</p></article>
        <article><strong>TCGA breast</strong><span>engine truth test</span><p>HER2/LumA contrast plus Basal/Her2/LumA multi-group inference.</p></article>
        <article><strong>TCGA breast · predictive outcome</strong><span>nested-CV test</span><p>Her2 vs LumA prediction: 5-fold outer nested CV, AUC 0.996; supervised multiblock score-target correlation r=0.915 across 105 tumours.</p></article>
        <article><strong>IntLIM NCI-60 + BRCA</strong><span>cross-omics truth test</span><p>Published condition-dependent gene–metabolite correlation changes.</p></article>
        <article><strong>AgingHFCD</strong><span>3-omics reference test</span><p>RNA, protein and metabolite effect directions checked against public reference results.</p></article>
        <article><strong>LRRK2 G2019S</strong><span>2-omics reference test</span><p>RNA/protein direction agreement plus RAB/endocytic biology.</p></article>
        <article><strong>STATegra</strong><span>time-course + replicates</span><p>Real Ikaros time course and 36-sample metabolomics replicate design; processed representations are checked for direction and scale consistency.</p></article>
        <article><strong>PaintOmics planted multi-omics</strong><span>known ground truth</span><p>RNA/protein convergence is checked against a planted molecular module and recorded pathway truth set.</p></article>
        <article><strong>Bioconductor missRows NCI-60</strong><span>partial-block integration test</span><p>60 subjects across the union of layers, 40 complete across RNA + protein; latent-only filling is used without inventing values for differential tests.</p></article>
        </div>
        <p class="note">{t('La suite automatisée est volontairement hétérogène : elle teste la vérité biologique, l’appariement des échantillons, l’inférence multi-groupes, la structure des réplicats et les statistiques inter-omiques, et pas seulement l’exécution du code.', 'The automated suite is intentionally heterogeneous: it tests biological truth, sample matching, multi-group inference, replicate structure and cross-omics statistics rather than only checking that code executes.')}</p>
      <a href="https://github.com/rberrah/rberrah.github.io/actions/workflows/multiomics-public-benchmark.yml" target="_blank" rel="noreferrer">{t('Ouvrir le benchmark public ↗', 'Open the public benchmark workflow ↗')}</a>
      </details>
    </div>
  </section>
</details>

{#if helpTooltip.visible}
  <div
    class="global-help-tooltip {helpTooltip.placement}"
    data-testid="global-help-tooltip"
    role="tooltip"
    style={'left:' + helpTooltip.left + 'px;top:' + helpTooltip.top + 'px'}
  >
    {helpTooltip.text}
  </div>
{/if}

<style>
  .simple-steps {
    display: flex;
    align-items: center;
    flex-wrap: wrap;
    gap: 1.2rem;
    margin: 0 0 var(--space-6);
    padding: 0 0 var(--space-6);
    border-bottom: 1px solid var(--border);
  }

  .simple-steps span {
    display: flex;
    align-items: center;
    gap: .5rem;
    color: var(--text-secondary);
    font-weight: 600;
    font-size: .86rem;
  }

  .simple-steps b {
    display: inline-grid;
    place-items: center;
    width: 1.7rem;
    height: 1.7rem;
    border-radius: 50%;
    border: 1px solid var(--border-strong);
    color: var(--text-primary);
  }

  .simple-hint {
    margin: .35rem 0 .9rem;
    color: var(--text-secondary);
    font-size: .9rem;
  }

  .simple-warning {
    padding: .7rem 1rem;
    border-inline-start: 3px solid var(--warning);
    background: var(--bg-secondary);
  }

  .simple-disclosure,
  .advanced-parameters,
  .advanced-results,
  .templates-details {
    margin: 1rem 0;
    border: 1px solid var(--border);
    border-radius: 12px;
    background: var(--bg-primary);
    padding: .75rem 1rem;
  }

  .simple-disclosure summary,
  .advanced-parameters summary,
  .advanced-results summary,
  .templates-details summary {
    display: flex;
    align-items: center;
    justify-content: space-between;
    flex-wrap: wrap;
    gap: .5rem;
    cursor: pointer;
    font-size: .92rem;
    font-weight: 650;
  }

  .simple-disclosure > :not(summary),
  .advanced-results > :not(summary),
  .advanced-parameters > :not(summary),
  .templates-details > :not(summary) {
    margin-top: 1rem;
  }

  .advanced-parameters > summary span,
  .methodology-details > summary span {
    color: var(--text-muted);
    font-size: .8rem;
    font-weight: 400;
  }

  .advanced-parameters > .panel {
    margin: 1rem 0;
  }

  .advanced-parameters .section-head p,
  .simple-study .section-head p {
    max-width: 70ch;
  }

  .templates-details {
    margin: 0;
    align-self: start;
    min-width: 0;
  }

  .inline-optional {
    min-width: 0;
  }

  .inline-optional label {
    display: grid;
    margin-top: .7rem;
  }

  .simple-file-values {
    padding: 1rem 0;
  }

  .simple-file-values > p {
    margin: .3rem 0 1rem;
    font-size: .87rem;
    color: var(--text-secondary);
  }

  .simple-file-values small {
    display: block;
    margin-top: .8rem;
    color: var(--text-muted);
  }


  .simple-qc-status {
    display: flex;
    align-items: center;
    flex-wrap: wrap;
    gap: .7rem;
    padding: .85rem 0;
    border-top: 1px solid var(--border);
    font-size: .85rem;
  }
  .simple-qc-status span {
    border: 1px solid var(--border);
    padding: .25rem .6rem;
    border-radius: 999px;
  }

  .simple-qc-status .qc-alert {
    border-color: var(--warning);
  }

  .export-detail {
    margin: 0;
  }

  .export-detail button {
    margin-top: .55rem;
  }

  .methodology-details {
    margin-top: 1.2rem;
  }

  @media (max-width:680px) {
    .simple-steps {
      gap: .6rem;
    }
    .simple-steps span {
      font-size: .76rem;
    }
    .simple-disclosure summary,
    .advanced-parameters summary,
    .advanced-results summary,
    .templates-details summary {
      display: block;
    }
    .simple-qc-status {
      align-items: flex-start;
    }
  }

  .hero {
    max-width: 920px;
    padding: var(--space-12) 0 var(--space-8);
  }

  .tool-back {
    display: inline-block;
    margin-bottom: var(--space-4);
    font-size: var(--text-sm);
  }

  .help-tip {
    position: relative;
    display: inline-grid;
    place-items: center;
    width: 1.05rem;
    height: 1.05rem;
    margin-left: 3px;
    border: 1px solid var(--border-strong);
    border-radius: 50%;
    font: 700 0.72rem/1 var(--font-sans);
    color: var(--text-secondary);
    cursor: help;
    vertical-align: middle;
  }

  .help-tip::after {
    display: none;
  }

  .global-help-tooltip {
    position: fixed;
    z-index: 10000;
    width: min(320px,76vw);
    padding: 9px 11px;
    border: 1px solid var(--border-strong);
    border-radius: 8px;
    background: var(--bg-primary);
    box-shadow: 0 10px 30px rgba(0,0,0,.22);
    color: var(--text-primary);
    font: 400 var(--text-xs)/1.45 var(--font-sans);
    text-align: left;
    white-space: normal;
    pointer-events: none;
  }

  .global-help-tooltip.above {
    transform: translate(-50%,-100%);
  }

  .global-help-tooltip.below {
    transform: translate(-50%,0);
  }

  .help-tip:focus-visible {
    outline: 2px solid var(--accent-pk);
    outline-offset: 2px;
  }

  h1 {
    font-size: clamp(2.4rem, 6vw, 4.8rem);
    line-height: .98;
    max-width: 14ch;
    margin: var(--space-3) 0 var(--space-6);
    letter-spacing: -.045em;
  }

  h2 {
    margin: 0;
    font-size: var(--text-2xl);
  }

  h3 {
    margin: 0 0 var(--space-2);
    font-size: var(--text-lg);
  }

  p {
    line-height: 1.65;
  }

  .lede {
    color: var(--text-secondary);
    font-size: var(--text-lg);
    max-width: 72ch;
  }

  .privacy {
    margin-top: var(--space-6);
    border-left: 3px solid var(--accent-ai);
    padding: var(--space-3) var(--space-4);
    background: var(--bg-secondary);
    color: var(--text-secondary);
    max-width: 76ch;
  }

  .privacy strong {
    color: var(--text-primary);
  }

  .workflow {
    display: grid;
    grid-template-columns: repeat(6, 1fr);
    border-block: 1px solid var(--border-strong);
    margin-bottom: var(--space-12);
  }

  .workflow div {
    padding: var(--space-4);
    border-right: 1px solid var(--border-subtle);
  }

  .workflow div:last-child {
    border-right: 0;
  }

  .workflow span {
    font-family: var(--font-mono);
    color: var(--accent-pk);
    display: block;
    margin-bottom: var(--space-2);
  }

  .workflow strong,
  .workflow small {
    display: block;
  }

  .workflow small {
    color: var(--text-muted);
    margin-top: 2px;
  }

  .panel,
  .results-preview {
    margin-top: var(--space-12);
    border-top: 1px solid var(--border-strong);
    padding-top: var(--space-6);
  }

  .section-head {
    display: flex;
    justify-content: space-between;
    align-items: end;
    gap: var(--space-6);
    margin-bottom: var(--space-6);
  }

  .section-head > p {
    max-width: 48ch;
    color: var(--text-secondary);
    margin: 0;
    font-size: var(--text-sm);
  }

  .form-grid {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: var(--space-4);
    max-width: 940px;
  }

  .form-grid label,
  .mapping-grid label {
    display: grid;
    gap: var(--space-2);
  }

  .form-grid label > span,
  .mapping-grid label > span {
    font-weight: 650;
  }

  small {
    color: var(--text-muted);
    font-weight: 400;
  }

  .wide {
    grid-column: 1 / -1;
  }

  input,
  select {
    width: 100%;
    box-sizing: border-box;
    padding: 11px 12px;
    color: var(--text-primary);
    background: var(--bg-primary);
    border: 1px solid var(--border-strong);
    border-radius: var(--radius);
    font: inherit;
  }

  input:focus,
  select:focus {
    outline: 2px solid var(--focus-ring);
    border-color: var(--accent-pk);
  }

  .method-card {
    max-width: 900px;
    margin-top: var(--space-6);
    padding: var(--space-5);
    background: var(--bg-secondary);
    border-left: 3px solid var(--accent-pd);
  }

  .method-card p {
    margin-bottom: 0;
    color: var(--text-secondary);
  }

  .method-tag {
    display: inline-block;
    margin-bottom: var(--space-2);
    font-family: var(--font-mono);
    font-size: var(--text-xs);
    color: var(--accent-pd);
  }

  .identity-grid,
  .validation {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
    margin-top: var(--space-6);
  }

  .identity-grid article,
  .validation article {
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-4);
  }

  .identity-grid code {
    display: block;
    color: var(--accent-pk);
    margin-bottom: var(--space-2);
  }

  .identity-grid strong,
  .identity-grid p {
    display: block;
  }

  .identity-grid p,
  .validation small {
    color: var(--text-secondary);
    margin-bottom: 0;
  }

  .omics-question-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
  }

  .omics-question-grid article {
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-4);
    background: var(--bg-secondary);
  }

  .omics-question-grid label {
    display: grid;
    gap: 6px;
    margin-top: var(--space-3);
  }

  .omics-question-grid label > span {
    font-weight: 650;
    font-size: var(--text-sm);
  }

  .contract {
    display: grid;
    grid-template-columns: 1.35fr .65fr;
    gap: var(--space-5);
  }

  .contract > div {
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-5);
    background: var(--bg-secondary);
  }

  pre {
    overflow-x: auto;
    padding: var(--space-3);
    background: var(--bg-primary);
    border: 1px solid var(--border-subtle);
    font-size: var(--text-xs);
  }

  .actions {
    display: flex;
    flex-wrap: wrap;
    gap: var(--space-2);
    margin-top: var(--space-4);
  }

  .demo-card {
    background: linear-gradient(135deg, color-mix(in srgb, var(--accent-ai) 10%, transparent), var(--bg-secondary)) !important;
  }

  .demo-links {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    margin-top: var(--space-4);
    font-size: var(--text-sm);
  }

  .dictionary,
  .aliases {
    margin-top: var(--space-5);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-4);
  }

  .dictionary summary,
  .aliases summary {
    cursor: pointer;
    font-weight: 700;
  }

  .dictionary-table {
    margin-top: var(--space-4);
  }

  .dictionary-table > div {
    display: grid;
    grid-template-columns: .8fr 1.2fr .8fr 2fr;
    gap: var(--space-3);
    padding: 8px 0;
    border-bottom: 1px solid var(--border-subtle);
    font-size: var(--text-sm);
  }

  .dictionary-table > div:first-child {
    color: var(--text-muted);
  }

  .aliases summary {
    cursor: pointer;
    font-weight: 700;
  }

  .alias-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-3);
    margin-top: var(--space-4);
  }

  .alias-grid article {
    background: var(--bg-secondary);
    padding: var(--space-3);
    border-radius: var(--radius);
  }

  .alias-grid article > span {
    margin-left: 8px;
    font-family: var(--font-mono);
    font-size: 10px;
    color: var(--text-muted);
    text-transform: uppercase;
  }

  .alias-grid p {
    font-size: var(--text-xs);
    color: var(--text-secondary);
    word-break: break-word;
    margin-bottom: 0;
  }

  .note {
    color: var(--text-secondary);
    font-size: var(--text-sm);
  }

  .uploads {
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: var(--space-4);
  }

  .uploads label {
    display: grid;
    gap: var(--space-2);
    padding: var(--space-5);
    border: 1px dashed var(--border-strong);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .uploads label.loaded {
    border-style: solid;
    border-color: var(--accent-pd);
  }

  .uploads label > span,
  .uploads label > small {
    color: var(--text-secondary);
  }

  .uploads input {
    background: var(--bg-primary);
  }

  .mapping {
    margin-top: var(--space-6);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-5);
  }

  .covariate-picker {
    margin-top: var(--space-5);
    padding-top: var(--space-4);
    border-top: 1px solid var(--border-subtle);
  }

  .covariate-picker h3 {
    margin-bottom: 4px;
  }

  .covariate-picker > div > p {
    margin-top: 0;
    color: var(--text-secondary);
  }

  .covariate-options {
    display: flex;
    flex-wrap: wrap;
    gap: 8px;
    margin: var(--space-3) 0;
  }

  .covariate-options label {
    display: flex;
    align-items: center;
    gap: 7px;
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
    padding: 6px 10px;
    background: var(--bg-primary);
  }

  .covariate-options input {
    width: auto;
    margin: 0;
  }

  .qc-result {
    margin-top: var(--space-5);
  }

  .qc-grid {
    display: grid;
    grid-template-columns: repeat(3,minmax(0,1fr));
    gap: var(--space-4);
    margin-top: var(--space-4);
  }

  .qc-card {
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-4);
    background: var(--bg-secondary);
    min-width: 0;
  }

  .qc-card-head {
    display: flex;
    justify-content: space-between;
    gap: var(--space-3);
    align-items: start;
  }

  .qc-card-head div {
    display: grid;
    gap: 2px;
  }

  .qc-card-head small {
    color: var(--text-secondary);
  }

  .qc-card-head > span {
    font-size: var(--text-xs);
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
    padding: 3px 7px;
  }

  .qc-card-head > span.qc-warning {
    border-color: var(--warning);
    color: var(--warning);
  }

  .qc-badges {
    display: flex;
    flex-wrap: wrap;
    justify-content: flex-end;
    gap: 5px;
  }

  .qc-badges > span {
    font-size: var(--text-xs);
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
    padding: 3px 7px;
  }

  .qc-badges > span.qc-warning {
    border-color: var(--warning);
    color: var(--warning);
  }

  .qc-badges > span.qc-screening {
    border-color: var(--accent-pk);
    color: var(--text-secondary);
  }

  .ms-qc-controls {
    display: grid;
    grid-template-columns: repeat(2,minmax(0,1fr));
    gap: 8px;
    padding: 10px;
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-primary);
  }

  .ms-qc-controls > strong {
    grid-column: 1/-1;
  }

  .ms-qc-controls label {
    display: grid;
    gap: 4px;
  }

  .ms-qc-summary {
    display: flex;
    flex-wrap: wrap;
    gap: 5px;
    margin: 10px 0;
  }

  .ms-qc-summary > strong,
  .ms-qc-summary > span {
    font-size: 10px;
    font-family: var(--font-mono);
    padding: 4px 6px;
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
  }

  .qc-metrics {
    display: grid;
    grid-template-columns: repeat(3,1fr);
    gap: 8px;
    margin: var(--space-4) 0;
  }

  .qc-metrics div {
    display: grid;
    gap: 2px;
  }

  .qc-metrics b {
    font-family: var(--font-mono);
  }

  .qc-metrics span {
    color: var(--text-secondary);
    font-size: var(--text-xs);
  }

  .qc-bars {
    display: grid;
    gap: 5px;
    max-height: 150px;
    overflow: auto;
  }

  .qc-bars > div {
    display: grid;
    grid-template-columns: minmax(50px,.8fr) 2fr;
    gap: 7px;
    align-items: center;
    font-size: var(--text-xs);
  }

  .qc-bars span {
    overflow: hidden;
    text-overflow: ellipsis;
    white-space: nowrap;
  }

  .qc-bars i {
    display: block;
    height: 6px;
    border-radius: 99px;
    background: var(--border-subtle);
    overflow: hidden;
  }

  .qc-bars em {
    display: block;
    height: 100%;
    background: var(--accent-pk);
    min-width: 1px;
  }

  .qc-pca {
    margin-top: var(--space-4);
  }

  .qc-pca-head {
    display: flex;
    justify-content: space-between;
    gap: var(--space-2);
    margin-bottom: 6px;
  }

  .qc-scatter {
    position: relative;
    height: 180px;
    border: 1px solid var(--border-subtle);
    background: var(--bg-primary);
    overflow: hidden;
  }

  .qc-scatter .axis-x {
    position: absolute;
    left: 0;
    right: 0;
    top: 50%;
    height: 1px;
    background: var(--border-subtle);
  }

  .qc-scatter .axis-y {
    position: absolute;
    top: 0;
    bottom: 0;
    left: 50%;
    width: 1px;
    background: var(--border-subtle);
  }

  .qc-dot {
    position: absolute;
    width: 9px;
    height: 9px;
    transform: translate(-50%,50%);
    border: 1px solid var(--bg-primary);
    border-radius: 50%;
    background: var(--accent-pd);
    padding: 0;
  }

  .qc-dot:focus-visible {
    outline: 2px solid var(--accent-pk);
    outline-offset: 2px;
  }

  .interpretation-box {
    margin-top: var(--space-5);
    padding: var(--space-5);
    border-left: 3px solid var(--accent-pd);
    background: var(--bg-secondary);
  }

  .interpretation-grid {
    display: grid;
    grid-template-columns: repeat(2, minmax(0,1fr));
    gap: var(--space-3);
    margin-top: var(--space-4);
  }

  .interpretation-grid article {
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-primary);
  }

  .interpretation-grid p {
    color: var(--text-secondary);
    font-size: var(--text-sm);
    margin-bottom: 0;
  }

  .reference-backend-box {
    margin-top: var(--space-4);
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-primary);
  }

  .reference-backend-box > div:first-child {
    display: flex;
    align-items: center;
    gap: 5px;
  }

  .backend-controls {
    display: grid;
    grid-template-columns: minmax(180px,.8fr) minmax(220px,1.4fr) auto;
    gap: 8px;
    margin: 10px 0 6px;
  }

  .btn-small {
    padding: 8px 11px;
  }

  .backend-packages {
    display: flex;
    flex-wrap: wrap;
    gap: 5px;
    margin-top: 8px;
  }

  .backend-packages span,
  .backend-method-grid article > span {
    font: 600 10px/1.2 var(--font-mono);
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
    padding: 4px 6px;
  }

  .backend-packages span.installed,
  .backend-method-grid article > span.installed {
    border-color: var(--accent-pd);
  }

  .backend-ok {
    color: var(--accent-pd);
  }

  .backend-bad {
    color: var(--accent-ai);
  }

  .backend-method-grid {
    display: grid;
    grid-template-columns: repeat(3,minmax(0,1fr));
    gap: 8px;
    margin-top: var(--space-4);
  }

  .backend-method-grid article {
    display: grid;
    grid-template-columns: 1fr auto;
    gap: 4px 8px;
    padding: 10px;
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-primary);
  }

  .backend-method-grid small {
    grid-column: 1/-1;
  }

  .backend-ms-detail {
    line-height: 1.55;
  }

  .predictive-result .api-summary {
    margin-top: var(--space-4);
  }

  .fold-grid {
    display: grid;
    gap: 5px;
    font-family: var(--font-mono);
    font-size: var(--text-xs);
    margin-top: var(--space-3);
  }

  .interpretation-glossary {
    margin-top: var(--space-4);
  }

  .mapping-head {
    display: flex;
    align-items: start;
    justify-content: space-between;
    gap: var(--space-4);
    margin-bottom: var(--space-4);
  }

  .mapping-head p {
    color: var(--text-secondary);
    margin: 0;
  }

  .mapping-head > span {
    font-family: var(--font-mono);
    font-size: var(--text-xs);
    border: 1px solid var(--border-strong);
    padding: 5px 8px;
    border-radius: 999px;
  }

  .mapping-head > span.ok {
    border-color: var(--accent-pd);
    color: var(--accent-pd);
  }

  .mapping-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
  }

  .validation {
    grid-template-columns: repeat(4, 1fr);
  }

  .validation article > span,
  .validation article > strong,
  .validation article > small {
    display: block;
  }

  .validation article > span {
    color: var(--text-secondary);
    font-size: var(--text-sm);
  }

  .validation article > strong {
    font-size: var(--text-2xl);
    margin: 4px 0;
  }

  .matrix-checks {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
    margin-top: var(--space-4);
  }

  .matrix-checks article {
    border-top: 2px solid var(--accent-pk);
    padding-top: var(--space-3);
  }

  .matrix-checks article > div {
    display: flex;
    justify-content: space-between;
    gap: 8px;
  }

  .matrix-checks article > div span {
    color: var(--text-muted);
    font-size: var(--text-xs);
  }

  .matrix-checks p {
    color: var(--text-secondary);
    margin-bottom: 4px;
  }

  .warning {
    color: var(--accent-ai) !important;
  }

  .error {
    color: var(--accent-ai);
    font-weight: 650;
  }

  .status {
    display: flex;
    gap: var(--space-3);
    align-items: baseline;
    margin-top: var(--space-5);
    padding: var(--space-3) var(--space-4);
    background: var(--bg-secondary);
    color: var(--text-secondary);
  }

  .status.ready strong {
    color: var(--accent-pd);
  }

  .run-box {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: var(--space-6);
    margin-top: var(--space-4);
    padding: var(--space-5);
    border: 1px solid var(--border-strong);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .run-box p {
    color: var(--text-secondary);
    max-width: 72ch;
  }

  .inline-check {
    display: flex;
    gap: 8px;
    align-items: center;
    font-size: var(--text-sm);
  }

  .inline-check input {
    width: auto;
  }

  .result-actions {
    display: flex;
    gap: 8px;
    flex-wrap: wrap;
  }

  .computed-summary {
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: var(--space-3);
  }

  .overlap-box {
    display: flex;
    justify-content: space-between;
    gap: var(--space-4);
    align-items: start;
    margin-top: var(--space-4);
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .overlap-pairs {
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    justify-content: flex-end;
  }

  .overlap-pairs span {
    font-size: var(--text-xs);
    padding: 5px 7px;
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
  }

  .identifier-resolution {
    margin-top: var(--space-5);
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
  }

  .mapping-result-table {
    font-size: var(--text-xs);
    overflow-x: auto;
  }

  .mapping-result-table > div {
    display: grid;
    grid-template-columns: 1fr 1.2fr .8fr .7fr;
    gap: 8px;
    padding: 7px 0;
    border-bottom: 1px solid var(--border-subtle);
    min-width: 650px;
  }

  .mapping-result-head {
    color: var(--text-muted);
  }

  .computed-summary article {
    padding: var(--space-3);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .computed-summary span,
  .computed-summary strong {
    display: block;
  }

  .computed-summary span {
    color: var(--text-muted);
    font-size: var(--text-xs);
  }

  .computed-summary strong {
    margin-top: 4px;
  }

  .actual-layer-grid {
    display: grid;
    grid-template-columns: repeat(3, minmax(0,1fr));
    gap: var(--space-4);
    margin-top: var(--space-5);
  }

  .actual-layer-grid > article {
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    padding: var(--space-4);
    background: var(--bg-secondary);
    min-width: 0;
  }

  .layer-title {
    display: flex;
    justify-content: space-between;
    gap: 8px;
    align-items: start;
  }

  .text-button {
    background: transparent;
    border: 0;
    color: var(--link);
    cursor: pointer;
    padding: 0;
    font: inherit;
    font-size: var(--text-sm);
  }

  .preprocess-list {
    display: flex;
    flex-wrap: wrap;
    gap: 5px;
    margin: var(--space-3) 0;
  }

  .preprocess-list span {
    font-size: 10px;
    font-family: var(--font-mono);
    padding: 3px 5px;
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
  }

  .selection-rule {
    color: var(--text-secondary);
    font-size: var(--text-xs);
  }

  .feature-table,
  .pathway-table {
    font-size: var(--text-xs);
    overflow-x: auto;
  }

  .feature-table > div {
    display: grid;
    grid-template-columns: 1.4fr .7fr .7fr .7fr;
    gap: 8px;
    padding: 7px 0;
    border-bottom: 1px solid var(--border-subtle);
    align-items: center;
  }

  .feature-head {
    color: var(--text-muted);
  }

  .feature-table .negative {
    color: var(--accent-ai);
  }

  .api-summary {
    display: flex;
    gap: var(--space-4);
    flex-wrap: wrap;
    margin: var(--space-3) 0;
  }

  .pathway-table > div {
    display: grid;
    grid-template-columns: minmax(220px,2fr) repeat(5,.65fr);
    gap: 8px;
    padding: 8px 0;
    border-bottom: 1px solid var(--border-subtle);
    align-items: center;
    min-width: 760px;
  }

  .pathway-head {
    color: var(--text-muted);
  }

  .cross-table {
    font-size: var(--text-xs);
    overflow-x: auto;
  }

  .cross-table > div {
    display: grid;
    grid-template-columns: minmax(280px,2.2fr) 1fr repeat(4,.65fr);
    gap: 8px;
    padding: 8px 0;
    border-bottom: 1px solid var(--border-subtle);
    align-items: center;
    min-width: 760px;
  }

  .cross-head {
    color: var(--text-muted);
  }

  .cross-table small {
    color: var(--text-muted);
    font-family: var(--font-mono);
  }

  .pattern-tag {
    font-family: var(--font-mono);
    font-size: 10px;
    text-transform: uppercase;
  }

  .api-error {
    padding: var(--space-4);
    border-left: 3px solid var(--accent-ai);
    background: var(--bg-primary);
  }

  .demo-results {
    scroll-margin-top: 90px;
  }

  .demo-story {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
  }

  .demo-story article {
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .demo-story ul {
    margin: 0;
    padding-left: 1.1rem;
    color: var(--text-secondary);
  }

  .demo-story li {
    margin-bottom: 6px;
  }

  .replicate-flow {
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    align-items: center;
    margin: var(--space-3) 0;
  }

  .replicate-flow code {
    padding: 4px 7px;
    border: 1px solid var(--border-strong);
    border-radius: 4px;
    background: var(--bg-primary);
  }

  .demo-omics-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
    margin-top: var(--space-5);
  }

  .demo-omics-grid > article {
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
  }

  .effect-row {
    display: grid;
    grid-template-columns: 1fr auto;
    gap: 4px 10px;
    align-items: center;
    margin: 12px 0;
  }

  .effect-row > span {
    font-family: var(--font-mono);
    color: var(--accent-pd);
  }

  .effect-row > span.down {
    color: var(--accent-ai);
  }

  .effect-row > span.stable {
    color: var(--text-muted);
  }

  .effect-row > i {
    grid-column: 1 / -1;
    display: block;
    height: 6px;
    border-radius: 99px;
    background: var(--accent-pd);
    min-width: 4px;
  }

  .effect-row > i.down {
    background: var(--accent-ai);
  }

  .effect-row > i.stable {
    background: var(--text-muted);
  }

  .integration-result {
    margin-top: var(--space-6);
    padding: var(--space-5);
    border: 1px solid var(--border-strong);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .integration-head {
    display: flex;
    justify-content: space-between;
    gap: var(--space-4);
    align-items: start;
  }

  .integration-head > span {
    font-family: var(--font-mono);
    font-size: var(--text-xs);
    color: var(--text-muted);
  }

  .module-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-4);
    margin-top: var(--space-4);
  }

  .module-grid article {
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-primary);
  }

  .module-grid article.primary {
    border-color: var(--accent-pd);
    border-width: 2px;
  }

  .module-strength {
    font-family: var(--font-mono);
    font-size: 10px;
    color: var(--accent-pk);
    text-transform: uppercase;
  }

  .module-evidence {
    display: flex;
    flex-wrap: wrap;
    gap: 5px;
    margin: var(--space-3) 0;
  }

  .module-evidence span {
    font-family: var(--font-mono);
    font-size: 10px;
    padding: 4px 6px;
    border: 1px solid var(--border-subtle);
    border-radius: 999px;
  }

  .mechanism-demo {
    margin-top: var(--space-5);
    padding: var(--space-5);
    border-left: 3px solid var(--accent-pd);
    background: var(--bg-secondary);
  }

  .mechanism-flow {
    display: flex;
    flex-wrap: wrap;
    gap: 10px;
    align-items: center;
  }

  .mechanism-flow > div {
    padding: var(--space-3);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-primary);
  }

  .mechanism-flow small,
  .mechanism-flow strong {
    display: block;
  }

  .mechanism-flow small {
    color: var(--text-muted);
  }

  .mechanism-flow > b {
    color: var(--accent-pk);
  }

  .evidence-layers {
    display: grid;
    grid-template-columns: repeat(4, 1fr);
    gap: var(--space-3);
    margin-top: var(--space-5);
  }

  .evidence-layers article {
    padding: var(--space-3);
    border-top: 2px solid var(--accent-pd);
    background: var(--bg-secondary);
  }

  .evidence-layers article.future {
    border-top-color: var(--text-muted);
  }

  .evidence-layers p {
    color: var(--text-secondary);
    font-size: var(--text-sm);
    margin-bottom: 0;
  }

  .plan {
    display: grid;
    grid-template-columns: .9fr 1.1fr;
    gap: var(--space-6);
    padding: var(--space-5);
    background: var(--bg-secondary);
    border-left: 3px solid var(--accent-pk);
  }

  .plan li {
    margin-bottom: var(--space-2);
    color: var(--text-secondary);
  }

  .database-bridge {
    margin-top: var(--space-6);
  }

  .section-head.compact {
    align-items: start;
    margin-bottom: var(--space-4);
  }

  .database-grid {
    display: grid;
    grid-template-columns: repeat(3, 1fr);
    gap: var(--space-3);
  }

  .database-grid article {
    padding: var(--space-4);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .database-grid article > span {
    font-family: var(--font-mono);
    font-size: 10px;
    color: var(--accent-pk);
    text-transform: uppercase;
  }

  .database-grid article p {
    color: var(--text-secondary);
    font-size: var(--text-sm);
  }

  .database-grid article a {
    font-size: var(--text-xs);
  }

  .resolution-flow {
    display: flex;
    flex-wrap: wrap;
    align-items: center;
    gap: 8px;
    margin-top: var(--space-4);
    padding: var(--space-4);
    background: var(--bg-secondary);
    border-left: 3px solid var(--accent-pd);
  }

  .resolution-flow code {
    border: 1px solid var(--border-strong);
    border-radius: 5px;
    padding: 4px 7px;
    background: var(--bg-primary);
  }

  .result-grid {
    display: grid;
    grid-template-columns: repeat(3, minmax(0, 1fr));
    gap: var(--space-4);
  }

  .result-grid article {
    min-height: 240px;
    padding: var(--space-5);
    border: 1px solid var(--border-subtle);
    border-radius: var(--radius);
    background: var(--bg-secondary);
  }

  .result-grid article p {
    color: var(--text-secondary);
  }

  .num {
    display: inline-block;
    margin-bottom: var(--space-5);
    font-family: var(--font-mono);
    color: var(--accent-pk);
  }

  .mini-factor {
    display: grid;
    gap: 6px;
    margin-top: var(--space-5);
  }

  .mini-factor i {
    display: block;
    height: 7px;
    background: var(--accent-pk);
    border-radius: 10px;
  }

  .mini-factor i:nth-child(2) {
    width: 72%;
    opacity: .65;
  }

  .mini-factor i:nth-child(3) {
    width: 44%;
    opacity: .4;
  }

  .layers,
  .evidence {
    display: flex;
    flex-wrap: wrap;
    gap: 6px;
    margin-top: var(--space-4);
  }

  .layers span,
  .evidence span {
    font-family: var(--font-mono);
    font-size: var(--text-xs);
    border: 1px solid var(--border-strong);
    padding: 4px 7px;
    border-radius: 999px;
  }

  .network-demo,
  .association {
    display: flex;
    align-items: center;
    justify-content: space-between;
    gap: 8px;
    margin-top: var(--space-6);
    font-family: var(--font-mono);
    font-size: var(--text-xs);
  }

  .network-demo b,
  .association span {
    border: 1px solid var(--border-strong);
    padding: 7px;
    border-radius: 5px;
    background: var(--bg-primary);
  }

  .network-demo em {
    font-style: normal;
    color: var(--accent-pk);
  }

  table {
    width: 100%;
    border-collapse: collapse;
    font-family: var(--font-mono);
    font-size: var(--text-xs);
  }

  th,
  td {
    padding: 5px;
    border-bottom: 1px solid var(--border-subtle);
    text-align: left;
  }

  .muted {
    color: var(--text-muted);
  }

  @media (max-width: 1000px) {
    .workflow, .result-grid { grid-template-columns: repeat(2, 1fr); }
    .contract, .plan { grid-template-columns: 1fr; }
    .alias-grid, .mapping-grid, .matrix-checks, .identity-grid, .omics-question-grid, .database-grid, .demo-story, .demo-omics-grid, .module-grid, .evidence-layers, .actual-layer-grid, .computed-summary, .interpretation-grid, .qc-grid { grid-template-columns: repeat(2, 1fr); }
    .validation { grid-template-columns: repeat(2, 1fr); }
  }

  @media (max-width: 640px) {
    .section-head, .mapping-head { align-items: start; flex-direction: column; }
    .form-grid, .uploads, .result-grid, .workflow, .alias-grid, .mapping-grid, .matrix-checks, .identity-grid, .validation, .omics-question-grid, .database-grid, .demo-story, .demo-omics-grid, .module-grid, .evidence-layers, .actual-layer-grid, .computed-summary, .interpretation-grid, .qc-grid { grid-template-columns: 1fr; }
    .run-box, .overlap-box { align-items: stretch; flex-direction: column; }
    .backend-controls, .backend-method-grid, .ms-qc-controls { grid-template-columns: 1fr; }
    .overlap-pairs { justify-content: flex-start; }
    .dictionary-table > div { grid-template-columns: 1fr; gap: 2px; padding: 12px 0; }
    .workflow div { border-right: 0; border-bottom: 1px solid var(--border-subtle); }
    .workflow div:last-child { border-bottom: 0; }
    .wide { grid-column: auto; }
    .status { align-items: start; flex-direction: column; }
  }
  /* Use the actual light/dark palette. --surface was never defined: in dark
     mode its #fff fallback inherited a near-white font (WCAG failure). */
  .scientific-assurance { margin:16px 0 24px; padding:16px 18px; background:var(--bg-secondary); border:1px solid var(--border-strong); border-inline-start:4px solid var(--accent-pk); color:var(--text-primary); border-radius:8px; }
  .scientific-assurance strong { display:block; margin-bottom:6px; color:var(--text-primary); }
  .scientific-assurance p { margin:4px 0 12px; color:var(--text-secondary); font-size:.92rem; }
  .scientific-assurance summary { color:var(--text-primary); font-weight:650; cursor:pointer; }
  .scientific-assurance li { padding:3px 0; color:var(--text-secondary); font-size:.89rem; }
  .quick-start { display:flex; align-items:center; justify-content:space-between; flex-wrap:wrap; gap:16px; margin:20px 0; padding:18px 20px; border:1px solid var(--border-strong); border-radius:14px; background:var(--bg-secondary); color:var(--text-primary); }
  .quick-start-copy { display:grid; gap:6px; flex:1 1 275px; min-width:0; }
  .quick-start strong { color:var(--text-primary); font-size:1rem; line-height:1.4; }
  .quick-start span { font-size:.91rem; line-height:1.55; color:var(--text-secondary); }
  .quick-start small { font-size:.8rem; line-height:1.5; color:var(--text-muted); }
  .quick-demo-actions { display:grid; grid-template-columns:repeat(3,minmax(0,1fr)); gap:8px; width:min(100%,485px); }
  .quick-demo-actions .btn { min-width:0; text-align:center; white-space:normal; overflow-wrap:anywhere; }
  .quick-demo-actions .quick-demo-single { font-size:.82rem; padding:10px 6px; }
  .quick-demo-actions .quick-demo-multi { grid-column:1 / -1; width:100%; min-height:48px; font-weight:700; }
  .optional-workflows { margin-top:16px; border:1px solid var(--border-strong); border-radius:12px; padding:12px 16px; }
  .optional-workflows > summary { cursor:pointer; font-weight:700; }
  .optional-workflows > .contract { margin-top:16px; }
  .annotation-error { color:var(--quiz-error-text); font-weight:700; }
  .ms-import-message { padding:10px 13px; border:1px solid var(--border-strong); border-radius:9px; color:var(--text-primary); background:var(--bg-secondary); font-size:.88rem; line-height:1.55; }
  @media (max-width:640px) {
    .quick-start { align-items:stretch; }
    .quick-demo-actions { width:100%; }
  }
</style>
