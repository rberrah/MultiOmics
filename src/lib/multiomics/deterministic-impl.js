// @ts-nocheck
// Public wrapper around the stable deterministic engine.
// The core remains isolated in deterministic-core.js; this layer adds
// conservative chemical resolution and study-feasibility/methodology diagnostics.

import * as core from './deterministic-core.js';
export * from './deterministic-core.js';

const WRAPPER_ENGINE_VERSION = '1.8.0';
const CHEBI_SEARCH_HOST = 'www.ebi.ac.uk';
const CHEBI_SEARCH_PATH = '/chebi/backend/api/public/es_search/';
const UNICHEM_BASE = 'https://www.ebi.ac.uk/unichem/rest/src_compound_id/';

// Current UniChem source identifiers used here: ChEBI = 7, HMDB = 18,
// PubChem Compounds = 22. KEGG is intentionally not routed through UniChem:
// it is no longer listed in the current UniChem source registry.
const UNICHEM_SOURCE = Object.freeze({
  chebi: 7,
  hmdb: 18,
  pubchem: 22
});

function cleanText(value) {
  return String(value ?? '').trim();
}

function traceKey(value) {
  return cleanText(value).toUpperCase();
}

function isKeggCompound(term, identifierType = 'unknown') {
  const value = cleanText(term);
  const type = cleanText(identifierType).toLowerCase();
  return /^C\d{5}$/i.test(value) && ['unknown','kegg','kegg_compound'].includes(type);
}

function declaredSource(term, identifierType = 'unknown') {
  const value = cleanText(term);
  const type = cleanText(identifierType).toLowerCase();
  if (/^HMDB\d+$/i.test(value)) return { name: 'hmdb', srcId: UNICHEM_SOURCE.hmdb };
  if ((type === 'pubchem' || type === 'pubchem_cid') && /^\d+$/.test(value)) {
    return { name: 'pubchem', srcId: UNICHEM_SOURCE.pubchem };
  }
  return null;
}

function canonicalChebi(value) {
  const text = cleanText(value);
  if (/^CHEBI:\d+$/i.test(text)) return text.toUpperCase();
  if (/^\d+$/.test(text)) return 'CHEBI:' + text;
  return null;
}

function uniqueChebiIds(json) {
  const rows = Array.isArray(json) ? json : [];
  const values = rows.flatMap((row) => {
    if (typeof row === 'string' || typeof row === 'number') return [row];
    const value = row?.src_compound_id ?? row?.srcCompoundId ?? row?.compound_id ?? row?.id;
    return Array.isArray(value) ? value : [value];
  });
  return [...new Set(values.map(canonicalChebi).filter(Boolean))];
}

function uniqueChebiSearchIds(json) {
  const rows = Array.isArray(json?.results) ? json.results : [];
  return [...new Set(rows
    .map((row) => row?._source || row?.source || row || null)
    .map((source) => canonicalChebi(source?.chebi_accession))
    .filter(Boolean))];
}

function syntheticChebiResponse(term, chebiIds) {
  const results = chebiIds.map((chebi) => ({
    _source: {
      chebi_accession: chebi,
      name: term,
      ascii_name: term,
      stars: null
    }
  }));
  return new Response(JSON.stringify({ results, total: results.length }), {
    status: 200,
    headers: { 'Content-Type': 'application/json' }
  });
}

function requestUrl(input) {
  try {
    if (typeof input === 'string' || input instanceof URL) return new URL(String(input));
    if (typeof Request !== 'undefined' && input instanceof Request) return new URL(input.url);
    if (input?.url) return new URL(String(input.url));
  } catch {
    return null;
  }
  return null;
}

function resolutionMethod(sourceName) {
  return `unichem_${sourceName}_to_chebi`;
}

function createStrictChemicalFetch(baseFetch, identifierType = 'unknown') {
  const trace = new Map();
  const safeFetch = typeof baseFetch === 'function'
    ? baseFetch
    : (...args) => globalThis.fetch(...args);

  const fetchFn = async (input, init) => {
    const url = requestUrl(input);
    if (!url || url.hostname !== CHEBI_SEARCH_HOST || url.pathname !== CHEBI_SEARCH_PATH) {
      return safeFetch(input, init);
    }

    const term = cleanText(url.searchParams.get('term'));
    const key = traceKey(term);

    // ChEBI 2.0 can search curated external cross-references. KEGG identifiers
    // are accepted only when one unique ChEBI accession is returned.
    if (isKeggCompound(term, identifierType)) {
      const method = 'chebi_exact_kegg_xref_search';
      try {
        const response = await safeFetch(input, init);
        if (!response.ok) {
          trace.set(key, {
            status: response.status === 404 ? 'unresolved' : 'api_error',
            method,
            candidates: [],
            error: response.status === 404 ? null : `HTTP ${response.status}`
          });
          return syntheticChebiResponse(term, []);
        }
        const ids = uniqueChebiSearchIds(await response.json());
        if (ids.length === 1) {
          trace.set(key, { status: 'resolved_external_id', method, candidates: ids, resolved: ids[0] });
          return syntheticChebiResponse(term, ids);
        }
        if (ids.length > 1) {
          trace.set(key, { status: 'ambiguous', method, candidates: ids, resolved: null });
          return syntheticChebiResponse(term, []);
        }
        trace.set(key, { status: 'unresolved', method, candidates: [], resolved: null });
        return syntheticChebiResponse(term, []);
      } catch (error) {
        trace.set(key, {
          status: 'api_error',
          method,
          candidates: [],
          resolved: null,
          error: error instanceof Error ? error.message : 'ChEBI KEGG cross-reference search failed'
        });
        return syntheticChebiResponse(term, []);
      }
    }

    const source = declaredSource(term, identifierType);
    if (!source) return safeFetch(input, init);

    const method = resolutionMethod(source.name);
    try {
      const mappingUrl = UNICHEM_BASE
        + encodeURIComponent(term)
        + '/' + source.srcId
        + '/' + UNICHEM_SOURCE.chebi;
      const response = await safeFetch(mappingUrl, { headers: { Accept: 'application/json' } });
      if (!response.ok) {
        trace.set(key, {
          status: response.status === 404 ? 'unresolved' : 'api_error',
          method,
          candidates: [],
          error: response.status === 404 ? null : `HTTP ${response.status}`
        });
        return syntheticChebiResponse(term, []);
      }

      const ids = uniqueChebiIds(await response.json());
      if (ids.length === 1) {
        trace.set(key, { status: 'resolved_external_id', method, candidates: ids, resolved: ids[0] });
        return syntheticChebiResponse(term, ids);
      }
      if (ids.length > 1) {
        trace.set(key, { status: 'ambiguous', method, candidates: ids, resolved: null });
        return syntheticChebiResponse(term, []);
      }
      trace.set(key, { status: 'unresolved', method, candidates: [], resolved: null });
      return syntheticChebiResponse(term, []);
    } catch (error) {
      trace.set(key, {
        status: 'api_error',
        method,
        candidates: [],
        resolved: null,
        error: error instanceof Error ? error.message : 'UniChem request failed'
      });
      // Do not fall back to a fuzzy name search after a failed external-ID lookup.
      return syntheticChebiResponse(term, []);
    }
  };

  return { fetchFn, trace };
}

function applyChemicalTrace(resolution, trace) {
  if (!resolution?.mappings || !(trace instanceof Map) || trace.size === 0) return resolution;
  const mappings = resolution.mappings.map((mapping) => {
    const traced = trace.get(traceKey(mapping.original));
    if (!traced) return mapping;
    return {
      ...mapping,
      resolved: traced.resolved || null,
      status: traced.status,
      method: traced.method,
      candidates: traced.candidates || [],
      ...(traced.error ? { error: traced.error } : {})
    };
  });
  return {
    ...resolution,
    mappings,
    resolvedCount: mappings.filter((item) => item.resolved).length,
    unresolvedCount: mappings.filter((item) => !item.resolved).length,
    ambiguousCount: mappings.filter((item) => item.status === 'ambiguous').length
  };
}

function mappedValue(row, mapping, key) {
  const source = cleanText(mapping?.[key]);
  return source && Object.prototype.hasOwnProperty.call(row || {}, source)
    ? cleanText(row[source])
    : '';
}

function uniqueNonEmpty(values) {
  return [...new Set(values.map(cleanText).filter(Boolean))];
}

function levelFromRatio(ratio, good = 0.8, caution = 0.5) {
  if (!Number.isFinite(ratio)) return 'not_applicable';
  if (ratio >= good) return 'good';
  if (ratio >= caution) return 'caution';
  return 'warning';
}

function canonicalBooleanEvent(value) {
  const text = cleanText(value).toLowerCase();
  if (!text) return null;
  if (['1','true','yes','y','event','dead','death','oui','événement','evenement'].includes(text)) return 1;
  if (['0','false','no','n','censored','censor','alive','non','censuré','censure'].includes(text)) return 0;
  const numeric = Number(text.replace(',', '.'));
  if (Number.isFinite(numeric) && (numeric === 0 || numeric === 1)) return numeric;
  return null;
}

function subjectsPerValue(rows, valueKey) {
  const output = {};
  for (const value of uniqueNonEmpty(rows.map((row) => row[valueKey]))) {
    output[value] = new Set(rows.filter((row) => row[valueKey] === value).map((row) => row.subject).filter(Boolean)).size;
  }
  return output;
}

function publishDiagnostics(diagnostics) {
  if (typeof window === 'undefined' || typeof window.dispatchEvent !== 'function') return;
  window.__PMX_MULTIOMICS_DIAGNOSTICS__ = diagnostics;
  try {
    window.dispatchEvent(new CustomEvent('pmx-multiomics-diagnostics', { detail: diagnostics }));
  } catch {
    // The analysis result still contains the diagnostics even if a legacy
    // browser does not expose CustomEvent as expected.
  }
}

function buildMethodEligibility({ protocol, canonicalRows, perLayer, overlap, conditions, groupCounts, blockers, warnings }) {
  const objective = protocol.objective || 'explore';
  const outcomeType = protocol.outcomeType || 'none';
  const commonSubjects = Number(overlap?.allMatched ?? 0);
  const singleOmic = Object.keys(perLayer).length === 1;
  const batchBlocked = Object.values(perLayer).some((item) => item.batchStatus === 'confounded');
  const batchIncomplete = Object.values(perLayer).some((item) => item.batchStatus === 'incomplete');

  const mofaReasons = [];
  if (singleOmic) mofaReasons.push('MOFA2 multi-omics factor integration requires at least two measured layers.');
  if (commonSubjects < 16) mofaReasons.push('MOFA2 factor analysis requires more than 15 shared samples for a defensible primary analysis.');
  if (batchBlocked) mofaReasons.push('Known technical series are confounded with the biological design.');
  if (batchIncomplete) mofaReasons.push('Technical-series annotation is incomplete.');
  const mofaStatus = singleOmic ? 'not_applicable' : mofaReasons.length ? 'not_recommended' : 'eligible';

  const categoricalTarget = objective === 'groups' || (objective === 'outcome' && ['binary','multiclass'].includes(outcomeType));
  let classCounts = groupCounts;
  if (objective === 'outcome' && ['binary','multiclass'].includes(outcomeType)) {
    classCounts = subjectsPerValue(canonicalRows, 'outcome');
  }
  const classSizes = Object.values(classCounts).filter((value) => value > 0);
  const smallestClass = classSizes.length ? Math.min(...classSizes) : null;
  const diabloReasons = [];
  if (singleOmic) diabloReasons.push('DIABLO multi-omics integration requires at least two measured layers.');
  if (!categoricalTarget) diabloReasons.push('DIABLO is only applicable to a supervised categorical target in this tool.');
  if (commonSubjects < 6) diabloReasons.push('Too few subjects are shared across all loaded omics blocks.');
  if (smallestClass != null && smallestClass < 3) diabloReasons.push('At least three subjects per class are required for stratified cross-validation.');
  if (batchBlocked) diabloReasons.push('Technical confounding must be resolved before supervised multi-omics discrimination.');
  const diabloStatus = !categoricalTarget || singleOmic ? 'not_applicable' : diabloReasons.length ? 'not_recommended' : 'eligible_with_internal_cv';

  return {
    mofa2: {
      status: mofaStatus,
      sharedSubjects: commonSubjects,
      reasons: mofaReasons,
      interpretation: 'Latent factors are exploratory covariance structures, not causal mechanisms or validated biomarkers.'
    },
    diablo: {
      status: diabloStatus,
      sharedSubjects: commonSubjects,
      smallestClass,
      classSizes: classCounts,
      reasons: diabloReasons,
      interpretation: 'A tuned DIABLO signature requires repeated internal CV and independent external validation before biomarker claims.'
    },
    global: {
      blockers: blockers.length,
      warnings: warnings.length
    }
  };
}

function buildPreAnalysisDiagnostics(args, result) {
  const rows = Array.isArray(args?.metadataRows) ? args.metadataRows : [];
  const mapping = args?.columnMapping || {};
  const protocol = args?.protocol || {};
  const biologicalRows = rows.filter((row) => {
    const sampleType = mappedValue(row, mapping, 'sample_type').toLowerCase();
    return !['blank','qc','pooled_qc','quality_control','pooled'].includes(sampleType);
  });

  const canonicalRows = biologicalRows.map((row) => ({
    subject: mappedValue(row, mapping, 'subject_id'),
    sample: mappedValue(row, mapping, 'sample_id'),
    assay: mappedValue(row, mapping, 'assay_id'),
    omic: mappedValue(row, mapping, 'omic').toLowerCase(),
    condition: mappedValue(row, mapping, 'condition'),
    timepoint: mappedValue(row, mapping, 'timepoint'),
    batch: mappedValue(row, mapping, 'batch'),
    outcome: mappedValue(row, mapping, 'outcome'),
    survivalTime: mappedValue(row, mapping, 'survival_time'),
    survivalEvent: mappedValue(row, mapping, 'survival_event')
  })).filter((row) => row.subject || row.sample || row.assay)
    .filter((row) => {
      // The sample sheet may describe other, non-uploaded modalities. Diagnostics
      // must use the same cohort and layers as the actual statistical analysis.
      const analysedLayers = result?.engine?.analysedLayers || Object.keys(result?.layers || {});
      return analysedLayers.includes(row.omic);
    });

  const conditions = uniqueNonEmpty(canonicalRows.map((row) => row.condition));
  const timepoints = uniqueNonEmpty(canonicalRows.map((row) => row.timepoint));
  const groupCounts = subjectsPerValue(canonicalRows, 'condition');
  const positiveGroupCounts = Object.values(groupCounts).filter((value) => value > 0);
  const minGroup = positiveGroupCounts.length ? Math.min(...positiveGroupCounts) : null;
  const maxGroup = positiveGroupCounts.length ? Math.max(...positiveGroupCounts) : null;
  const balanceRatio = minGroup != null && maxGroup ? minGroup / maxGroup : null;
  const groupBalanceStatus = conditions.length >= 2 ? levelFromRatio(balanceRatio) : 'not_applicable';

  const overlap = result?.metadataSummary?.overlap || {};
  const overlapRatio = overlap.totalSubjectsAnyLayer > 0
    ? overlap.allMatched / overlap.totalSubjectsAnyLayer
    : null;
  let overlapStatus = levelFromRatio(overlapRatio);
  if (protocol.partialOmicsExpected === 'yes' && overlapStatus === 'warning') overlapStatus = 'expected_partial';

  const perLayer = {};
  for (const [layer, layerResult] of Object.entries(result?.layers || {})) {
    const qc = layerResult?.qc || {};
    const medianMissing = Number.isFinite(qc.medianMissingFraction) ? qc.medianMissingFraction : null;
    const missingStatus = medianMissing == null
      ? 'not_available'
      : medianMissing <= 0.10 ? 'good' : medianMissing <= 0.30 ? 'caution' : 'warning';

    const layerRows = canonicalRows.filter((row) => row.omic === layer);
    const subjects = uniqueNonEmpty(layerRows.map((row) => row.subject));
    const observedTimes = uniqueNonEmpty(layerRows.map((row) => row.timepoint));
    let completeLongitudinalSubjects = null;
    let longitudinalCompleteness = null;
    let subjectsWithRepeatedMeasures = null;
    if (timepoints.length > 1 && subjects.length) {
      completeLongitudinalSubjects = subjects.filter((subject) => {
        const subjectTimes = new Set(layerRows.filter((row) => row.subject === subject).map((row) => row.timepoint).filter(Boolean));
        return timepoints.every((time) => subjectTimes.has(time));
      }).length;
      longitudinalCompleteness = completeLongitudinalSubjects / subjects.length;
      subjectsWithRepeatedMeasures = subjects.filter((subject) => {
        const subjectTimes = new Set(layerRows.filter((row) => row.subject === subject).map((row) => row.timepoint).filter(Boolean));
        return subjectTimes.size >= 2;
      }).length;
    }

    const batch = result?.metadataSummary?.batchAudit?.[layer] || null;
    perLayer[layer] = {
      subjects: subjects.length,
      observedTimepoints: observedTimes,
      medianMissingFraction: medianMissing,
      missingnessStatus: missingStatus,
      flaggedAssays: qc.outlierSamples?.length || 0,
      weakTechnicalReplicates: (qc.replicateCorrelations || []).filter((item) => item?.warning).length,
      completeLongitudinalSubjects,
      longitudinalCompleteness,
      subjectsWithRepeatedMeasures,
      batchStatus: batch?.status || 'not_available',
      batchCount: batch?.batchCount ?? batch?.batches?.length ?? 0,
      batchNote: batch?.note || ''
    };
  }

  const blockers = [];
  const warnings = [];
  const notes = [];
  const objective = protocol.objective || 'explore';
  const outcomeType = protocol.outcomeType || 'none';

  if (objective === 'groups' && conditions.length < 2) {
    blockers.push('Group comparison was requested but fewer than two biological conditions are represented.');
  }
  if (objective === 'groups' && minGroup != null && minGroup < 2) {
    blockers.push('At least one biological group has fewer than two independent subjects; feature-wise group inference is not estimable defensibly.');
  } else if (objective === 'groups' && minGroup != null && minGroup < 5) {
    warnings.push('At least one biological group has fewer than five independent subjects; treat inferential results as exploratory and emphasise effect sizes/uncertainty.');
  }

  if (conditions.length >= 2 && balanceRatio != null && balanceRatio < 0.5) {
    warnings.push('One biological group contains less than half as many unique subjects as the largest group.');
  } else if (conditions.length >= 2 && balanceRatio != null && balanceRatio < 0.8) {
    notes.push('Group sizes are moderately unbalanced; adjusted models remain possible but precision may differ between groups.');
  }

  if (overlapRatio != null && overlapRatio < 0.5 && protocol.partialOmicsExpected !== 'yes') {
    warnings.push('Fewer than half of subjects are represented in every loaded omics layer.');
  } else if (overlapRatio != null && overlapRatio < 0.8 && protocol.partialOmicsExpected !== 'yes') {
    notes.push('Not all subjects are represented in every omics layer; paired cross-omics analyses use the actual overlap only.');
  }

  for (const [layer, diag] of Object.entries(perLayer)) {
    if (diag.missingnessStatus === 'warning') warnings.push(`${layer}: median missingness exceeds 30%.`);
    else if (diag.missingnessStatus === 'caution') notes.push(`${layer}: median missingness is between 10% and 30%.`);
    if (diag.longitudinalCompleteness != null && diag.longitudinalCompleteness < 0.8) {
      notes.push(`${layer}: ${Math.round(100 * diag.longitudinalCompleteness)}% of subjects have all declared time points.`);
    }
    if ((objective === 'time' || protocol.longitudinal || protocol.designType === 'repeated') &&
        diag.subjectsWithRepeatedMeasures != null && diag.subjectsWithRepeatedMeasures < 2) {
      blockers.push(`${layer}: fewer than two subjects contain repeated time points, so a repeated-measures population effect is not defensibly estimable.`);
    }
  }

  const batchProblems = Object.entries(result?.metadataSummary?.batchAudit || {});
  const confounded = batchProblems.filter(([, audit]) => audit?.status === 'confounded');
  const incomplete = batchProblems.filter(([, audit]) => audit?.status === 'incomplete');
  if (confounded.length) {
    blockers.push('At least one technical series is completely confounded with the biological design; the corresponding biological effect cannot be separated statistically.');
  }
  if (incomplete.length) {
    warnings.push('Technical-series annotation is incomplete in at least one omics layer.');
  }

  let survival = null;
  if (objective === 'outcome' && outcomeType === 'survival') {
    const bySubject = new Map();
    for (const row of canonicalRows) {
      if (!row.subject) continue;
      const event = canonicalBooleanEvent(row.survivalEvent);
      const time = Number(cleanText(row.survivalTime).replace(',', '.'));
      if (!bySubject.has(row.subject)) bySubject.set(row.subject, { event: null, time: null });
      const entry = bySubject.get(row.subject);
      if (event != null) entry.event = event;
      if (Number.isFinite(time) && time >= 0) entry.time = time;
    }
    const usable = [...bySubject.values()].filter((entry) => entry.event != null && Number.isFinite(entry.time));
    const events = usable.filter((entry) => entry.event === 1).length;
    const censored = usable.filter((entry) => entry.event === 0).length;
    survival = { subjectsWithUsableSurvival: usable.length, events, censored };
    if (usable.length < 4 || events < 2) {
      blockers.push('Survival analysis requires more usable subjects and at least two observed events.');
    } else if (events < 10) {
      warnings.push('Fewer than ten observed survival events are available; Cox associations and internal prediction estimates may be highly unstable.');
    }
  }

  const smallestGroup = minGroup;
  const power = {
    status: 'not_calculated',
    smallestGroupSubjects: smallestGroup,
    note: 'Statistical power is not inferred from sample count alone. A defensible power calculation requires an expected effect size, variability and the intended model; simulation is recommended for complex multi-omics designs.'
  };

  const methodEligibility = buildMethodEligibility({
    protocol,
    canonicalRows,
    perLayer,
    overlap,
    conditions,
    groupCounts,
    blockers,
    warnings
  });

  const overallStatus = blockers.length
    ? 'blocked_or_redesign_required'
    : warnings.length
      ? 'review_required'
      : notes.length
        ? 'usable_with_cautions'
        : 'ready';

  return {
    status: overallStatus,
    plainLanguageTitle: 'Study check before interpretation',
    groupBalance: {
      conditions,
      subjectsPerGroup: groupCounts,
      smallestGroup: minGroup,
      largestGroup: maxGroup,
      smallestToLargestRatio: balanceRatio,
      status: groupBalanceStatus
    },
    omicsOverlap: {
      subjectsInEveryLoadedLayer: overlap.allMatched ?? null,
      subjectsInAnyLoadedLayer: overlap.totalSubjectsAnyLayer ?? null,
      fractionInEveryLoadedLayer: overlapRatio,
      pairwise: overlap.pairwise || [],
      status: overlapStatus,
      partialOmicsExpected: protocol.partialOmicsExpected === 'yes'
    },
    repeatedMeasures: {
      declaredTimepoints: timepoints,
      longitudinal: Boolean(protocol.longitudinal || objective === 'time' || protocol.designType === 'repeated'),
      perLayer: Object.fromEntries(Object.entries(perLayer).map(([layer, diag]) => [layer, {
        completeSubjects: diag.completeLongitudinalSubjects,
        fractionComplete: diag.longitudinalCompleteness,
        subjectsWithRepeatedMeasures: diag.subjectsWithRepeatedMeasures
      }]))
    },
    survival,
    perLayer,
    power,
    methodEligibility,
    blockers,
    warnings,
    notes,
    interpretation: {
      ready: 'No deterministic design/QC warning was detected. This does not prove adequate statistical power or external validity.',
      usable_with_cautions: 'The analysis can run, but the listed cautions should be considered when interpreting effect sizes, uncertainty and cross-omics integration.',
      review_required: 'At least one design/QC issue requires explicit review before biological conclusions are emphasised.',
      blocked_or_redesign_required: 'At least one deterministic design issue prevents a defensible strong interpretation. Correct the design/annotation or restrict the analysis to a descriptive route.'
    }[overallStatus]
  };
}

export async function resolveMetaboliteIdentifier(
  identifier,
  { fetchFn = globalThis.fetch, identifierType = 'unknown' } = {}
) {
  const strict = createStrictChemicalFetch(fetchFn, identifierType);
  const result = await core.resolveMetaboliteIdentifier(identifier, { fetchFn: strict.fetchFn });
  const traced = strict.trace.get(traceKey(identifier));
  if (!traced) return result;
  return {
    ...result,
    resolved: traced.resolved || null,
    status: traced.status,
    method: traced.method,
    candidates: traced.candidates || [],
    ...(traced.error ? { error: traced.error } : {})
  };
}

export async function resolveMetaboliteIdentifiers(
  identifiers,
  { fetchFn = globalThis.fetch, maxQueries = 30, identifierType = 'unknown' } = {}
) {
  const strict = createStrictChemicalFetch(fetchFn, identifierType);
  const result = await core.resolveMetaboliteIdentifiers(identifiers, {
    fetchFn: strict.fetchFn,
    maxQueries
  });
  return applyChemicalTrace(result, strict.trace);
}

let analysisQueue = Promise.resolve();

async function runWithStrictChemicalResolution(args) {
  const originalFetch = globalThis.fetch;
  if (typeof originalFetch !== 'function') {
    const result = await core.runDeterministicAnalysis(args);
    const diagnostics = buildPreAnalysisDiagnostics(args, result);
    result.preAnalysisDiagnostics = diagnostics;
    result.metadataSummary = { ...result.metadataSummary, preAnalysisDiagnostics: diagnostics };
    publishDiagnostics(diagnostics);
    if (result?.engine) {
      result.engine = {
        ...result.engine,
        version: WRAPPER_ENGINE_VERSION,
        feasibilityDiagnostics: 'deterministic blocking/caution audit for biological independence, group balance, batch confounding, cross-omics overlap, longitudinal support, missingness, survival events and advanced-method eligibility'
      };
    }
    return result;
  }

  const identifierType = args?.identifierTypes?.metabolomics || 'unknown';
  const strict = createStrictChemicalFetch(originalFetch.bind(globalThis), identifierType);
  globalThis.fetch = strict.fetchFn;
  try {
    const result = await core.runDeterministicAnalysis(args);
    if (result?.identifierResolution?.metabolomics) {
      result.identifierResolution.metabolomics = applyChemicalTrace(
        result.identifierResolution.metabolomics,
        strict.trace
      );
    }
    const diagnostics = buildPreAnalysisDiagnostics(args, result);
    result.preAnalysisDiagnostics = diagnostics;
    result.metadataSummary = { ...result.metadataSummary, preAnalysisDiagnostics: diagnostics };
    publishDiagnostics(diagnostics);
    if (result?.engine) {
      result.engine = {
        ...result.engine,
        version: WRAPPER_ENGINE_VERSION,
        chemicalIdentifierResolution: 'strict HMDB/PubChem UniChem mapping plus unique ChEBI KEGG cross-reference search; ambiguous mappings are never forced',
        feasibilityDiagnostics: 'deterministic blocking/caution audit for biological independence, group balance, batch confounding, cross-omics overlap, longitudinal support, missingness, survival events and advanced-method eligibility',
        interpretationPolicy: 'association, prediction and causality are reported as distinct claims; advanced multivariate methods are eligibility-gated'
      };
    }
    return result;
  } finally {
    globalThis.fetch = originalFetch;
  }
}

export function runDeterministicAnalysis(args) {
  const run = analysisQueue.then(
    () => runWithStrictChemicalResolution(args),
    () => runWithStrictChemicalResolution(args)
  );
  analysisQueue = run.catch(() => undefined);
  return run;
}
