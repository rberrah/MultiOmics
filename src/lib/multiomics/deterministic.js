// @ts-nocheck
// Public multi-omics entry point.
// The validated deterministic implementation lives in deterministic-impl.js.
// This thin browser-aware wrapper adds declared UI methodology context,
// external database provenance, cryptographic input provenance,
// chemical-mapping audit metadata and publishes the final result synchronously
// before it is returned to Svelte.

import * as impl from './deterministic-impl.js';
export * from './deterministic-impl.js';

const REACTOME_VERSION_URL = 'https://reactome.org/ContentService/data/database/version';
const REACTOME_ANALYSIS_SERVICE = 'https://reactome.org/AnalysisService/';
const CHEMICAL_MAPPING_POLICY = 'Canonicalise to ChEBI only after a canonical ChEBI input, one unique curated external-ID mapping, one unique KEGG cross-reference hit, or one exact ChEBI name match. Ambiguous mappings are never forced.';

function withBrowserMethodology(args) {
  if (typeof window === 'undefined') return args;
  const level = window.__PMX_METABOLOMICS_IDENTIFICATION_CONFIDENCE__;
  if (!level) return args;
  return {
    ...args,
    protocol: {
      ...(args?.protocol || {}),
      metabolomicsIdentificationConfidence: level
    }
  };
}

function cleanText(value) {
  return String(value ?? '').trim();
}

function inferMetaboliteSource(identifier, identifierType = 'unknown') {
  const value = cleanText(identifier);
  const declared = cleanText(identifierType).toLowerCase();
  if (/^CHEBI:\d+$/i.test(value) || declared === 'chebi') return 'ChEBI';
  if (/^HMDB\d+$/i.test(value) || declared === 'hmdb') return 'HMDB';
  if (/^C\d{5}$/i.test(value) || ['kegg','kegg_compound'].includes(declared)) return 'KEGG Compound';
  if ((['pubchem','pubchem_cid'].includes(declared) && /^\d+$/.test(value))) return 'PubChem Compound';
  if (/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/i.test(value) || ['inchikey','inchi_key'].includes(declared)) return 'InChIKey';
  if (['name','metabolite_name'].includes(declared)) return 'metabolite name';
  return 'unknown / name search';
}

function chemicalResolutionConfidence(mapping) {
  const status = mapping?.status;
  if (status === 'canonical') return 'canonical_input';
  if (status === 'resolved_external_id') return 'unique_external_mapping';
  if (status === 'resolved') return 'exact_unique_match';
  if (status === 'ambiguous') return 'ambiguous_not_used';
  if (status === 'api_error') return 'not_resolved_api_error';
  if (status === 'query_limit') return 'not_resolved_query_limit';
  if (status === 'empty') return 'not_applicable';
  return 'not_resolved';
}

function enrichChemicalMapping(mapping, identifierType = 'unknown') {
  if (!mapping || typeof mapping !== 'object') return mapping;
  return {
    ...mapping,
    sourceDatabase: inferMetaboliteSource(mapping.original, identifierType),
    canonicalDatabase: 'ChEBI',
    resolutionConfidence: chemicalResolutionConfidence(mapping),
    mappingUsedForIntegration: Boolean(mapping.resolved),
    ambiguityForced: false
  };
}

function countBy(items, keyFn) {
  const counts = {};
  for (const item of items) {
    const key = cleanText(keyFn(item)) || 'unknown';
    counts[key] = (counts[key] || 0) + 1;
  }
  return counts;
}

function summariseChemicalResolution(resolution, identifierType = 'unknown') {
  if (!resolution || typeof resolution !== 'object') return resolution;
  const mappings = Array.isArray(resolution.mappings)
    ? resolution.mappings.map((mapping) => enrichChemicalMapping(mapping, identifierType))
    : [];
  const total = mappings.length;
  const resolvedCount = mappings.filter((item) => item?.resolved).length;
  const ambiguousCount = mappings.filter((item) => item?.status === 'ambiguous').length;
  const apiErrorCount = mappings.filter((item) => item?.status === 'api_error').length;
  const queryLimitCount = mappings.filter((item) => item?.status === 'query_limit').length;
  const unresolvedCount = total - resolvedCount;
  const strictlyUnresolvedCount = mappings.filter((item) => item?.status === 'unresolved').length;
  return {
    ...resolution,
    mappings,
    resolvedCount,
    unresolvedCount,
    ambiguousCount,
    apiErrorCount,
    queryLimitCount,
    strictlyUnresolvedCount,
    coverageFraction: total ? resolvedCount / total : null,
    sourceTypeDeclared: identifierType,
    sourceDatabaseCounts: countBy(mappings, (item) => item?.sourceDatabase),
    statusCounts: countBy(mappings, (item) => item?.status),
    methodCounts: countBy(mappings, (item) => item?.method),
    canonicalDatabase: 'ChEBI',
    mappingPolicy: CHEMICAL_MAPPING_POLICY
  };
}

function attachChemicalResolutionAudit(result, args) {
  const raw = result?.identifierResolution?.metabolomics;
  if (!raw) return result;
  const identifierType = args?.identifierTypes?.metabolomics || 'unknown';
  const audited = summariseChemicalResolution(raw, identifierType);
  result.identifierResolution.metabolomics = audited;
  result.reproducibility = {
    ...(result.reproducibility || {}),
    chemicalIdentifierResolution: {
      metabolomics: {
        canonicalDatabase: audited.canonicalDatabase,
        sourceTypeDeclared: audited.sourceTypeDeclared,
        totalIdentifiers: audited.mappings?.length || 0,
        resolvedCount: audited.resolvedCount,
        unresolvedCount: audited.unresolvedCount,
        ambiguousCount: audited.ambiguousCount,
        apiErrorCount: audited.apiErrorCount,
        queryLimitCount: audited.queryLimitCount,
        coverageFraction: audited.coverageFraction,
        sourceDatabaseCounts: audited.sourceDatabaseCounts,
        statusCounts: audited.statusCounts,
        methodCounts: audited.methodCounts,
        mappingPolicy: audited.mappingPolicy
      }
    }
  };
  return result;
}

// Explicit wrappers intentionally override the star re-export above so callers
// receive the same mapping audit fields whether they resolve one identifier,
// a batch, or run the complete multi-omics pipeline.
export async function resolveMetaboliteIdentifier(identifier, options = {}) {
  const result = await impl.resolveMetaboliteIdentifier(identifier, options);
  return enrichChemicalMapping(result, options.identifierType || 'unknown');
}

export async function resolveMetaboliteIdentifiers(identifiers, options = {}) {
  const result = await impl.resolveMetaboliteIdentifiers(identifiers, options);
  return summariseChemicalResolution(result, options.identifierType || 'unknown');
}

function reactomeAnalysisTokens(result) {
  const tokens = [];
  const visit = (value) => {
    if (!value || typeof value !== 'object') return;
    if (typeof value.token === 'string' && value.token.trim()) tokens.push(value.token.trim());
    for (const child of Object.values(value)) visit(child);
  };
  visit(result?.reactome);
  return [...new Set(tokens)];
}

function stableCanonical(value) {
  if (Array.isArray(value)) return value.map(stableCanonical);
  if (value && typeof value === 'object') {
    return Object.fromEntries(
      Object.keys(value).sort().map((key) => [key, stableCanonical(value[key])])
    );
  }
  if (typeof value === 'number' && !Number.isFinite(value)) return null;
  return value;
}

function bytesToHex(buffer) {
  return [...new Uint8Array(buffer)].map((byte) => byte.toString(16).padStart(2, '0')).join('');
}

async function sha256Bytes(bytes) {
  if (!globalThis.crypto?.subtle) return null;
  const digest = await globalThis.crypto.subtle.digest('SHA-256', bytes);
  return bytesToHex(digest);
}

async function sha256Text(text) {
  return sha256Bytes(new TextEncoder().encode(String(text ?? '')));
}

async function sha256File(file) {
  if (!file || typeof file.arrayBuffer !== 'function') return null;
  return sha256Bytes(await file.arrayBuffer());
}

async function attachCryptographicInputProvenance(result, args) {
  const files = args?.files || {};
  const inputEntries = {};
  const keys = ['metadata', 'transcriptomics', 'proteomics', 'metabolomics'];
  let available = Boolean(globalThis.crypto?.subtle);

  for (const key of keys) {
    const file = files[key];
    if (!file) continue;
    let sha256 = null;
    try {
      sha256 = await sha256File(file);
    } catch {
      available = false;
    }
    inputEntries[key] = {
      name: file.name || key,
      size: Number.isFinite(file.size) ? file.size : null,
      type: file.type || null,
      lastModified: Number.isFinite(file.lastModified) ? file.lastModified : null,
      sha256
    };
  }

  let canonicalMetadataSha256 = null;
  try {
    const canonicalMetadata = JSON.stringify(stableCanonical(args?.metadataRows || []));
    canonicalMetadataSha256 = await sha256Text(canonicalMetadata);
  } catch {
    available = false;
  }

  const record = {
    status: available ? 'recorded' : 'unavailable',
    algorithm: 'SHA-256',
    files: inputEntries,
    canonicalMetadata: {
      rows: Array.isArray(args?.metadataRows) ? args.metadataRows.length : 0,
      sha256: canonicalMetadataSha256
    },
    note: available
      ? 'Cryptographic hashes are calculated locally in the browser from the exact uploaded bytes. Preserve them with immutable source files for manuscript or regulated-workflow archives.'
      : 'Web Crypto SHA-256 was unavailable; the older deterministic application fingerprint may still be present but is not a cryptographic integrity checksum.'
  };

  result.reproducibility = {
    ...(result.reproducibility || {}),
    cryptographicInputs: record
  };
  result.inputIntegrity = record;
  return result;
}

async function attachExternalDatabaseProvenance(result) {
  if (typeof window === 'undefined' || !result?.reactome || typeof globalThis.fetch !== 'function') return result;
  const capturedAt = new Date().toISOString();
  let reactome;
  try {
    const response = await globalThis.fetch(REACTOME_VERSION_URL, {
      method: 'GET',
      headers: { Accept: 'text/plain' },
      cache: 'no-store'
    });
    if (!response.ok) throw new Error(`HTTP ${response.status}`);
    const release = (await response.text()).trim();
    reactome = {
      database: 'Reactome',
      status: release ? 'recorded' : 'unavailable',
      release: release || null,
      versionEndpoint: REACTOME_VERSION_URL,
      analysisService: REACTOME_ANALYSIS_SERVICE,
      analysisTokens: reactomeAnalysisTokens(result),
      capturedAt
    };
  } catch (error) {
    reactome = {
      database: 'Reactome',
      status: 'unavailable',
      release: null,
      versionEndpoint: REACTOME_VERSION_URL,
      analysisService: REACTOME_ANALYSIS_SERVICE,
      analysisTokens: reactomeAnalysisTokens(result),
      capturedAt,
      error: error instanceof Error ? error.message : 'Reactome version request failed'
    };
  }

  result.externalDatabaseProvenance = {
    ...(result.externalDatabaseProvenance || {}),
    reactome,
    note: 'External database content can change independently of PMx Explain. Preserve database release identifiers alongside the analysis result whenever the service exposes them.'
  };

  // Keep the database release next to the pathway result as well as in the
  // generic reproducibility block. This ensures JSON/report exporters that
  // serialise either branch retain the exact external knowledge-base version.
  result.reactome.provenance = reactome;
  result.reproducibility = {
    ...(result.reproducibility || {}),
    externalDatabases: {
      ...(result.reproducibility?.externalDatabases || {}),
      reactome
    }
  };
  return result;
}

function publishFinalAnalysis(result) {
  if (typeof window === 'undefined' || typeof window.dispatchEvent !== 'function') return;
  window.__PMX_MULTIOMICS_ANALYSIS__ = result;
  try {
    // CustomEvent dispatch is synchronous: methodology guards listening to this
    // event modify the same object that is then returned and exported.
    window.dispatchEvent(new CustomEvent('pmx-multiomics-analysis', { detail: result }));
  } catch {
    // Exported results remain valid even on a legacy browser without CustomEvent.
  }
}

export async function runDeterministicAnalysis(args) {
  const preparedArgs = withBrowserMethodology(args);
  const result = await impl.runDeterministicAnalysis(preparedArgs);
  attachChemicalResolutionAudit(result, preparedArgs);
  await attachCryptographicInputProvenance(result, preparedArgs);
  await attachExternalDatabaseProvenance(result);
  publishFinalAnalysis(result);
  return result;
}
