// @ts-nocheck
// Deterministic, display-only data preparation for multi-omics figures.
// These helpers never perform inference. They reuse the uploaded matrices and
// the feature ranking already produced by the analysis engine.

import { parseDelimited } from './deterministic.js';
import { buildFocusedMetabolicNetwork } from './metabolic-network.js';

const LAYERS = ['transcriptomics', 'proteomics', 'metabolomics'];

function text(value) {
  return String(value ?? '').trim();
}

function key(value) {
  return text(value)
    .toLowerCase()
    .replace(/[α]/g, 'alpha')
    .replace(/[β]/g, 'beta')
    .replace(/[γ]/g, 'gamma')
    .replace(/[^a-z0-9]+/g, '');
}

function finite(value) {
  if (value === '' || value == null) return null;
  const parsed = Number(String(value).replace(',', '.'));
  return Number.isFinite(parsed) ? parsed : null;
}

function mean(values) {
  const x = values.filter(Number.isFinite);
  return x.length ? x.reduce((a, b) => a + b, 0) / x.length : NaN;
}

function variance(values) {
  const x = values.filter(Number.isFinite);
  if (x.length < 2) return NaN;
  const m = mean(x);
  return x.reduce((sum, value) => sum + (value - m) ** 2, 0) / (x.length - 1);
}

function median(values) {
  const x = values.filter(Number.isFinite).slice().sort((a, b) => a - b);
  if (!x.length) return NaN;
  const i = Math.floor(x.length / 2);
  return x.length % 2 ? x[i] : (x[i - 1] + x[i]) / 2;
}

function natural(values) {
  return [...new Set(values.filter(Boolean))].sort((a, b) =>
    String(a).localeCompare(String(b), undefined, { numeric: true, sensitivity: 'base' })
  );
}

function canonicalOmic(value) {
  const k = key(value);
  if (['transcriptomics', 'transcriptome', 'transcriptomique', 'rna', 'rnaseq', 'mrna', 'arn', 'geneexpression'].includes(k)) return 'transcriptomics';
  if (['proteomics', 'proteome', 'proteomique', 'protein', 'proteins', 'proteine', 'proteines', 'lfq'].includes(k)) return 'proteomics';
  if (['metabolomics', 'metabolome', 'metabolomique', 'metabolite', 'metabolites', 'met'].includes(k)) return 'metabolomics';
  return '';
}

function canonicalSampleType(value) {
  const k = key(value);
  if (['blank', 'blanc', 'solventblank'].includes(k)) return 'blank';
  if (['qc', 'pooledqc', 'qualitycontrol', 'pool'].includes(k)) return 'pooled_qc';
  return 'biological';
}

function mapped(row, mapping, name) {
  const column = mapping?.[name];
  return column ? text(row?.[column]) : '';
}

function canonicalMetadata(metadataRows, mapping, layer) {
  return metadataRows
    .filter((row) => canonicalOmic(mapped(row, mapping, 'omic')) === layer)
    .map((row) => ({
      subjectId: mapped(row, mapping, 'subject_id'),
      sampleId: mapped(row, mapping, 'sample_id'),
      assayId: mapped(row, mapping, 'assay_id'),
      condition: mapped(row, mapping, 'condition'),
      timepoint: mapped(row, mapping, 'timepoint'),
      batch: mapped(row, mapping, 'batch'),
      sampleType: canonicalSampleType(mapped(row, mapping, 'sample_type'))
    }))
    .filter((row) => row.assayId && row.sampleId && row.sampleType === 'biological');
}

function matrixFromParsed(parsed, expectedAssays) {
  if (!parsed?.headers?.length || parsed.headers.length < 2) return null;
  const first = parsed.headers[0];
  const columns = parsed.headers.slice(1);
  const rowIds = parsed.rows.map((row) => text(row[first])).filter(Boolean);
  const expected = new Set(expectedAssays);
  const columnHits = columns.filter((id) => expected.has(id)).length;
  const rowHits = rowIds.filter((id) => expected.has(id)).length;

  if (rowHits > columnHits) {
    const assays = parsed.rows.map((row) => text(row[first])).filter((id) => expected.has(id));
    const features = columns;
    const values = new Map(features.map((feature) => [feature, new Map()]));
    for (const row of parsed.rows) {
      const assay = text(row[first]);
      if (!expected.has(assay)) continue;
      for (const feature of features) values.get(feature).set(assay, finite(row[feature]));
    }
    return { features, assays, values, transposed: true };
  }

  const assays = columns.filter((id) => expected.has(id));
  const features = [];
  const values = new Map();
  for (const row of parsed.rows) {
    const feature = text(row[first]);
    if (!feature) continue;
    features.push(feature);
    values.set(feature, new Map(assays.map((assay) => [assay, finite(row[assay])])));
  }
  return { features, assays, values, transposed: false };
}

function minimumPositive(matrix) {
  let min = Infinity;
  for (const values of matrix.values.values()) {
    for (const value of values.values()) {
      if (Number.isFinite(value) && value > 0 && value < min) min = value;
    }
  }
  return Number.isFinite(min) ? min : 1;
}

function transformMatrix(matrix, layer, valueType) {
  const assays = matrix.assays;
  const processed = new Map();
  const minPositive = minimumPositive(matrix);
  const pseudo = Math.max(minPositive / 2, 1e-12);
  const allValues = [...matrix.values.values()].flatMap((values) => [...values.values()]);
  const hasNegative = allValues.some((value) => Number.isFinite(value) && value < 0);
  const isCounts = layer === 'transcriptomics' && valueType === 'raw_counts';
  const isSpectral = layer === 'proteomics' && valueType === 'spectral_count';
  const explicitlyLog = ['log_expression', 'log_intensity', 'log_abundance'].includes(valueType);
  const asSupplied = ['normalized', 'unknown'].includes(valueType) || (hasNegative && !isCounts && !isSpectral && !explicitlyLog);
  const logPositive = ['tpm', 'lfq_intensity', 'peak_area', 'concentration'].includes(valueType);
  const medianCenter = (layer === 'proteomics' && valueType === 'lfq_intensity') || (layer === 'metabolomics' && valueType === 'peak_area');
  let scale = 'as_supplied';
  const steps = [];

  const totals = new Map();
  if (isCounts || isSpectral) {
    scale = 'log2';
    for (const assay of assays) {
      let total = 0;
      for (const values of matrix.values.values()) {
        const value = values.get(assay);
        if (Number.isFinite(value) && value > 0) total += value;
      }
      totals.set(assay, total || 1);
    }
    steps.push(isCounts ? 'CPM + log2(CPM + 0.5)' : 'library-size normalization + log2');
  } else if (explicitlyLog) {
    scale = 'log2';
    steps.push('declared log scale');
  } else if (asSupplied) {
    scale = hasNegative ? 'transformed_unknown' : 'as_supplied';
    steps.push('values kept as supplied');
  } else if (logPositive) {
    scale = 'log2';
    steps.push('log2(value + half-minimum-positive pseudocount)');
  }

  for (const [feature, source] of matrix.values.entries()) {
    const next = new Map();
    for (const assay of assays) {
      const raw = source.get(assay);
      if (!Number.isFinite(raw)) {
        next.set(assay, null);
      } else if (isCounts || isSpectral) {
        next.set(assay, Math.log2(raw / totals.get(assay) * 1e6 + 0.5));
      } else if (explicitlyLog || asSupplied || !logPositive) {
        next.set(assay, raw);
      } else {
        next.set(assay, Math.log2(Math.max(0, raw) + pseudo));
      }
    }
    processed.set(feature, next);
  }

  if (medianCenter) {
    const assayMedians = new Map();
    for (const assay of assays) {
      assayMedians.set(assay, median([...processed.values()].map((values) => values.get(assay))));
    }
    const target = median([...assayMedians.values()]);
    for (const values of processed.values()) {
      for (const assay of assays) {
        const value = values.get(assay);
        if (Number.isFinite(value)) values.set(assay, value - assayMedians.get(assay) + target);
      }
    }
    steps.push('sample-wise median centering');
  }

  return { ...matrix, values: processed, scale, steps };
}

function aggregateReplicates(matrix, metadata) {
  const bySample = new Map();
  for (const row of metadata) {
    if (!bySample.has(row.sampleId)) bySample.set(row.sampleId, []);
    bySample.get(row.sampleId).push(row.assayId);
  }
  const sampleIds = [...bySample.keys()].sort((a, b) => {
    const ma = metadata.find((row) => row.sampleId === a);
    const mb = metadata.find((row) => row.sampleId === b);
    return [ma?.condition, ma?.timepoint, ma?.subjectId, a].join('|')
      .localeCompare([mb?.condition, mb?.timepoint, mb?.subjectId, b].join('|'), undefined, { numeric: true });
  });
  const sampleMeta = Object.fromEntries(sampleIds.map((sampleId) => [
    sampleId,
    metadata.find((row) => row.sampleId === sampleId)
  ]));
  const values = new Map();
  for (const feature of matrix.features) {
    const source = matrix.values.get(feature);
    const row = new Map();
    for (const sampleId of sampleIds) {
      const observed = bySample.get(sampleId).map((assay) => source.get(assay)).filter(Number.isFinite);
      row.set(sampleId, observed.length ? mean(observed) : null);
    }
    values.set(feature, row);
  }
  return { ...matrix, values, sampleIds, sampleMeta };
}

function topVariableFeatures(aggregated, limit) {
  return aggregated.features
    .map((feature) => {
      const values = aggregated.sampleIds.map((sample) => aggregated.values.get(feature)?.get(sample)).filter(Number.isFinite);
      return { feature, variance: variance(values) };
    })
    .filter((item) => Number.isFinite(item.variance) && item.variance > 0)
    .sort((a, b) => b.variance - a.variance || a.feature.localeCompare(b.feature))
    .slice(0, limit)
    .map((item) => item.feature);
}

function rankRows(layerResult, max = 30) {
  const rows = Array.isArray(layerResult?.rows) ? layerResult.rows : [];
  return rows.slice().sort((a, b) => {
    const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
    const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
    if (aq !== bq) return aq - bq;
    return Math.abs(Number(b.effect) || 0) - Math.abs(Number(a.effect) || 0);
  }).slice(0, max);
}

function buildHeatmap(aggregated, layerResult, maxFeatures = 30, maxSamples = 48) {
  const selected = Array.isArray(layerResult?.selected)
    ? layerResult.selected.map((row) => row.feature).filter((feature) => aggregated.features.includes(feature))
    : [];
  const features = [...new Set([...selected, ...topVariableFeatures(aggregated, maxFeatures)])].slice(0, maxFeatures);
  const samples = aggregated.sampleIds.slice(0, maxSamples);
  const stats = new Map((layerResult?.rows || []).map((row) => [row.feature, row]));
  const rows = [];

  for (const feature of features) {
    const source = aggregated.values.get(feature);
    const values = samples.map((sample) => source?.get(sample));
    const observed = values.filter(Number.isFinite);
    const center = mean(observed);
    const sd = Math.sqrt(variance(observed));
    if (!(sd > 0)) continue;
    const z = values.map((value) => Number.isFinite(value) ? (value - center) / sd : null);
    const stat = stats.get(feature);
    rows.push({
      feature,
      effect: Number.isFinite(stat?.effect) ? stat.effect : null,
      effectScale: stat?.effectScale || layerResult?.effectScale || aggregated.scale,
      qValue: Number.isFinite(stat?.qValue) ? stat.qValue : null,
      z,
      processedValues: values,
      center,
      sd
    });
  }

  return {
    scale: 'row_z_score',
    sourceScale: aggregated.scale,
    selectionRule: layerResult?.selectionRule || 'top variable features for descriptive visualization',
    samples: samples.map((sampleId) => ({
      sampleId,
      subjectId: aggregated.sampleMeta[sampleId]?.subjectId || '',
      condition: aggregated.sampleMeta[sampleId]?.condition || '',
      timepoint: aggregated.sampleMeta[sampleId]?.timepoint || '',
      batch: aggregated.sampleMeta[sampleId]?.batch || ''
    })),
    rows,
    samplesCapped: aggregated.sampleIds.length > samples.length,
    featuresCapped: features.length >= maxFeatures,
    method: 'Descriptive heatmap of display-transformed biological-sample values; technical replicates are averaged and rows z-scored. Display preprocessing is independent of the fitted model and may differ from its QC/batch adjustments. Statistical significance comes only from the feature models, not from the heatmap.',
    preprocessing: aggregated.steps
  };
}

// ChEBI names checked against EMBL-EBI ChEBI (stable identifiers). The mapping
// serves presentation only, except for the two explicitly curated carbon-map
// identities below; never guess an identity from an unknown numeric ID.
const KNOWN_METABOLITE_NAMES = Object.freeze({
  'chebi:16828': { fr: 'L-tryptophane', en: 'L-tryptophan' },
  'chebi:16946': { fr: 'L-kynurénine', en: 'L-kynurenine' },
  'chebi:24996': { fr: 'Lactate', en: 'Lactate' },
  'chebi:30031': { fr: 'Succinate', en: 'Succinate' }
});

/** Give known database IDs a readable name, preserving the original identifier. */
export function displayMetaboliteName(feature, language = 'fr') {
  const known = KNOWN_METABOLITE_NAMES[String(feature || '').trim().toLowerCase()];
  return known?.[language === 'en' ? 'en' : 'fr'] || String(feature || '');
}

function labelFor(feature, layer) {
  const known = layer === 'metabolomics' ? KNOWN_METABOLITE_NAMES[String(feature || '').trim().toLowerCase()] : null;
  return known ? { labelFr: known.fr, labelEn: known.en } : { labelFr: String(feature || ''), labelEn: String(feature || '') };
}

function effectEntries(layerResult, max = 28, layer = '') {
  return rankRows(layerResult, Math.max(max * 3, max))
    .filter((row) => Number.isFinite(row.effect))
    .map((row) => ({
      feature: row.feature,
      ...labelFor(row.feature, layer),
      effect: row.effect,
      qValue: Number.isFinite(row.qValue) ? row.qValue : null,
      effectScale: row.effectScale || layerResult?.effectScale || null
    }))
    .filter((row) => ['log2','as_supplied','transformed_unknown'].includes(row.effectScale))
    .slice(0, max);
}

function meanEffect(entries) {
  const x = entries.map((row) => row.effect).filter(Number.isFinite);
  return x.length ? mean(x) : null;
}

const CENTRAL_METABOLITES = [
  { id: 'glucose', label: 'Glucose', x: 90, y: 70, pathway: 'Glycolysis', aliases: ['glucose', 'dglucose'] },
  { id: 'g6p', label: 'G6P', x: 220, y: 70, pathway: 'Glycolysis / PPP', aliases: ['g6p', 'glucose6phosphate', 'glucose6p'] },
  { id: 'f6p', label: 'F6P', x: 350, y: 70, pathway: 'Glycolysis', aliases: ['f6p', 'fructose6phosphate', 'fructose6p'] },
  { id: 'fbp', label: 'FBP', x: 480, y: 70, pathway: 'Glycolysis', aliases: ['fbp', 'fructosebisphosphate', 'fructose16bisphosphate'] },
  { id: 'g3p', label: 'G3P', x: 610, y: 70, pathway: 'Glycolysis', aliases: ['g3p', 'glyceraldehyde3phosphate'] },
  { id: '3pg', label: '3-PG', x: 650, y: 125, pathway: 'Glycolysis', aliases: ['3pg', '3phosphoglycerate', '3phosphoglycericacid'] },
  { id: 'acetylcoa', label: 'AcCoA', x: 650, y: 250, pathway: 'TCA', aliases: ['acetylcoa', 'acetylcoenzymea'] },
  { id: 'pep', label: 'PEP', x: 740, y: 70, pathway: 'Glycolysis', aliases: ['pep', 'phosphoenolpyruvate'] },
  { id: 'pyruvate', label: 'PYR', x: 740, y: 180, pathway: 'Glycolysis', aliases: ['pyr', 'pyruvate', 'pyruvicacid'] },
  { id: 'lactate', label: 'LAC', x: 850, y: 180, pathway: 'Glycolysis', aliases: ['lac', 'lactate', 'llactate', 'lacticacid', 'CHEBI:24996'] },
  { id: 'alanine', label: 'ALA', x: 740, y: 290, pathway: 'Amino acids', aliases: ['ala', 'alanine', 'lalanine'] },
  { id: 'citrate', label: 'CIT', x: 560, y: 250, pathway: 'TCA', aliases: ['cit', 'citrate', 'citricacid'] },
  { id: 'akg', label: 'AKG', x: 430, y: 340, pathway: 'TCA', aliases: ['akg', 'alphaketoglutarate', '2oxoglutarate', 'ketoglutarate'] },
  { id: 'succinate', label: 'SUC', x: 500, y: 470, pathway: 'TCA', aliases: ['suc', 'succinate', 'succinicacid', 'CHEBI:30031'] },
  { id: 'fumarate', label: 'FUM', x: 650, y: 500, pathway: 'TCA', aliases: ['fum', 'fumarate', 'fumaricacid'] },
  { id: 'malate', label: 'MAL', x: 760, y: 430, pathway: 'TCA', aliases: ['mal', 'malate', 'lmalate', 'malicacid'] },
  { id: 'oaa', label: 'OAA', x: 710, y: 310, pathway: 'TCA', aliases: ['oaa', 'oxaloacetate', 'oxaloaceticacid'] },
  { id: 'glutamate', label: 'GLU', x: 290, y: 340, pathway: 'Amino acids', aliases: ['glu', 'glutamate', 'lglutamate', 'glutamicacid'] },
  { id: 'glutamine', label: 'GLN', x: 160, y: 340, pathway: 'Amino acids', aliases: ['gln', 'glutamine', 'lglutamine'] },
  { id: 'aspartate', label: 'ASP', x: 710, y: 550, pathway: 'MAS / amino acids', aliases: ['asp', 'aspartate', 'laspartate', 'asparticacid'] },
  { id: 'asparagine', label: 'ASN', x: 850, y: 550, pathway: 'Amino acids', aliases: ['asn', 'asparagine', 'lasparagine'] },
  { id: '6pg', label: '6PG', x: 220, y: 185, pathway: 'PPP', aliases: ['6pg', '6phosphogluconate', '6phosphogluconicacid'] },
  { id: 'r5p', label: 'R5P', x: 220, y: 285, pathway: 'PPP', aliases: ['r5p', 'ribose5phosphate'] },
  { id: 'prpp', label: 'PRPP', x: 90, y: 285, pathway: 'PPP', aliases: ['prpp', 'phosphoribosylpyrophosphate', '5phosphoribosyl1pyrophosphate'] },
  { id: 'serine', label: 'SER', x: 610, y: 180, pathway: 'Ser/Gly', aliases: ['ser', 'serine', 'lserine'] },
  { id: 'glycine', label: 'GLY', x: 500, y: 180, pathway: 'Ser/Gly', aliases: ['gly', 'glycine'] },
  { id: 'gly3p', label: 'Gly3P', x: 480, y: 150, pathway: 'Glycerol', aliases: ['gly3p', 'glycerol3phosphate'] },
  { id: '2hg', label: '2HG', x: 300, y: 455, pathway: 'TCA-related', aliases: ['2hg', '2hydroxyglutarate', '2hydroxyglutaricacid'] }
];

const CENTRAL_ENZYMES = [
  { id: 'hk', label: 'HK', x: 155, y: 48, aliases: ['hk1', 'hk2', 'hk3'] },
  { id: 'gpi', label: 'GPI', x: 285, y: 48, aliases: ['gpi'] },
  { id: 'pfk', label: 'PFK', x: 415, y: 48, aliases: ['pfkm', 'pfkl', 'pfkp'] },
  { id: 'aldo', label: 'ALDO', x: 545, y: 48, aliases: ['aldoa', 'aldob', 'aldoc'] },
  { id: 'gapdh', label: 'GAPDH', x: 675, y: 48, aliases: ['gapdh'] },
  { id: 'eno', label: 'ENO', x: 713, y: 102, aliases: ['eno1', 'eno2', 'eno3'] },
  { id: 'pkm', label: 'PKM', x: 760, y: 125, aliases: ['pkm', 'pklr'] },
  { id: 'ldha', label: 'LDHA', x: 815, y: 137, aliases: ['ldha'] },
  { id: 'ldhb', label: 'LDHB', x: 857, y: 210, aliases: ['ldhb'] },
  { id: 'pdh', label: 'PDH', x: 740, y: 235, aliases: ['pdha1', 'pdha2', 'pdhb'] },
  { id: 'cs', label: 'CS', x: 555, y: 295, aliases: ['cs'] },
  { id: 'idh', label: 'IDH', x: 490, y: 295, aliases: ['idh1', 'idh2', 'idh3a', 'idh3b', 'idh3g'] },
  { id: 'ogdh', label: 'OGDH', x: 455, y: 405, aliases: ['ogdh', 'ogdhl'] },
  { id: 'sdh', label: 'SDH', x: 575, y: 490, aliases: ['sdha', 'sdhb', 'sdhc', 'sdhd'] },
  { id: 'fh', label: 'FH', x: 705, y: 475, aliases: ['fh'] },
  { id: 'mdh', label: 'MDH', x: 760, y: 365, aliases: ['mdh1', 'mdh2'] },
  { id: 'got', label: 'GOT', x: 645, y: 400, aliases: ['got1', 'got2'] },
  { id: 'gls', label: 'GLS', x: 225, y: 315, aliases: ['gls', 'gls2'] },
  { id: 'glud', label: 'GLUD', x: 355, y: 315, aliases: ['glud1', 'glud2'] },
  { id: 'glul', label: 'GLUL', x: 225, y: 365, aliases: ['glul'] },
  { id: 'g6pd', label: 'G6PD', x: 190, y: 130, aliases: ['g6pd'] },
  { id: 'suclg', label: 'SUCLG', x: 475, y: 430, aliases: ['suclg1', 'suclg2', 'sucla2'] },
  { id: 'pgd', label: 'PGD', x: 190, y: 235, aliases: ['pgd'] },
  { id: 'phgdh', label: 'PHGDH', x: 585, y: 130, aliases: ['phgdh'] },
  { id: 'shmt', label: 'SHMT', x: 555, y: 160, aliases: ['shmt1', 'shmt2'] }
];

const CENTRAL_EDGES = [
  ['glucose','g6p'],['g6p','f6p'],['f6p','fbp'],['fbp','g3p'],['g3p','3pg'],['3pg','pep'],['pep','pyruvate'],
  ['pyruvate','lactate'],['pyruvate','alanine'],['pyruvate','acetylcoa'],['acetylcoa','citrate'],['citrate','akg'],['akg','succinate'],
  ['succinate','fumarate'],['fumarate','malate'],['malate','oaa'],['oaa','citrate'],['glutamine','glutamate'],
  ['glutamate','akg'],['oaa','aspartate'],['aspartate','asparagine'],['g6p','6pg'],['6pg','r5p'],['r5p','prpp'],
  ['3pg','serine'],['serine','glycine'],['akg','2hg']
];

function mappingAliases(result, layer) {
  const map = new Map();
  const rows = result?.identifierResolution?.[layer]?.mappings || [];
  for (const item of rows) {
    const aliases = [item?.original, item?.resolved, item?.label, item?.query].filter(Boolean);
    if (item?.original) map.set(item.original, aliases);
  }
  if (layer === 'metabolomics') for (const entry of result?.userMetaboliteAnnotations || []) {
    map.set(entry.feature, [...new Set([...(map.get(entry.feature)||[]),entry.feature,entry.chebi])]);
  }
  return map;
}

function featureKeys(feature, resolutionAliases) {
  return [...new Set([feature, ...(resolutionAliases.get(feature) || [])].map(key).filter(Boolean))];
}

function chooseEffect(rows, aliases, resolutionAliases, allowGeneFamily = false) {
  const wanted = new Set(aliases.map(key));
  const candidates = rows
    .filter((row) => featureKeys(row.feature, resolutionAliases).some((candidate) => wanted.has(candidate)))
    .filter((row) => Number.isFinite(row.effect));
  const byQ = (a, b) => {
    const aq = Number.isFinite(a.qValue) ? a.qValue : 1;
    const bq = Number.isFinite(b.qValue) ? b.qValue : 1;
    if (aq !== bq) return aq - bq;
    return Math.abs(Number(b.effect) || 0) - Math.abs(Number(a.effect) || 0);
  };
  candidates.sort(byQ);
  if (!candidates.length) return null;
  const valid = candidates.filter((row) =>
    ['log2','as_supplied','transformed_unknown'].includes(row.effectScale));
  if (!valid.length) {
    const best = candidates[0];
    return { feature: best.feature, effect: null, qValue: best.qValue ?? null,
      effectScale: best.effectScale || null, status: 'unavailable_effect_scale' };
  }

  // Gene families (e.g. HK1/HK2 or SDHA/SDHB) must not hide opposite
  // transcriptional directions behind the minimum q-value.
  const positives = valid.some((row) => row.effect > 0);
  const negatives = valid.some((row) => row.effect < 0);
  const scales = new Set(valid.map((row) => row.effectScale));
  const features = [...new Set(valid.map((row) => row.feature))];
  // Several MS peaks/adducts pointing to one ChEBI must never become one
  // arbitrarily selected biological measurement. Keep unresolved until explicit
  // quantification rules or replicate aggregation are specified.
  if (!allowGeneFamily && features.length > 1)
    return { feature: features.join(', '), effect:null, qValue:null,
      effectScale:scales.size===1 ? valid[0].effectScale : 'mixed',
      status:'multiple_features', matchedFeatures:features };
  if (features.length > 1 && (positives && negatives || scales.size > 1)) {
    return { feature: features.join(', '), effect: null, qValue: null,
      effectScale: scales.size === 1 ? valid[0].effectScale : 'mixed', status: 'discordant_isoforms', matchedFeatures: features };
  }
  const best = valid.slice().sort(byQ)[0];
  return {
    feature: best.feature,
    effect: best.effect,
    qValue: Number.isFinite(best.qValue) ? best.qValue : null,
    effectScale: best.effectScale,
    matchedFeatures: features,
    status: 'measured'
  };
}

function buildCentralCarbon(result) {
  const metaboliteRows = result?.layers?.metabolomics?.rows || [];
  const transcriptRows = result?.layers?.transcriptomics?.rows || [];
  const metaboliteAliases = mappingAliases(result, 'metabolomics');
  const transcriptAliases = mappingAliases(result, 'transcriptomics');

  const metabolites = CENTRAL_METABOLITES.map((node) => ({
    ...node,
    measurement: chooseEffect(metaboliteRows, node.aliases, metaboliteAliases)
  }));
  const enzymes = CENTRAL_ENZYMES.map((node) => ({
    ...node,
    measurement: chooseEffect(transcriptRows, node.aliases, transcriptAliases, true)
  }));
  const shownMetaboliteIds = new Set(metabolites.flatMap((node) =>
    node.measurement?.status === 'measured' ? node.measurement.matchedFeatures || [node.measurement.feature] : []));
  const shownTranscriptIds = new Set(enzymes.flatMap((node) =>
    node.measurement?.status === 'measured' ? node.measurement.matchedFeatures || [node.measurement.feature] : []));
  const offMap = (rows, shown, layer) => rows
    .filter((row) => Number.isFinite(row.effect) && row.effectScale === 'log2' && !shown.has(row.feature))
    .map((row) => ({ feature: row.feature, ...labelFor(row.feature, layer) }));
  return {
    metabolites,
    enzymes,
    edges: CENTRAL_EDGES,
    measuredMetabolites: metabolites.filter((node) => node.measurement?.status === 'measured').length,
    measuredTranscripts: enzymes.filter((node) => node.measurement?.status === 'measured').length,
    offMapMetabolites: offMap(metaboliteRows, shownMetaboliteIds, 'metabolomics'),
    offMapTranscripts: offMap(transcriptRows, shownTranscriptIds, 'transcriptomics'),
    method: 'Fixed schematic carbon-pathway map; colors use only model differential effects on a declared log2 scale. When multiple mapped features have opposite effects, the node is gray. Otherwise the lowest-q feature is a representative, not a pooled enzyme activity. Unmeasured or non-log2 effects also remain gray.',
    scope: 'Visualization aid only. The pathway layout is not used to compute enrichment, statistics or causality.'
  };
}


/**
 * Figure 4C-inspired, CURATED pathway selections. No over-representation or
 * pathway activity score is inferred by membership. Data are model log2 effects,
 * never 13C isotopologues or metabolic flux. Unresolved IDs are not guessed.
 */
const METABOLOGRAM_PATHWAYS = [
  { id: 'glycolysis', labelFr: 'Glycolyse', labelEn: 'Glycolysis',
    metabolites: ['glucose','g6p','f6p','fbp','g3p','3pg','pep','pyruvate','lactate'],
    enzymes: ['hk','gpi','pfk','aldo','gapdh','eno','pkm','ldha','ldhb'] },
  { id: 'tca', labelFr: 'Cycle de Krebs', labelEn: 'TCA cycle',
    metabolites: ['acetylcoa','citrate','akg','succinate','fumarate','malate','oaa'],
    enzymes: ['pdh','cs','idh','ogdh','suclg','sdh','fh','mdh'] },
  { id: 'aminoacids', labelFr: 'Acides aminés', labelEn: 'Amino acids',
    metabolites: ['alanine','serine','glycine','glutamate','glutamine','aspartate','asparagine'],
    enzymes: ['phgdh','shmt','gls','glud','glul','got'] },
  { id: 'ppp', labelFr: 'Voie des pentoses phosphates', labelEn: 'Pentose phosphate pathway',
    metabolites: ['g6p','6pg','r5p','prpp'],
    enzymes: ['g6pd','pgd'] }
];

function buildMetabologramPathways(central) {
  function entries(nodes, ids) {
    const used = new Set();
    return ids.map((id) => nodes.find((node) => node.id === id))
      .filter(Boolean)
      .map((node) => node.measurement)
      .filter((item) => item?.status === 'measured' && Number.isFinite(item.effect))
      // Multiple complexes may map to the same transcript: count each observed
      // feature at most once per pathway and omics layer.
      .filter((item) => !used.has(item.feature) && used.add(item.feature))
      .map((item) => ({ feature: item.feature, ...labelFor(item.feature, nodes === central.metabolites ? 'metabolomics' : 'transcriptomics'),
        effect: item.effect, qValue: item.qValue, effectScale: item.effectScale }));
  }
  return METABOLOGRAM_PATHWAYS.map((pathway) => {
    const metabolomics = entries(central.metabolites, pathway.metabolites);
    const transcriptomics = entries(central.enzymes, pathway.enzymes);
    return {
      id: pathway.id,
      labelFr: pathway.labelFr,
      labelEn: pathway.labelEn,
      metabolomics,
      transcriptomics,
      meanMetabolomicLog2Fc: meanEffect(metabolomics),
      meanTranscriptomicLog2Fc: meanEffect(transcriptomics),
      coverageMetabolites: metabolomics.length,
      coverageTranscripts: transcriptomics.length,
      method: 'Descriptive pathway selection from a fixed, reviewed list of central metabolites and enzyme genes. No pathway significance, isotope tracing, metabolic flux or causal interpretation is inferred.'
    };
  });
}

export async function buildMultiomicsVisualizationData({
  files,
  metadataRows,
  columnMapping,
  dataTypes,
  analysisResult,
  metaboliteAnnotations = /** @type {{feature:string,chebi:string,status:string,provenance:string}[]} */ ([]),
  maxFeatures = 30,
  maxSamples = 48
}) {
  const heatmaps = {};

  for (const layer of LAYERS) {
    if (!files?.[layer]) continue;
    const metadata = canonicalMetadata(metadataRows || [], columnMapping || {}, layer);
    if (!metadata.length) continue;
    const expectedAssays = metadata.map((row) => row.assayId);
    const parsed = parseDelimited(await files[layer].text());
    const matrix = matrixFromParsed(parsed, expectedAssays);
    if (!matrix?.features?.length || !matrix?.assays?.length) continue;
    const transformed = transformMatrix(matrix, layer, dataTypes?.[layer] || 'unknown');
    const aggregated = aggregateReplicates(transformed, metadata);
    heatmaps[layer] = buildHeatmap(
      aggregated,
      analysisResult?.layers?.[layer],
      maxFeatures,
      maxSamples
    );
  }

  // Pure visualization overlay. Statistical model objects remain untouched.
  const annotatedResult = metaboliteAnnotations.length
    ? { ...analysisResult, userMetaboliteAnnotations: metaboliteAnnotations } : analysisResult;
  const transcriptEntries = effectEntries(analysisResult?.layers?.transcriptomics, 28, 'transcriptomics');
  const metaboliteEntries = effectEntries(analysisResult?.layers?.metabolomics, 28, 'metabolomics');
  const centralCarbon = buildCentralCarbon(annotatedResult);

  return {
    heatmaps,
    metabologram: {
      transcriptomics: transcriptEntries,
      metabolomics: metaboliteEntries,
      meanTranscriptomicLog2Fc: meanEffect(transcriptEntries),
      meanMetabolomicLog2Fc: meanEffect(metaboliteEntries),
      method: 'Circular effect map ranked by fitted feature-model results. Effects retain their declared scales: log2 or native supplied units. Each layer uses its own color range; across-layer magnitude comparisons are invalid. Central semicircles are descriptive means, not inferential statistics.',
      requiresLog2Effect: false
    },
    metabologramPathways: buildMetabologramPathways(centralCarbon),
    centralCarbon,
    focusedMetabolicNetwork: buildFocusedMetabolicNetwork(annotatedResult, centralCarbon),
    methodologicalBoundary: [
      'Figures reuse the statistical results; they do not run additional hypothesis tests.',
      'Heatmap clustering is intentionally not used by default: deterministic sample order follows study annotations and feature order follows the analysis ranking.',
      'Heatmap row z-scores are descriptive and must not be interpreted as fold changes.',
      'Metabologram and pathway map retain effects on their declared scales; native supplied units are not converted to fold ratios. Color ranges are layer-specific and visual magnitudes must not be compared across different scales.',
      'The curated pathway panels are descriptive annotations, not enrichment or pathway activity tests. No metabolic flux is estimated from abundances or gene expression.',
      'The expanded network is an auditable schematic neighborhood, not a stoichiometrically validated reaction graph. Only explicitly curated ChEBI IDs are exact; name-only and unresolved matches are reported separately.',
      'An optional uploaded feature-to-ChEBI crosswalk is user-declared, not independently chemically verified; it affects figures only and must not be interpreted as MS/MS identification evidence.'
    ]
  };
}
