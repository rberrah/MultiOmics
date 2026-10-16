import fs from 'node:fs/promises';
import { File } from 'node:buffer';
import { parseDelimited, runDeterministicAnalysis, differentialCorrelationPair, layerOverlapSummary, resolveMetaboliteIdentifier, resolveMetaboliteIdentifiers, reactomeOverRepresentation } from '../src/lib/multiomics/deterministic.js';
import { buildMultiomicsVisualizationData } from '../src/lib/multiomics/visualization-data.js';

const root = process.argv[2] || 'tmp/public-benchmarks';

async function loadFile(path) {
  const text = await fs.readFile(path, 'utf8');
  return new File([text], path.split('/').pop(), { type: 'text/csv' });
}

async function loadDataset(dir, layers) {
  const metadata = await loadFile(`${dir}/metadata.csv`);
  const parsed = parseDelimited(await metadata.text());
  const files = { metadata, transcriptomics: null, proteomics: null, metabolomics: null };
  for (const layer of layers) files[layer] = await loadFile(`${dir}/${layer}.csv`);
  return {
    files,
    metadataRows: parsed.rows,
    mapping: {
      subject_id: 'subject_id',
      sample_id: 'sample_id',
      assay_id: 'assay_id',
      omic: 'omic',
      condition: 'condition',
      timepoint: 'timepoint',
      batch: 'batch',
      technical_replicate: 'technical_replicate',
      outcome: 'outcome'
    }
  };
}

function findRanks(rows, patterns) {
  const out = [];
  for (const pattern of patterns) {
    const regex = new RegExp(pattern, 'i');
    const hits = rows
      .map((row, index) => ({ ...row, rank: index + 1 }))
      .filter((row) => regex.test(row.feature));
    out.push({ pattern, hits: hits.slice(0, 10) });
  }
  return out;
}

function top(rows, n = 15) {
  return rows.slice(0, n).map((row) => ({
    feature: row.feature,
    foldRatio: row.foldRatio == null ? null : Number(row.foldRatio.toFixed(4)),
    effect: Number(row.effect.toFixed(4)),
    effectScale: row.effectScale,
    p: row.pValue,
    q: row.qValue
  }));
}

function pathwayMatches(consensus, patterns) {
  const regexes = patterns.map((p) => new RegExp(p, 'i'));
  return (consensus || []).filter((p) => regexes.some((r) => r.test(p.name))).slice(0, 25);
}

function requireTruth(condition, label, details = '') {
  if (!condition) {
    const suffix = details ? ` — ${details}` : '';
    throw new Error(`PUBLIC BENCHMARK FAILED: ${label}${suffix}`);
  }
  console.log(`TRUTH PASS: ${label}${details ? ` — ${details}` : ''}`);
}

const report = { generatedAt: new Date().toISOString(), benchmarks: {} };

// --- Nutrimouse ---
{
  const ds = await loadDataset(`${root}/nutrimouse`, ['transcriptomics', 'metabolomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'mouse',
      objective: 'groups',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'animal',
      groupCount: '2',
      sampleOverlap: 'same_specimen'
    },
    dataTypes: {
      transcriptomics: 'normalized',
      proteomics: 'log_intensity',
      metabolomics: 'concentration'
    },
    useReactome: true,
    resolveIdentifiers: false
  });

  const knownPatterns = [
    'Cyp3a11', 'Car$', 'Nr1i3', 'Lpin', 'Scd', 'Fads2', 'Eci', 'Ppara'
  ];

  const pathwayPatterns = [
    'fatty acid', 'lipid', 'peroxisom', 'xenobiotic', 'cytochrome p450', 'metabolism'
  ];

  report.benchmarks.nutrimouse = {
    truth: {
      design: 'WT vs PPARalpha-null mice under five diets',
      expectedBiology: 'PPARalpha-dependent lipid and xenobiotic metabolism'
    },
    metadata: result.metadataSummary,
    transcriptomics: {
      contrast: result.layers.transcriptomics.contrast,
      groupSizes: result.layers.transcriptomics.groupSizes,
      top: top(result.layers.transcriptomics.rows, 20),
      knownMarkerRanks: findRanks(result.layers.transcriptomics.rows, knownPatterns)
    },
    metabolomics: {
      contrast: result.layers.metabolomics.contrast,
      groupSizes: result.layers.metabolomics.groupSizes,
      top: top(result.layers.metabolomics.rows, 20)
    },
    reactomeError: result.reactomeError,
    pathwayMatches: pathwayMatches(result.reactome?.consensus, pathwayPatterns),
    topPathways: (result.reactome?.consensus || []).slice(0, 20).map((p) => ({
      id: p.id, name: p.name, fdr: p.fdr, supportingLayers: p.supportingLayers
    }))
  };

  const cyp3a11 = result.layers.transcriptomics.rows
    .map((row, index) => ({ ...row, rank: index + 1 }))
    .find((row) => /^CYP3A11$/i.test(row.feature));
  requireTruth(
    cyp3a11 && cyp3a11.rank <= 10 && cyp3a11.qValue <= 0.05,
    'Nutrimouse recovers CYP3A11 genotype signal',
    cyp3a11 ? `rank=${cyp3a11.rank}, q=${cyp3a11.qValue}` : 'CYP3A11 absent'
  );

  const linoleic = result.layers.metabolomics.rows.find((row) => /^C18\.2n\.6$/i.test(row.feature));
  requireTruth(
    linoleic && linoleic.foldRatio < 1 && linoleic.qValue <= 0.10,
    'Nutrimouse recovers elevated C18:2n-6 in PPARalpha-null mice',
    linoleic ? `WT/PPAR=${linoleic.foldRatio.toFixed(3)}, q=${linoleic.qValue}` : 'C18.2n.6 absent'
  );

  const expectedPathway = (result.reactome?.consensus || []).find((pathway) =>
    /(PPARalpha|PPARA|fatty acid metabolism|metabolism of lipids|triglyceride metabolism)/i.test(pathway.name)
  );
  requireTruth(
    Boolean(expectedPathway),
    'Nutrimouse Reactome mapping recovers PPAR/lipid biology',
    expectedPathway?.name || 'no expected pathway'
  );
}

// --- TCGA Her2 vs LumA ---
{
  const ds = await loadDataset(`${root}/tcga-her2-luma`, ['transcriptomics', 'proteomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'groups',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      groupCount: '2',
      sampleOverlap: 'same_specimen'
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'concentration'
    },
    useReactome: true,
    resolveIdentifiers: false
  });

  const her2Patterns = ['ERBB2', 'HER2', 'GRB7', 'STARD3', 'PGAP3'];
  const mrnaRanks = findRanks(result.layers.transcriptomics.rows, her2Patterns);
  const proteinRanks = findRanks(result.layers.proteomics.rows, her2Patterns);

  const allHer2Hits = [...mrnaRanks, ...proteinRanks].flatMap((x) => x.hits);
  const expectedDirectionHits = allHer2Hits.filter((x) => x.effect < 0); // result is LumA - Her2

  report.benchmarks.tcgaHer2LumA = {
    truth: {
      design: 'TCGA breast Her2-enriched vs Luminal A',
      expectedBiology: 'HER2/ERBB2-associated signal higher in Her2-enriched tumours'
    },
    metadata: result.metadataSummary,
    transcriptomics: {
      contrast: result.layers.transcriptomics.contrast,
      groupSizes: result.layers.transcriptomics.groupSizes,
      top: top(result.layers.transcriptomics.rows, 20),
      markerRanks: mrnaRanks
    },
    proteomics: {
      contrast: result.layers.proteomics.contrast,
      groupSizes: result.layers.proteomics.groupSizes,
      top: top(result.layers.proteomics.rows, 20),
      markerRanks: proteinRanks
    },
    markerHits: allHer2Hits.length,
    expectedDirectionHits: expectedDirectionHits.length,
    reactomeError: result.reactomeError,
    topPathways: (result.reactome?.consensus || []).slice(0, 20).map((p) => ({
      id: p.id, name: p.name, fdr: p.fdr, supportingLayers: p.supportingLayers
    }))
  };

  const her2 = result.layers.proteomics.rows
    .map((row, index) => ({ ...row, rank: index + 1 }))
    .find((row) => /^HER2$/i.test(row.feature));
  const her2Phospho = result.layers.proteomics.rows
    .map((row, index) => ({ ...row, rank: index + 1 }))
    .find((row) => /^HER2_pY1248$/i.test(row.feature));

  requireTruth(
    her2 && her2.rank <= 10 && her2.effect < 0 && her2.qValue <= 0.01,
    'TCGA recovers HER2 protein enrichment in Her2 subtype',
    her2 ? `rank=${her2.rank}, LumA/Her2=${her2.foldRatio.toFixed(3)}, q=${her2.qValue}` : 'HER2 absent'
  );
  requireTruth(
    her2Phospho && her2Phospho.rank <= 10 && her2Phospho.effect < 0 && her2Phospho.qValue <= 0.01,
    'TCGA recovers activated HER2-pY1248 enrichment in Her2 subtype',
    her2Phospho ? `rank=${her2Phospho.rank}, LumA/Her2=${her2Phospho.foldRatio.toFixed(3)}, q=${her2Phospho.qValue}` : 'HER2_pY1248 absent'
  );

  const erbb2Pathway = (result.reactome?.consensus || []).find((pathway) => /ERBB2|HER2/i.test(pathway.name));
  requireTruth(
    Boolean(erbb2Pathway),
    'TCGA Reactome result contains ERBB2/HER2 signalling',
    erbb2Pathway?.name || 'no ERBB2/HER2 pathway'
  );
}

// --- TCGA outcome prediction: same public biology, now testing supervised/outcome branches ---
{
  const ds = await loadDataset(`${root}/tcga-her2-luma`, ['transcriptomics', 'proteomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'outcome',
      outcomeType: 'binary',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      groupCount: '1',
      sampleOverlap: 'same_specimen',
      batchKnown: 'no',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'concentration'
    },
    useReactome: false,
    resolveIdentifiers: false
  });

  report.benchmarks.tcgaHer2LumAOutcome = {
    truth: {
      target: 'Her2 vs LumA subtype',
      rationale: 'HER2-enriched and Luminal A tumours should be separable from transcriptomic/proteomic profiles'
    },
    predictiveOutcome: result.predictiveOutcome,
    supervisedIntegration: result.supervisedIntegration
  };

  requireTruth(
    result.predictiveOutcome?.status === 'ok',
    'TCGA binary outcome branch completes nested cross-validation',
    result.predictiveOutcome?.reason || `folds=${result.predictiveOutcome?.outerFolds || 0}`
  );
  requireTruth(
    Number.isFinite(result.predictiveOutcome?.metrics?.auc) && result.predictiveOutcome.metrics.auc >= 0.80,
    'TCGA multi-omics prediction separates Her2 from LumA out of sample',
    `AUC=${Number(result.predictiveOutcome?.metrics?.auc || NaN).toFixed(3)}`
  );
  requireTruth(
    result.supervisedIntegration?.subjects >= 20 &&
      Number.isFinite(result.supervisedIntegration?.scoreTargetCorrelation) &&
      Math.abs(result.supervisedIntegration.scoreTargetCorrelation) >= 0.50,
    'TCGA supervised multiblock component is linked to subtype',
    `n=${result.supervisedIntegration?.subjects || 0}, r=${Number(result.supervisedIntegration?.scoreTargetCorrelation || NaN).toFixed(3)}`
  );
}

await fs.writeFile(`${root}/benchmark-report.json`, JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));


async function testPublishedPair(filename, reference, comparison, truth) {
  const parsed = parseDelimited(await fs.readFile(`${root}/${filename}`, 'utf8'));
  const clean = parsed.rows
    .map((row) => ({
      condition: row.condition,
      x: Number(row.gene),
      y: Number(row.metabolite)
    }))
    .filter((row) => Number.isFinite(row.x) && Number.isFinite(row.y));
  const ref = clean.filter((row) => row.condition === reference);
  const cmp = clean.filter((row) => row.condition === comparison);
  const stat = differentialCorrelationPair(
    ref.map((row) => row.x), ref.map((row) => row.y),
    cmp.map((row) => row.x), cmp.map((row) => row.y),
    'spearman'
  );
  truth(stat, ref.length, cmp.length);
  return { filename, reference, comparison, nReference: ref.length, nComparison: cmp.length, ...stat };
}

report.benchmarks.nci60IntLIM = {};
report.benchmarks.nci60IntLIM.FAM174B_malic = await testPublishedPair(
  'nci60_FAM174B_malic.csv',
  'Leukemia',
  'BPO',
  (stat, nRef, nCmp) => {
    requireTruth(nRef >= 5 && nCmp >= 10, 'NCI-60 IntLIM pair has expected group sizes', `Leukemia=${nRef}, BPO=${nCmp}`);
    requireTruth(
      stat.rReference < 0 && stat.rComparison > 0 && stat.deltaR > 0.5,
      'NCI-60 reproduces FAM174B–malic-acid correlation reversal',
      `rLeukemia=${stat.rReference.toFixed(3)}, rBPO=${stat.rComparison.toFixed(3)}, Δr=${stat.deltaR.toFixed(3)}`
    );
  }
);

report.benchmarks.nci60IntLIM.DNER_imidazole = await testPublishedPair(
  'nci60_DNER_imidazole.csv',
  'BPO',
  'Leukemia',
  (stat) => {
    requireTruth(
      stat.rReference < 0 && stat.rComparison > 0 && stat.deltaR > 0.5,
      'NCI-60 reproduces DNER–imidazolelactate correlation reversal',
      `rBPO=${stat.rReference.toFixed(3)}, rLeukemia=${stat.rComparison.toFixed(3)}, Δr=${stat.deltaR.toFixed(3)}`
    );
  }
);

report.benchmarks.brcaIntLIM = {};
report.benchmarks.brcaIntLIM.GPT2_2HG = await testPublishedPair(
  'brca_GPT2_2HG.csv',
  'NORMAL',
  'TUMOR',
  (stat, nRef, nCmp) => {
    requireTruth(nRef >= 40 && nCmp >= 50, 'BRCA IntLIM pair has substantial matched groups', `Normal=${nRef}, Tumor=${nCmp}`);
    requireTruth(
      stat.rComparison > stat.rReference && stat.deltaR > 0.25,
      'BRCA reproduces stronger GPT2–2-hydroxyglutarate coupling in tumour',
      `rNormal=${stat.rReference.toFixed(3)}, rTumor=${stat.rComparison.toFixed(3)}, Δr=${stat.deltaR.toFixed(3)}`
    );
  }
);

await fs.writeFile(`${root}/benchmark-report.json`, JSON.stringify(report, null, 2));
console.log('PUBLIC MULTI-OMICS TRUTH BENCHMARKS: PASS');


{
  const parsed = parseDelimited(await fs.readFile(`${root}/missrows_nci60_metadata.csv`, 'utf8'));
  const canonical = parsed.rows.map((row) => ({
    subjectId: row.subject_id,
    sampleId: row.sample_id,
    assayId: row.assay_id,
    omic: row.omic,
    condition: row.condition,
    timepoint: '',
    batch: '',
    technicalReplicate: '',
    outcome: ''
  }));
  const overlap = layerOverlapSummary(canonical, ['transcriptomics','proteomics']);
  report.benchmarks.missRowsNCI60 = overlap;
  requireTruth(overlap.layerSubjects.transcriptomics === 48, 'missRows preserves 48 transcriptomic individuals', String(overlap.layerSubjects.transcriptomics));
  requireTruth(overlap.layerSubjects.proteomics === 52, 'missRows preserves 52 proteomic individuals', String(overlap.layerSubjects.proteomics));
  requireTruth(overlap.allMatched === 40, 'missRows matches only the 40 genuinely shared individuals', `matched=${overlap.allMatched}`);
}

{
  const ds = await loadDataset(`${root}/missrows-nci60-partial`, ['transcriptomics','proteomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'explore',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'cell_line',
      groupCount: '1',
      sampleOverlap: 'partial',
      batchKnown: 'no',
      partialOmicsExpected: 'yes',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'normalized',
      proteomics: 'normalized',
      metabolomics: 'concentration'
    },
    useReactome: false,
    resolveIdentifiers: false
  });

  report.benchmarks.missRowsNCI60PartialIntegration = {
    metadata: result.metadataSummary,
    exploration: result.exploration
  };

  requireTruth(
    result.exploration?.subjects === 60 && result.exploration?.completeSubjects === 40,
    'missRows partial-block exploration keeps all 60 subjects without pretending all blocks are complete',
    `union=${result.exploration?.subjects}, complete=${result.exploration?.completeSubjects}`
  );
  requireTruth(
    result.exploration?.imputedCells > 0 &&
      /partial blocks/i.test(result.exploration?.missingDataPolicy || ''),
    'missRows uses explicit latent-only partial-block filling',
    `imputed latent cells=${result.exploration?.imputedCells || 0}`
  );
  requireTruth(
    result.metadataSummary?.overlap?.allMatched === 40 &&
      result.metadataSummary?.overlap?.totalSubjectsAnyLayer === 60,
    'missRows overlap audit remains exact inside the engine',
    `matched=${result.metadataSummary?.overlap?.allMatched}, any-layer=${result.metadataSummary?.overlap?.totalSubjectsAnyLayer}`
  );
}

// External identifier-resolution probe: useful operational evidence but an EBI outage must not invalidate biological truth benchmarks.
report.apiProbes = {};
for (const term of ['malic acid', 'C18.2n.6', '2-hydroxyglutarate']) {
  const mapping = await resolveMetaboliteIdentifier(term);
  report.apiProbes[term] = mapping;
  if (mapping.resolved) {
    console.log(`API PASS: ChEBI resolved ${term} -> ${mapping.resolved} (${mapping.label || mapping.query})`);
  } else {
    console.warn(`API WARN: ChEBI did not resolve ${term}: ${mapping.status}`);
  }
}



function signConcordance(engineRows, referenceRows) {
  const byFeature = new Map(engineRows.map((row) => [row.feature, row]));
  const pairs = referenceRows
    .filter((ref) => Number.isFinite(ref.officialLogFC) && Number.isFinite(ref.officialQ) && ref.officialQ <= 0.10)
    .map((ref) => ({ ref, row: byFeature.get(ref.feature) }))
    .filter((x) => x.row && Number.isFinite(x.row.effect) && x.row.effect !== 0 && x.ref.officialLogFC !== 0);
  const concordant = pairs.filter(({ref,row}) => Math.sign(row.effect) === -Math.sign(ref.officialLogFC)).length;
  return {
    n: pairs.length,
    concordant,
    fraction: pairs.length ? concordant/pairs.length : null
  };
}

{
  const ds = await loadDataset(`${root}/aging-hfcd-old-vs-young-cd`, ['transcriptomics','proteomics','metabolomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'mouse',
      objective: 'groups',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'animal',
      groupCount: '2',
      sampleOverlap: 'partial'
    },
    dataTypes: {
      transcriptomics: 'normalized',
      proteomics: 'normalized',
      metabolomics: 'normalized'
    },
    useReactome: true,
    resolveIdentifiers: false
  });

  const refParsed = parseDelimited(await fs.readFile(`${root}/aging-hfcd-old-vs-young-cd/reference.csv`, 'utf8'));
  const reference = refParsed.rows.map((row) => ({
    layer: row.layer,
    feature: row.feature,
    officialLogFC: Number(row.official_old_vs_young_logFC),
    officialQ: Number(row.official_q)
  }));

  const concordance = {};
  for (const layer of ['transcriptomics','proteomics','metabolomics']) {
    concordance[layer] = signConcordance(
      result.layers[layer].rows,
      reference.filter((row) => row.layer === layer)
    );
    requireTruth(
      concordance[layer].n >= 10 && concordance[layer].fraction >= 0.80,
      `AgingHFCD ${layer} agrees with reference Old-vs-Young effect directions`,
      `n=${concordance[layer].n}, concordance=${(100*concordance[layer].fraction).toFixed(1)}%`
    );
  }

  const ctsdRna = result.layers.transcriptomics.rows
    .map((row,index) => ({...row,rank:index+1}))
    .find((row) => /^Ctsd$/i.test(row.feature));
  const ctsdProtein = result.layers.proteomics.rows
    .map((row,index) => ({...row,rank:index+1}))
    .find((row) => /^Ctsd$/i.test(row.feature));

  // Engine contrast is Young_CD - Old_CD; published Ctsd aging association predicts a negative effect here.
  requireTruth(
    ctsdRna && ctsdRna.effect < 0 && ctsdRna.qValue <= 0.10,
    'AgingHFCD recovers higher Ctsd RNA in old control-diet mice',
    ctsdRna ? `rank=${ctsdRna.rank}, effect=${ctsdRna.effect.toFixed(3)}, q=${ctsdRna.qValue}` : 'Ctsd RNA absent'
  );
  requireTruth(
    ctsdProtein && ctsdProtein.effect < 0 && ctsdProtein.qValue <= 0.10,
    'AgingHFCD recovers higher CTSD protein in old control-diet mice',
    ctsdProtein ? `rank=${ctsdProtein.rank}, effect=${ctsdProtein.effect.toFixed(3)}, q=${ctsdProtein.qValue}` : 'CTSD protein absent'
  );

  requireTruth(
    result.metadataSummary.overlap.allMatched >= 20,
    'AgingHFCD has a non-trivial three-layer matched cohort',
    `all-three matched=${result.metadataSummary.overlap.allMatched}`
  );

  // Second, independent validation: run the actual Figure-4-inspired
  // visualization adapter on this PUBLIC 3-omics dataset rather than injecting
  // a synthetic result. A node must never be colored by a missing feature or
  // by a result on a non-log2 effect scale.
  const visuals = await buildMultiomicsVisualizationData({
    files: ds.files,
    metadataRows: ds.metadataRows,
    columnMapping: ds.mapping,
    dataTypes: { transcriptomics: 'normalized', proteomics: 'normalized', metabolomics: 'normalized' },
    analysisResult: result
  });
  const vizCoverage = {
    heatmaps: Object.fromEntries(Object.entries(visuals.heatmaps).map(([layer, view]) => [
      layer, { features: view.rows.length, biologicalSamples: view.samples.length }
    ])),
    metabologram: {
      geneSegments: visuals.metabologram.transcriptomics.length,
      metaboliteSegments: visuals.metabologram.metabolomics.length
    },
    focusedNetwork: {
      curatedMetabolites: visuals.focusedMetabolicNetwork.coverage.metabolites,
      curatedGenes: visuals.focusedMetabolicNetwork.coverage.genes,
      exactIds: visuals.focusedMetabolicNetwork.audit.filter((a) => a.status === 'verified_id').length,
      nameOnly: visuals.focusedMetabolicNetwork.audit.filter((a) => a.status === 'name_only').length,
      unknown: visuals.focusedMetabolicNetwork.audit.filter((a) => !a.id).length,
      displayedRegions: visuals.focusedMetabolicNetwork.regions.filter((r) => r.measuredCount > 0)
        .map((r) => ({ region: r.id, measuredCount: r.measuredCount }))
    }
  };
  const evidence = new Map(result.layers.transcriptomics.rows.map((row) => [row.feature, row]));
  requireTruth(
    visuals.metabologram.transcriptomics.every((item) =>
      evidence.has(item.feature) && evidence.get(item.feature).effect === item.effect),
    'AgingHFCD visualization exactly preserves the fitted transcript effects'
  );
  const metaboliteEvidence = new Map(result.layers.metabolomics.rows.map((row) => [row.feature, row]));
  requireTruth(
    visuals.metabologram.metabolomics.every((item) =>
      metaboliteEvidence.has(item.feature) && metaboliteEvidence.get(item.feature).effect === item.effect),
    'AgingHFCD visualization exactly preserves fitted metabolite effects'
  );
  requireTruth(
    Object.keys(visuals.heatmaps).length === 3
      && Object.values(visuals.heatmaps).every((view) => view.samples.length >= 10 && view.rows.length > 0),
    'AgingHFCD builds all three real sample-level heatmaps',
    JSON.stringify(vizCoverage.heatmaps)
  );
  requireTruth(
    visuals.focusedMetabolicNetwork.audit.every((item) => item.id ||
      !visuals.focusedMetabolicNetwork.nodes.some((node) =>
        node.measurement?.features?.includes(item.input))),
    'AgingHFCD does not color unmapped metabolite identifiers'
  );
  console.log('PUBLIC VISUALIZATION EVIDENCE AgingHFCD', JSON.stringify(vizCoverage));

  report.benchmarks.agingHFCD = {
    truth: 'Old_CD vs Young_CD; Ctsd is an aging-associated RNA/protein signal in the source publication',
    visualization: vizCoverage,
    metadata: result.metadataSummary,
    concordance,
    ctsdRna,
    ctsdProtein,
    crossOmics: {
      testedPairs: result.crossOmics.testedPairs,
      significantPairs: result.crossOmics.significantPairs,
      top: result.crossOmics.pairs.slice(0,20)
    },
    topPathways: (result.reactome?.consensus || []).slice(0,20).map((p) => ({
      id:p.id,name:p.name,fdr:p.fdr,supportingLayers:p.supportingLayers
    })),
    reactomeError: result.reactomeError
  };
}

await fs.writeFile(`${root}/benchmark-report.json`, JSON.stringify(report, null, 2));
console.log('PUBLIC MULTI-OMICS EXTENDED BENCHMARKS: PASS');


// --- STATegra public Ikaros time-course: precomputed Ikaros/control trajectories ---
{
  const paintRoot = process.env.PAINTOMICS_ROOT || 'external/PaintOmics';
  const dataRoot = `${paintRoot}/PaintomicsServer/src/examplefiles/datasets/08-stategra-multiomics/data`;

  const readTab = async (name) => parseDelimited(await fs.readFile(`${dataRoot}/${name}`, 'utf8'));
  const protein = await readTab('proteomics_values.tab');
  const metabolite = await readTab('metabolomics_values.tab');
  const relevantProteins = (await fs.readFile(`${dataRoot}/proteomics_relevant.tab`, 'utf8')).split(/\r?\n/).map(x=>x.trim()).filter(Boolean);
  const relevantMetabolites = (await fs.readFile(`${dataRoot}/metabolomics_relevant.tab`, 'utf8')).split(/\r?\n/).map(x=>x.trim()).filter(Boolean);

  const pName = protein.headers[0];
  const mName = metabolite.headers[0];
  const p24 = protein.headers[protein.headers.length-1];
  const m24 = metabolite.headers[metabolite.headers.length-1];

  const proteinRow = (name) => protein.rows.find((row) => String(row[pName]).toLowerCase() === name.toLowerCase());
  const metaboliteRow = (name) => metabolite.rows.find((row) => String(row[mName]).toLowerCase() === name.toLowerCase());
  const value24 = (row, col) => row ? Number(row[col]) : NaN;

  const ldha = proteinRow('Ldha');
  const hk2 = proteinRow('Hk2');
  const slc7a5 = proteinRow('Slc7a5');
  const ikzf1 = proteinRow('Ikzf1');
  const igll1 = proteinRow('Igll1');
  const lactate = metaboliteRow('Lactic acid');
  const pyruvate = metaboliteRow('Pyruvic acid');
  const malate = metaboliteRow('Malic acid');

  requireTruth(ldha && value24(ldha,p24) < -1, 'STATegra recovers late LDHA suppression', `24h log2 ratio=${value24(ldha,p24).toFixed(3)}`);
  requireTruth(hk2 && value24(hk2,p24) < -0.5, 'STATegra recovers late HK2 suppression', `24h log2 ratio=${value24(hk2,p24).toFixed(3)}`);
  requireTruth(slc7a5 && value24(slc7a5,p24) < -1, 'STATegra recovers late SLC7A5 suppression', `24h log2 ratio=${value24(slc7a5,p24).toFixed(3)}`);
  requireTruth(igll1 && value24(igll1,p24) < -2, 'STATegra recovers Igll1 suppression during differentiation', `24h log2 ratio=${value24(igll1,p24).toFixed(3)}`);
  requireTruth(ikzf1 && value24(ikzf1,p24) > 2, 'STATegra preserves strong IKZF1 protein signal', `24h log2 ratio=${value24(ikzf1,p24).toFixed(3)}`);
  requireTruth(lactate && value24(lactate,m24) < 0, 'STATegra recovers lower lactate at 24h', `24h log2 ratio=${value24(lactate,m24).toFixed(3)}`);
  requireTruth(pyruvate && value24(pyruvate,m24) < 0, 'STATegra recovers lower pyruvate at 24h', `24h log2 ratio=${value24(pyruvate,m24).toFixed(3)}`);
  requireTruth(malate && value24(malate,m24) < -0.5, 'STATegra recovers lower malate at 24h', `24h log2 ratio=${value24(malate,m24).toFixed(3)}`);

  // Deterministic knowledge-layer test anchored to the central-carbon biology
  // explicitly validated in the STATegra publication.
  const centralCarbonMetabolites = [
    'Lactic acid','Pyruvic acid','Malic acid','Citric acid',
    'Succinic acid','Fumaric acid','Alpha-ketoglutaric acid','Glucose'
  ];
  const metaboliteResolution = await resolveMetaboliteIdentifiers(centralCarbonMetabolites, { maxQueries: 12 });
  const resolvedMetabolites = metaboliteResolution.mappings.map(x=>x.resolved).filter(Boolean);
  const pathwayInput = ['LDHA','HK2','SLC7A5', ...resolvedMetabolites];
  let pathways = [];
  let reactomeError = null;
  try {
    const reactome = await reactomeOverRepresentation(pathwayInput, { projectToHuman: true, pageSize: 100 });
    pathways = reactome.pathways || [];
  } catch (error) {
    reactomeError = error instanceof Error ? error.message : String(error);
  }

  const glycolysisRelated = pathways.find((p) => /glycol|glucose|pyruvate|citric acid|tricarbox|TCA/i.test(p.name));
  if (!reactomeError) {
    requireTruth(
      Boolean(glycolysisRelated),
      'STATegra targeted knowledge integration recovers central-carbon metabolism',
      glycolysisRelated?.name || `top pathways=${pathways.slice(0,8).map(p=>p.name).join(' | ')}`
    );
  } else {
    console.warn(`API WARN: STATegra Reactome probe unavailable: ${reactomeError}`);
  }

  report.benchmarks.stategra = {
    truth: 'Ikaros-induced pre-B differentiation over 0–24 h; known late repression of glycolysis/TCA-related genes and metabolites',
    sourceFormat: 'published Ikaros-minus-control log2 ratios, not individual-level observations',
    protein24h: {
      Ldha: value24(ldha,p24),
      Hk2: value24(hk2,p24),
      Slc7a5: value24(slc7a5,p24),
      Ikzf1: value24(ikzf1,p24),
      Igll1: value24(igll1,p24)
    },
    metabolite24h: {
      lactate: value24(lactate,m24),
      pyruvate: value24(pyruvate,m24),
      malate: value24(malate,m24)
    },
    metaboliteResolution: {
      attempted: metaboliteResolution.mappings.length,
      resolved: metaboliteResolution.resolvedCount,
      unresolved: metaboliteResolution.unresolvedCount
    },
    centralCarbonPathway: glycolysisRelated ? {id:glycolysisRelated.id,name:glycolysisRelated.name,fdr:glycolysisRelated.fdr} : null,
    reactomeError
  };
}

// --- LRRK2 G2019S RNA + proteome at day 35 ---
{
  const ds = await loadDataset(`${root}/lrrk2-d35`, ['transcriptomics','proteomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'groups',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'cell_model',
      groupCount: '2',
      sampleOverlap: 'partial'
    },
    dataTypes: {
      transcriptomics: 'normalized',
      proteomics: 'normalized',
      metabolomics: 'normalized'
    },
    useReactome: true,
    resolveIdentifiers: false
  });

  const refParsed = parseDelimited(await fs.readFile(`${root}/lrrk2-d35/reference.csv`, 'utf8'));
  const reference = refParsed.rows.map((row) => ({
    layer: row.layer,
    feature: row.feature,
    officialLogFC: Number(row.official_logFC_G2019S_vs_Control),
    officialQ: Number(row.official_q)
  }));

  const concordance = {};
  for (const layer of ['transcriptomics','proteomics']) {
    const byFeature = new Map(result.layers[layer].rows.map((row) => [row.feature,row]));
    const pairs = reference
      .filter((row) => row.layer === layer && Number.isFinite(row.officialLogFC) && Number.isFinite(row.officialQ) && row.officialQ <= 0.10)
      .map((ref) => ({ref,row:byFeature.get(ref.feature)}))
      .filter((x) => x.row && Number.isFinite(x.row.effect) && x.row.effect !== 0 && x.ref.officialLogFC !== 0);
    const concordant = pairs.filter(({ref,row}) => Math.sign(row.effect) === Math.sign(ref.officialLogFC)).length;
    concordance[layer] = {
      n:pairs.length,
      concordant,
      fraction:pairs.length ? concordant/pairs.length : null
    };
    requireTruth(
      concordance[layer].n >= 20 && concordance[layer].fraction >= 0.80,
      `LRRK2 ${layer} agrees with published G2019S-vs-control effect directions`,
      `n=${concordance[layer].n}, concordance=${(100*concordance[layer].fraction).toFixed(1)}%`
    );
  }

  const rab29 = result.layers.transcriptomics.rows.find((row) => /^RAB29$/i.test(row.feature));
  requireTruth(
    rab29 && rab29.effect > 0 && rab29.qValue <= 0.10,
    'LRRK2 transcriptome recovers RAB29 up-regulation at day 35',
    rab29 ? `effect=${rab29.effect.toFixed(3)}, q=${rab29.qValue}` : 'RAB29 absent'
  );

  const rab25Protein = result.layers.proteomics.rows.find((row) => /^RAB25$/i.test(row.feature));
  requireTruth(
    rab25Protein && rab25Protein.effect > 0 && rab25Protein.qValue <= 0.10,
    'LRRK2 proteome recovers RAB25 up-regulation at day 35',
    rab25Protein ? `effect=${rab25Protein.effect.toFixed(3)}, q=${rab25Protein.qValue}` : 'RAB25 absent'
  );

  let targetedLrrkPathways = [];
  let targetedLrrkError = null;
  try {
    const targeted = await reactomeOverRepresentation(
      ['RAB29','RAB25','RAB10','RAB3A','RAB3B','CLTC','SYNJ1','SYNJ2','DNM1L','SH3GLB1','SH3GLB2'],
      { projectToHuman: true, pageSize: 100 }
    );
    targetedLrrkPathways = targeted.pathways || [];
  } catch (error) {
    targetedLrrkError = error instanceof Error ? error.message : String(error);
  }
  const endocyticPathway = targetedLrrkPathways.find((p) =>
    /endocyt|vesicle|RAB|clathrin|membrane trafficking|transport/i.test(p.name)
  );
  if (!targetedLrrkError) {
    requireTruth(
      Boolean(endocyticPathway),
      'LRRK2 targeted Reactome probe recovers RAB/endocytic trafficking biology',
      endocyticPathway?.name || `top pathways=${targetedLrrkPathways.slice(0,10).map(p=>p.name).join(' | ')}`
    );
  } else {
    console.warn(`API WARN: LRRK2 targeted Reactome probe unavailable: ${targetedLrrkError}`);
  }

  report.benchmarks.lrrk2 = {
    truth: 'LRRK2-G2019S dopaminergic-neuron model; published integrated RNA/protein analysis reports endocytic and RAB dysregulation',
    metadata: result.metadataSummary,
    concordance,
    rab29,
    rab25Protein,
    endocyticPathway: endocyticPathway ? {id:endocyticPathway.id,name:endocyticPathway.name,fdr:endocyticPathway.fdr} : null,
    targetedReactomeError: targetedLrrkError,
    reactomeError: result.reactomeError
  };
}

// --- TCGA breast 3-subtype multi-group inference ---
{
  const ds = await loadDataset(`${root}/tcga-three-subtypes`, ['transcriptomics','proteomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'groups',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      groupCount: '3',
      sampleOverlap: 'same_specimen'
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'normalized'
    },
    useReactome: true,
    resolveIdentifiers: false
  });

  requireTruth(
    result.layers.proteomics.mode === 'adjusted-multi-group-model',
    'TCGA three-subtype benchmark uses the explicit adjusted multi-group model',
    result.layers.proteomics.mode + ' · ' + result.layers.proteomics.inferenceMethod
  );

  const ranked = result.layers.proteomics.rows.map((row,index)=>({...row,rank:index+1}));
  const her2 = ranked.find((row)=>/^HER2$/i.test(row.feature));
  const er = ranked.find((row)=>/^ER-alpha$/i.test(row.feature));
  const pr = ranked.find((row)=>/^PR$/i.test(row.feature));

  requireTruth(
    her2 && her2.qValue <= 0.01 && her2.groupMeans?.Her2 > her2.groupMeans?.LumA && her2.groupMeans?.Her2 > her2.groupMeans?.Basal,
    'TCGA multi-group analysis identifies HER2 subtype-specific protein signal',
    her2 ? `rank=${her2.rank}, q=${her2.qValue}` : 'HER2 absent'
  );
  requireTruth(
    er && er.qValue <= 0.01 && er.groupMeans?.LumA > er.groupMeans?.Her2 && er.groupMeans?.LumA > er.groupMeans?.Basal,
    'TCGA multi-group analysis identifies LumA ER-alpha enrichment',
    er ? `rank=${er.rank}, q=${er.qValue}` : 'ER-alpha absent'
  );
  requireTruth(
    pr && pr.qValue <= 0.01 && pr.groupMeans?.LumA > pr.groupMeans?.Her2 && pr.groupMeans?.LumA > pr.groupMeans?.Basal,
    'TCGA multi-group analysis identifies LumA PR enrichment',
    pr ? `rank=${pr.rank}, q=${pr.qValue}` : 'PR absent'
  );

  report.benchmarks.tcgaThreeSubtypes = {
    truth: 'Basal / Her2 / LumA breast-cancer subtypes with established HER2 and hormone-receptor protein differences',
    metadata: result.metadataSummary,
    proteomics: {
      mode: result.layers.proteomics.mode,
      top: result.layers.proteomics.rows.slice(0,20),
      HER2: her2,
      ERalpha: er,
      PR: pr
    },
    reactomeError: result.reactomeError
  };
}

// --- PaintOmics planted multi-omics truth: independent RNA/protein layers should converge on the same planted module ---
{
  const paintRoot = process.env.PAINTOMICS_ROOT || 'external/PaintOmics';
  const datasetRoot = `${paintRoot}/PaintomicsServer/src/examplefiles/datasets/04-multiomics-integration`;
  const readList = async (path) => (await fs.readFile(path, 'utf8'))
    .split(/\r?\n/)
    .map((x) => x.trim())
    .filter((x) => x && !x.startsWith('#'));

  const geneRelevant = await readList(`${datasetRoot}/data/gene_expression_relevant.tab`);
  const proteinRelevant = await readList(`${datasetRoot}/data/proteomics_relevant.tab`);
  const metaboliteRelevant = await readList(`${datasetRoot}/data/metabolomics_relevant.tab`);
  const planted = new Set(await readList(`${datasetRoot}/expected/signal_features.txt`));
  const expectedPathways = await readList(`${datasetRoot}/expected/expected_pathways.txt`);

  const proteinSet = new Set(proteinRelevant);
  const sharedGeneProtein = geneRelevant.filter((id) => proteinSet.has(id));
  const sharedPlanted = sharedGeneProtein.filter((id) => planted.has(id));
  const plantedFraction = sharedGeneProtein.length ? sharedPlanted.length / sharedGeneProtein.length : 0;

  requireTruth(
    planted.size === 324,
    'PaintOmics planted multi-omics fixture exposes the recorded 324-feature ground truth',
    `planted=${planted.size}`
  );
  requireTruth(
    sharedGeneProtein.length >= 80 && plantedFraction >= 0.85,
    'PaintOmics RNA/protein relevant sets converge on the planted module',
    `shared=${sharedGeneProtein.length}, planted-shared=${sharedPlanted.length} (${(100*plantedFraction).toFixed(1)}%)`
  );
  requireTruth(
    metaboliteRelevant.length === 496 && expectedPathways.length === 8,
    'PaintOmics fixture preserves the recorded metabolite and pathway truth sets',
    `metabolites=${metaboliteRelevant.length}, pathways=${expectedPathways.length}`
  );

  report.benchmarks.paintOmicsPlanted = {
    truth: 'Simulated five-omics PaintOmics fixture with a shared planted signal; benchmark uses the three molecular layers supported by the current PMx engine.',
    plantedFeatures: planted.size,
    geneRelevant: geneRelevant.length,
    proteinRelevant: proteinRelevant.length,
    metaboliteRelevant: metaboliteRelevant.length,
    sharedGeneProtein: sharedGeneProtein.length,
    sharedGeneProteinPlanted: sharedPlanted.length,
    sharedGeneProteinPlantedFraction: plantedFraction,
    expectedPathways: expectedPathways.length
  };
}

// --- STATegra sample-level metabolomics: verify replicate design and direction against the public six-timepoint summary ---
{
  const paintRoot = process.env.PAINTOMICS_ROOT || 'external/PaintOmics';
  const replicateRoot = `${paintRoot}/PaintomicsServer/src/examplefiles/datasets/12-stategra-metabolomics-replicates/data`;
  const summaryRoot = `${paintRoot}/PaintomicsServer/src/examplefiles/datasets/08-stategra-multiomics/data`;

  const design = parseDelimited(await fs.readFile(`${replicateRoot}/experimental_design.tab`, 'utf8'));
  const replicateMatrix = parseDelimited(await fs.readFile(`${replicateRoot}/metabolomics_replicates.tab`, 'utf8'));
  const summaryMatrix = parseDelimited(await fs.readFile(`${summaryRoot}/metabolomics_values.tab`, 'utf8'));
  const relevantMetabolites = (await fs.readFile(`${replicateRoot}/metabolomics_relevant.tab`, 'utf8'))
    .split(/\r?\n/)
    .map((x) => x.trim())
    .filter((x) => x && !x.startsWith('#'));

  const sampleColumn = design.headers[0];
  const conditionColumn = design.headers[1];
  const conditionCounts = new Map();
  for (const row of design.rows) {
    const condition = row[conditionColumn];
    conditionCounts.set(condition, (conditionCounts.get(condition) || 0) + 1);
  }

  requireTruth(
    design.rows.length === 36 && conditionCounts.size === 12 && [...conditionCounts.values()].every((n) => n === 3),
    'STATegra replicate design preserves 12 condition-time cells with three biological replicates each',
    `samples=${design.rows.length}, cells=${conditionCounts.size}, replicates=${[...new Set(conditionCounts.values())].join('/')}`
  );
  requireTruth(
    replicateMatrix.rows.length === 58 && replicateMatrix.headers.length - 1 === 36 && relevantMetabolites.length === 34,
    'STATegra replicate matrix preserves the published 58-metabolite panel and 34-feature relevant set',
    `features=${replicateMatrix.rows.length}, sample-columns=${replicateMatrix.headers.length-1}, relevant=${relevantMetabolites.length}`
  );

  const repId = replicateMatrix.headers[0];
  const summaryId = summaryMatrix.headers[0];
  const summaryByName = new Map(summaryMatrix.rows.map((row) => [String(row[summaryId]).toLowerCase(), row]));
  const times = ['0','2','6','12','18','24'];
  const mean = (values) => values.reduce((sum, value) => sum + value, 0) / values.length;
  let comparable = 0;
  let concordant = 0;
  const scaleRatios = [];

  for (const row of replicateMatrix.rows) {
    const summaryRow = summaryByName.get(String(row[repId]).toLowerCase());
    if (!summaryRow) continue;
    for (const time of times) {
      const ctr = replicateMatrix.headers
        .filter((header) => header.startsWith(`Ctr_${time}H_`))
        .map((header) => Number(row[header]))
        .filter(Number.isFinite);
      const ik = replicateMatrix.headers
        .filter((header) => header.startsWith(`Ik_${time}H_`))
        .map((header) => Number(row[header]))
        .filter(Number.isFinite);
      const summaryColumn = `IKvsCtr_${time}h`;
      const published = Number(summaryRow[summaryColumn]);
      if (ctr.length !== 3 || ik.length !== 3 || !Number.isFinite(published)) continue;
      const derived = mean(ik) - mean(ctr);
      if (Math.abs(derived) < 1e-12 || Math.abs(published) < 1e-12) continue;
      comparable += 1;
      if (Math.sign(derived) === Math.sign(published)) concordant += 1;
      if (Math.abs(published) >= 0.05 && Math.sign(derived) === Math.sign(published)) {
        scaleRatios.push(Math.abs(derived / published));
      }
    }
  }

  scaleRatios.sort((a,b) => a-b);
  const medianScaleRatio = scaleRatios.length
    ? (scaleRatios.length % 2
        ? scaleRatios[(scaleRatios.length-1)/2]
        : (scaleRatios[scaleRatios.length/2-1] + scaleRatios[scaleRatios.length/2]) / 2)
    : null;
  const directionConcordance = comparable ? concordant / comparable : 0;

  requireTruth(
    comparable >= 300 && directionConcordance >= 0.95,
    'STATegra replicate-level Ikaros-minus-control trajectories agree in direction with the public multi-omics summary',
    `n=${comparable}, concordance=${(100*directionConcordance).toFixed(1)}%`
  );

  if (medianScaleRatio != null && Math.abs(medianScaleRatio - 1) > 0.20) {
    console.warn(
      `PUBLIC BENCHMARK NOTE: STATegra replicate and summary files agree in direction but not numeric scale (median |replicate contrast / summary contrast|=${medianScaleRatio.toFixed(3)}). The engine must therefore treat uploaded scale declarations explicitly rather than assuming these public files are numerically interchangeable.`
    );
  }

  report.benchmarks.stategraReplicates = {
    truth: 'Real STATegra metabolomics with 36 sample-level observations: control vs Ikaros, six time points, three biological replicates per cell.',
    sampleColumn,
    samples: design.rows.length,
    conditionTimeCells: conditionCounts.size,
    replicatesPerCell: [...new Set(conditionCounts.values())],
    metabolites: replicateMatrix.rows.length,
    relevantMetabolites: relevantMetabolites.length,
    comparableTrajectoryPoints: comparable,
    directionConcordance,
    medianAbsoluteScaleRatioVsSummary: medianScaleRatio,
    scaleNote: medianScaleRatio != null && Math.abs(medianScaleRatio - 1) > 0.20
      ? 'Direction is reproducible, but the two public processed representations are not on the same numeric scale; do not merge them without an explicit scale declaration.'
      : null
  };
}

// --- Public validation of the new unsupervised branch: TCGA breast without using subtype labels ---
{
  const ds = await loadDataset(root + '/tcga-three-subtypes', ['transcriptomics','proteomics']);
  const result = await runDeterministicAnalysis({
    ...ds,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'explore',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      groupCount: '1',
      sampleOverlap: 'same_specimen',
      batchKnown: 'no',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'normalized'
    },
    useReactome: false,
    resolveIdentifiers: false
  });

  const loadings = (result.exploration?.components || []).flatMap((component) =>
    (component.topLoadings || []).map((item) => ({ ...item, component: component.component }))
  );
  const breastMarkers = loadings.filter((item) => /HER2|ERBB2|ER-alpha|PR|GRB7/i.test(item.feature));
  requireTruth(
    result.exploration?.subjects >= 100 && result.exploration?.components?.length >= 1,
    'TCGA unsupervised multi-block branch computes shared latent structure without subtype labels',
    'subjects=' + result.exploration?.subjects + ', components=' + (result.exploration?.components?.length || 0)
  );
  requireTruth(
    breastMarkers.length >= 1,
    'TCGA unsupervised loadings contain established breast-subtype biology',
    breastMarkers.slice(0,5).map((item) => item.feature + '@PC' + item.component).join(', ') || 'no HER2/ER/PR/GRB7 marker among top loadings'
  );

  report.benchmarks.tcgaUnsupervised = {
    truth: 'Basal / Her2 / LumA samples analysed without giving subtype labels to the exploratory model.',
    subjects: result.exploration?.subjects,
    components: result.exploration?.components?.map((component) => ({
      component: component.component,
      explainedFraction: component.explainedFraction,
      topLoadings: component.topLoadings.slice(0,15)
    })),
    knownMarkerLoadings: breastMarkers
  };
}

// --- Public validation of the multiclass outcome branch: TCGA Basal / Her2 / LumA phenotype ---
{
  const ds = await loadDataset(root + '/tcga-three-subtypes', ['transcriptomics','proteomics']);
  const outcomeRows = ds.metadataRows.map((row) => ({ ...row, outcome: row.condition }));
  const result = await runDeterministicAnalysis({
    ...ds,
    metadataRows: outcomeRows,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'outcome',
      outcomeType: 'multiclass',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      groupCount: '1',
      sampleOverlap: 'same_specimen',
      batchKnown: 'no',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'normalized'
    },
    useReactome: false,
    resolveIdentifiers: false
  });

  const ranked = result.layers.proteomics.rows.map((row,index) => ({...row,rank:index+1}));
  const her2 = ranked.find((row) => /^HER2$/i.test(row.feature));
  const er = ranked.find((row) => /^ER-alpha$/i.test(row.feature));

  requireTruth(
    result.layers.proteomics.mode === 'outcome-multiclass',
    'TCGA three-subtype phenotype uses the multiclass outcome branch',
    result.layers.proteomics.inferenceMethod
  );
  requireTruth(
    her2 && her2.qValue <= 0.01 && er && er.qValue <= 0.01,
    'TCGA multiclass outcome model recovers HER2 and ER-alpha subtype biology',
    'HER2 rank=' + (her2?.rank ?? 'NA') + ', q=' + (her2?.qValue ?? 'NA') + '; ER-alpha rank=' + (er?.rank ?? 'NA') + ', q=' + (er?.qValue ?? 'NA')
  );

  report.benchmarks.tcgaMulticlassOutcome = {
    truth: 'Basal / Her2 / LumA subtype modelled as a multiclass phenotype.',
    HER2: her2,
    ERalpha: er,
    top: result.layers.proteomics.rows.slice(0,15)
  };
}

// --- Public validation of the outcome branch: TCGA Her2 vs LumA treated as a binary phenotype ---
{
  const ds = await loadDataset(root + '/tcga-her2-luma', ['transcriptomics','proteomics']);
  const outcomeRows = ds.metadataRows.map((row) => ({ ...row, outcome: row.condition }));
  const result = await runDeterministicAnalysis({
    ...ds,
    metadataRows: outcomeRows,
    columnMapping: ds.mapping,
    protocol: {
      organism: 'human',
      objective: 'outcome',
      outcomeType: 'binary',
      longitudinal: false,
      designType: 'independent',
      studySetting: 'clinical_observational',
      groupCount: '1',
      sampleOverlap: 'same_specimen',
      batchKnown: 'no',
      covariateColumns: []
    },
    dataTypes: {
      transcriptomics: 'log_expression',
      proteomics: 'log_intensity',
      metabolomics: 'normalized'
    },
    useReactome: false,
    resolveIdentifiers: false
  });

  const her2 = result.layers.proteomics.rows
    .map((row,index) => ({...row,rank:index+1}))
    .find((row) => /^HER2$/i.test(row.feature));
  requireTruth(
    result.layers.proteomics.mode === 'outcome-binary',
    'TCGA binary phenotype uses the logistic outcome branch',
    result.layers.proteomics.inferenceMethod
  );
  requireTruth(
    her2 && her2.effect < 0 && her2.qValue <= 0.10 && her2.exponentiatedEffect < 1,
    'TCGA logistic outcome model recovers HER2 as associated with Her2-vs-LumA phenotype',
    her2 ? 'rank=' + her2.rank + ', beta=' + her2.effect.toFixed(3) + ', OR=' + her2.exponentiatedEffect.toPrecision(3) + ', q=' + her2.qValue : 'HER2 absent'
  );

  report.benchmarks.tcgaBinaryOutcome = {
    truth: 'Her2-enriched vs Luminal A phenotype modelled as a binary outcome rather than as a group contrast.',
    HER2: her2,
    transcriptomicsTop: result.layers.transcriptomics.rows.slice(0,15),
    proteomicsTop: result.layers.proteomics.rows.slice(0,15)
  };
}

await fs.writeFile(`${root}/benchmark-report.json`, JSON.stringify(report, null, 2));
console.log('PUBLIC MULTI-OMICS BENCHMARK SUITE: PASS');
