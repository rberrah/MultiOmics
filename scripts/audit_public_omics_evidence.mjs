// A *scope-bound* qualification gate for real public-omics fixtures.
// This audits provenance and existing truth benchmarks, it does not invent
// a cross-cohort prospective replication or a clinical validation claim.
import fs from 'node:fs/promises';

const REQUIRED_SOURCES = Object.freeze({
  mixOmics: 'ef3e760526623d9e91e945e4b50be48d764efc40',
  NCI60_GeneMetabolite_Data: '4c140708c18c41910458414a626b1a51d3dd2a85',
  IntLIMVignettes: 'e1a03362dd113db9655ad3e4dabd62e1c47cd731',
  xOmicsShiny: 'db2065b87c6e7a3d76e24e272573dd74428ff191',
  PaintOmics: 'a166da377450977a6a95ec44f4ae703a31c99742',
  missRows: '5c5fd6c8940523bcf1511fafe9e7896a5b8a0a1e'
});

function assert(condition,message) {
  if(!condition)throw new Error('PUBLIC OMICS QUALIFICATION BLOCKED: '+message);
}
function plausibleSha(value) {return typeof value==='string' && /^[a-f0-9]{64}$/i.test(value);}
function plausibleGitSha(value) {return typeof value==='string' && /^[a-f0-9]{40}$/i.test(value);}
function finiteNumber(x){return typeof x==='number'&&Number.isFinite(x);}

export function auditPublicOmicsEvidence(report, provenance) {
  assert(report&&typeof report==='object'&&report.benchmarks&&
    typeof report.benchmarks==='object','No public truth benchmark report.');
  assert(provenance&&typeof provenance==='object'&&
    provenance.schema==='pmx-explain-public-multiomics-benchmark-provenance/v1',
    'Provenance manifest absent or wrong schema.');
  assert(plausibleGitSha(provenance.pmxExplain?.commit),
    'Current analysis commit SHA must be recorded.');
  for(const [source,commit] of Object.entries(REQUIRED_SOURCES)){
    const record=provenance.sourceRepositories?.[source];
    assert(record&&record.commit===commit,
      'Public dataset '+source+' has unknown/unreviewed git revision.');
  }
  const files=provenance.exportedInputFiles;
  assert(files&&typeof files==='object'&&Object.keys(files).length>=12,
    'Insufficient hashed exported public dataset inputs.');
  for(const [name,item] of Object.entries(files)){
    assert(!name.includes('..') && !name.startsWith('/'),
      'Unsafe relative exported file name.');
    assert(plausibleSha(item?.sha256) && Number.isInteger(item.bytes)&&item.bytes>0,
      'Missing SHA256 or byte count for '+name);
  }
  const benchmarks=report.benchmarks;
  const coverage=[
    {id:'nutrimouse',source:'mixOmics',layers:['transcriptomics','metabolomics'],
      verified:!!benchmarks.nutrimouse?.truth && !!benchmarks.nutrimouse?.transcriptomics &&
       !!benchmarks.nutrimouse?.metabolomics,
      boundary:'Public mouse RNA/lipid phenotype recovery; no independent model validation.'},
    {id:'tcgaHer2LumA',source:'mixOmics',layers:['transcriptomics','proteomics'],
      verified:!!benchmarks.tcgaHer2LumA?.truth &&
        Number(benchmarks.tcgaHer2LumA?.expectedDirectionHits)>=2,
      boundary:'Public breast-tumour subtype recovery; no held-out prospective confirmation.'},
    {id:'agingHFCD',source:'xOmicsShiny',
      layers:['transcriptomics','proteomics','metabolomics'],
      verified:['transcriptomics','proteomics','metabolomics'].every(layer=>{
        const x=benchmarks.agingHFCD?.concordance?.[layer];
        return Number.isInteger(x?.n)&&x.n>=10 &&
          finiteNumber(x.fraction)&&x.fraction>=.80;
      }),
      boundary:'Agreement with published processed-data effect directions, not raw RNA-seq statistical validation.'},
    {id:'stategraReplicates',source:'PaintOmics',layers:['metabolomics'],
      verified:Number(benchmarks.stategraReplicates?.samples)===36 &&
        Number(benchmarks.stategraReplicates?.conditionTimeCells)===12 &&
        finiteNumber(benchmarks.stategraReplicates?.directionConcordance) &&
        benchmarks.stategraReplicates.directionConcordance>=.95,
      boundary:'Replicated biological design and qualitative direction; public processed representations have scale differences.'},
    {id:'nci60IntLIM',source:'NCI60_GeneMetabolite_Data',
      layers:['transcriptomics','metabolomics'],
      verified:finiteNumber(benchmarks.nci60IntLIM?.FAM174B_malic?.deltaR) &&
        benchmarks.nci60IntLIM.FAM174B_malic.deltaR>.5,
      boundary:'Published gene/metabolite association recovery, not causality or independent patient replication.'}
  ];
  for(const x of coverage)assert(x.verified,
    'Required publicly sourced truth evidence missing or regression: '+x.id);
  const negative=report.negativeControls?.tcgaHer2LumALabelPermutation;
  assert(Array.isArray(negative?.permutations)&&negative.permutations.length>=5,
    'Subject-level negative controls were not executed.');
  assert(finiteNumber(negative?.summary?.meanAuc)&&
         negative.summary.meanAuc>=.35&&negative.summary.meanAuc<=.65,
    'Negative-label outcome benchmark fails its prespecified chance-range.');
  assert(finiteNumber(negative?.summary?.medianQ05Fraction)&&
         negative.summary.medianQ05Fraction<=.02,
    'Negative-label feature FDR screen fails the prespecified threshold.');
  const distinctSources=[...new Set(coverage.map(x=>x.source))];
  const modalities=[...new Set(coverage.flatMap(x=>x.layers))];
  assert(distinctSources.length>=4&&modalities.length===3,
    'Missing independently sourced public coverage for RNA/protein/metabolite layers.');
  return {
    schema:'pmx-multiomics-public-evidence-gate/v1',
    status:'public_reference_recovery_and_provenance_pass',
    datasets:coverage,
    distinct_source_repositories:distinctSources.length,
    modalities,
    provenance_commit:provenance.pmxExplain.commit,
    total_hashed_exported_files:Object.keys(files).length,
    negative_control_runs:negative.permutations.length,
    negative_control_precision:'exploratory_five_or_more_permutations_not_formal_calibration',
    independent_frozen_prediction_studies_verified_here:0,
    raw_rnaseq_confirmatory_studies_verified_here:0,
    clinical_or_universal_certification:false,
    limitations:'Public biological signal recovery uses published processed data and limited controls. Separate prospective or cross-cohort frozen-model validation and raw technology-specific benchmarking are still required.'
  };
}
if(process.argv[1]?.endsWith('audit_public_omics_evidence.mjs')){
  const root=process.argv[2]||'tmp/public-benchmarks';
  const report=JSON.parse(await fs.readFile(root+'/benchmark-report.json','utf8'));
  const provenance=JSON.parse(await fs.readFile(root+'/benchmark-provenance.json','utf8'));
  const audit=auditPublicOmicsEvidence(report,provenance);
  await fs.writeFile(root+'/qualification-scope.json',JSON.stringify(audit,null,2)+'\n');
  console.log('PUBLIC OMICS SCIENTIFIC SCOPE GATE: PASS (not clinical certification)');
}
