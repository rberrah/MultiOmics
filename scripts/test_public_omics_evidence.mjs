import assert from 'node:assert/strict';
import {auditPublicOmicsEvidence} from './audit_public_omics_evidence.mjs';
const revisions={
 mixOmics:'ef3e760526623d9e91e945e4b50be48d764efc40',
 NCI60_GeneMetabolite_Data:'4c140708c18c41910458414a626b1a51d3dd2a85',
 IntLIMVignettes:'e1a03362dd113db9655ad3e4dabd62e1c47cd731',
 xOmicsShiny:'db2065b87c6e7a3d76e24e272573dd74428ff191',
 PaintOmics:'a166da377450977a6a95ec44f4ae703a31c99742',
 missRows:'5c5fd6c8940523bcf1511fafe9e7896a5b8a0a1e'
};
const sha='a'.repeat(64);
const provenance={
 schema:'pmx-explain-public-multiomics-benchmark-provenance/v1',
 pmxExplain:{commit:'a'.repeat(40)},
 sourceRepositories:Object.fromEntries(Object.entries(revisions).map(([k,commit])=>[k,{commit}])),
 exportedInputFiles:Object.fromEntries(Array.from({length:13},(_,i)=>['input'+i+'.csv',{sha256:sha,bytes:1}]))
};
const bench={
 nutrimouse:{truth:'mouse',transcriptomics:{top:[]},metabolomics:{top:[]}},
 tcgaHer2LumA:{truth:'Her2',expectedDirectionHits:3},
 agingHFCD:{concordance:Object.fromEntries(['transcriptomics','proteomics','metabolomics'].map(x=>[x,{n:20,fraction:.9}]))},
 stategraReplicates:{samples:36,conditionTimeCells:12,directionConcordance:.97},
 nci60IntLIM:{FAM174B_malic:{deltaR:.8}}
};
const report={benchmarks:bench,negativeControls:{tcgaHer2LumALabelPermutation:{
  permutations:Array.from({length:8},(_,i)=>({seed:i})),
  summary:{meanAuc:.48,medianQ05Fraction:.01}
}}};
const output=auditPublicOmicsEvidence(report,provenance);
assert.equal(output.datasets.length,5);
assert.equal(output.distinct_source_repositories,4);
assert.equal(output.independent_frozen_prediction_studies_verified_here,0);
assert.equal(output.clinical_or_universal_certification,false);
assert.deepEqual(output.modalities.sort(),['metabolomics','proteomics','transcriptomics']);
const reject=(r=report,p=provenance)=>assert.throws(()=>auditPublicOmicsEvidence(r,p),/QUALIFICATION BLOCKED/);
const cloned=(x)=>JSON.parse(JSON.stringify(x));
let bad=cloned(report);delete bad.benchmarks.agingHFCD;reject(bad);
bad=cloned(report);bad.benchmarks.stategraReplicates.directionConcordance=.3;reject(bad);
bad=cloned(report);bad.negativeControls.tcgaHer2LumALabelPermutation.summary.meanAuc=.99;reject(bad);
bad=cloned(report);bad.negativeControls.tcgaHer2LumALabelPermutation.permutations=[];reject(bad);
let badP=cloned(provenance);badP.sourceRepositories.mixOmics.commit='b'.repeat(64);reject(report,badP);
badP=cloned(provenance);delete badP.exportedInputFiles['input0.csv'];delete badP.exportedInputFiles['input1.csv'];reject(report,badP);
badP=cloned(provenance);badP.pmxExplain.commit='unknown';reject(report,badP);
console.log('PUBLIC OMICS EVIDENCE SCOPE AND SOURCE-PINNING TEST: PASS');
