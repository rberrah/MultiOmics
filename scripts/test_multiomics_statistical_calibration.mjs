import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';

// Reproducible, independent simulation. This is a limited calibration screen,
// NOT evidence that every model / RNA-seq / MS / longitudinal design is valid.
function rng(seed) {
  let state = seed >>> 0;
  return () => ((state = (Math.imul(1664525, state) + 1013904223) >>> 0) + 0.5) / 4294967296;
}
function normal(random) {
  const u = Math.max(1e-12, random());
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * random());
}
function fixture(seed, planted = false, collinear = false) {
  const random = rng(seed);
  const ids = Array.from({length:24}, (_,i)=>'S'+String(i+1).padStart(2,'0'));
  const features = Array.from({length:24},(_,i)=>'P'+String(i+1).padStart(3,'0'));
  const metadataRows = ids.map((id,i)=>({
    subject_id:id, sample_id:id, assay_id:id, omic:'proteomics',
    condition:i<12?'control':'treated', timepoint:'T0', batch:'',
    covariate_group_proxy:i<12?'0':'1'
  }));
  const vals = features.map((feature,j)=>ids.map((_,i)=>
    8 + j * 0.1 + normal(random) + (planted && j<4 && i>=12 ? 2.2 : 0)));
  const csv = ['feature_id,'+ids.join(','),
    ...features.map((feature,j)=>[feature,...vals[j].map(v=>v.toFixed(10))].join(','))].join('\n');
  const files = {
    metadata:new File(['placeholder'],'samples.csv',{type:'text/csv'}),
    proteomics:new File([csv],'proteomics.csv',{type:'text/csv'}),
    transcriptomics:null,metabolomics:null
  };
  const columnMapping = {subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',omic:'omic',
    condition:'condition',timepoint:'timepoint',batch:'batch',covariate_group_proxy:'covariate_group_proxy'};
  const protocol = {organism:'human',objective:'groups',longitudinal:false,designType:'independent',
    studySetting:'clinical_observational',groupCount:'2',covariateColumns:collinear?['covariate_group_proxy']:[]};
  return { ids, features, metadataRows, vals, files, columnMapping, protocol };
}
async function analyze(data) {
  return runDeterministicAnalysis({
    files:data.files,metadataRows:data.metadataRows,columnMapping:data.columnMapping,
    protocol:data.protocol,
    dataTypes:{transcriptomics:'raw_counts',proteomics:'log_intensity',metabolomics:'peak_area'},
    useReactome:false,resolveIdentifiers:false
  });
}

if (process.argv.includes('--reference-input')) {
  const f = fixture(817341);
  console.log('subject_id,condition,feature,value');
  for (let j=0;j<f.features.length;j++) for (let i=0;i<f.ids.length;i++)
    console.log([f.ids[i], i<12?'control':'treated',f.features[j],f.vals[j][i].toFixed(10)].join(','));
} else if (process.argv.includes('--reference-results')) {
  const out = await analyze(fixture(817341));
  console.log('feature,effect,se,p_value,q_value');
  for (const row of out.layers.proteomics.rows)
    console.log([row.feature,row.effect,row.standardError,row.pValue,row.qValue].join(','));
} else {
  // Null with independent Gaussian noise and no batch/covariate confounding:
  // check nominal rejection, BH false discoveries, and invariant constraints.
  const repeats=60, qThreshold=0.10;
  let tests=0, rejections=0, bhFamilies=0, invalid=0, powerHits=0;
  for(let trial=0;trial<repeats;trial++) {
    const out=await analyze(fixture(10000+trial*13));
    const rows=out.layers.proteomics.rows;
    assert.equal(rows.length,24,'all 24 generated null features should be retained');
    tests+=rows.length;
    let anyBH=false;
    for(const row of rows) {
      assert(Number.isFinite(row.pValue) && row.pValue>=0 && row.pValue<=1,'valid null p');
      assert(Number.isFinite(row.qValue) && row.qValue>=0 && row.qValue<=1,'valid null BH q');
      if(row.pValue<0.05) rejections++;
      if(row.qValue<=qThreshold) anyBH=true;
      if(!Number.isFinite(row.effect) || !Number.isFinite(row.standardError)) invalid++;
    }
    if(anyBH)bhFamilies++;
    if(trial<30) {
      const planted=await analyze(fixture(50000+trial*17,true));
      const signal=new Set(['P001','P002','P003','P004']);
      for(const row of planted.layers.proteomics.rows) {
        if(signal.has(row.feature) && row.effect>0 && row.qValue<=qThreshold)powerHits++;
      }
    }
  }
  const nullRate=rejections/tests;
  const familyRate=bhFamilies/repeats;
  const power=powerHits/(30*4);
  // These are *predeclared tolerant smoke thresholds*, not acceptance of
  // scientific calibration at arbitrary sample sizes or distributions.
  assert.equal(invalid,0,'finite effects and standard errors');
  assert(nullRate>=0.02 && nullRate<=0.085,'grossly miscalibrated null p-value rate '+nullRate);
  assert(familyRate<=0.25,'grossly inflated BH family rejection '+familyRate);
  assert(power>=0.65,'loss of positive control detection '+power);

  // A perfectly aliased covariate must not produce falsely precise p-values
  // from the internal numerical ridge.
  const aliased=await analyze(fixture(817341,false,true));
  assert.equal(aliased.layers.proteomics.rows.length,0,
    'group inference must refuse an exactly collinear design');
  console.log(JSON.stringify({
    kind:'prespecified_limited_calibration',seed:10000,repeats,
    sampleSize:24,featuresPerDataset:24,nominalNullTests:tests,
    observedNullPBelow005:rejections,observedNullRate:nullRate,
    familiesWithAnyBHDiscovery:bhFamilies,observedFamilyRate:familyRate,
    plantedSignals:120,plantedBHHits:powerHits,plantedSensitivity:power,
    rankDeficiency:'refused',limitations:'One Gaussian independent two-group modality; not all distributions, study designs or R package equivalence.'
  },null,2));
}
