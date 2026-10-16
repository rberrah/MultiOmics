// Empirical calibration grid for two independent groups and one log-intensity
// modality. This is a software regression screen, NOT proof of nominal 5%
// calibration, FDR control or general validity on clinical cohorts.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import { File } from 'node:buffer';
import { runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';

function randomSource(seed) {
  let x = seed >>> 0;
  return () => ((x = (Math.imul(x,1664525)+1013904223)>>>0) + 0.5)/4294967296;
}
function normal(random) {
  const a = Math.max(1e-12,random());
  return Math.sqrt(-2*Math.log(a))*Math.cos(2*Math.PI*random());
}
function fixture(n,model,seed,planted=false) {
  const rand=randomSource(seed),features=12;
  const ids=Array.from({length:n},(_,i)=>'P'+String(i+1).padStart(3,'0'));
  const metadataRows=ids.map((id,i)=>({
    subject_id:id,sample_id:id,assay_id:id,omic:'proteomics',
    condition:i<n/2?'control':'treated',
    timepoint:'T0',batch:i%2?'B2':'B1',
    age:String(20+Math.floor(55*rand()))
  }));
  const lines=[['feature_id',...ids].join(',')];
  for(let j=0;j<features;j++) {
    const measurements=ids.map((id,i)=>{
      const treated=i>=n/2;
      let eps=normal(rand);
      if(model==='heteroscedastic_mcar') eps*=treated?1.9:0.75;
      const shift=planted && j<3 && treated?2.7:0;
      // A large technical-only shift must not masquerade as a group effect.
      const batch=metadataRows[i].batch==='B2'?(model==='batch_dominant_null'?4:0.5):0;
      const age=(Number(metadataRows[i].age)-47)*0.015;
      const missing=model==='heteroscedastic_mcar' && rand()<0.12;
      return missing?'':(8+j*0.25+batch+age+shift+eps).toFixed(8);
    });
    lines.push(['G'+String(j+1).padStart(2,'0'),...measurements].join(','));
  }
  return {
    files:{
      metadata:new File(['fixture'],'metadata.csv',{type:'text/csv'}),
      transcriptomics:null,
      proteomics:new File([lines.join('\n')],'proteomics.csv',{type:'text/csv'}),
      metabolomics:null
    },
    metadataRows,
    columnMapping:{
      subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',
      omic:'omic',condition:'condition',timepoint:'timepoint',batch:'batch',
      age:'age'
    },
    protocol:{
      objective:'groups',organism:'human',designType:'independent',
      longitudinal:false,studySetting:'synthetic_test',
      groupCount:'2',batchKnown:'yes',covariateColumns:['age']
    }
  };
}
// Wilson interval on independent *datasets* with any BH discovery.
// Feature-level p-values within a dataset are correlated and are NOT
// counted as independent Bernoulli replicates for this interval.
function wilson95(successes,n) {
  if(!n)return null;
  const z=1.959963984540054,p=successes/n,d=1+z*z/n;
  const center=(p+z*z/(2*n))/d;
  const margin=z*Math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d;
  return {lower:Math.max(0,center-margin),upper:Math.min(1,center+margin)};
}
async function analyse(input) {
  return runDeterministicAnalysis({
    ...input,dataTypes:{proteomics:'log_intensity'},
    resolveIdentifiers:false,useReactome:false
  });
}
const configs=[
  ...[16,40,80].flatMap(n=>[
    {n,model:'gaussian_balanced'},
    {n,model:'heteroscedastic_mcar'},
    {n,model:'batch_dominant_null'}
  ])
];
const nullTrials=20,signalTrials=6,results=[];
for(const {n,model} of configs) {
  let totalNullP=0,rawFalsePositives=0,empty=0,invalid=0;
  let qFamilies=0,signalEvaluated=0,signalDetected=0;
  for(let trial=0;trial<nullTrials;trial++) {
    const result=await analyse(fixture(n,model,90001+n*100+trial*17));
    const rows=result.layers.proteomics.rows;
    if(!rows.length)empty++;
    if(rows.some(x=>!Number.isFinite(x.pValue)||!Number.isFinite(x.qValue)))invalid++;
    let any=false;
    for(const row of rows) {
      totalNullP++;
      if(row.pValue<0.05)rawFalsePositives++;
      if(row.qValue<=0.05)any=true;
    }
    if(any)qFamilies++;
  }
  for(let trial=0;trial<signalTrials;trial++) {
    const result=await analyse(fixture(n,model,170001+n*100+trial*13,true));
    for(const row of result.layers.proteomics.rows) {
      if(!['G01','G02','G03'].includes(row.feature))continue;
      signalEvaluated++;
      if(row.qValue<=0.10 && row.effect>0)signalDetected++;
    }
  }
  assert.equal(empty,0,'all null scenarios must fit');
  assert.equal(invalid,0,'no invalid p/q-values in null scenarios');
  assert.ok(totalNullP>=nullTrials*9,'unexpected feature attrition');
  const typeIRate=rawFalsePositives/totalNullP;
  const bhFamilyRate=qFamilies/nullTrials;
  const sensitivity=signalDetected/Math.max(1,signalEvaluated);
  // Predeclared, conservative guardrails to detect GROSS regression only.
  // These are not equivalence tests to 5% or a validation of BH assumptions.
  assert.ok(typeIRate<=0.16,
    'large false-positive inflation in '+model+' N='+n+': '+typeIRate);
  assert.ok(bhFamilyRate<=0.40,
    'large false-discovery family rate in '+model+' N='+n+': '+bhFamilyRate);
  if(model==='batch_dominant_null')assert.ok(bhFamilyRate<=0.20,
    'batch-only negative control generated too many BH biological discoveries N='+n+': '+bhFamilyRate);
  if(n>=40)assert.ok(sensitivity>=0.40,
    'loss of planted effect detection in '+model+' N='+n+': '+sensitivity);
  results.push({
    n,model,
    nullTrials,validNullHypotheses:totalNullP,
    empiricalPBelow005:rawFalsePositives,
    rawTypeIRate:typeIRate,
    familiesWithAnyQ005:qFamilies,
    bhFamilyRate,
    bhFamilyWilson95: wilson95(qFamilies,nullTrials),
    plantedEffectEvaluated:signalEvaluated,
    plantedQ010Detection:signalDetected,
    plantedSensitivity:sensitivity
  });
}
const report={
  benchmark:'multiomics_independent_group_calibration_grid_v1',
  prespecifiedSeed:90001,
  outcomes:'two-group adjusted HC3 model with batch and age',
  nValues:[16,40,80],
  distributions:['gaussian_balanced','heteroscedastic_mcar','batch_dominant_null'],
  nominalAlpha:0.05,
  correction:'within-omic Benjamini-Hochberg',
  intervalMethod:'Wilson 95% on independent null datasets with any q<=0.05, NOT on pooled correlated features',
  guardrails:'gross regression: observed raw p<0.05 <=16%, null BH-family discovery <=40%, planted signal detection >=40% for N>=40',
  limitations:'180 all-null datasets and 54 planted datasets; 12 genes/proteins per dataset. Technical-only null uses group-balanced batches with a fixed +4-unit batch effect. Wilson intervals quantify uncertainty in the 20 independent null families per scenario, but cannot prove nominal 5% FDR. No confidence-interval equivalence test, formal FDR validation, other designs, covariate MNAR, publication approval or clinical external validation.',
  results
};
if (process.env.MULTIOMICS_CALIBRATION_REPORT_PATH) {
  const path=process.env.MULTIOMICS_CALIBRATION_REPORT_PATH;
  const dir=path.split('/').slice(0,-1).join('/');
  if(dir)await fs.mkdir(dir,{recursive:true});
  await fs.writeFile(path,JSON.stringify(report,null,2)+'\n');
}
console.log(JSON.stringify(report,null,2));
console.log('multiomics independent-group calibration grid: PASS');
