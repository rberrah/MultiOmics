// Independent synthetic stress study for a narrow two-independent-group,
// log-intensity proteomics design. Not a proof of nominal FDR or medical validity.
// MCAR/MAR/MNAR are known only to this simulator, NEVER inferred from observed data.
import fs from 'node:fs/promises';
import {File} from 'node:buffer';
import {runDeterministicAnalysis} from '../src/lib/multiomics/deterministic.js';
import {evaluateConfirmatoryReadiness,assessScientificAssurance} from '../src/lib/multiomics/scientific-assurance.js';

const SCENARIOS=['complete','mcar','mar_observed_age','strong_batch','mnar_opposite_tails'];
const SAMPLE_SIZES=[24,80];
const NULL_REPS=200;
const SIGNAL_REPS=25;
const FEATURE_COUNT=12;
const ALPHA=0.05;
const SIGNAL_Q=0.10;

function rng(seed) {
  let x=seed>>>0;
  return ()=>((x=(Math.imul(x,1664525)+1013904223)>>>0)+0.5)/4294967296;
}
function normal(rand) {
  return Math.sqrt(-2*Math.log(Math.max(1e-12,rand())))*Math.cos(2*Math.PI*rand());
}
function wilson(k,n) {
  if (!n) return null;
  const z=1.959963984540054,p=k/n,den=1+z*z/n;
  const mid=(p+z*z/(2*n))/den;
  const half=z*Math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den;
  return {lower:Math.max(0,mid-half),upper:Math.min(1,mid+half)};
}
function fixture({scenario,n,seed,planted}) {
  const rand=rng(seed);
  const ids=Array.from({length:n},(_,i)=>'P'+String(i+1).padStart(3,'0'));
  // Technical batches are balanced between groups (n is divisible by 4).
  const metadataRows=ids.map((id,i)=>({
    subject_id:id,sample_id:id,assay_id:id,omic:'proteomics',
    condition:i<n/2?'control':'treated',
    timepoint:'T0',batch:i%2?'B2':'B1',
    age:String(20+Math.floor(55*rand()))
  }));
  let missing=0,possible=0,missingCtrl=0,missingTreated=0;
  const csv=['feature_id,'+ids.join(',')];
  for(let j=0;j<FEATURE_COUNT;j++) {
    const row=[];
    for(let i=0;i<n;i++) {
      const meta=metadataRows[i],treated=i>=n/2;
      const eps=normal(rand);
      const batch=meta.batch==='B2'?(scenario==='strong_batch'?4:0.5):0;
      const age=(Number(meta.age)-47)*0.015;
      const truth=8+j*0.25+batch+age+eps+(planted && j<3 && treated?2.7:0);
      // MCAR: independent of all measured/unmeasured values.
      // MAR: depends on *observed* age, included as model covariate.
      // MNAR: opposite residual tails are censored within the two groups.
      // Missing rates are equal IN EXPECTATION, not forced to match by hand.
      let hide=false;
      if(scenario==='mcar')hide=rand()<0.12;
      if(scenario==='mar_observed_age')
        hide=rand()<(Number(meta.age)>=48?0.26:0.04);
      if(scenario==='mnar_opposite_tails')
        hide=(treated?eps>0:eps<0)&&rand()<0.46;
      possible++;
      if(hide) {
        missing++;
        if(treated)missingTreated++;else missingCtrl++;
      }
      row.push(hide?'':truth.toFixed(8));
    }
    csv.push(['F'+String(j+1).padStart(2,'0'),...row].join(','));
  }
  return {
    files:{
      metadata:new File(['benchmark'],'metadata.csv',{type:'text/csv'}),
      transcriptomics:null,
      proteomics:new File([csv.join('\n')],'proteomics.csv',{type:'text/csv'}),
      metabolomics:null
    },
    metadataRows,
    columnMapping:{
      subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',
      omic:'omic',condition:'condition',timepoint:'timepoint',
      batch:'batch',age:'age'
    },
    protocol:{
      objective:'groups',organism:'human',designType:'independent',
      longitudinal:false,studySetting:'synthetic_test',groupCount:'2',
      batchKnown:'yes',covariateColumns:['age']
    },
    simulation:{
      scenario,planted,n,missingFraction:missing/possible,
      missingByGroup:{control:missingCtrl/(n/2*FEATURE_COUNT),
        treated:missingTreated/(n/2*FEATURE_COUNT)}
    }
  };
}
async function analyse(f) {
  return runDeterministicAnalysis({
    files:f.files,metadataRows:f.metadataRows,
    columnMapping:f.columnMapping,protocol:f.protocol,
    dataTypes:{proteomics:'log_intensity'},
    resolveIdentifiers:false,useReactome:false
  });
}
const results=[];
const regressionFailures=[];
const gate=(condition,message)=>{ if(!condition)regressionFailures.push(message); };
for(const n of SAMPLE_SIZES)for(const scenario of SCENARIOS) {
  let valid=0,invalid=0,empty=0,rowsTested=0,rawReject=0,bhAny=0;
  let missingSum=0,imbalanceSum=0,missingReviewed=0,missingBlocked=0,exploratoryWarned=0;
  let signals=0,detected=0;
  for(let trial=0;trial<NULL_REPS;trial++){
    const seed=420001+100000*n+2009*SCENARIOS.indexOf(scenario)+17*trial;
    const f=fixture({n,scenario,seed,planted:false});
    const out=await analyse(f);
    const rows=out.layers?.proteomics?.rows||[];
    const bad=rows.some(row=>!Number.isFinite(row.pValue)||!Number.isFinite(row.qValue)
      ||row.pValue<0||row.pValue>1||row.qValue<0||row.qValue>1);
    if(!rows.length)empty++;
    if(bad)invalid++;else valid++;
    rowsTested+=rows.length;
    let any=false;
    for(const row of rows){
      if(row.pValue<ALPHA)rawReject++;
      if(row.qValue<=ALPHA)any=true;
    }
    if(any)bhAny++;
    missingSum+=f.simulation.missingFraction;
    imbalanceSum+=Math.abs(f.simulation.missingByGroup.control-f.simulation.missingByGroup.treated);
    if(scenario==='mnar_opposite_tails') {
      const exploratory=assessScientificAssurance(out);
      if(exploratory.notes.some(x=>x.code==='proteomics_missingness_not_ignorable'))
        exploratoryWarned++;
      // Granting a *hypothetical* successful R fit must not certify MNAR data.
      // The checker cannot identify the mechanism; it must at least flag
      // nonzero observed missingness and require independent review.
      const candidate={
        ...out,
        protocol:{...out.protocol,analysisIntent:'confirmatory'},
        referenceBackend:{status:'ok',methods:{proteomics_differential:{status:'ok'}}}
      };
      const readiness=evaluateConfirmatoryReadiness(candidate);
      const warns=readiness.checks.some(x=>
        x.code==='missingness_mechanism_proteomics'
        ||x.code==='high_missingness_proteomics'
        ||x.code==='missingness_imbalance_proteomics');
      if(warns)missingReviewed++;
      if(readiness.status==='blocked')missingBlocked++;
      gate(readiness.certified===false,'MNAR sensitivity must never be certified');
    }
  }
  for(let trial=0;trial<SIGNAL_REPS;trial++) {
    const f=fixture({n,scenario,seed:770001+3000*n+331*SCENARIOS.indexOf(scenario)+29*trial,planted:true});
    const out=await analyse(f);
    const hit=new Map((out.layers?.proteomics?.rows||[]).map(x=>[x.feature,x]));
    for(const name of ['F01','F02','F03']){
      signals++;
      if(hit.get(name)?.effect>0 && hit.get(name)?.qValue<=SIGNAL_Q)detected++;
    }
  }
  const familyRate=bhAny/NULL_REPS;
  const rawRate=rowsTested?rawReject/rowsTested:null;
  const meanMissing=missingSum/NULL_REPS;
  const meanImbalance=imbalanceSum/NULL_REPS;
  gate(empty===0,'null analysis lost all features in '+scenario+'/'+n);
  gate(invalid===0,'invalid p/q-values for '+scenario+'/'+n);
  gate(rowsTested>=NULL_REPS*9,'unexpected feature attrition for '+scenario+'/'+n);
  if(scenario==='complete'||scenario==='strong_batch')gate(meanMissing===0,'complete-data simulation has missingness '+scenario+'/'+n);
  else gate(meanMissing>0.04,'missingness generator did not run '+scenario);
  // These wide ceilings detect egregious regression, not a certification
  // that p-values or BH are calibrated at 5% in a specified target population.
  if(scenario!=='mnar_opposite_tails'){
    gate(rawRate<=0.11,'excessive raw null rejection '+scenario+'/'+n+': '+rawRate);
    gate(familyRate<=0.17,'excessive null BH family discovery '+scenario+'/'+n+': '+familyRate);
  }
  if(scenario==='mnar_opposite_tails'){
    gate(meanImbalance<0.10,'missingness should be comparable between groups in expectation');
    gate(exploratoryWarned===NULL_REPS,
      'Every MNAR exploratory output must warn about non-ignorable missingness');
    gate(missingReviewed>=Math.floor(NULL_REPS*0.90),
      'MNAR censoring must trigger visible missingness review in most studies');
  }
  const power=detected/signals;
  if(n===80 && scenario!=='mnar_opposite_tails')
    gate(power>=0.55,'positive control detection lost '+scenario+': '+power);
  results.push({
    n,scenario,nullDatasets:NULL_REPS,plantedDatasets:SIGNAL_REPS,
    featuresPerDataset:FEATURE_COUNT,validNullDatasets:valid,
    invalidNullDatasets:invalid,totalNullFeaturesTested:rowsTested,
    rawNullRejections:rawReject,rawNullRejectionRate:rawRate,
    nullFamiliesWithAnyBHDiscovery:bhAny,
    nullBHFamilyDiscoveryRate:familyRate,
    nullBHFamilyWilson95:wilson(bhAny,NULL_REPS),
    meanMissingFraction:meanMissing,
    meanAbsoluteGroupMissingnessDifference:meanImbalance,
    mnarWithMissingnessReview:scenario==='mnar_opposite_tails'?missingReviewed:null,
    mnarExploratoryReportsWarned:scenario==='mnar_opposite_tails'?exploratoryWarned:null,
    mnarBlockedAtSoftwareGate:scenario==='mnar_opposite_tails'?missingBlocked:null,
    plantedFeatures:signals,positiveDiscoveries:detected,
    plantedSensitivityAtQ010:power,
    interpretiveTier:scenario==='mnar_opposite_tails'
      ?'adversarial_biased_observation_do_not_assert_nominal_calibration'
      :'gross_regression_screen_not_formal_validation'
  });
}
const report={
  benchmark:'multiomics_missingness_stress_v1',
  seedPolicy:'LCG fixed per scenario, sample size, and replicate',
  nominalAlpha:ALPHA,qPlanted:SIGNAL_Q,
  sampleSizes:SAMPLE_SIZES,scenarios:SCENARIOS,
  nullDatasets:NULL_REPS*SAMPLE_SIZES.length*SCENARIOS.length,
  signalDatasets:SIGNAL_REPS*SAMPLE_SIZES.length*SCENARIOS.length,
  interval:'Wilson 95% on INDEPENDENT null simulated datasets, not pooled features',
  caveats:[
    'Only one simulated log-intensity omics layer and two independent groups.',
    'True MCAR/MAR/MNAR mechanisms are known to the generator, not the software.',
    'A global-null BH family discovery is not a general proof of FDR control.',
    'Test thresholds are operational gross-regression screens, not scientific equivalence margins; no automatic certification.',
    'Positive-control strength is fixed, and power is not transportable to real cohorts.',
    'No causal claims, clinical validation, statistical preregistration or publication certificate.'
  ],
  results,
  regressionFailures
};
if(process.env.MULTIOMICS_MISSINGNESS_REPORT_PATH){
  const file=process.env.MULTIOMICS_MISSINGNESS_REPORT_PATH;
  await fs.mkdir(file.split('/').slice(0,-1).join('/')||'.',{recursive:true});
  await fs.writeFile(file,JSON.stringify(report,null,2)+'\n');
}
console.log(JSON.stringify(report,null,2));
if(regressionFailures.length){
  console.error('multiomics missingness stress benchmark: FAIL',regressionFailures.join(' | '));
  process.exitCode=1;
}else console.log('multiomics missingness stress benchmark: PASS');
