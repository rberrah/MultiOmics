// Empirical mixed-truth FDR validation, NOT a universal scientific certificate.
// Unlike all-null FWER screens, every dataset here contains known true
// alternatives AND known null features. Independent datasets are the
// resampling unit; p-values from correlated proteins are not independent.
import fs from 'node:fs/promises';
import {File} from 'node:buffer';
import {runDeterministicAnalysis} from '../src/lib/multiomics/deterministic.js';
import {assessScientificAssurance} from '../src/lib/multiomics/scientific-assurance.js';

const MODELS=[
  'correlated_gaussian','heteroscedastic_correlated','mcar_correlated',
  'mar_age_correlated','heavy_tail_t3','mnar_opposite_tails'
];
const SIZES=[24,80],REPEATS=200,FEATURES=12;
const TRUE_EFFECTS=[2.0,1.4,0.9];
const Q=0.05;
function generator(seed){
  // Independent XORShift32 seed family, unlike prior LCG stress generator.
  let x=(seed>>>0)||0x9e3779b9;
  return ()=>{x^=x<<13;x^=x>>>17;x^=x<<5;return (x>>>0)/4294967296;};
}
function gauss(rng){
  return Math.sqrt(-2*Math.log(Math.max(rng(),1e-12)))*Math.cos(2*Math.PI*rng());
}
function t3(rng){
  const z=gauss(rng),a=gauss(rng),b=gauss(rng),c=gauss(rng);
  return z/Math.sqrt(Math.max(1e-10,(a*a+b*b+c*c)/3))/Math.sqrt(3);
}
function bootstrap95(values,seed){
  if(values.length<2)return null;
  const next=generator(seed),means=[];
  for(let b=0;b<800;b++){
    let sum=0;
    for(let i=0;i<values.length;i++)
      sum+=values[Math.floor(next()*values.length)];
    means.push(sum/values.length);
  }
  means.sort((a,b)=>a-b);
  return {lower:means[Math.floor(.025*(means.length-1))],
    upper:means[Math.ceil(.975*(means.length-1))],
    resamples:means.length,unit:'independent simulated dataset'};
}
function makeStudy(model,n,seed){
  const next=generator(seed);
  const ids=Array.from({length:n},(_,i)=>'P'+String(i+1).padStart(3,'0'));
  const meta=ids.map((id,i)=>({
    subject_id:id,sample_id:id,assay_id:id,omic:'proteomics',
    condition:i<n/2?'control':'treated',timepoint:'T0',
    batch:i%2?'B2':'B1',age:String(20+Math.floor(55*next()))
  }));
  const common=ids.map(()=>gauss(next));
  const data=['feature_id,'+ids.join(',')];
  let missing=0,total=0;
  for(let j=0;j<FEATURES;j++){
    const f='F'+String(j+1).padStart(2,'0');
    const entries=[f];
    for(let i=0;i<n;i++){
      const treated=i>=n/2;
      const independent=model==='heavy_tail_t3'?t3(next):gauss(next);
      const eps=0.48*common[i]+Math.sqrt(1-.48**2)*independent;
      const scale=model==='heteroscedastic_correlated'?(treated?1.7:.75):1;
      const signal=treated?(TRUE_EFFECTS[j]||0):0;
      const batch=meta[i].batch==='B2'?.65:0;
      const age=(Number(meta[i].age)-47)*.013;
      const y=8+j*.13+signal+batch+age+eps*scale;
      let hide=false;
      if(model==='mcar_correlated')hide=next()<.12;
      if(model==='mar_age_correlated')
        hide=next()<(Number(meta[i].age)>=48?.24:.04);
      if(model==='mnar_opposite_tails')
        hide=(treated?eps>0:eps<0)&&next()<.5;
      if(hide)missing++;
      total++;
      entries.push(hide?'':y.toFixed(9));
    }
    data.push(entries.join(','));
  }
  return {
    files:{metadata:new File(['synthetic'],'metadata.csv',{type:'text/csv'}),
      transcriptomics:null,proteomics:new File([data.join('\n')],'proteomics.csv',{type:'text/csv'}),
      metabolomics:null},
    metadataRows:meta,
    columnMapping:{subject_id:'subject_id',sample_id:'sample_id',
      assay_id:'assay_id',omic:'omic',condition:'condition',
      timepoint:'timepoint',batch:'batch',age:'age'},
    protocol:{objective:'groups',organism:'human',designType:'independent',
      longitudinal:false,groupCount:'2',studySetting:'synthetic_mixed_truth',
      batchKnown:'yes',covariateColumns:['age']},
    dataTypes:{proteomics:'log_intensity'},
    resolveIdentifiers:false,useReactome:false,
    simulation:{missingFraction:missing/total,seed,model,n}
  };
}
function average(x){return x.reduce((a,b)=>a+b,0)/x.length;}
const report={
  benchmark:'mixed_truth_fdr_v1',generatedBy:'independent xorshift32 and gaussian/t3 simulators',
  comparedPipeline:'production deterministic single-omic log intensity group HC3 + BH',
  seedPolicy:'base 810999, independent scenario x sample size x replicate XORShift streams',
  nominalQ:Q,sampleSizes:SIZES,scenarios:MODELS,replicatesPerCell:REPEATS,
  featuresPerDataset:FEATURES,trueSignalsPerDataset:TRUE_EFFECTS.length,
  plantedEffects:TRUE_EFFECTS,
  fdrDefinition:'mean false discoveries / max(1, total discoveries) over independent datasets',
  uncertainty:'800 fixed-seed percentile bootstrap resamples of independent dataset FDP values',
  results:[],regressionFailures:[]
};
for(const n of SIZES)for(const [modelIndex,model] of MODELS.entries()){
  const fdp=[],power=[],nullCoverage=[],trueCoverage=[];
  let anyDiscovery=0,erroneous=0,retained=0,missing=0;
  let zeroDiscovery=0,missingWarned=0,unreliable=0;
  for(let rep=0;rep<REPEATS;rep++){
    const seed=810999+n*100007+modelIndex*1000003+rep*73939;
    const study=makeStudy(model,n,seed);
    const result=await runDeterministicAnalysis(study);
    const rows=result?.layers?.proteomics?.rows||[];
    const index=new Map(rows.map(r=>[r.feature,r]));
    if(index.size!==rows.length || !rows.length)
      report.regressionFailures.push('Empty or duplicated features: '+model+'/'+n+'/'+rep);
    if(rows.some(r=>!Number.isFinite(r.pValue)||!Number.isFinite(r.qValue)
      ||r.pValue<0||r.pValue>1||r.qValue<0||r.qValue>1))
      report.regressionFailures.push('Invalid p/q: '+model+'/'+n+'/'+rep);
    const discoveries=rows.filter(r=>r.qValue<=Q);
    const falseCalls=discoveries.filter(r=>Number(r.feature.slice(1))>TRUE_EFFECTS.length);
    const trueCalls=discoveries.filter(r=>Number(r.feature.slice(1))<=TRUE_EFFECTS.length
      && r.effect>0);
    const f=discoveries.length?falseCalls.length/discoveries.length:0;
    fdp.push(f);
    power.push(trueCalls.length/TRUE_EFFECTS.length);
    if(discoveries.length)anyDiscovery++;else zeroDiscovery++;
    if(falseCalls.length)erroneous++;
    retained+=rows.length;
    missing+=study.simulation.missingFraction;
    nullCoverage.push(rows.filter(r=>Number(r.feature.slice(1))>TRUE_EFFECTS.length)
      .map(r=>Number.isFinite(r.ciLow)&&Number.isFinite(r.ciHigh)&&r.ciLow<=0&&r.ciHigh>=0?1:0));
    trueCoverage.push(TRUE_EFFECTS.map((effect,i)=>{
      const r=index.get('F'+String(i+1).padStart(2,'0'));
      return r&&r.ciLow<=effect&&r.ciHigh>=effect?1:0;
    }));
    const assurance=assessScientificAssurance(result);
    if(study.simulation.missingFraction>0){
      if(assurance.notes.some(note=>note.code==='proteomics_missingness_not_ignorable'))missingWarned++;
    }
    if(model==='mnar_opposite_tails' &&
      !assurance.notes.some(note=>note.code==='proteomics_missingness_not_ignorable'))
      unreliable++;
  }
  const nullCoverPerDataset=nullCoverage.map(x=>x.length?average(x):0);
  const trueCoverPerDataset=trueCoverage.map(average);
  const meanFdr=average(fdp),interval=bootstrap95(fdp,140009+n*101+modelIndex*17);
  const risk=model==='mnar_opposite_tails'
    ?'unidentifiable_MNAR_not_eligible_for_confirmatory_inference'
    : interval.lower>Q
      ?'possible_fdr_inflation_review'
      :'no_clear_inflation_in_this_simulation_only';
  // This is a failure detector for gross software regressions only; the
  // 5% hypothesis is evaluated in the REPORT and never 'certified'.
  if(model!=='mnar_opposite_tails' && meanFdr>0.25)
    report.regressionFailures.push('Gross mixed-truth FDR inflation: '+model+'/'+n+' = '+meanFdr);
  if(model==='mnar_opposite_tails' && unreliable)
    report.regressionFailures.push('No warning in '+unreliable+' MNAR studies N='+n);
  const record={
    n,model,datasets:REPEATS,mixedTruth:true,knownAlternatives:TRUE_EFFECTS.length,
    knownNulls:FEATURES-TRUE_EFFECTS.length,
    meanRealisedFDR:meanFdr,bootstrap95:interval,
    meanTrueSignalRecall:average(power),
    meanNull95CIcoverage:average(nullCoverPerDataset),
    meanTrueSignal95CIcoverage:average(trueCoverPerDataset),
    anyDiscoveryDatasets:anyDiscovery,
    datasetsWithFalseDiscovery:erroneous,
    zeroDiscoveryDatasets:zeroDiscovery,
    meanRetainedFeatures:retained/REPEATS,
    meanMissingFraction:missing/REPEATS,
    missingnessWarningDatasets:missingWarned,
    evidenceStatus:risk
  };
  report.results.push(record);
  console.log('CELL '+model+' N='+n+' FDR='+meanFdr.toFixed(4)+
    ' 95% CI=['+interval.lower.toFixed(4)+','+interval.upper.toFixed(4)+
    '] true_recall='+record.meanTrueSignalRecall.toFixed(3));
}
report.limitations=[
  'Only simulated proteomics log intensity and independent 2-group study designs; NOT multi-omics/clinical validation.',
  'HC3 uses finite-sample approximations; nominal 95% intervals may be miscalibrated.',
  'The measured FDR depends on simulated effect sizes, true-null fraction, feature correlations and mechanisms.',
  'Feature filtering changes the tested BH family; report mean retained counts, do not assume a fixed denominator.',
  'Bootstrap intervals are Monte-Carlo uncertainty across simulated datasets, not biological transportability.',
  'Only one generator family per scenario and 200 independent replicates; new simulation models remain needed.',
  'MNAR censoring is unidentifiable from observed values; no compatible numerical reference makes it valid.',
  'No automatic publication approval or proof that study-wide FDR is controlled.'
];
const dest=process.env.MULTIOMICS_MIXED_FDR_REPORT_PATH;
if(dest){
  await fs.mkdir(dest.split('/').slice(0,-1).join('/')||'.',{recursive:true});
  await fs.writeFile(dest,JSON.stringify(report,null,2)+'\n');
}
if(report.regressionFailures.length){
  console.error('Mixed-truth FDR validation failures: '+report.regressionFailures.join(' | '));
  process.exitCode=1;
}else console.log('multiomics mixed-truth FDR calibration: PASS (gross regression check only)');
