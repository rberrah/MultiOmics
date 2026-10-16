#!/usr/bin/env node
/**
 * Paired software-method verification: identical simulated subjects go through
 * the PUBLIC browser pipeline and an INDEPENDENT R lm/HC3 implementation.
 *
 * Modes:
 *   --generate <directory>: writes raw subject/feature observations, browser
 *      estimates and protocol, WITHOUT using R for simulation or estimates.
 *   --compare <directory>: consumes independently calculated R estimates,
 *      enforces predeclared tolerances and writes machine-readable diagnostics.
 *
 * MNAR outcomes are intentional negative methodological controls, NEVER an
 * endorsed confirmatory missing-data analysis.
 */
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import { File } from 'node:buffer';
import { runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';

const SPEC = Object.freeze({
  version: 'independent_hc3_reference_v1',
  seed: 20261009,
  subjects: [24, 64],
  scenarios: ['gaussian', 'heteroscedastic_mcar', 'mnar_left_censor', 'mnar_group_asymmetric'],
  nullReplicates: 24,
  effectReplicates: 12,
  features: 10,
  signalFeatures: 3,
  plantedEffect: 1.3,
  nominalAlpha: 0.05,
  rNumericalTolerance: {
    effect: 3e-4, standardError: 3e-4, pValue: 3e-4,
    qValue: 3e-4, ciLow: 1e-3, ciHigh: 1e-3
  }
});
const csv = values=>values.map(x=>String(x??'')).join(',');
function readCsv(text) {
  const [head,...rows]=text.trim().split(/\r?\n/);
  const names=head.split(',');
  return rows.filter(Boolean).map(line=>{
    const values=line.split(',').map(value=>value.replace(/^"|"$/g,'').replace(/""/g,'"'));
    return Object.fromEntries(names.map((k,i)=>[
      k.replace(/^"|"$/g,''),values[i]??''
    ]));
  });
}
function random(seed) {
  let state=seed>>>0;
  return ()=>((state=(Math.imul(state,1664525)+1013904223)>>>0)+0.5)/4294967296;
}
function gaussian(next) {
  const u=Math.max(next(),1e-12),v=next();
  return Math.sqrt(-2*Math.log(u))*Math.cos(2*Math.PI*v);
}

async function generate(directory) {
  await fs.mkdir(directory,{recursive:true});
  const data=[csv(['case_id','scenario','n','replicate','effect_planted','feature','subject_id','condition','batch','age','value'])];
  const estimates=[csv(['case_id','scenario','n','replicate','effect_planted','feature','effect','standardError','pValue','qValue','ciLow','ciHigh'])];
  const cases=[];
  for(const n of SPEC.subjects)for(const scenario of SPEC.scenarios) {
    for(const planted of [false,true]){
      const reps=planted?SPEC.effectReplicates:SPEC.nullReplicates;
      for(let rep=0;rep<reps;rep++) {
        const caseId=[scenario,n,planted?'signal':'null',rep].join('_');
        const seed=(SPEC.seed+ n*100000 +SPEC.scenarios.indexOf(scenario)*10000
          +(planted?5000:0)+rep*109)>>>0;
        const next=random(seed);
        const ids=Array.from({length:n},(_,i)=>'S'+String(i+1).padStart(3,'0'));
        const metadataRows=ids.map((id,i)=>({
          subject_id:id,sample_id:id,assay_id:id,omic:'proteomics',
          condition:i<n/2?'control':'treated',
          timepoint:'T0',batch:i%2===0?'B1':'B2',
          age:String(20+Math.floor(55*next()))
        }));
        const allObservations=[], matrix=[csv(['feature_id',...ids])];
        for(let j=0;j<SPEC.features;j++) {
          const feature='PROT_'+String(j+1).padStart(2,'0');
          const row=[feature];
          for(let i=0;i<n;i++){
            const treated=i>=n/2;
            const error=gaussian(next);
            const sd=scenario==='heteroscedastic_mcar'?(treated?1.8:0.7):1;
            const truth=planted&&j<SPEC.signalFeatures&&treated?SPEC.plantedEffect:0;
            const age=Number(metadataRows[i].age);
            const value=8+j*0.1+truth+(metadataRows[i].batch==='B2'?0.45:0)
              +(age-47)*0.014+sd*error;
            const missing=scenario==='heteroscedastic_mcar'
              ? next()<0.12
              : scenario==='mnar_left_censor'
                ? value<(7.3+j*0.1)
                : scenario==='mnar_group_asymmetric'
                  ? (treated?error>0.6:error< -0.6)
                  : false;
            const recorded=missing?'':value.toFixed(9);
            row.push(recorded);
            allObservations.push({
              case_id:caseId,scenario,n,replicate:rep,
              effect_planted:planted?1:0,feature,subject_id:ids[i],
              condition:metadataRows[i].condition,batch:metadataRows[i].batch,
              age,value:recorded
            });
          }
          matrix.push(csv(row));
        }
        const input={
          files:{
            metadata:new File(['fixture'],'subjects.csv',{type:'text/csv'}),
            transcriptomics:null,proteomics:new File([matrix.join('\n')],'protein.csv',{type:'text/csv'}),
            metabolomics:null
          },
          metadataRows,
          columnMapping:{
            subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',
            omic:'omic',condition:'condition',timepoint:'timepoint',batch:'batch',
            age:'age'
          },
          protocol:{
            organism:'human',objective:'groups',longitudinal:false,
            designType:'independent',batchKnown:'yes',groupCount:'2',
            studySetting:'synthetic_validation',covariateColumns:['age']
          },
          dataTypes:{proteomics:'log_intensity'},
          resolveIdentifiers:false,useReactome:false
        };
        const result=await runDeterministicAnalysis(input);
        const fitted=result.layers.proteomics.rows;
        // Retain feature attrition as a visible QC outcome; compare exactly the
        // actually tested hypotheses rather than inventing estimates for
        // features that did not pass the browser's own predeclared QC.
        assert.ok(fitted.length > 0 && fitted.length <= SPEC.features,
          'No estimable hypotheses remain for '+caseId);
        for(const p of fitted) {
          assert.ok(['effect','standardError','pValue','qValue','ciLow','ciHigh']
            .every(key=>Number.isFinite(p[key])),
            'Non-finite browser result: '+caseId+'/'+p.feature);
          estimates.push(csv([
            caseId,scenario,n,rep,planted?1:0,p.feature,
            p.effect,p.standardError,p.pValue,p.qValue,p.ciLow,p.ciHigh
          ]));
        }
        for(const observation of allObservations)data.push(csv(Object.values(observation)));
        cases.push({case_id:caseId,scenario,n,replicate:rep,planted,
          hypothesesSupplied:SPEC.features,hypothesesRetained:fitted.length,
          filteredFeatures:SPEC.features-fitted.length});
      }
    }
  }
  await Promise.all([
    fs.writeFile(path.join(directory,'observations.csv'),data.join('\n')+'\n'),
    fs.writeFile(path.join(directory,'browser-estimates.csv'),estimates.join('\n')+'\n'),
    fs.writeFile(path.join(directory,'protocol.json'),JSON.stringify({
      ...SPEC,cases:cases.length,
      dataGeneration:'independent seeded deterministic JS simulator',
      rFit:'R stats::lm + independent matrix HC3 variance estimator and Student-t inference',
      inferentialScope:'2 independent groups, featurewise adjusted Gaussian OLS; not publication approval',
      mnarPolicy:'MNAR is unidentifiable from observed data alone, independent R agreement does NOT rescue it',
      importantLimit:'Feature-wise model reference comparison conditions on browser-retained feature family; this does not validate upstream QC selection.',
    },null,2)+'\n')
  ]);
  console.log('Generated '+cases.length+' paired simulated datasets and '+(estimates.length-1)+' browser estimates');
}

function wilson(success,n,z=1.959963984540054){
  if(!n)return {lower:null,upper:null};
  const p=success/n,z2=z*z,den=1+z2/n;
  const mid=(p+z2/(2*n))/den,half=z*Math.sqrt((p*(1-p)+z2/(4*n))/n)/den;
  return {lower:Math.max(0,mid-half),upper:Math.min(1,mid+half)};
}
async function compare(directory) {
  const [bText,rText,pText]=await Promise.all([
    fs.readFile(path.join(directory,'browser-estimates.csv'),'utf8'),
    fs.readFile(path.join(directory,'r-estimates.csv'),'utf8'),
    fs.readFile(path.join(directory,'protocol.json'),'utf8')
  ]);
  const protocol=JSON.parse(pText),browser=readCsv(bText),reference=readCsv(rText);
  assert.equal(browser.length,reference.length,'Exact paired feature count is mandatory');
  const id=row=>row.case_id+'|'+row.feature;
  const rById=new Map(reference.map(row=>[id(row),row]));
  assert.equal(rById.size,reference.length,'Duplicated R feature rows');
  const metrics=Object.keys(SPEC.rNumericalTolerance),maximum=Object.fromEntries(metrics.map(x=>[x,0]));
  const scenarioMap=new Map(),failures=[];
  for(const row of browser) {
    const ref=rById.get(id(row));
    assert.ok(ref,'Missing independent R model: '+id(row));
    const key=[row.scenario,row.n,row.effect_planted].join('|');
    if(!scenarioMap.has(key))scenarioMap.set(key,[]);
    const differences={};
    for(const metric of metrics){
      const x=Number(row[metric]),y=Number(ref[metric]);
      if(!Number.isFinite(x)||!Number.isFinite(y)) {
        failures.push({feature:id(row),metric,message:'non-finite numerical output'});continue;
      }
      const delta=Math.abs(x-y);
      maximum[metric]=Math.max(maximum[metric],delta);
      differences[metric]=delta;
      if(delta>SPEC.rNumericalTolerance[metric])
        failures.push({feature:id(row),metric,absoluteDifference:delta,
          tolerance:SPEC.rNumericalTolerance[metric]});
    }
    scenarioMap.get(key).push({
      ...row,reference:ref,differences,
      trueEffect:row.effect_planted==='1' && Number(row.feature.replace('PROT_',''))<=SPEC.signalFeatures
        ? SPEC.plantedEffect:0
    });
  }
  const summaries=[];
  for(const [key,rows] of scenarioMap){
    const [scenario,nStr,plantedStr]=key.split('|');
    const isPlanted=plantedStr==='1',n=Number(nStr);
    const nullRows=rows.filter(x=>!isPlanted || x.trueEffect===0);
    const rateCount=nullRows.filter(x=>Number(x.pValue)<SPEC.nominalAlpha).length;
    const families=new Map();
    for(const row of rows){
      const id=row.case_id;
      if(!families.has(id))families.set(id,[]);
      families.get(id).push(row);
    }
    const familyFailures=[...families.values()].filter(values=>
      values.some(v=>v.trueEffect===0&&Number(v.qValue)<SPEC.nominalAlpha)).length;
    let totalDiscoveries=0,totalFalseDiscoveries=0;
    let sumFdp=0;
    for(const values of families.values()) {
      const called=values.filter(x=>Number(x.qValue)<=SPEC.nominalAlpha);
      const errors=called.filter(x=>x.trueEffect===0).length;
      totalDiscoveries+=called.length;
      totalFalseDiscoveries+=errors;
      sumFdp+=called.length?errors/called.length:0;
    }
    const coverage=rows.filter(x=>{
      const low=Number(x.ciLow),high=Number(x.ciHigh);
      return low<=x.trueEffect && x.trueEffect<=high;
    }).length;
    const signals=rows.filter(x=>x.trueEffect>0);
    const hits=signals.filter(x=>Number(x.qValue)<=0.10&&Number(x.effect)>0).length;
    const signalBias=signals.length?signals.reduce((s,x)=>s+Number(x.effect)-x.trueEffect,0)/signals.length:null;
    const mnar=scenario.startsWith('mnar');
    summaries.push({
      scenario,n,planted:isPlanted,datasets:families.size,hypotheses:rows.length,
      trueNullTests:nullRows.length,rawFalsePositives:rateCount,
      rawTypeIRate:nullRows.length?rateCount/nullRows.length:null,
      rawRateWilson95:wilson(rateCount,nullRows.length),
      familiesWithAnyFalseQ005:familyFailures,
      falseDiscoveryFamilyFraction:familyFailures/families.size,
      totalDiscoveries,totalFalseDiscoveries,
      meanRealisedFDP:sumFdp/families.size,
      intervalCoverage:coverage/rows.length,
      intervalCoverageWilson95:wilson(coverage,rows.length),
      knownSignals:signals.length,knownSignalsRecovered: hits,
      positiveSignalRecall:signals.length?hits/signals.length:null,
      positiveEffectBias:signalBias,
      evidenceStatus:mnar?'blocked_MNAR_nonidentifiable':'calibration_evidence_only',
      interpretation:mnar
        ? 'MNAR probabilities depend on unobserved abundances; reference agreement cannot identify unbiased biological effects. Sensitivity models needed.'
        : 'Limited empirical simulation: not an equivalence proof for 95% coverage, 5% size or study-wide FDR.'
    });
  }
  const report={
    version:SPEC.version,nominalAlpha:SPEC.nominalAlpha,
    simulations:protocol.cases,pairwiseModels:browser.length,
    featureAttrition:{
      requested:protocol.cases*SPEC.features,
      retained:browser.length,
      filtered:protocol.cases*SPEC.features-browser.length,
      interpretation:'The R comparison conditions on browser-retained features; prefiltering and post-selection calibration are not independently validated.'
    },
    comparison:'Separate JS application and R base::lm + independently implemented HC3',
    numericalTolerances:SPEC.rNumericalTolerance,
    maximumAbsoluteDeviation:maximum,verificationFailures:failures.length,
    failureExamples:failures.slice(0,12),
    numericalAgreement:failures.length===0?'pass':'fail',
    scientificCertification:false,
    observedDataMnarIdentifiable:false,
    summary:summaries,
    limitations:[
      'R checks the model and BH procedure on the same QC-retained features, not independent QC eligibility.',
      'Empirical size, FDR and coverage estimates carry Monte Carlo uncertainty; correlated feature tests are not independent.',
      'No validation of nonparametric/survival/longitudinal differential inference, MS chemistry, or patient cohorts.',
      'MNAR mechanisms are not learnable from these observations; default conclusion is blocked even if numerical models agree.'
    ]
  };
  await fs.writeFile(path.join(directory,'methodology-report.json'),JSON.stringify(report,null,2)+'\n');
  console.log(JSON.stringify({
    simulations:report.simulations,models:report.pairwiseModels,
    numericalAgreement:report.numericalAgreement,
    maxDeviation:maximum,
    scenarios:summaries.map(x=>({scenario:x.scenario,n:x.n,planted:x.planted,
      rawRate:x.rawTypeIRate,coverage:x.intervalCoverage,
      evidenceStatus:x.evidenceStatus}))
  },null,2));
  assert.equal(failures.length,0,'Independent R/JS OLS HC3 mismatch; see methodology-report.json');
}
const mode=process.argv[2],dir=process.argv[3]||'tmp/multiomics/methodology';
if(mode==='--generate')await generate(dir);
else if(mode==='--compare')await compare(dir);
else throw new Error('Usage: node scripts/validate_multiomics_methodology.mjs --generate|--compare <directory>');
