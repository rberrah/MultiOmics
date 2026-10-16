// Additional, predeclared simulation envelopes: adjusted heteroscedastic
// independent groups, technical-batch confounding, and repeated measures.
// These are diagnostic smoke tests, not comprehensive model validation.
import assert from 'node:assert/strict';
import { File } from 'node:buffer';
import { runDeterministicAnalysis } from '../src/lib/multiomics/deterministic.js';

const commonMapping={subject_id:'subject_id',sample_id:'sample_id',assay_id:'assay_id',
  omic:'omic',condition:'condition',timepoint:'timepoint',batch:'batch',age:'age'};
const valueTypes={transcriptomics:'raw_counts',proteomics:'log_intensity',metabolomics:'peak_area'};
const seeded=(seed)=>{let state=seed>>>0;return ()=>((state=(Math.imul(1664525,state)+1013904223)>>>0)+0.5)/4294967296};
const normal=(rng)=>Math.sqrt(-2*Math.log(Math.max(1e-15,rng()))) * Math.cos(2*Math.PI*rng());
const matrixFile=(ids,values)=>new File([[
  ['feature_id',...ids].join(','),
  ...values.map((valuesForFeature,j)=>['PROT_'+j,...valuesForFeature.map(v=>v.toFixed(8))].join(','))
].join('\n')],'proteomics.csv',{type:'text/csv'});

function independentFixture(seed,{planted=false,confounded=false}={}) {
  const rand=seeded(seed),subjects=40,features=16;
  const ids=Array.from({length:subjects},(_,i)=>'S'+i);
  const metadataRows=ids.map((id,i)=>{
    const condition=i<20?'control':'treatment';
    const batch=confounded?(i<20?'B1':'B2'):(i%2?'B2':'B1');
    return {subject_id:id,sample_id:id,assay_id:id,omic:'proteomics',
      condition,timepoint:'T0',batch,age:String(30+(i*17)%31)};
  });
  const rows=Array.from({length:features},(_,f)=>ids.map((_,i)=>{
    const m=metadataRows[i];
    const batch=m.batch==='B2'?1.7:0;
    const age=(Number(m.age)-45)*0.035;
    const variance=i<20?0.7:1.65;
    const effect=planted&&f<3&&i>=20?2.3:0;
    return 8+batch+age+effect+variance*normal(rand);
  }));
  return {files:{metadata:new File(['placeholder'],'samples.csv'),proteomics:matrixFile(ids,rows)},
    metadataRows,protocol:{objective:'groups',organism:'human',designType:'independent',
      longitudinal:false,covariateColumns:['age'],groupCount:'2'}};
}

function longitudinalFixture(seed,{planted=false}={}) {
  const rand=seeded(seed), subjects=24,features=8;
  const ids=[],metadataRows=[],values=Array.from({length:features},()=>[]);
  const randomIntercept=Array.from({length:subjects},()=>0.6*normal(rand));
  for(let i=0;i<subjects;i++) for(let visit=0;visit<2;visit++){
    const id='SUB'+i+'_T'+visit;
    const cond=i<12?'control':'treatment';
    ids.push(id);
    metadataRows.push({subject_id:'SUB'+i,sample_id:id,assay_id:id,omic:'proteomics',
      condition:cond,timepoint:visit?'T12':'T0',batch:'B'+(i%2),age:String(25+i)});
    for(let f=0;f<features;f++) {
      const effect=planted&&f<2&&i>=12&&visit===1?1.8:0;
      values[f].push(9+randomIntercept[i]+0.65*normal(rand)+effect);
    }
  }
  return {files:{metadata:new File(['placeholder'],'samples.csv'),proteomics:matrixFile(ids,values)},
    metadataRows,protocol:{objective:'time',organism:'human',designType:'repeated',
      longitudinal:true,covariateColumns:[],groupCount:'2'}};
}

const analyze=(fixture)=>runDeterministicAnalysis({
  files:fixture.files,metadataRows:fixture.metadataRows,columnMapping:commonMapping,
  protocol:fixture.protocol,dataTypes:valueTypes,resolveIdentifiers:false,useReactome:false
});

const threshold=0.05,qThreshold=0.10;
async function evaluateScenario(label,fixture,nullRuns,signalRuns,signals,pLimit,familyLimit) {
  const descriptive=label.includes('longitudinal');
  let pTotal=0,pSmall=0,families=0,hits=0,truthTotal=0,modelsEvaluated=0;
  for(let k=0;k<nullRuns;k++) {
    const result=await analyze(fixture(10001+101*k));
    const rows=result.layers.proteomics.rows;
    assert(rows.length>=(descriptive?7:14),
      label+': too few null models ('+rows.length+') error='+JSON.stringify(result.layers.proteomics.error)
      +' modes='+JSON.stringify({mode:result.layers.proteomics.mode,
        nSubjects:result.metadataSummary?.subjects,groupSizes:result.layers.proteomics.groupSizes,
        selectedFeatures:result.layers.proteomics.selected?.length}));
    modelsEvaluated+=rows.length;
    let any=false;
    for(const row of rows) {
      if(descriptive) {
        // Uncalibrated GLS previously had 10% p<.05 under null (12/120).
        // Refuse all inferential outputs, not just nominally significant ones.
        assert.equal(row.pValue,null,label+': unsafe longitudinal p-value');
        assert.equal(row.qValue,null,label+': unsafe longitudinal BH q');
        assert.equal(row.ciLow,null,label+': unsafe longitudinal confidence interval');
        assert.equal(row.ciHigh,null,label+': unsafe longitudinal confidence interval');
        assert.equal(row.inferentialStatus,'uncalibrated_reference_R_required');
      } else {
        if (!Number.isFinite(row.pValue)||!Number.isFinite(row.qValue)) throw Error(label+': invalid p/q');
        pTotal++;
        if(row.pValue<threshold)pSmall++;
        if(row.qValue<=qThreshold)any=true;
      }
    }
    if(any)families++;
  }
  for(let k=0;k<signalRuns;k++){
    const result=await analyze(fixture(20001+41*k,{planted:true}));
    const rows=result.layers.proteomics.rows;
    for(const row of rows){
      if(!signals.includes(row.feature))continue;
      truthTotal++;
      if(descriptive ? (row.effect>0.8 && row.pValue===null && row.qValue===null)
        : (row.qValue<=qThreshold && row.effect>0)) hits++;
    }
  }
  const nullRate=descriptive?null:pSmall/pTotal;
  const familyRate=descriptive?null:families/nullRuns;
  const power=hits/truthTotal;
  if(!descriptive) {
    assert(nullRate<=pLimit,`${label}: gross type-I inflation ${nullRate}`);
    assert(familyRate<=familyLimit,`${label}: gross BH family inflation ${familyRate}`);
  }
  assert(power>=0.55,`${label}: failed planted positive controls ${power}`);
  return {label,
    inferentialStatus:descriptive?'browser_p_q_ci_suppressed_pending_R_validation':'limited_null_calibration',
    nullRuns,modelsEvaluated,pTotal:descriptive?null:pTotal,
    pSmall:descriptive?null:pSmall,nullRate,families:descriptive?null:families,familyRate,
    signalRuns,truthTotal,hits,power};
}

const adjusted=await evaluateScenario('heteroscedastic_balanced_batch',
  independentFixture,24,8,['PROT_0','PROT_1','PROT_2'],0.15,0.36);
const repeated=await evaluateScenario('random_intercept_longitudinal',
  longitudinalFixture,15,6,['PROT_0','PROT_1'],0.20,0.40);
await assert.rejects(()=>analyze(independentFixture(9784,{confounded:true})),
  /confound|identif|batch/i,
  'perfectly confounded assay batch and condition must not produce p-values');
console.log(JSON.stringify({envelope:'multiomics_stress_v1',seeded:true,
  limitations:'Only two simulated gaussian designs. Not RNA negative-binomial, MNAR, Cox, missing-not-at-random, or real cohort external validation.',
  scenarios:[adjusted,repeated],fullyConfoundedBatch:'refused'},null,2));
