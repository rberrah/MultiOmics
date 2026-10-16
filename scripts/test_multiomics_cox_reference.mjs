// Strictly paired numerical reference for Cox proportional hazards coefficients.
// The test covers Breslow tied event times, continuous event times and censoring.
// Neither the test nor R equivalence establishes the PH assumption or validity
// with time-dependent exposures / informative censoring.
import assert from 'node:assert/strict';
import fs from 'node:fs/promises';
import path from 'node:path';
import { coxRegression } from '../src/lib/multiomics/deterministic-core.js';

function uniform(seed) {
  let x=seed>>>0;
  return ()=>((x=(Math.imul(x,1664525)+1013904223)>>>0)+0.5)/4294967296;
}
function gaussian(u) {
  return Math.sqrt(-2*Math.log(Math.max(u(),1e-12)))*Math.cos(2*Math.PI*u());
}
function parse(text) {
  const rows=text.trim().split(/\r?\n/);
  const head=rows.shift().split(',');
  return rows.map(row=>Object.fromEntries(row.split(',').map((value,i)=>[head[i],value])));
}
const dir=process.argv[3]||'tmp/multiomics/cox-methodology';
const mode=process.argv[2];
await fs.mkdir(dir,{recursive:true});
if(mode==='--generate') {
  const observations=[['case_id','sample_id','time','event','feature','batch','true_beta'].join(',')];
  const fits=[['case_id','n','events','ties','beta','se','p','tieMethod'].join(',')];
  for(const n of [48,80])for(const tied of [false,true])for(const betaTruth of [0,0.55]) {
    for(let rep=0;rep<10;rep++) {
      const id=[n,tied?'tied':'continuous',betaTruth?'signal':'null',rep].join('_');
      const u=uniform(20261009+n*1000+(tied?100:0)+(betaTruth?500:0)+rep*19);
      const design=[],times=[],events=[];
      for(let j=0;j<n;j++){
        const x=gaussian(u),batch=j%2;
        const hazard=0.065*Math.exp(betaTruth*x+0.38*batch);
        const eventTime=-Math.log(Math.max(u(),1e-12))/hazard;
        const censorTime=-Math.log(Math.max(u(),1e-12))/0.024;
        const event=eventTime<=censorTime?1:0;
        const rawTime=Math.min(eventTime,censorTime);
        const time=tied?Math.max(1,Math.ceil(rawTime)):rawTime;
        design.push([x,batch]);times.push(time);events.push(event);
        observations.push([id,'P'+(j+1),time,event,x,batch,betaTruth].join(','));
      }
      const fit=coxRegression(design,times,events,0);
      assert.ok(fit&&fit.converged,'Cox must converge on prespecified simulation '+id);
      assert.equal(fit.ties,'breslow');
      assert.ok(Number.isFinite(fit.se)&&fit.se>0);
      fits.push([id,n,events.reduce((a,b)=>a+b,0),tied?1:0,fit.beta[0],fit.se,fit.pValue,fit.ties].join(','));
    }
  }
  await fs.writeFile(path.join(dir,'survival-observations.csv'),observations.join('\n')+'\n');
  await fs.writeFile(path.join(dir,'browser-cox.csv'),fits.join('\n')+'\n');
  await fs.writeFile(path.join(dir,'protocol.json'),JSON.stringify({
    seed:20261009,n:[48,80],eventTime:['continuous','rounded_up_integer'],
    scenarios:['null_coefficient','known_log_hazard_ratio_0.55'],
    datasets:80,model:'Cox unstratified 2 covariate Breslow partial likelihood',
    censoring:'independent exponential',
    coefficientsRef:'R survival::coxph(ties=breslow)',
    nominalAlpha:0.05,
    numericalTolerance:{beta:0.001,standardError:0.001,pValue:0.001},
    scientificCertification:false,
    limitations:['Only independent right-censoring and time-constant PH simulation','No claim of PH, time-dependent exposures or transportability from numerical agreement']
  },null,2)+'\n');
  console.log('Generated 80 survival datasets, Breslow Cox fitted for each.');
} else if(mode==='--compare') {
  const [js,rs]=await Promise.all([
    fs.readFile(path.join(dir,'browser-cox.csv'),'utf8'),
    fs.readFile(path.join(dir,'r-cox.csv'),'utf8')
  ]);
  const a=parse(js),b=parse(rs),map=new Map(b.map(x=>[x.case_id,x]));
  assert.equal(a.length,80);assert.equal(b.length,80);assert.equal(map.size,80);
  const maxima={beta:0,standardError:0,pValue:0},disagreements=[],rows=[];
  for(const item of a) {
    const ref=map.get(item.case_id);
    assert.ok(ref,'R missing fit '+item.case_id);
    const deltas={
      beta:Math.abs(Number(item.beta)-Number(ref.beta)),
      standardError:Math.abs(Number(item.se)-Number(ref.se)),
      pValue:Math.abs(Number(item.p)-Number(ref.p))
    };
    for(const [key,delta] of Object.entries(deltas)) {
      assert.ok(Number.isFinite(delta));
      maxima[key]=Math.max(maxima[key],delta);
      if(delta>0.001)disagreements.push({id:item.case_id,key,delta});
    }
    rows.push({
      case_id:item.case_id,n:Number(item.n),
      tied:item.ties==='1',
      events:Number(item.events),
      rPHCheckP:Number(ref.ph_p),
      referenceConverged:ref.converged==='TRUE',
      ...deltas
    });
  }
  const report={
    version:'independent_cox_breslow_v1',
    pairs:a.length,
    maxAbsoluteDeviation:maxima,
    tolerance:0.001,
    mismatchedModels:disagreements.length,
    exampleDisagreements:disagreements.slice(0,10),
    result:disagreements.length?'fail':'pass',
    phDiagnostic:'R survival::cox.zph feature-level test for each simulated study (descriptive only, not automated proof of PH)',
    scientificCertification:false,
    rows
  };
  await fs.writeFile(path.join(dir,'cox-report.json'),JSON.stringify(report,null,2)+'\n');
  console.log(JSON.stringify({pairs:a.length,maxima,result:report.result}));
  assert.equal(disagreements.length,0,'Breslow Cox numerical R disagreement; inspect cox-report.json');
} else throw Error('Usage: node scripts/test_multiomics_cox_reference.mjs --generate|--compare <dir>');
