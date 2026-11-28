// Synthetic batch-design teaching model. No real dataset or clinical claim.
export function batchExample(crossed=true,batchEffect=2,biologicalEffect=0){
  if(typeof crossed!=='boolean'||![batchEffect,biologicalEffect].every(Number.isFinite)||batchEffect<0||batchEffect>4||biologicalEffect<0||biologicalEffect>2)throw new RangeError('Teaching domain');
  const noise=[-.3,-.1,.1,.3],points=[];
  for(const batch of [0,1])for(const group of [0,1]){
    if(!crossed&&batch!==group)continue;
    for(const e of noise)points.push({batch,group,value:batch*batchEffect+group*biologicalEffect+e});
  }
  const mean=(/** @type {{value:number}[]} */ rows)=>rows.reduce((s,p)=>s+p.value,0)/rows.length;
  const difference=mean(points.filter(p=>p.group===1))-mean(points.filter(p=>p.group===0));
  const within=crossed?[0,1].map(batch=>mean(points.filter(p=>p.batch===batch&&p.group===1))-mean(points.filter(p=>p.batch===batch&&p.group===0))):null;
  return {points,difference,within};
}
