// @ts-nocheck
// Publication-facing numeric tables are projected FROM the actual R reference
// backend's COMPLETE CSV, never inferred from browser statistics.
import {parseDelimited} from './deterministic.js';
const OMICS=['transcriptomics','proteomics','metabolomics'];
function numeric(value){
  if(value==null||String(value).trim()===''||String(value).toUpperCase()==='NA')return null;
  const n=Number(value);
  return Number.isFinite(n)?n:null;
}
function ready(method){
  return method?.status==='ok' && method?.referenceTableComplete===true &&
    typeof method?.referenceTableCsv==='string' && method.referenceTableCsv.trim().length>0;
}
function choose(layer,datum,backend){
  const keys=layer==='transcriptomics' && datum?.qc?.preprocessingAudit?.declaredValueType==='raw_counts'
    ? [layer+'_differential_deseq2',layer+'_differential_voom']
    : [layer+'_differential'];
  return keys.map(key=>({key,method:backend.methods?.[key]})).find(x=>ready(x.method))||null;
}
export function projectNativeRNumbers(result){
  if(!result || result.referenceBackend?.status!=='ok' ||
    result.protocol?.objective!=='groups')return result;
  const updated={...result,layers:{...result.layers},numericSourceByLayer:{},browserExploratoryLayers:{}};
  let nativeCount=0;
  for(const layer of OMICS){
    const original=result.layers?.[layer];
    if(!original)continue;
    const selected=choose(layer,original,result.referenceBackend);
    if(!selected){
      updated.numericSourceByLayer[layer]='browser_exploratory';
      continue;
    }
    const parsed=parseDelimited(selected.method.referenceTableCsv);
    const cols=parsed?.headers||[];
    const idKey=cols.includes('feature')?'feature':cols.includes('feature_id')?'feature_id':null;
    const effectKey=cols.includes('log2FoldChange')?'log2FoldChange':cols.includes('logFC')?'logFC':null;
    const pKey=cols.includes('pvalue')?'pvalue':cols.includes('P.Value')?'P.Value':
      cols.includes('PValue')?'PValue':null;
    const qKey=cols.includes('padj')?'padj':cols.includes('adj.P.Val')?'adj.P.Val':
      cols.includes('FDR')?'FDR':null;
    if(!idKey||!effectKey||!pKey||!qKey || !parsed.rows.length)
      throw new Error('Reference R table cannot be displayed faithfully: columns missing: '+layer);
    if(Number.isInteger(selected.method.referenceTableRows) &&
      parsed.rows.length!==selected.method.referenceTableRows)
      throw new Error('Incomplete R reference CSV for '+layer);
    const scale=selected.key.endsWith('_differential_deseq2')||
      selected.key.endsWith('_differential_voom')?'log2':original.effectScale;
    const seen=new Set();
    const rows=parsed.rows.map(row=>{
      const feature=String(row[idKey]||'');
      if(!feature||seen.has(feature))throw Error('Duplicated or missing R feature '+feature);
      seen.add(feature);
      const effect=numeric(row[effectKey]);
      const pValue=numeric(row[pKey]);
      const qValue=numeric(row[qKey]);
      return {
        feature,effect,effectScale:scale,
        foldRatio:scale==='log2'&&effect!==null?Math.pow(2,effect):null,
        pValue,qValue,
        model:selected.method.method||selected.key,
        originalReferenceColumns:row
      };
    });
    const ranked=rows.filter(r=>r.effect!==null)
      .sort((a,b)=>(a.qValue??1)-(b.qValue??1) ||
        Math.abs(b.effect)-Math.abs(a.effect) || a.feature.localeCompare(b.feature));
    const significant=ranked.filter(r=>r.qValue!==null&&r.qValue<=0.10);
    updated.browserExploratoryLayers[layer]=original;
    updated.layers[layer]={
      ...original,
      rows,selected:(significant.length>=10?significant:ranked).slice(0,50),
      mode:'reference-r',
      inferenceMethod:selected.method.method||selected.key,
      inferencePolicy:'Complete output from native R method; original browser exploratory fit remains separate.',
      selectionRule:'R native full table; q <= 0.10 exploratory ranking, otherwise effect ranking',
      effectScale:scale,
      contrast:original.contrast,
      steps:[...(original.steps||[]),'Native R results from '+selected.key],
      primaryReferenceMethodKey:selected.key,
      primaryReferenceTableRows:rows.length,
      numericalProvenance:'exact_native_r_reference_table'
    };
    updated.numericSourceByLayer[layer]='native_r_complete_table';
    nativeCount++;
  }
  updated.displayedReferenceMethodCount=nativeCount;
  return updated;
}
