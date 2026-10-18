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

function objectiveFamily(protocol){
  const objective=protocol?.objective;
  if(objective==='groups')return 'groups';
  if(objective==='time')return 'longitudinal';
  if(objective==='outcome' && protocol?.outcomeType==='survival')return 'survival';
  return null;
}

function candidateMethods(layer,datum,backend,protocol){
  const family=objectiveFamily(protocol);
  if(family==='groups'){
    const rawRna=layer==='transcriptomics' &&
      datum?.qc?.preprocessingAudit?.declaredValueType==='raw_counts';
    const keys=rawRna
      ? [layer+'_differential_deseq2',layer+'_differential_voom']
      : [layer+'_differential'];
    return keys.map(key=>({key,method:backend.methods?.[key],family}));
  }
  if(family==='longitudinal'){
    const key=layer+'_longitudinal';
    return [{key,method:backend.methods?.[key],family}];
  }
  if(family==='survival'){
    const key=layer+'_survival';
    return [{key,method:backend.methods?.[key],family}];
  }
  return [];
}

function choose(layer,datum,backend,protocol){
  return candidateMethods(layer,datum,backend,protocol).find(x=>ready(x.method))||null;
}

function columnsFor(family,cols){
  const idKey=cols.includes('feature')?'feature':cols.includes('feature_id')?'feature_id':null;
  if(family==='groups'){
    return {
      idKey,
      effectKey:cols.includes('log2FoldChange')?'log2FoldChange':cols.includes('logFC')?'logFC':null,
      pKey:cols.includes('pvalue')?'pvalue':cols.includes('P.Value')?'P.Value':
        cols.includes('PValue')?'PValue':null,
      qKey:cols.includes('padj')?'padj':cols.includes('adj.P.Val')?'adj.P.Val':
        cols.includes('FDR')?'FDR':null
    };
  }
  if(family==='longitudinal'){
    return {
      idKey,
      effectKey:cols.includes('estimate')?'estimate':null,
      pKey:cols.includes('p_value')?'p_value':null,
      qKey:cols.includes('q_bh')?'q_bh':null
    };
  }
  if(family==='survival'){
    return {
      idKey,
      effectKey:cols.includes('effect')?'effect':null,
      pKey:cols.includes('pValue')?'pValue':null,
      qKey:cols.includes('qValue')?'qValue':null
    };
  }
  return {idKey:null,effectKey:null,pKey:null,qKey:null};
}

export function projectNativeRNumbers(result){
  if(!result || result.referenceBackend?.status!=='ok')return result;
  const family=objectiveFamily(result.protocol);
  if(!family)return result;

  const updated={
    ...result,
    layers:{...result.layers},
    numericSourceByLayer:{},
    browserExploratoryLayers:{}
  };
  let nativeCount=0;

  for(const layer of OMICS){
    const original=result.layers?.[layer];
    if(!original)continue;
    const selected=choose(layer,original,result.referenceBackend,result.protocol);
    if(!selected){
      updated.numericSourceByLayer[layer]='browser_exploratory';
      continue;
    }

    const parsed=parseDelimited(selected.method.referenceTableCsv);
    const cols=parsed?.headers||[];
    const {idKey,effectKey,pKey,qKey}=columnsFor(selected.family,cols);
    if(!idKey||!effectKey||!pKey||!qKey||!parsed.rows.length)
      throw new Error('Reference R table cannot be displayed faithfully: columns missing: '+layer);
    if(Number.isInteger(selected.method.referenceTableRows) &&
      parsed.rows.length!==selected.method.referenceTableRows)
      throw new Error('Incomplete R reference CSV for '+layer);

    const groupLog2=selected.family==='groups' &&
      (selected.key.endsWith('_differential_deseq2')||
       selected.key.endsWith('_differential_voom'));
    const scale=selected.family==='survival'?'log_hazard':
      selected.family==='longitudinal'?'group_by_time_slope':
      groupLog2?'log2':original.effectScale;

    const seen=new Set();
    const rows=parsed.rows.map(row=>{
      const feature=String(row[idKey]||'');
      if(!feature||seen.has(feature))throw Error('Duplicated or missing R feature '+feature);
      seen.add(feature);
      const effect=numeric(row[effectKey]);
      const pValue=numeric(row[pKey]);
      const qValue=numeric(row[qKey]);
      const projected={
        feature,effect,effectScale:scale,
        foldRatio:scale==='log2'&&effect!==null?Math.pow(2,effect):null,
        pValue,qValue,
        model:selected.method.method||selected.key,
        originalReferenceColumns:row
      };
      if(selected.family==='survival'){
        projected.hazardRatio=numeric(row.hazardRatio);
        projected.proportionalHazardsP=numeric(row.proportionalHazardsP);
        projected.globalPhP=numeric(row.globalPhP);
        projected.nSubjects=numeric(row.nSubjects);
        projected.nEvents=numeric(row.nEvents);
        projected.status=String(row.status||'');
      } else if(selected.family==='longitudinal'){
        projected.term=String(row.term||'');
        projected.standardError=numeric(row.std_error);
        projected.degreesOfFreedom=numeric(row.df);
        projected.statistic=numeric(row.statistic);
        projected.nObservations=numeric(row.n_observations);
        projected.nSubjects=numeric(row.n_subjects);
        projected.randomEffectStructure=String(row.random_effect_structure||'');
        projected.status=String(row.status||'');
      }
      return projected;
    });

    const ranked=rows.filter(r=>r.effect!==null)
      .sort((a,b)=>(a.qValue??1)-(b.qValue??1) ||
        Math.abs(b.effect)-Math.abs(a.effect) || a.feature.localeCompare(b.feature));
    const significant=ranked.filter(r=>r.qValue!==null&&r.qValue<=0.10);

    updated.browserExploratoryLayers[layer]=original;
    updated.layers[layer]={
      ...original,
      rows,
      selected:(significant.length>=10?significant:ranked).slice(0,50),
      mode:'reference-r',
      inferenceMethod:selected.method.method||selected.key,
      inferencePolicy:'Complete output from native R method; original browser exploratory fit remains separate.',
      selectionRule:'R native full table; q <= 0.10 is an exploratory display/ranking threshold, not a universal confirmatory cutoff.',
      effectScale:scale,
      contrast:original.contrast,
      steps:[...(original.steps||[]),'Native R results from '+selected.key],
      primaryReferenceMethodKey:selected.key,
      primaryReferenceTableRows:rows.length,
      numericalProvenance:'exact_native_r_reference_table',
      nativeObjectiveFamily:selected.family
    };
    updated.numericSourceByLayer[layer]='native_r_complete_table';
    nativeCount++;
  }
  updated.displayedReferenceMethodCount=nativeCount;
  return updated;
}
