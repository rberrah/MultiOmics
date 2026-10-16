import assert from 'node:assert/strict';
import {projectNativeRNumbers} from '../src/lib/multiomics/r-primary-results.js';
const csv=['"feature","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"',
  '"ENSG_A",50,1.5,0.2,7.5,0.00001,0.00002',
  '"ENSG_B",20,-0.8,0.3,-2.67,0.012,0.04',
  '"ENSG_C",0,NA,NA,NA,NA,NA'].join('\n')+'\n';
const original={mode:'browser-groups',rows:[{feature:'ENSG_A',effect:999,qValue:0.9}],
  selected:[{feature:'ENSG_A'}],groupSizes:[4,4],contrast:'treated vs control',
  steps:['Browser exploratory'],effectScale:'log2',
  qc:{preprocessingAudit:{declaredValueType:'raw_counts'}}};
const base={protocol:{objective:'groups'},layers:{transcriptomics:original},
  referenceBackend:{status:'ok',methods:{transcriptomics_differential_deseq2:{
    status:'ok',method:'DESeq2',referenceTableComplete:true,referenceTableRows:3,
    referenceTableCsv:csv}}}};
const projected=projectNativeRNumbers(base);
assert.notEqual(projected,base);
assert.equal(projected.numericSourceByLayer.transcriptomics,'native_r_complete_table');
assert.equal(projected.displayedReferenceMethodCount,1);
assert.equal(projected.browserExploratoryLayers.transcriptomics.rows[0].effect,999);
assert.deepEqual(projected.layers.transcriptomics.rows.map(r=>r.feature),
 ['ENSG_A','ENSG_B','ENSG_C']);
assert.equal(projected.layers.transcriptomics.rows[0].effect,1.5);
assert.equal(projected.layers.transcriptomics.rows[0].foldRatio,Math.pow(2,1.5));
assert.equal(projected.layers.transcriptomics.rows[1].qValue,0.04);
assert.equal(projected.layers.transcriptomics.rows[2].qValue,null);
assert.equal(projected.layers.transcriptomics.rows[2].effect,null);
assert.equal(projected.layers.transcriptomics.rows[0].originalReferenceColumns.lfcSE,'0.2');
const incomplete=structuredClone(base);
incomplete.referenceBackend.methods.transcriptomics_differential_deseq2.referenceTableComplete=false;
const fallback=projectNativeRNumbers(incomplete);
assert.equal(fallback.numericSourceByLayer.transcriptomics,'browser_exploratory');
assert.equal(fallback.layers.transcriptomics.rows[0].effect,999);
const wrongRows=structuredClone(base);
wrongRows.referenceBackend.methods.transcriptomics_differential_deseq2.referenceTableRows=4;
assert.throws(()=>projectNativeRNumbers(wrongRows),/Incomplete R reference CSV/);
const repeated=structuredClone(base);
repeated.referenceBackend.methods.transcriptomics_differential_deseq2.referenceTableCsv =
  csv.replace('"ENSG_B"','"ENSG_A"');
assert.throws(()=>projectNativeRNumbers(repeated),/Duplicated or missing/);
console.log('PASS: site primary numeric values are those of complete R CSV; missingness and fallback explicit');
