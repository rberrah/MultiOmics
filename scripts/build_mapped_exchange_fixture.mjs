#!/usr/bin/env node
import fs from 'node:fs/promises';
import path from 'node:path';
import {webcrypto} from 'node:crypto';
import {buildDownstreamExchange} from '../src/lib/multiomics/downstream-exchange.js';
if(!globalThis.crypto)globalThis.crypto=webcrypto;
const out=process.argv[2];
if(!out)throw Error('Usage: node scripts/build_mapped_exchange_fixture.mjs OUTPUT_DIR');
const metadata=[
  'patient_id\tspecimen\trun\tmodality\ttissue_name\tvisit\tvisit_unit\tanchor',
  'P10\tS10\tA10\trna\tblood\t0\th\tfirst_dose',
  'P11\tS11\tA11\trna\tblood\t12\th\tfirst_dose'
].join('\n')+'\n';
const project={
  inputs:{
    metadata:{content:metadata,sha256:'d'.repeat(64)},
    transcriptomics:{content:'assay;ENSG_A;ENSG_B\nA10;10;\nA11;12;15\n'}
  },
  settings:{
    organism:'human',unitType:'participant',
    transcriptomicsValues:'normalized',transcriptomicsIdType:'ensembl_gene',
    columnMapping:{
      subject_id:'patient_id',sample_id:'specimen',assay_id:'run',omic:'modality',
      tissue:'tissue_name',collection_time:'visit',collection_time_unit:'visit_unit',
      time_origin:'anchor'
    }
  }
};
const exchange=await buildDownstreamExchange(project);
if(!exchange.summary.structural_exchange_valid||exchange.summary.pk_ready||
   !exchange.summary.pk_temporal_numeric_compatible||
   !exchange.summary.explicit_time_origin_complete)
  throw Error('Mapped fixture readiness boundary changed unexpectedly.');
await fs.rm(out,{recursive:true,force:true});
for(const entry of exchange.entries){
  const relative=entry.name.replace(/^exchange\//,'');
  const dest=path.join(out,relative);
  await fs.mkdir(path.dirname(dest),{recursive:true});
  await fs.writeFile(dest,entry.data,'utf8');
}
console.log(JSON.stringify(exchange.summary));
