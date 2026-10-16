// @ts-nocheck
// TRUE end-to-end research workflow: researcher-supplied CSV -> Svelte UI ->
// consented local Plumber R /run -> actual DESeq2 fit -> ZIP containing
// COMPLETE R rows -> independent Node/Python verification. No mocked /run,
// no fabricated statistics, no outside network, and no clinical claims.
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';
import crypto from 'node:crypto';
import {spawnSync} from 'node:child_process';
import os from 'node:os';
import path from 'node:path';
import {parseReferenceCsv} from '../../src/lib/multiomics/research-reference-csv.js';

const GENE_COUNT=172;
const SUBJECT_COUNT=28;
function csvInputs(){
  const subjects=Array.from({length:SUBJECT_COUNT},(_,i)=>'RESEARCH_S'+String(i+1).padStart(3,'0'));
  const meta=['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,sample_type'];
  for(let i=0;i<subjects.length;i++){
    const id=subjects[i];
    meta.push([id,id,id,'transcriptomics',i<14?'control':'treated','T0',
      i%2===0?'B1':'B2','1','biological'].join(','));
  }
  const matrix=['feature_id,'+subjects.join(',')];
  for(let gene=1;gene<=GENE_COUNT;gene++){
    const vals=subjects.map((_,i)=>{
      const biological=gene<=25&&i>=14?Math.round(48+gene%6):0;
      const labBatch=i%2===0?5:0;
      const noise=(gene*17+i*23+(i*gene)%13)%29;
      return 53+(gene%31)+biological+labBatch+noise;
    });
    matrix.push('GENE_'+String(gene).padStart(4,'0')+','+vals.join(','));
  }
  return {subjects,metadata:meta.join('\n')+'\n',matrix:matrix.join('\n')+'\n'};
}
function members(data){
  const b=Buffer.from(data),files={};
  let pos=0;
  while(pos+30<b.length&&b.readUInt32LE(pos)===0x04034b50){
    expect(b.readUInt16LE(pos+8)).toBe(0);
    const bytes=b.readUInt32LE(pos+18),len=b.readUInt16LE(pos+26),extra=b.readUInt16LE(pos+28);
    const name=b.subarray(pos+30,pos+30+len).toString('utf8');
    const content=b.subarray(pos+30+len+extra,pos+30+len+extra+bytes);
    expect(files[name]).toBeUndefined();
    files[name]=content;
    pos+=30+len+extra+bytes;
  }
  expect(b.readUInt32LE(pos)).toBe(0x02014b50);
  return files;
}
test('Scientific researcher -> actual R DESeq2 -> complete 172-gene reference CSV -> verifiable local ZIP',async({page})=>{
  // This scenario runs ONLY when CI starts the genuine R Plumber service.
  // General browser-only CI must not be required to install R or pretend to run it.
  test.skip(process.env.PMX_E2E_REAL_R !== '1', 'Requires actual local R backend; run multiomics-live-r-to-zip.yml');
  test.setTimeout(165_000);
  const errors=[],requests=[],responses=[];
  page.on('pageerror',e=>errors.push(e.message));
  page.on('request',r=>{if(r.url()==='http://127.0.0.1:8787/run')requests.push(r)});
  page.on('response',async r=>{
    if(r.url()==='http://127.0.0.1:8787/run'){
      const body=await r.json().catch(()=>null);
      responses.push({status:r.status(),body});
    }
  });
  const {metadata,matrix}=csvInputs();
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:18000}).toBe('true');
  // These CSVs are intentionally NOT the built-in website demo fixtures.
  await page.locator('.uploads input[type="file"]').nth(0).setInputFiles({
    name:'cohort_28_subjects.csv',mimeType:'text/csv',buffer:Buffer.from(metadata)
  });
  await page.locator('.uploads input[type="file"]').nth(1).setInputFiles({
    name:'raw_counts_172_genes.csv',mimeType:'text/csv',buffer:Buffer.from(matrix)
  });
  await page.getByLabel(/Objectif principal|Main objective/).selectOption('groups');
  // Expert-only fields are inside a collapsed disclosure. The beginner
  // independent-group defaults must already be correct; do not hang by
  // trying to interact with a hidden select.
  await expect(page.getByLabel(/Structure du design|Design structure/)).toHaveValue('independent');
  await expect(page.getByLabel(/Design longitudinal|Longitudinal design/)).toHaveValue('no');
  const toggle=page.getByTestId('multiomics-run-options');
  await toggle.locator('summary').click();
  // Explicit opt-in for the real R backend. Never let an E2E default turn it on.
  await page.getByLabel(/Mode backend R|R backend mode/).selectOption('required');
  await expect(page.getByTestId('multiomics-run')).toBeEnabled({timeout:20000});
  await page.getByTestId('multiomics-run').click();
  // Browser exploratory rows can appear BEFORE the awaited R POST finishes.
  // Wait for the RUN button to unlock, indicating completion of the WHOLE
  // analysis, before asserting complete R rows or downloading the ZIP.
  await expect(page.getByTestId('multiomics-run')).toBeEnabled({timeout:120_000});
  const runError=page.locator('section.panel:has(button[data-testid="multiomics-run"]) > p.error');
  if(await runError.isVisible())
    throw Error('Reference R research execution failed: '+await runError.innerText());
  const resultPanel=page.getByTestId('multiomics-results');
  await expect(resultPanel).toBeVisible({timeout:8000});
  await expect(page.getByTestId('multiomics-reference-preview-warning'))
    .toContainText('complet',{timeout:10000});
  expect(requests.length,'A real consented R /run request is required').toBe(1);
  const sent=requests[0].postDataJSON();
  expect(sent.matrices.transcriptomics).toBe(matrix);
  expect(sent.metadataCsv).toBe(metadata);
  expect(sent.dataTypes.transcriptomics).toBe('raw_counts');
  expect(sent.protocol.objective).toBe('groups');
  expect(responses.length).toBe(1);
  expect(responses[0].status).toBe(200);
  expect(responses[0].body?.status).toBe('ok');
  const fromBackend=responses[0].body.methods.transcriptomics_differential_deseq2;
  expect(fromBackend.status).toBe('ok');
  expect(fromBackend.referenceTableComplete).toBe(true);
  expect(fromBackend.referenceTableRows).toBeGreaterThan(100);
  expect(fromBackend.referenceTableRows).toBe(fromBackend.referenceTableExportedRows);
  expect(fromBackend.referenceTableOriginalFeatures).toBe(GENE_COUNT);
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:30000}),
    page.getByTestId('multiomics-research-bundle-download').click()
  ]);
  const content=await fs.readFile(await download.path());
  const entries=members(content);
  const manifest=JSON.parse(entries['provenance/manifest.json'].toString('utf8'));
  const archived=JSON.parse(entries['provenance/analysis_complete.json'].toString('utf8'));
  const project=JSON.parse(entries['source/projet_reproductible.json'].toString('utf8'));
  const method=manifest.r_reference_full_tables.find(x=>
    x.method_key==='transcriptomics_differential_deseq2');
  expect(method).toBeTruthy();
  expect(manifest.r_reference_full_feature_tables_exported).toBe(true);
  expect(manifest.confirmatory_certification).toBe(false);
  expect(manifest.clinical_validation).toBe(false);
  const name=method.name;
  const csv=entries[name].toString('utf8');
  const parsed=parseReferenceCsv(csv);
  expect(parsed.rows.length).toBe(method.rows);
  expect(parsed.rows.length).toBe(fromBackend.referenceTableRows);
  expect(parsed.headers).toEqual(method.columns);
  expect(parsed.headers).toContain('padj');
  expect(parsed.headers).toContain('log2FoldChange');
  expect(parsed.rows[0].feature).toMatch(/^GENE_\d{4}$/);
  expect(parsed.rows.some(x=>x.feature==='GENE_0172')).toBe(true);
  expect(new Set(parsed.rows.map(x=>x.feature)).size).toBe(parsed.rows.length);
  expect(entries[name].toString('utf8')).toBe(fromBackend.referenceTableCsv);
  expect(crypto.createHash('sha256').update(entries[name]).digest('hex'))
    .toBe(manifest.files[name].sha256);
  expect(project.inputs.transcriptomics.content).toBe(matrix);
  expect(project.inputs.metadata.content).toBe(metadata);
  expect(project.archivedResult).toBeNull();
  expect(archived.referenceBackend.methods.transcriptomics_differential_deseq2.referenceTableComplete)
    .toBe(true);
  expect(archived.referenceBackend.environment?.packages?.DESeq2).toBeTruthy();
  expect(entries['tableaux/transcriptomics.csv']).toBeTruthy();
  expect(entries['methods/methods.md'].toString('utf8')).toContain('DESeq2');
  // An INDEPENDENT Python3 verifier reads the same downloaded ZIP and all
  // original-file SHA256 values without extraction or network access.
  const scratch=await fs.mkdtemp(path.join(os.tmpdir(),'pmx-r-e2e-'));
  try{
    const zipfile=path.join(scratch,'study.zip');
    const verifier=path.join(scratch,'verifier_dossier.py');
    await fs.writeFile(zipfile,content);
    await fs.writeFile(verifier,entries['verification/verifier_dossier.py']);
    const check=spawnSync('python3',[verifier,zipfile],{encoding:'utf8'});
    expect(check.status,check.stderr).toBe(0);
    expect(check.stdout).toContain('INTEGRITY PASS');
  }finally{await fs.rm(scratch,{force:true,recursive:true})}
  expect(errors).toEqual([]);
});


function genuinePairedCsvInputs(){
  const assayIds=[];
  const metadata=['subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,sample_type'];
  for(let p=1;p<=4;p++){
    for(const condition of ['control','treated']){
      const subject='BIO_PAIR_'+String(p).padStart(2,'0');
      const assay=subject+'_'+condition;
      assayIds.push(assay);
      metadata.push([subject,assay,assay,'transcriptomics',condition,
         'T0','B1','1','biological'].join(','));
    }
  }
  const raw=['feature_id,'+assayIds.join(',')];
  for(let gene=1;gene<=180;gene++){
    const values=assayIds.map((assay,i)=>{
      const baseline=48+(gene%47);
      const repeatableNoise=(gene*17+i*23+(gene*i)%13)%36;
      const treatmentEffect=assay.endsWith('treated')&&gene<=24?48:0;
      return baseline+repeatableNoise+treatmentEffect;
    });
    raw.push('SIM_GENE_'+String(gene).padStart(4,'0')+','+values.join(','));
  }
  return {metadata:metadata.join('\n')+'\n',
          matrix:raw.join('\n')+'\n'};
}
test('Four genuinely paired biological units -> actual subject-blocked DESeq2 -> complete verified ZIP and independent review only',async({page})=>{
  test.skip(process.env.PMX_E2E_REAL_R!=='1',
    'This test requires the actual R Plumber server, never a mocked /run');
  test.setTimeout(180_000);
  const requests=[],responses=[],pageErrors=[];
  page.on('pageerror',e=>pageErrors.push(e.message));
  page.on('request',r=>{
    if(r.url()==='http://127.0.0.1:8787/run')requests.push(r);
  });
  page.on('response',async r=>{
    if(r.url()==='http://127.0.0.1:8787/run')
      responses.push({httpStatus:r.status(),body:await r.json().catch(()=>null)});
  });
  const input=genuinePairedCsvInputs();
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:18000}).toBe('true');
  await page.locator('.uploads input[type="file"]').nth(0).setInputFiles({
    name:'paired_subject_manifest.csv',mimeType:'text/csv',
    buffer:Buffer.from(input.metadata)});
  await page.locator('.uploads input[type="file"]').nth(1).setInputFiles({
    name:'paired_raw_counts.csv',mimeType:'text/csv',
    buffer:Buffer.from(input.matrix)});
  await page.getByLabel(/Objectif principal|Main objective/).selectOption('groups');
  await page.getByTestId('multiomics-simple-design').selectOption('paired');
  await expect(page.getByTestId('multiomics-paired-count')).toContainText('4',{timeout:20000});
  await expect(page.getByTestId('multiomics-paired-count')).toContainText('8');
  await page.getByTestId('multiomics-analysis-intent').locator('select').selectOption('confirmatory');
  const contrast=page.getByTestId('multiomics-rnaseq-contrast');
  await expect(contrast).toBeVisible();
  await contrast.locator('select').nth(0).selectOption('control');
  await contrast.locator('select').nth(1).selectOption('treated');
  await expect(page.getByTestId('multiomics-rnaseq-fdr-calibration-warning')).toBeVisible();
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-run')).toBeEnabled({timeout:150_000});
  const alert=page.locator('section.panel:has(button[data-testid="multiomics-run"]) > p.error');
  if(await alert.isVisible())
    throw Error('Paired genuine R run was blocked: '+await alert.innerText());
  await expect(page.getByTestId('multiomics-results')).toBeVisible();
  expect(requests.length).toBe(1);
  expect(responses.length).toBe(1);
  expect(responses[0].httpStatus).toBe(200);
  const submitted=requests[0].postDataJSON();
  expect(submitted.protocol.designType).toBe('paired');
  expect(submitted.protocol.analysisIntent).toBe('confirmatory');
  expect(submitted.protocol.rnaPrimaryContrast).toEqual(['condition','treated','control']);
  expect(submitted.dataTypes.transcriptomics).toBe('raw_counts');
  const methods=responses[0].body?.methods||{};
  const de=methods.transcriptomics_differential_deseq2;
  expect(de?.status).toBe('ok');
  const paired=de.pairedDesignContract;
  expect(paired?.pairing).toBe('same_subject_two_conditions');
  expect(paired?.biological_pairs).toBe(4);
  expect(paired?.independent_biological_subjects).toBe(4);
  expect(paired?.model_samples).toBe(8);
  expect(paired?.certified).toBe(false);
  expect(paired?.explicit_primary_contrast).toEqual(['condition','treated','control']);
  expect(de.referenceTableComplete).toBe(true);
  await expect(page.getByTestId('multiomics-confirmatory-readiness')).toBeVisible();
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:30000}),
    page.getByTestId('multiomics-research-bundle-download').click()
  ]);
  const contents=await fs.readFile(await download.path());
  const archive=members(contents);
  const exported=JSON.parse(archive['provenance/analysis_complete.json'].toString('utf8'));
  expect(exported.scientificAssurance?.confirmatoryReadiness?.status)
    .toBe('independent_review_required');
  expect(exported.scientificAssurance?.confirmatoryReadiness?.certified).toBe(false);
  const scientificChecks=exported.scientificAssurance.confirmatoryReadiness.checks;
  expect(scientificChecks.some(c=>c.code==='rnaseq_confirmatory_input_contract'&&c.status==='pass'))
    .toBe(true);
  expect(scientificChecks.some(c=>c.code==='rnaseq_null_fdr_calibration'&&c.status==='review'))
    .toBe(true);
  const manifest=JSON.parse(archive['provenance/manifest.json'].toString('utf8'));
  const method=manifest.r_reference_full_tables.find(x=>
    x.method_key==='transcriptomics_differential_deseq2');
  expect(method).toBeTruthy();
  const csv=archive[method.name].toString('utf8');
  const parsed=parseReferenceCsv(csv);
  expect(parsed.rows.length).toBe(de.referenceTableRows);
  expect(parsed.rows.some(x=>x.feature==='SIM_GENE_0180')).toBe(true);
  expect(new Set(parsed.rows.map(x=>x.feature)).size).toBe(parsed.rows.length);
  expect(crypto.createHash('sha256').update(archive[method.name]).digest('hex'))
    .toBe(manifest.files[method.name].sha256);
  expect(archive['methods/methods.md'].toString('utf8'))
    .toContain('Limites de calibration RNA-seq');
  expect(pageErrors).toEqual([]);
});
