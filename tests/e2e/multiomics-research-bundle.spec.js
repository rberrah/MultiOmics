// @ts-nocheck
import {test,expect} from '@playwright/test';
import fs from 'node:fs/promises';

function unzipStore(bytes){
  const archive=Buffer.from(bytes),result={};
  let offset=0;
  while(offset+30<=archive.length&&archive.readUInt32LE(offset)===0x04034b50){
    const method=archive.readUInt16LE(offset+8);
    const compressed=archive.readUInt32LE(offset+18);
    const n=archive.readUInt16LE(offset+26),extra=archive.readUInt16LE(offset+28);
    expect(method).toBe(0);
    const name=archive.subarray(offset+30,offset+30+n).toString('utf8');
    result[name]=archive.subarray(offset+30+n+extra,offset+30+n+extra+compressed).toString('utf8');
    offset+=30+n+extra+compressed;
  }
  expect(archive.readUInt32LE(offset)).toBe(0x02014b50);
  return result;
}
for(const [demo,expected] of [
  ['multiomics-quick-demo-rna',['transcriptomics']],
  ['multiomics-quick-demo-protein',['proteomics']],
  ['multiomics-quick-demo-metabolite',['metabolomics']],
  ['multiomics-quick-demo',['transcriptomics','proteomics','metabolomics']]
]){
  test('Scientific ZIP includes complete data and original sources: '+demo,async({page})=>{
    test.setTimeout(70_000);
    const errors=[];
    page.on('pageerror',e=>errors.push(e.message));
    await page.goto('/multiomics/tool?lang=fr');
    await expect.poll(()=>page.evaluate(()=>
      document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:12000}).toBe('true');
    await page.getByTestId(demo).click();
    await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
    const [download]=await Promise.all([
      page.waitForEvent('download',{timeout:30000}),
      page.getByTestId('multiomics-research-bundle-download').click()
    ]);
    expect(download.suggestedFilename()).toBe('multiomics_dossier_recherche_v1.zip');
    // Confirm browser completion while the page session is still active, before
    // parsing/validating the ZIP outside the browser process.
    await expect(page.getByTestId('multiomics-research-bundle-status')).toBeVisible({timeout:15000});
    const contents=unzipStore(await fs.readFile(await download.path()));
    const manifest=JSON.parse(contents['provenance/manifest.json']);
    const project=JSON.parse(contents['source/projet_reproductible.json']);
    const analysis=JSON.parse(contents['provenance/analysis_complete.json']);
    const audit=JSON.parse(contents['methods/methods_audit.json']);
    expect(manifest.confirmatory_certification).toBe(false);
    expect(manifest.sensitive_research_data).toBe(true);
    expect(manifest.downstream_exchange.status).toBe('included');
    expect(manifest.downstream_exchange.omics2pk_readiness.pk_dose_concentration_data_included).toBe(false);
    expect(manifest.downstream_exchange.omics2pk_readiness.status).toBe('source_exchange_ready_not_pk_ready');
    expect(audit.scientific_certification).toBe(false);
    expect(audit.clinical_validation).toBe(false);
    expect(project.archivedResult).toBeNull();
    expect(project.inputs.metadata.sha256).toMatch(/^[a-f\d]{64}$/);
    expect(Object.keys(analysis.layers).sort()).toEqual([...expected].sort());
    expect(Object.keys(manifest.browser_layer_rows).sort()).toEqual([...expected].sort());
    for(const layer of expected){
      const csv=contents['tableaux/'+layer+'.csv'];
      expect(csv).toBeTruthy();
      const lines=csv.replace(/^\uFEFF/,'').split('\r\n');
      expect(lines[0]).toContain('feature,effect');
      expect(manifest.browser_layer_rows[layer]).toBe(analysis.layers[layer].rows.length);
      // CSV has one entry per full JSON row for the provided demo datasets.
      expect(lines.length-2).toBe(analysis.layers[layer].rows.length);
    }
    expect(contents['methods/methods.md']).toContain('Méthodes exécutées');
    expect(errors).toEqual([]);
  });
}


// A completely user-supplied measurement file, rather than one of the built-in
// demo datasets. Assumption "one measurement column = one independent subject"
// is explicitly confirmed by the researcher; no group labels are invented.
test('External researcher imports own RNA matrix without metadata, exports complete provenance and never invents a contrast',async({page})=>{
  test.setTimeout(85_000);
  const errs=[];
  page.on('pageerror',e=>errs.push(e.message));
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  const rows=['feature_id,'+Array.from({length:8},(_,i)=>'S'+String(i+1).padStart(2,'0')).join(',')];
  for(let gene=1;gene<=45;gene++){
    rows.push('EXT_GENE_'+gene+','+Array.from({length:8},(_,i)=>
      String(25+gene*3+i*((gene%3)+1))).join(','));
  }
  const raw=rows.join('\n')+'\n';
  await page.locator('.uploads input[type="file"]').nth(1).setInputFiles({
    name:'independent_laboratory_expression.csv',mimeType:'text/csv',
    buffer:Buffer.from(raw,'utf8')
  });
  const helper=page.getByTestId('multiomics-single-sheet-helper');
  await expect(helper).toBeVisible({timeout:15000});
  await helper.getByRole('checkbox').check();
  await page.getByTestId('multiomics-generate-explore-metadata').click();
  await expect(page.getByTestId('multiomics-generated-sheet-notice')).toBeVisible();
  await expect(page.getByTestId('multiomics-generated-sheet-notice')).toContainText('aucun groupe inféré');
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:25000}),
    page.getByTestId('multiomics-research-bundle-download').click()
  ]);
  const contents=unzipStore(await fs.readFile(await download.path()));
  const project=JSON.parse(contents['source/projet_reproductible.json']);
  const details=JSON.parse(contents['provenance/analysis_complete.json']);
  const manifest=JSON.parse(contents['provenance/manifest.json']);
  const methods=JSON.parse(contents['methods/methods_audit.json']);
  expect(project.inputs.transcriptomics.filename).toBe('independent_laboratory_expression.csv');
  expect(project.inputs.transcriptomics.content).toBe(raw);
  expect(project.inputs.transcriptomics.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(project.inputs.metadata.content).toContain('subject_id');
  expect(details.protocol.objective).toBe('explore');
  expect(methods.objective).toBe('explore');
  expect(methods.clinical_validation).toBe(false);
  expect(methods.scientific_certification).toBe(false);
  expect(project.archivedResult).toBeNull();
  expect(manifest.browser_layer_rows.transcriptomics).toBe(details.layers.transcriptomics.rows.length);
  expect(contents['tableaux/transcriptomics.csv']).toContain('EXT_GENE_45');
  expect(contents['methods/methods.md']).toContain('Comprendre les colonnes');
  expect(errs).toEqual([]);
});


test('Researcher-provided long LC-MS AUC keeps original vendor rows and uncertain peak annotations in a reproducible ZIP',async({page})=>{
  test.setTimeout(100_000);
  const errors=[];
  page.on('pageerror',e=>errors.push(e.message));
  const cols=['feature_id','assay_id','peak_area','subject_id','sample_id',
              'condition','batch','sample_type','mz','retention_time'];
  const lines=[cols.join(',')];
  for(let sample=1;sample<=12;sample++){
    const sid='MS_SAMPLE_'+String(sample).padStart(2,'0');
    const condition=sample<=6?'control':'treated';
    for(let peak=1;peak<=25;peak++){
      const area=Number((1200+peak*41+sample*(peak%5+4)+
        (sample>6&&peak<=4?900:0)).toFixed(1));
      lines.push([ 'UNKNOWN_PEAK_'+peak,sid,area,
         'ANON_SUBJECT_'+sample,sid,condition,'RUN1','biological',
         (100+peak*.37).toFixed(4),(4+peak*.13).toFixed(4)].join(','));
    }
  }
  const vendorOriginal=lines.join('\n')+'\n';
  await page.goto('/multiomics/tool?lang=fr');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:15000}).toBe('true');
  await page.getByTestId('multiomics-ms-auc-upload').setInputFiles({
    name:'realistic_vendor_peak_auc.csv',mimeType:'text/csv',
    buffer:Buffer.from(vendorOriginal,'utf8')
  });
  await expect(page.getByTestId('multiomics-ms-auc-import-status'))
    .toBeVisible({timeout:18000});
  await page.getByTestId('multiomics-run').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:25000}),
    page.getByTestId('multiomics-research-bundle-download').click()
  ]);
  const zip=unzipStore(await fs.readFile(await download.path()));
  const source=JSON.parse(zip['source/projet_reproductible.json']);
  const output=JSON.parse(zip['provenance/analysis_complete.json']);
  const manifest=JSON.parse(zip['provenance/manifest.json']);
  expect(source.inputs.msOriginal.filename).toBe('realistic_vendor_peak_auc.csv');
  expect(source.inputs.msOriginal.content).toBe(vendorOriginal);
  expect(source.inputs.msOriginal.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(source.inputs.metabolomics.content).toContain('UNKNOWN_PEAK_25');
  expect(source.inputs.metadata.content).toContain('ANON_SUBJECT_12');
  expect(output.msImportProvenance.originalSource.importFormat).toBe('long_ms_auc');
  expect(output.engine.analysisMode).toBe('single_omic');
  expect(output.layers.metabolomics.rows.length).toBeGreaterThan(0);
  expect(manifest.browser_layer_rows.metabolomics).toBe(output.layers.metabolomics.rows.length);
  expect(zip['tableaux/metabolomics.csv']).toContain('UNKNOWN_PEAK_25');
  expect(manifest.confirmatory_certification).toBe(false);
  expect(manifest.clinical_validation).toBe(false);
  expect(errors).toEqual([]);
});


test('Direct Omics2PK download contains only strict source exchange, never invented PK outcomes',async({page})=>{
  test.setTimeout(70_000);
  await page.goto('/multiomics/tool?lang=en');
  await expect.poll(()=>page.evaluate(()=>
    document.documentElement.dataset.pmxMultiomicsHydrated),{timeout:12000}).toBe('true');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({timeout:30000});
  await expect(page.getByTestId('multiomics-omics2pk-export-panel')).toContainText('does not export group logFC/q-values');
  const [download]=await Promise.all([
    page.waitForEvent('download',{timeout:25000}),
    page.getByTestId('multiomics-omics2pk-export').click()
  ]);
  expect(download.suggestedFilename()).toBe('multiomics_omics2pk_exchange_v1.zip');
  const zip=unzipStore(await fs.readFile(await download.path()));
  const manifest=JSON.parse(zip['manifest.json']);
  expect(manifest.schema_version).toBe('1.0.0');
  expect(manifest.validation.status).toBe('exploratory');
  expect(manifest.files.samples.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(manifest.files.features.sha256).toMatch(/^[a-f0-9]{64}$/);
  expect(manifest.files.individual_signatures).toBeUndefined();
  expect(manifest.matrices.map(x=>x.omic).sort()).toEqual(
    ['metabolomics','proteomics','transcriptomics']
  );
  expect(zip['samples.csv']).toContain('timepoint');
  expect(zip['samples.csv']).toContain('T12');
  expect(zip['individual_signatures.csv']).toBeUndefined();
  await expect(page.getByTestId('multiomics-omics2pk-export-status')).toContainText(
    'PK dose/concentration data are not included'
  );
});
