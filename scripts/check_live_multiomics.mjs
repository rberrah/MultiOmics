// Live deployment verification (not local preview). Runs in GitHub Actions.
// An HTTP 200 can still be an unhydrated JS shell: exercise the actual UI.
import {chromium} from '@playwright/test';
import {createHash} from 'node:crypto';
import {readFile} from 'node:fs/promises';
const origin='https://rberrah.github.io';
const browser=await chromium.launch({headless:true});
try {
  const page=await browser.newPage({viewport:{width:1280,height:850}});
  const errors=[];
  const desktopNetworkFailures=[];
  page.on('pageerror',error=>errors.push(String(error?.message||error)));
  page.on('requestfailed',request=>desktopNetworkFailures.push({
    url:request.url(),reason:request.failure()?.errorText||'unknown'
  }));
  page.on('response',response=>{
    if(response.status()>=400)desktopNetworkFailures.push({
      url:response.url(),http_status:response.status()
    });
  });
  const presentation=await page.goto(origin+'/multiomics/',{
    waitUntil:'domcontentloaded',timeout:30000});
  if(presentation?.status()!==200)
    throw new Error('Presentation HTTP '+presentation?.status());
  await page.getByRole('heading',{level:1})
    .waitFor({state:'visible',timeout:30000});
  const tool=await page.goto(origin+'/multiomics/tool/',{
    waitUntil:'domcontentloaded',timeout:30000});
  if(tool?.status()!==200)
    throw new Error('Tool HTTP '+tool?.status());
  await page.getByTestId('multiomics-quick-demo-rna')
    .waitFor({state:'visible',timeout:30000});
  // Visible server-rendered HTML is not proof that the JS event handlers
  // are ready. Clicking SSR before hydration can silently do nothing.
  try {
    await page.waitForFunction(
      () => document.documentElement.dataset.pmxMultiomicsHydrated==='true',
      null,{timeout:30000}
    );
  } catch (error) {
    console.error('PUBLIC DESKTOP HYDRATION FAILED',JSON.stringify({
      url:page.url(),pageErrors:errors,
      networkFailures:desktopNetworkFailures.slice(0,20)
    },null,2));
    throw error;
  }
  console.log('PASS: real public page loads AND Svelte hydration completed');
  await page.getByTestId('multiomics-quick-demo-rna').click();
  try {
    await page.getByTestId('multiomics-results')
      .waitFor({state:'visible',timeout:30000});
  } catch (error) {
    const evidence=await page.evaluate(() => ({
      url:location.href,
      resultsPresent:!!document.querySelector('[data-testid="multiomics-results"]'),
      demoPresent:!!document.querySelector('[data-testid="multiomics-quick-demo-rna"]'),
      visibleAlerts:[...document.querySelectorAll('[role="alert"], .alert, .error')]
        .filter(el=>el.getClientRects().length)
        .map(el=>el.textContent?.trim().slice(0,350)).filter(Boolean).slice(0,8),
      mainText:document.querySelector('main')?.innerText?.slice(-1400) || ''
    }));
    console.error('PUBLIC DEMO DIAGNOSTICS',JSON.stringify({evidence,pageErrors:errors,networkFailures:desktopNetworkFailures.slice(0,20)},null,2));
    throw error;
  }
  const mode=await page.getByTestId('multiomics-analysis-mode').innerText();
  if(!/Omique unique|Single.omics/i.test(mode))
    throw new Error('Public one-omic demo returned an incomprehensible analysis mode');
  console.log('PASS: public one-omic beginner demo computes results');
  if(errors.length)
    throw new Error('Public desktop JS errors: '+errors.slice(0,5).join(' | '));

  const mobile=await browser.newPage({viewport:{width:390,height:844},
    isMobile:true,hasTouch:true});
  const mobileErrors=[];
  const mobileNetworkFailures=[];
  mobile.on('pageerror',error=>mobileErrors.push(String(error?.message||error)));
  mobile.on('requestfailed',request=>mobileNetworkFailures.push({
    url:request.url(),reason:request.failure()?.errorText||'unknown'
  }));
  mobile.on('response',response=>{
    if(response.status()>=400)mobileNetworkFailures.push({
      url:response.url(),http_status:response.status()
    });
  });
  await mobile.goto(origin+'/multiomics/tool/',{
    waitUntil:'domcontentloaded',timeout:30000});
  await mobile.getByTestId('multiomics-quick-demo-rna')
    .waitFor({state:'visible',timeout:30000});
  try {
    await mobile.waitForFunction(
      () => document.documentElement.dataset.pmxMultiomicsHydrated==='true',
      null,{timeout:30000}
    );
  } catch(error) {
    console.error('PUBLIC MOBILE HYDRATION FAILED',JSON.stringify({
      url:mobile.url(),pageErrors:mobileErrors,
      networkFailures:mobileNetworkFailures.slice(0,20)
    },null,2));
    throw error;
  }
  const overflow=await mobile.evaluate(()=>
    document.documentElement.scrollWidth-document.documentElement.clientWidth);
  if(overflow>12)
    throw new Error('Public mobile page horizontal overflow: '+overflow+'px');
  await mobile.getByTestId('multiomics-quick-demo-rna').click();
  try {
    await mobile.getByTestId('multiomics-results')
      .waitFor({state:'visible',timeout:30000});
  } catch (error) {
    const evidence=await mobile.evaluate(()=>({
      url:location.href,
      demoPresent:!!document.querySelector('[data-testid="multiomics-quick-demo-rna"]'),
      demoDisabled:document.querySelector('[data-testid="multiomics-quick-demo-rna"]')?.disabled,
      resultsPresent:!!document.querySelector('[data-testid="multiomics-results"]'),
      visibleAlerts:[...document.querySelectorAll('[role="alert"], .alert, .error, [role="status"]')]
        .filter(el=>el.getClientRects().length)
        .map(el=>el.textContent?.trim().slice(0,350)).filter(Boolean).slice(0,12),
      mainText:document.querySelector('main')?.innerText?.slice(-1400)||''
    }));
    console.error('PUBLIC MOBILE DEMO DIAGNOSTICS',JSON.stringify({
      evidence,pageErrors:mobileErrors,networkFailures:mobileNetworkFailures
    },null,2));
    throw error;
  }
  console.log('PASS: public mobile 390px beginner demo computes results');
  if(mobileErrors.length)
    throw new Error('Public mobile JS errors: '+mobileErrors.slice(0,5).join(' | '));

  // PR-grade production qualification: genuinely upload synthetically generated
  // vendor-style long MS peak-AUC CSV, compute in the PUBLIC browser, then
  // download and independently verify each original source and ZIP member.
  // Do not use any patient data, browser route mocks or fake reference R.
  const research=await browser.newPage({viewport:{width:1280,height:850},
    acceptDownloads:true});
  const researchErrors=[];
  const outbound=[];
  let auditStarted=false;
  research.on('pageerror',e=>researchErrors.push(String(e?.message||e)));
  research.on('request',req=>{
    const url=req.url();
    // Observe all network requests (including GET/beacons), not just POST:
    // identifiers could leak via URL query parameters. Start AFTER page load
    // so harmless browser/font/CDN initialization is not confused with a
    // researcher-controlled operation.
    if(auditStarted && /^https?:\/\//i.test(url) &&
       new URL(url).origin!==origin)
      outbound.push({method:req.method(),url});
  });
  const response=await research.goto(origin+'/multiomics/tool/',{
    waitUntil:'domcontentloaded',timeout:30000
  });
  if(response?.status()!==200)
    throw new Error('Live researcher tool unavailable (HTTP '+response?.status()+')');
  await research.waitForFunction(
    ()=>document.documentElement.dataset.pmxMultiomicsHydrated==='true',
    null,{timeout:30000});
  const columns=['feature_id','assay_id','peak_area','subject_id',
    'sample_id','condition','batch','sample_type','mz','retention_time'];
  const records=[columns.join(',')];
  for(let sample=1;sample<=12;sample++){
    const assay='PUBLIC_QA_MS_'+String(sample).padStart(2,'0');
    const condition=sample<=6?'control':'treated';
    for(let peak=1;peak<=25;peak++){
      const area=1200+peak*41+sample*(peak%5+4)+
        (sample>6 && peak<=4?900:0);
      records.push(['UNIDENTIFIED_QA_PEAK_'+peak,assay,area,
        'SYNTH_SUBJECT_'+sample,assay,condition,'RUN1','biological',
        (100+peak*.37).toFixed(4),(4+peak*.13).toFixed(4)].join(','));
    }
  }
  const raw=records.join('\n')+'\n';
  auditStarted=true; // Observe every external request after data is supplied.
  await research.getByTestId('multiomics-ms-auc-upload').setInputFiles({
    name:'synthetic_public_smoke_ms_auc.csv',mimeType:'text/csv',
    buffer:Buffer.from(raw,'utf8')
  });
  const status=research.getByTestId('multiomics-ms-auc-import-status');
  await status.waitFor({state:'visible',timeout:30000});
  await research.waitForFunction(()=>{
    const x=document.querySelector('[data-testid="multiomics-ms-auc-import-status"]');
    return x?.getAttribute('data-import-state')==='ready';
  },null,{timeout:30000});
  const label=await status.innerText();
  if(!/25/.test(label)||!/12/.test(label))
    throw new Error('Public MS AUC import did not preserve 25 features × 12 samples: '+label);
  await research.getByTestId('multiomics-run').click();
  await research.getByTestId('multiomics-results')
    .waitFor({state:'visible',timeout:30000});
  const [download]=await Promise.all([
    research.waitForEvent('download',{timeout:30000}),
    research.getByTestId('multiomics-research-bundle-download').click()
  ]);
  if(download.suggestedFilename()!=='multiomics_dossier_recherche_v1.zip')
    throw new Error('Unexpected researcher ZIP name');
  const zip=await readFile(await download.path());
  const members=new Map();
  let offset=0;
  while(offset+30<=zip.length && zip.readUInt32LE(offset)===0x04034b50){
    if(zip.readUInt16LE(offset+8)!==0)
      throw new Error('Unexpected compressed or encrypted research ZIP member');
    const size=zip.readUInt32LE(offset+18);
    const nameSize=zip.readUInt16LE(offset+26),extra=zip.readUInt16LE(offset+28);
    const name=zip.subarray(offset+30,offset+30+nameSize).toString('utf8');
    if(!name||members.has(name)||name.includes('..'))
      throw new Error('Duplicate or invalid research ZIP member '+name);
    const start=offset+30+nameSize+extra;
    if(start+size>zip.length)throw new Error('Truncated ZIP member '+name);
    members.set(name,zip.subarray(start,start+size));
    offset=start+size;
  }
  if(offset+4>zip.length||zip.readUInt32LE(offset)!==0x02014b50)
    throw new Error('Research ZIP central directory missing');
  const required=(path)=>{
    const bytes=members.get(path);
    if(!bytes)throw new Error('Missing research ZIP member '+path);
    return bytes;
  };
  const manifest=JSON.parse(required('provenance/manifest.json').toString('utf8'));
  const original=JSON.parse(required('source/projet_reproductible.json').toString('utf8'));
  const analysis=JSON.parse(required('provenance/analysis_complete.json').toString('utf8'));
  const methods=JSON.parse(required('methods/methods_audit.json').toString('utf8'));
  const source=original.inputs?.msOriginal;
  if(source?.filename!=='synthetic_public_smoke_ms_auc.csv' ||
     source.content!==raw)throw new Error('Research ZIP lost exact uploaded vendor input');
  const digest=(bytes)=>createHash('sha256').update(bytes).digest('hex');
  if(source.sha256!==digest(Buffer.from(raw,'utf8')))
    throw new Error('Original LC-MS CSV source SHA-256 mismatch');
  const files=manifest.files;
  if(!files||Object.keys(files).length!==members.size-1)
    throw new Error('Research dossier missing file-level integrity manifest');
  for(const [name,properties] of Object.entries(files)){
    const bytes=required(name);
    if(properties.sha256!==digest(bytes)||properties.bytes!==bytes.length)
      throw new Error('Research ZIP file SHA-256 or byte count mismatch: '+name);
  }
  if(analysis.msImportProvenance?.originalSource?.importFormat!=='long_ms_auc'||
     analysis.engine?.analysisMode!=='single_omic'||
     analysis.layers?.metabolomics?.rows?.length!==25||
     manifest.browser_layer_rows?.metabolomics!==25||
     !required('tableaux/metabolomics.csv').toString('utf8').includes('UNIDENTIFIED_QA_PEAK_25'))
    throw new Error('Public MS results lost original 25 peak features or provenance');
  if(manifest.confirmatory_certification!==false||
     manifest.clinical_validation!==false||
     methods.scientific_certification!==false||
     original.archivedResult!==null)
    throw new Error('Public researcher ZIP incorrectly certifies or freezes scientific findings');
  if(outbound.length || researchErrors.length)
    throw new Error('Public research flow attempted external data transfer or hit JS errors: '+
      JSON.stringify({outbound,researchErrors}));
  console.log('PASS: PUBLIC synthetic MS peak-AUC 300 rows / 25 features / 12 samples -> browser analysis -> full source-preserving researcher ZIP, SHA-256 per member; no external GET/POST requests after upload');
  await research.close();
} finally {
  await browser.close();
}
