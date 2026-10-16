import assert from 'node:assert/strict';
import {File} from 'node:buffer';
import {createHash} from 'node:crypto';
import {buildResearchBundle} from '../src/lib/multiomics/research-export-package.js';
import {readFileSync,writeFileSync,mkdtempSync,mkdirSync,rmSync,existsSync,copyFileSync} from 'node:fs';
import {join} from 'node:path';
import {tmpdir} from 'node:os';
import {execFileSync} from 'node:child_process';
import {parseDelimited} from '../src/lib/multiomics/deterministic.js';
import {projectNativeRNumbers} from '../src/lib/multiomics/r-primary-results.js';
import {buildMultiomicsVisualizationData} from '../src/lib/multiomics/visualization-data.js';

const evidence=process.env.PMX_PAIRED_AIRWAY_REPORT_DIR;
if(!evidence)throw Error('PMX_PAIRED_AIRWAY_REPORT_DIR required for real public Airway test');
const input=JSON.parse(readFileSync(join(evidence,'paired_airway_original_request.json'),'utf8'));
const native=JSON.parse(readFileSync(join(evidence,'paired_airway_original_response.json'),'utf8'));
assert.equal(native.status,'ok');
assert.equal(native.methods.transcriptomics_differential_deseq2.referenceTableComplete,true);
const stub={
  protocol:{objective:'groups',designType:'paired'},
  layers:{transcriptomics:{
    rows:[{feature:'ENSG00000103196',effect:-999,qValue:0.99}],
    selected:[],steps:['Browser exploratory'],groupSizes:[4,4],effectScale:'log2',
    contrast:'dexamethasone treated vs untreated',
    qc:{preprocessingAudit:{declaredValueType:'raw_counts'}}
  }},
  referenceBackend:native
};
const result=projectNativeRNumbers(stub);
const rows=result.layers.transcriptomics.rows;
const raw= parseDelimited(native.methods.transcriptomics_differential_deseq2.referenceTableCsv);
assert.equal(rows.length,raw.rows.length);
assert(rows.length>10000);
assert.equal(result.numericSourceByLayer.transcriptomics,'native_r_complete_table');
assert.equal(result.browserExploratoryLayers.transcriptomics.rows[0].effect,-999);
function value(s){
  if(s==null||s===''||s==='NA')return null;
  const n=Number(s);return Number.isFinite(n)?n:null;
}
function nearly(a,b){
  if(a===null||b===null)return a===b;
  return Math.abs(a-b)<=1e-10*Math.max(1,Math.abs(a),Math.abs(b));
}
for(let i=0;i<rows.length;i++){
  const x=rows[i], r=raw.rows[i];
  assert.equal(x.feature,r.feature);
  assert(nearly(x.effect,value(r.log2FoldChange)), 'log2FC mismatch '+x.feature);
  assert(nearly(x.pValue,value(r.pvalue)), 'p mismatch '+x.feature);
  assert(nearly(x.qValue,value(r.padj)), 'q mismatch '+x.feature);
}
const marker=rows.find(row=>row.feature==='ENSG00000103196');
assert(marker&&marker.effect>0&&marker.qValue<0.05);
const metadataRows=parseDelimited(input.metadataCsv).rows;
const visualizations=await buildMultiomicsVisualizationData({
  files:{transcriptomics:new File([input.matrices.transcriptomics],'airway_counts.csv')},
  metadataRows,
  columnMapping:input.columnMapping,
  dataTypes:input.dataTypes,
  analysisResult:result,
  maxFeatures:30,maxSamples:48
});
const heatmap=visualizations?.heatmaps?.transcriptomics;
assert(heatmap?.rows?.length>0 && heatmap.samples.length===8);
const root=mkdtempSync(join(tmpdir(),'pmx-public-airway-r-figures-'));
try{
  const json=join(root,'site_display_from_r.json');
  const out=join(root,'figures');
  writeFileSync(json,JSON.stringify({...result,visualizations}));
  execFileSync('Rscript',['multiomics-engine/plot_research_figures.R',json,out],
    {timeout:180000,stdio:'pipe'});
  const csv=parseDelimited(readFileSync(join(out,'heatmap_transcriptomics_plot_values.csv'),'utf8'));
  assert.equal(csv.rows.length,heatmap.rows.length);
  assert.deepEqual(csv.headers.slice(1),heatmap.samples.map(s=>s.sampleId));
  for(let i=0;i<heatmap.rows.length;i++){
    assert.equal(csv.rows[i].feature_id,heatmap.rows[i].feature);
    for(let j=0;j<heatmap.samples.length;j++){
      const expected=heatmap.rows[i].z[j]??null;
      const actual=value(csv.rows[i][csv.headers[j+1]]);
      assert(nearly(expected,actual),
        'R figure plotted value differs from JS site figure at '+i+','+j);
    }
  }
  const figureManifest=JSON.parse(readFileSync(join(out,'figure_manifest.json'),'utf8'));
  assert.equal(figureManifest.visual_contract.heatmap_scale,'row z-score');
  assert.deepEqual(figureManifest.visual_contract.heatmap_color_limits,[-2.5,2.5]);
  assert.equal(figureManifest.visual_contract.missing_values,
    'NA remains missing; never coerced to zero');
  assert.equal(figureManifest.visual_contract.no_statistical_refit,true);
  // Compare the actual color legend returned by the R figure generator to
  // the Svelte diverging palette formula, with ±2 RGB channels for interpolation.
  const legend=parseDelimited(readFileSync(join(out,'heatmap_color_legend.csv'),'utf8'));
  for(const swatch of legend.rows){
    const z=Number(swatch.z_score),ratio=Math.min(1,Math.abs(z)/2.5);
    const white=[247,247,247],target=z>=0?[204,62,62]:[53,103,178];
    const js=white.map((w,i)=>Math.round(w+(target[i]-w)*ratio));
    const hex=swatch.fill_hex.replace('#','');
    const rgb=[0,2,4].map(i=>parseInt(hex.slice(i,i+2),16));
    assert(rgb.every((n,i)=>Math.abs(n-js[i])<=2),
      'Public Airway heatmap JS vs R color difference at z='+z);
  }
  assert(figureManifest.figures.some(f=>f.kind==='heatmap'&&f.filename.endsWith('.pdf')));
  assert(figureManifest.figures.every(f=>f.comparison_to_browser==='scientifically_similar_not_pixel_identical'));
  assert(existsSync(join(out,'heatmap_transcriptomics.pdf')));
  // Real public study: assemble exactly the same research ZIP downloaded by users
  // and launch Rscript run_all.R from the extracted ZIP, with no file edits.
  const hash=(text)=>createHash('sha256').update(text).digest('hex');
  const portableProject={
    format:'pmx-multiomics-portable-project',schemaVersion:1,
    inputs:{
      metadata:{content:input.metadataCsv,sha256:hash(input.metadataCsv)},
      transcriptomics:{content:input.matrices.transcriptomics,
                       sha256:hash(input.matrices.transcriptomics)}
    },
    settings:{organism:'human',unitType:'culture',transcriptomicsValues:'raw_counts',
              transcriptomicsIdType:'ensembl_gene',columnMapping:input.columnMapping}
  };
  const researchResult={
    ...result,engine:{name:'R primary plus independent browser display',version:'test'},
    generatedAt:'2026-10-15T00:00:00Z',visualizations
  };
  const dossier=await buildResearchBundle({
    result:researchResult,portableProject,
    referenceReplay:{
      requestJson:readFileSync(join(evidence,'paired_airway_original_request.json'),'utf8'),
      responseJson:readFileSync(join(evidence,'paired_airway_original_response.json'),'utf8')
    }
  });
  assert.equal(dossier.manifest.r_source_replay.exact_runtime_source_verified,true,
    'Airway original R code must match the SHA-256 of the bundled implementation');
  const zipPath=join(root,'airway_research.zip');
  writeFileSync(zipPath,Buffer.from(await dossier.blob.arrayBuffer()));
  const extracted=join(root,'extracted_dossier');
  execFileSync('python3',['-c',
    'import sys,zipfile;zipfile.ZipFile(sys.argv[1]).extractall(sys.argv[2])',
    zipPath,extracted],{timeout:60000});
  assert(existsSync(join(extracted,'COMMENCER_ICI_FR_EN.md')));
  execFileSync('Rscript',['run_all.R'],{cwd:extracted,timeout:600000,stdio:'pipe'});
  const status=JSON.parse(readFileSync(join(extracted,'research_reproduction_status.json'),'utf8'));
  assert.equal(status.r_numerical_replay,'PASS_R_REFERENCE_ONLY',
    JSON.stringify(status.messages||[]));
  assert.equal(status.r_figures,'GENERATED_SIMILAR_PDF_SVG');
  assert(existsSync(join(extracted,'reproduction','reference_replayed.json')));
  assert(existsSync(join(extracted,'figures_r','output','heatmap_transcriptomics.pdf')));
  // Real biological cohort must also produce an explanatory R Markdown HTML
  // directly from the unmodified exported research ZIP.
  if(process.env.PMX_TEST_RMD==='1'){
    assert(existsSync(join(extracted,'ANALYSE_REPRODUCTIBLE.Rmd')));
    execFileSync('Rscript',['render_notebook.R'],{
      cwd:extracted,timeout:300000,stdio:'pipe'
    });
    const html=readFileSync(join(extracted,'analyse_reproductible.html'),'utf8');
    assert(html.length>3000 && html.includes('MultiOmics') && html.includes('PASS_R_REFERENCE_ONLY'));
    // The knitted HTML must actually display public-cohort values, not just
    // contain a generic "replay succeeded" statement.
    const nativeCsv=parseDelimited(
      readFileSync(join(extracted,'tableaux','transcriptomics.csv'),'utf8'));
    assert.equal(nativeCsv.rows.length,rows.length,
      'Knitted notebook must use the same full-feature native R table');
    assert(/Rows:\s*17199/.test(html.replace(/<[^>]+>/g,' ')),
      'HTML notebook must display its 17,199 full-feature rows');
    for(const entry of nativeCsv.rows.slice(0,6)){
      assert(html.includes(entry.feature),
        'Knitted HTML omits displayed native R feature '+entry.feature);
    }
    const numericAudit=JSON.parse(readFileSync(
      join(extracted,'reproduction','r_numeric_comparison_report.json'),'utf8'));
    assert.equal(numericAudit.result,'PASS',
      'Knitted public dataset must pass complete native R comparison');
    assert.equal(numericAudit.n_discrepancies,0);
    assert(html.length>3000 && existsSync(join(extracted,'research_reproduction_status.json')));
    const rerun=JSON.parse(readFileSync(join(extracted,'research_reproduction_status.json'),'utf8'));
    assert.equal(rerun.r_numerical_replay,'PASS_R_REFERENCE_ONLY');
    assert.equal(rerun.r_figures,'GENERATED_SIMILAR_PDF_SVG');
    // Persist ACTUAL public-cohort output so researchers can inspect the
    // knitted HTML and its numerical R replay audit in Actions artifacts.
    copyFileSync(join(extracted,'analyse_reproductible.html'),
      join(evidence,'airway_annotated_R_markdown_report.html'));
    copyFileSync(join(extracted,'reproduction','r_numeric_comparison_report.json'),
      join(evidence,'airway_native_R_numeric_comparison.json'));
    copyFileSync(join(extracted,'research_reproduction_status.json'),
      join(evidence,'airway_R_markdown_reproduction_status.json'));
    copyFileSync(join(extracted,'figures_r','output','figure_manifest.json'),
      join(evidence,'airway_R_figure_manifest.json'));
  }

    console.log('PASS: Airway/GSE52778',rows.length,
    'native R values reproduced feature-by-feature; browser & R figure plot data are identical for',
    heatmap.rows.length,'features ×',heatmap.samples.length,'samples');
}finally{rmSync(root,{recursive:true,force:true});}
