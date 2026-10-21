// @ts-nocheck — Playwright's browser-side DOM evaluation is runtime-validated.
import { test, expect } from '@playwright/test';

test('multi-omics one-click demo exposes guided heatmap, pathway and map views', async ({ page }) => {
  await page.route('https://reactome.org/AnalysisService/**', async (route) => {
    await route.fulfill({
      status: 200,
      contentType: 'application/json',
      body: JSON.stringify({ summary: { token: 'FIG4_VISUAL_TEST' }, pathways: [], pathwaysFound: 0 })
    });
  });
  await page.goto('/multiomics/tool');

  await expect(page.getByTestId('multiomics-quick-start')).toBeVisible();
  await page.getByTestId('multiomics-quick-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 30_000 });
  await expect(page.getByTestId('multiomics-visualizations')).toBeVisible();
  await expect(page.getByTestId('multiomics-heatmap')).toBeVisible();

  await page.locator('.figure-switcher button').nth(1).click();
  await expect(page.getByTestId('multiomics-metabologram')).toBeVisible();
  const legend = page.getByTestId('multiomics-readable-legend');
  await expect(legend).toBeVisible();
  await expect(legend).toContainText('L-tryptophane');
  await expect(legend).toContainText('L-kynurénine');
  await expect(legend).toContainText('Lactate');
  await expect(legend).toContainText('Succinate');
  await expect(legend).toContainText('CHEBI:24996');
  await expect(legend.locator('.legend-entry')).toHaveCount(9);
  // The portal supports light/dark schemes; the plot is always on white.
  // Never inherit white theme text onto the white annotation panel.
  for (const scheme of ['light', 'dark']) {
    await page.emulateMedia({ colorScheme: scheme });
    const contrast = await legend.evaluate((container) => {
      const color = (el, attribute) => {
        const raw = getComputedStyle(el)[attribute];
        const channels = raw.match(/\d+/g)?.slice(0, 3).map(Number) || [];
        return channels.map((x) => x / 255).map((x) => x <= .04045 ? x / 12.92 : ((x + .055) / 1.055) ** 2.4);
      };
      const luminance = (values) => values[0] * .2126 + values[1] * .7152 + values[2] * .0722;
      const ratio = (fg, bg) => {
        const a = luminance(color(fg, 'color')), b = luminance(color(bg, 'backgroundColor'));
        return (Math.max(a, b) + .05) / (Math.min(a, b) + .05);
      };
      const targets = [
        ...container.querySelectorAll('.legend-name strong, .legend-name small, .legend-effect')
      ];
      return {
        headingRatio: ratio(container.querySelector('.legend-section h4'), container.querySelector('.legend-section h4')),
        minTextRatio: Math.min(...targets.map((node) => ratio(node, node.closest('.legend-entry')))),
        surface: getComputedStyle(container).backgroundColor,
        count: targets.length
      };
    });
    expect(contrast.surface).toBe('rgb(255, 255, 255)');
    expect(contrast.headingRatio).toBeGreaterThanOrEqual(7);
    expect(contrast.minTextRatio).toBeGreaterThanOrEqual(7);
    expect(contrast.count).toBeGreaterThan(0);
  }
  await page.getByTestId('multiomics-pathway-select').selectOption('glycolysis');
  await expect(page.getByTestId('multiomics-pathway-coverage')).toContainText(/Glycolyse|Glycolysis/);

  await page.locator('.figure-switcher button').nth(2).click();
  const focused = page.getByTestId('multiomics-focused-network');
  await expect(focused).toBeVisible();
  await expect(page.getByTestId('multiomics-network-coverage')).toContainText(/variables métabolomiques appariées|metabolomics features mapped/i);
  await expect(page.getByTestId('multiomics-network-coverage')).toContainText(/transcrits projetés|transcripts projected/i);
  await expect(focused).toContainText(/automatisée mais non exhaustive|automated but not exhaustive/i);
  await expect(focused).toContainText(/Protéomique non projetée|Proteomics is not yet projected/i);
  await expect(focused).toContainText(/Tryptophane et kynurénine|Tryptophan and kynurenine/);
  await expect(focused).toContainText(/Glycolyse|Glycolysis/);
  await expect(focused).toContainText(/Cycle de Krebs|TCA cycle/);
  await expect(focused.getByTestId('multiomics-region-serotonin')).toHaveCount(0);
  await focused.getByTestId('multiomics-network-hops').selectOption('2');
  await expect(focused.getByTestId('multiomics-region-kynurenine')).toBeVisible();
  await focused.getByTestId('multiomics-network-region').selectOption('kynurenine');
  await expect(focused.getByTestId('multiomics-region-tca')).toHaveCount(0);
  await focused.getByTestId('multiomics-network-dictionary').locator('summary').click();
  await expect(focused.getByTestId('multiomics-network-dictionary')).toContainText('CHEBI:16828');
  await page.getByTestId('multiomics-central-carbon-details').locator('summary').first().click();
  await expect(page.getByTestId('multiomics-central-carbon-map')).toBeVisible();
  await expect(page.getByTestId('multiomics-central-carbon-map')).toContainText(/13C/);
  await expect(page.getByTestId('multiomics-off-map')).toContainText('L-tryptophane');
  await expect(page.getByTestId('multiomics-off-map')).toContainText('L-kynurénine');
  await expect(page.getByTestId('multiomics-central-carbon-map').locator('.map-summary')).toContainText('2');
  await expect(page.getByTestId('multiomics-central-carbon-map')).toContainText('LDHA');
  await expect(page.getByTestId('multiomics-central-carbon-map')).toContainText('LDHB');

  await page.locator('.figure-switcher button').first().click();
  await expect(page.getByTestId('multiomics-heatmap')).toBeVisible();

  // Inspect foreground/background of every figure—not just the metabologram
  // key—under both display color schemes. SVG annotations have explicit ink.
  for (const scheme of ['light','dark']) {
    await page.emulateMedia({colorScheme:scheme});
    for (const figure of ['multiomics-heatmap','multiomics-metabologram','multiomics-central-carbon-map']) {
      if (figure === 'multiomics-heatmap') await page.locator('.figure-switcher button').first().click();
      if (figure === 'multiomics-metabologram') await page.locator('.figure-switcher button').nth(1).click();
      if (figure === 'multiomics-central-carbon-map') {
        await page.locator('.figure-switcher button').nth(2).click();
        const disclosure=page.getByTestId('multiomics-central-carbon-details');
        if (await disclosure.getAttribute('open')===null) await disclosure.locator('summary').first().click();
      }
      const target=page.getByTestId(figure);
      await expect(target).toBeVisible();
      const measured=await target.evaluate((root)=>{
        const luminance=(rgb)=>{
          const c=(rgb.match(/\d+/g)||[]).slice(0,3).map(Number).map(v=>v/255)
            .map(v=>v<=.04045?v/12.92:((v+.055)/1.055)**2.4);
          return .2126*c[0]+.7152*c[1]+.0722*c[2];
        };
        const contrast=(a,b)=>(Math.max(a,b)+.05)/(Math.min(a,b)+.05);
        const bg=luminance(getComputedStyle(root).backgroundColor);
        const title=root.querySelector('.figure-title h3');
        const sample=root.querySelector('svg text.svg-note, svg text.pathway-label, svg text.half-label');
        return {background:getComputedStyle(root).backgroundColor,
          title:contrast(luminance(getComputedStyle(title).color),bg),
          svg:sample?contrast(luminance(getComputedStyle(sample).fill),1):null};
      });
      console.log('FIGURE_CONTRAST', scheme, figure, JSON.stringify(measured));
      expect(measured.background, figure + ' plot background (' + scheme + ')').toBe('rgb(255, 255, 255)');
      expect(measured.title, figure + ' heading contrast (' + scheme + ')').toBeGreaterThanOrEqual(7);
      if(measured.svg!==null)expect(measured.svg, figure + ' SVG annotation contrast (' + scheme + ')').toBeGreaterThanOrEqual(7);
    }
    await page.locator('.figure-switcher button').nth(2).click();
    const graph=page.getByTestId('multiomics-focused-network');
    const ink=await graph.evaluate(node=>({
      color:getComputedStyle(node).color,background:getComputedStyle(node).backgroundColor
    }));
    expect(ink.color).toBe('rgb(24, 43, 53)');
    expect(ink.background).toBe('rgb(255, 255, 255)');
  }
});

test('MS peak-AUC long CSV can be uploaded directly with group metadata', async ({page})=>{
  await page.goto('/multiomics/tool');
  await expect.poll(()=>page.evaluate(()=>document.documentElement.dataset.pmxMultiomicsHydrated), {timeout:15000}).toBe('true');
  const header='feature_id,assay_id,auc,subject_id,sample_id,condition';
  const rows=[header];
  for (let i=0;i<8;i++)for(const f of ['CHEBI:24996','CHEBI:30031'])
    rows.push([f,'INJ'+i,100+i*20+(f==='CHEBI:30031'?100:0),'SUBJ'+i,'S'+i,i<4?'control':'treated'].join(','));
  await page.getByTestId('multiomics-ms-auc-upload').setInputFiles({
    name:'ms_auc_export.csv',mimeType:'text/csv',buffer:Buffer.from(rows.join('\n'))
  });
  await expect(page.getByTestId('multiomics-ms-auc-import-status')).toHaveAttribute('data-import-state', 'ready', {timeout: 15_000});
  await expect(page.getByTestId('multiomics-ms-auc-import-status')).toContainText(/2 signaux MS|2 MS features/);
  await expect(page.getByTestId('multiomics-ms-auc-import-status')).toContainText(/8 injections|8 injections/);
  await expect(page.getByTestId('multiomics-simple-values')).toBeVisible();
  // A previous successful import must not survive a subsequently invalid CSV.
  const upload=page.getByTestId('multiomics-ms-auc-upload');
  await upload.setInputFiles({name:'invalid_ms_auc.csv',mimeType:'text/csv',
    buffer:Buffer.from(header+'\nCHEBI:24996,INJ0,not-a-number,SUBJ0,S0,control')});
  await expect(page.getByTestId('multiomics-ms-auc-import-status'))
    .toHaveAttribute('data-import-state','error');
  await expect(page.getByTestId('multiomics-ms-auc-import-status'))
    .toContainText(/MS AUC|import MS|MS import/i);
  await upload.setInputFiles({name:'recovered_ms_auc.csv',mimeType:'text/csv',
    buffer:Buffer.from(rows.join('\n'))});
  await expect(page.getByTestId('multiomics-ms-auc-import-status'))
    .toHaveAttribute('data-import-state','ready');
  await expect(page.getByTestId('multiomics-ms-auc-import-status'))
    .toContainText(/2 signaux MS|2 MS features/);
});


test('MS peak-AUC rapid replacement keeps only the newest import and its group metadata', async ({page})=>{
  await page.goto('/multiomics/tool');
  await expect.poll(()=>page.evaluate(()=>document.documentElement.dataset.pmxMultiomicsHydrated),
    {timeout:15000}).toBe('true');
  await page.evaluate(()=>{
    const input=document.querySelector('[data-testid="multiomics-ms-auc-upload"]');
    if(!(input instanceof HTMLInputElement))throw new Error('MS AUC uploader not found');
    const header='feature_id,assay_id,auc,subject_id,sample_id,condition';
    const replace=(name,body)=>{
      const data=new DataTransfer();
      data.items.add(new File([body],name,{type:'text/csv'}));
      input.files=data.files;
      input.dispatchEvent(new Event('change',{bubbles:true}));
    };
    // Both change events execute before either asynchronous File.text() finishes.
    // The old invalid import must never replace the newer, valid one.
    replace('obsolete.csv',header+'\nCHEBI:24996,OLD,invalid,OLD,OLD,control');
    const rows=[header];
    for(let i=0;i<8;i++)for(const feature of ['CHEBI:24996','CHEBI:30031'])
      rows.push([feature,'INJ'+i,100+i*5,'SUBJ'+i,'S'+i,i<4?'control':'treated'].join(','));
    replace('latest.csv',rows.join('\n'));
  });
  const status=page.getByTestId('multiomics-ms-auc-import-status');
  await expect(status).toHaveAttribute('data-import-state','ready',{timeout:15000});
  await expect(status).toContainText(/2 signaux MS|2 MS features/);
  await expect(status).toContainText(/8 injections/);
  await expect(page.getByTestId('multiomics-simple-values')).toBeVisible();
});

test('expert raw-MS and external validation workflows are disclosed on request', async ({ page }) => {
  await page.goto('/multiomics/tool');
  const expert = page.getByTestId('advanced-workflows');
  await expect(expert).toBeVisible();
  await expect(expert).not.toHaveAttribute('open', '');
  await expert.locator('summary').click();
  await expect(expert).toHaveAttribute('open', '');
  await expect(expert).toContainText('mzML');
  await expect(page.getByTestId('raw-ms-manifest-template')).toHaveAttribute('href', /raw_ms_manifest[.]csv$/);
});

test('simple mode keeps expert settings hidden until requested, without losing settings', async ({ page }) => {
  await page.goto('/multiomics/tool');
  const advanced = page.getByTestId('multiomics-advanced-settings');
  await expect(advanced).toBeVisible();
  await expect(advanced).not.toHaveAttribute('open', '');
  await expect(page.getByTestId('multiomics-simple-study')).toBeVisible();
  await expect(page.getByTestId('multiomics-simple-design')).toBeVisible();

  // Choosing repeated subjects is explicit; the expert form reflects the same design.
  await page.getByTestId('multiomics-simple-design').selectOption('repeated');
  await advanced.locator('summary').first().click();
  await expect(advanced).toHaveAttribute('open', '');
  await expect(page.getByLabel(/Structure du design|Design structure/)).toHaveValue('repeated');
  await expect(page.getByLabel(/Design longitudinal|Longitudinal design/)).toHaveValue('yes');
});

test('demo measurement scales are visible while technical tables start collapsed', async ({ page }) => {
  await page.route('https://reactome.org/AnalysisService/**', async route => route.fulfill({
    status: 200, contentType: 'application/json',
    body: JSON.stringify({ summary:{token:'SIMPLE_DEMO'},pathwaysFound:0,pathways:[] })
  }));
  await page.goto('/multiomics/tool');
  await page.getByTestId('multiomics-load-demo').click();
  await expect(page.getByTestId('multiomics-results')).toBeVisible({ timeout: 25000 });
  const types = page.getByTestId('multiomics-simple-values');
  await expect(types).toBeVisible();
  await expect(types.getByLabel(/Type de valeurs des gènes|Gene value type/)).toHaveValue('raw_counts');
  await expect(types.getByLabel(/Type de valeurs des protéines|Protein value type/)).toHaveValue('log_intensity');
  await expect(types.getByLabel(/Type de valeurs des métabolites|Metabolite value type/)).toHaveValue('peak_area');
  await expect(page.getByTestId('multiomics-visualizations')).toBeVisible();
  await expect(page.getByTestId('multiomics-interpretation')).toBeVisible();
  await expect(page.getByTestId('multiomics-quality-details')).not.toHaveAttribute('open', '');
  await expect(page.getByTestId('multiomics-detailed-results')).not.toHaveAttribute('open', '');
  await page.getByTestId('multiomics-quality-details').locator('summary').first().click();
  await expect(page.getByTestId('multiomics-qc')).toBeVisible();
});

test('MS peak annotations upload is explicit and rejects unknown chemical IDs', async ({ page }) => {
  await page.goto('/multiomics/tool');
  await page.waitForLoadState('networkidle');
  const upload=page.getByTestId('multiomics-ms-annotation-upload');
  await upload.setInputFiles({
    name:'peak_annotations.csv',
    mimeType:'text/csv',
    buffer:Buffer.from('feature_id,chebi_id\nInternalPeak43,CHEBI:24996\n')
  });
  await expect(page.getByTestId('multiomics-ms-annotation-count')).toContainText('1');
  await expect(page.getByTestId('multiomics-ms-annotation-error')).toHaveCount(0);
  await upload.setInputFiles({
    name:'bad_annotations.csv',mimeType:'text/csv',
    buffer:Buffer.from('feature_id,chebi_id\nInternalPeak43,NotKnown\n')
  });
  await expect(page.getByTestId('multiomics-ms-annotation-error')).toContainText('CHEBI');
  await expect(page.getByTestId('multiomics-ms-annotation-count')).toHaveCount(0);
});
