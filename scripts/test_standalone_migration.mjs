// Zero-dependency static migration preflight; does not certify statistics.
// Usage: node scripts/test_standalone_migration.mjs
import assert from 'node:assert/strict';
import { readFileSync, existsSync } from 'node:fs';

const source = (path) => {
  assert(existsSync(path), 'Missing required standalone source: '+path);
  return readFileSync(path,'utf8');
};
const config=source('svelte.config.js');
assert(/trailingSlash\s*:\s*['"]always['"]/.test(config),
  'GitHub Pages must receive directory index routes');
assert(config.includes("process.env.BASE_PATH"),
  'GitHub Pages project prefix must be configurable');
const packageFile=JSON.parse(source('package.json'));
assert.equal(packageFile.name,'multiomics');
assert(packageFile.scripts.check&&packageFile.scripts.build&&packageFile.scripts['test:multiomics']);
const html=source('src/app.html');
for(const name of [
  'guided-ui.js','plain-language-ui.js','validation-evidence-ui.js',
  'design-guide-ui.js','file-guide-ui.js','fdr-policy-ui.js'
]) {
  assert(html.includes('/multiomics/'+name),'Missing research UX module: '+name);
  source('static/multiomics/'+name);
}
for(const path of [
 'src/routes/+layout.svelte',
 'src/lib/stores/language.js',
 'src/routes/+page.svelte',
 'src/routes/multiomics/+page.svelte',
 'src/routes/multiomics/tool/+page.svelte',
 'src/routes/multiomics/learn/+page.svelte',
 'src/routes/multiomics/demo/+page.svelte',
 'multiomics-engine/server.R',
 'multiomics-engine/run_backend.R',
 'multiomics-engine/run_research_bundle.R',
 'multiomics-engine/verify_research_replay.R',
 'multiomics-engine/research_reproducibility_notebook.Rmd',
 'multiomics-engine/render_research_notebook.R',
 'tests/e2e/multiomics-methods-learning.spec.js'
]) source(path);
const page=source('src/routes/multiomics/tool/+page.svelte');
for(const name of [
  'multiomics-quick-demo-rna',
  'multiomics-quick-demo-protein',
  'multiomics-quick-demo-metabolite',
  'multiomics-quick-demo'
]) assert(page.includes('data-testid="'+name+'"'),'Quick demo missing: '+name);
assert(page.includes('quick-demo-multi { grid-column:1 / -1'),
  'Full-width multiomics call to action was lost');
assert(page.includes('exact_native_r_reference_table')||
       source('src/lib/multiomics/research-export-package.js').includes('exact_native_r_reference_table'),
       'Complete native R source table provenance missing');
assert(source('src/lib/multiomics/research-export-package.js').includes('ANALYSE_REPRODUCTIBLE.Rmd'),
  'Downloaded research ZIP must contain the annotated Rmd');
assert(source('src/routes/multiomics/learn/+page.svelte').includes('DIABLO'),
  'Beginner method guide missing');
assert(source('src/routes/multiomics/demo/+page.svelte').includes('8 sujets') &&
       source('src/routes/multiomics/demo/+page.svelte').includes('?demo=multi'),
  'Guided practical demo missing');
console.log('PASS: isolated repository migration static completeness and quick-demo UX');
console.log('NOTE: this does NOT replace npm build, browser tests or independent scientific R validation');
