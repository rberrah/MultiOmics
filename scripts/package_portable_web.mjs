// Build a shareable local-only static MultiOmics website from npm run build.
// The package contains rendered browser assets and a zero-dependency Node server.
// It does not claim to bundle Bioconductor or a hosted R service.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { cp, lstat, mkdir, readFile, readdir, rm, writeFile } from 'node:fs/promises';
import { resolve, join, relative, sep } from 'node:path';

const project = resolve('.');
const built = join(project, 'build');
const target = join(project, 'release', 'multiomics-web');
const site = join(target, 'site');
const required = ['index.html', 'multiomics/index.html',
  'multiomics/tool/index.html', 'multiomics/learn/index.html', 'multiomics/demo/index.html'];
for (const name of required) {
  const st = await lstat(join(built, name)).catch(() => null);
  assert(st?.isFile(), 'Missing portable site page: ' + name + '. Run npm run build with BASE_PATH unset.');
}
const homepage = await readFile(join(built, 'index.html'), 'utf8');
const tool = await readFile(join(built, 'multiomics/tool/index.html'), 'utf8');
assert(!homepage.includes('/pharmacometrie/_multiomics_app/') &&
       !tool.includes('/pharmacometrie/_multiomics_app/'),
       'This is a PMx-prefixed build, not a standalone one. Unset BASE_PATH and rebuild.');
assert(tool.includes('_multiomics_app/'), 'The compiled app bundle is missing.');
await rm(target, { recursive:true, force:true });
await mkdir(site, { recursive:true });
await cp(built, site, { recursive:true });
await cp(join(project,'scripts','portable_server.mjs'), join(target,'serve.mjs'));
for (const name of ['LICENSE', 'LICENSE-CONTENT.md', 'LICENSING.md', 'NOTICE', 'CITATION.cff'])
  await cp(join(project,name), join(target,name));
await writeFile(join(target, 'START_HERE.txt'),
'MultiOmics — standalone portable web edition\n' +
'=============================================\n' +
'1. Install Node.js 22 or newer if it is not installed.\n' +
'2. From this directory run: node serve.mjs\n' +
'3. Open: http://127.0.0.1:4173/multiomics/tool/\n' +
'4. Stop with Ctrl+C. No npm install is necessary to run this edition.\n\n' +
'This is a local, browser-first single/multi-omics exploratory tool.\n' +
'The package does NOT contain the reference R/Bioconductor backend.\n' +
'For reference R inference or to modify the code, download the full source\n' +
'from https://github.com/rberrah/MultiOmics and follow README.md.\n\n' +
'This static bundle does NOT transmit datasets to a backend by default;\n' +
'optional user-selected third-party lookups or external R backends need\n' +
'a separate privacy and security assessment. Do not treat its exploratory\n' +
'results as validated clinical or confirmatory inference.\n\n' +
'The canonical public URLs remain under\n' +
'https://rberrah.github.io/multiomics/\n' +
'Do not upload this root build over the PMx Explain deployment.\n' +
'The PMx site uses a separate, namespaced overlay built from the same source.\n\n' +
'Original source and licenses: https://github.com/rberrah/MultiOmics\n' +
'Legal notices and reuse conditions: read LICENSE, LICENSE-CONTENT.md,\n' +
'LICENSING.md and NOTICE in this folder. Noncommercial terms apply to\n' +
'newly licensed original code/content; earlier grants retain their terms.\n', 'utf8');

const manifest=[];
async function walk(dir) {
  for (const name of (await readdir(dir)).sort()) {
    const file=join(dir,name);
    const st=await lstat(file);
    assert(!st.isSymbolicLink(), 'Symlinks are not allowed in downloadable releases');
    if (st.isDirectory()) { await walk(file); continue; }
    assert(st.isFile(), 'Unexpected file kind in portable package');
    const pathname=relative(target,file).split(sep).join('/');
    const bytes=await readFile(file);
    manifest.push({path:pathname,bytes:bytes.length,
      sha256:createHash('sha256').update(bytes).digest('hex')});
  }
}
await walk(target);
manifest.sort((a,b)=>a.path.localeCompare(b.path,'en'));
assert(manifest.length>15,'Portable site is unexpectedly small');
await writeFile(join(target,'manifest.sha256.json'),JSON.stringify({
  schema:'multiomics.portable-web.v1',source:'rberrah/MultiOmics',
  basePath:'',includedBackend:'browser-only',files:manifest
},null,2)+'\n');
console.log('PASS: downloadable standalone browser package prepared:',manifest.length,'verified files');
console.log('Directory: release/multiomics-web; Node runtime only, no npm install after download.');
