/**
 * Read-only, SHA-256-verified handoff for https://rberrah.github.io/multiomics/.
 * Build with BASE_PATH unset: SvelteKit assets then resolve from root
 * /_multiomics_app/ rather than /pharmacometrie/_multiomics_app/.
 * This does not modify or publish PMx Explain.
 */
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { cp, lstat, mkdir, readdir, readFile, rm, writeFile } from 'node:fs/promises';
import { join, relative, resolve, sep } from 'node:path';

const root=resolve('.');
const build=join(root,'build');
const target=join(root,'short-url-overlay');
const site=join(target,'site');
const dirs=['multiomics','_multiomics_app'];
const paths=['multiomics/index.html','multiomics/tool/index.html',
  'multiomics/learn/index.html','multiomics/disciplines/index.html',
  'multiomics/resources/index.html'];
for(const p of paths){
  assert((await lstat(join(build,p))).isFile(),'Expected standalone route missing: '+p);
  const h=await readFile(join(build,p),'utf8');
  assert(h.includes('_multiomics_app/'),'Root bundle reference missing: '+p);
  assert(!h.includes('/pharmacometrie/_multiomics_app/'),
    'Refusing to publish legacy prefixed build to root short URL: '+p);
  assert(!h.includes('/MultiOmics/_multiomics_app/'),
    'Refusing project Pages-prefixed build to root short URL: '+p);
}
for(const dir of dirs)
  assert((await lstat(join(build,dir))).isDirectory(),'Missing '+dir);
await rm(target,{recursive:true,force:true});
await mkdir(site,{recursive:true});
for(const dir of dirs)
  await cp(join(build,dir),join(site,dir),{recursive:true});
const files=[];
async function walk(dir){
  for(const name of (await readdir(dir)).sort()){
    const p=join(dir,name);
    const s=await lstat(p);
    assert(!s.isSymbolicLink(),'Symlink forbidden in direct-URL payload: '+p);
    if(s.isDirectory()){await walk(p);continue;}
    assert(s.isFile(),'Unexpected direct-URL payload object');
    const path=relative(site,p).split(sep).join('/');
    assert(/^(?:multiomics|_multiomics_app)\//.test(path),
      'Root overlay violates exact namespace allowlist: '+path);
    const bytes=await readFile(p);
    files.push({path,bytes:bytes.length,
      sha256:createHash('sha256').update(bytes).digest('hex')});
  }
}
await walk(site);
assert(files.length>15,'Direct-URL static payload unexpectedly small');
await writeFile(join(target,'manifest.json'),JSON.stringify({
  schema:'multiomics.root-short-url-overlay.v1',
  source:'rberrah/MultiOmics',basePath:'',
  deployment:'NOT_DEPLOYED_BY_SOURCE_REPOSITORY',
  directories:dirs,
  routes:paths.map(p=>'/'+p.replace(/index\.html$/,'')),
  files
},null,2)+'\n');
console.log('PASS direct /multiomics/ root overlay:',files.length,
  'files, only /multiomics/ + /_multiomics_app/. No PMx code modified.');
