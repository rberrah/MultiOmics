import assert from 'node:assert/strict';
import {createHash} from 'node:crypto';
import {lstat,readFile,readdir} from 'node:fs/promises';
import {resolve,join,relative,sep} from 'node:path';
const root=resolve('short-url-overlay');
const x=JSON.parse(await readFile(join(root,'manifest.json'),'utf8'));
assert.equal(x.schema,'multiomics.root-short-url-overlay.v1');
assert.equal(x.basePath,'');
assert.deepEqual(x.directories,['multiomics','_multiomics_app']);
assert.deepEqual(x.routes,['/multiomics/','/multiomics/tool/','/multiomics/learn/',
  '/multiomics/ecosystem/','/multiomics/disciplines/','/multiomics/resources/']);
const paths=new Set();
for(const f of x.files){
  assert(/^(?:multiomics|_multiomics_app)\//.test(f.path));
  assert(!f.path.includes('..')&&!f.path.includes('\\'));
  assert(!paths.has(f.path));paths.add(f.path);
  const bytes=await readFile(join(root,'site',f.path));
  assert.equal(bytes.length,f.bytes);
  assert.equal(createHash('sha256').update(bytes).digest('hex'),f.sha256);
}
for(const p of ['multiomics/index.html','multiomics/tool/index.html',
'multiomics/learn/index.html','multiomics/ecosystem/index.html','multiomics/disciplines/index.html',
'multiomics/resources/index.html']){
  assert(paths.has(p));
  const h=await readFile(join(root,'site',p),'utf8');
  assert(h.includes('_multiomics_app/'));
  assert(!h.includes('/pharmacometrie/_multiomics_app/'));
}
assert([...paths].some(p=>p.startsWith('_multiomics_app/')));
async function walk(dir){
  for(const entry of await readdir(dir)){
    const p=join(dir,entry);
    const s=await lstat(p);
    assert(!s.isSymbolicLink());
    if(s.isDirectory()){await walk(p);continue;}
    const rel=relative(join(root,'site'),p).split(sep).join('/');
    assert(paths.has(rel),'Unmanifested file: '+rel);
  }
}
await walk(join(root,'site'));
console.log('PASS direct /multiomics/: no redirect, namespaced assets and sha256-verified files');
