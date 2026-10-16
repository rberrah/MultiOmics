// Validate a portable browser package independently from its source build.
import assert from 'node:assert/strict';
import { createHash } from 'node:crypto';
import { spawn } from 'node:child_process';
import { lstat, readFile, readdir } from 'node:fs/promises';
import { resolve, join, relative, sep } from 'node:path';

const root=resolve('release/multiomics-web');
const manifest=JSON.parse(await readFile(join(root,'manifest.sha256.json'),'utf8'));
assert.equal(manifest.schema,'multiomics.portable-web.v1');
assert.equal(manifest.basePath,'');
assert.equal(manifest.includedBackend,'browser-only');
assert(manifest.files.length>15);
const tracked=new Set();
for (const file of manifest.files) {
  assert(typeof file.path==='string' && !file.path.startsWith('/') &&
    !file.path.includes('\\') &&
    !file.path.split('/').includes('..') &&
    !tracked.has(file.path));
  tracked.add(file.path);
  const bytes=await readFile(join(root,file.path));
  assert.equal(bytes.length,file.bytes,'Incorrect size for '+file.path);
  assert.equal(createHash('sha256').update(bytes).digest('hex'),file.sha256,
    'Modified portable file: '+file.path);
}
for(const required of ['serve.mjs','START_HERE.txt','site/index.html',
 'site/multiomics/index.html','site/multiomics/tool/index.html',
 'site/multiomics/learn/index.html','LICENSE','LICENSE-CONTENT.md'])
  assert(tracked.has(required),'Missing portable file '+required);
async function scan(folder) {
  for(const name of await readdir(folder)) {
    const p=join(folder,name);
    const st=await lstat(p);
    assert(!st.isSymbolicLink());
    if(st.isDirectory()){await scan(p);continue;}
    const rel=relative(root,p).split(sep).join('/');
    assert(rel==='manifest.sha256.json'||tracked.has(rel),'Unmanifested packaged file: '+rel);
  }
}
await scan(root);

const port=4187, url='http://127.0.0.1:'+port;
const child=spawn(process.execPath,[join(root,'serve.mjs'),'--port',String(port)],
  {cwd:root,env:{...process.env,PORT:String(port)},stdio:['ignore','pipe','pipe']});
let output='';
child.stdout.on('data',x=>{output+=x.toString()});
child.stderr.on('data',x=>{output+=x.toString()});
async function wait() {
  for(let i=0;i<80;i++){
    if(child.exitCode!==null)throw Error('Portable server crashed: '+output);
    try{
      const r=await fetch(url+'/multiomics/tool/');
      if(r.status===200)return;
    }catch{}
    await new Promise(done=>setTimeout(done,100));
  }
  throw Error('Portable server did not become ready: '+output);
}
try {
  await wait();
  for(const path of ['/','/multiomics/','/multiomics/tool/','/multiomics/learn/',
     '/multiomics/README.txt']) {
    const r=await fetch(url+path);
    assert.equal(r.status,200,path+' did not load');
    const body=await r.text();
    assert(body.length>0,'Empty '+path);
  }
  const page=await (await fetch(url+'/multiomics/tool/')).text();
  assert(page.includes('_multiomics_app/'),'Portable JS references missing');
  const notFound=await fetch(url+'/nonexistent-path');
  assert.equal(notFound.status,404,'Nonexistent path should fail instead of returning the app shell');
  const post=await fetch(url+'/multiomics/tool/',{method:'POST'});
  assert.equal(post.status,405,'Portable server must be read-only');
  console.log('PASS: portable ZIP-ready package manifest, routes, JS assets, 404 and read-only localhost server');
} finally {
  child.kill('SIGTERM');
}
