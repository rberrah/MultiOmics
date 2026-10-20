import fs from 'node:fs';

const pkg=JSON.parse(fs.readFileSync('package.json','utf8'));
const lock=JSON.parse(fs.readFileSync('package-lock.json','utf8'));

if(lock.version!==pkg.version) throw new Error(`package-lock root version ${lock.version} != package.json ${pkg.version}`);
if(lock.packages?.['']?.version!==pkg.version) throw new Error(`package-lock packages[""].version ${lock.packages?.['']?.version} != package.json ${pkg.version}`);

const mismatches=[];
for(const [key,entry] of Object.entries(lock.packages||{})){
  if(!key.startsWith('node_modules/') || !entry || typeof entry!=='object') continue;
  const version=entry.version;
  const resolved=entry.resolved;
  if(typeof version!=='string' || typeof resolved!=='string') continue;
  let url;
  try{ url=new URL(resolved); }catch{ continue; }
  if(url.hostname!=='registry.npmjs.org') continue;
  const basename=decodeURIComponent(url.pathname.split('/').pop()||'');
  if(!basename.endsWith(`-${version}.tgz`)){
    mismatches.push({key,version,resolved,basename});
  }
}
if(mismatches.length){
  console.error('Package-lock metadata disagrees with resolved npm tarball versions:');
  for(const row of mismatches) console.error(JSON.stringify(row));
  process.exit(1);
}
console.log(`PASS: package-lock root version ${pkg.version} and ${Object.keys(lock.packages||{}).length-1} package entries are metadata-consistent`);
