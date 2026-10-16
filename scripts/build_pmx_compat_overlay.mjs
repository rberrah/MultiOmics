/**
 * Prepare a reproducible, allowlisted overlay for the existing PMx Explain URL.
 * This script does not publish or modify rberrah.github.io.
 * Build with BASE_PATH=/pharmacometrie before running.
 */
import { createHash } from 'node:crypto';
import { cp, lstat, mkdir, readFile, readdir, rm, writeFile } from 'node:fs/promises';
import { resolve, join, relative, sep } from 'node:path';

const project = resolve('.');
const built = join(project, 'build');
const output = join(project, 'compat-overlay');
const site = join(output, 'site', 'pharmacometrie');
const allowed = ['multiomics', '_multiomics_app'];
const required = [
  'multiomics/index.html',
  'multiomics/tool/index.html',
  'multiomics/learn/index.html'
];
const expectedBase = '/pharmacometrie';
function assert(condition, message) {
  if (!condition) throw new Error('PMx compatibility overlay: ' + message);
}
async function existsFile(path) {
  try { return (await lstat(path)).isFile(); } catch { return false; }
}
for (const path of required) {
  assert(await existsFile(join(built, path)), 'required prerendered page missing: ' + path);
  const html = await readFile(join(built, path), 'utf8');
  assert(html.includes(expectedBase + '/_multiomics_app/'),
    'page lacks namespaced standalone assets: ' + path);
  assert(!html.includes(expectedBase + '/_app/'),
    'PMx shared app assets would be overwritten: ' + path);
}
for (const dir of allowed) {
  const stat = await lstat(join(built, dir)).catch(() => null);
  assert(stat?.isDirectory(), 'missing isolated build directory: ' + dir);
}
await rm(output, { recursive: true, force: true });
await mkdir(site, { recursive: true });
for (const dir of allowed) await cp(join(built, dir), join(site, dir), { recursive: true });

const manifest = [];
async function walk(dir) {
  for (const name of (await readdir(dir)).sort()) {
    const file = join(dir, name);
    const stat = await lstat(file);
    assert(!stat.isSymbolicLink(), 'symlink forbidden in release payload: ' + file);
    if (stat.isDirectory()) { await walk(file); continue; }
    assert(stat.isFile(), 'unexpected non-file in release payload: ' + file);
    const path = relative(join(output, 'site'), file).split(sep).join('/');
    assert(path.startsWith('pharmacometrie/multiomics/') ||
      path.startsWith('pharmacometrie/_multiomics_app/'),
      'unapproved publication path: ' + path);
    const data = await readFile(file);
    manifest.push({ path, bytes: data.length,
      sha256: createHash('sha256').update(data).digest('hex') });
  }
}
await walk(site);
assert(manifest.length > required.length, 'release payload contains too few files');
const report = {
  schema: 'multiomics.compat-overlay.v1',
  source: 'rberrah/MultiOmics',
  urlBase: expectedBase,
  deployPolicy: 'APPEND_ONLY_NAMESPACED_PAYLOAD_NO_AUTOMATIC_PUBLISH',
  directories: allowed,
  routeUrls: required.map(path => expectedBase + '/' + path.replace(/index\.html$/, '')),
  files: manifest
};
await writeFile(join(output, 'manifest.json'), JSON.stringify(report, null, 2) + '\n');
console.log('Standalone MultiOmics PMx URL overlay prepared:', manifest.length, 'files.');
console.log('No PMx Explain code or public deployment has been modified.');
