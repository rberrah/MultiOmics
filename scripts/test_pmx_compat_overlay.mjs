import { createHash } from 'node:crypto';
import { readFile, stat } from 'node:fs/promises';
import { join, resolve } from 'node:path';
import assert from 'node:assert/strict';

const output = resolve('compat-overlay');
const manifest = JSON.parse(await readFile(join(output, 'manifest.json'), 'utf8'));
assert.equal(manifest.schema, 'multiomics.compat-overlay.v1');
assert.equal(manifest.urlBase, '/pharmacometrie');
assert.deepEqual(manifest.directories, ['multiomics', '_multiomics_app']);
const seen = new Set();
for (const item of manifest.files) {
  assert.match(item.path, /^pharmacometrie\/(?:multiomics|_multiomics_app)\//);
  assert(!seen.has(item.path), 'duplicate path: ' + item.path);
  seen.add(item.path);
  const p = join(output, 'site', item.path);
  const bytes = await readFile(p);
  assert.equal((await stat(p)).size, item.bytes);
  assert.equal(createHash('sha256').update(bytes).digest('hex'), item.sha256);
}
for (const page of ['multiomics/index.html', 'multiomics/tool/index.html', 'multiomics/learn/index.html']) {
  const path = 'pharmacometrie/' + page;
  assert(seen.has(path), 'required route is absent: ' + path);
  const html = await readFile(join(output, 'site', path), 'utf8');
  assert.match(html, /_multiomics_app\//);
  assert(!html.includes('/pharmacometrie/_app/'), 'shared PMx bundle reference in ' + path);
}
assert(seen.size > 3);
console.log('PASS: checked ' + seen.size + ' files, immutable hashes, historical URLs and asset namespace.');
