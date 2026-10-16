// Portable, read-only local server for the downloaded MultiOmics static website.
// Uses only the Node.js standard library; no npm install and no R needed.
import { createServer } from 'node:http';
import { readFile, lstat, realpath } from 'node:fs/promises';
import { dirname, extname, join, resolve, sep } from 'node:path';
import { fileURLToPath } from 'node:url';

const directory = dirname(fileURLToPath(import.meta.url));
const site = resolve(directory, 'site');
const siteReal = await realpath(site);
const args = process.argv.slice(2);
const idx = args.indexOf('--port');
if (idx >= 0 && (idx === args.length - 1 || !/^[0-9]+$/.test(args[idx + 1])))
  throw new Error('Usage: node serve.mjs [--port 4173]');
const port = idx >= 0 ? Number(args[idx + 1]) : Number(process.env.PORT || 4173);
if (!Number.isInteger(port) || port < 1 || port > 65535)
  throw new Error('Invalid TCP port');
const host = process.env.MULTIOMICS_HOST || '127.0.0.1';
if (!['127.0.0.1', 'localhost', '::1'].includes(host) &&
    process.env.MULTIOMICS_ALLOW_REMOTE !== 'yes') {
  throw new Error('Refusing to expose a research interface outside localhost. ' +
    'Remote serving requires explicit MULTIOMICS_ALLOW_REMOTE=yes and separate security review.');
}
const mime = {
  '.html':'text/html; charset=utf-8',
  '.js':'text/javascript; charset=utf-8',
  '.css':'text/css; charset=utf-8',
  '.json':'application/json',
  '.csv':'text/csv; charset=utf-8',
  '.txt':'text/plain; charset=utf-8',
  '.md':'text/markdown; charset=utf-8',
  '.svg':'image/svg+xml',
  '.woff':'font/woff',
  '.woff2':'font/woff2',
  '.ttf':'font/ttf',
  '.otf':'font/otf',
  '.png':'image/png',
  '.jpg':'image/jpeg',
  '.jpeg':'image/jpeg',
  '.webp':'image/webp',
  '.ico':'image/x-icon',
  '.webmanifest':'application/manifest+json',
  '.xml':'application/xml',
  '.wasm':'application/wasm'
};
const safeFile = (file) => file === site || file.startsWith(site + sep);
createServer(async (req, res) => {
  res.setHeader('Cache-Control', 'no-store');
  res.setHeader('X-Content-Type-Options', 'nosniff');
  if (!['GET', 'HEAD'].includes(req.method)) {
    res.writeHead(405, { Allow: 'GET, HEAD' }).end();
    return;
  }
  let pathname;
  try { pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname); }
  catch { res.writeHead(400).end('Invalid URL'); return; }
  if (pathname.includes('\0') || pathname.includes('\\')) {
    res.writeHead(400).end('Invalid path'); return;
  }
  let file = resolve(site, '.' + pathname);
  if (!safeFile(file)) { res.writeHead(403).end('Forbidden'); return; }
  try {
    const stat = await lstat(file);
    if (stat.isSymbolicLink()) { res.writeHead(403).end('Forbidden'); return; }
    if (stat.isDirectory()) file = join(file, 'index.html');
    const actual = await realpath(file);
    if (!safeFile(actual)) { res.writeHead(403).end('Forbidden'); return; }
    if (!(await lstat(actual)).isFile()) throw new Error('Not a file');
    const bytes = await readFile(actual);
    res.writeHead(200, { 'Content-Type': mime[extname(actual).toLowerCase()] || 'application/octet-stream' });
    res.end(req.method === 'HEAD' ? undefined : bytes);
  } catch (error) {
    if (error?.code === 'ENOENT' || error?.code === 'ENOTDIR' || error?.message === 'Not a file') {
      res.writeHead(404, { 'Content-Type': 'text/plain; charset=utf-8' }).end('Not found');
    } else {
      console.error('Local server error:', error.message);
      res.writeHead(500).end('Internal server error');
    }
  }
}).listen(port, host, () => {
  console.log('MultiOmics portable: http://' + host + ':' + port + '/multiomics/tool/');
  console.log('Browser exploratory analysis only; this server is not an R backend.');
});
