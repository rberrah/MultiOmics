import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import { execFileSync } from 'node:child_process';

const root = process.argv[2] || 'tmp/public-benchmarks';
const reportPath = path.join(root, 'benchmark-report.json');

const sources = [
  ['mixOmics', 'mixOmics-org/mixOmics', 'external/mixOmics'],
  ['NCI60_GeneMetabolite_Data', 'Mathelab/NCI60_GeneMetabolite_Data', 'external/NCI60_GeneMetabolite_Data'],
  ['IntLIMVignettes', 'ncats/IntLIMVignettes', 'external/IntLIMVignettes'],
  ['xOmicsShiny', 'GenomicsNX/xOmicsShiny_app', 'external/xOmicsShiny'],
  ['PaintOmics', 'ConesaLab/PaintOmics', 'external/PaintOmics'],
  ['missRows', 'bioc/missRows', 'external/missRows']
];

function gitHead(dir) {
  try {
    return execFileSync('git', ['-C', dir, 'rev-parse', 'HEAD'], { encoding: 'utf8' }).trim();
  } catch {
    return null;
  }
}

async function walk(dir) {
  const entries = await fs.readdir(dir, { withFileTypes: true });
  const files = [];
  for (const entry of entries) {
    const full = path.join(dir, entry.name);
    if (entry.isDirectory()) files.push(...await walk(full));
    else files.push(full);
  }
  return files;
}

async function sha256File(file) {
  const buffer = await fs.readFile(file);
  return {
    bytes: buffer.byteLength,
    sha256: crypto.createHash('sha256').update(buffer).digest('hex')
  };
}

const report = JSON.parse(await fs.readFile(reportPath, 'utf8'));
const exportedFiles = (await walk(root))
  .filter((file) => !file.endsWith('benchmark-report.json') && !file.endsWith('benchmark-provenance.json'))
  .sort();

const fileManifest = {};
for (const file of exportedFiles) {
  const rel = path.relative(root, file).replaceAll('\\', '/');
  fileManifest[rel] = await sha256File(file);
}

let rVersion = null;
try {
  rVersion = execFileSync('Rscript', ['--version'], { encoding: 'utf8', stdio: ['ignore', 'pipe', 'pipe'] }).trim();
} catch {
  try {
    rVersion = execFileSync('R', ['--version'], { encoding: 'utf8' }).split('\n')[0].trim();
  } catch {
    rVersion = null;
  }
}

const provenance = {
  schema: 'pmx-explain-public-multiomics-benchmark-provenance/v1',
  generatedAt: new Date().toISOString(),
  pmxExplain: {
    repository: process.env.GITHUB_REPOSITORY || 'rberrah/rberrah.github.io',
    commit: process.env.GITHUB_SHA || gitHead('.'),
    workflowRunId: process.env.GITHUB_RUN_ID || null,
    workflowRunAttempt: process.env.GITHUB_RUN_ATTEMPT || null
  },
  runtime: {
    node: process.version,
    platform: `${process.platform}-${process.arch}`,
    r: rVersion
  },
  sourceRepositories: Object.fromEntries(sources.map(([key, repository, dir]) => [key, {
    repository,
    commit: gitHead(dir)
  }])),
  exportedInputFiles: fileManifest,
  policy: {
    immutableSourceRefs: true,
    exportedFilesCryptographicallyHashed: true,
    hashAlgorithm: 'SHA-256',
    note: 'A benchmark result is reproducible only when source repository commits, PMx Explain commit and exported input hashes are retained together.'
  }
};

report.sourceProvenance = provenance;
await fs.writeFile(reportPath, JSON.stringify(report, null, 2) + '\n');
await fs.writeFile(path.join(root, 'benchmark-provenance.json'), JSON.stringify(provenance, null, 2) + '\n');

console.log(`Attached provenance for ${Object.keys(provenance.sourceRepositories).length} public sources and ${Object.keys(fileManifest).length} exported files.`);
