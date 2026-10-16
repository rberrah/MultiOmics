// @ts-nocheck
import {sveltekit} from '@sveltejs/kit/vite';
import {defineConfig} from 'vite';
import {execFileSync} from 'node:child_process';

// Embed the exact checkout SHA in the exported research exchange manifest.
// If no .git metadata is present (e.g. source archive), mark it as unknown.
function sourceCommit() {
  try {
    const sha=execFileSync('git',['rev-parse','HEAD'],{encoding:'utf8',timeout:3000}).trim();
    return /^[a-f0-9]{40}$/.test(sha)?sha:'unknown';
  } catch {
    return 'unknown';
  }
}
const commit=process.env.VITE_MULTIOMICS_SOURCE_SHA||sourceCommit();
export default defineConfig({
  plugins:[sveltekit()],
  define:{'import.meta.env.VITE_MULTIOMICS_SOURCE_SHA':JSON.stringify(commit)}
});
