// @ts-nocheck
import adapter from '@sveltejs/adapter-static';
import {sveltekit} from '@sveltejs/kit/vite';
import {defineConfig} from 'vite';

// SvelteKit 3: project configuration is passed to the Vite plugin.
export default defineConfig({
  plugins:[sveltekit({
    appDir:'_multiomics_app',
    adapter:adapter({fallback:'404.html'}),
    paths:{base:process.env.BASE_PATH||''},
    prerender:{entries:['*','/multiomics/','/multiomics/tool/','/multiomics/learn/'],handleHttpError:'fail'}
  })]
});
