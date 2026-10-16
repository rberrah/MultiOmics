import adapter from '@sveltejs/adapter-static';

/** @type {import('@sveltejs/kit').Config} */
const config = {
  trailingSlash: 'always',
  kit: {
    // Separate from PMx Explain's /_app assets when both builds share one domain.
    appDir: '_multiomics_app',
    adapter: adapter({fallback:'404.html'}),
    paths: {base:process.env.BASE_PATH||''},
    prerender:{entries:['*','/multiomics/','/multiomics/tool/','/multiomics/learn/'],handleHttpError:'fail'}
  }
};
export default config;
