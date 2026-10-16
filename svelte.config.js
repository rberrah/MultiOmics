import adapter from '@sveltejs/adapter-static';

/** @type {import('@sveltejs/kit').Config} */
const config = {
  kit: {
    adapter: adapter({fallback:'404.html'}),
    paths: {base:process.env.BASE_PATH||''},
    prerender:{entries:['*','/multiomics/','/multiomics/tool/','/multiomics/learn/'],handleHttpError:'fail'}
  }
};
export default config;
