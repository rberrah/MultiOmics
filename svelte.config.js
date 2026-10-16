import adapter from '@sveltejs/adapter-static';

/** @type {import('@sveltejs/kit').Config} */
const config = {
  trailingSlash: 'always',
  kit: {
    adapter: adapter({ fallback: '404.html' }),
    paths: { base: process.env.BASE_PATH || '' },
    prerender: { entries: ['*'] }
  }
};
export default config;
