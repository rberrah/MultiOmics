import {defineConfig,devices} from '@playwright/test';
const port=4173;
const url=process.env.LABS_E2E_URL||'http://127.0.0.1:'+port;
export default defineConfig({
 testDir:'./tests/e2e',
 testMatch:['multiomics*.spec.js'],
 timeout:45000,
 fullyParallel:false,
 retries:0,
 workers:1,
 reporter:process.env.CI?[['list'],['github']]:'list',
 use:{baseURL:url,trace:'off',...devices['Desktop Chrome']},
 projects:[{name:'chromium',use:{...devices['Desktop Chrome'],trace:'off'}}],
 webServer:process.env.LABS_E2E_URL?undefined:{
  command:process.env.LABS_E2E_PREBUILT
   ?'npm run preview -- --host 127.0.0.1 --port 4173'
   :'npm run build && npm run preview -- --host 127.0.0.1 --port 4173',
  url,timeout:180000,reuseExistingServer:!process.env.CI
 }
});
