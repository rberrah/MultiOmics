import { defineConfig, devices } from '@playwright/test';

export default defineConfig({
  testDir:'./tests/e2e',
  testMatch:['multiomics-portable.spec.js'],
  timeout:45000,
  fullyParallel:false,
  workers:1,
  retries:0,
  reporter:process.env.CI?[['list'],['github']]:'list',
  use:{
    ...devices['Desktop Chrome'],
    baseURL:'http://127.0.0.1:4186',
    trace:'retain-on-failure'
  },
  webServer:{
    command:'node release/multiomics-web/serve.mjs --port 4186',
    url:'http://127.0.0.1:4186/multiomics/tool/',
    timeout:30000,
    reuseExistingServer:!process.env.CI
  }
});
