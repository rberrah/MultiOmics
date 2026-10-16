import { defineConfig, devices } from '@playwright/test';
export default defineConfig({
  testDir: './tests/standalone',
  timeout: 45000,
  use: { baseURL: 'http://127.0.0.1:4173', trace: 'retain-on-failure' },
  projects: [ { name: 'desktop', use: { ...devices['Desktop Chrome'] } }, { name: 'mobile', use: { ...devices['iPhone 13'], browserName: 'chromium' } } ],
  webServer: { command: 'npm run dev -- --host 127.0.0.1 --port 4173', url: 'http://127.0.0.1:4173/multiomics/tool/', reuseExistingServer: !process.env.CI, timeout: 120000 }
});
