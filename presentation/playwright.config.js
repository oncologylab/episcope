const { defineConfig } = require("@playwright/test");

module.exports = defineConfig({
  testDir: "./tests",
  testMatch: "presentation.spec.js",
  timeout: 60000,
  expect: { timeout: 15000 },
  use: {
    baseURL: "http://127.0.0.1:4173/presentation/reveal/",
    viewport: { width: 1920, height: 1080 },
    trace: "retain-on-failure"
  },
  webServer: {
    command: "python3 -m http.server 4173 --bind 127.0.0.1 --directory ..",
    cwd: __dirname,
    url: "http://127.0.0.1:4173/presentation/reveal/",
    reuseExistingServer: true
  }
});
