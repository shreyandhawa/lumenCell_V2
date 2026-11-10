import { defineConfig, devices } from "@playwright/test";

// End-to-end: the real `lumencell launch` serving the built web app (web/dist),
// on a fresh synthetic dataset with a fresh cache and state, so every run goes
// through preprocessing, the root picker and a method run from scratch.
const PORT = 8799;
const DATA = ".lumencell/e2e";

export default defineConfig({
  testDir: "./e2e",
  timeout: 180_000,
  expect: { timeout: 60_000 },
  fullyParallel: false,
  workers: 1,
  reporter: [["list"]],
  use: {
    baseURL: `http://127.0.0.1:${PORT}`,
    trace: "retain-on-failure",
    viewport: { width: 1440, height: 900 },
  },
  projects: [{ name: "chromium", use: { ...devices["Desktop Chrome"], viewport: { width: 1440, height: 900 } } }],
  webServer: {
    cwd: "..",
    command: [
      `rm -rf ${DATA} && mkdir -p ${DATA}`,
      `uv run lumencell synth ${DATA}/e2e.h5ad --store-root --n-cells 1500 --n-genes 800`,
      `LUMENCELL_CACHE_DIR=${DATA}/cache LUMENCELL_DATA_DIR=${DATA}/state uv run lumencell launch ${DATA}/e2e.h5ad --no-browser --strict-port --port ${PORT}`,
    ].join(" && "),
    url: `http://127.0.0.1:${PORT}/api/v1/healthz`,
    timeout: 120_000,
    reuseExistingServer: false,
    stdout: "ignore",
    stderr: "pipe",
  },
});
