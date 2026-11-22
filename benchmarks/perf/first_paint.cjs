// First paint at 100k cells: time from navigation to the first WebGL draw of
// points (instanced draw with at least 1,000 instances), in headless Chromium on
// the real GPU, against the real `lumencell launch` serving the built app.
//
//   node benchmarks/perf/first_paint.cjs <data.h5ad> <root> [port]
//
// Run it with the server's cache and state dirs (the ones `timing_100k.py` used)
// set in the environment, and with nothing else heavy running on the machine.
// The timing run must be finished first: this measurement would skew it.
//
// Uses web/node_modules/playwright (run from the repo root; the script resolves it).

const { createRequire } = require("node:module");
const { spawn } = require("node:child_process");
const path = require("node:path");

const requireFromWeb = createRequire(path.resolve("web/package.json"));
const { chromium } = requireFromWeb("playwright");

const [data, root, portArg] = process.argv.slice(2);
const port = Number(portArg ?? 8790);
const url = `http://127.0.0.1:${port}/`;
const BUDGET_MS = 1500;
const RUNS = 3;

async function waitForHealth() {
  for (let i = 0; i < 600; i++) {
    try {
      const r = await fetch(`${url}api/v1/healthz`);
      if (r.ok) return;
    } catch {
      /* not up yet */
    }
    await new Promise((r) => setTimeout(r, 500));
  }
  throw new Error("server did not come up");
}

// Records the first instanced draw of at least 1,000 instances, after the
// navigation starts, plus the first contentful paint.
const INIT = `
  window.__firstPoints = null;
  const mark = (inst) => {
    if (window.__firstPoints === null && inst >= 1000) window.__firstPoints = performance.now();
  };
  for (const ctx of [WebGL2RenderingContext, WebGLRenderingContext]) {
    const p = ctx.prototype;
    const a = p.drawArraysInstanced, e = p.drawElementsInstanced;
    if (a) p.drawArraysInstanced = function (m, f, c, n) { mark(n); return a.call(this, m, f, c, n); };
    if (e) p.drawElementsInstanced = function (m, c, t, o, n) { mark(n); return e.call(this, m, c, t, o, n); };
  }
  window.__fcp = null;
  new PerformanceObserver((list) => {
    for (const en of list.getEntriesByName("first-contentful-paint")) window.__fcp = en.startTime;
  }).observe({ type: "paint", buffered: true });
`;

async function main() {
  const server = spawn(
    "lumencell",
    ["launch", data, "--root", root, "--no-browser", "--strict-port", "--port", String(port)],
    { stdio: ["ignore", "ignore", "inherit"] },
  );
  try {
    await waitForHealth();
    const browser = await chromium.launch({
      headless: true,
      args: ["--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"],
    });
    const results = [];
    for (let i = 0; i < RUNS; i++) {
      const context = await browser.newContext({ viewport: { width: 1440, height: 900 } });
      await context.addInitScript(INIT);
      const page = await context.newPage();
      const t0 = Date.now();
      await page.goto(url, { waitUntil: "commit" });
      await page.waitForFunction(() => window.__firstPoints !== null, null, { timeout: 180_000 });
      const m = await page.evaluate(() => ({ points: window.__firstPoints, fcp: window.__fcp }));
      results.push({ run: i + 1, first_points_ms: Math.round(m.points), fcp_ms: m.fcp && Math.round(m.fcp), wall_ms: Date.now() - t0 });
      await context.close();
    }
    await browser.close();
    const sorted = results.map((r) => r.first_points_ms).sort((a, b) => a - b);
    const median = sorted[Math.floor(sorted.length / 2)];
    console.log(JSON.stringify({ budget_ms: BUDGET_MS, median_first_points_ms: median, verdict: median < BUDGET_MS ? "pass" : "fail", runs: results }, null, 2));
  } finally {
    server.kill("SIGINT");
  }
}

main().catch((err) => {
  console.error(err);
  process.exitCode = 1;
});
