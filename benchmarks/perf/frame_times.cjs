// Frame times while the plot moves: the 2-D <-> 3-D morph, auto-rotate and an orbit drag, with 3-D
// at rest as the baseline. Headless Chromium on the real GPU (the renderer is printed, so a software
// fallback shows), against the real `lumencell launch` serving the built app, at 1440x900 and device
// pixel ratio 2 (a Retina laptop). A frame time is the gap between two requestAnimationFrame
// callbacks. While the plot moves, only the frames that drew the cells (`drawn`) are timed; at rest,
// all of them. `main_ms_per_frame` is the page's main-thread busy time (Chrome's TaskDuration) per
// frame, Playwright's own waits included, so it overstates: near the frame time, the main thread is
// the limit; far below it, the GPU is. Swap in use is recorded, since paging skews both.
//
//   node benchmarks/perf/frame_times.cjs <data.h5ad> <root> [port]
//
// Run it like first_paint.cjs: from the repo root, with LUMENCELL_CACHE_DIR and LUMENCELL_DATA_DIR
// set to the dirs the dataset was prepared with, and nothing else heavy running. It waits for the
// run that --root starts, so the trajectory curves are drawn too.

const { createRequire } = require("node:module");
const { execSync, spawn } = require("node:child_process");
const path = require("node:path");
const readline = require("node:readline");

const requireFromWeb = createRequire(path.resolve("web/package.json"));
const { chromium } = requireFromWeb("playwright");

const [data, root, portArg] = process.argv.slice(2);
if (!data || !root) {
  console.error("usage: node benchmarks/perf/frame_times.cjs <data.h5ad> <root> [port]");
  process.exit(2);
}
const port = Number(portArg ?? 8791);
const url = `http://127.0.0.1:${port}/`;
const REPEATS = 3;
const FRAME_MS = 1000 / 60;
// A gap this long missed at least one 60 Hz vsync.
const MISSED_MS = 1.5 * FRAME_MS;
const sleep = (ms) => new Promise((r) => setTimeout(r, ms));

async function getJson(p) {
  const r = await fetch(url + p);
  if (!r.ok) throw new Error(`${p}: HTTP ${r.status}`);
  return r.json();
}

async function waitForHealth(server) {
  for (let i = 0; i < 600; i++) {
    if (server.exitCode !== null) throw new Error(`the server exited with status ${server.exitCode}`);
    try {
      if ((await fetch(`${url}api/v1/healthz`)).ok) return;
    } catch {
      /* not up yet */
    }
    await sleep(500);
  }
  throw new Error("server did not come up");
}

// The run --root asked for has finished, and no run is still going (it would compete with the
// browser for the CPU). Run ids hash their inputs, so a relaunch names the run it already has.
async function waitForRun(server, said) {
  for (let i = 0; i < 3600; i++) {
    if (server.exitCode !== null) throw new Error(`the server exited with status ${server.exitCode}`);
    if (said.rootError) throw new Error(said.rootError);
    if (said.runId) {
      const { runs } = await getJson("api/v1/trajectories");
      const ours = runs.find((r) => r.run_id === said.runId);
      if (ours?.state === "failed" || ours?.state === "cancelled") throw new Error(`run ${ours.run_id} ${ours.state}: ${ours.error ?? ""}`);
      if (ours?.state === "done" && !runs.some((r) => r.state === "queued" || r.state === "running")) return;
    }
    await sleep(1000);
  }
  throw new Error("no run finished within an hour");
}

// Counts instanced draws of at least 1,000 instances (the cells), and records rAF gaps on demand,
// each with whether the cells were drawn in it.
const INIT = `
  window.__points = 0;
  for (const ctx of [WebGL2RenderingContext, WebGLRenderingContext]) {
    const p = ctx.prototype;
    const a = p.drawArraysInstanced, e = p.drawElementsInstanced;
    if (a) p.drawArraysInstanced = function (m, f, c, n) { if (n >= 1000) window.__points++; return a.call(this, m, f, c, n); };
    if (e) p.drawElementsInstanced = function (m, c, t, o, n) { if (n >= 1000) window.__points++; return e.call(this, m, c, t, o, n); };
  }
  window.__record = () => {
    const gaps = [], drew = [];
    let last = null, seen = window.__points, on = true;
    requestAnimationFrame(function tick(now) {
      if (last !== null) {
        gaps.push(now - last);
        drew.push(window.__points > seen);
      }
      last = now;
      seen = window.__points;
      if (on) requestAnimationFrame(tick);
    });
    return () => { on = false; return { gaps, drew }; };
  };
`;

function stats({ gaps, drew, mainMs }) {
  const drawn = drew.filter(Boolean).length;
  const s = (drawn ? gaps.filter((_, i) => drew[i]) : [...gaps]).sort((a, b) => a - b);
  const q = (p) => Math.round(s[Math.min(s.length - 1, Math.floor(p * s.length))] * 10) / 10;
  return {
    frames: gaps.length,
    drawn,
    p50_ms: q(0.5),
    p95_ms: q(0.95),
    max_ms: q(1),
    missed: s.filter((g) => g > MISSED_MS).length,
    main_ms_per_frame: Math.round((mainMs / gaps.length) * 10) / 10,
  };
}

async function main() {
  const server = spawn(
    process.env.LUMENCELL_BIN ?? "lumencell",
    ["launch", data, "--root", root, "--no-browser", "--strict-port", "--port", String(port)],
    { stdio: ["ignore", "pipe", "pipe"] },
  );
  // launch prints which run --root started ("running paga_dpt from ... (run <id>)"), or why it did not.
  const said = { runId: null, rootError: null };
  readline.createInterface({ input: server.stdout }).on("line", (line) => {
    said.runId ??= /\(run ([^)\s]+)\)$/.exec(line)?.[1] ?? null;
  });
  readline.createInterface({ input: server.stderr }).on("line", (line) => {
    console.error(line);
    if (/--root was not applied|could not start/.test(line)) said.rootError = line;
  });
  let browser = null;
  try {
    await waitForHealth(server);
    await waitForRun(server, said);
    browser = await chromium.launch({
      headless: true,
      args: ["--use-angle=metal", "--enable-gpu", "--ignore-gpu-blocklist"],
    });
    const context = await browser.newContext({
      viewport: { width: 1440, height: 900 },
      deviceScaleFactor: 2,
      reducedMotion: "no-preference",
    });
    await context.addInitScript(INIT);
    const page = await context.newPage();
    await page.goto(url);
    const view = (name) => page.getByRole("group", { name: "View" }).getByRole("button", { name, exact: true });
    const autoRotate = page.getByRole("button", { name: "Auto-rotate", exact: true });
    // A click that did not take would time a still plot.
    const expectPressed = async (button, want) => {
      if ((await button.getAttribute("aria-pressed")) !== String(want)) throw new Error("a click on the plot's controls did not take");
    };
    // The run's curves are drawn (2-D caption), and 3-D has opened.
    await page.getByText("curves are projections of the tree").waitFor({ timeout: 600_000 });
    await page.waitForFunction(
      () => [...document.querySelectorAll('[aria-label="View"] button')].some((b) => b.textContent === "3-D" && !b.disabled),
      null,
      { timeout: 600_000 },
    );
    const gpu = await page.evaluate(() => {
      const gl = document.createElement("canvas").getContext("webgl2");
      const ext = gl.getExtension("WEBGL_debug_renderer_info");
      return ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER);
    });
    const { n_obs } = await getJson("api/v1/dataset");

    const cdp = await context.newCDPSession(page);
    await cdp.send("Performance.enable");
    const busy = async () => (await cdp.send("Performance.getMetrics")).metrics.find((m) => m.name === "TaskDuration").value * 1000;
    const measure = async (act) => {
      await page.evaluate(() => (window.__stop = window.__record()));
      const before = await busy();
      await act();
      const r = await page.evaluate(() => window.__stop());
      return { ...r, mainMs: (await busy()) - before };
    };
    const settleIn = async (name) => {
      await view(name).click();
      await expectPressed(view(name), true);
      await sleep(1000);
    };
    // The morph lasts 0.6 s; 0.7 s covers it and little rest.
    const morphTo = (name) => async () => {
      await view(name).click();
      await expectPressed(view(name), true);
      await sleep(700);
    };
    const scenarios = {
      rest_3d: async () => {
        await settleIn("3-D");
        return measure(() => sleep(2000));
      },
      morph_to_3d: async () => {
        await settleIn("2-D");
        return measure(morphTo("3-D"));
      },
      morph_to_2d: async () => {
        await settleIn("3-D");
        return measure(morphTo("2-D"));
      },
      auto_rotate: async () => {
        await settleIn("3-D");
        const r = await measure(async () => {
          await autoRotate.click();
          await expectPressed(autoRotate, true);
          await sleep(3000);
        });
        await autoRotate.click();
        await expectPressed(autoRotate, false);
        return r;
      },
      // 120 mouse moves, inside the stage (left of where the evidence sheet opens). Playwright delivers
      // about one move every two frames, so about half the frames draw.
      orbit_drag: async () => {
        await settleIn("3-D");
        if ((await page.evaluate(() => document.elementFromPoint(500, 400)?.tagName)) !== "CANVAS")
          throw new Error("orbit_drag: the drag would not start on the plot");
        const before = await page.evaluate(() => location.hash);
        const r = await measure(async () => {
          await page.mouse.move(500, 400);
          await page.mouse.down();
          for (let k = 1; k <= 120; k++) {
            await page.mouse.move(500 + 2 * k, 400 + k / 2);
            await sleep(16);
          }
          await page.mouse.up();
        });
        await sleep(400); // the view reaches the URL at most every 300 ms
        if ((await page.evaluate(() => location.hash)) === before) throw new Error("orbit_drag: the drag did not move the camera");
        return r;
      },
    };
    // A moving scenario in which the cells drew on few frames measured mostly a still plot.
    const run = async (name) => {
      const r = await scenarios[name]();
      if (name !== "rest_3d" && r.drew.filter(Boolean).length < 0.3 * r.gaps.length)
        throw new Error(`${name}: the cells drew on too few frames to time`);
      return r;
    };
    // One untimed pass compiles the shaders and loads the 3-D curves.
    for (const name of Object.keys(scenarios)) await run(name);
    const results = {};
    for (const name of Object.keys(scenarios)) {
      const all = { gaps: [], drew: [], mainMs: 0 };
      for (let i = 0; i < REPEATS; i++) {
        const r = await run(name);
        all.gaps.push(...r.gaps);
        all.drew.push(...r.drew);
        all.mainMs += r.mainMs;
      }
      results[name] = stats(all);
    }
    console.log(
      JSON.stringify(
        {
          date: new Date().toISOString().slice(0, 10),
          n_cells: n_obs,
          gpu,
          software_gpu: /swiftshader|llvmpipe|software/i.test(gpu),
          swap: process.platform === "darwin" ? execSync("sysctl -n vm.swapusage").toString().trim() : null,
          viewport: "1440x900 @2x",
          repeats: REPEATS,
          missed_means_gap_over_ms: Math.round(MISSED_MS * 10) / 10,
          scenarios: results,
        },
        null,
        2,
      ),
    );
  } finally {
    await browser?.close();
    server.kill("SIGINT");
  }
}

main().catch((err) => {
  console.error(err);
  process.exitCode = 1;
});
