/// <reference lib="dom" />
import fs from "node:fs";
import { expect, type Page, test } from "@playwright/test";

test.use({ reducedMotion: "reduce" });

/** Collect console errors and uncaught exceptions; the test fails if there are any. */
function watchErrors(page: Page): string[] {
  const errors: string[] = [];
  page.on("console", (msg) => {
    if (msg.type() === "error") errors.push(msg.text());
  });
  page.on("pageerror", (err) => errors.push(err.message));
  return errors;
}

/** Share of the plot's pixels (between two fractions of its width) that are not the dark stage: 0 means nothing was drawn. */
function drawnShare(page: Page, from = 0, to = 1): Promise<number> {
  return page
    .locator("canvas")
    .first()
    .evaluate((gl: HTMLCanvasElement, [from, to]) => {
      const copy = document.createElement("canvas");
      copy.width = gl.width;
      copy.height = gl.height;
      const ctx = copy.getContext("2d")!;
      ctx.drawImage(gl, 0, 0);
      const x0 = Math.round(from * copy.width);
      const { data } = ctx.getImageData(x0, 0, Math.max(1, Math.round(to * copy.width) - x0), copy.height);
      let drawn = 0;
      for (let i = 0; i < data.length; i += 4) if (data[i]! + data[i + 1]! + data[i + 2]! > 30) drawn++;
      return drawn / (data.length / 4);
    }, [from, to]);
}

test("prepare, choose a root, run, and inspect the lineages", async ({ page }) => {
  const errors = watchErrors(page);
  await page.goto("/");

  // 1. Preprocessing runs in the background; the Prepare step lists its stages, then the root picker opens.
  await expect(page.getByRole("heading", { level: 1 })).toContainText("lumenCell");
  await expect(page.getByText("Where does the trajectory start?")).toBeVisible({ timeout: 120_000 });
  await expect(page.locator("canvas").first()).toBeVisible();

  // 2. Nothing runs until a root is chosen and pseudotime waits for a run, but 3-D opens after preparation.
  await expect(page.getByRole("button", { name: /Run PAGA/ })).toHaveCount(0);
  await expect(page.getByRole("combobox", { name: "Colour by" })).toHaveCount(0);
  await expect(page.getByText("Pseudotime appears after a run")).toBeVisible();
  await page.getByRole("button", { name: "3-D", exact: true }).click();
  await expect(page.getByText(/DC1–DC3 · 3 of \d+ diffusion components/)).toBeVisible();
  await page.getByRole("button", { name: "2-D", exact: true }).click();

  // 3. Use the root stored in the file (the first suggestion).
  const suggestions = page.getByRole("list", { name: "Suggested roots" });
  await expect(suggestions).toContainText("Stored in the file");
  await suggestions.getByRole("button").first().click();
  await expect(page.getByRole("region", { name: "Pick a root" })).toContainText("Chosen");

  // 4. Run the method; the run finishes in the worker and two lineages appear.
  await page.getByRole("button", { name: /Run PAGA \+ DPT/ }).click();
  const lineages = page.getByRole("list", { name: "Lineages" });
  await expect(lineages.getByRole("listitem")).toHaveCount(2, { timeout: 120_000 });
  await expect(page.getByRole("combobox", { name: "Colour by" })).toHaveValue("pseudotime");

  // 5. The Evidence sheet states the verdict in plain words (never "supported") and where the result came from.
  const evidence = page.getByRole("complementary", { name: "Evidence" });
  await expect(evidence.getByRole("region", { name: "Verdict" })).toContainText(
    /Not shown to be a tree|Clones do not back these lineages|Not supported: cells left without a pseudotime/,
  );
  await expect(evidence.getByRole("region", { name: "Provenance" })).toContainText("X_diffmap");
  await expect(evidence.getByRole("region", { name: "Provenance" })).toContainText("experimental");
  // Each lineage's curve fidelity, on the 2-D embedding and on DC1-DC3.
  await expect(evidence.getByRole("region", { name: "Curve fidelity" })).toContainText(/DC1–DC3 -?\d\.\d\d/);

  // 6. Focus a lineage.
  const first = lineages.getByRole("button").first();
  await first.click();
  await expect(first).toHaveAttribute("aria-pressed", "true");
  await expect(first).toContainText("focused");

  // 7. The 3-D view: DC1-DC3 with the gimbal, then back to 2-D.
  await page.getByRole("button", { name: "3-D", exact: true }).click();
  await expect(page.getByRole("button", { name: "3-D", exact: true })).toHaveAttribute("aria-pressed", "true");
  await expect(page.getByText(/DC1–DC3 · 3 of \d+ diffusion components/)).toBeVisible();
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC3" })).toBeVisible();
  // The cells are drawn (the depth shader compiled and ran), and the fog can be turned off.
  await expect.poll(() => drawnShare(page)).toBeGreaterThan(0.001);
  await expect(page.getByText(/far cells fade/)).toBeVisible();
  await page.getByRole("button", { name: "Display options" }).click();
  await page.getByRole("checkbox", { name: "Fade far cells in depth" }).uncheck();
  await expect(page.getByText(/far cells fade/)).toHaveCount(0);
  // Other diffusion components: z = DC4 relabels the caption and the gimbal; DC3 again restores them.
  await page.getByRole("combobox", { name: "Z axis" }).selectOption("4");
  await expect(page.getByText(/DC1, DC2, DC4 · 3 of \d+ diffusion components/)).toBeVisible();
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC4" })).toBeVisible();
  await page.getByRole("combobox", { name: "Z axis" }).selectOption("3");
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC3" })).toBeVisible();
  // Graphs drawn over the view: PAGA and the neighbour graph load (and draw without a console error).
  const graphs = Promise.all([page.waitForResponse(/\/graph\/paga\?embedding=/), page.waitForResponse(/\/graph\/knn/)]);
  await page.getByRole("checkbox", { name: "PAGA cluster graph" }).check();
  await page.getByRole("checkbox", { name: "Neighbour graph" }).check();
  for (const res of await graphs) expect(res.status()).toBe(200);
  await page.getByRole("checkbox", { name: "PAGA cluster graph" }).uncheck();
  await page.getByRole("checkbox", { name: "Neighbour graph" }).uncheck();
  await page.getByRole("button", { name: "Display options" }).click();
  // Keyboard: Shift+arrow orbits, R resets. Auto-rotate stays off under reduced motion.
  const gimbal = () =>
    page
      .getByRole("img", { name: "Orientation of DC1, DC2 and DC3" })
      .locator("line")
      .evaluateAll((lines) => lines.map((l) => l.getAttribute("x2")).join());
  const tilted = await gimbal();
  await page.getByRole("application", { name: /^Plot/ }).focus();
  await page.keyboard.press("Shift+ArrowLeft");
  await expect.poll(gimbal).not.toBe(tilted);
  await page.keyboard.press("r");
  await expect.poll(gimbal).toBe(tilted);
  await expect(page.getByRole("button", { name: "Auto-rotate" })).toBeDisabled();
  // Saved views: save one, turn away, go back to it. Playing a tour is off under reduced motion.
  await page.getByRole("button", { name: "Camera tour" }).click();
  await page.getByRole("button", { name: "Save this view" }).click();
  const saved = await gimbal();
  await page.getByRole("application", { name: /^Plot/ }).focus();
  await page.keyboard.press("Shift+ArrowLeft");
  await expect.poll(gimbal).not.toBe(saved);
  await page.getByRole("button", { name: "Go to view 1" }).click();
  await expect.poll(gimbal).toBe(saved);
  await expect(page.getByRole("button", { name: "Play tour" })).toBeDisabled();
  await page.getByRole("button", { name: "Camera tour" }).click();
  // Lasso: a drag around the middle of the plot selects cells; Esc clears them.
  await page.getByRole("button", { name: "Lasso select" }).click();
  const selection = page.getByRole("group", { name: "Selection" });
  await expect(selection).toContainText("No cells selected");
  const [cx, cy] = [680, 440];
  await page.mouse.move(cx + 260, cy);
  await page.mouse.down();
  for (let a = 0; a <= 2 * Math.PI + 0.01; a += Math.PI / 12) await page.mouse.move(cx + 260 * Math.cos(a), cy + 220 * Math.sin(a));
  await page.mouse.up();
  await expect(selection).toContainText(/[\d,]+ cells selected/);
  await page.keyboard.press("Escape");
  await expect(selection).toContainText("No cells selected");
  await page.getByRole("button", { name: "Lasso select" }).click();

  // Split: the 2-D embedding and DC1-DC3 side by side, each drawn in its own pane.
  await page.getByRole("button", { name: "Split", exact: true }).click();
  await expect(page.getByText(/a 2-D display projection/)).toBeVisible();
  await expect(page.getByText(/DC1–DC3 · 3 of \d+ diffusion components/)).toBeVisible();
  const stage = (await page.locator("canvas").first().boundingBox())!;
  const half = (stage.width - 304 - 384) / 2; // between the rail and the open Evidence sheet
  const [leftX, rightX, midY] = [stage.x + 304 + half / 2, stage.x + 304 + 1.5 * half, stage.y + (stage.height - 80) / 2];
  await expect.poll(() => drawnShare(page, 304 / stage.width, (304 + half) / stage.width)).toBeGreaterThan(0.001);
  await expect.poll(() => drawnShare(page, (304 + half) / stage.width, (304 + 2 * half) / stage.width)).toBeGreaterThan(0.001);
  // Keys move the pane the pointer was over last: Shift+arrow over the 2-D pane leaves the 3-D pane alone.
  const sideBySide = await gimbal();
  await page.getByRole("application", { name: /^Plot/ }).focus();
  await page.mouse.move(leftX, midY);
  await page.keyboard.press("Shift+ArrowLeft");
  await page.waitForTimeout(400);
  expect(await gimbal()).toBe(sideBySide);
  await page.mouse.move(rightX, midY);
  await page.keyboard.press("Shift+ArrowLeft");
  await expect.poll(gimbal).not.toBe(sideBySide);
  // A lasso in the 2-D pane selects cells (both panes show them).
  await page.getByRole("button", { name: "Lasso select" }).click();
  await page.mouse.move(leftX + 110, midY);
  await page.mouse.down();
  for (let a = 0; a <= 2 * Math.PI + 0.01; a += Math.PI / 12) await page.mouse.move(leftX + 110 * Math.cos(a), midY + 150 * Math.sin(a));
  await page.mouse.up();
  await expect(selection).toContainText(/[\d,]+ cells selected/);
  await page.keyboard.press("Escape");
  await page.getByRole("button", { name: "Lasso select" }).click();
  await page.getByRole("button", { name: "2-D", exact: true }).click();
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC3" })).toHaveCount(0);

  // 8. Scrub pseudotime with the keyboard; the gene trends open on the same cursor.
  const scrub = page.getByRole("slider", { name: "Pseudotime cursor" });
  await scrub.focus();
  for (let i = 0; i < 40; i++) await scrub.press("ArrowRight");
  await expect(page.getByText("pt 0.40")).toBeVisible();
  const trends = page.getByRole("region", { name: "Gene trends" });
  await expect(trends.getByRole("list", { name: "Genes" }).getByRole("listitem")).toHaveCount(4);

  // 9. Play, then pause.
  await page.getByRole("button", { name: "Play", exact: true }).click();
  await page.getByRole("button", { name: "Pause", exact: true }).click();
  await expect(page.getByRole("button", { name: "Play", exact: true })).toBeVisible();

  // 10. Hide and show the Evidence sheet; switch to the light theme.
  await evidence.getByRole("button", { name: "Hide evidence" }).click();
  await expect(evidence).toHaveCount(0);
  await page.getByRole("button", { name: "Show evidence" }).click();
  await expect(page.getByRole("complementary", { name: "Evidence" })).toBeVisible();
  await page.getByRole("button", { name: "Switch to light theme" }).click();
  await expect(page.locator("html")).toHaveAttribute("data-theme", "light");

  // 11. Export PNG: a real PNG of the plot, with its legend and provenance drawn in.
  const download = page.waitForEvent("download");
  await page.getByRole("button", { name: "Export PNG" }).click();
  const png = fs.readFileSync((await (await download).path())!);
  expect(png.subarray(0, 8).toString("hex")).toBe("89504e470d0a1a0a");

  // 12. The view is in the address: after a reload the run, focus, cursor, 3-D view and axes are back.
  await page.getByRole("button", { name: "3-D", exact: true }).click();
  await page.getByRole("button", { name: "Display options" }).click();
  await page.getByRole("combobox", { name: "Z axis" }).selectOption("4");
  await page.getByRole("button", { name: "Display options" }).click();
  const cursorText = await page.getByText(/^pt \d\.\d\d$/).innerText();
  await expect.poll(() => page.evaluate(() => location.hash)).toContain("axes=1%2C2%2C4");
  await page.reload();
  await expect(page.getByRole("list", { name: "Lineages" }).getByRole("listitem")).toHaveCount(2);
  await expect(page.getByRole("button", { name: "3-D", exact: true })).toHaveAttribute("aria-pressed", "true");
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC4" })).toBeVisible();
  await expect(page.getByRole("list", { name: "Lineages" }).getByRole("button").first()).toHaveAttribute("aria-pressed", "true");
  await expect(page.getByText(cursorText, { exact: true })).toBeVisible();

  expect(errors, errors.join("\n")).toEqual([]);
});
