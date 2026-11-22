/// <reference lib="dom" />
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

/** Share of the plot's pixels that are not the dark stage: 0 means nothing was drawn. */
function drawnShare(page: Page): Promise<number> {
  return page
    .locator("canvas")
    .first()
    .evaluate((gl: HTMLCanvasElement) => {
      const copy = document.createElement("canvas");
      copy.width = gl.width;
      copy.height = gl.height;
      const ctx = copy.getContext("2d")!;
      ctx.drawImage(gl, 0, 0);
      const { data } = ctx.getImageData(0, 0, copy.width, copy.height);
      let drawn = 0;
      for (let i = 0; i < data.length; i += 4) if (data[i]! + data[i + 1]! + data[i + 2]! > 30) drawn++;
      return drawn / (data.length / 4);
    });
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
  await expect(page.getByText(/DC1–DC3 · the first 3 diffusion components/)).toBeVisible();
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

  // 11. The run is kept: after a reload it is open again.
  await page.reload();
  await expect(page.getByRole("list", { name: "Lineages" }).getByRole("listitem")).toHaveCount(2);

  expect(errors, errors.join("\n")).toEqual([]);
});
