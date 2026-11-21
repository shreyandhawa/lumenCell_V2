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

test("prepare, choose a root, run, and inspect the lineages", async ({ page }) => {
  const errors = watchErrors(page);
  await page.goto("/");

  // 1. Preprocessing runs in the background; the Prepare step lists its stages, then the root picker opens.
  await expect(page.getByRole("heading", { level: 1 })).toContainText("lumenCell");
  await expect(page.getByText("Where does the trajectory start?")).toBeVisible({ timeout: 120_000 });
  await expect(page.locator("canvas").first()).toBeVisible();

  // 2. Nothing runs until a root is chosen; pseudotime and 3-D wait for a run.
  await expect(page.getByRole("button", { name: /Run PAGA/ })).toHaveCount(0);
  await expect(page.getByRole("combobox", { name: "Colour by" })).toHaveCount(0);
  await expect(page.getByText("Pseudotime appears after a run")).toBeVisible();
  await expect(page.getByRole("button", { name: "3-D" })).toBeDisabled();

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

  // 6. Focus a lineage.
  const first = lineages.getByRole("button").first();
  await first.click();
  await expect(first).toHaveAttribute("aria-pressed", "true");
  await expect(first).toContainText("focused");

  // 7. The 3-D view: DC1-DC3 with the gimbal, then back to 2-D.
  await page.getByRole("button", { name: "3-D" }).click();
  await expect(page.getByRole("button", { name: "3-D" })).toHaveAttribute("aria-pressed", "true");
  await expect(page.getByText(/DC1–DC3 · 3 of \d+ diffusion components/)).toBeVisible();
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC3" })).toBeVisible();
  await page.getByRole("button", { name: "2-D" }).click();
  await expect(page.getByRole("img", { name: "Orientation of DC1, DC2 and DC3" })).toHaveCount(0);

  // 8. Scrub pseudotime with the keyboard; the gene trends open on the same cursor.
  const scrub = page.getByRole("slider", { name: "Pseudotime cursor" });
  await scrub.focus();
  for (let i = 0; i < 40; i++) await scrub.press("ArrowRight");
  await expect(page.getByText("pt 0.40")).toBeVisible();
  const trends = page.getByRole("region", { name: "Gene trends" });
  await expect(trends.getByRole("list", { name: "Genes" }).getByRole("listitem")).toHaveCount(4);

  // 9. Play, then pause.
  await page.getByRole("button", { name: "Play" }).click();
  await page.getByRole("button", { name: "Pause" }).click();
  await expect(page.getByRole("button", { name: "Play" })).toBeVisible();

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
