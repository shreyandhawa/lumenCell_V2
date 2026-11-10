import { expect, type Page, test } from "@playwright/test";

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

  // 1. Preprocessing runs in the background; the panel shows its stages, then the root picker.
  await expect(page.getByRole("heading", { level: 1 })).toContainText("lumenCell");
  await expect(page.getByText("Where does the trajectory start?")).toBeVisible({ timeout: 120_000 });
  await expect(page.locator("canvas").first()).toBeVisible();

  // 2. Nothing runs until a root is chosen; pseudotime colouring is unavailable.
  const colourBy = page.getByRole("combobox", { name: "Colour by" });
  await expect(colourBy.locator('option[value="pseudotime"]')).toBeDisabled();
  await expect(page.getByRole("button", { name: /Run PAGA/ })).toHaveCount(0);

  // 3. Use the root stored in the file (the first suggestion).
  const suggestions = page.getByLabel("Suggested roots");
  await expect(suggestions).toContainText("stored in file");
  await suggestions.getByRole("button", { name: "Start here" }).first().click();
  await expect(page.getByRole("region", { name: "Root" })).toContainText("the root stored in the file");

  // 4. Run the method; the run finishes in the worker and two lineages appear.
  await page.getByRole("button", { name: /Run PAGA \+ DPT/ }).click();
  const lineages = page.getByRole("list", { name: "Lineages" });
  await expect(lineages.getByRole("listitem")).toHaveCount(2, { timeout: 120_000 });
  await expect(colourBy).toHaveValue("pseudotime");
  await expect(page.getByLabel("Legend").or(page.getByText("pseudotime", { exact: true })).first()).toBeVisible();

  // 5. Rename a lineage; the name sticks after a reload.
  const first = lineages.getByRole("listitem").first();
  await first.hover();
  await first.getByRole("button", { name: /^Rename/ }).click();
  const input = page.getByRole("textbox", { name: /^Rename/ });
  await input.fill("Branch one");
  await input.press("Enter");
  await expect(lineages).toContainText("Branch one");

  // 6. Narrow the pseudotime window with the keyboard.
  const from = page.getByRole("slider", { name: "Pseudotime from" });
  await from.focus();
  for (let i = 0; i < 40; i++) await from.press("ArrowRight");
  const window_ = page.getByRole("region", { name: "Pseudotime window" });
  await expect(window_).toContainText("0.40");
  await window_.getByRole("button", { name: "Reset" }).click();
  await expect(window_).toContainText("0.00");

  // 7. The trust panel states where the result came from.
  const trust = page.getByRole("region", { name: "How far to trust this" });
  await expect(trust).toContainText("X_diffmap");
  await expect(trust).toContainText("experimental");

  // 8. Genes along pseudotime: rank the top temporal genes, then plot one.
  const genes = page.getByRole("region", { name: "Genes along pseudotime" });
  await genes.getByRole("tab", { name: "Top temporal genes" }).click();
  const top = genes.getByRole("list", { name: "Top temporal genes" });
  await expect(top.getByRole("listitem")).toHaveCount(30);
  const firstGene = (await top.getByRole("button").first().textContent())!.trim();
  await top.getByRole("button").first().click();
  await expect(colourBy).toHaveValue(`gene:${firstGene}`);
  await genes.getByRole("tab", { name: "Trends" }).click();
  await expect(genes.getByRole("figure", { name: `${firstGene} along pseudotime` })).toBeVisible();

  await page.reload();
  await expect(page.getByRole("list", { name: "Lineages" })).toContainText("Branch one");

  expect(errors, errors.join("\n")).toEqual([]);
});
