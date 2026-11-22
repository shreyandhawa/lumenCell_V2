import { render, screen, within } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import type * as T from "@/api/generated";
import { ReferenceResults } from "@/components/trajectory/ReferenceResults";

const source = "benchmarks/results/real-2026-11-20.json @ a7d69c4";

const row = (over: Partial<T.ReferenceRow>): T.ReferenceRow => ({
  dataset: "Pancreas (mouse, E15.5)",
  metric: "Branch accuracy",
  gate: "branch_accuracy",
  score: 0.6306,
  threshold: "≥ 0.8",
  passed: false,
  detail: null,
  source,
  ...over,
});

describe("ReferenceResults", () => {
  it("says these are reference results and shows each dataset with its source once", () => {
    render(
      <ReferenceResults
        rows={[
          row({}),
          row({ metric: "Lineage recall", gate: "lineage_recall", score: 0.5, threshold: "≥ 1.0" }),
          row({ dataset: "Bone marrow (human CD34+)", score: 0.9426, threshold: "≥ 0.75", passed: true }),
        ]}
      />,
    );
    expect(screen.getByText(/not a score for your dataset/)).toBeInTheDocument();
    const pancreas = screen.getByRole("list", { name: "Pancreas (mouse, E15.5): reference results" });
    expect(within(pancreas).getAllByRole("listitem")).toHaveLength(2);
    expect(screen.getAllByText(source)).toHaveLength(2);
    expect(screen.getByRole("heading", { name: "Bone marrow (human CD34+)" })).toBeInTheDocument();
  });

  it("marks a failed gate as failed, with its threshold and the failing cases", () => {
    render(
      <ReferenceResults
        rows={[row({ detail: "worst case; fails 3 of 12 cases", score: 0.8125, threshold: "≥ 0.9", gate: "branch_accuracy" })]}
      />,
    );
    const item = screen.getByRole("listitem");
    expect(within(item).getByText("Fail")).toBeInTheDocument();
    expect(item).toHaveTextContent("needs ≥ 0.9");
    expect(item).toHaveTextContent("worst case; fails 3 of 12 cases");
    expect(within(item).getByText("0.812")).toBeInTheDocument();
  });

  it("never shows an unmeasured gate as a pass, and says why", () => {
    render(
      <ReferenceResults
        rows={[
          row({
            dataset: "LCMV CD8+ T cells (mouse)",
            score: null,
            passed: null,
            detail: "lineage gates: not measured (no gated lineage, declared 2026-11-20)",
          }),
        ]}
      />,
    );
    const item = screen.getByRole("listitem");
    expect(item).toHaveTextContent("Not measured");
    expect(item).not.toHaveTextContent("Pass");
    expect(item).toHaveTextContent("lineage gates: not measured");
  });

  it("rounds a score down, so a value under its bound never reads as the bound", () => {
    render(<ReferenceResults rows={[row({ score: 0.7996, threshold: "≥ 0.8", passed: false })]} />);
    expect(screen.getByText("0.799")).toBeInTheDocument();
  });

  it("says when a method has no published results", () => {
    render(<ReferenceResults rows={[]} />);
    expect(screen.getByText("This method has no published benchmark results.")).toBeInTheDocument();
  });
});
