import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";
import { checkParams, fieldError, initialValues, ParamForm, type ParamSchema } from "@/components/trajectory/ParamForm";

// The shapes the engine's methods produce (see engine/methods/*.py), trimmed to one of each kind.
const schema: ParamSchema = {
  properties: {
    n_dcs: { title: "N Dcs", description: "Diffusion components used by DPT.", type: "integer", default: 10, minimum: 2, maximum: 14 },
    min_connectivity: { title: "Min Connectivity", type: "number", default: 0.05, minimum: 0, maximum: 1 },
    ppt_lambda: { title: "Ppt Lambda", type: "number", default: 1.0, exclusiveMinimum: 0 },
    n_eigs: {
      title: "N Eigs",
      description: "Empty = Palantir's eigengap rule.",
      anyOf: [{ type: "integer", minimum: 2, maximum: 30 }, { type: "null" }],
      default: null,
    },
    terminal_detection: {
      title: "Terminal Detection",
      type: "string",
      enum: ["paga_leaves", "palantir"],
      default: "paga_leaves",
    },
    cluster_key: { title: "Cluster Key", type: "string", default: "lumencell_leiden" },
  },
};

describe("initialValues", () => {
  it("starts each field at its schema default, empty when the default is null", () => {
    expect(initialValues(schema)).toEqual({
      n_dcs: "10",
      min_connectivity: "0.05",
      ppt_lambda: "1",
      n_eigs: "",
      terminal_detection: "paga_leaves",
      cluster_key: "lumencell_leiden",
    });
  });
});

describe("fieldError", () => {
  const n_dcs = schema.properties!.n_dcs!;

  it("accepts values inside the range and rejects values outside it", () => {
    expect(fieldError(n_dcs, "2")).toBeNull();
    expect(fieldError(n_dcs, "14")).toBeNull();
    expect(fieldError(n_dcs, "1")).toBe("Must be at least 2.");
    expect(fieldError(n_dcs, "15")).toBe("Must be at most 14.");
  });

  it("rejects fractions for integers and text that is not a number", () => {
    expect(fieldError(n_dcs, "8.5")).toBe("Enter a whole number.");
    expect(fieldError(n_dcs, "eight")).toBe("Enter a number.");
    expect(fieldError(n_dcs, "")).toBe("Enter a value.");
  });

  it("treats an exclusive minimum as strict", () => {
    const lambda = schema.properties!.ppt_lambda!;
    expect(fieldError(lambda, "0")).toBe("Must be more than 0.");
    expect(fieldError(lambda, "0.5")).toBeNull();
  });

  it("allows an empty optional value, and only listed choices", () => {
    expect(fieldError(schema.properties!.n_eigs!, "")).toBeNull();
    expect(fieldError(schema.properties!.terminal_detection!, "palantir")).toBeNull();
    expect(fieldError(schema.properties!.terminal_detection!, "other")).toBe("Choose one of the listed values.");
  });
});

describe("checkParams", () => {
  it("sends typed values, null for an empty optional field, and reports the fields that cannot be sent", () => {
    const ok = checkParams(schema, { ...initialValues(schema), n_dcs: "8" });
    expect(ok.errors).toEqual({});
    expect(ok.params).toEqual({
      n_dcs: 8,
      min_connectivity: 0.05,
      ppt_lambda: 1,
      n_eigs: null,
      terminal_detection: "paga_leaves",
      cluster_key: "lumencell_leiden",
    });

    const bad = checkParams(schema, { ...initialValues(schema), n_dcs: "99" });
    expect(bad.errors).toEqual({ n_dcs: "Must be at most 14." });
  });
});

describe("ParamForm", () => {
  it("shows the defaults and the range, and reports each edit with the field's name", async () => {
    const onChange = vi.fn();
    const user = userEvent.setup();
    render(<ParamForm schema={schema} values={initialValues(schema)} errors={{}} onChange={onChange} />);

    const dcs = screen.getByRole("spinbutton", { name: "N Dcs" });
    expect(dcs).toHaveValue(10);
    expect(dcs).toHaveAttribute("min", "2");
    expect(dcs).toHaveAttribute("max", "14");
    expect(screen.getByText("Diffusion components used by DPT. · 2 to 14")).toBeInTheDocument();

    await user.type(dcs, "1");
    expect(onChange).toHaveBeenLastCalledWith("n_dcs", "101");

    const method = screen.getByRole("combobox", { name: "Terminal Detection" });
    await user.selectOptions(method, "palantir");
    expect(onChange).toHaveBeenLastCalledWith("terminal_detection", "palantir");
  });

  it("marks a field with its error", () => {
    render(
      <ParamForm
        schema={schema}
        values={{ ...initialValues(schema), n_dcs: "99" }}
        errors={{ n_dcs: "Must be at most 14." }}
        onChange={() => {}}
      />,
    );
    expect(screen.getByRole("spinbutton", { name: "N Dcs" })).toHaveAttribute("aria-invalid", "true");
    expect(screen.getByText("Must be at most 14.")).toBeInTheDocument();
  });
});
