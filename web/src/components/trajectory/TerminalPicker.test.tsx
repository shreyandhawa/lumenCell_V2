import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { describe, expect, it, vi } from "vitest";
import { TerminalPicker } from "@/components/trajectory/TerminalPicker";

describe("TerminalPicker", () => {
  it("offers the dataset's labels in their own order and reports the ticked ones in that order", async () => {
    const onChange = vi.fn();
    const user = userEvent.setup();
    render(<TerminalPicker column="lumencell_leiden" labels={["0", "1", "2"]} value={["2"]} onChange={onChange} />);

    expect(screen.getByRole("group", { name: "Terminal states, labels of lumencell_leiden" })).toBeInTheDocument();
    expect(screen.getByRole("checkbox", { name: "2" })).toBeChecked();
    await user.click(screen.getByRole("checkbox", { name: "0" }));
    expect(onChange).toHaveBeenLastCalledWith(["0", "2"]);
    await user.click(screen.getByRole("checkbox", { name: "2" }));
    expect(onChange).toHaveBeenLastCalledWith([]);
  });

  it("says the method finds the terminals when the column has no labels", () => {
    render(<TerminalPicker column="batch" labels={null} value={[]} onChange={() => {}} />);
    expect(screen.getByText(/batch has no labels in this dataset/)).toBeInTheDocument();
    expect(screen.queryByRole("checkbox")).not.toBeInTheDocument();
  });
});
