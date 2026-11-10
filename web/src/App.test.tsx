import { render, screen } from "@testing-library/react";
import { describe, expect, it } from "vitest";
import App from "@/App";

describe("App", () => {
  it("renders the app name and the empty state", () => {
    render(<App />);
    expect(screen.getByRole("heading", { level: 1 })).toHaveTextContent("lumenCell");
    expect(screen.getByText("No dataset loaded.")).toBeInTheDocument();
  });
});
