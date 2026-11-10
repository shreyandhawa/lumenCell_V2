import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { render, screen, waitFor, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { http, HttpResponse } from "msw";
import { setupServer } from "msw/node";
import { afterAll, afterEach, beforeAll, beforeEach, describe, expect, it, vi } from "vitest";
import type * as T from "@/api/generated";
import App from "@/App";
import { useUi } from "@/state/ui";
import * as fx from "@/test/fixtures";

// deck.gl needs WebGL, which jsdom lacks; the plot is covered by the Playwright test.
vi.mock("@/components/Scatter", () => ({ Scatter: () => <div data-testid="scatter" /> }));

const api = (path: string) => `*/api/v1${path}`;

interface State {
  dataset: T.DatasetInfo;
  root: T.ResolvedRoot | null;
  runs: T.RunSummary[];
  putRoot: unknown[];
}
let state: State;

const server = setupServer(
  http.get(api("/dataset"), () => HttpResponse.json(state.dataset)),
  http.get(api("/root"), () => HttpResponse.json({ root: state.root })),
  http.put(api("/root"), async ({ request }) => {
    state.putRoot.push(await request.json());
    state.root = fx.resolvedRoot;
    return HttpResponse.json(fx.resolvedRoot);
  }),
  http.get(api("/root/suggestions"), () => HttpResponse.json(fx.suggestions)),
  http.get(api("/methods"), () => HttpResponse.json(fx.methods)),
  http.get(api("/trajectories"), () => HttpResponse.json({ runs: state.runs })),
  http.get(api("/trajectories/run1"), () => HttpResponse.json(fx.run)),
  http.get(api("/trajectories/run1/curves"), () =>
    HttpResponse.json({ run_id: "run1", embedding: "umap", curve_source: "binned", curves: [], milestones: [], warnings: [] }),
  ),
);

beforeAll(() => server.listen({ onUnhandledFrame: "error" }));
afterEach(() => server.resetHandlers());
afterAll(() => server.close());
beforeEach(() => {
  state = { dataset: fx.dataset("ready"), root: null, runs: [], putRoot: [] };
  useUi.setState({ embedding: null, colorBy: null, runId: null, ptRange: [0, 1], focusLineage: null, selection: null });
});

function renderApp() {
  const client = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={client}>
      <App />
    </QueryClientProvider>,
  );
}

describe("App", () => {
  it("shows preprocessing progress until the dataset is prepared", async () => {
    state.dataset = fx.dataset("running");
    renderApp();
    const stages = await screen.findByRole("list", { name: "Preprocessing stages" });
    expect(within(stages).getByText("PCA")).toBeInTheDocument();
    expect(screen.getByRole("status")).toHaveTextContent("Preparing: PCA");
    expect(screen.queryByText("Where does the trajectory start?")).not.toBeInTheDocument();
  });

  it("asks for a root and never picks one itself", async () => {
    renderApp();
    expect(await screen.findByText("Where does the trajectory start?")).toBeInTheDocument();
    expect(await screen.findByText(/earliest day \(day 0\)/)).toBeInTheDocument();
    expect(state.putRoot).toEqual([]);
    expect(screen.queryByRole("button", { name: /Run PAGA/ })).not.toBeInTheDocument();
  });

  it("keeps pseudotime and lineage colouring disabled until there is a finished run", async () => {
    const user = userEvent.setup();
    renderApp();
    await user.click(await screen.findByRole("combobox", { name: "Colour by" }));
    const option = await screen.findByRole("option", { name: /Pseudotime/ });
    expect(option).toHaveAttribute("aria-disabled", "true");
    expect(option).toHaveTextContent("after a run");
  });

  it("saves the root the user picks, then offers methods", async () => {
    const user = userEvent.setup();
    renderApp();
    await user.click(await screen.findByRole("button", { name: /Start here/ }));
    await waitFor(() => expect(state.putRoot).toEqual([{ root: { kind: "cluster", key: "cell_type", value: "stem" } }]));
    expect(await screen.findByRole("button", { name: /Run PAGA \+ DPT/ })).toBeInTheDocument();
    expect(screen.getByText(/Experimental/)).toBeInTheDocument();
  });

  it("shows lineages and the trust panel for a finished run", async () => {
    state.root = fx.resolvedRoot;
    state.runs = [{ run_id: "run1", method: "paga_dpt", state: "done", created_at: "2026-11-10T00:00:00Z", n_lineages: 1, error: null }];
    renderApp();
    expect(await screen.findByText("Erythroid")).toBeInTheDocument();
    const trust = screen.getByRole("region", { name: "How far to trust this" });
    expect(within(trust).getByText("X_diffmap · 10 dims")).toBeInTheDocument();
    expect(within(trust).getByText("2 cells are not connected.")).toBeInTheDocument();
    await userEvent.setup().click(screen.getByRole("combobox", { name: "Colour by" }));
    expect(await screen.findByRole("option", { name: "Pseudotime" })).not.toHaveAttribute("aria-disabled", "true");
  });

  it("explains how to start the server when the API is unreachable", async () => {
    server.use(http.get(api("/dataset"), () => HttpResponse.error()));
    renderApp();
    expect(await screen.findByRole("alert")).toHaveTextContent("lumencell launch");
  });
});
