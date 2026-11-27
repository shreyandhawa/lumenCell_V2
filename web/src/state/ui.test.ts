import { afterEach, describe, expect, it, vi } from "vitest";
import { initialTheme, useUi } from "@/state/ui";

const systemPrefers = (light: boolean) =>
  vi.stubGlobal("matchMedia", (q: string) => ({ matches: light && q.includes("light") }) as MediaQueryList);

afterEach(() => {
  localStorage.clear();
  vi.unstubAllGlobals();
});

describe("theme", () => {
  it("starts from the saved choice, else the system's, else dark", () => {
    systemPrefers(true);
    expect(initialTheme()).toBe("light");
    localStorage.setItem("lumencell-theme", "dark");
    expect(initialTheme()).toBe("dark"); // the saved choice beats the system
    localStorage.clear();
    systemPrefers(false);
    expect(initialTheme()).toBe("dark");
  });

  it("sets data-theme and remembers the choice", () => {
    useUi.getState().setTheme("light");
    expect(document.documentElement.dataset.theme).toBe("light");
    expect(localStorage.getItem("lumencell-theme")).toBe("light");
    expect(initialTheme()).toBe("light");
  });
});
