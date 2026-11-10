export default function App() {
  return (
    <div className="flex h-full flex-col bg-surface text-ink">
      <header className="flex h-10 shrink-0 items-center border-b border-line px-3">
        <h1 className="text-sm font-semibold tracking-tight">
          lumen<span className="text-accent">Cell</span>
        </h1>
      </header>
      <main className="flex flex-1 items-center justify-center text-sm text-ink-muted">
        No dataset loaded.
      </main>
    </div>
  );
}
