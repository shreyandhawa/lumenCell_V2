"""Command-line entry point (`lumencell`)."""

import signal
import threading
from pathlib import Path
from types import FrameType
from typing import TYPE_CHECKING, Annotated

import typer

from lumencell import __version__

if TYPE_CHECKING:
    from lumencell.data.loader import Dataset
    from lumencell.prep.pipeline import PrepResult
    from lumencell.prep.stages import PrepParams

app = typer.Typer(
    name="lumencell",
    help="Explore pseudotime and lineage trajectories in single-cell data.",
    no_args_is_help=True,
    add_completion=False,
)


@app.callback()
def main() -> None:
    """Explore pseudotime and lineage trajectories in single-cell data."""


@app.command()
def version() -> None:
    """Print the lumenCell version."""
    typer.echo(__version__)


DataArg = Annotated[Path, typer.Argument(exists=True, dir_okay=False, help="An .h5ad file.")]
XIsOpt = Annotated[
    str | None, typer.Option("--x-is", help="Override what X holds: counts or lognorm (default: detect).")
]
UmapOpt = Annotated[bool, typer.Option(help="Compute a UMAP even if the file has an embedding.")]
RootOpt = Annotated[
    str | None,
    typer.Option(help="Where the trajectory starts: cell:<name>, cluster:<column>=<value> or file."),
]


def _prep_params(recompute_umap: bool) -> "PrepParams":
    from lumencell.prep.stages import PrepParams

    return PrepParams(recompute_umap=recompute_umap)


def _prepared(data: Path, x_is: str | None, recompute_umap: bool) -> "tuple[Dataset, PrepResult, PrepParams]":
    from lumencell.data.loader import DatasetError, load
    from lumencell.prep.cache import PrepCache

    try:
        ds = load(data, x_is)  # type: ignore[arg-type]
    except DatasetError as exc:
        raise typer.BadParameter(str(exc), param_hint="DATA") from exc
    params = _prep_params(recompute_umap)

    def report(stage: str, i: int, n: int) -> None:
        if stage != "done":
            typer.echo(f"  [{i + 1}/{n}] {stage}", err=True)

    typer.echo(
        f"preparing {data.name} ({ds.n_obs:,} cells x {ds.n_vars:,} genes, X holds {ds.x_kind})", err=True
    )
    prep = PrepCache(ds).run(params, progress=report)
    return ds, prep, params


@app.command()
def prepare(data: DataArg, x_is: XIsOpt = None, recompute_umap: UmapOpt = False) -> None:
    """Run preprocessing ahead of time and cache it (launch reuses the cache)."""
    _ds, prep, _params = _prepared(data, x_is, recompute_umap)
    for r in prep.stages:
        typer.echo(f"{r.name:10s} {r.state:8s} {r.seconds:6.2f}s")


@app.command()
def roots(data: DataArg, x_is: XIsOpt = None) -> None:
    """Print root suggestions with their evidence. Nothing is applied."""
    from lumencell.engine.context import InferenceContext
    from lumencell.engine.roots import format_root, suggest_roots

    ds, prep, params = _prepared(data, x_is, False)
    suggestions = suggest_roots(InferenceContext.from_prep(ds, prep, params))
    if not suggestions:
        typer.echo("no suggestions; pick the starting cluster you know, e.g. --root cluster:<column>=<value>")
    for s in suggestions:
        typer.echo(f"{format_root(s.spec):40s} [{s.source}] {s.evidence}")


@app.command()
def run(
    data: DataArg,
    root: Annotated[str, typer.Option(help="cell:<name>, cluster:<column>=<value> or file.")],
    method: str = "paga_dpt",
    terminals: Annotated[str | None, typer.Option(help="Comma-separated terminal cluster labels.")] = None,
    x_is: XIsOpt = None,
) -> None:
    """Run a trajectory method from the terminal and save the run (launch shows it)."""
    from lumencell.engine.context import InferenceContext
    from lumencell.engine.registry import UnknownMethodError, run_method
    from lumencell.engine.roots import RootError, format_root, parse_root
    from lumencell.engine.store import RunStore, run_id
    from lumencell.prep.cache import dataset_dir, prep_key

    ds, prep, params = _prepared(data, x_is, False)
    ctx = InferenceContext.from_prep(ds, prep, params)
    terms = [t.strip() for t in terminals.split(",")] if terminals else None
    try:
        spec = parse_root(root)
        result = run_method(ctx, method, spec, terminals=terms)
    except (RootError, UnknownMethodError, ValueError) as exc:
        raise typer.BadParameter(str(exc)) from exc
    rid = run_id(method, result.params, format_root(spec), terms, prep_key(ds, params))
    RunStore(dataset_dir(ds.fingerprint)).save(rid, result, terminals=terms)
    typer.echo(f"run {rid}: {len(result.lineages)} lineages, {result.n_unassigned} cells unassigned")
    for lin in result.lineages:
        typer.echo(f"  {lin.name}: {lin.n_cells:,} cells")
    for w in result.warnings:
        typer.echo(f"  {w.severity}: {w.message}")


@app.command()
def export(
    data: DataArg,
    run_id_: Annotated[str, typer.Argument(metavar="RUN_ID", help="From `lumencell run` or the web app.")],
    out: Annotated[Path, typer.Argument(help="Output .h5ad (a new file; the input is never modified).")],
    x_is: XIsOpt = None,
) -> None:
    """Write a run's pseudotime and lineages into a copy of the dataset."""
    from lumencell.engine.store import RunStore
    from lumencell.export import export_run
    from lumencell.prep.cache import dataset_dir

    ds, prep, _params = _prepared(data, x_is, False)
    store = RunStore(dataset_dir(ds.fingerprint))
    if not store.exists(run_id_):
        raise typer.BadParameter(f"no run {run_id_!r} for this dataset", param_hint="RUN_ID")
    try:
        export_run(
            data, out, store.load(run_id_), run_id_, prep.artifacts["leiden"], prep.artifacts.get("umap")
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="OUT") from exc
    typer.echo(f"wrote {out}")


def _exit_on_signal(signum: int, _frame: FrameType | None) -> None:
    raise SystemExit(128 + signum)


@app.command()
def launch(
    data: DataArg,
    root: RootOpt = None,
    method: str = "paga_dpt",
    port: int = 8765,
    host: Annotated[str, typer.Option(help="Interface to bind. Keep the default: the API has no auth.")] = (
        "127.0.0.1"
    ),
    browser: Annotated[bool, typer.Option(help="Open the browser when the server is ready.")] = True,
    strict_port: Annotated[
        bool, typer.Option(help="Fail if the port is taken instead of trying the next.")
    ] = False,
    x_is: XIsOpt = None,
    recompute_umap: UmapOpt = False,
) -> None:
    """Start the local server and open the explorer."""
    import uvicorn

    from lumencell.api.convert import to_api_spec
    from lumencell.api.errors import ApiError
    from lumencell.api.models import TrajectoryRequest
    from lumencell.api.session import Session
    from lumencell.data.loader import DatasetError
    from lumencell.engine.roots import RootError, parse_root
    from lumencell.jobs.manager import JobManager
    from lumencell.server import LOCAL_HOSTS, build_app, free_port, open_when_ready

    if root is not None:
        try:
            parse_root(root)
        except RootError as exc:
            raise typer.BadParameter(str(exc), param_hint="--root") from exc
    try:
        jobs = JobManager()
    except ValueError as exc:  # a malformed LUMENCELL_JOB_TIME_LIMITS
        typer.echo(f"lumencell: {exc}", err=True)
        raise typer.Exit(2) from None
    jobs.start()
    try:
        session = Session(data, jobs, x_is=x_is, prep_params=_prep_params(recompute_umap))
    except DatasetError as exc:
        jobs.shutdown()
        raise typer.BadParameter(str(exc), param_hint="DATA") from exc

    if not session.prep_ready():
        prep_job = session.start_prep()
        typer.echo(f"preparing {data.name} in the background (job {prep_job.job_id})")

    if root is not None:

        def run_when_prepared() -> None:
            job = jobs.find(f"prepare:{session.prep_key}")
            if job is not None:
                job = jobs.wait(job.job_id, timeout=24 * 3600)
                if job.state != "done":
                    reason = f": {job.error}" if job.error else ""
                    typer.echo(f"preparing {data.name} {job.state}{reason}; --root was not applied", err=True)
                    return
            try:
                spec = to_api_spec(parse_root(root))
                session.set_root(spec)
                created, _ = session.runs.create_run(TrajectoryRequest(method=method, root=spec))
            except ApiError as exc:
                # e.g. a cluster name that is not in the file. The server keeps
                # running; the root can still be picked in the app.
                typer.echo(f"could not start {method} from --root {root}: {exc.message}", err=True)
                return
            typer.echo(f"running {method} from {root} (run {created.run_id})")

        threading.Thread(target=run_when_prepared, daemon=True).start()

    try:
        port = free_port(host, port, attempts=1 if strict_port else 20)
    except RuntimeError:
        jobs.shutdown()
        raise typer.BadParameter(f"port {port} is already in use", param_hint="--port") from None
    url = f"http://{host}:{port}"
    trusted = LOCAL_HOSTS if host in LOCAL_HOSTS else [host, *LOCAL_HOSTS]
    typer.echo(f"lumenCell is running at {url}  (Ctrl+C to stop)")
    if browser:
        open_when_ready(url)
    # uvicorn catches SIGTERM, shuts down, then raises it again under the handler that was in
    # place before it started. Python's default handler would end the process there, skipping
    # the `finally` and leaving the job worker running; SystemExit unwinds instead (exit code
    # 143, as a shell reports SIGTERM). SIGINT already unwinds, as KeyboardInterrupt.
    previous = signal.signal(signal.SIGTERM, _exit_on_signal)
    try:
        uvicorn.run(build_app(session, trusted), host=host, port=port, log_level="warning")
    finally:
        signal.signal(signal.SIGTERM, previous or signal.SIG_DFL)  # a second SIGTERM ends it at once
        jobs.shutdown()


cache_app = typer.Typer(help="See what the cache holds, and clean it up.", no_args_is_help=True)
app.add_typer(cache_app, name="cache")


@cache_app.command("list")
def cache_list() -> None:
    """List the datasets in the cache, with their size and when each was last used."""
    import time
    from datetime import datetime

    from lumencell.prep.cache import default_cache_root
    from lumencell.prep.housekeeping import format_ago, format_size, list_cache

    root = default_cache_root()
    entries = list_cache(root)
    if not entries:
        typer.echo(f"the cache at {root} is empty")
        return
    now = time.time()
    typer.echo(f"cache at {root}")
    for e in entries:
        when = datetime.fromtimestamp(e.last_used).strftime("%Y-%m-%d")
        where = "(a symbolic link, never pruned)" if e.is_link else e.source or ""
        typer.echo(
            f"  {e.name or e.fingerprint:24s} {format_size(e.size):>9s}  "
            f"last used {when} ({format_ago(e.last_used, now)})  {e.fingerprint}  {where}".rstrip()
        )
    typer.echo(f"  {len(entries)} datasets, {format_size(sum(e.size for e in entries))} in all")


@cache_app.command("prune")
def cache_prune(
    older_than: Annotated[
        str | None, typer.Option(help="Remove datasets not used for this long: 12h, 30d or 2w.")
    ] = None,
    everything: Annotated[bool, typer.Option("--all", help="Remove every dataset in the cache.")] = False,
    dry_run: Annotated[bool, typer.Option(help="Only say what would be removed.")] = False,
) -> None:
    """Remove datasets from the cache, and leftovers of interrupted jobs.

    Only dataset folders and what is inside them are ever removed. Downloaded benchmark
    data and the saved roots (the state directory) are never touched.
    """
    from lumencell.prep.cache import default_cache_root
    from lumencell.prep.housekeeping import (
        Removal,
        UnsafePath,
        apply_prune,
        check_prune,
        format_size,
        list_cache,
        parse_age,
        plan_prune,
    )

    if (older_than is not None) == everything:
        raise typer.BadParameter("give either --older-than or --all", param_hint="--older-than / --all")
    try:
        age = parse_age(older_than) if older_than is not None else None
    except ValueError as exc:
        raise typer.BadParameter(str(exc), param_hint="--older-than") from exc
    root = default_cache_root()
    entries = list_cache(root)
    if not entries:
        typer.echo(f"the cache at {root} is empty")
        return
    for e in entries:
        if e.is_link:
            typer.echo(
                f"skipped {e.path}: a symbolic link, which may lead outside the cache; it is never removed"
            )
    removals = plan_prune(root, older_than=age, everything=everything)
    try:
        if dry_run:
            check_prune(root, removals)
            outcomes: list[tuple[Removal, str | None]] = [(r, None) for r in removals]
        else:
            outcomes = apply_prune(root, removals)
    except UnsafePath as exc:
        typer.echo(f"lumencell: {exc}; nothing was removed", err=True)
        raise typer.Exit(1) from None
    if not outcomes:
        typer.echo("nothing to remove")
        return
    verb = "would remove" if dry_run else "removed"
    for r, problem in outcomes:
        if problem is None:
            typer.echo(f"{verb} {r.reason} ({format_size(r.size)}): {r.path}")
        else:
            typer.echo(
                f"{problem.split(':')[0]} {r.reason}: {r.path} ({problem})", err=problem != "already gone"
            )
    freed = sum(r.size for r, problem in outcomes if problem is None)
    typer.echo(f"{'would free' if dry_run else 'freed'} {format_size(freed)}")
    if any(problem not in (None, "already gone") for _r, problem in outcomes):
        raise typer.Exit(1)


bench_app = typer.Typer(help="Accuracy benchmarks against known trajectories.", no_args_is_help=True)
app.add_typer(bench_app, name="bench")


@bench_app.command("fetch")
def bench_fetch(
    dataset: Annotated[str, typer.Argument(help="Dataset id, or 'all'.")],
    yes: Annotated[bool, typer.Option("--yes", "-y", help="Download without asking.")] = False,
) -> None:
    """Download a real benchmark dataset (asks first; shows source and size)."""
    from lumencell.bench.datasets import DATASETS, FetchError, fetch, is_cached, local_path, sha256_of

    ids = list(DATASETS) if dataset == "all" else [dataset]
    unknown = [i for i in ids if i not in DATASETS]
    if unknown:
        raise typer.BadParameter(f"unknown dataset {unknown[0]!r}; choose from {', '.join(DATASETS)} or all")
    for i in ids:
        ds = DATASETS[i]
        if is_cached(ds):
            typer.echo(f"{ds.id}: already at {local_path(ds)}")
            continue
        typer.echo(f"{ds.id}: {ds.title} ({ds.citation})")
        typer.echo(f"  from {ds.url}\n  {ds.size / 1e6:,.1f} MB -> {local_path(ds)}")
        if not yes and not typer.confirm("  download it?", default=False):
            typer.echo("  skipped")
            continue
        try:
            path = fetch(ds)
        except FetchError as exc:
            raise typer.BadParameter(str(exc)) from exc
        typer.echo(f"  done, sha256 {sha256_of(path)}")


@bench_app.command("real")
def bench_real(
    methods: Annotated[
        str, typer.Option(help="Comma-separated method names.")
    ] = "paga_dpt,palantir,scfates_ppt",
    datasets: Annotated[str, typer.Option(help="Comma-separated dataset ids; default all downloaded.")] = "",
    out: Annotated[
        Path | None, typer.Option(help="Results JSON (default benchmarks/results/real-<date>.json).")
    ] = None,
) -> None:
    """Run methods on the downloaded real datasets and check the release gates."""
    import json
    import subprocess
    from datetime import UTC, datetime

    from lumencell.bench.datasets import DATASETS, is_cached
    from lumencell.bench.gates import REPO, check_gates, load_thresholds
    from lumencell.bench.runner import run_real_case
    from lumencell.engine.registry import METHODS

    cfg = load_thresholds()["real"]
    ids = [d for d in (datasets.split(",") if datasets else DATASETS) if d]
    missing = [d for d in ids if d not in DATASETS or not is_cached(DATASETS[d])]
    if missing:
        raise typer.BadParameter(f"not downloaded: {', '.join(missing)} (lumencell bench fetch ...)")
    rows: dict[str, dict[str, object]] = {}
    failures = []
    for d in ids:
        for method in methods.split(","):
            if not METHODS[method].available().ok:
                typer.echo(f"{d:24s} {method:12s} skipped (not installed)")
                continue
            row = run_real_case(d, method).row()
            rows.setdefault(method, {})[d] = row
            failed = check_gates(row, cfg[d]["gates"])
            failures += [(method, f) for f in failed]
            verdict = "PASS" if not failed else "FAIL"
            typer.echo(
                f"{d:24s} {method:12s} sp={row['pt_spearman_global']:.3f} "
                f"sp_lin={row['pt_spearman_lineage']:.3f} branch={row['branch_accuracy']:.3f} "
                f"recall={row['lineage_recall']:.2f} terminal={row['terminal_share']:.3f} "
                f"{verdict} ({row['runtime_seconds']}s)"
            )
    for method, f in failures:
        typer.echo(f"gate     {method}: {f}")
    sha = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True
    ).stdout.strip()
    date = datetime.now(UTC).date().isoformat()
    target = out or REPO / "benchmarks" / "results" / f"real-{date}.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(
        json.dumps({"generated": date, "commit": sha, "methods": rows}, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    typer.echo(f"wrote {target}")


@bench_app.command("report")
def bench_report(
    readme: Annotated[bool, typer.Option(help="Rewrite the results section of README.md.")] = False,
) -> None:
    """Markdown tables of the latest synthetic baseline and real-data results."""
    from lumencell.bench.report import render, update_readme

    text = render()
    if readme:
        update_readme(text)
        typer.echo("README.md updated")
    else:
        typer.echo(text)


@bench_app.command("synth")
def bench_synth(
    methods: Annotated[str, typer.Option(help="Comma-separated method names.")] = "paga_dpt",
    topologies: Annotated[
        str | None, typer.Option(help="Comma-separated; default from thresholds.toml.")
    ] = None,
    seeds: Annotated[str | None, typer.Option(help="Comma-separated; default from thresholds.toml.")] = None,
    out: Annotated[Path | None, typer.Option(help="Write the results JSON here.")] = None,
    write_baseline: Annotated[bool, typer.Option(help="Replace benchmarks/baseline.json.")] = False,
) -> None:
    """Run methods on the synthetic matrix; report gates and regressions against the baseline."""
    import json
    import subprocess
    import tempfile
    from datetime import UTC, datetime

    from lumencell.bench.gates import BASELINE, check_gates, check_regression, load_baseline, load_thresholds
    from lumencell.bench.runner import run_synthetic_case

    cfg = load_thresholds()
    syn = cfg["synthetic"]
    topo_list = topologies.split(",") if topologies else list(syn["topologies"])
    seed_list = [int(s) for s in seeds.split(",")] if seeds else list(syn["seeds"])
    baseline = load_baseline()
    rows: dict[str, dict[str, dict[str, object]]] = {}
    gate_failures, regressions = [], []
    with tempfile.TemporaryDirectory() as tmp:
        for method in methods.split(","):
            rows[method] = {}
            for topo in topo_list:
                for seed in seed_list:
                    res = run_synthetic_case(
                        Path(tmp), topo, method=method, preset=syn["preset"], seed=seed,
                        n_cells=syn["n_cells"], n_genes=syn["n_genes"],
                    )  # fmt: skip
                    row = res.row()
                    rows[method][res.case] = row
                    gate_failures += [(method, f) for f in check_gates(row, syn["gates"])]
                    base = baseline.get(method, {}).get(res.case)
                    if base and not write_baseline:
                        regressions += [
                            (method, f) for f in check_regression(row, base, cfg["regression"]["max_drop"])
                        ]
                    typer.echo(
                        f"{method:10s} {res.case:28s} sp={row['pt_spearman_global']:.3f} "
                        f"sp_lin={row['pt_spearman_lineage']:.3f} branch={row['branch_accuracy']:.3f} "
                        f"recall={row['lineage_recall']:.2f} spur={row['n_spurious']} "
                        f"rootX={row['root_pt_excess']:.3f} ({row['runtime_seconds']}s)"
                    )

    for method, f in gate_failures:
        typer.echo(f"gate     {method}: {f}")
    for method, f in regressions:
        typer.echo(f"REGRESSION {method}: {f}", err=True)

    try:
        sha = subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True
        ).stdout.strip()
    except OSError:
        sha = ""
    payload = {"generated": datetime.now(UTC).date().isoformat(), "commit": sha, "methods": rows}
    if out:
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        typer.echo(f"wrote {out}")
    if write_baseline:
        BASELINE.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        typer.echo(f"wrote {BASELINE}")
    if regressions:
        raise typer.Exit(1)


@app.command()
def synth(
    output: Annotated[Path, typer.Argument(help="Where to write the .h5ad file.")],
    topology: Annotated[
        str, typer.Option(help="linear, bifurcation, trifurcation or binary2.")
    ] = "bifurcation",
    preset: Annotated[str, typer.Option(help="easy, realistic or hard.")] = "realistic",
    n_cells: Annotated[int, typer.Option(min=10)] = 3000,
    n_genes: Annotated[int, typer.Option(min=10)] = 2000,
    seed: int = 0,
    store_root: Annotated[bool, typer.Option(help="Also store the true root in uns['iroot'].")] = False,
) -> None:
    """Write a synthetic dataset with a known trajectory (raw counts, no embedding)."""
    from lumencell.bench.synth import PRESETS, TOPOLOGIES, generate

    if topology not in TOPOLOGIES:
        raise typer.BadParameter(f"choose one of {', '.join(sorted(TOPOLOGIES))}", param_hint="--topology")
    if preset not in PRESETS:
        raise typer.BadParameter(f"choose one of {', '.join(PRESETS)}", param_hint="--preset")
    adata = generate(topology, preset, n_cells=n_cells, n_genes=n_genes, seed=seed, store_root=store_root)
    output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output)
    typer.echo(
        f"wrote {output}: {adata.n_obs} cells x {adata.n_vars} genes, {topology} ({preset}), seed {seed}"
    )
