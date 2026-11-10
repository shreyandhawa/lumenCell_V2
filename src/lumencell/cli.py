"""Command-line entry point (`lumencell`)."""

from pathlib import Path
from typing import Annotated

import typer

from lumencell import __version__

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


bench_app = typer.Typer(help="Accuracy benchmarks against known trajectories.", no_args_is_help=True)
app.add_typer(bench_app, name="bench")


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
