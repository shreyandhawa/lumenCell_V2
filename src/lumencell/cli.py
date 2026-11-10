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
