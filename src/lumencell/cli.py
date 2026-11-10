"""Command-line entry point (`lumencell`)."""

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
