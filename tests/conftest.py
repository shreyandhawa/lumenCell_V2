from dataclasses import dataclass
from pathlib import Path

import pytest

from lumencell.bench.synth import generate
from lumencell.data.loader import Dataset, load
from lumencell.engine.context import InferenceContext
from lumencell.prep.pipeline import PrepResult, run_pipeline
from lumencell.prep.stages import PrepParams


@dataclass(frozen=True)
class Prepared:
    ds: Dataset
    prep: PrepResult
    ctx: InferenceContext


def prepare(
    path: Path, topology: str = "bifurcation", n_cells: int = 1500, seed: int = 0, **kw: object
) -> Prepared:
    adata = generate(topology, "realistic", n_cells=n_cells, n_genes=1000, seed=seed, **kw)  # type: ignore[arg-type]
    adata.write_h5ad(path)
    ds = load(path)
    params = PrepParams()
    prep = run_pipeline(ds, params)
    return Prepared(ds, prep, InferenceContext.from_prep(ds, prep, params))


@pytest.fixture(scope="session")
def bifurcation(tmp_path_factory: pytest.TempPathFactory) -> Prepared:
    """A prepared 1500-cell synthetic bifurcation, shared by engine tests."""
    return prepare(tmp_path_factory.mktemp("bif") / "bif.h5ad", store_root=True)
