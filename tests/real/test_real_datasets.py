"""Checks on the downloaded real datasets; skipped when they are not cached."""

import anndata as ad
import pytest

from lumencell.bench.datasets import DATASETS, is_cached, local_path, sha256_of, truth_for

pytestmark = [
    pytest.mark.real_data,
    # These files were written by 2019-era anndata; reading them is fine.
    pytest.mark.filterwarnings("ignore::anndata.OldFormatWarning"),
    pytest.mark.filterwarnings("ignore:Moving element from .uns:FutureWarning"),
]


@pytest.mark.parametrize("dataset_id", sorted(DATASETS))
def test_cached_dataset_matches_its_checksum_and_truth(dataset_id: str) -> None:
    ds = DATASETS[dataset_id]
    if not is_cached(ds):
        pytest.skip(f"run: lumencell bench fetch {dataset_id}")
    assert sha256_of(local_path(ds)) == ds.sha256
    obs = ad.read_h5ad(local_path(ds), backed="r").obs
    truth = truth_for(ds, obs)
    assert truth.root.sum() > 0
    for name, mask in truth.lineages.items():
        assert mask.sum() > 50, f"{dataset_id}/{name} has almost no cells"
