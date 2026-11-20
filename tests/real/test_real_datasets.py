"""Checks on the downloaded real datasets; skipped when they are not cached."""

import anndata as ad
import pytest

from lumencell.bench.datasets import DATASETS, is_cached, local_path, sha256_of, source_path, truth_for

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
    for s in ds.sources:
        assert sha256_of(source_path(s)) == s.sha256
    obs = ad.read_h5ad(local_path(ds), backed="r").obs
    truth = truth_for(ds, obs)
    assert truth.root.sum() > 0
    for name, mask in truth.lineages.items():
        assert mask.sum() > 50, f"{dataset_id}/{name} has almost no cells"


def test_case_study_files_match_their_checksums() -> None:
    from lumencell.bench.datasets import CASE_STUDIES

    for cs in CASE_STUDIES.values():
        if not all(source_path(s).exists() for s in cs.sources):
            pytest.skip(f"run: lumencell bench fetch {cs.id}")
        for s in cs.sources:
            assert source_path(s).stat().st_size == s.size
            assert sha256_of(source_path(s)) == s.sha256
