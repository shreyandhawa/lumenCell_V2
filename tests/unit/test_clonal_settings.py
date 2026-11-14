"""Pure units of the clonal settings: the checks, the summary, the candidates and the column kind."""

import numpy as np
import pandas as pd
import pytest

from lumencell.api.dataset_view import MAX_CATEGORIES_LISTED, column_kind
from lumencell.api.models import ClonalSettings
from lumencell.trust.clonal import CANDIDATE, check_settings, clonal_candidates, settings_summary
from lumencell.trust.constants import MAX_DONORS

S = ClonalSettings(clone_key="clone", donor_key="donor")


def _obs(clone: list[object], donor: list[object]) -> pd.DataFrame:
    return pd.DataFrame({"clone": pd.Series(clone, dtype=object), "donor": pd.Series(donor, dtype=object)})


def _text(problems: list[str], needle: str) -> bool:
    return any(needle in p for p in problems)


# ---- summary --------------------------------------------------------------------------


def test_summary_pins_every_statistic() -> None:
    sizes = [1, 1, 1, 2, 3, 10]
    clone = [f"c{i}" for i, n in enumerate(sizes) for _ in range(n)]
    obs = _obs(clone, ["d0"] * len(clone))
    got = settings_summary(obs, S)
    assert got.n_cells_labelled == sum(sizes)
    assert got.n_clones == 6
    assert got.n_clones_2plus == 3
    assert got.size_median == 1.5
    assert got.size_p90 == pytest.approx(float(np.percentile(sizes, 90)))
    assert got.size_p99 == pytest.approx(float(np.percentile(sizes, 99)))
    assert got.size_p90 != got.size_p99 and got.size_max == 10
    assert got.n_donors == 1


def test_summary_identity_is_donor_and_label() -> None:
    obs = _obs(["a", "a", "a", "a", "b"], ["d0", "d0", "d1", "d1", "d1"])
    got = settings_summary(obs, S)
    assert got.n_clones == 3  # (d0,a), (d1,a), (d1,b)
    assert got.n_clones_2plus == 2
    assert got.size_max == 2 and got.n_donors == 2


def test_summary_leaves_out_cells_without_a_label_or_a_donor() -> None:
    obs = _obs(
        ["a", "a", "b", None, "c", np.nan],
        ["d0", "d0", np.nan, "d9", "d1", "d2"],
    )
    got = settings_summary(obs, S)
    # (d0,a) twice and (d1,c) once; "b" has no donor; the cells with no label do not count.
    assert got.n_cells_labelled == 3 and got.n_clones == 2 and got.n_clones_2plus == 1
    assert got.n_donors == 4  # d0, d1, d2 and d9 are present; the missing donor is not a donor


def test_summary_refuses_what_check_settings_refuses() -> None:
    with pytest.raises(ValueError, match="no cell has both a clone label and a donor"):
        settings_summary(_obs(["a", "b"], [np.nan, np.nan]), S)


# ---- check_settings -------------------------------------------------------------------


def test_a_valid_obs_has_no_problems() -> None:
    assert check_settings(_obs(["a", "b"], ["d0", "d1"]), S) == []


def test_strata_rules() -> None:
    obs = _obs(["a", "b"], ["d0", "d1"]).assign(t=["x", "y"], u=["x", "y"])
    ok = ClonalSettings(clone_key="clone", donor_key="donor", strata_keys=["t", "u"])
    assert check_settings(obs, ok) == []
    dup = ClonalSettings(clone_key="clone", donor_key="donor", strata_keys=["t", "t"])
    assert _text(check_settings(obs, dup), "distinct")
    on_donor = ClonalSettings(clone_key="clone", donor_key="donor", strata_keys=["donor"])
    assert _text(check_settings(obs, on_donor), "strata column 'donor'")
    on_clone = ClonalSettings(clone_key="clone", donor_key="donor", strata_keys=["clone"])
    assert _text(check_settings(obs, on_clone), "strata column 'clone'")
    single = ClonalSettings(clone_key="clone", single_donor=True, strata_keys=["donor"])
    assert check_settings(obs, single) == []


def test_the_clone_column_must_hold_labels_and_some_cell_must_have_one() -> None:
    all_nan = _obs([np.nan, None], ["d0", "d1"])
    assert _text(check_settings(all_nan, S), "no clone label")
    cat = all_nan.assign(clone=pd.Categorical([None, None], categories=["a", "b"]))
    assert _text(check_settings(cat, S), "no clone label")
    ints = _obs([1, 2], ["d0", "d1"])
    assert _text(check_settings(ints, S), "must hold text labels")
    mixed = _obs(["a", 2], ["d0", "d1"])
    assert _text(check_settings(mixed, S), "must hold text labels")
    numeric = pd.DataFrame({"clone": [1.0, 2.0], "donor": ["d0", "d1"]})
    assert _text(check_settings(numeric, S), "must hold text labels")
    string = _obs(["a", "b"], ["d0", "d1"]).astype({"clone": "string"})
    assert check_settings(string, S) == []
    assert check_settings(_obs(["a", "b"], ["d0", "d1"]).astype({"clone": "category"}), S) == []


@pytest.mark.parametrize(
    "n, ok",
    [(0, False), (1, True), (MAX_DONORS, True), (MAX_DONORS + 1, False)],
)
def test_donor_bounds(n: int, ok: bool) -> None:
    cells = max(n, 1)
    donors: list[object] = [f"d{i}" for i in range(n)] + [np.nan] * (cells - n)
    obs = _obs(["a"] * cells, donors)
    problems = check_settings(obs, S)
    if ok:
        assert problems == []
    else:
        assert _text(problems, f"{n} distinct values")


def test_missing_donors_are_not_counted_as_a_donor() -> None:
    donors: list[object] = [f"d{i}" for i in range(MAX_DONORS)] + [np.nan, None, "", "nan"]
    obs = _obs(["a"] * len(donors), donors)
    assert check_settings(obs, S) == []


def test_unknown_and_equal_columns() -> None:
    obs = _obs(["a"], ["d0"])
    assert _text(check_settings(obs, ClonalSettings(clone_key="zz", donor_key="donor")), "'zz'")
    assert _text(check_settings(obs, ClonalSettings(clone_key="clone", donor_key="zz")), "'zz'")
    assert _text(check_settings(obs, ClonalSettings(clone_key="clone", donor_key="clone")), "different")


def test_labelled_cells_must_have_a_donor() -> None:
    obs = _obs(["a", None], [np.nan, "d0"])
    assert _text(check_settings(obs, S), "no cell has both")


# ---- candidates -----------------------------------------------------------------------


def test_candidates_by_name_in_obs_order() -> None:
    columns = [
        "raw_clonotype_id",
        "donor",
        "xclone_id",
        "AIRR:Clonotype",
        "clone_id2",
        "my_clone_id",
        "airr:raw_clonotype_id",
        "CLONE_ID",
        "clonotype",
        "clone",
    ]
    assert clonal_candidates(columns) == [
        "raw_clonotype_id",
        "AIRR:Clonotype",
        "airr:raw_clonotype_id",
        "CLONE_ID",
        "clonotype",
    ]
    assert CANDIDATE.match("clone_id") and not CANDIDATE.match("clone_id ")
    assert clonal_candidates([]) == []


# ---- column_kind ----------------------------------------------------------------------


def test_column_kind_rules() -> None:
    assert column_kind(pd.Series([1, 2, 3]), ["1", "2", "3"]) == "continuous"
    assert column_kind(pd.Series([1.5, 2.5]), []) == "continuous"
    assert column_kind(pd.Series([True, False]), ["False", "True"]) == "boolean"
    assert column_kind(pd.Series([True, None], dtype="boolean"), ["True"]) == "boolean"
    few = pd.Series(["a", "b"])
    assert column_kind(few, [str(i) for i in range(MAX_CATEGORIES_LISTED)]) == "categorical"
    assert column_kind(few, [str(i) for i in range(MAX_CATEGORIES_LISTED + 1)]) == "string"
    # A string column is never boolean, even one that spells it out.
    assert column_kind(pd.Series(["True", "False"]), ["False", "True"]) == "categorical"
    cat = pd.Series(
        pd.Categorical(["a", "b"], categories=["a", "b", *(str(i) for i in range(MAX_CATEGORIES_LISTED + 5))])
    )
    assert column_kind(cat, ["a", "b"]) == "categorical"
    assert column_kind(pd.Series(pd.Categorical([True, False])), ["False", "True"]) == "categorical"
