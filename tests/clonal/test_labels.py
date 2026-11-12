import numpy as np
import pandas as pd

from lumencell.clonal.labels import encode_clones, is_missing


def objects(values: list[object]) -> np.ndarray:
    return np.array(values, dtype=object)


def test_missing_values_are_unlabelled() -> None:
    labels = objects([None, np.nan, "", "nan", "None", pd.NA, "c1"])
    assert is_missing(labels).tolist() == [True] * 6 + [False]
    codes = encode_clones(labels, None)
    assert (codes.clone[:6] == -1).all()
    assert codes.clone[6] >= 0


def test_clone_identity_is_donor_and_label() -> None:
    codes = encode_clones(objects(["a", "a", "a", "b"]), objects(["d1", "d2", "d1", "d1"]))
    assert codes.clone[0] == codes.clone[2] != codes.clone[1]
    assert len(set(codes.clone.tolist())) == 3
    assert (codes.clone >= 0).all()


def test_missing_donor_drops_the_cell_and_counts_it() -> None:
    codes = encode_clones(objects(["a", "a", "a"]), objects([None, "d1", "d1"]))
    assert codes.clone[0] == -1
    assert codes.donor[0] == -1
    assert codes.n_missing_donor == 1
    assert codes.clone[1] == codes.clone[2] >= 0


def test_an_unlabelled_cell_without_a_donor_is_not_counted_as_missing_a_donor() -> None:
    codes = encode_clones(objects([None, "a", "a"]), objects([None, "d1", "d1"]))
    assert codes.n_missing_donor == 0


def test_labels_are_compared_as_whole_strings() -> None:
    codes = encode_clones(objects(["a", "a\x00", 1, "1"]), None)
    assert codes.clone[0] != codes.clone[1]
    assert codes.clone[2] == codes.clone[3]  # a label is its string


def test_unused_categories_are_not_clones() -> None:
    labels = pd.Series(pd.Categorical(["a", "a"], categories=["a", "b", "zz"])).to_numpy(dtype=object)
    codes = encode_clones(labels, None)
    assert set(codes.clone.tolist()) == {0}


def test_strata_combine_with_donor_and_missing_is_its_own_level() -> None:
    codes = encode_clones(
        objects(["a", "a", "a", "a"]),
        objects(["d1", "d1", "d1", "d2"]),
        [objects(["t1", "t1", None, "t1"])],
    )
    assert codes.stratum is not None
    assert codes.stratum[0] == codes.stratum[1]
    assert codes.stratum[2] not in (codes.stratum[0], codes.stratum[1])
    assert codes.stratum[3] != codes.stratum[0]


def test_missing_stratum_values_share_one_level_apart_from_every_string() -> None:
    # is_missing flags None and the string "None" alike, so both are the one missing level ...
    codes = encode_clones(objects(["a", "a"]), None, [objects([None, "None"])])
    assert codes.stratum is not None and codes.stratum[0] == codes.stratum[1]
    # ... and a missing value never shares a level with a present string.
    codes = encode_clones(objects(["a", "a"]), None, [objects([None, "x"])])
    assert codes.stratum is not None and codes.stratum[0] != codes.stratum[1]


def test_no_strata_means_no_stratum_codes() -> None:
    assert encode_clones(objects(["a", "a"]), None).stratum is None


def test_no_donor_column_means_one_donor() -> None:
    codes = encode_clones(objects(["a", "b", None]), None)
    assert codes.donor.tolist() == [0, 0, 0]
    assert codes.n_missing_donor == 0


def test_codes_are_int64() -> None:
    codes = encode_clones(objects(["a", None]), objects(["d1", "d1"]), [objects(["t", "t"])])
    assert codes.stratum is not None
    assert codes.clone.dtype == codes.donor.dtype == codes.stratum.dtype == np.int64
