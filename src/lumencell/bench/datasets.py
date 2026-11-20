"""Real datasets with a known trajectory, for the release benchmark.

They are never downloaded implicitly: ``lumencell bench fetch <id>`` shows the
source and size and asks first. Files go to the user cache directory
(``LUMENCELL_DATASETS_DIR`` overrides it), are checked against the size and
SHA-256 recorded here, and are never modified.

Truth comes from published annotations, not from any method:

* ``order`` is a stage ordinal along each expected lineage (pancreas, bone
  marrow) or real embryo time (gastrulation);
* each expected lineage is the set of cells whose annotation lies on it;
* the root is the annotated starting population, and it is also the root the
  methods are given, the way a user would name it.

Methods never see these annotations: they cluster with Leiden. The
annotation column is only used to name lineages and, in the benchmark, to say
where the trajectory starts.
"""

import re
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

from lumencell.bench.builders import sheih2020
from lumencell.bench.download import (
    USER_AGENT,
    Build,
    FetchError,
    Source,
    datasets_dir,
    fetch,
    is_cached,
    local_path,
    sha256_of,
    source_path,
)
from lumencell.bench.truth import Truth

__all__ = [
    "CASE_STUDIES",
    "DATASETS",
    "USER_AGENT",
    "Build",
    "CaseStudyDataset",
    "FetchError",
    "RealDataset",
    "Source",
    "datasets_dir",
    "fetch",
    "is_cached",
    "local_path",
    "sha256_of",
    "source_path",
    "truth_for",
]


@dataclass(frozen=True)
class RealDataset:
    id: str
    title: str
    citation: str
    sources: tuple[Source, ...]
    licence: str  # what the page hosting the files states, or that it states none
    label_key: str  # annotation column holding the stages
    root_label: str
    # Expected lineages: stage labels in order. Gated lineages count towards the
    # gates; the others are reported only.
    lineages: dict[str, list[str]]
    gated: list[str]
    build: Build | None = None  # turns the sources into one derived .h5ad
    time_key: str | None = None  # obs column with real time, used as the order when present
    notes: str = ""
    extra: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.build is None and len(self.sources) != 1:
            raise ValueError(f"{self.id}: a dataset without a build step has exactly one source")
        names = [src.filename for src in self.sources]
        if not names or len(set(names)) != len(names):
            raise ValueError(f"{self.id}: sources need unique filenames")
        if not self.licence.strip():
            raise ValueError(f"{self.id}: a licence is required")

    @property
    def root_spec(self) -> str:
        return f"cluster:{self.label_key}={self.root_label}"


DATASETS: dict[str, RealDataset] = {
    d.id: d
    for d in (
        RealDataset(
            id="pancreas",
            title="Mouse pancreatic endocrinogenesis, E15.5",
            citation="Bastidas-Ponce et al. 2019, Development 146:dev173849",
            sources=(
                Source(
                    url="https://raw.githubusercontent.com/theislab/scvelo_notebooks/master/data/Pancreas/endocrinogenesis_day15.h5ad",
                    filename="pancreas_endocrinogenesis_day15.h5ad",
                    size=52_477_898,
                    sha256="9e3e459eca00ba06b496ec80def32941b5b2889918720e3e7aa6ffb811fbe7c6",
                ),
            ),
            licence=(
                "no licence stated on https://github.com/theislab/scvelo_notebooks, read 2026-11-19; "
                "cite Bastidas-Ponce et al. 2019"
            ),
            label_key="clusters",
            root_label="Ductal",
            lineages={
                "Alpha": ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Alpha"],
                "Beta": ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Beta"],
                "Delta": ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Delta"],
                "Epsilon": ["Ductal", "Ngn3 low EP", "Ngn3 high EP", "Pre-endocrine", "Epsilon"],
            },
            gated=["Alpha", "Beta"],
            notes="Tuning set: defaults may be adjusted on this dataset.",
        ),
        RealDataset(
            id="bonemarrow",
            title="Human CD34+ bone marrow",
            citation="Setty et al. 2019, Nature Biotechnology 37:451-460",
            sources=(
                Source(
                    url="https://ndownloader.figshare.com/files/27686835",
                    filename="human_cd34_bone_marrow.h5ad",
                    size=169_575_532,
                    sha256="12222efe4a9dd5916fa279860a2a4fdec383bd2e0db0249a1cd18c549c4a4c2c",
                ),
            ),
            licence=(
                "CC BY 4.0, as stated on "
                "https://figshare.com/articles/dataset/human_cd34_bone_marrow_h5ad/14465559 "
                "(figshare file 27686835), read 2026-11-19; cite Setty et al. 2019"
            ),
            label_key="clusters",
            root_label="HSC_1",
            lineages={
                "Ery": ["HSC_1", "Ery_1", "Ery_2"],
                "Mono": ["HSC_1", "Mono_1", "Mono_2"],
            },
            gated=["Ery", "Mono"],
            notes="Held out: never used to choose defaults.",
        ),
        RealDataset(
            id="gastrulation_erythroid",
            title="Mouse gastrulation, erythroid lineage",
            citation="Pijuan-Sala et al. 2019, Nature 566:490-495",
            sources=(
                Source(
                    url="https://ndownloader.figshare.com/files/27686871",
                    filename="erythroid_lineage.h5ad",
                    size=485_759_710,
                    sha256="0af49504137dc8dc930a19e3534243a121da1ddddbe28ebc5fab0470294d60d4",
                ),
            ),
            licence=(
                "CC BY 4.0, as stated on "
                "https://figshare.com/articles/dataset/Gastrulation_erythroid_lineage/14465562 "
                "(figshare file 27686871), read 2026-11-19; cite Pijuan-Sala et al. 2019"
            ),
            label_key="celltype",
            root_label="Blood progenitors 1",
            lineages={
                "Erythroid": [
                    "Blood progenitors 1",
                    "Blood progenitors 2",
                    "Erythroid1",
                    "Erythroid2",
                    "Erythroid3",
                ]
            },
            gated=["Erythroid"],
            time_key="stage",
            notes="Held out. Order is real embryo time (E7.0-E8.5), not annotation.",
        ),
    )
}


@dataclass(frozen=True)
class CaseStudyDataset:
    """A case study: fetched and built like a benchmark dataset, but not in ``DATASETS``.

    It has no truth. ``bench fetch all`` skips it; it is fetched by its id.
    """

    id: str
    title: str
    citation: str
    sources: tuple[Source, ...]
    licence: str
    build: Build
    notes: str = ""


# build.params equal the declaration's chains, cdr3, filters and max_per_chain (thresholds.toml,
# [case_study.sheih2020]); tests/unit/test_sheih2020.py pins that.
CASE_STUDIES: dict[str, CaseStudyDataset] = {
    "sheih2020": CaseStudyDataset(
        id="sheih2020",
        title="CD19 CAR-T cells in four patients (CLL-1, CLL-2, NHL-6, NHL-7), scRNA with TCR clonotypes",
        citation="Sheih et al. 2020, Nat Commun 11:219, doi:10.1038/s41467-019-13880-1",
        sources=sheih2020.SOURCES,
        licence=sheih2020.LICENCE,
        build=Build(
            builder=sheih2020.build,
            version=sheih2020.BUILD_VERSION,
            params={"chains": "paired", "cdr3": "nt", "filters": [], "max_per_chain": 2},
        ),
        notes="Case study (evaluation only): the declaration is [case_study.sheih2020].",
    ),
}


# ---- truth --------------------------------------------------------------------------------

_NUMBER = re.compile(r"[-+]?\d*\.?\d+")


def truth_for(ds: RealDataset, obs: pd.DataFrame) -> Truth:
    labels = obs[ds.label_key].astype(str).to_numpy()
    lineages = {name: np.isin(labels, stages) for name, stages in ds.lineages.items() if name in ds.gated}
    if ds.time_key is not None:
        order = np.array(
            [float(m.group()) if (m := _NUMBER.search(str(v))) else np.nan for v in obs[ds.time_key]],
            dtype=np.float64,
        )
    else:
        order = np.full(len(labels), np.nan)
        for stages in ds.lineages.values():
            for i, stage in enumerate(stages):
                order[labels == stage] = i
    root = labels == ds.root_label
    terminals = {name: labels == ds.lineages[name][-1] for name in ds.gated}
    return Truth(order=order, lineages=lineages, root=root, meta={"dataset": ds.id}, terminals=terminals)
