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

import hashlib
import os
import re
import ssl
import tempfile
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

import certifi
import numpy as np
import pandas as pd
import platformdirs

from lumencell.bench.truth import Truth

USER_AGENT = "lumencell-bench"


@dataclass(frozen=True)
class RealDataset:
    id: str
    title: str
    citation: str
    url: str
    filename: str
    size: int
    sha256: str | None
    label_key: str  # annotation column holding the stages
    root_label: str
    # Expected lineages: stage labels in order. Gated lineages count towards the
    # gates; the others are reported only.
    lineages: dict[str, list[str]]
    gated: list[str]
    time_key: str | None = None  # obs column with real time, used as the order when present
    notes: str = ""
    extra: dict[str, str] = field(default_factory=dict)

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
            url="https://raw.githubusercontent.com/theislab/scvelo_notebooks/master/data/Pancreas/endocrinogenesis_day15.h5ad",
            filename="pancreas_endocrinogenesis_day15.h5ad",
            size=52_477_898,
            sha256="9e3e459eca00ba06b496ec80def32941b5b2889918720e3e7aa6ffb811fbe7c6",
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
            url="https://ndownloader.figshare.com/files/27686835",
            filename="human_cd34_bone_marrow.h5ad",
            size=169_575_532,
            sha256="12222efe4a9dd5916fa279860a2a4fdec383bd2e0db0249a1cd18c549c4a4c2c",
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
            url="https://ndownloader.figshare.com/files/27686871",
            filename="erythroid_lineage.h5ad",
            size=485_759_710,
            sha256="0af49504137dc8dc930a19e3534243a121da1ddddbe28ebc5fab0470294d60d4",
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


def datasets_dir() -> Path:
    override = os.environ.get("LUMENCELL_DATASETS_DIR")
    return Path(override) if override else Path(platformdirs.user_cache_dir("lumencell")) / "datasets"


def local_path(ds: RealDataset) -> Path:
    return datasets_dir() / ds.filename


def is_cached(ds: RealDataset) -> bool:
    p = local_path(ds)
    return p.exists() and p.stat().st_size == ds.size


def sha256_of(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class FetchError(RuntimeError):
    pass


def fetch(ds: RealDataset, progress: Callable[[int, int], None] | None = None) -> Path:
    """Download ``ds`` (the caller has already asked the user). Returns the local path."""
    target = local_path(ds)
    if is_cached(ds):
        return target
    target.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(ds.url, headers={"User-Agent": USER_AGENT})
    fd, tmp_name = tempfile.mkstemp(prefix=f".{ds.filename}.", dir=target.parent)
    tmp = Path(tmp_name)
    try:
        h = hashlib.sha256()
        done = 0
        # certifi's CA bundle: some Python builds (e.g. python.org on macOS) ship without CA certificates.
        context = ssl.create_default_context(cafile=certifi.where())
        with os.fdopen(fd, "wb") as out, urllib.request.urlopen(request, timeout=60, context=context) as res:
            while chunk := res.read(1 << 20):
                out.write(chunk)
                h.update(chunk)
                done += len(chunk)
                if progress:
                    progress(done, ds.size)
        if done != ds.size:
            raise FetchError(f"{ds.id}: expected {ds.size:,} bytes, got {done:,}")
        digest = h.hexdigest()
        if ds.sha256 and digest != ds.sha256:
            raise FetchError(f"{ds.id}: checksum mismatch ({digest} != {ds.sha256})")
        tmp.rename(target)
        return target
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise


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
