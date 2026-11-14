"""Saving and loading trajectory runs.

A run lives in ``<dataset cache dir>/runs/<run_id>/``::

    result.json     graph, lineages, root, parameters, warnings, provenance
    cells.lcb       every per-cell array, already in the wire format, so the API
                    serves GET /trajectories/{id}/cells straight from disk
    weights.npz     the milestone weight matrix used to draw on embeddings
    curve_nodes.npz cell weights of a fitted curve's nodes (only for methods that fit one)

The run id hashes everything that determines the result, so asking for the same
run twice finds the existing one instead of recomputing it.
"""

import hashlib
import json
import logging
import os
import re
import shutil
import sys
import tempfile
from collections.abc import Sequence
from dataclasses import asdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np
from scipy import sparse

from lumencell import __version__
from lumencell.api import binary
from lumencell.clonal.isolation import ISOLATION_VERSION
from lumencell.engine.schema import (
    Divergence,
    Lineage,
    Milestone,
    Notice,
    ResolvedRoot,
    Segment,
    TrajectoryResult,
)

log = logging.getLogger(__name__)

# 2: result.json records the prep key, the clone-isolation version and the obs keys used.
STORE_VERSION = 2

# A trust file is named <kind>-<16 hex>; the name is checked before it reaches a path.
TRUST_NAME = re.compile(r"[a-z][a-z_]{0,31}-[0-9a-f]{16}")


def run_id(
    method: str,
    params: dict[str, Any],
    root_spec: str,
    terminals: list[str] | None,
    prep_key: str,
    clonal_columns: Sequence[str] = (),
) -> str:
    payload: dict[str, Any] = {
        "store": STORE_VERSION,
        "lumencell": __version__,
        "method": method,
        "params": params,
        "root": root_spec,
        "terminals": sorted(terminals) if terminals else None,
        "prep": prep_key,
    }
    # Naming clone or donor columns can change the detected label key, and with it the lineage
    # names, so the run depends on them. Left out when empty: every earlier run id is unchanged.
    if clonal_columns:
        payload["clonal_columns"] = sorted(clonal_columns)
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def cells_payload(result: TrajectoryResult) -> bytes:
    """The LCB1 buffer served by GET /trajectories/{id}/cells."""
    arrays: dict[str, tuple[str, Any]] = {
        "pseudotime": ("f32", result.pseudotime),
        "membership": ("u8", result.membership.astype(np.uint8)),
        "primary_lineage": ("i16", result.primary_lineage),
        "segment": ("i32", result.segment),
        "progress": ("f32", result.progress),
    }
    if result.fate is not None:
        arrays["fate"] = ("f32", result.fate)
    for name, values in result.extras.items():
        arrays[f"extras.{name}"] = ("f32", values)
    return binary.encode(arrays, {"n_obs": result.n_obs, "n_lineages": len(result.lineages)})


def _graph_json(
    result: TrajectoryResult,
    rid: str,
    created_at: str,
    terminals: list[str] | None,
    prep_key: str | None,
) -> dict[str, Any]:
    root = asdict(result.root)
    return {
        "store_version": STORE_VERSION,
        "run_id": rid,
        "created_at": created_at,
        # Provenance: whether the terminals were named by the user or found by the method.
        "terminals_source": "user" if terminals else "auto",
        # Provenance: the prep the run was computed on, how clone columns were kept out of
        # inference, and which obs columns steered it.
        "prep_key": prep_key,
        "inference_isolation": ISOLATION_VERSION,
        "obs_keys": result.obs_keys,
        "method": result.method,
        "method_version": result.method_version,
        "lumencell_version": __version__,
        "topology_kind": result.topology_kind,
        "milestones": [asdict(m) for m in result.milestones],
        "segments": [asdict(s) for s in result.segments],
        "lineages": [asdict(lin) for lin in result.lineages],
        "divergence": [asdict(d) for d in result.divergence],
        "root": root,
        "inference_rep": result.inference_rep,
        "inference_dims": result.inference_dims,
        "params": result.params,
        "warnings": [asdict(w) for w in result.warnings],
        "runtime_seconds": result.runtime_seconds,
        "n_obs": result.n_obs,
        "n_unassigned": result.n_unassigned,
        "has_fate": result.fate is not None,
        "extras": sorted(result.extras),
        "curve_paths": result.curve_paths,
    }


class RunStore:
    def __init__(self, dataset_cache_dir: Path) -> None:
        self.dir = dataset_cache_dir / "runs"

    def path(self, rid: str) -> Path:
        return self.dir / rid

    def exists(self, rid: str) -> bool:
        """A completed run written by this store version. Any other is absent, so it is recomputed."""
        if not (self.path(rid) / ".complete").exists():
            return False
        try:
            return bool(self.graph(rid).get("store_version") == STORE_VERSION)
        except (OSError, ValueError, AttributeError):
            return False

    def save(
        self,
        rid: str,
        result: TrajectoryResult,
        terminals: list[str] | None = None,
        prep_key: str | None = None,
    ) -> None:
        """Save a run. ``terminals`` are the terminal labels the run was asked for, if any;
        ``prep_key`` is the prep cache key the run was computed on."""
        self.dir.mkdir(parents=True, exist_ok=True)
        final = self.path(rid)
        tmp = Path(tempfile.mkdtemp(prefix=f".{rid}-", dir=self.dir))
        try:
            created = datetime.now(UTC).isoformat(timespec="seconds")
            (tmp / "result.json").write_text(
                json.dumps(_graph_json(result, rid, created, terminals, prep_key), indent=2),
                encoding="utf-8",
            )
            (tmp / "cells.lcb").write_bytes(cells_payload(result))
            sparse.save_npz(
                tmp / "weights.npz", sparse.csr_matrix(result.milestone_weights), compressed=False
            )
            if result.curve_nodes is not None:
                sparse.save_npz(
                    tmp / "curve_nodes.npz", sparse.csr_matrix(result.curve_nodes), compressed=False
                )
            (tmp / ".complete").write_text(created, encoding="utf-8")
            old: Path | None = None
            if final.exists():  # set the old run aside, so a failure below can put it back
                old = Path(tempfile.mkdtemp(prefix=f".{rid}-", dir=self.dir))
                old.rmdir()
                final.rename(old)
            try:
                if old is not None and (old / "trust").is_dir():  # trust outlives a re-save
                    (old / "trust").rename(tmp / "trust")
                tmp.rename(final)
            except BaseException:
                if old is not None:
                    if (tmp / "trust").is_dir() and not (old / "trust").exists():
                        (tmp / "trust").rename(old / "trust")
                    if not final.exists():
                        old.rename(final)
                raise
            if old is not None:
                shutil.rmtree(old, ignore_errors=True)
        except BaseException:
            shutil.rmtree(tmp, ignore_errors=True)
            raise

    def trust_path(self, rid: str, name: str) -> Path:
        if not TRUST_NAME.fullmatch(name):
            raise ValueError(f"not a trust file name: {name!r}")
        return self.path(rid) / "trust" / f"{name}.json"

    def write_trust(self, rid: str, name: str, payload: dict[str, Any]) -> Path:
        """Write one trust result into a complete run, by atomic rename."""
        target = self.trust_path(rid, name)
        if not (self.path(rid) / ".complete").exists():
            raise FileNotFoundError(f"no complete run {rid!r}")
        text = json.dumps(payload, indent=2, allow_nan=False)
        target.parent.mkdir(exist_ok=True)
        fd, tmp = tempfile.mkstemp(prefix=f".{name}-", suffix=".tmp", dir=target.parent)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as f:
                f.write(text)
            os.replace(tmp, target)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise
        return target

    def read_trust(self, rid: str, name: str) -> dict[str, Any] | None:
        """A trust result, or None when missing, unreadable or not a JSON object."""
        path = self.trust_path(rid, name)
        try:
            data = json.loads(path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            return None
        except (OSError, ValueError) as e:
            print(f"lumencell: ignoring unreadable trust file {path}: {e}", file=sys.stderr)
            return None
        if not isinstance(data, dict):
            print(f"lumencell: ignoring trust file {path}: not a JSON object", file=sys.stderr)
            return None
        return data

    def graph(self, rid: str) -> dict[str, Any]:
        data: dict[str, Any] = json.loads((self.path(rid) / "result.json").read_text(encoding="utf-8"))
        return data

    def cells_bytes(self, rid: str) -> bytes:
        return (self.path(rid) / "cells.lcb").read_bytes()

    def update_graph(self, rid: str, graph: dict[str, Any]) -> None:
        p = self.path(rid) / "result.json"
        tmp = p.with_suffix(".json.tmp")
        tmp.write_text(json.dumps(graph, indent=2), encoding="utf-8")
        tmp.replace(p)

    def load(self, rid: str) -> TrajectoryResult:
        if not self.exists(rid):
            raise FileNotFoundError(f"no run {rid!r} saved by store version {STORE_VERSION}")
        g = self.graph(rid)
        cells = binary.decode(self.cells_bytes(rid)).arrays
        fate = np.array(cells["fate"]) if "fate" in cells else None
        extras = {k.removeprefix("extras."): np.array(v) for k, v in cells.items() if k.startswith("extras.")}
        root = dict(g["root"])
        root["warnings"] = tuple(Notice(**w) for w in root.get("warnings", ()))
        return TrajectoryResult(
            method=g["method"],
            method_version=g["method_version"],
            topology_kind=g["topology_kind"],
            milestones=[Milestone(**m) for m in g["milestones"]],
            segments=[Segment(**s) for s in g["segments"]],
            lineages=[Lineage(**lin) for lin in g["lineages"]],
            divergence=[Divergence(**d) for d in g["divergence"]],
            pseudotime=np.array(cells["pseudotime"]),
            membership=np.array(cells["membership"]).astype(bool),
            segment=np.array(cells["segment"]),
            progress=np.array(cells["progress"]),
            root=ResolvedRoot(**root),
            inference_rep=g["inference_rep"],
            inference_dims=g["inference_dims"],
            params=g["params"],
            milestone_weights=sparse.csr_matrix(sparse.load_npz(self.path(rid) / "weights.npz")),
            fate=fate,
            extras=extras,
            warnings=[Notice(**w) for w in g["warnings"]],
            runtime_seconds=g["runtime_seconds"],
            curve_nodes=(
                sparse.csr_matrix(sparse.load_npz(self.path(rid) / "curve_nodes.npz"))
                if (self.path(rid) / "curve_nodes.npz").exists()
                else None
            ),
            curve_paths=g.get("curve_paths"),
            obs_keys=list(g["obs_keys"]),
        )

    def list(self) -> list[dict[str, Any]]:
        """Summaries of every completed run this version can read, newest first.

        A run folder that is damaged, or was written by another store version,
        is left out rather than failing the whole list.
        """
        if not self.dir.exists():
            return []
        out = []
        other_version = 0
        for d in self.dir.iterdir():
            if d.name.startswith(".") or not (d / ".complete").exists():
                continue
            try:
                g = self.graph(d.name)
                if g.get("store_version") != STORE_VERSION:
                    other_version += 1
                    continue
                datetime.fromisoformat(g["created_at"])
                summary = {
                    "run_id": str(g["run_id"]),
                    "method": str(g["method"]),
                    "created_at": str(g["created_at"]),
                    "n_lineages": len(g["lineages"]),
                }
            except (OSError, ValueError, KeyError, TypeError, AttributeError):
                continue
            out.append(summary)
        if other_version:
            log.warning(
                "%d saved run(s) in %s were written by another store version and are not listed; "
                "running them again recomputes them",
                other_version,
                self.dir,
            )
        return sorted(out, key=lambda r: r["created_at"], reverse=True)
