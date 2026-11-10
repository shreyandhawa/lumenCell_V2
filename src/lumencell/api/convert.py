"""Conversions between the engine's types and the API's."""

from typing import Any

from lumencell.api import models
from lumencell.engine.roots import CellRoot, ClusterRoot, FileRoot, RootSpec, parse_root
from lumencell.engine.schema import ResolvedRoot

ApiRootSpec = models.CellRoot | models.ClusterRoot | models.FileRoot


def to_api_spec(spec: RootSpec) -> models.CellRoot | models.ClusterRoot | models.FileRoot:
    if isinstance(spec, CellRoot):
        return models.CellRoot(kind="cell", obs_name=spec.obs_name)
    if isinstance(spec, ClusterRoot):
        return models.ClusterRoot(kind="cluster", key=spec.key, value=spec.value)
    return models.FileRoot(kind="file")


def from_api_spec(spec: models.CellRoot | models.ClusterRoot | models.FileRoot) -> RootSpec:
    if isinstance(spec, models.CellRoot):
        return CellRoot(spec.obs_name)
    if isinstance(spec, models.ClusterRoot):
        return ClusterRoot(spec.key, spec.value)
    return FileRoot()


def to_api_root(root: ResolvedRoot | dict[str, Any]) -> models.ResolvedRoot:
    r = root if isinstance(root, dict) else root.__dict__
    warnings = [w if isinstance(w, dict) else w.__dict__ for w in r.get("warnings", ())]
    return models.ResolvedRoot(
        spec=to_api_spec(parse_root(r["spec"])),
        cell_index=r["cell_index"],
        obs_name=r["obs_name"],
        cluster_key=r.get("cluster_key"),
        cluster_label=r.get("cluster_label"),
        resolved_by=r["resolved_by"],
        warnings=[models.Notice(**w) for w in warnings],
    )
