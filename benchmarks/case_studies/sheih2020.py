"""The Sheih 2020 case-study numbers: every method on each patient alone, plus the pooled diagnostic.

Wiring only: the declaration is [case_study.sheih2020] in benchmarks/thresholds.toml, the data is the
derived file `lumencell bench fetch sheih2020` builds, and every number comes from
lumencell.bench.case_study. Takes about 15 minutes on the development machine (4 patients x 3
methods with stability on each, then the pooled run; the 2026-11-20 run took 14).
Run as: uv run python -I benchmarks/case_studies/sheih2020.py [--out PATH] [--workdir DIR]
"""

import argparse
import subprocess
import sys
import tempfile
from datetime import UTC, datetime
from pathlib import Path

from lumencell.bench.case_study import failed_methods, load_declaration, run_units, write_report
from lumencell.bench.datasets import CASE_STUDIES
from lumencell.bench.download import is_cached, local_path
from lumencell.bench.gates import REPO


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", type=Path, default=None, help="results file (default: dated, under results/)")
    ap.add_argument(
        "--workdir", type=Path, default=None, help="unit files and prep caches (default: a temp dir)"
    )
    args = ap.parse_args()

    item = CASE_STUDIES["sheih2020"]
    if not is_cached(item):
        raise SystemExit("sheih2020 is not built; run: uv run lumencell bench fetch sheih2020")
    decl = load_declaration("sheih2020")
    commit = subprocess.run(
        ["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, cwd=REPO, check=True
    ).stdout.strip()
    date = datetime.now(UTC).date().isoformat()
    out = args.out or REPO / "benchmarks" / "case_studies" / "results" / f"sheih2020-{date}.json"

    if args.workdir is not None:
        report = run_units(local_path(item), decl, args.workdir)
    else:
        with tempfile.TemporaryDirectory(prefix="lumencell-case-") as tmp:
            report = run_units(local_path(item), decl, Path(tmp))
    write_report(report, out, commit)
    print(f"wrote {out}")
    failed = failed_methods(report)
    for unit, method, exc in failed:
        print(f"FAILED: {unit} {method} ({exc})", file=sys.stderr)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
