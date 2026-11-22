"""Time Palantir and scFates on 100k synthetic cells, once, and record the gates.

Follows the recipe in docs/backend-robustness-2.md (Issue 2): `lumencell synth`,
`lumencell prepare`, then `lumencell run` for each method, with wall time and
peak RSS from /usr/bin/time -l. A run that passes the cap is killed and recorded
as a miss. Budgets are fixed here and are not changed to fit results.

    caffeinate -i env PYTHONPATH=src LUMENCELL_BIN=<venv>/bin/lumencell \
        python benchmarks/perf/timing_100k.py
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import platform
import re
import signal
import subprocess
import time
from pathlib import Path

BUDGET_S = {"palantir": 600, "scfates_ppt": 600}
CAP_S = 3600  # longer than this is recorded as a miss; the run is not repeated


def timed(cmd: list[str], env: dict[str, str], cap: float | None = None) -> dict[str, object]:
    """Run `cmd` under /usr/bin/time -l; wall seconds, peak RSS in MB, status."""
    start = time.monotonic()
    proc = subprocess.Popen(
        ["/usr/bin/time", "-l", *cmd],
        env=env,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        start_new_session=True,
    )
    status: object
    try:
        out, err = proc.communicate(timeout=cap)
        status = proc.returncode
    except subprocess.TimeoutExpired:
        os.killpg(proc.pid, signal.SIGKILL)
        out, err = proc.communicate()
        status = f"killed at the {int(cap or 0)} s cap"
    wall = time.monotonic() - start
    rss = re.search(r"(\d+)\s+maximum resident set size", err or "")
    return {
        "wall_s": round(wall, 1),
        "peak_rss_mb": round(int(rss.group(1)) / 2**20) if rss else None,
        "status": status,
        "stdout": out or "",
        "stderr_tail": "\n".join((err or "").splitlines()[-8:]),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--workdir", default=".lumencell/perf", help="Dataset and caches (kept out of git).")
    parser.add_argument("--n-cells", type=int, default=100_000)
    parser.add_argument("--methods", default="palantir,scfates_ppt")
    args = parser.parse_args()

    today = datetime.date.today().isoformat()
    work = Path(args.workdir).resolve()
    work.mkdir(parents=True, exist_ok=True)
    data = work / f"synth{args.n_cells // 1000}k.h5ad"
    env = {
        **os.environ,
        "LUMENCELL_CACHE_DIR": str(work / "cache"),
        "LUMENCELL_DATA_DIR": str(work / "data"),
    }
    cli = [os.environ.get("LUMENCELL_BIN", "lumencell")]
    out = Path("benchmarks/results") / f"perf-100k-{today}.json"
    prior = json.loads(out.read_text()) if out.exists() else {}
    result: dict[str, object] = {
        "date": today,
        "n_cells": args.n_cells,
        "machine": f"{platform.platform()}, {os.cpu_count()} logical CPUs",
        "budgets_s": BUDGET_S,
        "cap_s": CAP_S,
        "root": prior.get("root"),
        "steps": prior.get("steps", {}),
        "methods": prior.get("methods", {}),
    }
    steps: dict[str, object] = result["steps"]  # type: ignore[assignment]
    methods: dict[str, object] = result["methods"]  # type: ignore[assignment]

    # A rerun of some methods keeps the earlier synth, prepare and root steps.
    if "prepare" not in steps:
        if not data.exists():
            subprocess.run([*cli, "synth", str(data), "--n-cells", str(args.n_cells)], env=env, check=True)
        steps["prepare"] = timed([*cli, "prepare", str(data)], env)
        print("prepare:", {k: v for k, v in steps["prepare"].items() if k != "stdout"}, flush=True)  # type: ignore[union-attr]
        roots = timed([*cli, "roots", str(data)], env)
        steps["roots"] = {k: v for k, v in roots.items() if k != "stdout"}
        lines = [ln for ln in str(roots["stdout"]).splitlines() if ln.strip()]
        result["root"] = None if not lines or lines[0].startswith("no suggestions") else lines[0].split()[0]
    root = result["root"]
    if root is None:
        print("no root suggestion; nothing to time", flush=True)
    else:
        root = str(root)
        for method in args.methods.split(","):
            r = timed([*cli, "run", str(data), "--root", root, "--method", method], env, cap=CAP_S)
            methods[method] = {k: v for k, v in r.items() if k != "stdout"}
            wall = float(r["wall_s"])  # type: ignore[arg-type]
            verdict = "pass" if r["status"] == 0 and wall < BUDGET_S[method] else "fail"
            methods[method]["gate"] = verdict  # type: ignore[index]
            print(f"{method}: {verdict}, wall {wall} s, status {r['status']}", flush=True)

    out.write_text(json.dumps(result, indent=2) + "\n")
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
