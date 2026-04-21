#!/usr/bin/env python3
"""Expand grids.yaml sweep entries into one-per-line CSV
(sweep_name,nx,nz,lin_solver,n_threads) on stdout.

Uses a deliberately tiny parser so the harness has no external YAML
dependency; the structure of grids.yaml is fixed (see schema_version).

Invocation: `GRIDS_FILE=... python3 _parse_grids.py`
"""
from __future__ import annotations

import os
import re
import sys


def parse(path: str) -> list[str]:
    with open(path) as f:
        txt = f.read()
    rows: list[str] = []
    current: dict | None = None

    def emit() -> None:
        if not current:
            return
        skips = set()
        for sk in current.get("skip", []):
            skips.add((sk.get("lin_solver"), sk.get("nx"), sk.get("nz")))
        for g in current.get("grids", []):
            for s in current.get("lin_solvers", []):
                for t in current.get("n_threads", []):
                    if (s, g[0], g[1]) in skips:
                        continue
                    rows.append(f"{current['name']},{g[0]},{g[1]},{s},{t}")

    # Track whether we are currently parsing a `skip:` block (list of dicts
    # under the current sweep, same indent as `grids:`).
    in_skip = False
    for line in txt.splitlines():
        if re.match(r"^  - name:", line):
            emit()
            in_skip = False
            current = {
                "name": line.split(":", 1)[1].strip(),
                "grids": [],
                "lin_solvers": [],
                "n_threads": [],
                "skip": [],
            }
        elif current is None:
            continue
        elif re.match(r"^    skip:\s*$", line):
            in_skip = True
        elif re.match(r"^    \w+:", line):
            # Any other `    key:` at sweep-level closes the skip block.
            in_skip = False
            if re.match(r"^    lin_solvers:", line):
                current["lin_solvers"] = [int(x) for x in re.findall(r"\d+", line.split(":", 1)[1])]
            elif re.match(r"^    n_threads:", line):
                current["n_threads"] = [int(x) for x in re.findall(r"\d+", line.split(":", 1)[1])]
        elif in_skip and re.match(r"^      - \{\s*lin_solver:", line):
            m = re.search(r"lin_solver:\s*(\d+).*nx:\s*(\d+).*nz:\s*(\d+)", line)
            if m:
                current["skip"].append({
                    "lin_solver": int(m.group(1)),
                    "nx": int(m.group(2)),
                    "nz": int(m.group(3)),
                })
        elif re.match(r"^      - \{\s*nx:", line):
            m = re.search(r"nx:\s*(\d+).*nz:\s*(\d+)", line)
            if m:
                current["grids"].append((int(m.group(1)), int(m.group(2))))
    emit()
    return rows


if __name__ == "__main__":
    grids_file = os.environ.get("GRIDS_FILE")
    if not grids_file:
        print("GRIDS_FILE env var is required", file=sys.stderr)
        sys.exit(2)
    for row in parse(grids_file):
        print(row)
