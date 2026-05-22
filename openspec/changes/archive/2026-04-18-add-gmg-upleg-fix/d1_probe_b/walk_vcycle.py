#!/usr/bin/env python3
"""add-gmg-upleg-fix §4.2 — walk a V-cycle dump directory and print
per-event residual norms covering BOTH downleg and upleg.

The §3.1 `UplegAmplificationBound` fixture only instruments upleg stages
(prolongate, post_smooth, coarse_solve); §4 D1 probing on the 81² dump
showed per-stage upleg ratios all ≤ 1.4× yet whole-V-cycle ≈ 8.72×, so
the amplification must originate on the downleg. This script fills that
blind spot: it orders all snapshots from the first V-cycle by sequence
number, computes ‖r‖₂ = sqrt(Σres_u² + Σres_v² + Σres_p²), and prints a
linear trace with level, stage, phase, and the ratio to the previous
event's residual.

Usage:
    ./walk_vcycle.py <dump_dir>

Example:
    ./walk_vcycle.py cmake-build/TESTS/SolViUpleg81/vcycle_dump
"""
from __future__ import annotations

import os
import re
import sys

try:
    import h5py
except ImportError:
    sys.stderr.write("This script needs h5py. pip install h5py\n")
    sys.exit(2)

RX = re.compile(r"^level_(\d+)_(.+)_(pre|post)_(\d+)\.h5$")


def residual_norm(path: str) -> float:
    with h5py.File(path, "r") as f:
        ru = f["Fields/res_u"][()]
        rv = f["Fields/res_v"][()]
        rp = f["Fields/res_p"][()]
    return float(((ru * ru).sum() + (rv * rv).sum() + (rp * rp).sum()) ** 0.5)


def main() -> int:
    if len(sys.argv) != 2:
        sys.stderr.write(f"usage: {sys.argv[0]} <dump_dir>\n")
        return 2

    dump_dir = sys.argv[1]
    rows = []
    for name in os.listdir(dump_dir):
        m = RX.match(name)
        if not m:
            continue
        rows.append({
            "level": int(m.group(1)),
            "step":  m.group(2),
            "phase": m.group(3),
            "seq":   int(m.group(4)),
            "path":  os.path.join(dump_dir, name),
        })
    rows.sort(key=lambda r: r["seq"])

    if not rows:
        sys.stderr.write(f"no V-cycle snapshots found in {dump_dir}\n")
        return 1

    first_start = next(
        (r["seq"] for r in rows if r["level"] == 0 and r["step"] == "pre_smooth" and r["phase"] == "pre"),
        None,
    )
    first_end = next(
        (r["seq"] for r in rows
         if r["level"] == 0 and r["step"] == "post_smooth" and r["phase"] == "post"
         and first_start is not None and r["seq"] > first_start),
        None,
    )
    if first_start is None or first_end is None:
        sys.stderr.write("could not locate first V-cycle boundaries\n")
        return 1

    cycle = [r for r in rows if first_start <= r["seq"] <= first_end]

    print(f"dump_dir={dump_dir}")
    print(f"first V-cycle window: seq {first_start}..{first_end} ({len(cycle)} events)")
    print()
    print(f"{'seq':>4} {'lvl':>3} {'stage':<14} {'phase':<5} {'||r||':>12} {'ratio':>10} {'cum':>10}")
    print("-" * 72)

    prev = None
    entry = None
    for r in cycle:
        n = residual_norm(r["path"])
        ratio_prev = (n / prev) if prev and prev > 0 else float("nan")
        if entry is None:
            entry = n
        cum = (n / entry) if entry > 0 else float("nan")
        marker = ""
        if ratio_prev == ratio_prev and ratio_prev > 2.0:
            marker = "  <-- AMPLIFY > 2x"
        elif ratio_prev == ratio_prev and ratio_prev < 0.5:
            marker = "  <-- damp"
        print(f"{r['seq']:>4} {r['level']:>3} {r['step']:<14} {r['phase']:<5} "
              f"{n:>12.4e} {ratio_prev:>10.3f} {cum:>10.3f}{marker}")
        prev = n

    print()
    if entry and entry > 0:
        print(f"whole V-cycle ratio (exit / entry): {prev / entry:.4e} "
              f"(bound = 4.0; FAIL if > 4.0)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
