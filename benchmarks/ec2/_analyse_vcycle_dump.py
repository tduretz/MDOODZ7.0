#!/usr/bin/env python3
"""Analyse V-cycle HDF5 dump tree from gmg_dump_vcycle=1.

Produces two diagnostics for the `add-gmg-stokes-defence` stall probe:

1. Per-V-cycle level-0 residual norm history (ρ per V-cycle).
2. Per-level residual reduction inside the FIRST V-cycle only: for each
   level, ||res|| BEFORE restrict vs. AFTER prolongate. A level whose
   coarse correction is "blind" leaves ||res|| roughly unchanged.

Usage: python3 _analyse_vcycle_dump.py <vcycle_dump_dir> <out_prefix>
"""
from __future__ import annotations
import sys
import os
import re
from pathlib import Path

import h5py
import numpy as np


_FNAME_RE = re.compile(r"level_(\d+)_(\w+?)_(pre|post)_(\d+)\.h5$")


def _load_res_norm(path: Path) -> float:
    with h5py.File(path, "r") as f:
        ru = np.asarray(f["/Fields/res_u"])
        rv = np.asarray(f["/Fields/res_v"])
        rp = np.asarray(f["/Fields/res_p"])
    return float(np.sqrt(np.sum(ru * ru) + np.sum(rv * rv) + np.sum(rp * rp)))


def _scan(dump_dir: Path) -> list[dict]:
    rows: list[dict] = []
    for p in sorted(dump_dir.glob("*.h5")):
        m = _FNAME_RE.match(p.name)
        if not m:
            continue
        rows.append(
            {
                "level": int(m.group(1)),
                "step": m.group(2),
                "phase": m.group(3),
                "seq": int(m.group(4)),
                "path": p,
            }
        )
    rows.sort(key=lambda r: r["seq"])
    return rows


def _vcycle_boundaries(rows: list[dict]) -> list[tuple[int, int]]:
    """Return (start_seq, end_seq) for each V-cycle.

    Each V-cycle begins with `level_0 pre_smooth pre` and ends with
    `level_0 post_smooth post`. Seqs are global over the whole run.
    """
    starts = [r["seq"] for r in rows
              if r["level"] == 0 and r["step"] == "pre_smooth" and r["phase"] == "pre"]
    ends = [r["seq"] for r in rows
            if r["level"] == 0 and r["step"] == "post_smooth" and r["phase"] == "post"]
    return list(zip(starts, ends))


def vcycle_history(rows: list[dict], out_csv: Path) -> None:
    by_seq = {r["seq"]: r for r in rows}
    vcycles = _vcycle_boundaries(rows)
    lines = ["vcycle,r_start,r_end,rho"]
    for i, (s, e) in enumerate(vcycles):
        r_start = _load_res_norm(by_seq[s]["path"])
        r_end = _load_res_norm(by_seq[e]["path"])
        rho = (r_end / r_start) if r_start > 0 else float("nan")
        lines.append(f"{i},{r_start:.6e},{r_end:.6e},{rho:.6e}")
    out_csv.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_csv}  ({len(vcycles)} V-cycles)")


def first_vcycle_per_level(rows: list[dict], out_csv: Path) -> None:
    """Inside the first V-cycle, report per-level:
    ||res_before_restrict|| and ||res_after_prolongate||.

    On the downleg, for level k>0, `restrict pre` at level k-1 carries
    the residual that will be restricted; after coarse work and
    `prolongate post` at level k-1, the residual is what remains.
    So the "damage ratio" for the coarse correction through level k
    (taking k's contribution plus everything below) is
      ||res(level k-1, prolongate, post)|| / ||res(level k-1, restrict, pre)||.
    For k == n_levels - 1 the coarse_solve replaces restrict/prolongate.
    """
    by_seq = {r["seq"]: r for r in rows}
    vcycles = _vcycle_boundaries(rows)
    if not vcycles:
        out_csv.write_text("no vcycles detected\n")
        return
    seq_start, seq_end = vcycles[0]
    first = [r for r in rows if seq_start <= r["seq"] <= seq_end]

    n_levels = 1 + max(r["level"] for r in first)

    lines = ["level,event,res_before,res_after,ratio"]

    # For levels 0 .. n_levels-2: downleg restrict (pre/post) plus upleg
    # prolongate (pre/post). The pair we want is
    #   before = restrict.pre at level k (that's the residual at k just
    #            before we descend from k to k+1)
    #   after  = prolongate.post at level k (residual at k after the
    #            coarse correction from k+1 has been prolonged back and
    #            the post-smoother applied).
    for k in range(n_levels - 1):
        b = next((r for r in first
                  if r["level"] == k and r["step"] == "restrict" and r["phase"] == "pre"), None)
        a = next((r for r in first
                  if r["level"] == k and r["step"] == "prolongate" and r["phase"] == "post"), None)
        if b is None or a is None:
            lines.append(f"{k},MISSING,,,")
            continue
        rb = _load_res_norm(b["path"])
        ra = _load_res_norm(a["path"])
        ratio = (ra / rb) if rb > 0 else float("nan")
        lines.append(f"{k},coarse_correction_{k}_to_{k + 1},{rb:.6e},{ra:.6e},{ratio:.6e}")

    # Deepest level: coarse_solve.pre vs coarse_solve.post
    k = n_levels - 1
    b = next((r for r in first
              if r["level"] == k and r["step"] == "coarse_solve" and r["phase"] == "pre"), None)
    a = next((r for r in first
              if r["level"] == k and r["step"] == "coarse_solve" and r["phase"] == "post"), None)
    if b is not None and a is not None:
        rb = _load_res_norm(b["path"])
        ra = _load_res_norm(a["path"])
        ratio = (ra / rb) if rb > 0 else float("nan")
        lines.append(f"{k},coarse_solve,{rb:.6e},{ra:.6e},{ratio:.6e}")

    # Also: at each level on the downleg, pre_smooth before/after to see
    # whether the smoother is actually doing work per-level.
    for k in range(n_levels):
        b = next((r for r in first
                  if r["level"] == k and r["step"] == "pre_smooth" and r["phase"] == "pre"), None)
        a = next((r for r in first
                  if r["level"] == k and r["step"] == "pre_smooth" and r["phase"] == "post"), None)
        if b is None or a is None:
            continue
        rb = _load_res_norm(b["path"])
        ra = _load_res_norm(a["path"])
        ratio = (ra / rb) if rb > 0 else float("nan")
        lines.append(f"{k},pre_smooth,{rb:.6e},{ra:.6e},{ratio:.6e}")

    out_csv.write_text("\n".join(lines) + "\n")
    print(f"wrote {out_csv}  (n_levels={n_levels})")


def main() -> int:
    if len(sys.argv) < 3:
        print("usage: _analyse_vcycle_dump.py <vcycle_dump_dir> <out_prefix>",
              file=sys.stderr)
        return 2
    dump_dir = Path(sys.argv[1])
    out_prefix = Path(sys.argv[2])
    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    rows = _scan(dump_dir)
    print(f"scanned {len(rows)} snapshots from {dump_dir}")
    vcycle_history(rows, out_prefix.with_suffix(".vcycle_history.csv"))
    first_vcycle_per_level(rows, out_prefix.with_suffix(".first_vcycle_per_level.csv"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
