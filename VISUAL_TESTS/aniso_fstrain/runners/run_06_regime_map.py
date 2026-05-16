#!/usr/bin/env python3
"""Test 6 — 2D T × ε̇ regime map (42 cells).

For each (T, ε̇) cell, generate a per-cell .txt by patching the template
and run the unified AnisoFstrainBox binary. Collect the final δ at γ=5
for the heatmap.

Sweep grid: T ∈ {500, 700, 900, 1100, 1300, 1500} K
            ε̇ ∈ {1e-16, …, 1e-10} s^-1
NT=200 steps; dt is computed per-cell so total strain γ̇·Nt·dt = γ_total.
"""
from __future__ import annotations

import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from _lib import (BIN, RUNS_DIR, Cell, default_subs, ensure_build, render_template, run_cell)

TS = [500, 700, 900, 1100, 1300, 1500]                       # K
EDOTS = [1e-16, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10]    # s^-1
NT = 200
GAMMA_TOTAL = 5.0
TEST_NAME = "06_regime_map"


def build_subs(T_K: float, edot_SI: float) -> dict:
    # γ_total = γ̇·Nt·dt = 2·ε̇·Nt·dt → dt = γ_total / (2·ε̇·Nt) in SI seconds
    # (the .txt parser divides dt by scaling.t, with scaling.t = L/V = 1e13).
    dt_SI = GAMMA_TOTAL / (2.0 * edot_SI * NT)
    return default_subs(
        WRITER_STEP=NT,  # only output first & last; final δ at last step
        SCALE_ETA="1e+22", SCALE_L="10000.0", SCALE_V="1e-09",
        SCALE_T="1000.0",
        NX=11, NZ=11, NT=NT, DT=f"{dt_SI:.6e}",
        MECHANICAL=1,
        BKG_STRAIN_RATE=f"{edot_SI:.6e}",
        BKG_TEMPERATURE=f"{T_K:.1f}",
        PWLV=40, ETA0="1e22",
        ANISO_ANGLE=0,
        ANISO_FACTOR="4.0",
        ANI_FSTRAIN=3,
        ANI_RELAX_EPS_MAX="-1.0",
    )


def run_one(T_K: float, edot_SI: float, test_dir: Path) -> dict:
    name = f"cell_T{T_K:04.0f}_E{abs(int(round(np.log10(edot_SI)))):02d}"
    subs = build_subs(T_K, edot_SI)
    res = run_cell(test_dir / name, subs, threads=2, timeout=300)
    # Read final δ from last HDF5 file
    final_delta = np.nan
    h5dir = res["output_dir"]
    h5s = sorted(h5dir.glob("Output*.gzip.h5")) if h5dir.exists() else []
    if h5s:
        try:
            with h5py.File(h5s[-1], "r") as f:
                final_delta = float(f["Centers/aniso_delta"][:].mean())
        except Exception:
            pass
    return {
        "name":   name,
        "T":      T_K,
        "edot":   edot_SI,
        "rc":     res["rc"],
        "wall":   res["wall"],
        "delta":  final_delta,
        "nsteps": len(h5s),
    }


def main():
    ensure_build()
    test_dir = RUNS_DIR / TEST_NAME
    test_dir.mkdir(parents=True, exist_ok=True)
    tasks = [(T, e) for T in TS for e in EDOTS]
    print(f"[{TEST_NAME}] running {len(tasks)} cells (4 workers, OMP=2 each) ...",
          flush=True)
    results: list[dict] = []
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=4) as ex:
        futs = {ex.submit(run_one, T, e, test_dir): (T, e) for T, e in tasks}
        done = 0
        for fut in as_completed(futs):
            T, e = futs[fut]
            try:
                r = fut.result()
            except Exception as ex2:
                r = {"name": f"cell_T{T:04.0f}_E??", "T": T, "edot": e,
                     "rc": -2, "wall": 0.0, "delta": float("nan"), "nsteps": 0,
                     "exception": repr(ex2)}
            results.append(r)
            done += 1
            print(f"  [{done}/{len(tasks)}] T={r['T']} K, ε̇={r['edot']:.0e}: "
                  f"rc={r['rc']} wall={r['wall']:.1f}s δ_final={r['delta']:.3f} "
                  f"steps={r['nsteps']}", flush=True)
    print(f"[{TEST_NAME}] total wall: {time.time()-t0:.1f}s", flush=True)
    out_json = test_dir / "results.json"
    out_json.write_text(json.dumps(results, indent=2))
    print(f"results saved to {out_json}")


if __name__ == "__main__":
    main()
