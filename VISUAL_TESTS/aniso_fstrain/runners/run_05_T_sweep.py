#!/usr/bin/env python3
"""Test 5 — T-sweep pure δ-decay (Arrhenius staircase).

Six cells at T ∈ {500, 700, 900, 1100, 1300, 1500} K, mechanical=1 with
tiny BSR=1e-20 (F effectively frozen), ani_relax_eps_max=-1, δ_init=4.
Total time ~1.27 Myr (Nt=400, dt=1e11 s).

Uses SI scales (eta=1e22, L=1e4, V=1e-9, T=1000 — so all .txt values are
SI; MDOODZ divides through internally).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _lib import Cell, default_subs, run_sweep

TEMPS = [500, 700, 900, 1100, 1300, 1500]


def build_cells():
    cells = []
    for T in TEMPS:
        subs = default_subs(
            WRITER_STEP=1,
            SCALE_ETA="1e+22", SCALE_L="10000.0", SCALE_V="1e-09",
            SCALE_T="1000.0",
            NX=11, NZ=11, NT=400, DT="1e11",
            MECHANICAL=1,
            BKG_STRAIN_RATE="1e-20",
            BKG_TEMPERATURE=f"{T:.1f}",
            PWLV=0, ETA0="1e0",
            ANISO_ANGLE=0,
            ANISO_FACTOR="4.0",
            ANI_FSTRAIN=3,
            ANI_RELAX_EPS_MAX="-1.0",
        )
        cells.append(Cell(name=f"T{T:04d}", subs=subs, meta={"T_K": T}))
    return cells


def main():
    run_sweep("05_T_sweep", build_cells(),
              max_workers=3, threads_per_cell=4, timeout=900)


if __name__ == "__main__":
    main()
