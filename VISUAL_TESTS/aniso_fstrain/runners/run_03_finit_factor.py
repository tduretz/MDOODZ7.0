#!/usr/bin/env python3
"""Test 3 — δ round-trip across aniso_factor.

Nine cells with aniso_factor ∈ {1, 2, 3, 4, 6, 8, 10, 12, 14}, all else
fixed: aniso_angle=0, ani_fstrain=3, Nt=2, mechanical=0, BSR=1e-15.
The plot extracts Centers/aniso_delta at step 0 and confirms it lies on
the y=x identity (Hansen round-trip is exact).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _lib import Cell, default_subs, run_sweep

FACTORS = [1, 2, 3, 4, 6, 8, 10, 12, 14]


def build_cells():
    cells = []
    for f in FACTORS:
        subs = default_subs(
            WRITER_STEP=1,
            SCALE_ETA="1.0", SCALE_L="1.0", SCALE_V="1.0", SCALE_T="1",
            NX=21, NZ=21, NT=2, DT="0.001",
            MECHANICAL=0,
            BKG_STRAIN_RATE="1e-15",
            BKG_TEMPERATURE="300.0",
            PWLV=0, ETA0="1e0",
            ANISO_ANGLE=0,
            ANISO_FACTOR=f"{f:.1f}",
            ANI_FSTRAIN=3,
            ANI_RELAX_EPS_MAX="-1.0",
        )
        cells.append(Cell(name=f"f{f:02d}", subs=subs,
                          meta={"aniso_factor": f}))
    return cells


def main():
    run_sweep("03_finit_factor", build_cells(),
              max_workers=4, threads_per_cell=2, timeout=300)


if __name__ == "__main__":
    main()
