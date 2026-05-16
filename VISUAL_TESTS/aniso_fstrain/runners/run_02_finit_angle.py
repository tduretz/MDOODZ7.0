#!/usr/bin/env python3
"""Test 2 — F-init geometry across aniso_angle.

Six cells with aniso_angle ∈ {0, 30, 60, 80, 90, 120}°, all else fixed:
aniso_factor=4 (→ γ_eff≈1.03, FS_AR≈2.68), ani_fstrain=3 inheritance ON,
Nt=2 (so we can read step-0 outputs), mechanical=0, BSR=1e-15. Plot
shows that F-major aligns with the foliation (aniso_angle−90°).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _lib import Cell, default_subs, run_sweep

ANGLES = [0, 30, 60, 80, 90, 120]


def build_cells():
    cells = []
    for ang in ANGLES:
        subs = default_subs(
            WRITER_STEP=1,
            SCALE_ETA="1.0", SCALE_L="1.0", SCALE_V="1.0", SCALE_T="1",
            NX=21, NZ=21, NT=2, DT="0.001",
            MECHANICAL=0,
            BKG_STRAIN_RATE="1e-15",
            BKG_TEMPERATURE="300.0",
            PWLV=0, ETA0="1e0",
            ANISO_ANGLE=ang,
            ANISO_FACTOR="4.0",
            ANI_FSTRAIN=3,
            ANI_RELAX_EPS_MAX="-1.0",
        )
        cells.append(Cell(name=f"a{ang:03d}", subs=subs,
                          meta={"aniso_angle": ang}))
    return cells


def main():
    run_sweep("02_finit_angle", build_cells(),
              max_workers=4, threads_per_cell=2, timeout=300)


if __name__ == "__main__":
    main()
