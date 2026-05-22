#!/usr/bin/env python3
"""Test 4 — inheritance reset-strain sweep.

Four cells with aniso_factor ∈ {1, 2, 4, 8} initial fabric, γ̇=2 simple
shear to γ=20 (Nt=21, dt=0.5), aniso_angle=80°, ani_fstrain=3,
ani_relax_eps_max=-1 (gate OFF to isolate F-init effect).
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _lib import Cell, default_subs, run_sweep

FACTORS = [1, 2, 4, 8]


def build_cells():
    cells = []
    for f in FACTORS:
        subs = default_subs(
            WRITER_STEP=1,
            SCALE_ETA="1.0", SCALE_L="1.0", SCALE_V="1.0", SCALE_T="1",
            NX=21, NZ=21, NT=21, DT="0.5",
            MECHANICAL=1,
            BKG_STRAIN_RATE="1.0",
            BKG_TEMPERATURE="300.0",
            PWLV=0, ETA0="1e0",
            ANISO_ANGLE=80,
            ANISO_FACTOR=f"{f:.1f}",
            ANI_FSTRAIN=3,
            ANI_RELAX_EPS_MAX="-1.0",
        )
        cells.append(Cell(name=f"f{f}", subs=subs,
                          meta={"aniso_factor": f}))
    return cells


def main():
    run_sweep("04_inheritance", build_cells(),
              max_workers=4, threads_per_cell=2, timeout=600)


if __name__ == "__main__":
    main()
