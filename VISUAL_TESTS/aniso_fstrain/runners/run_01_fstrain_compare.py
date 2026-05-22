#!/usr/bin/env python3
"""Test 1 — ani_fstrain 1 vs 2 vs 3 comparison.

Three cells (F1/F2/F3) sharing all kinematics; the only sweep dim is
`ani_fstrain ∈ {1, 2, 3}`. Olivine simple shear γ̇=2·BSR=2, Nt=17,
dt=0.5 → γ_max≈16, T=300K, aniso_factor=1 (no inheritance).
Cold T + default strain-rate gate ON → F3 should overlay F2 exactly.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
from _lib import Cell, default_subs, run_sweep

CASES = {
    "f1": 1,
    "f2": 2,
    "f3": 3,
}


def build_cells():
    cells = []
    for name, fmode in CASES.items():
        subs = default_subs(
            WRITER_STEP=1,
            SCALE_ETA="1e0", SCALE_L="1e0", SCALE_V="1e0", SCALE_T="1",
            NX=21, NZ=21, NT=17, DT="0.5",
            MECHANICAL=1,
            BKG_STRAIN_RATE="1.0",
            BKG_TEMPERATURE="300.0",
            PWLV=0, ETA0="1e0",
            ANISO_ANGLE=45,
            ANISO_FACTOR="1.0",
            ANI_FSTRAIN=fmode,
            # Test 1 sister .txt files DID NOT set ani_relax_eps_max →
            # MDOODZ default 1e-13 s^-1 (gate ON). We make that explicit.
            ANI_RELAX_EPS_MAX="1e-13",
        )
        cells.append(Cell(name=name, subs=subs, meta={"ani_fstrain": fmode}))
    return cells


def main():
    run_sweep("01_fstrain_compare", build_cells(),
              max_workers=3, threads_per_cell=4, timeout=600)


if __name__ == "__main__":
    main()
