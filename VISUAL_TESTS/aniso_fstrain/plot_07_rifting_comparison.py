#!/usr/bin/env python3
"""
plot_07_rifting_comparison.py — regenerate the 3-panel ~5 Ma comparison
figure used in `CALIBRATION_AND_DYNAMIC_FSTRAIN.md` §3.

Reads one HDF5 snapshot from each of three RiftingRoman runs
(const-aniso, ani_fstrain=3 init1, ani_fstrain=3 v3 with F-init) and
plots the cell-centered δ field (`/Centers/ani_fac`) on a shared color
scale 1 → 14.132 (the Hansen 38-point-fit ceiling derived in §1.2 of
the summary), overlaid with thermal isotherms (200–1400 °C) and the
phase-2 (lithospheric-mantle) boundary.

Provenance / field definitions / lithospheric-mantle statistics:
see `notes/07_rifting_comparison_provenance.md`.

Typical use:

    python3 plot_07_rifting_comparison.py \
        --const  cmake-exec/RiftingRoman/run_constaniso/Output00650.gzip.h5 \
        --init1  cmake-exec/RiftingRoman/run_aniso3_init1/Output00700.gzip.h5 \
        --v3     cmake-exec/RiftingRoman/run_aniso3_v3/Output00650.gzip.h5 \
        --out    VISUAL_TESTS/img/aniso_fstrain/07_rifting_comparison_5Ma.png \
        --print-stats
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import NamedTuple

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

# Hansen 38-point combined-fit ceiling — see §1.2 of the summary markdown.
DELTA_MAX = 14.132

# Isotherm levels in °C; converted to K for /Centers/T which is in K.
ISOTHERMS_C = (200, 400, 600, 800, 1000, 1200, 1400)

# Phase index for lithospheric mantle in the RiftingRoman setup.
PHASE_LITH_MANTLE = 2

# Seconds per Myr (using 365.25 d).
SEC_PER_MYR = 365.25 * 86400.0 * 1.0e6


class Snapshot(NamedTuple):
    path: Path
    xc: np.ndarray         # (nx-1,) cell-center x [m]
    zc: np.ndarray         # (nz-1,) cell-center z [m]
    delta: np.ndarray      # (nz-1, nx-1) δ = ani_fac
    T_C: np.ndarray        # (nz-1, nx-1) temperature in °C
    compo: np.ndarray      # (nz-1, nx-1) phase from /VizGrid/compo
    time_s: float          # /Model/Params[0]


def load_snapshot(path: Path) -> Snapshot:
    with h5py.File(path, "r") as h:
        xc = h["/Model/xc_coord"][:].astype(np.float64)
        zc = h["/Model/zc_coord"][:].astype(np.float64)
        nx_c = xc.size
        nz_c = zc.size

        ani = h["/Centers/ani_fac"][:].astype(np.float64).reshape(nz_c, nx_c)
        T_K = h["/Centers/T"][:].astype(np.float64).reshape(nz_c, nx_c)
        # /VizGrid/compo is on the same (nz_c, nx_c) cell-center grid for
        # the standard RiftingRoman writer (verified against /Centers shape).
        compo = h["/VizGrid/compo"][:].reshape(nz_c, nx_c)

        time_s = float(h["/Model/Params"][0])

    return Snapshot(
        path=path,
        xc=xc,
        zc=zc,
        delta=ani,
        T_C=T_K - 273.15,
        compo=compo.astype(np.int32),
        time_s=time_s,
    )


def lith_mantle_stats(snap: Snapshot) -> dict:
    mask = snap.compo == PHASE_LITH_MANTLE
    if not mask.any():
        return dict(n=0, delta_min=float("nan"), delta_max=float("nan"),
                    median=float("nan"), frac_gt4_pct=float("nan"))
    sel = snap.delta[mask]
    return dict(
        n=int(mask.sum()),
        delta_min=float(sel.min()),
        delta_max=float(sel.max()),
        median=float(np.median(sel)),
        frac_gt4_pct=float((sel > 4.0).mean() * 100.0),
    )


def panel(ax, snap: Snapshot, title: str, *, norm: Normalize, cmap: str) -> None:
    x_km = snap.xc / 1.0e3
    z_km = snap.zc / 1.0e3

    pcm = ax.pcolormesh(
        x_km, z_km, snap.delta,
        norm=norm, cmap=cmap, shading="nearest",
    )

    # Isotherms.
    ax.contour(
        x_km, z_km, snap.T_C,
        levels=list(ISOTHERMS_C),
        colors="black", linewidths=0.5, alpha=0.7,
    )
    # Phase-2 boundary.
    ax.contour(
        x_km, z_km, (snap.compo == PHASE_LITH_MANTLE).astype(float),
        levels=[0.5], colors="red", linewidths=0.7,
    )

    ax.set_xlabel("x [km]")
    ax.set_title(f"{title}\nt = {snap.time_s / SEC_PER_MYR:.2f} Ma")
    ax.set_aspect("equal")
    return pcm


def make_figure(const: Snapshot, init1: Snapshot, v3: Snapshot, out: Path) -> None:
    fig, axes = plt.subplots(
        1, 3, figsize=(15, 4.2), sharey=True,
        gridspec_kw=dict(wspace=0.08),
    )
    norm = Normalize(vmin=1.0, vmax=DELTA_MAX)
    cmap = "viridis"

    pcm = panel(axes[0], const, "const-aniso (δ=4, ani_fstrain=0)", norm=norm, cmap=cmap)
    panel(axes[1], init1, "ani_fstrain=3 init1 (δ₀=1, F=I)", norm=norm, cmap=cmap)
    panel(axes[2], v3,    "ani_fstrain=3 v3 (δ₀=4 via F-init)", norm=norm, cmap=cmap)

    axes[0].set_ylabel("z [km]")

    cax = fig.add_axes([0.92, 0.18, 0.012, 0.68])
    cb = fig.colorbar(pcm, cax=cax)
    cb.set_label(r"$\delta = (\eta_n / \eta_s)^{1/n}$")
    cb.set_ticks([1, 2, 4, 6, 8, 10, 12, DELTA_MAX])

    fig.suptitle(
        r"Olivine viscous anisotropy $\delta$ at ~5 Ma — three rifting setups",
        y=1.02, fontsize=12,
    )

    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=160, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--const", required=True, type=Path,
                   help="HDF5 snapshot from run_constaniso (~5 Ma)")
    p.add_argument("--init1", required=True, type=Path,
                   help="HDF5 snapshot from run_aniso3_init1 (~5 Ma)")
    p.add_argument("--v3", required=True, type=Path,
                   help="HDF5 snapshot from run_aniso3_v3 (~5 Ma)")
    p.add_argument("--out", required=True, type=Path,
                   help="Output PNG path")
    p.add_argument("--print-stats", action="store_true",
                   help="Print per-panel lith-mantle δ statistics to stdout")
    args = p.parse_args()

    const = load_snapshot(args.const)
    init1 = load_snapshot(args.init1)
    v3 = load_snapshot(args.v3)

    if args.print_stats:
        for name, snap in [("const", const), ("init1", init1), ("v3", v3)]:
            s = lith_mantle_stats(snap)
            print(f"[{name}] t={snap.time_s/SEC_PER_MYR:.3f} Ma  "
                  f"δ_range=[{s['delta_min']:.2f}, {s['delta_max']:.2f}]  "
                  f"median={s['median']:.2f}  "
                  f"% > 4 = {s['frac_gt4_pct']:.1f}%  "
                  f"(N_cells={s['n']})")

    make_figure(const, init1, v3, args.out)
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
