#!/usr/bin/env python3
"""
Visualise rotation advection test results.

Reads HDF5 output from all 3 reseed modes and plots a 4-panel figure:
  Panel 1: Analytical (initial composition)
  Panel 2: Mode 0 (legacy) final composition
  Panel 3: Mode 1 (current) final composition
  Panel 4: Mode 2 (v2) final composition

Usage:
    python3 plot_rotation.py [build_dir]

    build_dir  — path to the cmake build directory containing TESTS/
                 (default: cmake-build)
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import h5py
except ImportError:
    print("ERROR: h5py not installed. Run: pip install h5py")
    sys.exit(1)


def read_compo_hr(h5path):
    """Read VizGrid/compo_hr and Model/Params from an HDF5 output file."""
    with h5py.File(h5path, "r") as f:
        compo = np.array(f["VizGrid"]["compo_hr"], dtype=np.float64)
        params = np.array(f["Model"]["Params"])
        Nx = int(params[3])
        Nz = int(params[4])
        Lx = params[1]
        Lz = params[2]
    # compo_hr is on a 2×-resolution grid
    nxhr = 2 * (Nx - 1)
    nzhr = 2 * (Nz - 1)
    compo2d = compo.reshape((nzhr, nxhr))
    extent = [-Lx / 2, Lx / 2, -Lz / 2, Lz / 2]
    return compo2d, extent


def compute_l2(init, final):
    """RMS difference between two composition fields."""
    diff = final.astype(np.float64) - init.astype(np.float64)
    return np.sqrt(np.mean(diff ** 2))


def compute_mismatch_pct(init, final):
    """Percentage of cells that differ."""
    return 100.0 * np.sum(init != final) / init.size


def main():
    build_dir = sys.argv[1] if len(sys.argv) > 1 else "cmake-build"
    tests_dir = os.path.join(build_dir, "TESTS")

    # ---- Figure 1: Rigid rotation (analytical comparison) ----
    rotation_modes = [
        ("RotationMode0", "Mode 0 (legacy)"),
        ("RotationMode1", "Mode 1 (current)"),
        ("RotationMode2", "Mode 2 (v2)"),
    ]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)

    init_path = os.path.join(tests_dir, "RotationMode0", "Output00000.gzip.h5")
    if not os.path.exists(init_path):
        print(f"ERROR: {init_path} not found. Run the test first.")
        sys.exit(1)
    compo_init, extent = read_compo_hr(init_path)
    axes[0].imshow(compo_init, origin="lower", extent=extent, cmap="RdBu_r",
                   vmin=-0.5, vmax=1.5, aspect="equal")
    axes[0].set_title("Analytical\n(initial)", fontsize=11)
    axes[0].set_xlabel("x")
    axes[0].set_ylabel("z")

    for i, (folder, label) in enumerate(rotation_modes):
        final_path = os.path.join(tests_dir, folder, "Output01000.gzip.h5")
        if not os.path.exists(final_path):
            print(f"WARNING: {final_path} not found — skipping {label}")
            axes[i + 1].text(0.5, 0.5, "No data", transform=axes[i + 1].transAxes,
                             ha="center", va="center", fontsize=14)
            axes[i + 1].set_title(label, fontsize=11)
            continue

        compo_final, _ = read_compo_hr(final_path)
        l2 = compute_l2(compo_init, compo_final)
        mismatch = compute_mismatch_pct(compo_init, compo_final)

        axes[i + 1].imshow(compo_final, origin="lower", extent=extent, cmap="RdBu_r",
                           vmin=-0.5, vmax=1.5, aspect="equal")
        axes[i + 1].set_title(f"{label}\nL2={l2:.4f}  mis={mismatch:.1f}%", fontsize=11)
        axes[i + 1].set_xlabel("x")

    for ax in axes[1:]:
        ax.set_yticklabels([])

    fig.suptitle("Rigid-Body Rotation — 1 Full Revolution", fontsize=13, fontweight="bold")

    out_path = os.path.join(tests_dir, "rotation_advection_comparison.png")
    fig.savefig(out_path, dpi=150)
    print(f"Saved: {out_path}")
    plt.close(fig)

    # ---- Figure 2: Vortex (strong shear, inter-mode comparison) ----
    vortex_modes = [
        ("VortexMode0", "Mode 0 (legacy)"),
        ("VortexMode1", "Mode 1 (current)"),
        ("VortexMode2", "Mode 2 (v2)"),
    ]

    vortex_init_path = os.path.join(tests_dir, "VortexMode0", "Output00000.gzip.h5")
    if not os.path.exists(vortex_init_path):
        print("Vortex outputs not found — skipping vortex figure.")
        return

    compo_v_init, extent_v = read_compo_hr(vortex_init_path)

    fig2, axes2 = plt.subplots(1, 4, figsize=(16, 4), constrained_layout=True)

    axes2[0].imshow(compo_v_init, origin="lower", extent=extent_v, cmap="RdBu_r",
                    vmin=-0.5, vmax=1.5, aspect="equal")
    axes2[0].set_title("Initial\n(t = 0)", fontsize=11)
    axes2[0].set_xlabel("x")
    axes2[0].set_ylabel("z")

    vortex_finals = {}
    for i, (folder, label) in enumerate(vortex_modes):
        final_path = os.path.join(tests_dir, folder, "Output00500.gzip.h5")
        if not os.path.exists(final_path):
            print(f"WARNING: {final_path} not found — skipping {label}")
            axes2[i + 1].text(0.5, 0.5, "No data", transform=axes2[i + 1].transAxes,
                              ha="center", va="center", fontsize=14)
            axes2[i + 1].set_title(label, fontsize=11)
            continue

        compo_final, _ = read_compo_hr(final_path)
        vortex_finals[i] = compo_final

        axes2[i + 1].imshow(compo_final, origin="lower", extent=extent_v, cmap="RdBu_r",
                            vmin=-0.5, vmax=1.5, aspect="equal")
        axes2[i + 1].set_title(label, fontsize=11)
        axes2[i + 1].set_xlabel("x")

    # Compute inter-mode L2 (the actual quality metric) and print summary
    mode_names = ["Mode 0", "Mode 1", "Mode 2"]
    pairs = [(0, 1), (0, 2), (1, 2)]
    print("\n  Vortex inter-mode L2:")
    for a, b in pairs:
        if a in vortex_finals and b in vortex_finals:
            l2 = compute_l2(vortex_finals[a], vortex_finals[b])
            print(f"    {mode_names[a]} vs {mode_names[b]}: L2 = {l2:.4f}")

    for ax in axes2[1:]:
        ax.set_yticklabels([])

    fig2.suptitle("Vortex Advection — Strong Shear (Nx_part=4, min_part_cell=4)",
                  fontsize=13, fontweight="bold")

    out_path2 = os.path.join(tests_dir, "vortex_advection_comparison.png")
    fig2.savefig(out_path2, dpi=150)
    print(f"Saved: {out_path2}")
    plt.close(fig2)


if __name__ == "__main__":
    main()
