#!/usr/bin/env python3
"""Extract Nusselt number and Vrms from Blankenbach HDF5 output.

Mirrors the C++ BlankenBenchTests::NusseltAndVrms computation exactly.

Usage:
  python3 extract_blankenbach.py <hdf5_file>
  python3 extract_blankenbach.py BlankenBench/Output12100.gzip.h5
  python3 extract_blankenbach.py BlankenBench/Output*.gzip.h5   # time series
"""
import sys
import glob
import numpy as np
import h5py

# Blankenbach Case 1a physical parameters
T_top_K  = 273.15    # K
DeltaT   = 1000.0    # K
H        = 1.0e5     # m
rho      = 3300.0    # kg/m³
k_cond   = 3.3       # W/m/K
Cp       = 1250.0    # J/kg/K
kappa    = k_cond / (rho * Cp)  # thermal diffusivity

# Published reference (Blankenbach et al. 1989, Case 1a)
Nu_ref   = 4.884409
Vrms_ref = 42.864947


def extract(path):
    with h5py.File(path, "r") as f:
        params = f["Model/Params"][:]
        Nx = int(params[3])
        Nz = int(params[4])
        ncx, ncz = Nx - 1, Nz - 1

        T  = np.array(f["Centers/T"], dtype=np.float64)
        Vx = np.array(f["VxNodes/Vx"], dtype=np.float64)
        Vz = np.array(f["VzNodes/Vz"], dtype=np.float64)

        Nb_part = int(f["Model/Nb_part"][0]) if "Model/Nb_part" in f else -1

    T_bot_K = T_top_K + DeltaT  # 1273.15 K
    dz = H / ncz

    # --- Nusselt number (top-surface heat flux) ---
    T_top_row = T[(ncz - 1) * ncx : ncz * ncx]
    dTdz_top = (T_top_K - T_top_row) / (dz / 2.0)
    Nu_top = -(H / DeltaT) * np.mean(dTdz_top)

    # --- Nusselt number (bottom-surface heat flux) ---
    T_bot_row = T[0 : ncx]
    dTdz_bot = (T_bot_row - T_bot_K) / (dz / 2.0)
    Nu_bot = -(H / DeltaT) * np.mean(dTdz_bot)

    # --- Vrms (interpolate staggered → centres) ---
    sumV2 = 0.0
    for j in range(ncz):
        for i in range(ncx):
            vx_c = 0.5 * (Vx[j * Nx + i] + Vx[j * Nx + (i + 1)])
            vz_c = 0.5 * (Vz[j * (Nx + 1) + i] + Vz[(j + 1) * (Nx + 1) + i])
            sumV2 += vx_c * vx_c + vz_c * vz_c
    Vrms_SI = np.sqrt(sumV2 / (ncx * ncz))
    Vrms_nd = Vrms_SI * H / kappa

    return Nu_top, Nu_bot, Vrms_nd, Nb_part


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    # Expand globs
    paths = []
    for arg in sys.argv[1:]:
        expanded = sorted(glob.glob(arg))
        paths.extend(expanded if expanded else [arg])

    if len(paths) == 1:
        Nu_top, Nu_bot, Vrms, Nb = extract(paths[0])
        print(f"File: {paths[0]}")
        print(f"  Nu_top = {Nu_top:.6f}  (ref: {Nu_ref}  err: {abs(Nu_top - Nu_ref) / Nu_ref * 100:.2f}%)")
        print(f"  Nu_bot = {Nu_bot:.6f}  (ref: {Nu_ref}  err: {abs(Nu_bot - Nu_ref) / Nu_ref * 100:.2f}%)")
        print(f"  Vrms   = {Vrms:.6f}  (ref: {Vrms_ref}  err: {abs(Vrms - Vrms_ref) / Vrms_ref * 100:.2f}%)")
        if Nb >= 0:
            print(f"  Nb_part = {Nb}")
        if abs(Nu_top - Nu_bot) / max(Nu_ref, 1e-10) > 0.1:
            print(f"  WARNING: Nu_top != Nu_bot — thermal field NOT at steady state")
    else:
        # Time series
        print(f"{'step':>8s}  {'Nu_top':>10s}  {'Nu_bot':>10s}  {'Vrms':>10s}  {'Vrms_err%':>9s}  {'Nb_part':>8s}")
        print("-" * 68)
        for p in paths:
            import os, re
            m = re.search(r"Output(\d+)", os.path.basename(p))
            step = m.group(1) if m else "?"
            Nu_top, Nu_bot, Vrms, Nb = extract(p)
            vrms_err = abs(Vrms - Vrms_ref) / Vrms_ref * 100
            nb_str = str(Nb) if Nb >= 0 else "N/A"
            print(f"{step:>8s}  {Nu_top:10.6f}  {Nu_bot:10.6f}  {Vrms:10.6f}  {vrms_err:8.2f}%  {nb_str:>8s}")


if __name__ == "__main__":
    main()
