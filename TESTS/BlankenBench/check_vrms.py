#!/usr/bin/env python3
"""Quick check: compare buggy vs correct Vrms/Nu from existing HDF5 output."""
import sys, h5py, numpy as np

h5path = sys.argv[1] if len(sys.argv) > 1 else "Output26000.gzip.h5"
f = h5py.File(h5path, 'r')

params = f['Model']['Params'][:]
Nx = int(params[3]); Nz = int(params[4])
ncx = Nx - 1; ncz = Nz - 1
print(f"Grid: Nx={Nx}, Nz={Nz}, ncx={ncx}, ncz={ncz}")

T  = f['Centers']['T'][:].reshape(ncz, ncx)
Vx = f['VxNodes']['Vx'][:].reshape(Nz+1, Nx)
Vz = f['VzNodes']['Vz'][:].reshape(Nz, Nx+1)
f.close()

H = 1e5; DeltaT = 1000.0; T_top_K = 273.15
rho = 3300.0; k = 3.3; Cp = 1250.0
kappa = k / (rho * Cp)
dz = H / ncz

# --- Nusselt number (first-order one-sided FD) ---
sumDtDz = 0.0
for i in range(ncx):
    T_centre = T[ncz-1, i]
    dTdz = (T_top_K - T_centre) / (dz / 2.0)
    sumDtDz += dTdz
Nu_1st = -(H / DeltaT) * (sumDtDz / ncx)

# --- Nusselt number (second-order one-sided FD) ---
sumDtDz2 = 0.0
for i in range(ncx):
    T_m1 = T[ncz-1, i]   # top cell centre
    T_m2 = T[ncz-2, i]   # second from top
    # f'(0) = [9*f(h) - f(3h) - 8*f(0)] / (3 * 2h)  where h=dz/2
    dTdz2 = (9*T_m1 - T_m2 - 8*T_top_K) / (3 * dz)
    sumDtDz2 += dTdz2
Nu_2nd = -(H / DeltaT) * (sumDtDz2 / ncx)

# --- Nusselt (bottom boundary, first-order) ---
sumDtDz_bot = 0.0
for i in range(ncx):
    T_centre = T[0, i]
    T_bot_K = T_top_K + DeltaT
    dTdz = (T_centre - T_bot_K) / (dz / 2.0)
    sumDtDz_bot += dTdz
Nu_bot_1st = (H / DeltaT) * (sumDtDz_bot / ncx)

# --- Nusselt (bottom boundary, second-order) ---
sumDtDz_bot2 = 0.0
T_bot_K = T_top_K + DeltaT
for i in range(ncx):
    T_m1 = T[0, i]   # bottom cell centre
    T_m2 = T[1, i]   # second from bottom
    dTdz2 = (9*T_m1 - T_m2 - 8*T_bot_K) / (3 * dz)
    sumDtDz_bot2 += dTdz2
Nu_bot_2nd = (H / DeltaT) * (sumDtDz_bot2 / ncx)

# --- Vrms BUGGY (test code: missing j+1 and i+1 offsets) ---
sumV2_bug = 0.0
for j in range(ncz):
    for i in range(ncx):
        vx_c = 0.5 * (Vx[j, i] + Vx[j, i+1])
        vz_c = 0.5 * (Vz[j, i] + Vz[j+1, i])
        sumV2_bug += vx_c**2 + vz_c**2
Vrms_bug = np.sqrt(sumV2_bug / (ncx * ncz)) * H / kappa

# --- Vrms CORRECT (extract_fields.cpp: j+1 for Vx rows, i+1 for Vz cols) ---
sumV2_fix = 0.0
for j in range(ncz):
    for i in range(ncx):
        vx_c = 0.5 * (Vx[j+1, i] + Vx[j+1, i+1])
        vz_c = 0.5 * (Vz[j, i+1] + Vz[j+1, i+1])
        sumV2_fix += vx_c**2 + vz_c**2
Vrms_fix = np.sqrt(sumV2_fix / (ncx * ncz)) * H / kappa

print(f"\n{'='*55}")
print(f"Blankenbach Case 1a — Reference: Nu=4.884409, Vrms=42.864947")
print(f"{'='*55}")
print(f"  Nu  (top, 1st-order):   {Nu_1st:.6f}  ({(Nu_1st/4.884409-1)*100:+.1f}%)")
print(f"  Nu  (top, 2nd-order):   {Nu_2nd:.6f}  ({(Nu_2nd/4.884409-1)*100:+.1f}%)")
print(f"  Nu  (bot, 1st-order):   {Nu_bot_1st:.6f}  ({(Nu_bot_1st/4.884409-1)*100:+.1f}%)")
print(f"  Nu  (bot, 2nd-order):   {Nu_bot_2nd:.6f}  ({(Nu_bot_2nd/4.884409-1)*100:+.1f}%)")
print(f"  Vrms (BUGGY indexing):  {Vrms_bug:.6f}  ({(Vrms_bug/42.864947-1)*100:+.1f}%)")
print(f"  Vrms (FIXED indexing):  {Vrms_fix:.6f}  ({(Vrms_fix/42.864947-1)*100:+.1f}%)")
print(f"{'='*55}")
