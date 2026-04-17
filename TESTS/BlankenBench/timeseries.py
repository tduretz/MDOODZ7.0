#!/usr/bin/env python3
"""Time series of Nu and Vrms across all Blankenbach output files."""
import h5py, numpy as np, glob

files = sorted(glob.glob('BlankenBench/Output*.gzip.h5'))
print(f"Found {len(files)} output files")
print(f"First: {files[0]}, Last: {files[-1]}")
print()

H = 1e5; DeltaT = 1000.0; T_top_K = 273.15; T_bot_K = T_top_K + DeltaT
rho = 3300.0; k = 3.3; Cp = 1250.0; kappa = k / (rho * Cp)

# Sample every ~10 files for time series
indices = list(range(0, len(files), max(1, len(files)//20))) + [len(files)-1]
indices = sorted(set(indices))

hdr = "{:>8s} {:>10s} {:>10s} {:>10s} {:>10s}".format(
    "Step", "Nu_top", "Nu_bot", "Vrms_bug", "Vrms_fix")
print(hdr)
print("-" * len(hdr))

for idx in indices:
    fp = files[idx]
    step = int(fp.split('Output')[1].split('.')[0])
    f = h5py.File(fp, 'r')
    params = f['Model']['Params'][:]
    Nx = int(params[3]); Nz = int(params[4])
    ncx = Nx - 1; ncz = Nz - 1
    dz = H / ncz

    T = f['Centers']['T'][:].reshape(ncz, ncx)
    Vx = f['VxNodes']['Vx'][:].reshape(Nz + 1, Nx)
    Vz = f['VzNodes']['Vz'][:].reshape(Nz, Nx + 1)
    f.close()

    # Nu top (1st order)
    sumDtDz = sum((T_top_K - T[ncz - 1, i]) / (dz / 2.0) for i in range(ncx))
    Nu_top = -(H / DeltaT) * (sumDtDz / ncx)

    # Nu bottom (1st order)
    sumDtDz_b = sum((T[0, i] - T_bot_K) / (dz / 2.0) for i in range(ncx))
    Nu_bot = (H / DeltaT) * (sumDtDz_b / ncx)

    # Vrms buggy & fixed
    sv2b = 0.0; sv2f = 0.0
    for j in range(ncz):
        for i in range(ncx):
            vxb = 0.5 * (Vx[j, i] + Vx[j, i + 1])
            vzb = 0.5 * (Vz[j, i] + Vz[j + 1, i])
            sv2b += vxb**2 + vzb**2
            vxf = 0.5 * (Vx[j + 1, i] + Vx[j + 1, i + 1])
            vzf = 0.5 * (Vz[j, i + 1] + Vz[j + 1, i + 1])
            sv2f += vxf**2 + vzf**2
    Vrms_b = np.sqrt(sv2b / (ncx * ncz)) * H / kappa
    Vrms_f = np.sqrt(sv2f / (ncx * ncz)) * H / kappa

    print("{:>8d} {:>10.4f} {:>10.4f} {:>10.4f} {:>10.4f}".format(
        step, Nu_top, Nu_bot, Vrms_b, Vrms_f))
