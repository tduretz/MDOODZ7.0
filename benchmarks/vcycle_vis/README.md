# V-cycle animation generator

Produces the V-cycle residual animation embedded in the
`add-gmg-stokes-defence` defence document (§10.10).

## Pipeline

```
  TESTS/VcycleAnimGenerate.cpp          (produces dumps)
             │  ctest -L experimental -R VcycleAnimGenerate
             ▼
  build/TESTS/SolViRes41_gmg_vcycle_anim/vcycle_dump/*.h5
             │  python3 benchmarks/vcycle_vis/make_gif.py <dump_dir>
             ▼
  openspec/changes/add-gmg-stokes-defence/figs/
     ├── vcycle_sol_cx_41.gif          (animated, 2 fps)
     └── vcycle_still_{01..06}.png     (printable stills)
```

## Regenerating the figures

```bash
# from repository root
cmake -S . -B build -DTEST=ON
cmake --build build --target VcycleAnimGenerate

cd build/TESTS
./VcycleAnimGenerate   # runs the 41x41 SolCx analog with gmg_dump_vcycle = 1
cd ../..

python3 benchmarks/vcycle_vis/make_gif.py \
    build/TESTS/SolViRes41_gmg_vcycle_anim/vcycle_dump \
    --output-dir openspec/changes/add-gmg-stokes-defence/figs
```

## Dependencies

`requirements.txt` lists the Python packages. All are stock scientific
Python:

* `h5py`    — load MDOODZ HDF5 dumps
* `numpy`   — FFT + radial binning
* `matplotlib` — 2-panel figure rendering
* `Pillow`  — animated-GIF assembly

## Dump schema

Each HDF5 file is produced by `gmg_vcycle_dump_snapshot` in
`MDLIB/MultigridStokes.c`:

| group     | dataset  | type / shape               | meaning                      |
|-----------|----------|----------------------------|------------------------------|
| `/Level`  | `dims`   | int[4] (Nx, Nz, Ncx, Ncz)  | face + centre counts         |
| `/Level`  | `spacing`| dbl[2] (dx, dz)            | cell size at this MG level   |
| `/Level`  | `level`  | int                        | 0 = fine, increases coarsely |
| `/Fields` | `Vx`     | dbl[Nx·Ncz]                | x-velocity iterate           |
| `/Fields` | `Vz`     | dbl[Ncx·Nz]                | z-velocity iterate           |
| `/Fields` | `P`      | dbl[Ncx·Ncz]               | pressure iterate             |
| `/Fields` | `res_u`  | dbl[Nx·Ncz]                | momentum-x residual          |
| `/Fields` | `res_v`  | dbl[Ncx·Nz]                | momentum-z residual          |
| `/Fields` | `res_p`  | dbl[Ncx·Ncz]               | continuity residual          |
| `/Meta`   | `step`   | c-str                      | operator name                |
| `/Meta`   | `phase`  | c-str                      | `pre` (before) / `post`      |
| `/Meta`   | `seq`    | int                        | global traversal order       |

Operator names correspond to the V-cycle in `Vcycle()`:
`pre_smooth`, `restrict`, `coarse_solve`, `prolongate`, `post_smooth`.

## Design notes

The animation is intentionally built from a single V-cycle of the very
first FGMRES iteration of the very first Picard iteration of step 1.
The `gmg_dump_vcycle` instrumentation is one-shot (`fired` latches at
finish), so running the fixture within a single process guarantees
exactly one V-cycle's worth of snapshots in the output tree.

The radial-power-spectrum panel is the defence's §10.3 "why multigrid"
money plot: the smoother attacks high wavenumbers, the coarse-grid
correction kills the low wavenumbers that the smoother can't touch. In
the animation this shows up as the spectrum "collapsing" from a broad
log-log line toward a narrow low-k tail as the V-cycle progresses.
