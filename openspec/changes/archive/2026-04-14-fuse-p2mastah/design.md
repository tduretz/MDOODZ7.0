## Context

After Experiment 6 (`fix-accumulated-strain-omp`), the interp phase accounts for 29% of wall time (2.56s at 16 threads on c5ad.4xlarge, 1001×801) and is bandwidth-saturated — flat from 8 to 16 threads. The root cause is redundant particle-coordinate reads: 24 individual P2Mastah calls per timestep (hot path) each loop over all ~12M particles, re-reading `particles.x[k]` and `particles.z[k]` every time.

The existing P2Mastah function (mode 2, atomic scatter) already eliminates per-call allocation and thread-local buffers. The remaining optimization opportunity is fusing multiple calls that share a centroid type and stencil width into a single pass over the particle array, computing grid indices and bilinear weights once per particle and scattering N fields simultaneously.

**Call-site inventory (hot path — per timestep):**
- Lines 539–594 (pre-solve interpolation): 22 calls across cent/vert/vxnodes/vznodes
- Lines 1094–1095 (chemical diffusion): 2 calls (vxnodes/vznodes singletons)
- Total: 24 calls per timestep

**Init phase** (lines 214–350): ~29 calls, runs once — not targeted for fusion.

## Goals

- Reduce `interp_s` from ~2.56s to <1.5s at 16 threads (1001×801 on c5ad.4xlarge)
- Maintain full backward compatibility: `interp_mode = 0/1/2` unchanged; mode 3 is opt-in
- Fused output bit-identical to mode 2 (same atomic scatter math)
- CI: single BlankenBench with all optimizations (mode 3 + thermal_solver=1)

## Non-Goals

- Fusing init-phase P2Mastah calls (runs once, negligible impact)
- Mixed-stencil fusion (combining stencil=1 and stencil=9 fields in one batch)
- Changing CountPartCell / phase-percentage (prop=1) paths
- Fusing across different centroid types

## Decisions

### D1: Group by (centroid, stencil) — not just centroid

Fields can only share a WM (weight accumulator) array when they use the same centroid type AND stencil width. Stencil=1 uses `dx_itp = dx/2` (1-cell), stencil=9 uses `dx_itp = 3*dx/2` (9-cell), producing different weights. Sharing WM across stencil widths would produce incorrect normalization.

**Consequence**: The pre-solve hot path splits into batches by `(centroid, stencil)`:

| Batch | Centroid | Stencil | Fields | Count |
|-------|----------|---------|--------|-------|
| 1 | cent (1) | 1 | Cp, Qr, T0_n, divth0_n | 4 |
| 2 | cent (1) | interp | p0_n, phi0_n, [sxxd0, szzd0]_elastic, [d1_n, d2_n, angle_n]_aniso, X0_n, noise_n, d0_n | 6–10 |
| 3 | vert (0) | interp | p0_s, [sxz0]_elastic, [d1_s, d2_s, angle_s]_aniso, X0_s, noise_s | 3–7 |
| — | vxnodes | 1 | kx | 1 (singleton) |
| — | vznodes | 1 | kz | 1 (singleton) |

**Alternative rejected**: Mixed-stencil fusion (compute both weight sets per particle). Adds significant complexity with marginal benefit — only 4 stencil=1 fields would join the cent-interp mega-batch.

### D2: Dependency barriers require sub-batching cent-interp and vert-interp

Intervening grid operations between P2Mastah calls create dependency barriers:

```
Batch 1 (cent, s=1):  Cp, Qr, T0_n, divth0_n          [4 fields]
  ── ArrayEqualArray(mesh.T, mesh.T0_n) ──
Batch 2a (cent, interp): p0_n, phi0_n, [sxxd0, szzd0]  [2–4 fields]
Batch 2b (vert, interp): p0_s, [sxz0]                   [1–2 fields]
  ── InterpCentroidsToVerticesDouble × 2, InterpVerticesToCentroidsDouble, ShearModCompExpGrid ──
Batch 3a (cent, interp): [d1_n, d2_n, angle_n]          [0–3 fields]
Batch 3b (vert, interp): [d1_s, d2_s, angle_s]          [0–3 fields]
  ── FiniteStrainAspectRatio ──
Batch 4a (cent, interp): X0_n, noise_n, d0_n            [3 fields]
Batch 4b (vert, interp): X0_s, noise_s                  [2 fields]
  ── ArrayEqualArray(mesh.d_n, mesh.d0_n) ──
Singletons: kx (vxnodes, s=1), kz (vznodes, s=1)       [2 calls]
```

**Result**: 24 individual calls → 7 fused batches + 2 singletons = 9 particle loops instead of 24. A **2.7× reduction** in particle-coordinate memory traffic.

Calls within each sub-batch are reordered from their original interleaved cent/vert order to group by centroid. This is safe because all P2Mastah calls within the same barrier region write to independent destination arrays.

### D3: P2MastahField descriptor struct

```c
typedef struct {
    double *src;       // source: mat_prop array or particles.field (NULL for flag<0)
    double *dst;       // destination grid array (e.g. mesh.Cp)
    char   *BCtype;    // boundary condition type array (e.g. mesh.BCp.type)
    int     flag;      // 0=mat_prop, 1=particle field, -1/-2/-3=anisotropy
    int     avg;       // 0=arithmetic, 1=harmonic, 2=geometric
    int     prop;      // 0=field interpolation (always 0 in fused batches)
    int     stencil;   // 1 or 9 (must match all fields in batch)
} P2MastahField;
```

The `stencil` field is included per-descriptor for self-documentation but must be uniform within a batch. `prop` must be 0 (phase percentage calls excluded per spec). `avg` is 0 in all current call sites but included for forward compatibility.

### D4: P2Mastah_Fused function signature and architecture

```c
void P2Mastah_Fused(
    params       *model,
    markers       particles,
    grid         *mesh,
    int           centroid,
    int           interp_stencil,   // shared stencil for this batch
    P2MastahField *fields,
    int           nfields,
    InterpBufPool *pool
);
```

**Architecture:**

1. **Setup** (same as P2Mastah): Derive Nx, Nz, X_vect, Z_vect from `centroid`; compute `dx_itp`, `dz_itp` from `interp_stencil`.

2. **Buffer init**: Get `pool_idx` from centroid. Memset shared `WM = pool->WM[pool_idx]` and `nfields` BMWM arrays from `pool->BMWM_fused[0..nfields-1]`.

3. **Scatter loop** (`#pragma omp parallel for`): For each active particle k:
   - Compute `ip`, `jp` (grid column/row) — **once**
   - Compute `dxm`, `dzm` (distances to node) — **once**
   - For each field f in `0..nfields-1`:
     - Compute `mark_val` from `fields[f].flag` and field's `src` array
     - Accumulate `mark_val * weight` into `pool->BMWM_fused[f][kp]` via `omp atomic`
   - Accumulate `weight` into `WM[kp]` via `omp atomic` — **once per stencil node**
   - For stencil=9: repeat for 8 neighbors (N, S, E, W, NE, NW, SE, SW)

4. **Normalization** (`#pragma omp parallel for`): For each node i:
   - For each field f: `fields[f].dst[i] = pool->BMWM_fused[f][i] / WM[i]`
   - Apply harmonic/geometric inversion if `fields[f].avg != 0`

**Key optimization**: WM is accumulated once per stencil node per particle (not per-field), reducing atomic operations from `2N` to `N+1` per stencil node (N BMWM atomics + 1 WM atomic, vs 2N in separate calls).

### D5: Shared BMWM pool — flat array reused across centroids

```c
// Additions to InterpBufPool for mode 3:
int     max_fused_fields;      // 16 (sufficient for largest batch: 10 fields)
int     max_centroid_size;     // max(sizes[0..3])
double *BMWM_fused[16];       // [field_idx] → double[max_centroid_size]
```

Fused batches for different centroids are called sequentially, never concurrently, so the BMWM arrays are shared across all centroid types. Each array is sized to `max_centroid_size` — the max of all four centroid grid sizes.

**Memory cost** at 1001×801: `max_centroid_size` ≈ 803K doubles. 16 arrays × 6.4MB = **102MB** additional beyond existing mode 2 buffers (~50MB). Total mode 3 memory: ~152MB. Acceptable for the target c5ad.4xlarge (32GB RAM).

**Init** (`InterpBufPoolInit`, mode 3 branch):
```c
pool->max_fused_fields = 16;
pool->max_centroid_size = 0;
for (int c = 0; c < 4; c++) {
    if (pool->sizes[c] > pool->max_centroid_size)
        pool->max_centroid_size = pool->sizes[c];
}
for (int f = 0; f < pool->max_fused_fields; f++) {
    pool->BMWM_fused[f] = DoodzCalloc(pool->max_centroid_size, sizeof(double));
}
```

**Overflow**: If `nfields > pool->max_fused_fields`, log an error and fall back to individual mode 2 calls.

**Alternative rejected**: Pre-allocate per-centroid BMWM slots (4 centroids × 16 slots = 410MB). Wasteful since batches don't overlap across centroids.

### D6: Mode 3 dispatch pattern at call sites

At each fusible call group in `Main_DOODZ.c`:

```c
if (pool && pool->interp_mode == 3) {
    P2MastahField fields[16];
    int nf = 0;
    fields[nf++] = (P2MastahField){ .src = input.materials.Cp, .dst = mesh.Cp,
                                     .BCtype = mesh.BCp.type, .flag = 0, .avg = 0,
                                     .prop = 0, .stencil = 1 };
    fields[nf++] = (P2MastahField){ .src = input.materials.Qr, .dst = mesh.Qr,
                                     .BCtype = mesh.BCp.type, .flag = 0, .avg = 0,
                                     .prop = 0, .stencil = 1 };
    // ... more fields ...
    P2Mastah_Fused(&input.model, particles, &mesh, cent, 1, fields, nf, pool);
} else {
    P2Mastah(..., pool);  // existing calls unchanged
    P2Mastah(..., pool);
}
```

For conditional fields (elastic, anisotropy), build the array dynamically:

```c
if (pool && pool->interp_mode == 3) {
    P2MastahField fields[16];
    int nf = 0;
    fields[nf++] = (P2MastahField){ particles.P, mesh.p0_n, mesh.BCp.type, 1, 0, 0, input.model.interp_stencil };
    fields[nf++] = (P2MastahField){ particles.phi, mesh.phi0_n, mesh.BCp.type, 1, 0, 0, input.model.interp_stencil };
    if (input.model.elastic == 1) {
        fields[nf++] = (P2MastahField){ particles.sxxd, mesh.sxxd0, mesh.BCp.type, 1, 0, 0, input.model.interp_stencil };
        fields[nf++] = (P2MastahField){ particles.szzd, mesh.szzd0, mesh.BCp.type, 1, 0, 0, input.model.interp_stencil };
    }
    P2Mastah_Fused(&input.model, particles, &mesh, cent, input.model.interp_stencil, fields, nf, pool);
} else {
    P2Mastah(...p0_n..., pool);
    if (input.model.elastic == 1) { P2Mastah(...sxxd0..., pool); P2Mastah(...szzd0..., pool); }
    // ...
}
```

### D7: Singleton calls fall back to P2Mastah mode 2

Call sites with only 1 field for a centroid (kx on vxnodes, kz on vznodes, kc_x, kc_z) use the existing `P2Mastah` with `interp_mode == 2` semantics. When `interp_mode == 3`, the pool's `interp_mode` field is 3, but `P2Mastah` treats `imode >= 2` identically (atomic scatter path). No changes needed in `P2Mastah` for this.

### D8: CI — single BlankenBench with all optimizations

One new BlankenBench test fixture: `interp_mode = 3`, `thermal_solver = 1`. No separate mode-3-only variants. The existing mode 1/2 test fixtures in CI remain unchanged.

The BlankenBench `.txt` file gets the optimized parameters via `sed` in the CMake test definition, matching the existing pattern for mode 1/2 tests.

## Risks and Trade-offs

### R1: Increased atomic contention per stencil node
With N fields per particle, each stencil node gets `N+1` atomic operations (N BMWM + 1 WM) vs `2` in mode 2. For the largest batch (10 fields, stencil=9), this is `11 × 9 = 99` atomics per particle per batch vs `2 × 9 = 18` per individual call — but the total across all fields is `99` vs `10 × 18 = 180`. Net reduction: 45%. The risk is L1 cache line contention at hot grid nodes. Mitigation: benchmarking will reveal if contention negates the bandwidth savings.

### R2: Memory overhead (102MB for BMWM pool)
The 16 pre-allocated BMWM arrays add ~102MB at 1001×801. This is 0.3% of available RAM on the target machine (32GB). For smaller grids (e.g., 201×201), the overhead drops to ~5MB. Acceptable.

### R3: Code complexity at call sites
Each fusible group gets an `if (mode 3) { ... } else { ... }` block. The pre-solve section expands from ~22 lines to ~60–80 lines. Mitigation: the `else` branch preserves the exact original code, making diff review straightforward. A future refactoring could extract helper functions, but this is out of scope.

### R4: Call-site reordering correctness
Reordering P2Mastah calls within barrier regions is safe because all calls write to independent destination arrays and read from unchanged source arrays. Each barrier region has been manually verified (see D2). Risk: future code changes could add dependencies between reordered calls. Mitigation: comments at each barrier point documenting the dependency.

### R5: Stencil validation overhead
`P2Mastah_Fused` must validate that all descriptors share the same stencil value. This is a one-time loop over `nfields` (max 10) at the start of each call — negligible cost.
