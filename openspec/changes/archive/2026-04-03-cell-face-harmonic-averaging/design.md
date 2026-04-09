## Context

MDOODZ's Stokes solver uses a staggered-grid FD scheme where the constitutive tensor components (`D11_n`, `D22_n` at cell centres, `D33_s` at vertices) are evaluated at their respective grid points and then read by the stencil assembly in `StokesAssemblyDecoupled.c`. For the x-momentum equation, the stencil reads `D11_n` at iPrE and iPrW (the pressure cells east and west of the Vx node), and `D33_s` at ixyN and ixyS (the vertex points north and south). For z-momentum, it reads `D22_n` at iPrS/iPrN and `D33_s` at ixyE/ixyW.

Currently there is no averaging between adjacent D values — each stencil coefficient uses the D value from its specific grid point. When a viscosity discontinuity falls between two grid points, this produces first-order truncation error at the interface cells.

The `eta_average` parameter (0/1/2) controls an unrelated operation: how marker viscosities are mixed to grid nodes *within* each cell during P2G interpolation.

## Goals / Non-Goals

**Goals:**
- Add a `cell_avg` model parameter (0=direct, 1=harmonic) that controls cell-face averaging of the constitutive tensor when assembling the Stokes stencil
- Default value `cell_avg=0` preserves all existing behaviour — zero regression risk
- Verify the effect on SolVi L2 error and convergence order at 3+ resolutions
- Ensure Newton convergence is not degraded by maintaining consistency between stencil and Jacobian

**Non-Goals:**
- Geometric averaging mode (only harmonic, which is the theoretically optimal choice for viscosity jumps)
- Modifying the anisotropic D-tensor components beyond D11/D22/D33 (D12, D13, D14 etc. are small correction terms; averaging them is a future extension)
- Changing how `eta_average` works (marker-to-node phase mixing is a separate concern)
- Achieving second-order convergence (that would require immersed-boundary methods)

## Decisions

### 1. Average at the D-tensor read, not upstream

**Decision**: Insert harmonic averaging in `StokesAssemblyDecoupled.c` at the point where `D11E`, `D11W`, `D33N`, `D33S` (and the z-momentum equivalents) are read from the mesh arrays — not in `StokesRoutines.c` where D is computed.

**Rationale**: 
- The stencil for Vx at point (i+½, j) needs the normal viscosity at the east face (between cell centres i and i+1) and at the west face (between i-1 and i). The D11 values at these centres are already computed. The cell-face value is `D11_face = 2·D11_i·D11_{i+1} / (D11_i + D11_{i+1})`. 
- Doing it at read-time means: (a) no new mesh arrays needed, (b) the change is localised to 2 functions, (c) `D11_n`/`D22_n`/`D33_s` remain point-values, available for other uses (stress evaluation, output).

**Alternative considered**: Pre-compute averaged arrays (e.g. `D11_n_face_x[]`). Rejected because it doubles memory for D arrays and requires extra loops.

### 2. Harmonic mean formula with safety for near-zero values

**Decision**: Use `D_face = 2·D_a·D_b / (D_a + D_b)` with a guard: if `D_a + D_b < 1e-30` then `D_face = 0.0`.

**Rationale**: Standard harmonic mean. The guard prevents division by zero when both sides are air cells (D=0). The `inE`/`inW` flags already zero out contributions from air cells, so the guard is a safety net.

### 3. Average only D11, D22, D33 (principal diagonal terms)

**Decision**: Apply harmonic averaging only to `D11_n`, `D22_n`, and `D33_s` reads. The cross-terms (`D12_n`, `D13_n`, `D14_n`, `D23_n`, `D24_n`, `D31_s`, `D32_s`, `D34_s`) are left as direct reads.

**Rationale**:
- In the isotropic case, D12=D13=D14=0, so averaging is irrelevant for cross-terms.
- In the anisotropic case, the cross-terms are proportional to `aniS·d1·eta` or `aniS·d2·eta` — they're corrections to the principal viscosity. Harmonic-averaging these could produce non-physical results (they can be negative). Starting with diagonal-only is safe and captures the dominant effect.
- Literature benchmarks (Deubelbeiss & Kaus 2008) tested harmonic averaging of the scalar viscosity η, which corresponds to our D11=D22=2η, D33=2η in the isotropic case.

### 4. Gate on `model.cell_avg` in the stencil functions

**Decision**: Wrap the averaging in a simple conditional:
```c
double D11E, D11W;
if (model.cell_avg == 1) {
    // Harmonic mean of D11 at east face (between iPrW and iPrE) 
    // and its further-east neighbour is not needed — iPrE and iPrW
    // ARE the two cells sharing the face where the Vx node sits.
    D11E = mesh->D11_n[iPrE]; // keep as-is: iPrE is already the cell-centre east of the face
    D11W = mesh->D11_n[iPrW]; // keep as-is
    // The harmonic average applies when these values enter ADJACENT 
    // flux differences — see detail below
} else {
    D11E = mesh->D11_n[iPrE];
    D11W = mesh->D11_n[iPrW];
}
```

**Correction — deeper analysis**: Actually, looking at the stencil more carefully, `D11E` and `D11W` are already the viscosity values at the two cell centres flanking the Vx node. The finite-difference for `∂(τ_xx)/∂x` at the Vx node uses `(D11E·ε̇_E - D11W·ε̇_W) / dx`. The "cell face" is where the Vx node sits, and D11E/D11W are the two sides. 

The actual averaging opportunity is: instead of using `D11E` directly in the `uE` coefficient and `D11W` directly in `uW`, we could replace both with a single harmonic-averaged face value. But this changes the stencil structure — the 5-point stencil coefficients are derived symbolically and use D11E and D11W independently.

**Revised approach**: The cleanest intervention is to average the D values **before** they enter the symbolic stencil expressions:
```c
if (model.cell_avg == 1 && (mesh->D11_n[iPrE] + mesh->D11_n[iPrW]) > 1e-30) {
    double harm = 2.0 * mesh->D11_n[iPrE] * mesh->D11_n[iPrW] / 
                       (mesh->D11_n[iPrE] + mesh->D11_n[iPrW]);
    D11E = harm;
    D11W = harm;
} else {
    D11E = mesh->D11_n[iPrE];
    D11W = mesh->D11_n[iPrW];
}
```
Setting both D11E=D11W=harmonic_mean effectively gives a constant viscosity across the face, which is exactly what harmonic face averaging does. Similarly for D33N/D33S in x-momentum, and D22S/D22N, D33E/D33W in z-momentum.

### 5. Jacobian consistency

**Decision**: The Jacobian in `FD_Jacobian.c` does NOT directly read D-tensor arrays — it recomputes viscosity via `PhaseRheologyLoop_v1` with perturbed strain rates. The D tensors are then recomputed in `StokesRoutines.c`. The stencil assembly function is called again with the perturbed D values.

**Key insight**: The Jacobian calls the same `Xmomentum_InnerNodesDecoupled` and `Zmomentum_InnerNodesDecoupled` functions. So if we modify those functions to do harmonic averaging when `model.cell_avg=1`, the Jacobian automatically gets the consistent linearisation. No separate Jacobian modification is needed.

### 6. Parameter naming and reading

**Decision**: Use `cell_avg` in the `.txt` file, read in `InputOutput.c` alongside `eta_average`:
```c
model.cell_avg = ReadInt2(fin, "cell_avg", 0);
```
Add `int cell_avg` field to the `params` struct in `mdoodz-private.h`.

**Alternative considered**: Extending `eta_average` to include cell-face modes (e.g., values 3, 4). Rejected because `eta_average` controls marker-to-node mixing (a different operation) and overloading it would be confusing.

## Risks / Trade-offs

**[Risk] Setting D11E=D11W=harmonic homogenises the viscosity across the control volume, potentially over-smoothing in regions with gradual viscosity gradients**
→ Mitigation: the harmonic mean approaches the direct value when D11E ≈ D11W. Over-smoothing only occurs at sharp jumps, which is exactly where it helps. Verify with SolVi that uniform-viscosity regions are unaffected.

**[Risk] Anisotropic stencil has cross-terms (D13, D31, etc.) that couple Vx to Vz through off-diagonal viscosity — these are not averaged**
→ Mitigation: cross-terms are zero in the isotropic case and small corrections in the anisotropic case. Start isotropic-only testing. Future work can extend to cross-terms if needed.

**[Risk] Newton convergence could degrade if the harmonic mean creates steeper viscosity gradients in the linearised system**
→ Mitigation: Monitor Newton iteration count in SolVi tests. The Jacobian inherits the same averaging since it calls the same stencil functions.

**[Risk] Boundary cells near air (phase=-1) may have D=0 on one side**
→ Mitigation: the `inE`/`inW`/`inN`/`inS` flags already zero out air-cell contributions in the stencil. The guard `D_a + D_b > 1e-30` prevents division by zero.
