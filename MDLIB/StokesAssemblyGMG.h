// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2023-2026 MDOODZ Developers team
//
// All rights reserved. This file is part of MDOODZ.
// =========================================================================
//
// StokesAssemblyGMG — pure read-only matrix-vector product accessors that
// mirror the coefficient formulas used by the sparse-matrix assembly in
// StokesAssemblyDecoupled.c. Option A per design D11: these accessors exist
// so the GMG V-cycle's fine-level operator can reproduce the exact same
// linear operator that CHOLMOD solves, without touching the assembly code
// that is shared by the other lin_solver paths.
//
// Correctness contract: for every BC tag, every rheology branch (Picard in
// Phase 1, Newton-Jacobian in Phase 2), and every anisotropy / compressibility
// setting, the matvec output here SHALL agree with A·x built via
// BuildStokesOperatorDecoupled to within 1e-12 relative tolerance. This is
// enforced by TESTS/StokesMatvecEquivalence.cpp (task 15.9–15.13).
//
// Non-goals: no assembly, no smoother, no transfer. Just A·x evaluation that
// reads the mesh's D-tensor entries, BC tags, and the input (u, v, p).

#pragma once

#include "mdoodz-private.h"

#ifdef __cplusplus
extern "C" {
#endif

// Evaluate the MDOODZ Stokes operator on full-grid velocity and pressure
// fields. The output arrays (ru, rv, rp) MUST be sized identically to the
// mesh's velocity/pressure arrays:
//   ru: Nx * (Nz+1)                (cell-face Vx layout, same as u)
//   rv: (Nx+1) * Nz                (cell-face Vz layout, same as v)
//   rp: (Nx-1) * (Nz-1)            (cell-centred P  layout, same as p)
//
// At cells that do NOT carry an equation (BC types 0, 30, 31, 11, 13, -12
// for velocity; 0, 30, 31 for pressure) the output is set to 0.0 — those
// rows have no corresponding entry in the assembled sparse matrix and so
// do not contribute to A·x.
//
// At each equation cell the output equals celvol · Σ coeff_j · x_j, where
// coeff_j are the per-neighbour stencil coefficients computed identically
// to Xmomentum_InnerNodesDecoupled / Zmomentum_InnerNodesDecoupled /
// Continuity_InnerNodesDecoupled (Picard variant for Phase 1; the Newton
// Jacobian variant is added in Phase 2). Neighbours whose BC type is 0,
// 31, 11, or 13 are skipped, mirroring AddCoeff3's behaviour (those values
// go into the assembly RHS `bbc`, not into A).
void ApplyStokesOperatorMDOODZ(grid *mesh, params model,
                               const double *u, const double *v, const double *p,
                               double *ru, double *rv, double *rp);

// Assemble the 5x5 Vanka block at pressure cell (i, j) — the local
// coupling between Vx_L, Vx_R, Vz_B, Vz_T, P_C, with the block's rows
// drawn from the same MDOODZ coefficient formulas as
// `ApplyStokesOperatorMDOODZ` above (design D11: both the matvec and
// the Vanka block dispatch on the fine level must be bit-identical to
// `StokesAssemblyDecoupled.c`).
//
// DOF order inside the block: [Vx_L, Vx_R, Vz_B, Vz_T, P_C] at MDOODZ-
// padded full-grid positions:
//   Vx_L = u[i   + (j+1)*Nx]        Vx_R = u[i+1 + (j+1)*Nx]
//   Vz_B = v[i+1 +  j   *(Nx+1)]    Vz_T = v[i+1 + (j+1)*(Nx+1)]
//   P_C  = p[i   +  j   *(Nx-1)]
//
// Inputs:
//   mesh, model: as for ApplyStokesOperatorMDOODZ.
//   u, v, p    : the current iterate in MDOODZ-padded layout. Used to
//                evaluate `out_b[r]` (the contribution from outside-
//                block neighbours that the caller subtracts from its
//                rhs to form the block's local rhs).
//   i, j       : pressure-cell coordinates, i ∈ [0, Nx-1), j ∈ [0, Nz-1).
//
// Outputs:
//   out_A[5][5]: dense 5×5 block in RAW (no-celvol) scaling, matching
//                the `L->rhs_u` conventions produced by
//                `MultigridStokes.c::BuildLevelRHSFromCompressed`.
//                Callers that want celvol-scaled entries can multiply
//                after the fact.
//   out_b[5]   : for real-equation rows, Σ coeff · state over
//                outside-block neighbours (contributions the caller
//                subtracts from rhs_u/rhs_v/rhs_p at the block DOFs).
//                For rows whose DOF is NOT a real equation (Dirichlet /
//                Neumann / inactive tag in {0, 11, 13, 30, 31, -12}),
//                the row is the identity and `out_b[r]` is the pinned
//                BC value `mesh->BC*.val[dof]`; caller uses this as the
//                block rhs directly.
//
// Returns 0 on success, non-zero if (i, j) is out of range or inputs
// are NULL (caller should skip the cell in that case).
//
// Phase 1 returns the Picard stencil block. Phase 2 (task 15.21) will
// extend this to return the *symmetric Picard part* of the 5×5 block
// in Newton mode (per design D7 — Vanka smoother stability).
int StokesCellBlockMDOODZ(grid *mesh, params model,
                          const double *u, const double *v, const double *p,
                          int i, int j,
                          double out_A[5][5], double out_b[5]);

#ifdef __cplusplus
}
#endif
