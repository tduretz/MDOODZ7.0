// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2023-2026 MDOODZ Developers team
//
// All rights reserved. This file is part of MDOODZ.
//
// MDOODZ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License version 3
// as published by the Free Software Foundation.
// =========================================================================
//
// Geometric-multigrid Stokes solver for the add-gmg-stokes-solver change:
//   - Variable-viscosity Picard Stokes matvec on a staggered grid
//   - Vanka 5x5 block smoother (red-black coloured)
//   - Coarsest-level CHOLMOD direct solve (assembled CSC)
//   - Recursive V-cycle driver
//   - FGMRES outer iteration with V-cycle preconditioning
//   - Adapter to MDOODZ's lin_solver = 3 dispatch (best-effort; see note
//     in MultigridStokes.c on stencil fidelity vs. StokesAssemblyDecoupled.c)

#ifndef MDOODZ_MULTIGRID_STOKES_H
#define MDOODZ_MULTIGRID_STOKES_H

#include "mdoodz-private.h"
#include "MultigridLevels.h"

#ifdef __cplusplus
extern "C" {
#endif

// --- Stokes stencil matvec -----------------------------------------------
//
// Applies the staggered-grid Picard Stokes operator to the level's current
// (Vx, Vz, P) state and stores the result in (out_u, out_v, out_p).
// The three output buffers must have the same shape as Vx, Vz, P.
// BC handling: Dirichlet Vx = Vz = 0 on the outer boundary (rows held as
// identities with zero RHS). Pressure gauge: equation row at cell (0, 0)
// replaced by P(0, 0) = 0.
void StokesApplyA(const MultigridLevel *L,
                  const double *Vx, const double *Vz, const double *P,
                  double *out_u, double *out_v, double *out_p);

// Compute residual r = f - A*x into the level's res_u/res_v/res_p,
// where f is rhs_u/rhs_v/rhs_p and x is Vx/Vz/P.
void StokesResidual(MultigridLevel *L);

// Scale solution/RHS per the pressure gauge before exporting.
double StokesResidualNorm(const MultigridLevel *L);

// --- Vanka smoother (Section 4) ------------------------------------------
//
// One red-black Vanka sweep (updates Vx/Vz/P in place using the current
// RHS). nu = number of full red-black sweeps to perform.
void VankaSweep(MultigridLevel *L, int nu);

// Single-cell 5x5 assemble + solve, exposed for unit tests.
// Returns 0 on success, non-zero on singular block.
int  VankaBlockAssembleSolve(const MultigridLevel *L,
                             const double *Vx, const double *Vz, const double *P,
                             const double *rhs_u, const double *rhs_v, const double *rhs_p,
                             int i, int j,
                             double block[5][5], double rhs[5], double sol[5]);

// --- V-cycle driver (Section 7) ------------------------------------------

// Run one V-cycle starting at level k (0 = finest). Uses the stored rhs_*
// at level k and updates Vx/Vz/P at level k in place with an additive
// correction from the recursion.
// nu_pre/nu_post taken from H; pre-smoothing sweeps at the current level,
// then restriction, recursion, prolongation, post-smoothing.
void Vcycle(MultigridHierarchy *H, int k);

// --- Coarse-level direct solve (Section 6) -------------------------------

// Build the assembled CSC matrix for the coarsest level (n_levels - 1) and
// factorise with CHOLMOD. Cached on the hierarchy (invalidated via
// H->coarse_factor_valid).
int  CoarseAssembleAndFactor(MultigridHierarchy *H);

// Solve the coarsest level: rhs_u/rhs_v/rhs_p → Vx/Vz/P at the coarsest level.
int  CoarseSolve(MultigridHierarchy *H);

// --- FGMRES outer iteration (Section 8) ----------------------------------

typedef struct _fgmres_stats {
    int    iterations;
    int    restarts;
    double initial_residual;  // ‖r_0‖ = ‖b - A*x_0‖ (gauges the relative-reduction tolerance)
    double final_residual;    // ‖r_k‖ at exit; compare to tol via final_residual / initial_residual
    int    fell_back_to_cholmod;
} FGMRESStats;

// Run FGMRES with V-cycle preconditioning on the finest level. Uses the
// level's rhs_u/rhs_v/rhs_p as RHS, Vx/Vz/P as initial guess AND result.
// Returns 0 on convergence within max_restarts * restart iterations, 1
// if it hit the iteration budget (caller decides whether to fall back).
int  FGMRES_GMG(MultigridHierarchy *H, double tol, int restart, int max_restarts,
                FGMRESStats *out_stats);

// --- MDOODZ <-> GMG bridge (Section 10, task 10) -------------------------
//
// Populate the fine GMG level from an MDOODZ mesh:
//   - viscosity fields (etan <- mesh->eta_n, etas <- mesh->eta_s);
//   - velocity/pressure initial guess (Vx/Vz/P <- mesh->u_in/v_in/p_in);
//   - active masks via BuildActiveMask;
//   - boundary velocities for Dirichlet faces (outer ring of Vx/Vz).
// Returns 0 on success, non-zero if shapes don't match.
int PopulateLevelFromMesh(MultigridLevel *L, const grid *mesh);

// Build the level's RHS from mesh body forces and continuity source.
// rhs_u[i] <- mesh->roger_x[i_in_mesh] for interior rows, or the Dirichlet
// boundary value at the outer ring. rhs_p[i] <- mesh->rhs_p[i].
int BuildMeshStokesRHS(MultigridLevel *L, const grid *mesh);

// Copy the level's current Vx/Vz/P back into the mesh's u_in/v_in/p_in
// arrays. Only overwrites entries whose DOF was "solved" (non-30 tag);
// entries with fixed boundary values are left untouched so the caller
// still sees their own BC values.
int UnpackLevelToMesh(const MultigridLevel *L, grid *mesh);

// --- MDOODZ dispatch (Section 10) ----------------------------------------
//
// Top-level adapter called from StokesRoutines.c when lin_solver = 3.
// Extracts the coupled Stokes state from the MDOODZ mesh, runs the GMG
// solver, and copies the resulting velocity/pressure back into the caller-
// provided RHS/solution vectors.
//
// Returns 0 on success, non-zero on failure (caller should fall back).
int  SolveStokesGMG(SparseMat *StokesA, SparseMat *StokesB,
                    SparseMat *StokesC, SparseMat *StokesD,
                    SparseMat *Stokes,  DirectSolver *PardisoStokes,
                    double *rhs_u, double *rhs_p, double *x,
                    params model, grid *mesh, scale scaling);

// --- V-cycle dump (test-only helpers) ------------------------------------
//
// These three helpers expose the otherwise file-static V-cycle dumper
// so unit tests can arm / finish / reset independently of
// `SolveStokesGMG`. Production code only uses the dump through
// `gmg_dump_vcycle = 1` in the input file, which arms internally via
// `SolveStokesGMG`. See design D9 of add-gmg-stokes-defence.

// Reset the one-shot guard. Production never calls this; tests call it
// between arms to exercise multiple dumps per process.
void GmgVcycleDumpReset(void);

// Arm the dumper with an explicit output directory (created recursively
// if needed). If the one-shot guard has already fired, this is a no-op
// and returns 0. Returns 1 when armed, 0 otherwise.
int  GmgVcycleDumpArmForTest(const char *output_dir);

// Mark the current recorded cycle as finished and latch the one-shot
// guard. Safe to call when not armed. Returns the number of snapshots
// written during the just-closed recording (0 if not armed).
int  GmgVcycleDumpFinishForTest(void);

#ifdef __cplusplus
}
#endif

#endif // MDOODZ_MULTIGRID_STOKES_H
