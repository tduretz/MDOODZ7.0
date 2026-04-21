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
//
// MDOODZ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MDOODZ. If not, see <https://www.gnu.org/licenses/>.
// =========================================================================
//
// Geometric-multigrid hierarchy + transfer operators for the
// add-gmg-stokes-solver change.
//
// The GMG module operates on a self-contained staggered-grid data
// structure (MultigridLevel / MultigridHierarchy) with a standard
// variable-viscosity Picard Stokes discretization à la Gerya 2010.
// MDOODZ integration lives in MultigridStokes.c (adapter layer at the
// lin_solver = 3 dispatch seam in StokesRoutines.c).
//
// Grid layout per level (MAC staggered, 2D):
//   Nx  = number of nodes along x (vertex grid)
//   Nz  = number of nodes along z (vertex grid)
//   Ncx = Nx - 1 (cell count along x)
//   Ncz = Nz - 1 (cell count along z)
//
//   Vx:   sized (Nx    ) * (Ncz   )   -- horizontal faces between vertices along z
//   Vz:   sized (Ncx   ) * (Nz    )   -- vertical   faces between vertices along x
//   P:    sized (Ncx   ) * (Ncz   )   -- cell-centred
//   etan: sized (Ncx   ) * (Ncz   )   -- cell-centred viscosity
//   etas: sized (Nx    ) * (Nz    )   -- vertex      viscosity
//
// Index convention: row-major, j (z) slowest, i (x) fastest.
//   idx_P    (i, j) = i + j * Ncx
//   idx_Vx   (i, j) = i + j * Nx           // i in [0, Nx), j in [0, Ncz)
//   idx_Vz   (i, j) = i + j * Ncx          // i in [0, Ncx), j in [0, Nz)
//   idx_etan (i, j) = i + j * Ncx
//   idx_etas (i, j) = i + j * Nx

#ifndef MDOODZ_MULTIGRID_LEVELS_H
#define MDOODZ_MULTIGRID_LEVELS_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

// Coarsest-level DOF floor (coarsen no further once we are below this).
#ifndef MDOODZ_GMG_COARSE_DOF_FLOOR
#define MDOODZ_GMG_COARSE_DOF_FLOOR 5000
#endif
// Minimum per-side cell count before we bail out into the coarse solve.
#ifndef MDOODZ_GMG_MIN_SIDE
#define MDOODZ_GMG_MIN_SIDE 16
#endif
// Viscosity ratio above which restriction switches from arithmetic to harmonic.
#ifndef MDOODZ_GMG_HARMONIC_RATIO
#define MDOODZ_GMG_HARMONIC_RATIO 10.0
#endif

// Inline index helpers (header-only so unit tests can use them directly).
static inline size_t MdIdx_P   (int i, int j, int Ncx) { return (size_t)i + (size_t)j * (size_t)Ncx; }
static inline size_t MdIdx_Vx  (int i, int j, int Nx ) { return (size_t)i + (size_t)j * (size_t)Nx;  }
static inline size_t MdIdx_Vz  (int i, int j, int Ncx) { return (size_t)i + (size_t)j * (size_t)Ncx; }
static inline size_t MdIdx_etan(int i, int j, int Ncx) { return (size_t)i + (size_t)j * (size_t)Ncx; }
static inline size_t MdIdx_etas(int i, int j, int Nx ) { return (size_t)i + (size_t)j * (size_t)Nx;  }

// One level of the GMG hierarchy.
typedef struct _multigrid_level {
    int    Nx, Nz;                  // vertex grid
    int    Ncx, Ncz;                // cell grid (Nx-1, Nz-1)
    double dx, dz;                  // grid spacing
    int    level;                   // 0 = finest

    // Field buffers owned by the level
    double *Vx;                     // [Nx  * Ncz]
    double *Vz;                     // [Ncx * Nz ]
    double *P;                      // [Ncx * Ncz]
    double *rhs_u;                  // momentum-x RHS, same shape as Vx
    double *rhs_v;                  // momentum-z RHS, same shape as Vz
    double *rhs_p;                  // continuity   RHS, same shape as P
    // Work buffers used by the V-cycle (residuals, corrections)
    double *res_u, *res_v, *res_p;  // residual
    double *cor_u, *cor_v, *cor_p;  // correction

    // Viscosity fields (restricted from the fine level above)
    double *etan;                   // [Ncx * Ncz]
    double *etas;                   // [Nx  * Nz ]

    // Active mask (tag-30 handling). 1 = active, 0 = inactive.
    // One mask per staggered variable so Vx/Vz face rows and P cell rows
    // can be independently masked.
    char   *act_P;                  // [Ncx * Ncz]
    char   *act_Vx;                 // [Nx  * Ncz]
    char   *act_Vz;                 // [Ncx * Nz ]

    // MDOODZ stencil bridge (design D11, tasks 15.7/15.8). Originally a
    // single `use_mdoodz_matvec` toggle gated both the outer matvec and
    // the Vanka 5×5 block assembly; add-gmg-upleg-fix §4 localised a
    // divergent L0 pre-smoother to the Vanka bridge path and split the
    // flag in two so the outer matvec can keep the MDOODZ bit-identical
    // operator (for FGMRES fidelity) while the Vanka block falls back to
    // the textbook path (which damps correctly on L1+).
    //
    // `use_mdoodz_matvec_outer` — when non-zero, `StokesApplyA` on this
    //   level dispatches to `StokesApplyA_MDOODZ_bridge`. Set only on
    //   level 0 from `SolveStokesGMG`.
    // `use_mdoodz_matvec_vanka` — when non-zero, `VankaBlockAssembleSolve`
    //   on this level dispatches to `VankaBlockAssembleSolve_MDOODZ_bridge`
    //   (which uses `StokesCellBlockMDOODZ`). **Currently left zero**
    //   everywhere pending a root-cause fix in the bridge 5×5 block; see
    //   openspec/changes/add-gmg-upleg-fix/STATUS.md "Finding 3" and
    //   "Fix candidate F" for the rationale.
    // Coarse levels and stand-alone unit tests leave both flags zero and
    // use the textbook Picard Stokes stencil (D4 explicitly allows
    // coarse-level Picard-on-restricted-viscosity rediscretisation).
    int           use_mdoodz_matvec_outer;
    int           use_mdoodz_matvec_vanka;
    void         *mesh_ref;         // `grid*`; kept opaque to avoid circular include
    void         *model_ref;        // `params*`; kept opaque to avoid circular include

    // Scratch buffers for the MDOODZ-padded state used by the fine-level
    // Vanka bridge (task 15.8). `StokesCellBlockMDOODZ` expects (u, v, p)
    // in MDOODZ's Nx*(Nz+1)/(Nx+1)*Nz/Ncx*Ncz padded layout, not the
    // strict-interior layout stored in L->Vx/Vz/P. `VankaSweep` packs
    // L->Vx/Vz/P into these buffers once at the start of each sweep and
    // keeps them in sync with each cell's per-DOF update — that way the
    // Gauss-Seidel freshness of the textbook sweep is preserved without
    // repacking the entire grid per cell. Allocated lazily from
    // `VankaSweep` when `use_mdoodz_matvec_vanka` is non-zero; freed by
    // `FreeMultigridHierarchy`. NULL otherwise.
    double       *md_u_pad;         // [Nx   * (Nz + 1)]
    double       *md_v_pad;         // [(Nx + 1) * Nz ]
    double       *md_p_pad;         // [Ncx  * Ncz]  (same as P layout)
} MultigridLevel;

// Hierarchy container. Level 0 is finest; levels[n_levels - 1] is the
// coarsest, which is handed to UMFPACK for a direct factor/solve.
// (The Stokes saddle-point is indefinite, so CHOLMOD is not suitable;
// UMFPACK's general-unsymmetric LU handles it cleanly.)
typedef struct _multigrid_hierarchy {
    int             n_levels;
    MultigridLevel *levels;         // [n_levels]

    // Coarsest-level UMFPACK state (lazy-initialised on first coarse solve,
    // invalidated when viscosity is restricted afresh or masks change).
    // Stored as opaque void* so the header needn't include <umfpack.h>.
    void           *coarse_umf_numeric;    // void* from umfpack_di_numeric
    int             coarse_factor_valid;

    // Coarsest-level assembled CSR/CSC matrix (kept so we can re-factorise
    // cheaply when viscosity changes structurally). Stored in "UMFPACK-CSR"
    // form: Ap = row pointers of length n+1, Aj = column indices of length nnz,
    // Ax = values. UMFPACK is then called with UMFPACK_Aat (A-transpose-solve)
    // to interpret this CSR layout as a CSC of A^T (the standard MDOODZ trick).
    int            *coarse_Ap;             // row pointers (size n+1)
    int            *coarse_Aj;             // column indices (size nnz)
    double         *coarse_Ax;             // values (size nnz)
    int             coarse_n;              // matrix dimension
    int             coarse_nnz;
    // Map from (dof_kind, local_index) to global coarse DOF index, and
    // the inverse — filled in when the coarse system is assembled.
    // kind 0 = Vx, 1 = Vz, 2 = P. -1 means "not a DOF" (inactive or on
    // the outer boundary, which is held as an identity row).
    int            *coarse_eqn_Vx;         // [Nx_c  * Ncz_c]
    int            *coarse_eqn_Vz;         // [Ncx_c * Nz_c ]
    int            *coarse_eqn_P;          // [Ncx_c * Ncz_c]

    // Tuning (copied from params at construction time).
    int    nu_pre, nu_post;         // Vanka sweeps per level

    // Newton/Jacobian flag (task 9.3). When set, downstream hooks may
    // switch to a Jacobian discretization. The current solver uses Picard
    // for the smoother and coarse operator on every call; the flag is
    // plumbed end-to-end so Section 9 extensions can specialise without a
    // second dispatch layer.
    int    use_jacobian;
} MultigridHierarchy;

// --- Hierarchy lifetime --------------------------------------------------

// Compute the number of levels for a given fine-grid size, honouring
// the gmg_levels override (0 = auto).
int  MultigridComputeDefaultLevelCount(int Nx, int Nz, int gmg_levels_override);

// Allocate a complete hierarchy from a fine-level (Nx, Nz, dx, dz) and
// a gmg_levels override. Each coarser level has ceil(N/2) cells per side.
// Returns 0 on success, non-zero on error.
int  MultigridHierarchyAllocate(MultigridHierarchy *H, int Nx_fine, int Nz_fine,
                                double dx_fine, double dz_fine,
                                int gmg_levels_override);

// Free all level buffers, coarse CSC storage, and CHOLMOD state.
void MultigridHierarchyFree(MultigridHierarchy *H);

// --- Transfer operators --------------------------------------------------
// These operate between two adjacent levels; "fine" is the finer level,
// "coarse" is the coarser level (one step deeper in the hierarchy).

// Full-weighting restriction of a velocity component (Vx or Vz).
// component = 0 -> Vx, component = 1 -> Vz.
void RestrictVelocity(const MultigridLevel *fine, MultigridLevel *coarse,
                      const double *fine_V, double *coarse_V, int component);

// Full-weighting restriction of the pressure field.
void RestrictPressure(const MultigridLevel *fine, MultigridLevel *coarse,
                      const double *fine_P, double *coarse_P);

// Bilinear prolongation of a velocity component (Vx or Vz) with *addition*
// to the destination (additive correction step of a V-cycle).
void ProlongateVelocityAdd(const MultigridLevel *coarse, const MultigridLevel *fine,
                           const double *coarse_V, double *fine_V, int component);

// Piecewise-constant (injection) prolongation of pressure, additive.
void ProlongatePressureAdd(const MultigridLevel *coarse, const MultigridLevel *fine,
                           const double *coarse_P, double *fine_P);

// Restrict viscosity fields (cell-centred etan + vertex-based etas) from
// fine to coarse level. Uses harmonic mean where the local ratio exceeds
// MDOODZ_GMG_HARMONIC_RATIO, arithmetic mean otherwise. Inactive children
// are skipped.
void RestrictViscosity(const MultigridLevel *fine, MultigridLevel *coarse);

// Restrict the active mask: a coarse cell is active iff at least one of
// its fine children is active.
void RestrictActiveMask(const MultigridLevel *fine, MultigridLevel *coarse);

// --- Active-mask construction (task 5.1) --------------------------------
//
// Build the per-level active masks from MDOODZ boundary-condition tag
// arrays (BCu.type for Vx, BCv.type for Vz, BCp.type for P). A DOF is
// considered active (mask = 1) unless it carries the "air" tag 30, which
// is the only MDOODZ convention meaning "do not compute". Dirichlet and
// Neumann boundary tags (0, 11, 13, 31, 2, -2, -12) remain "active" here
// because their row is still part of the GMG operator (as an identity
// row at the Dirichlet face, or as a stencil row for Neumann); the Vanka
// block handles each of them via the boundary-row logic in
// VankaBlockAssembleSolve.
//
// The input arrays must match MDOODZ's mesh layout:
//   BCu_type: length Nx * (Nz + 1)   (indexed as [k + l * Nx], l in [0, Nz])
//   BCv_type: length (Nx + 1) * Nz
//   BCp_type: length (Nx - 1) * (Nz - 1)
// Our fine level stores:
//   act_Vx:   [Nx    * (Nz - 1)]
//   act_Vz:   [(Nx-1)* Nz      ]
//   act_P:    [(Nx-1)*(Nz - 1) ]
// We map MDOODZ's inner Vx rows l = 1 .. Nz-1 onto our j = 0 .. Ncz-1,
// and MDOODZ's inner Vz columns k = 1 .. Nx-1 onto our i = 0 .. Ncx-1.
// BC entries outside those inner ranges encode ghost values and do not
// feed directly into our active-mask arrays; they are read by the
// population helpers when computing Dirichlet boundary values.
void BuildActiveMask(MultigridLevel       *L,
                     const signed char    *BCu_type,
                     const signed char    *BCv_type,
                     const signed char    *BCp_type);

// --- Utility -------------------------------------------------------------

// Fill a whole level with zeros (fields + work buffers; mask untouched).
void LevelZeroFields(MultigridLevel *L);

// 2-norm over active entries of a vector with given mask.
double MaskedL2(const double *x, const char *mask, size_t n);

#ifdef __cplusplus
}
#endif

#endif // MDOODZ_MULTIGRID_LEVELS_H
