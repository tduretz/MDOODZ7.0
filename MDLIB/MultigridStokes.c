// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2023-2026 MDOODZ Developers team
//
// All rights reserved. This file is part of MDOODZ.
// =========================================================================
//
// GMG Stokes solver implementation.
//
// Discretization. Variable-viscosity Picard Stokes on a MAC staggered grid
// (Gerya 2010, chap 14). For cell (i, j) ∈ [0, Ncx) × [0, Ncz):
//
//   mom_x at (i, j+0.5) (Vx(i, j) in array terms, i ∈ [1, Nx-1) interior):
//     ∂/∂x (2 η_n ∂Vx/∂x)  +  ∂/∂z (η_s (∂Vx/∂z + ∂Vz/∂x))  -  ∂P/∂x  =  fx
//   mom_z at (i+0.5, j):
//     ∂/∂z (2 η_n ∂Vz/∂z)  +  ∂/∂x (η_s (∂Vx/∂z + ∂Vz/∂x))  -  ∂P/∂z  =  fz
//   continuity at (i+0.5, j+0.5):
//     ∂Vx/∂x + ∂Vz/∂z = fp
//
// Boundaries: Dirichlet on the outer ring of Vx (i in {0, Nx-1}) and Vz
// (j in {0, Nz-1}); those rows are held as identities with the level's
// rhs_u/rhs_v at the boundary giving the Dirichlet value. Inactive cells
// (tag 30) are likewise identity rows tied to the current iterate.
// Pressure null-space: we do *not* write an explicit gauge row; instead
// the V-cycle projects the constant-pressure component out of both RHS
// and iterate at every level (ProjectPressureMean below) to keep the
// Galerkin-like transfer consistent between fine and coarse.
//
// Coarsest-level solve. We assemble the coarsest stencil into CSR form
// (by evaluating StokesApplyA against local unit vectors), factor it
// with UMFPACK (which handles indefinite systems — CHOLMOD is SPD-only),
// cache the numeric factor on the hierarchy, and reuse it for every
// V-cycle until the hierarchy is freed or the factor is invalidated.

#include "MultigridStokes.h"
#include "mdoodz-log.h"
#include "StokesAssemblyGMG.h"
#include "mdoodz-private.h"   // CreateOutputHDF5 / AddGroupToHDF5 / AddFieldToGroup

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>

#include "umfpack.h"

// Forward declaration: dispatcher into the MDOODZ-stencil matvec for the
// fine level (design D11 / tasks 15.7). Defined below StokesApplyA.
static void StokesApplyA_MDOODZ_bridge(const MultigridLevel *L,
                                       const double *Vx, const double *Vz, const double *P,
                                       double *out_u, double *out_v, double *out_p);

// ---------- utility ------------------------------------------------------

// Inner-face index ranges:
//   Interior Vx rows for mom-x:  i ∈ [1, Nx-1),  j ∈ [0, Ncz)
//   Interior Vz rows for mom-z:  i ∈ [0, Ncx),   j ∈ [1, Nz-1)
// (face positions at the outer boundary are held by Dirichlet rows)

// Read helpers with shape.
static inline double gVx(const double *V, int i, int j, int Nx ) { return V[(size_t)i + (size_t)j * Nx]; }
static inline double gVz(const double *V, int i, int j, int Ncx) { return V[(size_t)i + (size_t)j * Ncx]; }
static inline double gP (const double *P, int i, int j, int Ncx) { return P[(size_t)i + (size_t)j * Ncx]; }

// ---------- Stokes matvec (Section 7.2) ---------------------------------

void StokesApplyA(const MultigridLevel *L,
                  const double *Vx, const double *Vz, const double *P,
                  double *out_u, double *out_v, double *out_p) {
    // Stencil bridge dispatch (design D11, task 15.7). On the fine
    // level of a real SolveStokesGMG run we route the matvec through
    // ApplyStokesOperatorMDOODZ so the fine-level operator is bit-
    // identical to StokesAssemblyDecoupled.c (enforced by the golden
    // cross-check in TESTS/StokesMatvecEquivalence.cpp). Coarse levels
    // and manufactured-RHS unit tests leave the flag at its default
    // zero and continue on the textbook Picard Stokes stencil.
    if (L->use_mdoodz_matvec_outer) {
        StokesApplyA_MDOODZ_bridge(L, Vx, Vz, P, out_u, out_v, out_p);
        return;
    }

    const int Nx = L->Nx,  Nz = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const double dx = L->dx, dz = L->dz;
    const double idx = 1.0 / dx, idz = 1.0 / dz;
    const double idx2 = idx * idx, idz2 = idz * idz;

    const double *etan = L->etan;
    const double *etas = L->etas;

    // --- momentum-x ---
    // Rows for Vx(i, j): identity on boundary (i = 0 or i = Nx-1); interior
    // otherwise. Also identity when inactive.
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Nx; ++i) {
            const size_t idx_u = MdIdx_Vx(i, j, Nx);
            if (i == 0 || i == Nx - 1 || !L->act_Vx[idx_u]) {
                out_u[idx_u] = gVx(Vx, i, j, Nx);
                continue;
            }
            // Stresses around (i, j+0.5):
            const double etan_L = etan[MdIdx_P(i - 1, j, Ncx)];
            const double etan_R = etan[MdIdx_P(i,     j, Ncx)];
            // etas at the face's top (i, j+1) and bottom (i, j) vertex
            // positions in vertex grid (Nx × Nz).
            const double etas_T = etas[MdIdx_etas(i, j + 1, Nx)];
            const double etas_B = etas[MdIdx_etas(i, j,     Nx)];

            // sigma_xx right = 2 etan_R (Vx(i+1, j) - Vx(i, j)) / dx
            // sigma_xx left  = 2 etan_L (Vx(i, j) - Vx(i-1, j)) / dx
            const double sxx_R = 2.0 * etan_R * (gVx(Vx, i + 1, j, Nx) - gVx(Vx, i, j, Nx)) * idx;
            const double sxx_L = 2.0 * etan_L * (gVx(Vx, i,     j, Nx) - gVx(Vx, i - 1, j, Nx)) * idx;

            // sigma_xz top (at (i, j+1)) and bot (at (i, j))
            // Vz indices: Vz(i-1, j) and Vz(i, j) are valid because i ∈ [1, Nx-1)
            //   and so i-1 ∈ [0, Ncx-1), i ∈ [1, Ncx] — but Vz indices in x go [0, Ncx).
            //   When i == Nx - 1 we've already bailed out above, so i < Nx - 1 ⇒ i ≤ Ncx - 1. OK.
            const double Vz_right_T = (i   < Ncx) ? gVz(Vz, i,     j + 1, Ncx) : 0.0;
            const double Vz_left_T  = (i-1 >= 0 ) ? gVz(Vz, i - 1, j + 1, Ncx) : 0.0;
            const double Vz_right_B = (i   < Ncx) ? gVz(Vz, i,     j,     Ncx) : 0.0;
            const double Vz_left_B  = (i-1 >= 0 ) ? gVz(Vz, i - 1, j,     Ncx) : 0.0;

            // (j+1) for Vx top may go out of bounds when j == Ncz - 1. Clamp
            // to the Dirichlet ghost (= 0) by restricting bounds.
            const double Vx_top = (j + 1 < Ncz) ? gVx(Vx, i, j + 1, Nx) : 0.0;
            const double Vx_bot = (j - 1 >= 0 ) ? gVx(Vx, i, j - 1, Nx) : 0.0;

            const double sxz_T = etas_T * ((Vx_top - gVx(Vx, i, j, Nx)) * idz
                                            + (Vz_right_T - Vz_left_T) * idx);
            const double sxz_B = etas_B * ((gVx(Vx, i, j, Nx) - Vx_bot) * idz
                                            + (Vz_right_B - Vz_left_B) * idx);

            // pressure: P(i, j) - P(i-1, j)
            const double PR = gP(P, i,     j, Ncx);
            const double PL = gP(P, i - 1, j, Ncx);

            out_u[idx_u] = (sxx_R - sxx_L) * idx + (sxz_T - sxz_B) * idz - (PR - PL) * idx;
            (void)idx2; (void)idz2; // (computed coefficients inline)
        }
    }

    // --- momentum-z ---
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const size_t idx_v = MdIdx_Vz(i, j, Ncx);
            if (j == 0 || j == Nz - 1 || !L->act_Vz[idx_v]) {
                out_v[idx_v] = gVz(Vz, i, j, Ncx);
                continue;
            }
            const double etan_B = etan[MdIdx_P(i, j - 1, Ncx)];
            const double etan_T = etan[MdIdx_P(i, j,     Ncx)];
            const double etas_R = etas[MdIdx_etas(i + 1, j, Nx)];
            const double etas_L = etas[MdIdx_etas(i,     j, Nx)];

            const double szz_T = 2.0 * etan_T * (gVz(Vz, i, j + 1, Ncx) - gVz(Vz, i, j, Ncx)) * idz;
            const double szz_B = 2.0 * etan_B * (gVz(Vz, i, j,     Ncx) - gVz(Vz, i, j - 1, Ncx)) * idz;

            const double Vx_top_R = (j   < Ncz) ? gVx(Vx, i + 1, j,     Nx) : 0.0;
            const double Vx_bot_R = (j-1 >= 0 ) ? gVx(Vx, i + 1, j - 1, Nx) : 0.0;
            const double Vx_top_L = (j   < Ncz) ? gVx(Vx, i,     j,     Nx) : 0.0;
            const double Vx_bot_L = (j-1 >= 0 ) ? gVx(Vx, i,     j - 1, Nx) : 0.0;

            const double Vz_right = (i + 1 < Ncx) ? gVz(Vz, i + 1, j, Ncx) : 0.0;
            const double Vz_left  = (i - 1 >= 0 ) ? gVz(Vz, i - 1, j, Ncx) : 0.0;

            const double sxz_R = etas_R * ((Vx_top_R - Vx_bot_R) * idz
                                            + (Vz_right - gVz(Vz, i, j, Ncx)) * idx);
            const double sxz_L = etas_L * ((Vx_top_L - Vx_bot_L) * idz
                                            + (gVz(Vz, i, j, Ncx) - Vz_left) * idx);

            const double PT = gP(P, i, j,     Ncx);
            const double PB = gP(P, i, j - 1, Ncx);

            out_v[idx_v] = (szz_T - szz_B) * idz + (sxz_R - sxz_L) * idx - (PT - PB) * idz;
        }
    }

    // --- continuity ---
    //
    // No pressure-gauge identity row here. The Stokes system is singular
    // with a 1-D null space {constant pressure}; the multigrid driver
    // projects this null space out of the RHS and the iterate once per
    // V-cycle (ProjectPressureMean below). This keeps Galerkin consistency
    // between fine and coarse operators, which an identity row would
    // break because coarse cell (0,0) is a different physical location
    // than fine cell (0,0).
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const size_t idx_p = MdIdx_P(i, j, Ncx);
            if (!L->act_P[idx_p]) { out_p[idx_p] = gP(P, i, j, Ncx); continue; }
            const double dudx = (gVx(Vx, i + 1, j, Nx) - gVx(Vx, i, j, Nx)) * idx;
            const double dwdz = (gVz(Vz, i, j + 1, Ncx) - gVz(Vz, i, j, Ncx)) * idz;
            out_p[idx_p] = dudx + dwdz;
        }
    }
}

// ---------------------------------------------------------------------------
// MDOODZ-stencil matvec bridge (design D11 / task 15.7).
//
// This translates the GMG-level staggered arrays (Vx: Nx×Ncz, Vz: Ncx×Nz,
// P: Ncx×Ncz) into the MDOODZ-mesh padded layout (u_in: Nx×(Nz+1),
// v_in: (Nx+1)×Nz, p_in: Ncx×Ncz), calls the golden-test-validated
// `ApplyStokesOperatorMDOODZ` on the fine-level mesh, and writes the
// result back into GMG layout. Scaling: `ApplyStokesOperatorMDOODZ`
// returns `celvol · Σ coef · x` (matching the assembled sparse matrix),
// while GMG expects raw PDE scaling, so we divide by celvol on output.
// Identity rows (boundary/inactive DOFs) are held against the input
// iterate exactly as the textbook matvec does so `StokesResidual` sees
// a zero residual there (since `BuildMeshStokesRHS` stores the iterate
// itself at those locations).
//
// Ghost rows of u_tmp / columns of v_tmp are left at zero. The MDOODZ
// stencil skips any neighbour tagged 0, 11, 13, 31 via `gmg_skip_neighbour`,
// which for a well-formed MDOODZ mesh covers every ghost face, so their
// contribution to A·x is zero regardless of the fill value. This keeps
// the bridge stateless with respect to mesh->u_in / mesh->v_in ghost
// Dirichlet values (which are already folded into `mesh->roger_x/z`
// through `bbc` at assembly time).
static void StokesApplyA_MDOODZ_bridge(const MultigridLevel *L,
                                       const double *Vx, const double *Vz, const double *P,
                                       double *out_u, double *out_v, double *out_p) {
    grid   *mesh  = (grid   *)L->mesh_ref;
    params *model = (params *)L->model_ref;
    if (mesh == NULL || model == NULL) {
        LOG_ERR("GMG: bridge dispatch invoked without mesh/model refs; "
                "falling back to zero output");
        const size_t nVx = (size_t)L->Nx  * L->Ncz;
        const size_t nVz = (size_t)L->Ncx * L->Nz;
        const size_t nP  = (size_t)L->Ncx * L->Ncz;
        memset(out_u, 0, nVx * sizeof(double));
        memset(out_v, 0, nVz * sizeof(double));
        memset(out_p, 0, nP  * sizeof(double));
        return;
    }

    const int Nx  = L->Nx,  Nz  = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const int Nzvx = Nz + 1;
    const int Nxvz = Nx + 1;
    const double celvol = L->dx * L->dz;
    const double inv_celvol = 1.0 / celvol;

    const size_t nu_full = (size_t)Nx   * (size_t)Nzvx;
    const size_t nv_full = (size_t)Nxvz * (size_t)Nz;
    const size_t np_full = (size_t)Ncx  * (size_t)Ncz;

    double *u_tmp  = (double *)calloc(nu_full, sizeof(double));
    double *v_tmp  = (double *)calloc(nv_full, sizeof(double));
    double *p_tmp  = (double *)calloc(np_full, sizeof(double));
    double *ru_tmp = (double *)calloc(nu_full, sizeof(double));
    double *rv_tmp = (double *)calloc(nv_full, sizeof(double));
    double *rp_tmp = (double *)calloc(np_full, sizeof(double));
    if (!u_tmp || !v_tmp || !p_tmp || !ru_tmp || !rv_tmp || !rp_tmp) {
        LOG_ERR("GMG: bridge allocation failed (Nx=%d Nz=%d)", Nx, Nz);
        free(u_tmp); free(v_tmp); free(p_tmp);
        free(ru_tmp); free(rv_tmp); free(rp_tmp);
        return;
    }

    // --- pack: GMG layout -> MDOODZ padded layout ---
    for (int j = 0; j < Ncz; ++j) {
        const int l = j + 1;              // GMG j∈[0,Ncz) maps to MDOODZ l∈[1,Nz]
        for (int i = 0; i < Nx; ++i) {
            u_tmp[(size_t)i + (size_t)l * Nx] = Vx[MdIdx_Vx(i, j, Nx)];
        }
    }
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const int k = i + 1;          // GMG i∈[0,Ncx) maps to MDOODZ k∈[1,Nx]
            v_tmp[(size_t)k + (size_t)j * Nxvz] = Vz[MdIdx_Vz(i, j, Ncx)];
        }
    }
    memcpy(p_tmp, P, np_full * sizeof(double));

    // --- MDOODZ matvec ---
    ApplyStokesOperatorMDOODZ(mesh, *model, u_tmp, v_tmp, p_tmp, ru_tmp, rv_tmp, rp_tmp);

    // --- unpack: MDOODZ padded layout -> GMG layout, descaled by celvol.
    //            Boundary/inactive rows held as identity against the input.
    //
    // Sign convention (design D11 + D4 reconciliation, task 15.8
    // follow-up). MDOODZ's assembled matrix rows have the form
    // `-∇·σ + ∇p = f` — i.e. the Vx_C / Vz_C diagonal is POSITIVE. The
    // textbook matvec in `StokesApplyA` above instead discretises
    // `+∇·σ - ∇p = f` with NEGATIVE diagonal. Both are valid —
    // they just differ in sign on momentum rows. When the bridge is
    // engaged on the fine level but coarse levels keep the textbook
    // stencil (per D4), the V-cycle correction (P·A_c^{-1}·R·r_f) is
    // only stable if A_c and A_f share a sign. The cheapest
    // reconciliation is to negate the bridge's momentum output here
    // so the fine level speaks textbook convention too; the
    // caller-supplied rhs_u / rhs_v are negated in
    // `BuildLevelRHSFromCompressed` to match. Continuity (out_p) has
    // the SAME sign in both conventions (`∇·u = 0`) and is not
    // flipped. Identity rows (bnd / inactive) are also unflipped —
    // they pin the DOF value and are sign-invariant.
    for (int j = 0; j < Ncz; ++j) {
        const int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            const size_t idx_u = MdIdx_Vx(i, j, Nx);
            if (i == 0 || i == Nx - 1 || !L->act_Vx[idx_u]) {
                out_u[idx_u] = Vx[idx_u];
            } else {
                out_u[idx_u] = -ru_tmp[(size_t)i + (size_t)l * Nx] * inv_celvol;
            }
        }
    }
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const size_t idx_v = MdIdx_Vz(i, j, Ncx);
            if (j == 0 || j == Nz - 1 || !L->act_Vz[idx_v]) {
                out_v[idx_v] = Vz[idx_v];
            } else {
                const int k = i + 1;
                out_v[idx_v] = -rv_tmp[(size_t)k + (size_t)j * Nxvz] * inv_celvol;
            }
        }
    }
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const size_t idx_p = MdIdx_P(i, j, Ncx);
            if (!L->act_P[idx_p]) {
                out_p[idx_p] = P[idx_p];
            } else {
                out_p[idx_p] = rp_tmp[idx_p] * inv_celvol;
            }
        }
    }

    free(u_tmp); free(v_tmp); free(p_tmp);
    free(ru_tmp); free(rv_tmp); free(rp_tmp);
}

// Remove the constant-pressure null-space component from an array sized
// Ncx*Ncz using the active mask as the support. Idempotent.
static void project_pressure_mean(double *Pfield, const char *mask, int Ncx, int Ncz) {
    double sum = 0.0;
    long count = 0;
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            if (mask != NULL && !mask[MdIdx_P(i, j, Ncx)]) continue;
            sum += Pfield[MdIdx_P(i, j, Ncx)];
            ++count;
        }
    }
    if (count == 0) return;
    double mean = sum / (double)count;
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            if (mask != NULL && !mask[MdIdx_P(i, j, Ncx)]) continue;
            Pfield[MdIdx_P(i, j, Ncx)] -= mean;
        }
    }
}

void StokesResidual(MultigridLevel *L) {
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    StokesApplyA(L, L->Vx, L->Vz, L->P, L->res_u, L->res_v, L->res_p);
    for (size_t k = 0; k < nVx; ++k) L->res_u[k] = L->rhs_u[k] - L->res_u[k];
    for (size_t k = 0; k < nVz; ++k) L->res_v[k] = L->rhs_v[k] - L->res_v[k];
    for (size_t k = 0; k < nP ; ++k) L->res_p[k] = L->rhs_p[k] - L->res_p[k];
}

double StokesResidualNorm(const MultigridLevel *L) {
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    double s = 0.0;
    for (size_t k = 0; k < nVx; ++k) s += L->res_u[k] * L->res_u[k];
    for (size_t k = 0; k < nVz; ++k) s += L->res_v[k] * L->res_v[k];
    for (size_t k = 0; k < nP ; ++k) s += L->res_p[k] * L->res_p[k];
    return sqrt(s);
}

// ---------- Vanka 5x5 block (Section 4) ---------------------------------
//
// DOF order per block: [Vx_L, Vx_R, Vz_B, Vz_T, P_C]
// For cell (i, j):
//   Vx_L = Vx(i,   j)
//   Vx_R = Vx(i+1, j)
//   Vz_B = Vz(i, j)
//   Vz_T = Vz(i, j+1)
//   P_C  = P(i, j)

// Solve a general 5x5 system A*x = b by Gaussian elimination with partial
// pivoting. A is modified in place. Returns 0 on success, non-zero on
// singular matrix.
static int solve5x5(double A[5][5], double b[5], double x[5]) {
    int piv[5] = {0, 1, 2, 3, 4};
    for (int k = 0; k < 5; ++k) {
        // Pivot selection
        int maxrow = k;
        double maxval = fabs(A[piv[k]][k]);
        for (int r = k + 1; r < 5; ++r) {
            double v = fabs(A[piv[r]][k]);
            if (v > maxval) { maxval = v; maxrow = r; }
        }
        if (maxval < 1e-300) return -1;
        if (maxrow != k) { int tmp = piv[k]; piv[k] = piv[maxrow]; piv[maxrow] = tmp; }
        // Eliminate
        double diag = A[piv[k]][k];
        for (int r = k + 1; r < 5; ++r) {
            double f = A[piv[r]][k] / diag;
            A[piv[r]][k] = f;
            for (int c = k + 1; c < 5; ++c) {
                A[piv[r]][c] -= f * A[piv[k]][c];
            }
            b[piv[r]] -= f * b[piv[k]];
        }
    }
    // Back-substitute
    for (int k = 4; k >= 0; --k) {
        double s = b[piv[k]];
        for (int c = k + 1; c < 5; ++c) s -= A[piv[k]][c] * x[c];
        x[k] = s / A[piv[k]][k];
    }
    return 0;
}

// Fine-level MDOODZ stencil bridge for the 5×5 Vanka block (task 15.8 /
// design D11). Dispatches into `StokesCellBlockMDOODZ` so that the Vanka
// smoother at the fine level consumes the SAME coefficient formulas
// used by the matvec (`StokesApplyA_MDOODZ_bridge` → `xmom/zmom/cont_
// fill_coeffs`) — eliminating any drift between the operator that
// FGMRES/Richardson sees and the preconditioner's block solve. The
// block is built in raw-PDE scaling (no celvol) to match the level's
// `rhs_u/rhs_v/rhs_p` produced by `BuildLevelRHSFromCompressed`.
//
// Preconditions (enforced by guards):
//   L->use_mdoodz_matvec_vanka is non-zero
//   L->mesh_ref, L->model_ref, L->md_u_pad, L->md_v_pad, L->md_p_pad
//   are all non-NULL (allocated by `VankaSweep` on entry).
static int VankaBlockAssembleSolve_MDOODZ_bridge(
        const MultigridLevel *L,
        const double *rhs_u, const double *rhs_v, const double *rhs_p,
        int i, int j,
        double block[5][5], double rhs[5], double sol[5]) {
    grid   *mesh  = (grid   *)L->mesh_ref;
    params *model = (params *)L->model_ref;
    if (mesh == NULL || model == NULL) return -1;
    if (L->md_u_pad == NULL || L->md_v_pad == NULL || L->md_p_pad == NULL) return -1;

    const int Nx = L->Nx,  Nz = L->Nz, Ncx = L->Ncx;
    if (i < 0 || i >= L->Ncx || j < 0 || j >= L->Ncz) return -1;
    if (!L->act_P [MdIdx_P (i,     j, Ncx)]) return -1;
    if (!L->act_Vx[MdIdx_Vx(i,     j, Nx )]) return -1;
    if (!L->act_Vx[MdIdx_Vx(i + 1, j, Nx )]) return -1;
    if (!L->act_Vz[MdIdx_Vz(i, j,     Ncx)]) return -1;
    if (!L->act_Vz[MdIdx_Vz(i, j + 1, Ncx)]) return -1;

    double out_A[5][5];
    double out_b[5];
    int rc = StokesCellBlockMDOODZ(mesh, *model,
                                   L->md_u_pad, L->md_v_pad, L->md_p_pad,
                                   i, j, out_A, out_b);
    if (rc != 0) return rc;

    for (int r = 0; r < 5; ++r)
        for (int c = 0; c < 5; ++c)
            block[r][c] = out_A[r][c];

    // Classify rows: real-equation rows consume the caller's rhs minus
    // the external contribution; identity rows (Dirichlet/Neumann/
    // inactive DOF per MDOODZ BC tags) are pinned to the BC value
    // returned in `out_b`. This mirrors `StokesCellBlockMDOODZ`'s own
    // identity-row encoding, and guarantees `block * sol = rhs` yields
    // `sol[r] = bc_val` for the identity rows.
    const int Nxvz = Nx + 1;
    const int iVxL = i       + (j + 1) * Nx;
    const int iVxR = (i + 1) + (j + 1) * Nx;
    const int iVzB = (i + 1) +  j      * Nxvz;
    const int iVzT = (i + 1) + (j + 1) * Nxvz;
    const int iPC  = i       +  j      * Ncx;
    const signed char tVxL = mesh->BCu.type[iVxL];
    const signed char tVxR = mesh->BCu.type[iVxR];
    const signed char tVzB = mesh->BCv.type[iVzB];
    const signed char tVzT = mesh->BCv.type[iVzT];
    const signed char tPC  = mesh->BCp.type[iPC];
    const int real[5] = {
        (tVxL == -1 || tVxL == 2 || tVxL == -2) ? 1 : 0,
        (tVxR == -1 || tVxR == 2 || tVxR == -2) ? 1 : 0,
        (tVzB == -1 || tVzB == 2 || tVzB == -2) ? 1 : 0,
        (tVzT == -1 || tVzT == 2 || tVzT == -2) ? 1 : 0,
        (tPC  == -1) ? 1 : 0,
    };
    // rhs_u / rhs_v arrive in TEXTBOOK convention (negated momentum
    // rows — see `BuildLevelRHSFromCompressed` and the sign comment in
    // `StokesApplyA_MDOODZ_bridge`). `out_A` / `out_b` from
    // `StokesCellBlockMDOODZ` are in MDOODZ convention (positive Vx
    // diagonal, coupling carries positive contributions). To solve
    // `A_text · x = rhs_text - coupling_text · x_outside` using the
    // MDOODZ block (`A_mdoodz = -A_text`, `coupling_mdoodz = -coupling_text`):
    //   A_mdoodz · x = -rhs_text - out_b
    // for momentum real rows (r < 4). Continuity (r = 4) has identical
    // sign in both conventions and is untouched.
    const double rhs_dof[5] = {
        -rhs_u[MdIdx_Vx(i,     j, Nx)],
        -rhs_u[MdIdx_Vx(i + 1, j, Nx)],
        -rhs_v[MdIdx_Vz(i, j,     Ncx)],
        -rhs_v[MdIdx_Vz(i, j + 1, Ncx)],
         rhs_p[MdIdx_P (i,     j, Ncx)],
    };
    for (int r = 0; r < 5; ++r) {
        rhs[r] = real[r] ? (rhs_dof[r] - out_b[r]) : out_b[r];
    }
    return solve5x5(block, rhs, sol);
}

// Build a 5x5 Vanka block for cell (i, j). The block matrix holds the
// coupling between the five local DOFs; contributions from neighbouring
// DOFs are evaluated against the current (Vx, Vz, P) state and subtracted
// from the RHS so that solving block * delta = rhs_in_block updates the
// five DOFs in a consistent Gauss-Seidel-like step.
//
// Returns 0 on success, non-zero if cell is near a boundary where the
// block cannot be built (caller should skip).
int VankaBlockAssembleSolve(const MultigridLevel *L,
                            const double *Vx, const double *Vz, const double *P,
                            const double *rhs_u, const double *rhs_v, const double *rhs_p,
                            int i, int j,
                            double block[5][5], double rhs[5], double sol[5]) {
    // MDOODZ stencil bridge dispatch (design D11, task 15.8). Gated
    // by `use_mdoodz_matvec_vanka` (was `use_mdoodz_matvec` pre
    // add-gmg-upleg-fix; the split is explained in MultigridLevels.h).
    // Currently SolveStokesGMG leaves this flag at zero even on L0
    // because the bridge 5×5 block caused a 7.55× per-sweep L0
    // pre-smooth amplification (STATUS.md Finding 3). The dispatch is
    // kept here so the bridge path can be re-enabled once its root
    // cause is fixed in a follow-up change.
    if (L->use_mdoodz_matvec_vanka && L->md_u_pad != NULL) {
        (void)Vx; (void)Vz; (void)P;
        return VankaBlockAssembleSolve_MDOODZ_bridge(L, rhs_u, rhs_v, rhs_p,
                                                     i, j, block, rhs, sol);
    }

    const int Nx = L->Nx,  Nz = L->Nz;
    const int Ncx = L->Ncx;
    const double dx = L->dx, dz = L->dz;
    const double idx = 1.0 / dx, idz = 1.0 / dz;
    const double idx2 = idx * idx, idz2 = idz * idz;
    const double ixz  = idx * idz;
    const double *etan = L->etan;
    const double *etas = L->etas;

    // Process every active cell; boundary DOF rows are rewritten as
    // identity at the end so that Dirichlet faces are naturally held.
    if (i < 0 || i >= L->Ncx || j < 0 || j >= L->Ncz) return -1;
    if (!L->act_P[MdIdx_P(i, j, Ncx)]) return -1;
    // Skip if any local face is inactive (keeps the block well-posed).
    if (!L->act_Vx[MdIdx_Vx(i,     j, Nx )]) return -1;
    if (!L->act_Vx[MdIdx_Vx(i + 1, j, Nx )]) return -1;
    if (!L->act_Vz[MdIdx_Vz(i, j,     Ncx)]) return -1;
    if (!L->act_Vz[MdIdx_Vz(i, j + 1, Ncx)]) return -1;

    // Tag which of the five local DOFs is a Dirichlet (boundary) face.
    // Vx is boundary at i == 0 or i == Nx - 1; Vz at j == 0 or j == Nz - 1.
    const int bnd_VxL = (i     == 0)      || (i     == Nx - 1);
    const int bnd_VxR = (i + 1 == 0)      || (i + 1 == Nx - 1);
    const int bnd_VzB = (j     == 0)      || (j     == Nz - 1);
    const int bnd_VzT = (j + 1 == 0)      || (j + 1 == Nz - 1);

    // Viscosities used in the block
    const double etn_c     = etan[MdIdx_P(i,     j,     Ncx)]; // (i, j)
    const double etn_W     = etan[MdIdx_P(i - 1, j,     Ncx)]; // (i-1, j)   - left cell (for Vx_L's -σxx_L)
    const double etn_E     = etan[MdIdx_P(i + 1 < Ncx ? i + 1 : i, j, Ncx)]; // (i+1, j) if valid
    const double etn_S     = etan[MdIdx_P(i,     j - 1, Ncx)]; // (i, j-1)   - bottom cell
    const double etn_N     = etan[MdIdx_P(i,     j + 1 < L->Ncz ? j + 1 : j, Ncx)]; // (i, j+1)
    // etas at vertices around the cell
    const double es_SW     = etas[MdIdx_etas(i,     j,     Nx)]; // vertex (i, j)
    const double es_SE     = etas[MdIdx_etas(i + 1, j,     Nx)]; // vertex (i+1, j)
    const double es_NW     = etas[MdIdx_etas(i,     j + 1, Nx)]; // vertex (i, j+1)
    const double es_NE     = etas[MdIdx_etas(i + 1, j + 1, Nx)]; // vertex (i+1, j+1)

    // Zero block/rhs
    for (int r = 0; r < 5; ++r) { rhs[r] = 0.0; for (int c = 0; c < 5; ++c) block[r][c] = 0.0; }

    // ---- Row 0: mom-x at Vx_L = Vx(i, j) ----
    // Local coefficients (same derivation as in the header comment):
    //   Vx_L: -(2 etn_W + 2 etn_c)/dx² - (es_NW + es_SW)/dz²
    //   Vx_R:  2 etn_c / dx²
    //   Vz_B: -es_SW / (dx dz)
    //   Vz_T: +es_NW / (dx dz)
    //   P_C : -1 / dx
    block[0][0] = -(2.0 * etn_W + 2.0 * etn_c) * idx2 - (es_NW + es_SW) * idz2;
    block[0][1] =  2.0 * etn_c * idx2;
    block[0][2] = -es_SW * ixz;
    block[0][3] = +es_NW * ixz;
    block[0][4] = -idx;
    // RHS = rhs_u(i, j) - (contributions from neighbours outside the block):
    //   - outside Vx : Vx(i-1, j), Vx(i, j-1), Vx(i, j+1)
    //   - outside Vz : Vz(i-1, j), Vz(i-1, j+1)
    //   - outside P  : P(i-1, j)
    double neigh = 0.0;
    // Vx(i-1, j)
    neigh += (2.0 * etn_W * idx2) * gVx(Vx, i - 1, j, Nx);
    // Vx(i, j-1) and Vx(i, j+1) (may be at y-boundary — use value directly; Dirichlet ghost = 0 implicit)
    double Vx_below = (j - 1 >= 0      ) ? gVx(Vx, i, j - 1, Nx) : 0.0;
    double Vx_above = (j + 1 <  L->Ncz ) ? gVx(Vx, i, j + 1, Nx) : 0.0;
    neigh += (es_SW * idz2) * Vx_below;
    neigh += (es_NW * idz2) * Vx_above;
    // Vz(i-1, j) and Vz(i-1, j+1)
    double Vz_im1_j   = (i - 1 >= 0) ? gVz(Vz, i - 1, j,     Ncx) : 0.0;
    double Vz_im1_jp1 = (i - 1 >= 0) ? gVz(Vz, i - 1, j + 1, Ncx) : 0.0;
    neigh += (+es_SW * ixz) * Vz_im1_j;
    neigh += (-es_NW * ixz) * Vz_im1_jp1;
    // P(i-1, j) — from -∂P/∂x: -(P(i,j) - P(i-1,j))/dx → P(i-1,j) coef +1/dx
    neigh += (+idx) * gP(P, i - 1, j, Ncx);
    rhs[0] = rhs_u[MdIdx_Vx(i, j, Nx)] - neigh;

    // ---- Row 1: mom-x at Vx_R = Vx(i+1, j) ----
    //   Vx_L: +2 etn_c / dx²
    //   Vx_R: -(2 etn_c + 2 etn_E)/dx² - (es_NE + es_SE)/dz²
    //   Vz_B: +es_SE / (dx dz)
    //   Vz_T: -es_NE / (dx dz)
    //   P_C : +1 / dx
    // (NB: etn_E equals etan[i+1, j] in principle but we protect by clamping.)
    block[1][0] = 2.0 * etn_c * idx2;
    block[1][1] = -(2.0 * etn_c + 2.0 * etn_E) * idx2 - (es_NE + es_SE) * idz2;
    block[1][2] = +es_SE * ixz;
    block[1][3] = -es_NE * ixz;
    block[1][4] = +idx;
    // Outside contributions for Vx_R row:
    //   Vx(i+2, j) : coef 2 etn_E / dx²
    //   Vx(i+1, j-1), Vx(i+1, j+1) : coefs es_SE/dz², es_NE/dz²
    //   Vz(i+1, j), Vz(i+1, j+1) : coefs -es_SE/(dxdz), +es_NE/(dxdz)
    //   P(i+1, j) : coef -1/dx
    neigh = 0.0;
    double Vx_R_R = (i + 2 < Nx) ? gVx(Vx, i + 2, j, Nx) : 0.0;
    neigh += (2.0 * etn_E * idx2) * Vx_R_R;
    double Vx_Rbelow = (j - 1 >= 0      ) ? gVx(Vx, i + 1, j - 1, Nx) : 0.0;
    double Vx_Rabove = (j + 1 <  L->Ncz ) ? gVx(Vx, i + 1, j + 1, Nx) : 0.0;
    neigh += (es_SE * idz2) * Vx_Rbelow;
    neigh += (es_NE * idz2) * Vx_Rabove;
    double Vz_ip1_j   = (i + 1 < Ncx) ? gVz(Vz, i + 1, j,     Ncx) : 0.0;
    double Vz_ip1_jp1 = (i + 1 < Ncx) ? gVz(Vz, i + 1, j + 1, Ncx) : 0.0;
    neigh += (-es_SE * ixz) * Vz_ip1_j;
    neigh += (+es_NE * ixz) * Vz_ip1_jp1;
    double P_ip1 = (i + 1 < Ncx) ? gP(P, i + 1, j, Ncx) : 0.0;
    neigh += (-idx) * P_ip1;
    rhs[1] = rhs_u[MdIdx_Vx(i + 1, j, Nx)] - neigh;

    // ---- Row 2: mom-z at Vz_B = Vz(i, j) ----
    //   Vx_L: -es_SW / (dx dz)
    //   Vx_R: +es_SE / (dx dz)
    //   Vz_B: -(2 etn_c + 2 etn_S)/dz² - (es_SE + es_SW)/dx²
    //   Vz_T:  2 etn_c / dz²
    //   P_C : -1 / dz
    block[2][0] = -es_SW * ixz;
    block[2][1] = +es_SE * ixz;
    block[2][2] = -(2.0 * etn_c + 2.0 * etn_S) * idz2 - (es_SE + es_SW) * idx2;
    block[2][3] =  2.0 * etn_c * idz2;
    block[2][4] = -idz;
    // Outside contributions:
    //   Vz(i, j-1) : coef 2 etn_S / dz²
    //   Vz(i-1, j), Vz(i+1, j) : coefs es_SW/dx², es_SE/dx²
    //   Vx(i, j-1) : coef +es_SW / (dx dz)
    //   Vx(i+1, j-1) : coef -es_SE / (dx dz)
    //   P(i, j-1) : coef +1/dz
    neigh = 0.0;
    double Vz_B_B = (j - 1 >= 0) ? gVz(Vz, i, j - 1, Ncx) : 0.0;
    neigh += (2.0 * etn_S * idz2) * Vz_B_B;
    double Vz_W = (i - 1 >= 0)   ? gVz(Vz, i - 1, j, Ncx) : 0.0;
    double Vz_E = (i + 1 < Ncx)  ? gVz(Vz, i + 1, j, Ncx) : 0.0;
    neigh += (es_SW * idx2) * Vz_W;
    neigh += (es_SE * idx2) * Vz_E;
    double Vx_Bleft  = (j - 1 >= 0) ? gVx(Vx, i,     j - 1, Nx) : 0.0;
    double Vx_Bright = (j - 1 >= 0) ? gVx(Vx, i + 1, j - 1, Nx) : 0.0;
    neigh += (+es_SW * ixz) * Vx_Bleft;
    neigh += (-es_SE * ixz) * Vx_Bright;
    double P_jm1 = (j - 1 >= 0) ? gP(P, i, j - 1, Ncx) : 0.0;
    neigh += (+idz) * P_jm1;
    rhs[2] = rhs_v[MdIdx_Vz(i, j, Ncx)] - neigh;

    // ---- Row 3: mom-z at Vz_T = Vz(i, j+1) ----
    //   Vx_L: +es_NW / (dx dz)
    //   Vx_R: -es_NE / (dx dz)
    //   Vz_B:  2 etn_c / dz²
    //   Vz_T: -(2 etn_N + 2 etn_c)/dz² - (es_NE + es_NW)/dx²
    //   P_C : +1 / dz
    block[3][0] = +es_NW * ixz;
    block[3][1] = -es_NE * ixz;
    block[3][2] =  2.0 * etn_c * idz2;
    block[3][3] = -(2.0 * etn_N + 2.0 * etn_c) * idz2 - (es_NE + es_NW) * idx2;
    block[3][4] = +idz;
    // Outside contributions:
    //   Vz(i, j+2) : coef 2 etn_N / dz²
    //   Vz(i-1, j+1), Vz(i+1, j+1) : coefs es_NW/dx², es_NE/dx²
    //   Vx(i, j+1), Vx(i+1, j+1) : coefs -es_NW / (dx dz), +es_NE / (dx dz)
    //   P(i, j+1) : coef -1/dz
    neigh = 0.0;
    double Vz_T_T = (j + 2 < Nz) ? gVz(Vz, i, j + 2, Ncx) : 0.0;
    neigh += (2.0 * etn_N * idz2) * Vz_T_T;
    double Vz_TW = (i - 1 >= 0)  ? gVz(Vz, i - 1, j + 1, Ncx) : 0.0;
    double Vz_TE = (i + 1 < Ncx) ? gVz(Vz, i + 1, j + 1, Ncx) : 0.0;
    neigh += (es_NW * idx2) * Vz_TW;
    neigh += (es_NE * idx2) * Vz_TE;
    double Vx_Tleft  = (j + 1 < L->Ncz) ? gVx(Vx, i,     j + 1, Nx) : 0.0;
    double Vx_Tright = (j + 1 < L->Ncz) ? gVx(Vx, i + 1, j + 1, Nx) : 0.0;
    neigh += (-es_NW * ixz) * Vx_Tleft;
    neigh += (+es_NE * ixz) * Vx_Tright;
    double P_jp1 = (j + 1 < L->Ncz) ? gP(P, i, j + 1, Ncx) : 0.0;
    neigh += (-idz) * P_jp1;
    rhs[3] = rhs_v[MdIdx_Vz(i, j + 1, Ncx)] - neigh;

    // ---- Row 4: continuity at P_C ----
    //   Vx_L: -1/dx, Vx_R: +1/dx, Vz_B: -1/dz, Vz_T: +1/dz, P_C: 0
    block[4][0] = -idx;
    block[4][1] = +idx;
    block[4][2] = -idz;
    block[4][3] = +idz;
    block[4][4] =  0.0;
    rhs[4] = rhs_p[MdIdx_P(i, j, Ncx)]; // continuity has no "outside" DOFs in its stencil

    // For any DOF that sits on a Dirichlet (boundary) face, replace its
    // equation with an identity row: row_r = e_r, rhs_r = rhs_* at that
    // boundary index. This matches StokesApplyA (identity rows on the
    // outer ring), so the block agrees with the global operator and the
    // sweep preserves boundary values instead of drifting from them.
    // Additionally, the identity row decouples the boundary DOF from the
    // other four unknowns in the 5x5 so partial-pivoting always succeeds.
    const int bnd_row[5] = { bnd_VxL, bnd_VxR, bnd_VzB, bnd_VzT, 0 };
    const double bnd_rhs[5] = {
        rhs_u[MdIdx_Vx(i,     j, Nx )],
        rhs_u[MdIdx_Vx(i + 1, j, Nx )],
        rhs_v[MdIdx_Vz(i, j,     Ncx)],
        rhs_v[MdIdx_Vz(i, j + 1, Ncx)],
        0.0,
    };
    for (int r = 0; r < 5; ++r) {
        if (!bnd_row[r]) continue;
        for (int c = 0; c < 5; ++c) block[r][c] = (c == r) ? 1.0 : 0.0;
        rhs[r] = bnd_rhs[r];
        // Remove this boundary DOF's contribution from the other rows so
        // that their unknowns really are the interior four. We add back
        // block[r'][r] * rhs[r] to rhs[r'], then zero the column.
        for (int rp = 0; rp < 5; ++rp) {
            if (rp == r) continue;
            rhs[rp] -= block[rp][r] * rhs[r];
            block[rp][r] = 0.0;
        }
    }

    // Zero pressure block diagonal is a saddle-point — solve5x5 with
    // partial pivoting handles it fine as long as the velocity block is
    // non-singular.
    return solve5x5(block, rhs, sol);
}

// Allocate and populate the MDOODZ-padded scratch buffers on a
// `MultigridLevel`. Called lazily from `VankaSweep` when the level has
// `use_mdoodz_matvec_vanka` set. Returns 0 on success, non-zero if allocation
// failed. Populating is done here once; `MdoodzPadFromLevel` repacks
// per-sweep, and `MdoodzPadUpdateCell` keeps the buffers in lock-step
// with per-cell under-relaxed updates so Gauss-Seidel freshness is
// preserved.
static int EnsureMdoodzPadScratch(MultigridLevel *L) {
    if (L->md_u_pad != NULL && L->md_v_pad != NULL && L->md_p_pad != NULL) return 0;
    const int Nx = L->Nx, Nz = L->Nz, Ncx = L->Ncx, Ncz = L->Ncz;
    const int Nzvx = Nz + 1, Nxvz = Nx + 1;
    const size_t nu_pad = (size_t)Nx   * (size_t)Nzvx;
    const size_t nv_pad = (size_t)Nxvz * (size_t)Nz;
    const size_t np_pad = (size_t)Ncx  * (size_t)Ncz;
    double *u_pad = (double *)calloc(nu_pad, sizeof(double));
    double *v_pad = (double *)calloc(nv_pad, sizeof(double));
    double *p_pad = (double *)calloc(np_pad, sizeof(double));
    if (u_pad == NULL || v_pad == NULL || p_pad == NULL) {
        free(u_pad); free(v_pad); free(p_pad);
        LOG_ERR("GMG: failed to allocate Vanka bridge scratch (Nx=%d Nz=%d)", Nx, Nz);
        return -1;
    }
    L->md_u_pad = u_pad; L->md_v_pad = v_pad; L->md_p_pad = p_pad;
    return 0;
}

// Repack L->Vx/Vz/P (strict interior) → L->md_{u,v,p}_pad (MDOODZ padded).
// The identity rows written by `BuildLevelRHSFromCompressed` guarantee
// the ghost rows/cols see the Dirichlet value in L->Vx/Vz/P, so a plain
// copy is sufficient and the bridge block-builder sees the correct
// boundary values when it walks MDOODZ-padded neighbour indices.
static void MdoodzPadFromLevel(const MultigridLevel *L) {
    const int Nx = L->Nx, Nz = L->Nz, Ncx = L->Ncx, Ncz = L->Ncz;
    const int Nxvz = Nx + 1;
    for (int j = 0; j < Ncz; ++j) {
        const int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            L->md_u_pad[(size_t)i + (size_t)l * Nx] = L->Vx[MdIdx_Vx(i, j, Nx)];
        }
    }
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const int k = i + 1;
            L->md_v_pad[(size_t)k + (size_t)j * Nxvz] = L->Vz[MdIdx_Vz(i, j, Ncx)];
        }
    }
    memcpy(L->md_p_pad, L->P, (size_t)Ncx * Ncz * sizeof(double));
}

// Propagate the 5 DOFs updated by a Vanka cell solve into the MDOODZ
// padded scratch. Cheap (5 writes) and keeps the per-sweep repack's
// work consistent with Gauss-Seidel colouring.
static inline void MdoodzPadUpdateCell(MultigridLevel *L, int i, int j) {
    const int Nx = L->Nx, Ncx = L->Ncx, Nxvz = Nx + 1;
    L->md_u_pad[(size_t)i     + (size_t)(j + 1) * Nx]   = L->Vx[MdIdx_Vx(i,     j,     Nx )];
    L->md_u_pad[(size_t)(i+1) + (size_t)(j + 1) * Nx]   = L->Vx[MdIdx_Vx(i + 1, j,     Nx )];
    L->md_v_pad[(size_t)(i+1) + (size_t) j      * Nxvz] = L->Vz[MdIdx_Vz(i,     j,     Ncx)];
    L->md_v_pad[(size_t)(i+1) + (size_t)(j + 1) * Nxvz] = L->Vz[MdIdx_Vz(i,     j + 1, Ncx)];
    L->md_p_pad[MdIdx_P(i, j, Ncx)]                     = L->P [MdIdx_P (i,     j,     Ncx)];
}

void VankaSweep(MultigridLevel *L, int nu) {
    const int Nx = L->Nx;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    // Under-relaxation for the saddle-point block. Without damping, the
    // overlapping block-GS iteration can diverge on a pure Stokes system.
    // ω ≈ 0.6 is a robust default (Vanka 1986, John & Tobiska 2000).
    const double omega = 0.6;

    // MDOODZ bridge scratch for the fine-level Vanka path (task 15.8).
    // Allocated once per level and repacked at the start of each sweep.
    const int use_bridge = (L->use_mdoodz_matvec_vanka != 0);
    if (use_bridge) {
        if (EnsureMdoodzPadScratch(L) != 0) return;
    }

    for (int sweep = 0; sweep < nu; ++sweep) {
        if (use_bridge) MdoodzPadFromLevel(L);
        for (int colour = 0; colour < 2; ++colour) {
            // Parallelisable over cells of a single colour.
            for (int j = 0; j < Ncz; ++j) {
                for (int i = 0; i < Ncx; ++i) {
                    if (((i + j) & 1) != colour) continue;
                    double block[5][5], rhs[5], sol[5];
                    int rc = VankaBlockAssembleSolve(L, L->Vx, L->Vz, L->P,
                                                    L->rhs_u, L->rhs_v, L->rhs_p,
                                                    i, j, block, rhs, sol);
                    if (rc != 0) continue;
                    // sol holds the new candidate values (not corrections),
                    // because we built the RHS as f - A_outside * x_outside.
                    // Apply under-relaxation: new = (1-omega)*old + omega*sol.
                    L->Vx[MdIdx_Vx(i,     j, Nx)]  = (1.0 - omega) * L->Vx[MdIdx_Vx(i,     j, Nx)]  + omega * sol[0];
                    L->Vx[MdIdx_Vx(i + 1, j, Nx)]  = (1.0 - omega) * L->Vx[MdIdx_Vx(i + 1, j, Nx)]  + omega * sol[1];
                    L->Vz[MdIdx_Vz(i, j,     Ncx)] = (1.0 - omega) * L->Vz[MdIdx_Vz(i, j,     Ncx)] + omega * sol[2];
                    L->Vz[MdIdx_Vz(i, j + 1, Ncx)] = (1.0 - omega) * L->Vz[MdIdx_Vz(i, j + 1, Ncx)] + omega * sol[3];
                    L->P [MdIdx_P (i,     j, Ncx)] = (1.0 - omega) * L->P [MdIdx_P (i,     j, Ncx)] + omega * sol[4];
                    if (use_bridge) MdoodzPadUpdateCell(L, i, j);
                }
            }
        }
    }
}

// ---------- Coarsest-level solve (Section 6) ----------------------------
//
// We assemble the coarsest-level stencil into a CSR sparse matrix
// (obtained by evaluating StokesApplyA against local unit vectors), and
// factor it with UMFPACK. The Stokes saddle-point is indefinite, so
// CHOLMOD is not applicable; UMFPACK's general unsymmetric LU handles it
// cleanly. The factor is cached on the hierarchy and re-used by every
// V-cycle until invalidated (e.g. when viscosities change structurally).

// Build the coarse DOF numbering. "DOF" == one scalar unknown in the
// assembled system. For simplicity we include every Vx/Vz/P entry
// (both interior and outer-ring boundary); boundary and inactive rows
// get collapsed to identity rows during the stencil extraction below.
// This keeps indexing trivial (DOF index is a linear function of array
// position) at the modest cost of a few dozen identity rows per level —
// dwarfed by the O(n^2) stencil extraction itself.
static void build_coarse_dof_numbering(MultigridHierarchy *H) {
    MultigridLevel *C = &H->levels[H->n_levels - 1];
    const int Nx = C->Nx, Nz = C->Nz;
    const int Ncx = C->Ncx, Ncz = C->Ncz;
    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    const size_t nP  = (size_t)Ncx * Ncz;
    const int    n   = (int)(nVx + nVz + nP);

    H->coarse_n   = n;
    free(H->coarse_eqn_Vx);
    free(H->coarse_eqn_Vz);
    free(H->coarse_eqn_P);
    H->coarse_eqn_Vx = (int *)malloc(nVx * sizeof(int));
    H->coarse_eqn_Vz = (int *)malloc(nVz * sizeof(int));
    H->coarse_eqn_P  = (int *)malloc(nP  * sizeof(int));

    int eq = 0;
    for (size_t k = 0; k < nVx; ++k) H->coarse_eqn_Vx[k] = eq++;
    for (size_t k = 0; k < nVz; ++k) H->coarse_eqn_Vz[k] = eq++;
    for (size_t k = 0; k < nP;  ++k) H->coarse_eqn_P [k] = eq++;
}

// Extract row `eq` of the coarsest-level operator into `row_buf`, where
// row_buf is a length-n dense scratch that has already been zeroed.
// Returns 0 if the row is the DOF's stencil row, non-zero if the row is
// an identity (boundary or inactive) row. When identity, only row_buf[eq]
// is set.
//
// We extract a row by a double-matvec trick: given basis vector e_eq we
// compute y = A e_eq. The k-th component of y is A[k, eq]; the eq-th
// component of y is A[eq, eq]. To get the eq-th ROW we'd need A^T; since
// our operator is symmetric in the velocity block and antisymmetric
// between B and C, we instead build the matrix column-by-column (easier)
// and convert into CSR form at the end.
//
// (UMFPACK interprets a CSR of A as a CSC of A^T, and UMFPACK_Aat tells
// it to solve A^T x = b — i.e. A x = b in our CSR interpretation.)
static void extract_coarse_column(MultigridHierarchy *H, int eq, double *col_buf) {
    MultigridLevel *C = &H->levels[H->n_levels - 1];
    const int Nx = C->Nx, Nz = C->Nz;
    const int Ncx = C->Ncx, Ncz = C->Ncz;
    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    const size_t nP  = (size_t)Ncx * Ncz;

    // Use res_u/res_v/res_p as unit-vector scratch; cor_* as output.
    double *ux = C->res_u;
    double *uz = C->res_v;
    double *up = C->res_p;
    double *yx = C->cor_u;
    double *yz = C->cor_v;
    double *yp = C->cor_p;

    memset(ux, 0, nVx * sizeof(double));
    memset(uz, 0, nVz * sizeof(double));
    memset(up, 0, nP  * sizeof(double));

    if ((size_t)eq < nVx) {
        ux[eq] = 1.0;
    } else if ((size_t)eq < nVx + nVz) {
        uz[eq - nVx] = 1.0;
    } else {
        up[eq - nVx - nVz] = 1.0;
    }

    StokesApplyA(C, ux, uz, up, yx, yz, yp);

    // Pack into col_buf: identity-row patching is deferred to the caller
    // when it assembles the row-wise matrix, because a given column can
    // legitimately hold off-diagonal stencil entries for boundary rows
    // that WE want to overwrite with identity.
    for (size_t k = 0; k < nVx; ++k) col_buf[k]             = yx[k];
    for (size_t k = 0; k < nVz; ++k) col_buf[nVx + k]       = yz[k];
    for (size_t k = 0; k < nP;  ++k) col_buf[nVx + nVz + k] = yp[k];
}

// Identify boundary / inactive "identity" rows on the coarsest level.
// Returns 1 iff the row equation `eq` is an identity row in StokesApplyA.
static int coarse_row_is_identity(const MultigridLevel *C, int eq) {
    const int Nx = C->Nx, Nz = C->Nz;
    const int Ncx = C->Ncx, Ncz = C->Ncz;
    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    if ((size_t)eq < nVx) {
        // Vx DOF at (i, j)
        int j = eq / Nx;
        int i = eq - j * Nx;
        if (i == 0 || i == Nx - 1) return 1;
        return !C->act_Vx[(size_t)eq];
    }
    int k = eq - (int)nVx;
    if ((size_t)k < nVz) {
        int j = k / Ncx;
        int i = k - j * Ncx;
        if (j == 0 || j == Nz - 1) return 1;
        return !C->act_Vz[(size_t)k];
    }
    int p = k - (int)nVz;
    return !C->act_P[(size_t)p];
}

int CoarseAssembleAndFactor(MultigridHierarchy *H) {
    if (H == NULL || H->n_levels < 1) return -1;
    if (H->coarse_factor_valid) return 0;

    MultigridLevel *C = &H->levels[H->n_levels - 1];
    build_coarse_dof_numbering(H);
    const int n = H->coarse_n;

    // Dense workspace (n up to a few thousand DOFs at the coarsest level).
    double *col = (double *)calloc((size_t)n, sizeof(double));
    if (col == NULL) return -1;

    // First pass: count nnz per row.
    int *row_nnz = (int *)calloc((size_t)n, sizeof(int));
    if (row_nnz == NULL) { free(col); return -1; }

    for (int eq_col = 0; eq_col < n; ++eq_col) {
        extract_coarse_column(H, eq_col, col);
        for (int r = 0; r < n; ++r) {
            double v = col[r];
            if (v == 0.0) continue;
            // Skip writing a stencil entry if the destination row is
            // going to be an identity row.
            if (coarse_row_is_identity(C, r) && r != eq_col) continue;
            ++row_nnz[r];
        }
        // Ensure every identity row has a diagonal entry (might miss
        // otherwise if the unit vector happens to land on a zero).
    }
    // Guarantee at least one entry per identity row (the diagonal).
    for (int r = 0; r < n; ++r) {
        if (coarse_row_is_identity(C, r) && row_nnz[r] == 0) row_nnz[r] = 1;
    }

    // Build CSR row pointers.
    free(H->coarse_Ap);
    free(H->coarse_Aj);
    free(H->coarse_Ax);
    H->coarse_Ap = (int *)malloc((size_t)(n + 1) * sizeof(int));
    if (H->coarse_Ap == NULL) { free(col); free(row_nnz); return -1; }
    H->coarse_Ap[0] = 0;
    for (int r = 0; r < n; ++r) H->coarse_Ap[r + 1] = H->coarse_Ap[r] + row_nnz[r];
    int nnz = H->coarse_Ap[n];
    H->coarse_nnz = nnz;
    H->coarse_Aj  = (int *)malloc((size_t)nnz * sizeof(int));
    H->coarse_Ax  = (double *)malloc((size_t)nnz * sizeof(double));
    if (H->coarse_Aj == NULL || H->coarse_Ax == NULL) {
        free(col); free(row_nnz);
        free(H->coarse_Ap); free(H->coarse_Aj); free(H->coarse_Ax);
        H->coarse_Ap = NULL; H->coarse_Aj = NULL; H->coarse_Ax = NULL;
        return -1;
    }
    // Running insertion pointer per row.
    int *row_pos = (int *)calloc((size_t)n, sizeof(int));
    if (row_pos == NULL) { free(col); free(row_nnz); return -1; }

    // Second pass: fill entries column by column, respecting identity
    // rewriting.
    for (int eq_col = 0; eq_col < n; ++eq_col) {
        extract_coarse_column(H, eq_col, col);
        for (int r = 0; r < n; ++r) {
            double v = col[r];
            int is_id = coarse_row_is_identity(C, r);
            if (is_id && r != eq_col) continue;
            if (is_id && r == eq_col) v = 1.0;  // identity-row diagonal
            if (!is_id && v == 0.0)   continue;
            int idx = H->coarse_Ap[r] + row_pos[r]++;
            H->coarse_Aj[idx] = eq_col;
            H->coarse_Ax[idx] = v;
        }
    }
    // Patch identity rows that got no stencil entry from the columns
    // (happens when the row's DOF never couples back into any column's
    // stencil — e.g. an inactive Vx whose neighbours are all inactive
    // too). Write an identity diagonal.
    for (int r = 0; r < n; ++r) {
        if (!coarse_row_is_identity(C, r)) continue;
        if (row_pos[r] != 0) continue;
        if (row_nnz[r] < 1)  continue;
        int idx = H->coarse_Ap[r];
        H->coarse_Aj[idx] = r;
        H->coarse_Ax[idx] = 1.0;
        row_pos[r] = 1;
    }
    free(col); free(row_nnz); free(row_pos);

    // Factor.
    void *Symbolic = NULL;
    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    umfpack_di_defaults(Control);
    Control[UMFPACK_STRATEGY] = UMFPACK_STRATEGY_UNSYMMETRIC;
    Control[UMFPACK_IRSTEP]   = 0;

    int status = umfpack_di_symbolic(n, n, H->coarse_Ap, H->coarse_Aj, H->coarse_Ax,
                                     &Symbolic, Control, Info);
    if (status < 0) {
        LOG_ERR("GMG coarse: umfpack_di_symbolic failed (status=%d)", status);
        return -1;
    }
    void *Numeric = NULL;
    status = umfpack_di_numeric(H->coarse_Ap, H->coarse_Aj, H->coarse_Ax,
                                Symbolic, &Numeric, Control, Info);
    umfpack_di_free_symbolic(&Symbolic);
    if (status < 0) {
        LOG_ERR("GMG coarse: umfpack_di_numeric failed (status=%d)", status);
        return -1;
    }
    if (H->coarse_umf_numeric != NULL) {
        void *old = H->coarse_umf_numeric;
        umfpack_di_free_numeric(&old);
    }
    H->coarse_umf_numeric  = Numeric;
    H->coarse_factor_valid = 1;
    return 0;
}

int CoarseSolve(MultigridHierarchy *H) {
    if (H == NULL || H->n_levels < 1) return -1;
    MultigridLevel *C = &H->levels[H->n_levels - 1];

    // If UMFPACK factor is not ready, (re-)build it.
    if (!H->coarse_factor_valid) {
        if (CoarseAssembleAndFactor(H) != 0) {
            // Fall back to the iterative pseudo-solve path.
            LOG_WARN("GMG coarse: UMFPACK assembly/factor failed; using iterative Vanka.");
            const int    max_iters = 5000;
            const double rtol      = 1e-10;
            StokesResidual(C);
            double r0 = StokesResidualNorm(C);
            if (r0 < 1e-300) return 0;
            double r  = r0;
            int    it = 0, stall_count = 0;
            while (it < max_iters && r > rtol * r0) {
                VankaSweep(C, 1);
                StokesResidual(C);
                double rn = StokesResidualNorm(C);
                if (rn > r * 0.9999) ++stall_count; else stall_count = 0;
                if (stall_count > 100 && it > 200) break;
                r = rn;
                ++it;
            }
            return 0;
        }
    }

    // Pack RHS into a flat vector in (Vx | Vz | P) layout.
    const int n = H->coarse_n;
    const int Nx = C->Nx, Nz = C->Nz;
    const int Ncx = C->Ncx, Ncz = C->Ncz;
    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    const size_t nP  = (size_t)Ncx * Ncz;

    double *b = (double *)malloc((size_t)n * sizeof(double));
    double *x = (double *)malloc((size_t)n * sizeof(double));
    if (b == NULL || x == NULL) { free(b); free(x); return -1; }
    // Copy RHS; for identity rows use the current boundary/fixed value
    // in Vx/Vz/P so that iteration after solve preserves them.
    for (size_t k = 0; k < nVx; ++k) {
        int eq = (int)k;
        b[eq] = coarse_row_is_identity(C, eq) ? C->Vx[k] : C->rhs_u[k];
    }
    for (size_t k = 0; k < nVz; ++k) {
        int eq = (int)(nVx + k);
        b[eq] = coarse_row_is_identity(C, eq) ? C->Vz[k] : C->rhs_v[k];
    }
    for (size_t k = 0; k < nP;  ++k) {
        int eq = (int)(nVx + nVz + k);
        b[eq] = coarse_row_is_identity(C, eq) ? C->P[k]  : C->rhs_p[k];
    }

    double Control[UMFPACK_CONTROL], Info[UMFPACK_INFO];
    umfpack_di_defaults(Control);
    Control[UMFPACK_IRSTEP] = 0;

    // Our matrix is stored in CSR; UMFPACK wants CSC. Solving UMFPACK_Aat
    // with a CSR layout is algebraically the same as A x = b in CSR.
    int status = umfpack_di_solve(UMFPACK_Aat, H->coarse_Ap, H->coarse_Aj, H->coarse_Ax,
                                  x, b, H->coarse_umf_numeric, Control, Info);
    if (status < 0) {
        LOG_WARN("GMG coarse: umfpack_di_solve failed (status=%d); retaining previous iterate.",
                 status);
        free(b); free(x);
        return -1;
    }

    // Unpack solution.
    memcpy(C->Vx, x,                nVx * sizeof(double));
    memcpy(C->Vz, x + nVx,          nVz * sizeof(double));
    memcpy(C->P , x + nVx + nVz,    nP  * sizeof(double));
    free(b); free(x);
    return 0;
}

// ---------- V-cycle driver (Section 7) ----------------------------------

// Forward decl for the pressure-mean projector (defined earlier in this file).
static void project_pressure_mean(double *Pfield, const char *mask, int Ncx, int Ncz);

// ---- V-cycle dump (gmg_dump_vcycle = 1) --------------------------------
//
// Design D5/D9 of add-gmg-stokes-defence: on the first V-cycle of the first
// FGMRES iteration of the first Picard iteration of the first time step,
// emit per-level HDF5 snapshots of Vx/Vz/P/res_u/res_v/res_p immediately
// before and after each V-cycle operator (pre-smooth, restrict, coarse
// solve, prolongate, post-smooth). After that one cycle, the static
// `fired` flag pins the dumper shut for the rest of the process lifetime
// so a long run does not flood disk.
//
// Filenames carry a monotonically increasing step counter so the post-
// processor (`benchmarks/vcycle_vis/make_gif.py`) can reconstruct the
// traversal order deterministically even on deep hierarchies:
//
//   ${writer_subfolder}/vcycle_dump/level_{N}_{step}_{phase}_{seq}.h5
//
// where `seq` is zero-padded 3-digit order.

static struct {
    int   active;            // dumper is recording this V-cycle
    int   fired;              // one-shot: set after the first recorded cycle
    char  output_dir[1024];
    int   step_counter;
} g_vcycle_dumper = { 0, 0, {0}, 0 };

static void gmg_dump_mkdir_p(const char *path) {
    // Cheap mkdir -p: walk the path, mkdir each component. Silently
    // ignores EEXIST. Good enough for ${writer_subfolder}/vcycle_dump/.
    if (!path || !*path) return;
    char buf[1024];
    size_t n = strlen(path);
    if (n >= sizeof(buf)) return;
    memcpy(buf, path, n + 1);
    for (size_t i = 1; i < n; ++i) {
        if (buf[i] == '/') {
            buf[i] = 0;
#ifdef _WIN32
            mkdir(buf);
#else
            mkdir(buf, 0775);
#endif
            buf[i] = '/';
        }
    }
#ifdef _WIN32
    mkdir(buf);
#else
    mkdir(buf, 0775);
#endif
}

// Arm the dumper if the caller asked for it and it has not fired yet.
// Called at SolveStokesGMG entry; returns 1 if the current invocation
// will record, 0 otherwise.
static int gmg_vcycle_dump_arm(const char *writer_subfolder,
                               int gmg_dump_vcycle) {
    if (gmg_dump_vcycle != 1 || g_vcycle_dumper.fired) {
        g_vcycle_dumper.active = 0;
        return 0;
    }
    const char *base = (writer_subfolder && *writer_subfolder) ? writer_subfolder : "./";
    // Trim trailing slash from base for clean concatenation.
    char trimmed[900];
    size_t nb = strlen(base);
    if (nb >= sizeof(trimmed)) {
        LOG_WARN("GMG dump: writer_subfolder too long (%zu chars); dump disabled", nb);
        return 0;
    }
    memcpy(trimmed, base, nb + 1);
    if (nb > 0 && trimmed[nb - 1] == '/') trimmed[nb - 1] = 0;
    int written = snprintf(g_vcycle_dumper.output_dir,
                           sizeof(g_vcycle_dumper.output_dir),
                           "%s/vcycle_dump", trimmed);
    if (written < 0 || (size_t)written >= sizeof(g_vcycle_dumper.output_dir)) {
        LOG_WARN("GMG dump: output path overflow; dump disabled");
        return 0;
    }
    gmg_dump_mkdir_p(g_vcycle_dumper.output_dir);
    g_vcycle_dumper.step_counter = 0;
    g_vcycle_dumper.active       = 1;
    LOG_INFO("GMG V-cycle dump armed: output directory = %s",
             g_vcycle_dumper.output_dir);
    return 1;
}

// Disarm + record that a dump has fired. Called after the first V-cycle
// on the fine level returns.
static void gmg_vcycle_dump_finish(void) {
    if (!g_vcycle_dumper.active) return;
    LOG_INFO("GMG V-cycle dump complete: %d snapshots written to %s",
             g_vcycle_dumper.step_counter, g_vcycle_dumper.output_dir);
    g_vcycle_dumper.active = 0;
    g_vcycle_dumper.fired  = 1;
}

// Reset the one-shot guard. Exposed for tests that need to exercise the
// dump more than once per process (unit tests live outside production
// paths). Production never calls this.
void GmgVcycleDumpReset(void) {
    g_vcycle_dumper.active = 0;
    g_vcycle_dumper.fired  = 0;
    g_vcycle_dumper.step_counter = 0;
    g_vcycle_dumper.output_dir[0] = 0;
}

// Test-only: arm the dumper with an explicit output directory. Production
// always arms via gmg_vcycle_dump_arm(writer_subfolder,gmg_dump_vcycle).
int GmgVcycleDumpArmForTest(const char *output_dir) {
    if (g_vcycle_dumper.fired) {
        g_vcycle_dumper.active = 0;
        return 0;
    }
    if (output_dir == NULL || *output_dir == 0) return 0;
    size_t n = strlen(output_dir);
    if (n >= sizeof(g_vcycle_dumper.output_dir)) return 0;
    memcpy(g_vcycle_dumper.output_dir, output_dir, n + 1);
    gmg_dump_mkdir_p(g_vcycle_dumper.output_dir);
    g_vcycle_dumper.step_counter = 0;
    g_vcycle_dumper.active = 1;
    return 1;
}

// Test-only: close the current recording and latch the one-shot guard.
// Returns the snapshot count written during this recording.
int GmgVcycleDumpFinishForTest(void) {
    int count = g_vcycle_dumper.step_counter;
    if (!g_vcycle_dumper.active) return 0;
    g_vcycle_dumper.active = 0;
    g_vcycle_dumper.fired  = 1;
    return count;
}

// Emit one snapshot. `step_name` is one of
//   "pre_smooth" / "restrict" / "coarse_solve" / "prolongate" / "post_smooth"
// and `phase` is "pre" or "post". Snapshot contains:
//   /Level/{Nx, Nz, Ncx, Ncz, dx, dz, level}
//   /Fields/{Vx, Vz, P, res_u, res_v, res_p}
//   /Meta/{step, phase, seq}
static void gmg_vcycle_dump_snapshot(MultigridLevel *L,
                                     const char *step_name,
                                     const char *phase) {
    if (!g_vcycle_dumper.active || L == NULL) return;

    // Refresh residual on this level so Fields/res_* matches Fields/V*.
    // StokesResidual = rhs - A * iterate; always well-defined at every
    // point in the V-cycle.
    StokesResidual(L);

    const int seq = g_vcycle_dumper.step_counter++;
    char filename[1100];
    int written = snprintf(filename, sizeof(filename),
                           "%s/level_%d_%s_%s_%03d.h5",
                           g_vcycle_dumper.output_dir, L->level,
                           step_name, phase, seq);
    if (written < 0 || (size_t)written >= sizeof(filename)) {
        LOG_WARN("GMG dump: filename overflow at seq=%d", seq);
        return;
    }

    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;

    CreateOutputHDF5(filename);
    AddGroupToHDF5(filename, "Level");
    AddGroupToHDF5(filename, "Fields");
    AddGroupToHDF5(filename, "Meta");

    int    level_dims[4] = { L->Nx, L->Nz, L->Ncx, L->Ncz };
    double level_spacing[2] = { L->dx, L->dz };
    int    level_index = L->level;
    AddFieldToGroup(filename, "Level", "dims",     'i', 4, level_dims,    1);
    AddFieldToGroup(filename, "Level", "spacing",  'd', 2, level_spacing, 1);
    AddFieldToGroup(filename, "Level", "level",    'i', 1, &level_index,  1);

    AddFieldToGroup(filename, "Fields", "Vx",      'd', (int)nVx, L->Vx,    1);
    AddFieldToGroup(filename, "Fields", "Vz",      'd', (int)nVz, L->Vz,    1);
    AddFieldToGroup(filename, "Fields", "P",       'd', (int)nP,  L->P,     1);
    AddFieldToGroup(filename, "Fields", "res_u",   'd', (int)nVx, L->res_u, 1);
    AddFieldToGroup(filename, "Fields", "res_v",   'd', (int)nVz, L->res_v, 1);
    AddFieldToGroup(filename, "Fields", "res_p",   'd', (int)nP,  L->res_p, 1);

    AddFieldToGroup(filename, "Meta", "step",  'c',
                    (int)(strlen(step_name) + 1), (void *)step_name, 1);
    AddFieldToGroup(filename, "Meta", "phase", 'c',
                    (int)(strlen(phase) + 1), (void *)phase, 1);
    AddFieldToGroup(filename, "Meta", "seq",   'i', 1, (void *)&seq, 1);
}

void Vcycle(MultigridHierarchy *H, int k) {
    if (H == NULL || k < 0 || k >= H->n_levels) return;
    MultigridLevel *L = &H->levels[k];

    if (k == H->n_levels - 1) {
        // Coarsest: solve. Project pressure null-space out of the RHS to
        // ensure solvability.
        gmg_vcycle_dump_snapshot(L, "coarse_solve", "pre");
        project_pressure_mean(L->rhs_p, L->act_P, L->Ncx, L->Ncz);
        CoarseSolve(H);
        project_pressure_mean(L->P, L->act_P, L->Ncx, L->Ncz);
        gmg_vcycle_dump_snapshot(L, "coarse_solve", "post");
        return;
    }

    // Pre-smooth, compute residual, restrict to coarser grid.
    gmg_vcycle_dump_snapshot(L, "pre_smooth", "pre");
    VankaSweep(L, H->nu_pre);
    gmg_vcycle_dump_snapshot(L, "pre_smooth", "post");

    StokesResidual(L);

    MultigridLevel *C = &H->levels[k + 1];
    gmg_vcycle_dump_snapshot(L, "restrict", "pre");
    RestrictVelocity(L, C, L->res_u, C->rhs_u, 0);
    RestrictVelocity(L, C, L->res_v, C->rhs_v, 1);
    RestrictPressure(L, C, L->res_p, C->rhs_p);

    const size_t nVx_c = (size_t)C->Nx  * C->Ncz;
    const size_t nVz_c = (size_t)C->Ncx * C->Nz;
    const size_t nP_c  = (size_t)C->Ncx * C->Ncz;
    memset(C->Vx, 0, nVx_c * sizeof(double));
    memset(C->Vz, 0, nVz_c * sizeof(double));
    memset(C->P,  0, nP_c  * sizeof(double));
    gmg_vcycle_dump_snapshot(C, "restrict", "post");

    Vcycle(H, k + 1);

    gmg_vcycle_dump_snapshot(L, "prolongate", "pre");
    ProlongateVelocityAdd(C, L, C->Vx, L->Vx, 0);
    ProlongateVelocityAdd(C, L, C->Vz, L->Vz, 1);
    ProlongatePressureAdd(C, L, C->P,  L->P );
    gmg_vcycle_dump_snapshot(L, "prolongate", "post");

    // Post-smooth and project pressure null-space at every level so
    // intermediate states stay in the gauge used by StokesApplyA.
    gmg_vcycle_dump_snapshot(L, "post_smooth", "pre");
    VankaSweep(L, H->nu_post);
    project_pressure_mean(L->P, L->act_P, L->Ncx, L->Ncz);
    gmg_vcycle_dump_snapshot(L, "post_smooth", "post");
}

// ---------- Flatten / scatter helpers (for FGMRES) ----------------------
//
// FGMRES operates on a single concatenated vector [Vx | Vz | P]. We pack
// and unpack between this layout and the level's arrays as needed.

static void level_pack(const MultigridLevel *L, const double *Vx, const double *Vz,
                       const double *P, double *flat) {
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    memcpy(flat,                 Vx, nVx * sizeof(double));
    memcpy(flat + nVx,           Vz, nVz * sizeof(double));
    memcpy(flat + nVx + nVz,     P,  nP  * sizeof(double));
}

static void level_unpack(const MultigridLevel *L, const double *flat,
                         double *Vx, double *Vz, double *P) {
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    memcpy(Vx, flat,                 nVx * sizeof(double));
    memcpy(Vz, flat + nVx,           nVz * sizeof(double));
    memcpy(P,  flat + nVx + nVz,     nP  * sizeof(double));
}

static size_t level_total_dof(const MultigridLevel *L) {
    return (size_t)L->Nx * L->Ncz + (size_t)L->Ncx * L->Nz + (size_t)L->Ncx * L->Ncz;
}

// dot product
static double vec_dot(const double *a, const double *b, size_t n) {
    double s = 0.0;
    for (size_t k = 0; k < n; ++k) s += a[k] * b[k];
    return s;
}
static double vec_norm(const double *a, size_t n) { return sqrt(vec_dot(a, a, n)); }

// ---------- FGMRES (Section 8) ------------------------------------------

// Apply A to a flattened vector x, writing into y.
static void apply_A_flat(MultigridLevel *L, const double *x_flat, double *y_flat) {
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    // We need scratch with same shape as Vx/Vz/P; use res_* as temp.
    level_unpack(L, x_flat, L->res_u, L->res_v, L->res_p);
    StokesApplyA(L, L->res_u, L->res_v, L->res_p, L->cor_u, L->cor_v, L->cor_p);
    memcpy(y_flat,              L->cor_u, nVx * sizeof(double));
    memcpy(y_flat + nVx,        L->cor_v, nVz * sizeof(double));
    memcpy(y_flat + nVx + nVz,  L->cor_p, nP  * sizeof(double));
}

// Apply the V-cycle preconditioner: given r_flat, produce z_flat ≈ A^{-1} r.
static void apply_precond_flat(MultigridHierarchy *H, const double *r_flat, double *z_flat) {
    MultigridLevel *L = &H->levels[0];
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    // Prepare level: zero guess, rhs = r.
    memset(L->Vx, 0, nVx * sizeof(double));
    memset(L->Vz, 0, nVz * sizeof(double));
    memset(L->P , 0, nP  * sizeof(double));
    memcpy(L->rhs_u, r_flat,              nVx * sizeof(double));
    memcpy(L->rhs_v, r_flat + nVx,        nVz * sizeof(double));
    memcpy(L->rhs_p, r_flat + nVx + nVz,  nP  * sizeof(double));
    Vcycle(H, 0);
    memcpy(z_flat,              L->Vx, nVx * sizeof(double));
    memcpy(z_flat + nVx,        L->Vz, nVz * sizeof(double));
    memcpy(z_flat + nVx + nVz,  L->P,  nP  * sizeof(double));
}

int FGMRES_GMG(MultigridHierarchy *H, double tol, int restart, int max_restarts,
               FGMRESStats *out_stats) {
    if (H == NULL || H->n_levels < 1) return -1;
    MultigridLevel *L = &H->levels[0];
    const size_t n = level_total_dof(L);
    if (restart < 1) restart = 30;
    if (max_restarts < 1) max_restarts = 30;

    // Snapshot the RHS (needed repeatedly).
    double *b = (double *)malloc(n * sizeof(double));
    double *x = (double *)malloc(n * sizeof(double));
    double *r = (double *)malloc(n * sizeof(double));
    double *w = (double *)malloc(n * sizeof(double));
    if (b == NULL || x == NULL || r == NULL || w == NULL) {
        free(b); free(x); free(r); free(w);
        return -1;
    }
    level_pack(L, L->rhs_u, L->rhs_v, L->rhs_p, b);
    level_pack(L, L->Vx,    L->Vz,    L->P,     x);

    // FGMRES workspace: V[0..m], Z[0..m-1], H (restart+1 x restart), cs, sn, g.
    const int m = restart;
    double **V  = (double **)malloc((size_t)(m + 1) * sizeof(double *));
    double **Z  = (double **)malloc((size_t)m       * sizeof(double *));
    if (V == NULL || Z == NULL) { free(b); free(x); free(r); free(w); free(V); free(Z); return -1; }
    for (int k = 0; k <= m; ++k) V[k] = (double *)calloc(n, sizeof(double));
    for (int k = 0; k <  m; ++k) Z[k] = (double *)calloc(n, sizeof(double));

    double *Hm = (double *)calloc((size_t)(m + 1) * (size_t)m, sizeof(double));
    double *cs = (double *)calloc((size_t)m,       sizeof(double));
    double *sn = (double *)calloc((size_t)m,       sizeof(double));
    double *g  = (double *)calloc((size_t)(m + 1), sizeof(double));

    // Relative-reduction tolerance (add-gmg-upleg-fix D4): the predicate
    // tests ‖r_k‖ / ‖r_0‖ ≤ tol, where ‖r_0‖ is the initial residual of
    // this FGMRES invocation (computed lazily on the first restart so x_0
    // can be any caller-provided warm start). Log messages downstream
    // report the same quantity so "converged" is always numerically
    // consistent with `gmg_fgmres_tol`.
    double r0_norm   = 0.0;
    double target    = 0.0;
    int    have_r0   = 0;
    int total_iter = 0;
    int converged = 0;
    double rnorm = 0.0;

    for (int restart_idx = 0; restart_idx < max_restarts && !converged; ++restart_idx) {
        // r = b - A*x
        apply_A_flat(L, x, w);
        for (size_t k = 0; k < n; ++k) r[k] = b[k] - w[k];
        rnorm = vec_norm(r, n);
        if (!have_r0) {
            r0_norm = rnorm;
            // If the initial residual is already zero (x_0 is exact), the
            // ratio r/r0 is formally undefined; the predicate degenerates
            // to testing the absolute residual against tol, which is the
            // correct behaviour for this edge case.
            target  = tol * (r0_norm > 0.0 ? r0_norm : 1.0);
            have_r0 = 1;
        }
        if (rnorm <= target) { converged = 1; break; }
        // V[0] = r / rnorm
        const double inv = 1.0 / rnorm;
        for (size_t k = 0; k < n; ++k) V[0][k] = r[k] * inv;
        // g = [rnorm, 0, 0, ..., 0]
        memset(g, 0, (size_t)(m + 1) * sizeof(double));
        g[0] = rnorm;

        int j_final = m;
        for (int j = 0; j < m; ++j) {
            // Z[j] = M^{-1} V[j]
            apply_precond_flat(H, V[j], Z[j]);
            // w = A Z[j]
            apply_A_flat(L, Z[j], w);
            // Modified Gram-Schmidt against V[0..j]
            for (int i = 0; i <= j; ++i) {
                double h_ij = vec_dot(w, V[i], n);
                Hm[(size_t)i * m + j] = h_ij;
                for (size_t k = 0; k < n; ++k) w[k] -= h_ij * V[i][k];
            }
            double h_jp1_j = vec_norm(w, n);
            Hm[(size_t)(j + 1) * m + j] = h_jp1_j;
            if (h_jp1_j > 1e-300) {
                double inv2 = 1.0 / h_jp1_j;
                for (size_t k = 0; k < n; ++k) V[j + 1][k] = w[k] * inv2;
            } else {
                // Happy breakdown: solution lies in current subspace.
                j_final = j + 1;
            }
            // Apply previous Givens rotations to column j of Hm.
            for (int i = 0; i < j; ++i) {
                double h1 = Hm[(size_t)i       * m + j];
                double h2 = Hm[(size_t)(i + 1) * m + j];
                Hm[(size_t)i       * m + j] =  cs[i] * h1 + sn[i] * h2;
                Hm[(size_t)(i + 1) * m + j] = -sn[i] * h1 + cs[i] * h2;
            }
            // Compute new Givens rotation for (Hm[j, j], Hm[j+1, j]).
            double h_jj  = Hm[(size_t)j       * m + j];
            double h_jpj = Hm[(size_t)(j + 1) * m + j];
            double nu    = hypot(h_jj, h_jpj);
            if (nu < 1e-300) nu = 1.0;
            cs[j] = h_jj  / nu;
            sn[j] = h_jpj / nu;
            Hm[(size_t)j       * m + j] = nu;
            Hm[(size_t)(j + 1) * m + j] = 0.0;
            // Update g.
            double g_j = g[j];
            g[j    ] =  cs[j] * g_j;
            g[j + 1] = -sn[j] * g_j;

            ++total_iter;
            rnorm = fabs(g[j + 1]);
            if (rnorm <= target) { j_final = j + 1; converged = 1; break; }
            if (h_jp1_j <= 1e-300) break;
        }

        // Solve upper-triangular Hm(0..j_final, 0..j_final-1) y = g(0..j_final-1)
        int jf = converged ? j_final : j_final;
        if (!converged && jf == m) jf = m;
        double y[1024];
        if (jf > 1024) jf = 1024;
        for (int i = jf - 1; i >= 0; --i) {
            double s = g[i];
            for (int k = i + 1; k < jf; ++k) s -= Hm[(size_t)i * m + k] * y[k];
            y[i] = s / Hm[(size_t)i * m + i];
        }
        // x = x + sum_i y[i] * Z[i]
        for (int i = 0; i < jf; ++i) {
            double yi = y[i];
            for (size_t k = 0; k < n; ++k) x[k] += yi * Z[i][k];
        }
        const double rel = (r0_norm > 0.0) ? (rnorm / r0_norm) : rnorm;
        LOG_INFO("GMG-FGMRES: restart %d, total_iter %d, final_res_rel %.3e (tol %.3e)",
                 restart_idx, total_iter, rel, tol);
    }

    // Store result back into level Vx/Vz/P.
    level_unpack(L, x, L->Vx, L->Vz, L->P);

    if (out_stats != NULL) {
        out_stats->iterations           = total_iter;
        out_stats->restarts             = (total_iter + m - 1) / m;
        out_stats->initial_residual     = r0_norm;
        out_stats->final_residual       = rnorm;
        out_stats->fell_back_to_cholmod = 0;
    }

    for (int k = 0; k <= m; ++k) free(V[k]);
    for (int k = 0; k <  m; ++k) free(Z[k]);
    free(V); free(Z);
    free(Hm); free(cs); free(sn); free(g);
    free(b); free(x); free(r); free(w);

    return converged ? 0 : 1;
}

// ---------- Mesh bridge (Section 10, task 10) ---------------------------
//
// The bridge is deliberately narrow: it moves viscosity, the current
// velocity/pressure iterate, and boundary values between MDOODZ's
// ghost-padded Vx/Vz arrays (shapes Nx*(Nz+1) and (Nx+1)*Nz) and our
// interior-only Vx/Vz arrays (Nx*Ncz, Ncx*Nz). The index maps are
// documented in MultigridLevels.h (BuildActiveMask).
//
// The rediscretization uses Gerya 2010's variable-viscosity Picard
// stencil against mesh->eta_n (cell centres) and mesh->eta_s (vertices);
// this matches MDOODZ's assembled operator exactly in the isotropic,
// incompressible, pure-Dirichlet regime (SolVi etc.) and serves as a
// robust preconditioner everywhere else. MDOODZ's more elaborate BC
// tags (2, 11, 13, -2, -12) are mapped as follows:
//   - active_Vx/Vz/P set by BuildActiveMask (tag 30 only → inactive)
//   - Dirichlet-valued boundary faces pull their fixed value from
//     mesh->BCu.val / mesh->BCv.val (or the caller-provided initial
//     mesh->u_in / mesh->v_in, which MDOODZ has already populated with
//     the boundary values on its own outer ring).
// For tags 2/13 (Neumann) we currently leave the stencil as-is, which
// is correct for homogeneous Neumann but under-constraints in general.
// SolveStokesGMG returns a non-zero code in that regime so the caller
// falls back to the direct solve; see the design doc for the extension
// path.

int PopulateLevelFromMesh(MultigridLevel *L, const grid *mesh) {
    if (L == NULL || mesh == NULL) return -1;
    if (L->Nx != mesh->Nx || L->Nz != mesh->Nz) {
        LOG_ERR("GMG: mesh/level shape mismatch (mesh %dx%d, level %dx%d)",
                mesh->Nx, mesh->Nz, L->Nx, L->Nz);
        return -1;
    }
    const int Nx  = L->Nx, Nz = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const int nzvx = Nz + 1;
    const int nxvz = Nx + 1;

    // --- viscosity ---
    // mesh->eta_n has the same (Ncx, Ncz) layout as our etan.
    // mesh->eta_s has the same (Nx, Nz) layout as our etas.
    memcpy(L->etan, mesh->eta_n, (size_t)Ncx * Ncz * sizeof(double));
    memcpy(L->etas, mesh->eta_s, (size_t)Nx  * Nz  * sizeof(double));

    // --- Vx iterate / boundary values ---
    // Our Vx layout is Nx*Ncz with j in [0, Ncz). MDOODZ's mesh->u_in has
    // layout Nx*(Nz+1) with l in [0, Nz]. Map l = j + 1.
    for (int j = 0; j < Ncz; ++j) {
        int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            L->Vx[MdIdx_Vx(i, j, Nx)] = mesh->u_in[(size_t)i + (size_t)l * Nx];
        }
    }

    // --- Vz iterate / boundary values ---
    // Our Vz layout is Ncx*Nz with i in [0, Ncx). MDOODZ's mesh->v_in has
    // layout (Nx+1)*Nz with k in [0, Nx]. Map k = i + 1.
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            int k = i + 1;
            L->Vz[MdIdx_Vz(i, j, Ncx)] = mesh->v_in[(size_t)k + (size_t)j * nxvz];
        }
    }
    (void)nzvx;

    // --- Pressure iterate ---
    memcpy(L->P, mesh->p_in, (size_t)Ncx * Ncz * sizeof(double));

    // --- Active mask ---
    BuildActiveMask(L, mesh->BCu.type, mesh->BCv.type, mesh->BCp.type);
    return 0;
}

int BuildMeshStokesRHS(MultigridLevel *L, const grid *mesh) {
    if (L == NULL || mesh == NULL) return -1;
    const int Nx  = L->Nx, Nz = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const int nxvz = Nx + 1;
    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    const size_t nP  = (size_t)Ncx * Ncz;

    memset(L->rhs_u, 0, nVx * sizeof(double));
    memset(L->rhs_v, 0, nVz * sizeof(double));
    memset(L->rhs_p, 0, nP  * sizeof(double));

    // Interior (solved) rows: body forces and continuity source.
    // MDOODZ loads mesh->roger_x / mesh->roger_z / mesh->rhs_p with the
    // cell/face-weighted body-force contributions that become the RHS
    // for solved DOFs. For outer-ring / inactive faces we want the
    // identity row's RHS to equal the current iterate (so that
    // StokesResidual = 0 at those rows — see StokesApplyA's identity
    // branches).
    for (int j = 0; j < Ncz; ++j) {
        int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            size_t m  = (size_t)i + (size_t)l * Nx;           // mesh index
            size_t ix = MdIdx_Vx(i, j, Nx);
            if (i == 0 || i == Nx - 1 || !L->act_Vx[ix]) {
                L->rhs_u[ix] = L->Vx[ix];
            } else {
                L->rhs_u[ix] = mesh->roger_x[m];
            }
        }
    }
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            int k = i + 1;
            size_t m  = (size_t)k + (size_t)j * nxvz;
            size_t ix = MdIdx_Vz(i, j, Ncx);
            if (j == 0 || j == Nz - 1 || !L->act_Vz[ix]) {
                L->rhs_v[ix] = L->Vz[ix];
            } else {
                L->rhs_v[ix] = mesh->roger_z[m];
            }
        }
    }
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            size_t m  = MdIdx_P(i, j, Ncx);
            if (!L->act_P[m]) {
                L->rhs_p[m] = L->P[m];
            } else {
                L->rhs_p[m] = mesh->rhs_p[m];
            }
        }
    }
    return 0;
}

int UnpackLevelToMesh(const MultigridLevel *L, grid *mesh) {
    if (L == NULL || mesh == NULL) return -1;
    if (L->Nx != mesh->Nx || L->Nz != mesh->Nz) return -1;
    const int Nx  = L->Nx, Nz = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const int nxvz = Nx + 1;

    // Write Vx back; leave inactive / boundary faces untouched so caller
    // sees their MDOODZ BC values as they were.
    for (int j = 0; j < Ncz; ++j) {
        int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            size_t ix = MdIdx_Vx(i, j, Nx);
            if (!L->act_Vx[ix]) continue;
            if (i == 0 || i == Nx - 1) continue;
            mesh->u_in[(size_t)i + (size_t)l * Nx] = L->Vx[ix];
        }
    }
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            size_t ix = MdIdx_Vz(i, j, Ncx);
            if (!L->act_Vz[ix]) continue;
            if (j == 0 || j == Nz - 1) continue;
            int k = i + 1;
            mesh->v_in[(size_t)k + (size_t)j * nxvz] = L->Vz[ix];
        }
    }
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            size_t ix = MdIdx_P(i, j, Ncx);
            if (!L->act_P[ix]) continue;
            mesh->p_in[ix] = L->P[ix];
        }
    }
    return 0;
}

// Restrict fine-level viscosity & masks all the way down the hierarchy.
static void propagate_hierarchy_coefficients(MultigridHierarchy *H) {
    for (int k = 0; k < H->n_levels - 1; ++k) {
        RestrictViscosity(&H->levels[k], &H->levels[k + 1]);
        RestrictActiveMask(&H->levels[k], &H->levels[k + 1]);
    }
}

// ---------- Compressed-equation-space bridge (Section 15 / D11) ---------
//
// MDOODZ's direct solvers operate on the compressed equation vector
// produced by `BuildStokesOperatorDecoupled` — momentum equations
// numbered 0…neq_mom-1 by walking the full u_in/v_in arrays and keeping
// only cells with BCu/BCv.type ∈ {-1, 2, -2} (or {-12} for periodic
// ghosts that mirror their interior twin), then continuity equations
// numbered neq_mom…neq-1 over pressure cells with BCp.type == -1.
// StokesA->b / StokesA->F are sized `neq_mom`, StokesC->b /
// StokesC->F are sized `neq_cont`, and BOTH are celvol-scaled (see the
// `*= celvol` lines in StokesAssemblyDecoupled.c:164-166 / 1121).
//
// These helpers convert between that compressed form and the GMG
// level's full-array layout so SolveStokesGMG can be used
// interchangeably with DirectStokesDecoupled from the nonlinear driver.
// The bridge matvec lives in raw PDE scaling (see
// `StokesApplyA_MDOODZ_bridge` above), so we divide by celvol on the
// way in.
static int BuildLevelRHSFromCompressed(MultigridLevel *L,
                                       const SparseMat *Stokes,
                                       const double *rhs_u,
                                       const double *rhs_p) {
    if (L == NULL || Stokes == NULL || rhs_u == NULL || rhs_p == NULL) return -1;
    const int Nx  = L->Nx,  Nz  = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const int nxvz = Nx + 1;
    const double celvol = L->dx * L->dz;
    const double inv_celvol = 1.0 / celvol;

    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    const size_t nP  = (size_t)Ncx * Ncz;
    memset(L->rhs_u, 0, nVx * sizeof(double));
    memset(L->rhs_v, 0, nVz * sizeof(double));
    memset(L->rhs_p, 0, nP  * sizeof(double));

    // --- Vx: full-array index m = i + l*Nx with l = j + 1 in MDOODZ padded layout.
    for (int j = 0; j < Ncz; ++j) {
        const int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            const size_t m  = (size_t)i + (size_t)l * Nx;
            const size_t ix = MdIdx_Vx(i, j, Nx);
            const int eqn   = Stokes->eqn_u[m];
            if (i == 0 || i == Nx - 1 || !L->act_Vx[ix] || eqn < 0 || eqn >= Stokes->neq_mom) {
                // Identity row: A·x = x. We keep rhs equal to the
                // iterate so `StokesResidual` sees a zero residual there
                // — identical handling to the textbook adapter test path.
                L->rhs_u[ix] = L->Vx[ix];
            } else {
                // Sign flip: the bridge matvec negates momentum rows to
                // match the textbook's `+∇·σ - ∇p = f` convention, so
                // the rhs must also flip for the residual `rhs - A·x`
                // to reduce as the solve progresses. See the comment
                // block in `StokesApplyA_MDOODZ_bridge` for the
                // full derivation.
                L->rhs_u[ix] = -rhs_u[eqn] * inv_celvol;
            }
        }
    }

    // --- Vz: full-array index m = k + j*nxvz with k = i + 1.
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const int k    = i + 1;
            const size_t m = (size_t)k + (size_t)j * nxvz;
            const size_t ix = MdIdx_Vz(i, j, Ncx);
            const int eqn  = Stokes->eqn_v[m];
            if (j == 0 || j == Nz - 1 || !L->act_Vz[ix] || eqn < 0 || eqn >= Stokes->neq_mom) {
                L->rhs_v[ix] = L->Vz[ix];
            } else {
                L->rhs_v[ix] = -rhs_u[eqn] * inv_celvol;
            }
        }
    }

    // --- P: eqn_p maps into [neq_mom, neq); caller's rhs_p[k] holds
    //       continuity equation k = (global eqn) - neq_mom.
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const size_t m   = MdIdx_P(i, j, Ncx);
            const int eqn    = Stokes->eqn_p[m];
            if (!L->act_P[m] || eqn < Stokes->neq_mom || eqn >= Stokes->neq) {
                L->rhs_p[m] = L->P[m];
            } else {
                L->rhs_p[m] = rhs_p[eqn - Stokes->neq_mom] * inv_celvol;
            }
        }
    }
    return 0;
}

// Inverse of BuildLevelRHSFromCompressed: extract the solution back into
// the caller's compressed vector. Writes only equations that exist
// (eqn >= 0); does NOT touch other entries.
static int ExtractLevelToCompressed(const MultigridLevel *L,
                                    const SparseMat *Stokes,
                                    double *sol) {
    if (L == NULL || Stokes == NULL || sol == NULL) return -1;
    const int Nx  = L->Nx,  Nz  = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;
    const int nxvz = Nx + 1;

    for (int j = 0; j < Ncz; ++j) {
        const int l = j + 1;
        for (int i = 0; i < Nx; ++i) {
            const size_t m  = (size_t)i + (size_t)l * Nx;
            const size_t ix = MdIdx_Vx(i, j, Nx);
            const int eqn = Stokes->eqn_u[m];
            if (eqn >= 0 && eqn < Stokes->neq_mom) sol[eqn] = L->Vx[ix];
        }
    }
    for (int j = 0; j < Nz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const int k = i + 1;
            const size_t m  = (size_t)k + (size_t)j * nxvz;
            const size_t ix = MdIdx_Vz(i, j, Ncx);
            const int eqn = Stokes->eqn_v[m];
            if (eqn >= 0 && eqn < Stokes->neq_mom) sol[eqn] = L->Vz[ix];
        }
    }
    for (int j = 0; j < Ncz; ++j) {
        for (int i = 0; i < Ncx; ++i) {
            const size_t m = MdIdx_P(i, j, Ncx);
            const int eqn = Stokes->eqn_p[m];
            if (eqn >= Stokes->neq_mom && eqn < Stokes->neq) sol[eqn] = L->P[m];
        }
    }
    return 0;
}

// ---------- MDOODZ dispatch (Section 10) ---------------------------------

int SolveStokesGMG(SparseMat *StokesA, SparseMat *StokesB,
                   SparseMat *StokesC, SparseMat *StokesD,
                   SparseMat *Stokes,  DirectSolver *PardisoStokes,
                   double *rhs_u, double *rhs_p, double *x,
                   params model, grid *mesh, scale scaling) {
    (void)StokesA; (void)StokesB; (void)StokesD;
    (void)PardisoStokes; (void)scaling;

    if (mesh == NULL) return -1;
    (void)StokesC;
    const int Nx = mesh->Nx, Nz = mesh->Nz;
    if (Nx < MDOODZ_GMG_MIN_SIDE || Nz < MDOODZ_GMG_MIN_SIDE) {
        LOG_WARN("GMG: grid %dx%d below minimum (%d); falling back to CHOLMOD.",
                 Nx, Nz, MDOODZ_GMG_MIN_SIDE);
        return -2;
    }
    // Newton mode is now supported end-to-end (tasks 15.19-15.21 /
    // design D11 Phase 2): the fine-level matvec routes through the
    // Newton Jacobian branch of `ApplyStokesOperatorMDOODZ`, while the
    // Vanka block returns the symmetric Picard part of the smoother
    // (per D7). Coarse levels retain their Picard rediscretization;
    // this is consistent with the D7 symmetric-smoothing rule across
    // the hierarchy.
    // Compressed-space operation requires the equation maps owned by
    // `Stokes` plus the caller-provided rhs and output buffers; without
    // them we cannot participate in the standard MDOODZ
    // solver-swap-in contract (DirectStokesDecoupled,
    // KSPStokesDecoupled, KillerSolver), so cleanly request fallback
    // instead of silently producing a wrong answer.
    if (Stokes == NULL || Stokes->eqn_u == NULL || Stokes->eqn_v == NULL ||
        Stokes->eqn_p == NULL || rhs_u == NULL || rhs_p == NULL || x == NULL) {
        LOG_WARN("GMG: missing Stokes equation maps or compressed rhs/sol; "
                 "falling back to CHOLMOD.");
        return -1;
    }

    // Keep a mutable copy of `model` with stack lifetime so we can hand
    // its address to the stencil bridge (D11 / task 15.7) without using
    // a static (which would not be thread-safe). `model` arrives by
    // value so this is already an independent copy.
    params model_local = model;

    // Allocate a hierarchy sized to the fine mesh.
    MultigridHierarchy H;
    const double dx_fine = mesh->dx, dz_fine = mesh->dz;
    if (MultigridHierarchyAllocate(&H, Nx, Nz, dx_fine, dz_fine, model.gmg_levels) != 0) {
        LOG_ERR("GMG: hierarchy allocation failed");
        return -1;
    }
    H.nu_pre  = (model.gmg_nu_pre  >= 1) ? model.gmg_nu_pre  : 2;
    H.nu_post = (model.gmg_nu_post >= 1) ? model.gmg_nu_post : 2;
    H.use_jacobian = (model.Newton != 0) ? 1 : 0;

    // Populate viscosity / BC masks from mesh and restrict down.
    // Velocities / pressure start at zero — the caller's rhs already
    // carries any Dirichlet contributions via AddCoeff3's bbc folding.
    if (PopulateLevelFromMesh(&H.levels[0], mesh) != 0) {
        LOG_ERR("GMG: mesh->level population failed");
        MultigridHierarchyFree(&H);
        return -1;
    }
    {
        MultigridLevel *L0 = &H.levels[0];
        memset(L0->Vx, 0, (size_t)L0->Nx  * L0->Ncz * sizeof(double));
        memset(L0->Vz, 0, (size_t)L0->Ncx * L0->Nz  * sizeof(double));
        memset(L0->P , 0, (size_t)L0->Ncx * L0->Ncz * sizeof(double));
    }

    // Engage the MDOODZ stencil bridge on the fine level only (D11 /
    // task 15.7). `StokesApplyA` routes through
    // `ApplyStokesOperatorMDOODZ` for H.levels[0] so the FGMRES outer
    // operator matches StokesAssemblyDecoupled.c bit-identically.
    // Coarse levels keep textbook Picard rediscretisation on the
    // restricted viscosity (design D4). The hierarchy borrows
    // mesh/model_local by pointer; model_local is a stack copy that
    // outlives the V-cycle (both live on SolveStokesGMG's frame).
    //
    // The Vanka 5×5 block builder stays on the TEXTBOOK path at L0:
    // add-gmg-upleg-fix §4 localised a 7.55× per-sweep L0 pre-smooth
    // amplification to `VankaBlockAssembleSolve_MDOODZ_bridge`, which
    // disappears (ratio 0.174×, UplegAmplificationBound passes) once
    // the bridge Vanka path is bypassed. Outer FGMRES fidelity is
    // preserved because `StokesApplyA` still dispatches through the
    // bridge. A follow-up change is tracked to fix the bridge 5×5
    // block in-place; until then the two dispatches are kept split.
    H.levels[0].use_mdoodz_matvec_outer = 1;
    H.levels[0].use_mdoodz_matvec_vanka = 0;
    H.levels[0].mesh_ref                = (void *)mesh;
    H.levels[0].model_ref               = (void *)&model_local;

    // Build RHS in compressed equation space from the caller's vectors.
    // IMPORTANT: this must come AFTER PopulateLevelFromMesh (for the
    // active masks) and AFTER the L0->Vx/Vz/P zero-init (for the
    // identity rows that fall back to the iterate value).
    if (BuildLevelRHSFromCompressed(&H.levels[0], Stokes, rhs_u, rhs_p) != 0) {
        LOG_ERR("GMG: rhs translation (compressed -> full) failed");
        MultigridHierarchyFree(&H);
        return -1;
    }

    propagate_hierarchy_coefficients(&H);

    // Arm the V-cycle dumper if requested and it has not fired yet this
    // process. Exactly the first V-cycle after arming records; the
    // finish call below marks `fired` so no subsequent solver invocation
    // dumps again (see design D9: first V-cycle of first FGMRES of first
    // Picard of first time step).
    const int vcycle_dump_armed = gmg_vcycle_dump_arm(model.writer_subfolder,
                                                      model.gmg_dump_vcycle);

    // Run the solver. `gmg_standalone` toggles between bare V-cycles
    // (diagnostic) and FGMRES preconditioned by V-cycles (production).
    const double tol    = (model.gmg_fgmres_tol > 0.0) ? model.gmg_fgmres_tol : 1e-8;
    const int restart   = (model.gmg_fgmres_restart > 0) ? model.gmg_fgmres_restart : 30;
    const int max_restarts = (model.gmg_fgmres_max_restarts >= 1) ? model.gmg_fgmres_max_restarts : 20;

    // add-gmg-upleg-fix §1.3: echo the active GMG configuration once per
    // Stokes solve so STATUS receipts, perf dumps, and ctest logs all see
    // which knobs were in effect (including the new max_restarts).
    LOG_INFO("GMG-FGMRES config: levels=%d (0=auto), nu_pre=%d, nu_post=%d, "
             "restart=%d, max_restarts=%d, tol=%.3e, standalone=%d, dump=%d",
             model.gmg_levels, H.nu_pre, H.nu_post,
             restart, max_restarts, tol,
             model.gmg_standalone, model.gmg_dump_vcycle);

    FGMRESStats stats = {0};
    int rc = 0;
    if (model.gmg_standalone == 1) {
        // Bare V-cycles until the residual is below tol or 100 sweeps.
        const int max_sweeps = 100;
        int sweep = 0;
        StokesResidual(&H.levels[0]);
        double r0 = StokesResidualNorm(&H.levels[0]);
        double r  = r0;
        while (sweep < max_sweeps && r > tol * (r0 > 0.0 ? r0 : 1.0)) {
            Vcycle(&H, 0);
            StokesResidual(&H.levels[0]);
            r = StokesResidualNorm(&H.levels[0]);
            ++sweep;
        }
        stats.iterations       = sweep;
        stats.initial_residual = r0;
        stats.final_residual   = r;
        rc = (r > tol * (r0 > 0.0 ? r0 : 1.0)) ? 1 : 0;
    } else {
        rc = FGMRES_GMG(&H, tol, restart, max_restarts, &stats);
    }

    // Close out the one-shot dumper. It has already recorded the first
    // (and only) V-cycle of this solve; finish logs the total count and
    // flips the process-wide `fired` flag so no subsequent invocation
    // dumps again.
    if (vcycle_dump_armed) gmg_vcycle_dump_finish();

    // add-gmg-upleg-fix D4/§1.6-1.7: the predicate tests ‖r_k‖/‖r_0‖ ≤ tol;
    // emit `final_res_rel` (the same quantity) in both success and failure
    // paths so "converged" in the log always means "rel ≤ tol" numerically.
    const double final_res_rel = (stats.initial_residual > 0.0)
                                     ? (stats.final_residual / stats.initial_residual)
                                     : stats.final_residual;

    if (rc != 0) {
        LOG_WARN("GMG-FGMRES did not converge: iters=%d, final_res_rel=%.3e, tol=%.3e "
                 "(‖r_0‖=%.3e, ‖r_k‖=%.3e)",
                 stats.iterations, final_res_rel, tol,
                 stats.initial_residual, stats.final_residual);
        MultigridHierarchyFree(&H);
        return -4;
    }

    LOG_INFO("GMG-FGMRES converged: iters=%d, restarts=%d, final_res_rel=%.3e, tol=%.3e",
             stats.iterations, stats.restarts, final_res_rel, tol);

    // Copy the result from the GMG level's full-array layout back into
    // the caller's compressed-equation-space `x` vector. The caller
    // (DirectStokesDecoupled style) is responsible for feeding `x`
    // into line-search / ExtractSolutions — we deliberately do NOT
    // write mesh->u_in/v_in/p_in here because:
    //   (a) in defect-correction mode `x` is the Newton delta and
    //       must not be applied blindly to mesh->u_in;
    //   (b) the caller's line search may scale the increment by a
    //       non-unit alpha before committing it to the mesh.
    if (ExtractLevelToCompressed(&H.levels[0], Stokes, x) != 0) {
        LOG_ERR("GMG: compressed solution extraction failed");
        MultigridHierarchyFree(&H);
        return -1;
    }

    MultigridHierarchyFree(&H);
    return 0;
}
