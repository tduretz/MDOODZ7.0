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
// GMG hierarchy + transfer operators. See MultigridLevels.h for layout.

#include "MultigridLevels.h"
#include "mdoodz-log.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "umfpack.h"

// --- Small helpers -------------------------------------------------------

static int imin(int a, int b) { return a < b ? a : b; }
static int imax(int a, int b) { return a > b ? a : b; }

static void *xcalloc(size_t n, size_t s) {
    void *p = calloc(n, s);
    if (p == NULL && n * s > 0) {
        LOG_ERR("GMG: calloc failed for %zu bytes", n * s);
        exit(1);
    }
    return p;
}

// Fill a char array with 1 (all active).
static void fill_active(char *m, size_t n) {
    if (m == NULL) return;
    for (size_t k = 0; k < n; ++k) m[k] = 1;
}

// --- Level count policy (task 2.1) ---------------------------------------

int MultigridComputeDefaultLevelCount(int Nx, int Nz, int gmg_levels_override) {
    if (Nx < 3 || Nz < 3) return 1;

    const int min_side = imin(Nx, Nz);
    // Absolute upper bound on levels: we cannot halve below ~2 cells per side.
    int max_possible = 1;
    int s = min_side;
    while (s > 2) { s = (s + 1) / 2; ++max_possible; }

    if (gmg_levels_override > 0) {
        int clamped = imax(2, imin(gmg_levels_override, max_possible));
        return clamped;
    }

    // Auto-policy: keep halving while the CURRENT level's min side stays
    // >= MIN_SIDE; once a coarsened level falls below MIN_SIDE it becomes
    // the coarsest (CHOLMOD-handled). The 5000-DOF clause is an invariant
    // the output must satisfy, not a separate trigger. For huge grids
    // (DOF at MIN_SIDE still > FLOOR) the min_side rule is sufficient
    // because CHOLMOD handles any size at the coarsest level; for small
    // enough grids, we always land at DOF << FLOOR.
    int levels = 1;
    int Nx_k = Nx, Nz_k = Nz;
    while (levels < max_possible) {
        if (imin(Nx_k, Nz_k) < MDOODZ_GMG_MIN_SIDE) break; // current already coarsest
        Nx_k = (Nx_k + 1) / 2;
        Nz_k = (Nz_k + 1) / 2;
        ++levels;
    }
    return imax(2, levels);
}

// --- Hierarchy allocation / free (tasks 2.2, 2.3) ------------------------

static int LevelAlloc(MultigridLevel *L, int level, int Nx, int Nz, double dx, double dz) {
    if (Nx < 3 || Nz < 3) {
        LOG_ERR("GMG: level %d too small (Nx=%d, Nz=%d); need >= 3 per side", level, Nx, Nz);
        return -1;
    }
    L->level = level;
    L->Nx  = Nx;      L->Nz  = Nz;
    L->Ncx = Nx - 1;  L->Ncz = Nz - 1;
    L->dx  = dx;      L->dz  = dz;

    const size_t nVx = (size_t)Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    const size_t nES = (size_t)Nx  * Nz;

    L->Vx    = (double *)xcalloc(nVx, sizeof(double));
    L->Vz    = (double *)xcalloc(nVz, sizeof(double));
    L->P     = (double *)xcalloc(nP,  sizeof(double));
    L->rhs_u = (double *)xcalloc(nVx, sizeof(double));
    L->rhs_v = (double *)xcalloc(nVz, sizeof(double));
    L->rhs_p = (double *)xcalloc(nP,  sizeof(double));
    L->res_u = (double *)xcalloc(nVx, sizeof(double));
    L->res_v = (double *)xcalloc(nVz, sizeof(double));
    L->res_p = (double *)xcalloc(nP,  sizeof(double));
    L->cor_u = (double *)xcalloc(nVx, sizeof(double));
    L->cor_v = (double *)xcalloc(nVz, sizeof(double));
    L->cor_p = (double *)xcalloc(nP,  sizeof(double));

    L->etan  = (double *)xcalloc(nP,  sizeof(double));
    L->etas  = (double *)xcalloc(nES, sizeof(double));

    L->act_P  = (char *)xcalloc(nP,  sizeof(char));
    L->act_Vx = (char *)xcalloc(nVx, sizeof(char));
    L->act_Vz = (char *)xcalloc(nVz, sizeof(char));
    // Default: everything active.
    fill_active(L->act_P,  nP);
    fill_active(L->act_Vx, nVx);
    fill_active(L->act_Vz, nVz);

    return 0;
}

static void LevelFree(MultigridLevel *L) {
    if (L == NULL) return;
    free(L->Vx);    free(L->Vz);    free(L->P);
    free(L->rhs_u); free(L->rhs_v); free(L->rhs_p);
    free(L->res_u); free(L->res_v); free(L->res_p);
    free(L->cor_u); free(L->cor_v); free(L->cor_p);
    free(L->etan);  free(L->etas);
    free(L->act_P); free(L->act_Vx); free(L->act_Vz);
    free(L->md_u_pad); free(L->md_v_pad); free(L->md_p_pad);
    memset(L, 0, sizeof(*L));
}

int MultigridHierarchyAllocate(MultigridHierarchy *H, int Nx_fine, int Nz_fine,
                               double dx_fine, double dz_fine,
                               int gmg_levels_override) {
    if (H == NULL) return -1;
    memset(H, 0, sizeof(*H));

    const int n_levels = MultigridComputeDefaultLevelCount(Nx_fine, Nz_fine, gmg_levels_override);
    if (n_levels < 1) {
        LOG_ERR("GMG: invalid level count %d", n_levels);
        return -1;
    }
    H->n_levels = n_levels;
    H->levels   = (MultigridLevel *)xcalloc((size_t)n_levels, sizeof(MultigridLevel));

    int Nx = Nx_fine, Nz = Nz_fine;
    double dx = dx_fine, dz = dz_fine;
    for (int k = 0; k < n_levels; ++k) {
        if (LevelAlloc(&H->levels[k], k, Nx, Nz, dx, dz) != 0) {
            MultigridHierarchyFree(H);
            return -1;
        }
        // Prepare dimensions for the next (coarser) level.
        if (k < n_levels - 1) {
            int Nx_next = (Nx + 1) / 2;
            int Nz_next = (Nz + 1) / 2;
            // Coarse grid still spans the same physical domain:
            // fine domain length = (Nx-1)*dx; coarse spacing = fine_len / (Nx_next-1).
            dx = (double)(Nx - 1) * dx / (double)(Nx_next - 1);
            dz = (double)(Nz - 1) * dz / (double)(Nz_next - 1);
            Nx = Nx_next;
            Nz = Nz_next;
        }
    }
    H->coarse_umf_numeric  = NULL;
    H->coarse_factor_valid = 0;
    H->use_jacobian        = 0;
    return 0;
}

void MultigridHierarchyFree(MultigridHierarchy *H) {
    if (H == NULL) return;
    if (H->levels != NULL) {
        for (int k = 0; k < H->n_levels; ++k) LevelFree(&H->levels[k]);
        free(H->levels);
    }
    free(H->coarse_Ap); free(H->coarse_Aj); free(H->coarse_Ax);
    free(H->coarse_eqn_Vx); free(H->coarse_eqn_Vz); free(H->coarse_eqn_P);
    if (H->coarse_umf_numeric != NULL) {
        void *n = H->coarse_umf_numeric;
        umfpack_di_free_numeric(&n);
        H->coarse_umf_numeric = NULL;
    }
    memset(H, 0, sizeof(*H));
}

// --- Transfer operators --------------------------------------------------
//
// Staggering reminder. At the fine level:
//   Vx_f(i, j), i ∈ [0, Nx_f),   j ∈ [0, Ncz_f).
//   Vz_f(i, j), i ∈ [0, Ncx_f),  j ∈ [0, Nz_f).
//   P_f (i, j), i ∈ [0, Ncx_f),  j ∈ [0, Ncz_f).
//
// Under 2x coarsening with ceil-rounding, coarse index k maps to fine
// index 2k (vertex-aligned for Vx's i-axis, Vz's j-axis, and vertex grid)
// or to fine cell (2k, 2k+1) (cell-aligned for P, Vx's j-axis, Vz's i-axis).
//
// For the transfer operators we use boundary-safe full-weighting with
// mirrored indices at the edges. That's standard for a staggered Stokes
// V-cycle and preserves constants exactly.

// Clamp index to valid range.
static int clampi(int i, int lo, int hi) { return i < lo ? lo : (i > hi ? hi : i); }

// Read V with index clamping.
static double V_clamped(const double *V, int i, int j, int Ni, int Nj) {
    int ii = clampi(i, 0, Ni - 1);
    int jj = clampi(j, 0, Nj - 1);
    return V[(size_t)ii + (size_t)jj * Ni];
}

// ---- Vx restriction (vertex-aligned in x, cell-aligned in z) -----------
//
// Coarse Vx has shape Nx_c x Ncz_c with Nx_c = (Nx_f + 1) / 2,
// Ncz_c = Ncz_f / 2 (integer floor) — or more precisely, sized from
// (Nx_c, Ncz_c) = (coarse->Nx, coarse->Ncz).
//
// For a coarse Vx node at (I, J):
//   x-axis: vertex-aligned, fine index i_f = 2*I (clamped)
//   z-axis: cell-aligned,   fine indices j_f ∈ {2J, 2J+1} (averaged)
//
// Full-weighting stencil (symmetric about (i_f, j_f_centre)):
//   (1/16) sum over 3x3 nbhd with weights [1 2 1 / 2 4 2 / 1 2 1]
//   where the z-centre is (2J + 2J+1)/2, i.e. between two fine rows.
//   We evaluate it as a 4x2 window averaging the two z-rows symmetrically.

static void RestrictVx(const MultigridLevel *F, const MultigridLevel *C,
                       const double *Vf, double *Vc) {
    const int Nxf  = F->Nx,  Ncz_f = F->Ncz;
    const int Nxc  = C->Nx,  Ncz_c = C->Ncz;

    for (int J = 0; J < Ncz_c; ++J) {
        // Cell-aligned z: two fine rows (2J, 2J+1).
        int jf0 = clampi(2 * J,     0, Ncz_f - 1);
        int jf1 = clampi(2 * J + 1, 0, Ncz_f - 1);
        for (int I = 0; I < Nxc; ++I) {
            int If = clampi(2 * I, 0, Nxf - 1);
            // 3x2 window in x (vertex-centred) times 2 fine z-rows.
            // Weights on x: 1/4, 1/2, 1/4 (Laplacian row)
            // Weights on z: 1/2, 1/2 (arithmetic mean of two cell rows)
            double sum = 0.0;
            const double wx[3] = {0.25, 0.5, 0.25};
            const int    dx[3] = {-1, 0, +1};
            const double wz[2] = {0.5, 0.5};
            const int    jrows[2] = {jf0, jf1};
            for (int kz = 0; kz < 2; ++kz) {
                for (int kx = 0; kx < 3; ++kx) {
                    sum += wx[kx] * wz[kz] *
                           V_clamped(Vf, If + dx[kx], jrows[kz], Nxf, Ncz_f);
                }
            }
            Vc[MdIdx_Vx(I, J, Nxc)] = sum;
        }
    }
}

// Vz restriction (cell-aligned in x, vertex-aligned in z): symmetric dual.
static void RestrictVz(const MultigridLevel *F, const MultigridLevel *C,
                       const double *Vf, double *Vc) {
    const int Ncx_f = F->Ncx, Nzf = F->Nz;
    const int Ncx_c = C->Ncx, Nzc = C->Nz;
    for (int J = 0; J < Nzc; ++J) {
        int Jf = clampi(2 * J, 0, Nzf - 1);
        for (int I = 0; I < Ncx_c; ++I) {
            int if0 = clampi(2 * I,     0, Ncx_f - 1);
            int if1 = clampi(2 * I + 1, 0, Ncx_f - 1);
            double sum = 0.0;
            const double wz[3] = {0.25, 0.5, 0.25};
            const int    dz[3] = {-1, 0, +1};
            const double wx[2] = {0.5, 0.5};
            const int    icols[2] = {if0, if1};
            for (int kx = 0; kx < 2; ++kx) {
                for (int kz = 0; kz < 3; ++kz) {
                    sum += wx[kx] * wz[kz] *
                           V_clamped(Vf, icols[kx], Jf + dz[kz], Ncx_f, Nzf);
                }
            }
            Vc[MdIdx_Vz(I, J, Ncx_c)] = sum;
        }
    }
}

void RestrictVelocity(const MultigridLevel *fine, MultigridLevel *coarse,
                      const double *fine_V, double *coarse_V, int component) {
    if (component == 0) RestrictVx(fine, coarse, fine_V, coarse_V);
    else                RestrictVz(fine, coarse, fine_V, coarse_V);

    // Dirichlet rows are held as identity by StokesApplyA and the Vanka
    // block; the discrete residual at those rows is identically zero on
    // any consistent iterate. If the restriction window happens to pull
    // in neighbouring interior values (because it is clamped near the
    // boundary), the resulting coarse RHS at the boundary index would
    // otherwise fossilise a non-zero residual that Vanka cannot remove
    // (boundary DOFs keep their identity row in the 5x5 block). Zero
    // those entries explicitly so the boundary remains a "don't care".
    const int Nxc = coarse->Nx, Nzc = coarse->Nz;
    const int Ncx_c = coarse->Ncx, Ncz_c = coarse->Ncz;
    if (component == 0) {
        for (int J = 0; J < Ncz_c; ++J) {
            coarse_V[MdIdx_Vx(0,       J, Nxc)] = 0.0;
            coarse_V[MdIdx_Vx(Nxc - 1, J, Nxc)] = 0.0;
        }
    } else {
        for (int I = 0; I < Ncx_c; ++I) {
            coarse_V[MdIdx_Vz(I, 0,       Ncx_c)] = 0.0;
            coarse_V[MdIdx_Vz(I, Nzc - 1, Ncx_c)] = 0.0;
        }
    }
}

// Pressure restriction: cell-centred -> cell-centred, full-weighting over
// the 2x2 block of fine children that make up one coarse cell. If the
// fine grid has an odd count the boundary children are handled by
// clamping (which preserves constants exactly).
void RestrictPressure(const MultigridLevel *fine, MultigridLevel *coarse,
                      const double *fine_P, double *coarse_P) {
    const int Ncx_f = fine->Ncx, Ncz_f = fine->Ncz;
    const int Ncx_c = coarse->Ncx, Ncz_c = coarse->Ncz;
    for (int J = 0; J < Ncz_c; ++J) {
        for (int I = 0; I < Ncx_c; ++I) {
            int i0 = clampi(2 * I,     0, Ncx_f - 1);
            int i1 = clampi(2 * I + 1, 0, Ncx_f - 1);
            int j0 = clampi(2 * J,     0, Ncz_f - 1);
            int j1 = clampi(2 * J + 1, 0, Ncz_f - 1);
            double s = 0.0;
            s += fine_P[MdIdx_P(i0, j0, Ncx_f)];
            s += fine_P[MdIdx_P(i1, j0, Ncx_f)];
            s += fine_P[MdIdx_P(i0, j1, Ncx_f)];
            s += fine_P[MdIdx_P(i1, j1, Ncx_f)];
            coarse_P[MdIdx_P(I, J, Ncx_c)] = 0.25 * s;
        }
    }
}

// ---- Prolongation (coarse -> fine, additive) ---------------------------

static void ProlongateVxAdd(const MultigridLevel *C, const MultigridLevel *F,
                            const double *Vc, double *Vf) {
    const int Nxc = C->Nx, Ncz_c = C->Ncz;
    const int Nxf = F->Nx, Ncz_f = F->Ncz;
    // Fine node (if, jf) -> fractional coarse coords.
    //  x-axis (vertex): Ic = if / 2 (integer), frac = (if % 2) * 0.5
    //  z-axis (cell):   Jc = jf / 2, frac = (jf % 2) * 0.5
    // Bilinear weights derive from these fractional coords.
    for (int jf = 0; jf < Ncz_f; ++jf) {
        double zf = (jf + 0.5) / 2.0 - 0.5;      // coarse-cell coordinate
        int    Jc = (int)floor(zf);
        double tz = zf - (double)Jc;
        int    J0 = clampi(Jc,     0, Ncz_c - 1);
        int    J1 = clampi(Jc + 1, 0, Ncz_c - 1);
        for (int if_ = 0; if_ < Nxf; ++if_) {
            double xf = if_ * 0.5;              // coarse-vertex coordinate
            int    Ic = (int)floor(xf);
            double tx = xf - (double)Ic;
            int    I0 = clampi(Ic,     0, Nxc - 1);
            int    I1 = clampi(Ic + 1, 0, Nxc - 1);
            double v = (1.0 - tx) * (1.0 - tz) * Vc[MdIdx_Vx(I0, J0, Nxc)]
                     + (      tx) * (1.0 - tz) * Vc[MdIdx_Vx(I1, J0, Nxc)]
                     + (1.0 - tx) * (      tz) * Vc[MdIdx_Vx(I0, J1, Nxc)]
                     + (      tx) * (      tz) * Vc[MdIdx_Vx(I1, J1, Nxc)];
            Vf[MdIdx_Vx(if_, jf, Nxf)] += v;
        }
    }
}

static void ProlongateVzAdd(const MultigridLevel *C, const MultigridLevel *F,
                            const double *Vc, double *Vf) {
    const int Ncx_c = C->Ncx, Nzc = C->Nz;
    const int Ncx_f = F->Ncx, Nzf = F->Nz;
    for (int jf = 0; jf < Nzf; ++jf) {
        double zf = jf * 0.5;
        int    Jc = (int)floor(zf);
        double tz = zf - (double)Jc;
        int    J0 = clampi(Jc,     0, Nzc - 1);
        int    J1 = clampi(Jc + 1, 0, Nzc - 1);
        for (int if_ = 0; if_ < Ncx_f; ++if_) {
            double xf = (if_ + 0.5) / 2.0 - 0.5;
            int    Ic = (int)floor(xf);
            double tx = xf - (double)Ic;
            int    I0 = clampi(Ic,     0, Ncx_c - 1);
            int    I1 = clampi(Ic + 1, 0, Ncx_c - 1);
            double v = (1.0 - tx) * (1.0 - tz) * Vc[MdIdx_Vz(I0, J0, Ncx_c)]
                     + (      tx) * (1.0 - tz) * Vc[MdIdx_Vz(I1, J0, Ncx_c)]
                     + (1.0 - tx) * (      tz) * Vc[MdIdx_Vz(I0, J1, Ncx_c)]
                     + (      tx) * (      tz) * Vc[MdIdx_Vz(I1, J1, Ncx_c)];
            Vf[MdIdx_Vz(if_, jf, Ncx_f)] += v;
        }
    }
}

void ProlongateVelocityAdd(const MultigridLevel *coarse, const MultigridLevel *fine,
                           const double *coarse_V, double *fine_V, int component) {
    if (component == 0) ProlongateVxAdd(coarse, fine, coarse_V, fine_V);
    else                ProlongateVzAdd(coarse, fine, coarse_V, fine_V);
}

// Pressure prolongation: piecewise-constant injection (additive).
void ProlongatePressureAdd(const MultigridLevel *coarse, const MultigridLevel *fine,
                           const double *coarse_P, double *fine_P) {
    const int Ncx_c = coarse->Ncx, Ncz_c = coarse->Ncz;
    const int Ncx_f = fine->Ncx,   Ncz_f = fine->Ncz;
    for (int jf = 0; jf < Ncz_f; ++jf) {
        int Jc = clampi(jf / 2, 0, Ncz_c - 1);
        for (int if_ = 0; if_ < Ncx_f; ++if_) {
            int Ic = clampi(if_ / 2, 0, Ncx_c - 1);
            fine_P[MdIdx_P(if_, jf, Ncx_f)] += coarse_P[MdIdx_P(Ic, Jc, Ncx_c)];
        }
    }
}

// --- Viscosity restriction (harmonic over contrasts, arithmetic else) ---
// Helper: average of up to 4 positive viscosity samples, harmonic when
// the ratio max/min exceeds HARMONIC_RATIO, else arithmetic. n > 0.
static double viscosity_average(const double *samples, const char *mask_samples, int n) {
    double sum_a  = 0.0;
    double sum_hi = 0.0; // sum of 1/eta
    int    count  = 0;
    double mn = 0.0, mx = 0.0;
    for (int k = 0; k < n; ++k) {
        if (mask_samples != NULL && !mask_samples[k]) continue;
        double e = samples[k];
        if (!(e > 0.0)) continue; // skip zero / negative
        if (count == 0) { mn = mx = e; }
        else { if (e < mn) mn = e; if (e > mx) mx = e; }
        sum_a  += e;
        sum_hi += 1.0 / e;
        ++count;
    }
    if (count == 0) return 0.0;
    double arith = sum_a / (double)count;
    if (mn > 0.0 && mx / mn > MDOODZ_GMG_HARMONIC_RATIO) {
        return (double)count / sum_hi;  // harmonic
    }
    return arith;
}

void RestrictViscosity(const MultigridLevel *fine, MultigridLevel *coarse) {
    const int Ncx_f = fine->Ncx, Ncz_f = fine->Ncz;
    const int Ncx_c = coarse->Ncx, Ncz_c = coarse->Ncz;

    // etan: cell-centred, 2x2 block of fine cells -> 1 coarse cell.
    for (int J = 0; J < Ncz_c; ++J) {
        for (int I = 0; I < Ncx_c; ++I) {
            int i0 = clampi(2 * I,     0, Ncx_f - 1);
            int i1 = clampi(2 * I + 1, 0, Ncx_f - 1);
            int j0 = clampi(2 * J,     0, Ncz_f - 1);
            int j1 = clampi(2 * J + 1, 0, Ncz_f - 1);
            double s[4] = {
                fine->etan[MdIdx_etan(i0, j0, Ncx_f)],
                fine->etan[MdIdx_etan(i1, j0, Ncx_f)],
                fine->etan[MdIdx_etan(i0, j1, Ncx_f)],
                fine->etan[MdIdx_etan(i1, j1, Ncx_f)],
            };
            char m[4] = {
                fine->act_P[MdIdx_P(i0, j0, Ncx_f)],
                fine->act_P[MdIdx_P(i1, j0, Ncx_f)],
                fine->act_P[MdIdx_P(i0, j1, Ncx_f)],
                fine->act_P[MdIdx_P(i1, j1, Ncx_f)],
            };
            coarse->etan[MdIdx_etan(I, J, Ncx_c)] = viscosity_average(s, m, 4);
        }
    }

    // etas: vertex-based, Nx_c * Nz_c output. Each coarse vertex is co-located
    // with the fine vertex at (2*I, 2*J) (clamped). Use a 3x3 full-weighting
    // average with harmonic-or-arithmetic rule.
    const int Nxf = fine->Nx, Nzf = fine->Nz;
    const int Nxc = coarse->Nx, Nzc = coarse->Nz;
    for (int J = 0; J < Nzc; ++J) {
        for (int I = 0; I < Nxc; ++I) {
            int If = clampi(2 * I, 0, Nxf - 1);
            int Jf = clampi(2 * J, 0, Nzf - 1);
            double samples[9];
            int n = 0;
            for (int dj = -1; dj <= 1; ++dj) {
                for (int di = -1; di <= 1; ++di) {
                    int ii = clampi(If + di, 0, Nxf - 1);
                    int jj = clampi(Jf + dj, 0, Nzf - 1);
                    samples[n++] = fine->etas[MdIdx_etas(ii, jj, Nxf)];
                }
            }
            coarse->etas[MdIdx_etas(I, J, Nxc)] = viscosity_average(samples, NULL, n);
        }
    }
}

// --- Active mask restriction (task 5.2) ---------------------------------

void RestrictActiveMask(const MultigridLevel *fine, MultigridLevel *coarse) {
    const int Ncx_f = fine->Ncx, Ncz_f = fine->Ncz;
    const int Ncx_c = coarse->Ncx, Ncz_c = coarse->Ncz;
    for (int J = 0; J < Ncz_c; ++J) {
        for (int I = 0; I < Ncx_c; ++I) {
            int i0 = clampi(2 * I,     0, Ncx_f - 1);
            int i1 = clampi(2 * I + 1, 0, Ncx_f - 1);
            int j0 = clampi(2 * J,     0, Ncz_f - 1);
            int j1 = clampi(2 * J + 1, 0, Ncz_f - 1);
            int any_active =
                  fine->act_P[MdIdx_P(i0, j0, Ncx_f)]
                | fine->act_P[MdIdx_P(i1, j0, Ncx_f)]
                | fine->act_P[MdIdx_P(i0, j1, Ncx_f)]
                | fine->act_P[MdIdx_P(i1, j1, Ncx_f)];
            coarse->act_P[MdIdx_P(I, J, Ncx_c)] = (char)(any_active ? 1 : 0);
        }
    }
    // Vx mask: coarse Vx at (I, J) active iff any fine Vx in the 2x1 z-cell
    // slab above/below is active.
    const int Nxf = fine->Nx,  Nxc = coarse->Nx;
    for (int J = 0; J < coarse->Ncz; ++J) {
        int j0 = clampi(2 * J,     0, fine->Ncz - 1);
        int j1 = clampi(2 * J + 1, 0, fine->Ncz - 1);
        for (int I = 0; I < Nxc; ++I) {
            int If = clampi(2 * I, 0, Nxf - 1);
            int any_active = fine->act_Vx[MdIdx_Vx(If, j0, Nxf)]
                           | fine->act_Vx[MdIdx_Vx(If, j1, Nxf)];
            coarse->act_Vx[MdIdx_Vx(I, J, Nxc)] = (char)(any_active ? 1 : 0);
        }
    }
    // Vz mask: coarse Vz at (I, J) active iff any fine Vz in the 1x2 x-slab is active.
    const int Ncx_fz = fine->Ncx;
    for (int J = 0; J < coarse->Nz; ++J) {
        int Jf = clampi(2 * J, 0, fine->Nz - 1);
        for (int I = 0; I < coarse->Ncx; ++I) {
            int i0 = clampi(2 * I,     0, Ncx_fz - 1);
            int i1 = clampi(2 * I + 1, 0, Ncx_fz - 1);
            int any_active = fine->act_Vz[MdIdx_Vz(i0, Jf, Ncx_fz)]
                           | fine->act_Vz[MdIdx_Vz(i1, Jf, Ncx_fz)];
            coarse->act_Vz[MdIdx_Vz(I, J, coarse->Ncx)] = (char)(any_active ? 1 : 0);
        }
    }
}

// --- Utilities ----------------------------------------------------------

void LevelZeroFields(MultigridLevel *L) {
    if (L == NULL) return;
    const size_t nVx = (size_t)L->Nx  * L->Ncz;
    const size_t nVz = (size_t)L->Ncx * L->Nz;
    const size_t nP  = (size_t)L->Ncx * L->Ncz;
    memset(L->Vx,    0, nVx * sizeof(double));
    memset(L->Vz,    0, nVz * sizeof(double));
    memset(L->P,     0, nP  * sizeof(double));
    memset(L->rhs_u, 0, nVx * sizeof(double));
    memset(L->rhs_v, 0, nVz * sizeof(double));
    memset(L->rhs_p, 0, nP  * sizeof(double));
    memset(L->res_u, 0, nVx * sizeof(double));
    memset(L->res_v, 0, nVz * sizeof(double));
    memset(L->res_p, 0, nP  * sizeof(double));
    memset(L->cor_u, 0, nVx * sizeof(double));
    memset(L->cor_v, 0, nVz * sizeof(double));
    memset(L->cor_p, 0, nP  * sizeof(double));
}

double MaskedL2(const double *x, const char *mask, size_t n) {
    double s = 0.0;
    for (size_t k = 0; k < n; ++k) {
        if (mask != NULL && !mask[k]) continue;
        double v = x[k];
        s += v * v;
    }
    return sqrt(s);
}

// --- Active mask from MDOODZ BC tags (task 5.1) --------------------------

// Only tag 30 means "do not compute" for the purposes of our discrete
// operator; everything else (Dirichlet, Neumann, periodic, interior) is
// handled by a live row in the 5x5 Vanka block and StokesApplyA.
static inline char tag_active(signed char t) { return (t == 30) ? 0 : 1; }

void BuildActiveMask(MultigridLevel       *L,
                     const signed char    *BCu_type,
                     const signed char    *BCv_type,
                     const signed char    *BCp_type) {
    if (L == NULL) return;
    const int Nx  = L->Nx,  Nz  = L->Nz;
    const int Ncx = L->Ncx, Ncz = L->Ncz;

    // Default to all-active (e.g. if caller passes NULL BC arrays for a
    // pure box problem with no air).
    const size_t nVx = (size_t)Nx  * Ncz;
    const size_t nVz = (size_t)Ncx * Nz;
    const size_t nP  = (size_t)Ncx * Ncz;
    for (size_t k = 0; k < nVx; ++k) L->act_Vx[k] = 1;
    for (size_t k = 0; k < nVz; ++k) L->act_Vz[k] = 1;
    for (size_t k = 0; k < nP;  ++k) L->act_P [k] = 1;

    // --- Vx mask ---
    // MDOODZ Vx is laid out as Nx * (Nz+1) indexed by k + l*Nx with
    // l in [0, Nz]. Our fine Vx is Nx * Ncz indexed by i + j*Nx with
    // j in [0, Ncz) = [0, Nz-1). The mapping is l = j + 1: our j=0 is
    // MDOODZ's first inner Vx row, our j=Ncz-1 is the last one below
    // the top ghost row.
    if (BCu_type != NULL) {
        for (int j = 0; j < Ncz; ++j) {
            int l = j + 1;  // skip south ghost
            for (int i = 0; i < Nx; ++i) {
                signed char t = BCu_type[(size_t)i + (size_t)l * (size_t)Nx];
                L->act_Vx[MdIdx_Vx(i, j, Nx)] = tag_active(t);
            }
        }
    }

    // --- Vz mask ---
    // MDOODZ Vz layout: (Nx+1) * Nz indexed by k + l*(Nx+1). Our fine Vz
    // is Ncx * Nz indexed by i + j*Ncx. The mapping is k = i + 1 to skip
    // the west ghost column.
    if (BCv_type != NULL) {
        const int nxvz = Nx + 1;
        for (int j = 0; j < Nz; ++j) {
            for (int i = 0; i < Ncx; ++i) {
                int k = i + 1;  // skip west ghost
                signed char t = BCv_type[(size_t)k + (size_t)j * (size_t)nxvz];
                L->act_Vz[MdIdx_Vz(i, j, Ncx)] = tag_active(t);
            }
        }
    }

    // --- P mask ---
    // MDOODZ P layout matches ours (Ncx * Ncz, indexed i + j*Ncx).
    if (BCp_type != NULL) {
        for (int j = 0; j < Ncz; ++j) {
            for (int i = 0; i < Ncx; ++i) {
                signed char t = BCp_type[(size_t)i + (size_t)j * (size_t)Ncx];
                L->act_P[MdIdx_P(i, j, Ncx)] = tag_active(t);
            }
        }
    }
}
