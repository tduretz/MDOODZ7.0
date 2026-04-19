// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2023-2026 MDOODZ Developers team
//
// All rights reserved. This file is part of MDOODZ.
// =========================================================================
//
// Implementation of ApplyStokesOperatorMDOODZ / StokesCellBlockMDOODZ.
//
// The coefficient formulas below are COPIED verbatim (with only mechanical
// edits to drop the Assemble-side AddCoeff3 calls) from
// Xmomentum_InnerNodesDecoupled, Zmomentum_InnerNodesDecoupled and
// Continuity_InnerNodesDecoupled in StokesAssemblyDecoupled.c (Picard
// variant — Phase 1). Any future change to a coefficient formula over
// there MUST be mirrored here, enforced by
// TESTS/StokesMatvecEquivalence.cpp (task 15.9–15.13).
//
// Matvec semantics:
//   row[eqn] = celvol * Σ_{j : active} coeff_j * x_j
// where "active" skips neighbours whose BC type is in {0, 31, 11, 13},
// matching AddCoeff3's behaviour (those values are sent to the Dirichlet
// contribution `bbc`, not into A).
//
// 5×5 Vanka block (task 15.6) re-uses the same coefficient fillers so the
// block builder and the row matvec can never drift from each other: both
// consume `XmomCoeffs` / `ZmomCoeffs` / `ContCoeffs` filled by the same
// helper, they only differ in how the struct's neighbour slots are
// projected (block → 5×5 diagonal sub-matrix + external RHS, matvec →
// single scalar dot product).

#include "StokesAssemblyGMG.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <string.h>

// Skip predicate used to mirror AddCoeff3: neighbours whose BC type is
// 0 (disabled), 31 (Neumann-like, pressure), 11 (Dirichlet velocity), or
// 13 (Neumann velocity) do NOT contribute to A·x — their contribution
// goes to StokesA->bbc during assembly. Periodic (-2) and -12 tags must
// still be honoured; Dirichlet of type 2 (stress BC on a boundary Vx/Vz
// face) is kept as an equation and participates normally.
static inline bool gmg_skip_neighbour(char bc_type) {
    return (bc_type == 0) || (bc_type == 31) || (bc_type == 11) || (bc_type == 13);
}

// Drift canary (task 15.13). A deliberate, compile-flag-gated one-line
// perturbation of a single coefficient. Build with
//   -DMDOODZ_GMG_DRIFT_CANARY
// to activate it; the golden cross-check test in
// TESTS/StokesMatvecEquivalence.cpp MUST then fail. This guarantees the
// test actually catches drift (an always-green test is not a test).
#ifdef MDOODZ_GMG_DRIFT_CANARY
#  define GMG_CANARY_PERTURB(x) ((x) * 1.0000001)
#else
#  define GMG_CANARY_PERTURB(x) (x)
#endif

// ---------------------------------------------------------------------------
// Coefficient bundles. These mirror the (coefficient, index, Set-flag)
// tuples used by the Assemble branch in StokesAssemblyDecoupled.c, but
// without any celvol multiplication — the caller scales up once at the
// end (matvec: `row*celvol`; block builder: callers that want the MDOODZ-
// matrix scaling also multiply, callers that want raw-PDE scaling do
// not). Keeping coefficients celvol-free here lets us route the same
// data to either target cleanly.
// ---------------------------------------------------------------------------

typedef struct {
    // Coefficients (celvol-free).
    double uS, uN, uW, uE, uC;
    double vSW, vSE, vNW, vNE;
    double pW, pE;
    // Full-grid indices these coefficients address.
    int iVxS, iVxW, iVxC, iVxE, iVxN;
    int iVzSW, iVzSE, iVzNW, iVzNE;
    int iPrW, iPrE;
    // Set flags — identical to the Assemble branch's SetVx*/SetVz*/SetPr*
    // gates. `false` means "this neighbour is skipped entirely"
    // (neighbouring cell is tag 30, or we're on the outer ring without a
    // periodic wrap). `true` means "this neighbour contributes unless
    // further suppressed by `gmg_skip_neighbour` on its BC tag".
    bool SetVxS, SetVxW, SetVxC, SetVxE, SetVxN;
    bool SetVzSW, SetVzSE, SetVzNW, SetVzNE;
    bool SetPrW, SetPrE;
} XmomCoeffs;

typedef struct {
    double vS, vN, vW, vE, vC;
    double uSW, uSE, uNW, uNE;
    double pS, pN;
    int iVzS, iVzW, iVzC, iVzE, iVzN;
    int iVxSW, iVxSE, iVxNW, iVxNE;
    int iPrS, iPrN;
    bool SetVzS, SetVzW, SetVzC, SetVzE, SetVzN;
    bool SetVxSW, SetVxSE, SetVxNW, SetVxNE;
    bool SetPrS, SetPrN;
} ZmomCoeffs;

typedef struct {
    // Continuity coefficients — all four neighbours are inside the block
    // for the cell centred at the current P; we keep the struct shape
    // homogeneous with the Xmom/Zmom ones so the block builder can treat
    // them uniformly.
    double uW, uE, vS, vN;
    int iVxW, iVxE, iVzS, iVzN;
    bool SetVxW, SetVxE, SetVzS, SetVzN;
} ContCoeffs;

// ---------------------------------------------------------------------------
// X-momentum coefficient extractor. Mirrors Xmomentum_InnerNodesDecoupled
// (Picard branch) line-for-line except that (a) coefficients are stored
// in a bundle instead of fed into AddCoeff3, and (b) no celvol.
// ---------------------------------------------------------------------------
static void xmom_fill_coeffs(grid *mesh, params model, int ix, int iz,
                             int stab, int comp, double om,
                             double one_dx, double one_dz,
                             int c1, int c2, int c3, int nx, int nz,
                             XmomCoeffs *C) {
    const int nxvz = nx + 1, ncx = nx - 1;
    (void)iz;
    const double dx = mesh->dx, dz = mesh->dz;
    double OOP = 1.0, uC_corr = 0.0;
    if (model.out_of_plane == 1) OOP = 3.0 / 2.0;
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    double SxxBCW = 0., SxxBCE = 0.;
    if (ix == 0    && BCu_t[c1] == 2) SxxBCW = 1.0;
    if (ix == nx-1 && BCu_t[c1] == 2) SxxBCE = 1.0;

    int iVxC  = c1;
    int iVxW  = iVxC - 1;
    int iVxE  = iVxC + 1;
    int iVxS  = iVxC - nx;
    int iVxN  = iVxC + nx;
    int iVzSW = c3 - nxvz;
    int iVzSE = c3 - nxvz + 1;
    int iVzNW = c3;
    int iVzNE = c3 + 1;
    int iPrW  = c2;
    int iPrE  = c2 + 1;
    int ixyN  = c1;
    int ixyS  = c1 - nx;

    if (BCu_t[iVxC] == -2) {
        iVxW  = c1 + nx - 2;
        iPrW  = c2 + ncx;
        iVzSW = c3 - 2;
        iVzNW = c3 + nxvz - 2;
    }

    const double D11E = mesh->D11_n[iPrE];
    const double D11W = mesh->D11_n[iPrW];
    const double D33N = mesh->D33_s[ixyN];
    const double D33S = mesh->D33_s[ixyS];

    double inE = 0.0, inW = 0.0;
    if (BCp_t[iPrW] == -1) inW = 1.0;
    if (BCp_t[iPrE] == -1) inE = 1.0;

    double inS = 1.0, inN = 1.0, RegS = 1.0, RegN = 1.0;
    double NeuS = 0.0, NeuN = 0.0, DirS = 0.0, DirN = 0.0;
    if (BCg_t[ixyS] == 30) inS = 0.0;
    if (BCg_t[ixyN] == 30) inN = 0.0;
    if (BCu_t[iVxS] == 13) { NeuS = 1.0; RegS = 0.0; }
    if (BCu_t[iVxN] == 13) { NeuN = 1.0; RegN = 0.0; }
    if (BCu_t[iVxS] == 11) { DirS = 1.0; RegS = 0.0; }
    if (BCu_t[iVxN] == 11) { DirN = 1.0; RegN = 0.0; }
    (void)NeuS; (void)NeuN;

    double DirSW = 0.0, DirNW = 0.0;
    if (BCv_t[c3 - nxvz] == 11) DirSW = 1.0;
    if (BCv_t[c3]        == 11) DirNW = 1.0;
    double DirSE = 0.0, DirNE = 0.0;
    if (BCv_t[c3 - nxvz + 1] == 11) DirSE = 1.0;
    if (BCv_t[c3 + 1]        == 11) DirNE = 1.0;

    double uW = GMG_CANARY_PERTURB((1.0/3.0)*D11W*inW*(1 - SxxBCW)*(OOP*comp*inW - 3)/pow(dx, 2));
    double uC = -1.0/6.0*(SxxBCE*(2*D11W*pow(dz, 2)*inW*(OOP*comp*inW - 3) - 3*pow(dx, 2)*(D33N*inN*(2*DirN + RegN) + D33S*inS*(2*DirS + RegS))) + SxxBCW*(2*D11E*pow(dz, 2)*inE*(OOP*comp*inE - 3) - 3*pow(dx, 2)*(D33N*inN*(2*DirN + RegN) + D33S*inS*(2*DirS + RegS))) + 2*(3*pow(dx, 2)*(D33N*inN*(2*DirN + RegN) + D33S*inS*(2*DirS + RegS)) - pow(dz, 2)*(D11E*inE*(OOP*comp*inE - 3) + D11W*inW*(OOP*comp*inW - 3)))*(SxxBCE + SxxBCW - 1))/(pow(dx, 2)*pow(dz, 2));
    double uE = (1.0/3.0)*D11E*inE*(1 - SxxBCE)*(OOP*comp*inE - 3)/pow(dx, 2);
    double uS = (1.0/2.0)*D33S*RegS*inS*(SxxBCE + SxxBCW - 2)/pow(dz, 2);
    double uN = (1.0/2.0)*D33N*RegN*inN*(SxxBCE + SxxBCW - 2)/pow(dz, 2);
    double vSW = (1.0/6.0)*(-3*D33S*SxxBCW*inS*(2*DirSE*SxxBCE - SxxBCE - SxxBCW + 1) + SxxBCE*(2*D11W*OOP*comp*pow(inW, 2) - 3*D33S*inS*(2*DirSE*SxxBCE - SxxBCE - SxxBCW + 1)) - 2*(D11W*OOP*comp*pow(inW, 2) - 3*D33S*inS*(2*DirSE*SxxBCE - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    double vSE = (1.0/6.0)*(3*D33S*SxxBCE*inS*(2*DirSW*SxxBCW - SxxBCE - SxxBCW + 1) - SxxBCW*(2*D11E*OOP*comp*pow(inE, 2) - 3*D33S*inS*(2*DirSW*SxxBCW - SxxBCE - SxxBCW + 1)) + 2*(D11E*OOP*comp*pow(inE, 2) - 3*D33S*inS*(2*DirSW*SxxBCW - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    double vNW = (1.0/6.0)*(3*D33N*SxxBCW*inN*(2*DirNE*SxxBCE - SxxBCE - SxxBCW + 1) - SxxBCE*(2*D11W*OOP*comp*pow(inW, 2) - 3*D33N*inN*(2*DirNE*SxxBCE - SxxBCE - SxxBCW + 1)) + 2*(D11W*OOP*comp*pow(inW, 2) - 3*D33N*inN*(2*DirNE*SxxBCE - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    double vNE = (1.0/6.0)*(-3*D33N*SxxBCE*inN*(2*DirNW*SxxBCW - SxxBCE - SxxBCW + 1) + SxxBCW*(2*D11E*OOP*comp*pow(inE, 2) - 3*D33N*inN*(2*DirNW*SxxBCW - SxxBCE - SxxBCW + 1)) - 2*(D11E*OOP*comp*pow(inE, 2) - 3*D33N*inN*(2*DirNW*SxxBCW - SxxBCE - SxxBCW + 1))*(SxxBCE + SxxBCW - 1))/(dx*dz);
    double pW = (SxxBCW - 1)/dx;
    double pE = (1 - SxxBCE)/dx;

    if (stab == 1) {
        double drhodx = (mesh->rho_n[c2 + 1] - mesh->rho_n[c2]) * one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        if (uC + uC_corr < 0.0) uC_corr = 0.0;
        uC += uC_corr;
    }
    (void)one_dz;
    (void)nz;

    const bool SetVxS  = BCu_t[iVxS] != 30;
    const bool SetVxW  = BCu_t[iVxW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVxC  = BCu_t[iVxC] != 30;
    const bool SetVxE  = BCu_t[iVxE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetVxN  = BCu_t[iVxN] != 30;
    const bool SetVzSW = BCv_t[iVzSW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVzSE = BCv_t[iVzSE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetVzNW = BCv_t[iVzNW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVzNE = BCv_t[iVzNE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetPrW  = BCp_t[iPrW]  != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetPrE  = BCp_t[iPrE]  != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);

    C->uS = uS; C->uN = uN; C->uW = uW; C->uE = uE; C->uC = uC;
    C->vSW = vSW; C->vSE = vSE; C->vNW = vNW; C->vNE = vNE;
    C->pW = pW; C->pE = pE;
    C->iVxS = iVxS; C->iVxW = iVxW; C->iVxC = iVxC; C->iVxE = iVxE; C->iVxN = iVxN;
    C->iVzSW = iVzSW; C->iVzSE = iVzSE; C->iVzNW = iVzNW; C->iVzNE = iVzNE;
    C->iPrW = iPrW; C->iPrE = iPrE;
    C->SetVxS = SetVxS; C->SetVxW = SetVxW; C->SetVxC = SetVxC;
    C->SetVxE = SetVxE; C->SetVxN = SetVxN;
    C->SetVzSW = SetVzSW; C->SetVzSE = SetVzSE;
    C->SetVzNW = SetVzNW; C->SetVzNE = SetVzNE;
    C->SetPrW = SetPrW; C->SetPrE = SetPrE;
}

static double xmom_row_matvec(grid *mesh, params model, int ix, int iz,
                              int stab, int comp, double om,
                              double one_dx, double one_dz, double celvol,
                              int c1, int c2, int c3, int nx, int nz,
                              const double *u, const double *v, const double *p) {
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCp_t = mesh->BCp.type;
    XmomCoeffs C;
    xmom_fill_coeffs(mesh, model, ix, iz, stab, comp, om, one_dx, one_dz,
                     c1, c2, c3, nx, nz, &C);

    double row = 0.0;
    if (C.SetVxS  && !gmg_skip_neighbour(BCu_t[C.iVxS]))  row += C.uS  * u[C.iVxS];
    if (C.SetVxW  && !gmg_skip_neighbour(BCu_t[C.iVxW]))  row += C.uW  * u[C.iVxW];
    if (C.SetVxC  && !gmg_skip_neighbour(BCu_t[C.iVxC]))  row += C.uC  * u[C.iVxC];
    if (C.SetVxE  && !gmg_skip_neighbour(BCu_t[C.iVxE]))  row += C.uE  * u[C.iVxE];
    if (C.SetVxN  && !gmg_skip_neighbour(BCu_t[C.iVxN]))  row += C.uN  * u[C.iVxN];
    if (C.SetVzSW && !gmg_skip_neighbour(BCv_t[C.iVzSW])) row += C.vSW * v[C.iVzSW];
    if (C.SetVzSE && !gmg_skip_neighbour(BCv_t[C.iVzSE])) row += C.vSE * v[C.iVzSE];
    if (C.SetVzNW && !gmg_skip_neighbour(BCv_t[C.iVzNW])) row += C.vNW * v[C.iVzNW];
    if (C.SetVzNE && !gmg_skip_neighbour(BCv_t[C.iVzNE])) row += C.vNE * v[C.iVzNE];
    if (C.SetPrW  && !gmg_skip_neighbour(BCp_t[C.iPrW]))  row += C.pW  * p[C.iPrW];
    if (C.SetPrE  && !gmg_skip_neighbour(BCp_t[C.iPrE]))  row += C.pE  * p[C.iPrE];

    return row * celvol;
}

// ---------------------------------------------------------------------------
// Newton X-momentum Jacobian (task 15.19 / design D11 Phase 2).
//
// Ports `Xjacobian_InnerNodesDecoupled3` (StokesAssemblyDecoupled.c:219)
// line-for-line, dropping the AddCoeff3 assembly calls and celvol scaling.
// This is the full Newton Jacobian row of the x-momentum equation — a
// 9-Vx / 12-Vz / 6-P stencil (vs. Picard's 5-Vx / 4-Vz / 2-P) with the
// anisotropic D-tensor (D11..D14 on normal-stress cells, D31..D34 on
// shear-stress vertices) coupling every component to every neighbour.
//
// Correctness contract: for every Newton-mode configuration, this
// matvec SHALL agree with `BuildJacobianOperatorDecoupled`'s assembled
// sparse Jacobian to within 1e-12 relative tolerance. Enforced by the
// `JacobianMatvecMatches*` family in `TESTS/StokesMatvecEquivalence.cpp`.
// ---------------------------------------------------------------------------

typedef struct {
    double uC, uS, uN, uW, uE;
    double uSW, uSE, uNW, uNE;
    double vSW, vSE, vNW, vNE;
    double vSSW, vSSE, vNNW, vNNE;
    double vSWW, vSEE, vNWW, vNEE;
    double pW, pE;
    double pSW, pSE, pNW, pNE;
    int iVxC, iVxS, iVxN, iVxW, iVxE;
    int iVxSW, iVxSE, iVxNW, iVxNE;
    int iVzSW, iVzSE, iVzNW, iVzNE;
    int iVzSSW, iVzSSE, iVzNNW, iVzNNE;
    int iVzSWW, iVzSEE, iVzNWW, iVzNEE;
    int iPrW, iPrE;
    int iPrSW, iPrSE, iPrNW, iPrNE;
    bool SetVxS, SetVxW, SetVxC, SetVxE, SetVxN;
    bool SetVxSW, SetVxSE, SetVxNW, SetVxNE;
    bool SetVzSW, SetVzSE, SetVzNW, SetVzNE;
    bool SetVzSSW, SetVzSSE, SetVzNNW, SetVzNNE;
    bool SetVzSWW, SetVzSEE, SetVzNWW, SetVzNEE;
    bool SetPrW, SetPrE;
    bool SetPrSW, SetPrSE, SetPrNW, SetPrNE;
} XmomNewtonCoeffs;

static void xmom_fill_coeffs_newton(grid *mesh, params model, int ix, int iz,
                                    int stab, int comp, double om,
                                    double one_dx, double one_dz,
                                    int c1, int c2, int c3, int nx, int nz,
                                    XmomNewtonCoeffs *C) {
    // `Xjacobian_InnerNodesDecoupled3` uses `ix == k` and `iz == l`.
    const int k = ix, l = iz;
    const int nxvz = nx + 1, ncx = nx - 1, nzvx = nz + 1;
    const double dx = mesh->dx, dz = mesh->dz;
    double out_of_plane = 1.0, uC_corr = 0.0;
    if (model.out_of_plane == 1) out_of_plane = 3.0 / 2.0;
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    double SxxBCW = 0., SxxBCE = 0.;
    if (ix == 0      && BCu_t[c1] == 2) SxxBCW = 1.0;
    if (ix == nx - 1 && BCu_t[c1] == 2) SxxBCE = 1.0;

    int iVxC   = c1;
    int iVxW   = iVxC - 1;
    int iVxE   = iVxC + 1;
    int iVxS   = iVxC - nx;
    int iVxN   = iVxC + nx;
    int iVxSW  = iVxS - 1;
    int iVxSE  = iVxS + 1;
    int iVxNW  = iVxN - 1;
    int iVxNE  = iVxN + 1;
    int iVzSW  = c3 - nxvz;
    int iVzSE  = c3 - nxvz + 1;
    int iVzNW  = c3;
    int iVzNE  = c3 + 1;
    int iVzSWW = iVzSW - 1;
    int iVzSEE = iVzSE + 1;
    int iVzNWW = iVzNW - 1;
    int iVzNEE = iVzNE + 1;
    int iVzNNW = iVzNW + nxvz, iVzSSW = iVzSW - nxvz;
    int iVzNNE = iVzNE + nxvz, iVzSSE = iVzSE - nxvz;
    int iPrW   = c2;
    int iPrE   = c2 + 1;
    int ixyN   = c1;
    int ixyS   = c1 - nx;
    int ixySW  = ixyS - 1;
    int ixySE  = ixyS + 1;
    int ixyNW  = ixyN - 1;
    int ixyNE  = ixyN + 1;
    int iPrSW  = iPrW - ncx;
    int iPrSE  = iPrE - ncx;
    int iPrNW  = iPrW + ncx;
    int iPrNE  = iPrE + ncx;

    if (l > 1)       if (BCv_t[iVzSWW] == -12) iVzSWW = iVzSW + (nx - 2);
    if (l < nz - 1)  if (BCv_t[iVzNWW] == -12) iVzNWW = iVzNW + (nx - 2);
    if (l > 1)       if (BCv_t[iVzSEE] == -12) iVzSEE = iVzSE - (nx - 2);
    if (l < nz - 1)  if (BCv_t[iVzNEE] == -12) iVzNEE = iVzNE - (nx - 2);

    if (BCu_t[iVxC] == -2) {
        iVxW   = c1 + nx - 2;      iPrW = c2 + ncx; iVzSW = c3 - 2; iVzNW = c3 + nxvz - 2;
        iPrNW  = iPrW + ncx;       iPrSW = iPrW - ncx;
        iVxSW  = iVxW - nx;        iVxNW = iVxW + nx;
        iVzNNW = iVzNW + nxvz;     iVzSSW = iVzSW - nxvz;
        iVzSWW = iVzSW - 1;        iVzNWW = iVzNW - 1;
        ixySW  = ixyS + nx;
        ixyNW  = ixyN + nx;
    }

    if (k == 0)      { ixySW = ixyS;  ixyNW = ixyN; }
    if (k == nx - 1) { ixySE = ixyS;  ixyNE = ixyN; }

    // D-tensor entries (four cell centres + two vertices).
    double D11E = mesh->D11_n[iPrE];
    double D12E = mesh->D12_n[iPrE];
    double D13E = mesh->D13_n[iPrE];
    double D14E = mesh->D14_n[iPrE];

    double D11W = mesh->D11_n[iPrW];
    double D12W = mesh->D12_n[iPrW];
    double D13W = mesh->D13_n[iPrW];
    double D14W = mesh->D14_n[iPrW];

    double D31N = mesh->D31_s[ixyN];
    double D32N = mesh->D32_s[ixyN];
    double D33N = mesh->D33_s[ixyN];
    double D34N = mesh->D34_s[ixyN];

    double D31S = mesh->D31_s[ixyS];
    double D32S = mesh->D32_s[ixyS];
    double D33S = mesh->D33_s[ixyS];
    double D34S = mesh->D34_s[ixyS];

    double inE = 0.0, inW = 0.0, inS = 0.0, inN = 0.0, inSv = 0.0, inNv = 0.0;
    if (BCp_t[iPrW] == -1) inW = 1.0;
    if (BCp_t[iPrE] == -1) inE = 1.0;
    if (BCu_t[iVxS] != 30 && BCu_t[iVxS] != 13) inS = 1.0;
    if (BCu_t[iVxN] != 30 && BCu_t[iVxN] != 13) inN = 1.0;
    if (BCg_t[ixyS] != 30) inSv = 1.0;
    if (BCg_t[ixyN] != 30) inNv = 1.0;

    double wE = 0.0, wW = 0.0, wS = 0.0, wN = 0.0;
    double inSWc = 0.0, inSEc = 0.0, inNWc = 0.0, inNEc = 0.0;
    double inSWv = 0.0, inSEv = 0.0, inNWv = 0.0, inNEv = 0.0;

    if (l > 1) {
        if (BCp_t[iPrSW] == -1) inSWc = 1.0;
        if (BCp_t[iPrSE] == -1) inSEc = 1.0;
    }
    if (l < nzvx - 2) {
        if (BCp_t[iPrNW] == -1) inNWc = 1.0;
        if (BCp_t[iPrNE] == -1) inNEc = 1.0;
    }
    if ((k > 0) || (k == 0 && (BCv_t[iVzSW] == -1 || BCv_t[iVzSW] == 0))) {
        if (BCg_t[ixySW] != 30) inSWv = 1.0;
        if (BCg_t[ixyNW] != 30) inNWv = 1.0;
    }
    if ((k < nx - 1) || (k == nx - 1 && (BCv_t[iVzSE] == -1 || BCv_t[iVzSE] == 0))) {
        if (BCg_t[ixySE] != 30) inSEv = 1.0;
        if (BCg_t[ixyNE] != 30) inNEv = 1.0;
    }

    double inSWW = 0.0, inNWW = 0.0;
    if (BCv_t[iVzSWW] == -1) inSWW = 1.0;
    if (BCv_t[iVzNWW] == -1) inNWW = 1.0;
    double inSEE = 0.0, inNEE = 0.0;
    if (BCv_t[iVzSEE] == -1) inSEE = 1.0;
    if (BCv_t[iVzNEE] == -1) inNEE = 1.0;

    wE = inN + inS + inNEv + inSEv;
    wW = inN + inS + inNWv + inSWv;
    wS = inW + inE + inSWc + inSEc;
    wN = inW + inE + inNWc + inNEc;
    if (wW > 1.0) wW = 1.0 / wW;
    if (wE > 1.0) wE = 1.0 / wE;
    if (wS > 1.0) wS = 1.0 / wS;
    if (wN > 1.0) wN = 1.0 / wN;

    // Coefficient formulas — verbatim from Xjacobian_InnerNodesDecoupled3.
    double uW = (1.0/3.0)*(-0.75*D13W*dx*(inN*inNWv - inS*inSWv) + dx*inW*(-D31N*wN*(comp*out_of_plane - 3) + D31S*wS*(comp*out_of_plane - 3) - D32N*comp*out_of_plane*wN + D32S*comp*out_of_plane*wS) + dz*inW*(D11W*(comp*out_of_plane - 3) + D12W*comp*out_of_plane))/(pow(dx, 2)*dz);
    double uC = (1.0/3.0)*(dx*(3*dx*(D33N*pow(inN, 2)*inNv + D33S*pow(inS, 2)*inSv) + dz*(inE - inW)*(-D31N*wN*(comp*out_of_plane - 3) + D31S*wS*(comp*out_of_plane - 3) - D32N*comp*out_of_plane*wN + D32S*comp*out_of_plane*wS)) - dz*(0.75*dx*(-D13E + D13W)*(pow(inN, 2)*inNv - pow(inS, 2)*inSv) + dz*(D11E*inE*(comp*out_of_plane - 3) + D11W*inW*(comp*out_of_plane - 3) + D12E*comp*inE*out_of_plane + D12W*comp*inW*out_of_plane)))/(pow(dx, 2)*pow(dz, 2));
    double uE = (1.0/3.0)*(0.75*D13E*dx*(inN*inNEv - inS*inSEv) + dx*inE*(D31N*wN*(comp*out_of_plane - 3) - D31S*wS*(comp*out_of_plane - 3) + D32N*comp*out_of_plane*wN - D32S*comp*out_of_plane*wS) + dz*inE*(D11E*(comp*out_of_plane - 3) + D12E*comp*out_of_plane))/(pow(dx, 2)*dz);
    double uS = (-D33S*dx*pow(inS, 2)*inSv + 0.25*dz*pow(inS, 2)*inSv*(D13E - D13W) + (1.0/3.0)*dz*wS*(inSEc - inSWc)*(D31S*(comp*out_of_plane - 3) + D32S*comp*out_of_plane))/(dx*pow(dz, 2));
    double uN = (-D33N*dx*pow(inN, 2)*inNv + 0.25*dz*pow(inN, 2)*inNv*(-D13E + D13W) - 1.0/3.0*dz*wN*(inNEc - inNWc)*(D31N*(comp*out_of_plane - 3) + D32N*comp*out_of_plane))/(dx*pow(dz, 2));
    double vSW = (1.0/3.0)*(-dx*(3*D33S*dz*inS*inSv + dx*(D31N*comp*inW*out_of_plane*wN + D31S*comp*out_of_plane*wS*(inSWc - inW) + D32N*inW*wN*(comp*out_of_plane - 3) + D32S*wS*(inSWc - inW)*(comp*out_of_plane - 3))) + dz*(dx*inW*(D11W*comp*out_of_plane + D12W*(comp*out_of_plane - 3)) + 0.75*dz*(D13E*inS*inSv - D13W*(inS*inSv - inSWW*inSWv))))/(pow(dx, 2)*pow(dz, 2));
    double vSE = (1.0/3.0)*(dx*(3*D33S*dz*inS*inSv + dx*(-D31N*comp*inE*out_of_plane*wN + D31S*comp*out_of_plane*wS*(inE - inSEc) - D32N*inE*wN*(comp*out_of_plane - 3) + D32S*wS*(inE - inSEc)*(comp*out_of_plane - 3))) - dz*(dx*inE*(D11E*comp*out_of_plane + D12E*(comp*out_of_plane - 3)) + 0.75*dz*(D13E*(inS*inSv - inSEE*inSEv) - D13W*inS*inSv)))/(pow(dx, 2)*pow(dz, 2));
    double vNW = (1.0/3.0)*(dx*(3*D33N*dz*inN*inNv - dx*(D31N*comp*out_of_plane*wN*(inNWc - inW) + D31S*comp*inW*out_of_plane*wS + D32N*wN*(inNWc - inW)*(comp*out_of_plane - 3) + D32S*inW*wS*(comp*out_of_plane - 3))) - dz*(dx*inW*(D11W*comp*out_of_plane + D12W*(comp*out_of_plane - 3)) + 0.75*dz*(-D13E*inN*inNv + D13W*(inN*inNv - inNWW*inNWv))))/(pow(dx, 2)*pow(dz, 2));
    double vNE = (1.0/3.0)*(-dx*(3*D33N*dz*inN*inNv + dx*(-D31N*comp*out_of_plane*wN*(inE - inNEc) + D31S*comp*inE*out_of_plane*wS - D32N*wN*(inE - inNEc)*(comp*out_of_plane - 3) + D32S*inE*wS*(comp*out_of_plane - 3))) + dz*(dx*inE*(D11E*comp*out_of_plane + D12E*(comp*out_of_plane - 3)) + 0.75*dz*(-D13E*(inN*inNv - inNEE*inNEv) + D13W*inN*inNv)))/(pow(dx, 2)*pow(dz, 2));
    double uSW = (1.0/3.0)*(-0.75*D13W*inS*inSWv + D31S*inSWc*wS*(comp*out_of_plane - 3) + D32S*comp*inSWc*out_of_plane*wS)/(dx*dz);
    double uSE = (1.0/3.0)*(0.75*D13E*inS*inSEv - D31S*inSEc*wS*(comp*out_of_plane - 3) - D32S*comp*inSEc*out_of_plane*wS)/(dx*dz);
    double uNW = (1.0/3.0)*(0.75*D13W*inN*inNWv - D31N*inNWc*wN*(comp*out_of_plane - 3) - D32N*comp*inNWc*out_of_plane*wN)/(dx*dz);
    double uNE = (1.0/3.0)*(-0.75*D13E*inN*inNEv + D31N*inNEc*wN*(comp*out_of_plane - 3) + D32N*comp*inNEc*out_of_plane*wN)/(dx*dz);
    double vSWW = -0.25*D13W*inSWW*inSWv/pow(dx, 2);
    double vSEE = -0.25*D13E*inSEE*inSEv/pow(dx, 2);
    double vNWW = -0.25*D13W*inNWW*inNWv/pow(dx, 2);
    double vNEE = -0.25*D13E*inNEE*inNEv/pow(dx, 2);
    double vSSW = (1.0/3.0)*inSWc*wS*(D31S*comp*out_of_plane + D32S*(comp*out_of_plane - 3))/pow(dz, 2);
    double vSSE = (1.0/3.0)*inSEc*wS*(D31S*comp*out_of_plane + D32S*(comp*out_of_plane - 3))/pow(dz, 2);
    double vNNW = (1.0/3.0)*inNWc*wN*(D31N*comp*out_of_plane + D32N*(comp*out_of_plane - 3))/pow(dz, 2);
    double vNNE = (1.0/3.0)*inNEc*wN*(D31N*comp*out_of_plane + D32N*(comp*out_of_plane - 3))/pow(dz, 2);
    double pW   = -inW*one_dx + inW*(D14W*dz + dx*(-D34N*inN*wN + D34S*inS*wS))/(dx*dz);
    double pE   =  inE*one_dx + inE*(-D14E*dz + dx*(-D34N*inN*wN + D34S*inS*wS))/(dx*dz);
    double pSW  =  D34S*inS*inSWc*wS/dz;
    double pSE  =  D34S*inS*inSEc*wS/dz;
    double pNW  = -D34N*inN*inNWc*wN/dz;
    double pNE  = -D34N*inN*inNEc*wN/dz;

    // Stabilisation with density gradients (same as Picard path).
    if (stab == 1) {
        double drhodx = (mesh->rho_n[c2 + 1] - mesh->rho_n[c2]) * one_dx;
        uC_corr = 1.00 * om * model.dt * mesh->gx[c1] * drhodx;
        if (uC + uC_corr > 0.0) uC += uC_corr;
    }

    // Non-conforming Dirichlet corrections — mirror
    // Xjacobian_InnerNodesDecoupled3 lines 420-429 exactly. These fold
    // the contributions of Dirichlet-neighbour stencil nodes into the
    // diagonal so that `bbc` sees their value as a fixed state rather
    // than as a row entry, preserving the same matvec output as the
    // assembled matrix.
    if (BCu_t[iVxS]   == 11) uC  -=  uS;
    if (BCu_t[iVxN]   == 11) uC  -=  uN;
    if (BCu_t[iVxSW]  == 11) uW  -=  uSW;
    if (BCu_t[iVxSE]  == 11) uE  -=  uSE;
    if (BCu_t[iVxNW]  == 11) uW  -=  uNW;
    if (BCu_t[iVxNE]  == 11) uE  -=  uNE;
    if (BCv_t[iVzNWW] == 11) vNW -= vNWW;
    if (BCv_t[iVzNEE] == 11) vNE -= vNEE;
    if (BCv_t[iVzSWW] == 11) vSW -= vSWW;
    if (BCv_t[iVzSEE] == 11) vSE -= vSEE;

    (void)one_dz; (void)nz;

    // Set flags — mirror Xjacobian_InnerNodesDecoupled3's Assemble-side
    // SetVx/SetVz/SetPr gates verbatim. Two extra sets of guards (the
    // `Newton==1 && l>1` / `Newton==1 && l<nz-1` pressure-corner
    // gates, and the `(k > 1 || (k∈{0,1} && -1))` / symmetric east/
    // north patterns for the double-offset Vz neighbours) carry over
    // the out-of-stencil handling that `Xjacobian3` applies.
    const bool SetVxSW = BCu_t[iVxSW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVxS  = BCu_t[iVxS]  != 30;
    const bool SetVxSE = BCu_t[iVxSE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetVxW  = BCu_t[iVxW]  != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVxC  = BCu_t[iVxC]  != 30;
    const bool SetVxE  = BCu_t[iVxE]  != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetVxNW = BCu_t[iVxNW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVxN  = BCu_t[iVxN]  != 30;
    const bool SetVxNE = BCu_t[iVxNE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);

    const bool SetVzSW = BCv_t[iVzSW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVzSE = BCv_t[iVzSE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetVzNW = BCv_t[iVzNW] != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetVzNE = BCv_t[iVzNE] != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);

    const bool SetVzSSW = (l > 1)      && (BCv_t[iVzSSW] != 30);
    const bool SetVzSSE = (l > 1)      && (BCv_t[iVzSSE] != 30);
    const bool SetVzNNW = (l < nz - 1) && (BCv_t[iVzNNW] != 30);
    const bool SetVzNNE = (l < nz - 1) && (BCv_t[iVzNNE] != 30);

    const bool SetVzSWW = ((k > 1) || ((k == 0 || k == 1) && BCv_t[iVzSWW] == -1)) && BCv_t[iVzSWW] != 30;
    const bool SetVzSEE = ((k < nx - 2) || ((k == nx - 2 || k == nx - 1) && BCv_t[iVzSEE] == -1)) && BCv_t[iVzSEE] != 30;
    const bool SetVzNWW = ((k > 1) || ((k == 0 || k == 1) && BCv_t[iVzNWW] == -1)) && BCv_t[iVzNWW] != 30;
    const bool SetVzNEE = ((k < nx - 2) || ((k == nx - 2 || k == nx - 1) && BCv_t[iVzNEE] == -1)) && BCv_t[iVzNEE] != 30;

    const bool SetPrW  = BCp_t[iPrW]  != 30 && (ix > 0      || BCu_t[iVxC] == -2);
    const bool SetPrE  = BCp_t[iPrE]  != 30 && (ix < nx - 1 || BCu_t[iVxC] == -2);
    const bool SetPrSW = (l > 1)      && (BCp_t[iPrSW] != 30);
    const bool SetPrSE = (l > 1)      && (BCp_t[iPrSE] != 30);
    const bool SetPrNW = (l < nz - 1) && (BCp_t[iPrNW] != 30);
    const bool SetPrNE = (l < nz - 1) && (BCp_t[iPrNE] != 30);

    // Pack.
    C->uC = uC; C->uS = uS; C->uN = uN; C->uW = uW; C->uE = uE;
    C->uSW = uSW; C->uSE = uSE; C->uNW = uNW; C->uNE = uNE;
    C->vSW = vSW; C->vSE = vSE; C->vNW = vNW; C->vNE = vNE;
    C->vSSW = vSSW; C->vSSE = vSSE; C->vNNW = vNNW; C->vNNE = vNNE;
    C->vSWW = vSWW; C->vSEE = vSEE; C->vNWW = vNWW; C->vNEE = vNEE;
    C->pW = pW; C->pE = pE;
    C->pSW = pSW; C->pSE = pSE; C->pNW = pNW; C->pNE = pNE;

    C->iVxC = iVxC; C->iVxS = iVxS; C->iVxN = iVxN; C->iVxW = iVxW; C->iVxE = iVxE;
    C->iVxSW = iVxSW; C->iVxSE = iVxSE; C->iVxNW = iVxNW; C->iVxNE = iVxNE;
    C->iVzSW = iVzSW; C->iVzSE = iVzSE; C->iVzNW = iVzNW; C->iVzNE = iVzNE;
    C->iVzSSW = iVzSSW; C->iVzSSE = iVzSSE; C->iVzNNW = iVzNNW; C->iVzNNE = iVzNNE;
    C->iVzSWW = iVzSWW; C->iVzSEE = iVzSEE; C->iVzNWW = iVzNWW; C->iVzNEE = iVzNEE;
    C->iPrW = iPrW; C->iPrE = iPrE;
    C->iPrSW = iPrSW; C->iPrSE = iPrSE; C->iPrNW = iPrNW; C->iPrNE = iPrNE;

    C->SetVxS = SetVxS; C->SetVxW = SetVxW; C->SetVxC = SetVxC;
    C->SetVxE = SetVxE; C->SetVxN = SetVxN;
    C->SetVxSW = SetVxSW; C->SetVxSE = SetVxSE;
    C->SetVxNW = SetVxNW; C->SetVxNE = SetVxNE;
    C->SetVzSW = SetVzSW; C->SetVzSE = SetVzSE;
    C->SetVzNW = SetVzNW; C->SetVzNE = SetVzNE;
    C->SetVzSSW = SetVzSSW; C->SetVzSSE = SetVzSSE;
    C->SetVzNNW = SetVzNNW; C->SetVzNNE = SetVzNNE;
    C->SetVzSWW = SetVzSWW; C->SetVzSEE = SetVzSEE;
    C->SetVzNWW = SetVzNWW; C->SetVzNEE = SetVzNEE;
    C->SetPrW = SetPrW; C->SetPrE = SetPrE;
    C->SetPrSW = SetPrSW; C->SetPrSE = SetPrSE;
    C->SetPrNW = SetPrNW; C->SetPrNE = SetPrNE;
}

static double xmom_row_matvec_newton(grid *mesh, params model, int ix, int iz,
                                     int stab, int comp, double om,
                                     double one_dx, double one_dz, double celvol,
                                     int c1, int c2, int c3, int nx, int nz,
                                     const double *u, const double *v, const double *p) {
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCp_t = mesh->BCp.type;
    XmomNewtonCoeffs C;
    xmom_fill_coeffs_newton(mesh, model, ix, iz, stab, comp, om, one_dx, one_dz,
                            c1, c2, c3, nx, nz, &C);

    double row = 0.0;
    // 9 Vx neighbours.
    if (C.SetVxSW && !gmg_skip_neighbour(BCu_t[C.iVxSW])) row += C.uSW * u[C.iVxSW];
    if (C.SetVxS  && !gmg_skip_neighbour(BCu_t[C.iVxS]))  row += C.uS  * u[C.iVxS];
    if (C.SetVxSE && !gmg_skip_neighbour(BCu_t[C.iVxSE])) row += C.uSE * u[C.iVxSE];
    if (C.SetVxW  && !gmg_skip_neighbour(BCu_t[C.iVxW]))  row += C.uW  * u[C.iVxW];
    if (C.SetVxC  && !gmg_skip_neighbour(BCu_t[C.iVxC]))  row += C.uC  * u[C.iVxC];
    if (C.SetVxE  && !gmg_skip_neighbour(BCu_t[C.iVxE]))  row += C.uE  * u[C.iVxE];
    if (C.SetVxNW && !gmg_skip_neighbour(BCu_t[C.iVxNW])) row += C.uNW * u[C.iVxNW];
    if (C.SetVxN  && !gmg_skip_neighbour(BCu_t[C.iVxN]))  row += C.uN  * u[C.iVxN];
    if (C.SetVxNE && !gmg_skip_neighbour(BCu_t[C.iVxNE])) row += C.uNE * u[C.iVxNE];
    // 12 Vz neighbours.
    if (C.SetVzSSW && !gmg_skip_neighbour(BCv_t[C.iVzSSW])) row += C.vSSW * v[C.iVzSSW];
    if (C.SetVzSSE && !gmg_skip_neighbour(BCv_t[C.iVzSSE])) row += C.vSSE * v[C.iVzSSE];
    if (C.SetVzSWW && !gmg_skip_neighbour(BCv_t[C.iVzSWW])) row += C.vSWW * v[C.iVzSWW];
    if (C.SetVzSW  && !gmg_skip_neighbour(BCv_t[C.iVzSW]))  row += C.vSW  * v[C.iVzSW];
    if (C.SetVzSE  && !gmg_skip_neighbour(BCv_t[C.iVzSE]))  row += C.vSE  * v[C.iVzSE];
    if (C.SetVzSEE && !gmg_skip_neighbour(BCv_t[C.iVzSEE])) row += C.vSEE * v[C.iVzSEE];
    if (C.SetVzNWW && !gmg_skip_neighbour(BCv_t[C.iVzNWW])) row += C.vNWW * v[C.iVzNWW];
    if (C.SetVzNW  && !gmg_skip_neighbour(BCv_t[C.iVzNW]))  row += C.vNW  * v[C.iVzNW];
    if (C.SetVzNE  && !gmg_skip_neighbour(BCv_t[C.iVzNE]))  row += C.vNE  * v[C.iVzNE];
    if (C.SetVzNEE && !gmg_skip_neighbour(BCv_t[C.iVzNEE])) row += C.vNEE * v[C.iVzNEE];
    if (C.SetVzNNW && !gmg_skip_neighbour(BCv_t[C.iVzNNW])) row += C.vNNW * v[C.iVzNNW];
    if (C.SetVzNNE && !gmg_skip_neighbour(BCv_t[C.iVzNNE])) row += C.vNNE * v[C.iVzNNE];
    // 6 P neighbours.
    if (C.SetPrSW && !gmg_skip_neighbour(BCp_t[C.iPrSW])) row += C.pSW * p[C.iPrSW];
    if (C.SetPrSE && !gmg_skip_neighbour(BCp_t[C.iPrSE])) row += C.pSE * p[C.iPrSE];
    if (C.SetPrW  && !gmg_skip_neighbour(BCp_t[C.iPrW]))  row += C.pW  * p[C.iPrW];
    if (C.SetPrE  && !gmg_skip_neighbour(BCp_t[C.iPrE]))  row += C.pE  * p[C.iPrE];
    if (C.SetPrNW && !gmg_skip_neighbour(BCp_t[C.iPrNW])) row += C.pNW * p[C.iPrNW];
    if (C.SetPrNE && !gmg_skip_neighbour(BCp_t[C.iPrNE])) row += C.pNE * p[C.iPrNE];

    return row * celvol;
}

// ---------------------------------------------------------------------------
// Newton Z-momentum Jacobian (task 15.20 / design D11 Phase 2).
//
// Ports `Zjacobian_InnerNodesDecoupled3` (StokesAssemblyDecoupled.c:732)
// line-for-line, dropping the AddCoeff3 assembly calls and celvol scaling.
// The Z-momentum Newton row is a 9-Vz / 12-Vx / 6-P stencil.
// ---------------------------------------------------------------------------

typedef struct {
    double vC, vS, vN, vW, vE;
    double vSW, vSE, vNW, vNE;
    double uSW, uSE, uNW, uNE;
    double uSSW, uSSE, uNNW, uNNE;
    double uSWW, uSEE, uNWW, uNEE;
    double pS, pN;
    double pSW, pSE, pNW, pNE;
    int iVzC, iVzS, iVzN, iVzW, iVzE;
    int iVzSW, iVzSE, iVzNW, iVzNE;
    int iVxSW, iVxSE, iVxNW, iVxNE;
    int iVxSSW, iVxSSE, iVxNNW, iVxNNE;
    int iVxSWW, iVxSEE, iVxNWW, iVxNEE;
    int iPrS, iPrN;
    int iPrSW, iPrSE, iPrNW, iPrNE;
    bool SetVzS, SetVzW, SetVzC, SetVzE, SetVzN;
    bool SetVzSW, SetVzSE, SetVzNW, SetVzNE;
    bool SetVxSW, SetVxSE, SetVxNW, SetVxNE;
    bool SetVxSSW, SetVxSSE, SetVxNNW, SetVxNNE;
    bool SetVxSWW, SetVxSEE, SetVxNWW, SetVxNEE;
    bool SetPrS, SetPrN;
    bool SetPrSW, SetPrSE, SetPrNW, SetPrNE;
} ZmomNewtonCoeffs;

static void zmom_fill_coeffs_newton(grid *mesh, params model, int ix, int iz,
                                    int stab, int comp, double om,
                                    double one_dx, double one_dz,
                                    int c1, int c2, int c3, int nx, int nz,
                                    ZmomNewtonCoeffs *C) {
    // `Zjacobian_InnerNodesDecoupled3` uses `ix == k`, `iz == l`.
    const int k = ix, l = iz;
    const int nxvz = nx + 1, ncx = nx - 1;
    const double dx = mesh->dx, dz = mesh->dz;
    double out_of_plane = 1.0, vC_corr = 0.0;
    if (model.out_of_plane == 1) out_of_plane = 3.0 / 2.0;
    const int periodix = model.periodic_x;
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    double SzzBCS = 0., SzzBCN = 0.;
    if (iz == 0    && BCv_t[c3] == 2) SzzBCS = 1.0;
    if (iz == nz-1 && BCv_t[c3] == 2) SzzBCN = 1.0;

    int iVzC   = c3;
    int iVzW   = iVzC - 1,   iVzE   = iVzC + 1;
    int iVzS   = c3 - nxvz,  iVzN   = c3 + nxvz;
    int iVzSW  = iVzS - 1,   iVzSE  = iVzS + 1;
    int iVzNW  = iVzN - 1,   iVzNE  = iVzN + 1;

    int iVxNW  = c1 + nx - 1;
    int iVxNE  = c1 + nx;
    int iVxSW  = c1 - 1;
    int iVxSE  = c1;
    int iVxNWW = iVxNW - 1;
    int iVxNEE = iVxNE + 1;
    int iVxSWW = iVxSW - 1;
    int iVxSEE = iVxSE + 1;
    int iVxNNW = iVxNW + nx;
    int iVxNNE = iVxNE + nx;
    int iVxSSW = iVxSW - nx;
    int iVxSSE = iVxSE - nx;

    int iPrS  = c2;
    int iPrSW = iPrS - 1;
    int iPrSE = iPrS + 1;
    int iPrN  = c2 + ncx;
    int iPrNW = iPrN - 1;
    int iPrNE = iPrN + 1;

    int ixyW  = c1 - 1;
    int ixySW = ixyW - nx;
    int ixyNW = ixyW + nx;
    int ixyE  = c1;
    int ixySE = ixyE - nx;
    int ixyNE = ixyE + nx;

    if (k == 1 && BCv_t[iVzW] == -12) {
        iVxNWW = iVxNW + (nx - 2);
        iVxSWW = iVxSW + (nx - 2);
        iPrNW  = iPrN  + (ncx - 1);
        iPrSW  = iPrS  + (ncx - 1);
        iVzW   = iVzC  + nxvz - 3;
        iVzSW  = iVzW  - nxvz;
        iVzNW  = iVzW  + nxvz;
    }
    if (k == nx - 1 && BCv_t[iVzE] == -12) {
        iVxNEE = iVxNE - (nx - 2);
        iVxSEE = iVxSE - (nx - 2);
        iPrNE  = iPrN  - (ncx - 1);
        iPrSE  = iPrS  - (ncx - 1);
        iVzE   = iVzC  - nxvz + 3;
        iVzNE  = iVzE  + nxvz;
        iVzSE  = iVzE  - nxvz;
    }

    // D-tensor entries.
    double D31W = mesh->D31_s[ixyW];
    double D32W = mesh->D32_s[ixyW];
    double D33W = mesh->D33_s[ixyW];
    double D34W = mesh->D34_s[ixyW];
    double D31E = mesh->D31_s[ixyE];
    double D32E = mesh->D32_s[ixyE];
    double D33E = mesh->D33_s[ixyE];
    double D34E = mesh->D34_s[ixyE];
    double D21S = mesh->D21_n[iPrS];
    double D22S = mesh->D22_n[iPrS];
    double D23S = mesh->D23_n[iPrS];
    double D24S = mesh->D24_n[iPrS];
    double D21N = mesh->D21_n[iPrN];
    double D22N = mesh->D22_n[iPrN];
    double D23N = mesh->D23_n[iPrN];
    double D24N = mesh->D24_n[iPrN];

    double inS = 0.0, inN = 0.0, inW = 0.0, inE = 0.0, inWv = 0.0, inEv = 0.0;
    if (BCp_t[iPrS] == -1) inS = 1.0;
    if (BCp_t[iPrN] == -1) inN = 1.0;
    if (BCv_t[iVzW] != 30 && BCv_t[iVzW] != 13) inW = 1.0;
    if (BCv_t[iVzE] != 30 && BCv_t[iVzE] != 13) inE = 1.0;
    if (BCg_t[ixyW] != 30) inWv = 1.0;
    if (BCg_t[ixyE] != 30) inEv = 1.0;

    double wS = 0.0, wN = 0.0, wE = 0.0, wW = 0.0;
    double inSWc = 0.0, inSEc = 0.0, inNWc = 0.0, inNEc = 0.0;
    double inSWv = 0.0, inSEv = 0.0, inNWv = 0.0, inNEv = 0.0;

    if ((k > 1) || (k == 1 && BCv_t[iVzW] == -1)) {
        if (BCp_t[iPrSW] == -1) inSWc = 1.0;
        if (BCp_t[iPrNW] == -1) inNWc = 1.0;
    }
    if ((k < nxvz - 2) || (k == nxvz - 2 && BCv_t[iVzE] == -1)) {
        if (BCp_t[iPrSE] == -1) inSEc = 1.0;
        if (BCp_t[iPrNE] == -1) inNEc = 1.0;
    }
    if (l > 0) {
        if (BCg_t[ixySW] != 30) inSWv = 1.0;
        if (BCg_t[ixySE] != 30) inSEv = 1.0;
    }
    if (l < nz - 1) {
        if (BCg_t[ixyNW] != 30) inNWv = 1.0;
        if (BCg_t[ixyNE] != 30) inNEv = 1.0;
    }

    double inSSW = 0.0, inSSE = 0.0, inNNW = 0.0, inNNE = 0.0;
    if (BCu_t[iVxSSW] == -1 || periodix == 1) inSSW = 1.0;
    if (BCu_t[iVxSSE] == -1 || periodix == 1) inSSE = 1.0;
    if (BCu_t[iVxNNW] == -1 || periodix == 1) inNNW = 1.0;
    if (BCu_t[iVxNNE] == -1 || periodix == 1) inNNE = 1.0;

    wE = inN + inS + inNEc + inSEc;
    wW = inN + inS + inNWc + inSWc;
    wS = inW + inE + inSWv + inSEv;
    wN = inW + inE + inNWv + inNEv;
    if (wW > 1.0) wW = 1.0 / wW;
    if (wE > 1.0) wE = 1.0 / wE;
    if (wS > 1.0) wS = 1.0 / wS;
    if (wN > 1.0) wN = 1.0 / wN;

    // FD coefficients — verbatim from Zjacobian_InnerNodesDecoupled3.
    double vW  = (-D33W*dz*pow(inW, 2)*inWv + 0.25*dx*pow(inW, 2)*inWv*(D23N*inN - D23S*inS) + (1.0/3.0)*dx*wW*(inNWc - inSWc)*(D31W*comp*out_of_plane + D32W*(comp*out_of_plane - 3)))/(pow(dx, 2)*dz);
    double vC  = (1.0/3.0)*(-dx*(dx*(D21N*comp*inN*out_of_plane + D21S*comp*inS*out_of_plane + D22N*inN*(comp*out_of_plane - 3) + D22S*inS*(comp*out_of_plane - 3)) + 0.75*dz*(-D23N*inN + D23S*inS)*(pow(inE, 2)*inEv - pow(inW, 2)*inWv)) + dz*(dx*(inN - inS)*(-D31E*comp*out_of_plane*wE + D31W*comp*out_of_plane*wW - D32E*wE*(comp*out_of_plane - 3) + D32W*wW*(comp*out_of_plane - 3)) + 3*dz*(D33E*pow(inE, 2)*inEv + D33W*pow(inW, 2)*inWv)))/(pow(dx, 2)*pow(dz, 2));
    double vE  = (-D33E*dz*pow(inE, 2)*inEv + 0.25*dx*pow(inE, 2)*inEv*(-D23N*inN + D23S*inS) - 1.0/3.0*dx*wE*(inNEc - inSEc)*(D31E*comp*out_of_plane + D32E*(comp*out_of_plane - 3)))/(pow(dx, 2)*dz);
    double vS  = (1.0/3.0)*inS*(-0.75*D23S*dz*(inE*inSEv - inSWv*inW) + dx*(D21S*comp*out_of_plane + D22S*(comp*out_of_plane - 3)) + dz*(-D31E*comp*out_of_plane*wE + D31W*comp*out_of_plane*wW - D32E*wE*(comp*out_of_plane - 3) + D32W*wW*(comp*out_of_plane - 3)))/(dx*pow(dz, 2));
    double vN  = (1.0/3.0)*inN*(0.75*D23N*dz*(inE*inNEv - inNWv*inW) + dx*(D21N*comp*out_of_plane + D22N*(comp*out_of_plane - 3)) + dz*(D31E*comp*out_of_plane*wE - D31W*comp*out_of_plane*wW + D32E*wE*(comp*out_of_plane - 3) - D32W*wW*(comp*out_of_plane - 3)))/(dx*pow(dz, 2));
    double uSW = (1.0/3.0)*(dx*(0.75*dx*(D23N*inN*inW*inWv + D23S*inS*(inSSW*inSWv - inW*inWv)) + dz*inS*(D21S*(comp*out_of_plane - 3) + D22S*comp*out_of_plane)) - dz*(3*D33W*dx*inW*inWv + dz*(D31E*inS*wE*(comp*out_of_plane - 3) - D31W*wW*(inS - inSWc)*(comp*out_of_plane - 3) + D32E*comp*inS*out_of_plane*wE - D32W*comp*out_of_plane*wW*(inS - inSWc))))/(pow(dx, 2)*pow(dz, 2));
    double uSE = (1.0/3.0)*(-dx*(0.75*dx*(-D23N*inE*inEv*inN + D23S*inS*(inE*inEv - inSEv*inSSE)) + dz*inS*(D21S*(comp*out_of_plane - 3) + D22S*comp*out_of_plane)) + dz*(3*D33E*dx*inE*inEv + dz*(D31E*wE*(inS - inSEc)*(comp*out_of_plane - 3) - D31W*inS*wW*(comp*out_of_plane - 3) + D32E*comp*out_of_plane*wE*(inS - inSEc) - D32W*comp*inS*out_of_plane*wW)))/(pow(dx, 2)*pow(dz, 2));
    double uNW = (1.0/3.0)*(dx*(0.75*dx*(D23N*inN*(inNNW*inNWv - inW*inWv) + D23S*inS*inW*inWv) - dz*inN*(D21N*(comp*out_of_plane - 3) + D22N*comp*out_of_plane)) + dz*(3*D33W*dx*inW*inWv + dz*(-D31E*inN*wE*(comp*out_of_plane - 3) + D31W*wW*(inN - inNWc)*(comp*out_of_plane - 3) - D32E*comp*inN*out_of_plane*wE + D32W*comp*out_of_plane*wW*(inN - inNWc))))/(pow(dx, 2)*pow(dz, 2));
    double uNE = (1.0/3.0)*(dx*(0.75*dx*(-D23N*inN*(inE*inEv - inNEv*inNNE) + D23S*inE*inEv*inS) + dz*inN*(D21N*(comp*out_of_plane - 3) + D22N*comp*out_of_plane)) - dz*(3*D33E*dx*inE*inEv + dz*(-D31E*wE*(inN - inNEc)*(comp*out_of_plane - 3) + D31W*inN*wW*(comp*out_of_plane - 3) - D32E*comp*out_of_plane*wE*(inN - inNEc) + D32W*comp*inN*out_of_plane*wW)))/(pow(dx, 2)*pow(dz, 2));
    double vSW = (1.0/3.0)*(-0.75*D23S*inS*inSWv*inW + D31W*comp*inSWc*out_of_plane*wW + D32W*inSWc*wW*(comp*out_of_plane - 3))/(dx*dz);
    double vSE = (1.0/3.0)*(0.75*D23S*inE*inS*inSEv - D31E*comp*inSEc*out_of_plane*wE - D32E*inSEc*wE*(comp*out_of_plane - 3))/(dx*dz);
    double vNW = (1.0/3.0)*(0.75*D23N*inN*inNWv*inW - D31W*comp*inNWc*out_of_plane*wW - D32W*inNWc*wW*(comp*out_of_plane - 3))/(dx*dz);
    double vNE = (1.0/3.0)*(-0.75*D23N*inE*inN*inNEv + D31E*comp*inNEc*out_of_plane*wE + D32E*inNEc*wE*(comp*out_of_plane - 3))/(dx*dz);
    double uSWW = (1.0/3.0)*inSWc*wW*(D31W*(comp*out_of_plane - 3) + D32W*comp*out_of_plane)/pow(dx, 2);
    double uSEE = (1.0/3.0)*inSEc*wE*(D31E*(comp*out_of_plane - 3) + D32E*comp*out_of_plane)/pow(dx, 2);
    double uNWW = (1.0/3.0)*inNWc*wW*(D31W*(comp*out_of_plane - 3) + D32W*comp*out_of_plane)/pow(dx, 2);
    double uNEE = (1.0/3.0)*inNEc*wE*(D31E*(comp*out_of_plane - 3) + D32E*comp*out_of_plane)/pow(dx, 2);
    double uSSW = -0.25*D23S*inS*inSSW*inSWv/pow(dz, 2);
    double uSSE = -0.25*D23S*inS*inSEv*inSSE/pow(dz, 2);
    double uNNW = -0.25*D23N*inN*inNNW*inNWv/pow(dz, 2);
    double uNNE = -0.25*D23N*inN*inNEv*inNNE/pow(dz, 2);
    double pS  = -inS*one_dz + (inS*(D24S*dx + dz*(-D34E*inE*wE + D34W*inW*wW))/(dx*dz));
    double pN  =  inN*one_dz + (inN*(-D24N*dx + dz*(-D34E*inE*wE + D34W*inW*wW))/(dx*dz));
    double pSW = (D34W*inSWc*inW*wW/dx);
    double pSE = (-D34E*inE*inSEc*wE/dx);
    double pNW = (D34W*inNWc*inW*wW/dx);
    double pNE = (-D34E*inE*inNEc*wE/dx);

    // Stabilisation with density gradients — mirrors Zjacobian branch
    // exactly (which does NOT reflect Picard's min(vC+corr,vW,vE) guard).
    if (stab == 1) {
        double drhodz = (mesh->rho_n[c2 + ncx] - mesh->rho_n[c2]) * one_dz;
        vC_corr = 1.00 * om * model.dt * mesh->gz[c3] * drhodz;
        if (vC + vC_corr < 0.0) vC_corr = 0.0;
        vC += vC_corr;
    }

    // Non-conforming Dirichlet corrections — mirror Zjacobian
    // StokesAssemblyDecoupled.c:934-943 exactly.
    if (BCv_t[iVzW]   == 11) vC  -= vW;
    if (BCv_t[iVzE]   == 11) vC  -= vE;
    if (BCv_t[iVzSW]  == 11) vW  -= vSW;
    if (BCv_t[iVzSE]  == 11) vE  -= vSE;
    if (BCv_t[iVzNW]  == 11) vW  -= vNW;
    if (BCv_t[iVzNE]  == 11) vE  -= vNE;
    if (BCu_t[iVxNNW] == 11) uNW -= uNNW;
    if (BCu_t[iVxNNE] == 11) uNE -= uNNE;
    if (BCu_t[iVxSSW] == 11) uSW -= uSSW;
    if (BCu_t[iVxSSE] == 11) uSE -= uSSE;

    (void)one_dx;

    // Set flags — mirror Zjacobian_InnerNodesDecoupled3 assembly branch
    // (StokesAssemblyDecoupled.c:948-1053) for SetVx*/SetVz*/SetPr*,
    // SetVxSSW/E, SetVxNNW/E, SetVxSWW/SEE/NWW/NEE (guarded by the
    // `(k>1)|| ...` k-edge patterns and `Newton==1 && ...` pressure
    // corners).
    const bool SetVzSW = BCv_t[iVzSW] != 30 && iz > 0;
    const bool SetVzS  = BCv_t[iVzS]  != 30 && iz > 0;
    const bool SetVzSE = BCv_t[iVzSE] != 30 && iz > 0;
    const bool SetVzW  = BCv_t[iVzW]  != 30;
    const bool SetVzC  = BCv_t[iVzC]  != 30;
    const bool SetVzE  = BCv_t[iVzE]  != 30;
    const bool SetVzNW = BCv_t[iVzNW] != 30 && iz < nz - 1;
    const bool SetVzN  = BCv_t[iVzN]  != 30 && iz < nz - 1;
    const bool SetVzNE = BCv_t[iVzNE] != 30 && iz < nz - 1;

    const bool SetVxSW = BCu_t[iVxSW] != 30 && iz > 0;
    const bool SetVxSE = BCu_t[iVxSE] != 30 && iz > 0;
    const bool SetVxNW = BCu_t[iVxNW] != 30 && iz < nz - 1;
    const bool SetVxNE = BCu_t[iVxNE] != 30 && iz < nz - 1;

    // SSW/SSE/NNW/NNE — unconditional BC != 30 in the assembler.
    const bool SetVxSSW = BCu_t[iVxSSW] != 30;
    const bool SetVxSSE = BCu_t[iVxSSE] != 30;
    const bool SetVxNNW = BCu_t[iVxNNW] != 30;
    const bool SetVxNNE = BCu_t[iVxNNE] != 30;

    // SWW/SEE/NWW/NEE — guarded by k-edge patterns.
    const bool SetVxSWW = ((k > 1) || (k == 1 && BCu_t[iVxSWW] == -1)) && BCu_t[iVxSWW] != 30;
    const bool SetVxSEE = ((k < nx - 1) || (k == nx - 1 && BCu_t[iVxSEE] == -1)) && BCu_t[iVxSEE] != 30;
    const bool SetVxNWW = ((k > 1) || (k == 1 && BCu_t[iVxNWW] == -1)) && BCu_t[iVxNWW] != 30;
    const bool SetVxNEE = ((k < nx - 1) || (k == nx - 1 && BCu_t[iVxNEE] == -1)) && BCu_t[iVxNEE] != 30;

    // Pressure gates — pS & pN require inS & inN (both be-interior); the
    // assembler writes them only when `BCp_t[iPrS]!=30 && BCp_t[iPrN]!=30`.
    const bool SetPrS = (BCp_t[iPrS] != 30) && (BCp_t[iPrN] != 30);
    const bool SetPrN = (BCp_t[iPrS] != 30) && (BCp_t[iPrN] != 30);

    // Pressure corners — `Newton==1 && (k>1 || (k==1 && periodix))` or
    // symmetric east pattern.
    const bool SetPrSW = ((k > 1) || (k == 1 && periodix == 1)) && (BCp_t[iPrSW] != 30);
    const bool SetPrSE = ((k < nx - 1) || (k == nx - 1 && periodix == 1)) && (BCp_t[iPrSE] != 30);
    const bool SetPrNW = ((k > 1) || (k == 1 && periodix == 1)) && (BCp_t[iPrNW] != 30);
    const bool SetPrNE = ((k < nx - 1) || (k == nx - 1 && periodix == 1)) && (BCp_t[iPrNE] != 30);

    // Pack.
    C->vC = vC; C->vS = vS; C->vN = vN; C->vW = vW; C->vE = vE;
    C->vSW = vSW; C->vSE = vSE; C->vNW = vNW; C->vNE = vNE;
    C->uSW = uSW; C->uSE = uSE; C->uNW = uNW; C->uNE = uNE;
    C->uSSW = uSSW; C->uSSE = uSSE; C->uNNW = uNNW; C->uNNE = uNNE;
    C->uSWW = uSWW; C->uSEE = uSEE; C->uNWW = uNWW; C->uNEE = uNEE;
    C->pS = pS; C->pN = pN;
    C->pSW = pSW; C->pSE = pSE; C->pNW = pNW; C->pNE = pNE;

    C->iVzC = iVzC; C->iVzS = iVzS; C->iVzN = iVzN; C->iVzW = iVzW; C->iVzE = iVzE;
    C->iVzSW = iVzSW; C->iVzSE = iVzSE; C->iVzNW = iVzNW; C->iVzNE = iVzNE;
    C->iVxSW = iVxSW; C->iVxSE = iVxSE; C->iVxNW = iVxNW; C->iVxNE = iVxNE;
    C->iVxSSW = iVxSSW; C->iVxSSE = iVxSSE; C->iVxNNW = iVxNNW; C->iVxNNE = iVxNNE;
    C->iVxSWW = iVxSWW; C->iVxSEE = iVxSEE; C->iVxNWW = iVxNWW; C->iVxNEE = iVxNEE;
    C->iPrS = iPrS; C->iPrN = iPrN;
    C->iPrSW = iPrSW; C->iPrSE = iPrSE; C->iPrNW = iPrNW; C->iPrNE = iPrNE;

    C->SetVzS = SetVzS; C->SetVzW = SetVzW; C->SetVzC = SetVzC;
    C->SetVzE = SetVzE; C->SetVzN = SetVzN;
    C->SetVzSW = SetVzSW; C->SetVzSE = SetVzSE;
    C->SetVzNW = SetVzNW; C->SetVzNE = SetVzNE;
    C->SetVxSW = SetVxSW; C->SetVxSE = SetVxSE;
    C->SetVxNW = SetVxNW; C->SetVxNE = SetVxNE;
    C->SetVxSSW = SetVxSSW; C->SetVxSSE = SetVxSSE;
    C->SetVxNNW = SetVxNNW; C->SetVxNNE = SetVxNNE;
    C->SetVxSWW = SetVxSWW; C->SetVxSEE = SetVxSEE;
    C->SetVxNWW = SetVxNWW; C->SetVxNEE = SetVxNEE;
    C->SetPrS = SetPrS; C->SetPrN = SetPrN;
    C->SetPrSW = SetPrSW; C->SetPrSE = SetPrSE;
    C->SetPrNW = SetPrNW; C->SetPrNE = SetPrNE;
}

static double zmom_row_matvec_newton(grid *mesh, params model, int ix, int iz,
                                     int stab, int comp, double om,
                                     double one_dx, double one_dz, double celvol,
                                     int c1, int c2, int c3, int nx, int nz,
                                     const double *u, const double *v, const double *p) {
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCp_t = mesh->BCp.type;
    ZmomNewtonCoeffs C;
    zmom_fill_coeffs_newton(mesh, model, ix, iz, stab, comp, om, one_dx, one_dz,
                            c1, c2, c3, nx, nz, &C);

    double row = 0.0;
    // 12 Vx neighbours (SSW/SSE, SWW/SW/SE/SEE, NWW/NW/NE/NEE, NNW/NNE).
    if (C.SetVxSSW && !gmg_skip_neighbour(BCu_t[C.iVxSSW])) row += C.uSSW * u[C.iVxSSW];
    if (C.SetVxSSE && !gmg_skip_neighbour(BCu_t[C.iVxSSE])) row += C.uSSE * u[C.iVxSSE];
    if (C.SetVxSWW && !gmg_skip_neighbour(BCu_t[C.iVxSWW])) row += C.uSWW * u[C.iVxSWW];
    if (C.SetVxSW  && !gmg_skip_neighbour(BCu_t[C.iVxSW]))  row += C.uSW  * u[C.iVxSW];
    if (C.SetVxSE  && !gmg_skip_neighbour(BCu_t[C.iVxSE]))  row += C.uSE  * u[C.iVxSE];
    if (C.SetVxSEE && !gmg_skip_neighbour(BCu_t[C.iVxSEE])) row += C.uSEE * u[C.iVxSEE];
    if (C.SetVxNWW && !gmg_skip_neighbour(BCu_t[C.iVxNWW])) row += C.uNWW * u[C.iVxNWW];
    if (C.SetVxNW  && !gmg_skip_neighbour(BCu_t[C.iVxNW]))  row += C.uNW  * u[C.iVxNW];
    if (C.SetVxNE  && !gmg_skip_neighbour(BCu_t[C.iVxNE]))  row += C.uNE  * u[C.iVxNE];
    if (C.SetVxNEE && !gmg_skip_neighbour(BCu_t[C.iVxNEE])) row += C.uNEE * u[C.iVxNEE];
    if (C.SetVxNNW && !gmg_skip_neighbour(BCu_t[C.iVxNNW])) row += C.uNNW * u[C.iVxNNW];
    if (C.SetVxNNE && !gmg_skip_neighbour(BCu_t[C.iVxNNE])) row += C.uNNE * u[C.iVxNNE];
    // 9 Vz neighbours.
    if (C.SetVzSW && !gmg_skip_neighbour(BCv_t[C.iVzSW])) row += C.vSW * v[C.iVzSW];
    if (C.SetVzS  && !gmg_skip_neighbour(BCv_t[C.iVzS]))  row += C.vS  * v[C.iVzS];
    if (C.SetVzSE && !gmg_skip_neighbour(BCv_t[C.iVzSE])) row += C.vSE * v[C.iVzSE];
    if (C.SetVzW  && !gmg_skip_neighbour(BCv_t[C.iVzW]))  row += C.vW  * v[C.iVzW];
    if (C.SetVzC  && !gmg_skip_neighbour(BCv_t[C.iVzC]))  row += C.vC  * v[C.iVzC];
    if (C.SetVzE  && !gmg_skip_neighbour(BCv_t[C.iVzE]))  row += C.vE  * v[C.iVzE];
    if (C.SetVzNW && !gmg_skip_neighbour(BCv_t[C.iVzNW])) row += C.vNW * v[C.iVzNW];
    if (C.SetVzN  && !gmg_skip_neighbour(BCv_t[C.iVzN]))  row += C.vN  * v[C.iVzN];
    if (C.SetVzNE && !gmg_skip_neighbour(BCv_t[C.iVzNE])) row += C.vNE * v[C.iVzNE];
    // 6 P neighbours.
    if (C.SetPrSW && !gmg_skip_neighbour(BCp_t[C.iPrSW])) row += C.pSW * p[C.iPrSW];
    if (C.SetPrS  && !gmg_skip_neighbour(BCp_t[C.iPrS]))  row += C.pS  * p[C.iPrS];
    if (C.SetPrSE && !gmg_skip_neighbour(BCp_t[C.iPrSE])) row += C.pSE * p[C.iPrSE];
    if (C.SetPrNW && !gmg_skip_neighbour(BCp_t[C.iPrNW])) row += C.pNW * p[C.iPrNW];
    if (C.SetPrN  && !gmg_skip_neighbour(BCp_t[C.iPrN]))  row += C.pN  * p[C.iPrN];
    if (C.SetPrNE && !gmg_skip_neighbour(BCp_t[C.iPrNE])) row += C.pNE * p[C.iPrNE];

    return row * celvol;
}

// ---------------------------------------------------------------------------
// Z-momentum coefficient extractor (mirrors Zmomentum_InnerNodesDecoupled).
// ---------------------------------------------------------------------------
static void zmom_fill_coeffs(grid *mesh, params model, int ix, int iz,
                             int stab, int comp, double om,
                             double one_dx, double one_dz,
                             int c1, int c2, int c3, int nx, int nz,
                             ZmomCoeffs *C) {
    const int nxvz = nx + 1, ncx = nx - 1;
    const double dx = mesh->dx, dz = mesh->dz;
    double OOP = 1.0, vC_corr = 0.0;
    if (model.out_of_plane == 1) OOP = 3.0 / 2.0;
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCg_t = mesh->BCg.type, *BCp_t = mesh->BCp.type;

    double SzzBCS = 0., SzzBCN = 0.;
    if (iz == 0    && BCv_t[c3] == 2) SzzBCS = 1.0;
    if (iz == nz-1 && BCv_t[c3] == 2) SzzBCN = 1.0;

    int iVzC  = c3;
    int iVzW  = iVzC - 1;
    int iVzE  = iVzC + 1;
    int iVzS  = iVzC - nxvz;
    int iVzN  = iVzC + nxvz;
    int iVxNW = c1 + nx - 1;
    int iVxNE = c1 + nx;
    int iVxSW = c1 - 1;
    int iVxSE = c1;
    int iPrS  = c2;
    int iPrN  = c2 + ncx;
    int ixyE  = c1;
    int ixyW  = c1 - 1;

    if (BCv_t[iVzW] == -12 && ix == 1)        iVzW = c3 + nxvz - 3;
    if (BCv_t[iVzE] == -12 && ix == nxvz - 2) iVzE = c3 - nxvz + 3;

    const double D22S = mesh->D22_n[iPrS];
    const double D22N = mesh->D22_n[iPrN];
    const double D33E = mesh->D33_s[ixyE];
    const double D33W = mesh->D33_s[ixyW];

    double inN = 0.0, inS = 0.0;
    if (BCp_t[iPrS] == -1) inS = 1.0;
    if (BCp_t[iPrN] == -1) inN = 1.0;
    (void)inS; (void)inN;

    double inW = 1.0, inE = 1.0, RegW = 1., RegE = 1.;
    double NeuW = 0., NeuE = 0., DirW = 0., DirE = 0.;
    if (BCg_t[ixyW] == 30) inW = 0.0;
    if (BCg_t[ixyE] == 30) inE = 0.0;
    if (BCv_t[iVzW] == 13) { NeuW = 1.0; RegW = 0.0; }
    if (BCv_t[iVzE] == 13) { NeuE = 1.0; RegE = 0.0; }
    if (BCv_t[iVzW] == 11) { DirW = 1.0; RegW = 0.0; }
    if (BCv_t[iVzE] == 11) { DirE = 1.0; RegE = 0.0; }
    (void)NeuW; (void)NeuE;

    double DirSW = 0., DirNW = 0.;
    if (mesh->BCu.type[c1 - 1]    == 11) DirSW = 1.0;
    if (mesh->BCu.type[c1 + nx-1] == 11) DirNW = 1.0;
    double DirSE = 0., DirNE = 0.;
    if (mesh->BCu.type[c1]        == 11) DirSE = 1.0;
    if (mesh->BCu.type[c1 + nx]   == 11) DirNE = 1.0;

    double vW = (1.0/2.0)*D33W*RegW*inW*(SzzBCN + SzzBCS - 2)/pow(dx, 2);
    double vC = (1.0/6.0)*(-SzzBCN*(2*D22S*pow(dx, 2)*(OOP*comp - 3) - 3*pow(dz, 2)*(D33E*inE*(2*DirE + RegE) + D33W*inW*(2*DirW + RegW))) - SzzBCS*(2*D22N*pow(dx, 2)*(OOP*comp - 3) - 3*pow(dz, 2)*(D33E*inE*(2*DirE + RegE) + D33W*inW*(2*DirW + RegW))) + 2*(pow(dx, 2)*(D22N + D22S)*(OOP*comp - 3) - 3*pow(dz, 2)*(D33E*inE*(2*DirE + RegE) + D33W*inW*(2*DirW + RegW)))*(SzzBCN + SzzBCS - 1))/(pow(dx, 2)*pow(dz, 2));
    double vE = (1.0/2.0)*D33E*RegE*inE*(SzzBCN + SzzBCS - 2)/pow(dx, 2);
    double vS = (1.0/3.0)*D22S*(1 - SzzBCS)*(OOP*comp - 3)/pow(dz, 2);
    double vN = (1.0/3.0)*D22N*(1 - SzzBCN)*(OOP*comp - 3)/pow(dz, 2);
    double uSW = (1.0/6.0)*(-3*D33W*SzzBCS*inW*(2*DirNW*SzzBCN - SzzBCN - SzzBCS + 1) + SzzBCN*(2*D22S*OOP*comp - 3*D33W*inW*(2*DirNW*SzzBCN - SzzBCN - SzzBCS + 1)) - 2*(D22S*OOP*comp - 3*D33W*inW*(2*DirNW*SzzBCN - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    double uSE = (1.0/6.0)*(3*D33E*SzzBCS*inE*(2*DirNE*SzzBCN - SzzBCN - SzzBCS + 1) - SzzBCN*(2*D22S*OOP*comp - 3*D33E*inE*(2*DirNE*SzzBCN - SzzBCN - SzzBCS + 1)) + 2*(D22S*OOP*comp - 3*D33E*inE*(2*DirNE*SzzBCN - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    double uNW = (1.0/6.0)*(3*D33W*SzzBCN*inW*(2*DirSW*SzzBCS - SzzBCN - SzzBCS + 1) - SzzBCS*(2*D22N*OOP*comp - 3*D33W*inW*(2*DirSW*SzzBCS - SzzBCN - SzzBCS + 1)) + 2*(D22N*OOP*comp - 3*D33W*inW*(2*DirSW*SzzBCS - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    double uNE = (1.0/6.0)*(-3*D33E*SzzBCN*inE*(2*DirSE*SzzBCS - SzzBCN - SzzBCS + 1) + SzzBCS*(2*D22N*OOP*comp - 3*D33E*inE*(2*DirSE*SzzBCS - SzzBCN - SzzBCS + 1)) - 2*(D22N*OOP*comp - 3*D33E*inE*(2*DirSE*SzzBCS - SzzBCN - SzzBCS + 1))*(SzzBCN + SzzBCS - 1))/(dx*dz);
    double pS = (SzzBCS - 1)/dz;
    double pN = (1 - SzzBCN)/dz;

    if (stab == 1) {
        double drhodz = (mesh->rho_n[c2 + ncx] - mesh->rho_n[c2]) * one_dz;
        vC_corr = 1.00 * om * model.dt * mesh->gz[c3] * drhodz;
        if ((vC + vC_corr) < 0.0 || (vC + vC_corr) < vW || (vC + vC_corr) < vE) vC_corr = 0.0;
        vC += vC_corr;
    }
    (void)one_dx;
    (void)ncx;

    const bool SetVzS  = BCv_t[iVzS] != 30 && iz > 0;
    const bool SetVzW  = BCv_t[iVzW] != 30;
    const bool SetVzC  = BCv_t[iVzC] != 30;
    const bool SetVzE  = BCv_t[iVzE] != 30;
    const bool SetVzN  = BCv_t[iVzN] != 30 && iz < nz - 1;
    const bool SetVxSW = BCu_t[iVxSW] != 30 && iz > 0;
    const bool SetVxSE = BCu_t[iVxSE] != 30 && iz > 0;
    const bool SetVxNW = BCu_t[iVxNW] != 30 && iz < nz - 1;
    const bool SetVxNE = BCu_t[iVxNE] != 30 && iz < nz - 1;
    const bool SetPrS  = BCp_t[iPrS] != 30 && iz > 0;
    const bool SetPrN  = BCp_t[iPrN] != 30 && iz < nz - 1;

    C->vS = vS; C->vN = vN; C->vW = vW; C->vE = vE; C->vC = vC;
    C->uSW = uSW; C->uSE = uSE; C->uNW = uNW; C->uNE = uNE;
    C->pS = pS; C->pN = pN;
    C->iVzS = iVzS; C->iVzW = iVzW; C->iVzC = iVzC; C->iVzE = iVzE; C->iVzN = iVzN;
    C->iVxSW = iVxSW; C->iVxSE = iVxSE; C->iVxNW = iVxNW; C->iVxNE = iVxNE;
    C->iPrS = iPrS; C->iPrN = iPrN;
    C->SetVzS = SetVzS; C->SetVzW = SetVzW; C->SetVzC = SetVzC;
    C->SetVzE = SetVzE; C->SetVzN = SetVzN;
    C->SetVxSW = SetVxSW; C->SetVxSE = SetVxSE;
    C->SetVxNW = SetVxNW; C->SetVxNE = SetVxNE;
    C->SetPrS = SetPrS; C->SetPrN = SetPrN;
}

static double zmom_row_matvec(grid *mesh, params model, int ix, int iz,
                              int stab, int comp, double om,
                              double one_dx, double one_dz, double celvol,
                              int c1, int c2, int c3, int nx, int nz,
                              const double *u, const double *v, const double *p) {
    signed char *BCv_t = mesh->BCv.type, *BCu_t = mesh->BCu.type,
                *BCp_t = mesh->BCp.type;
    ZmomCoeffs C;
    zmom_fill_coeffs(mesh, model, ix, iz, stab, comp, om, one_dx, one_dz,
                     c1, c2, c3, nx, nz, &C);

    double row = 0.0;
    if (C.SetVxSW && !gmg_skip_neighbour(BCu_t[C.iVxSW])) row += C.uSW * u[C.iVxSW];
    if (C.SetVxSE && !gmg_skip_neighbour(BCu_t[C.iVxSE])) row += C.uSE * u[C.iVxSE];
    if (C.SetVxNW && !gmg_skip_neighbour(BCu_t[C.iVxNW])) row += C.uNW * u[C.iVxNW];
    if (C.SetVxNE && !gmg_skip_neighbour(BCu_t[C.iVxNE])) row += C.uNE * u[C.iVxNE];
    if (C.SetVzS  && !gmg_skip_neighbour(BCv_t[C.iVzS]))  row += C.vS  * v[C.iVzS];
    if (C.SetVzW  && !gmg_skip_neighbour(BCv_t[C.iVzW]))  row += C.vW  * v[C.iVzW];
    if (C.SetVzC  && !gmg_skip_neighbour(BCv_t[C.iVzC]))  row += C.vC  * v[C.iVzC];
    if (C.SetVzE  && !gmg_skip_neighbour(BCv_t[C.iVzE]))  row += C.vE  * v[C.iVzE];
    if (C.SetVzN  && !gmg_skip_neighbour(BCv_t[C.iVzN]))  row += C.vN  * v[C.iVzN];
    if (C.SetPrS  && !gmg_skip_neighbour(BCp_t[C.iPrS]))  row += C.pS  * p[C.iPrS];
    if (C.SetPrN  && !gmg_skip_neighbour(BCp_t[C.iPrN]))  row += C.pN  * p[C.iPrN];

    return row * celvol;
}

// ---------------------------------------------------------------------------
// Continuity coefficient extractor (mirrors Continuity_InnerNodesDecoupled).
// ---------------------------------------------------------------------------
static void cont_fill_coeffs(grid *mesh, params model, int comp,
                             double one_dx, double one_dz,
                             int c1, int c3, int nxvz,
                             ContCoeffs *C) {
    (void)comp;
    double oop_fact = 1.0;
    const double adv_rho = 1.0;
    if (model.out_of_plane == 1) oop_fact = 3.0 / 2.0;

    // Note: in compressible mode the assembler writes pc = β/dt into
    // StokesC->F (the residual), NOT into StokesD (the p-p block).
    // StokesD is therefore zero-matrix on every row; the β/dt term
    // couples only through the nonlinear residual, not through A·x.
    // See StokesAssemblyDecoupled.c:1120-1127 for the Assemble path.

    const int iVxW = c1;
    const int iVxE = c1 + 1;
    const int iVzS = c3;
    const int iVzN = c3 + nxvz;

    const bool SetVxW = (mesh->BCu.type[iVxW] != 13);
    const bool SetVxE = (mesh->BCu.type[iVxE] != 13);
    const bool SetVzS = (mesh->BCv.type[iVzS] != 13);
    const bool SetVzN = (mesh->BCv.type[iVzN] != 13);

    C->uW = SetVxW ? (-adv_rho * one_dx * oop_fact) : 0.0;
    C->uE = SetVxE ? ( adv_rho * one_dx * oop_fact) : 0.0;
    C->vS = SetVzS ? (-adv_rho * one_dz * oop_fact) : 0.0;
    C->vN = SetVzN ? ( adv_rho * one_dz * oop_fact) : 0.0;

    C->iVxW = iVxW; C->iVxE = iVxE; C->iVzS = iVzS; C->iVzN = iVzN;
    C->SetVxW = SetVxW; C->SetVxE = SetVxE;
    C->SetVzS = SetVzS; C->SetVzN = SetVzN;
}

static double cont_row_matvec(grid *mesh, params model, int comp,
                              double one_dx, double one_dz, double celvol,
                              int c1, int c2, int c3, int nx, int nxvz,
                              const double *u, const double *v, const double *p) {
    (void)p; (void)c2; (void)nx;
    ContCoeffs C;
    cont_fill_coeffs(mesh, model, comp, one_dx, one_dz, c1, c3, nxvz, &C);

    double row = 0.0;
    if (C.SetVxW && !gmg_skip_neighbour(mesh->BCu.type[C.iVxW])) row += C.uW * u[C.iVxW];
    if (C.SetVxE && !gmg_skip_neighbour(mesh->BCu.type[C.iVxE])) row += C.uE * u[C.iVxE];
    if (C.SetVzS && !gmg_skip_neighbour(mesh->BCv.type[C.iVzS])) row += C.vS * v[C.iVzS];
    if (C.SetVzN && !gmg_skip_neighbour(mesh->BCv.type[C.iVzN])) row += C.vN * v[C.iVzN];

    return row * celvol;
}

// ---------------------------------------------------------------------------
// Top-level dispatch (mirrors BuildStokesOperatorDecoupled's iterator).
// ---------------------------------------------------------------------------
void ApplyStokesOperatorMDOODZ(grid *mesh, params model,
                               const double *u, const double *v, const double *p,
                               double *ru, double *rv, double *rp) {
    const int nx = mesh->Nx, nz = mesh->Nz;
    const int nxvz = nx + 1, nzvx = nz + 1, ncx = nx - 1, ncz = nz - 1;
    const double celvol = mesh->dx * mesh->dz;
    const double one_dx = 1.0 / mesh->dx;
    const double one_dz = 1.0 / mesh->dz;

    int comp = 0, stab = 0;
    if (model.free_surface_stab > 0) stab = 1;
    // Note: the Jacobian path in BuildJacobianOperatorDecoupled forces
    // `comp = 1` unconditionally (see StokesAssemblyDecoupled.c:1822).
    // The Picard path (BuildStokesOperatorDecoupled) sets it from
    // `model.compressible`. Mirror the dispatching assembler exactly.
    if (model.Newton == 1) {
        comp = 1;
    } else {
        if (model.compressible == 1) comp = 1;
    }
    const double theta = model.free_surface_stab;

    if (ru) memset(ru, 0, (size_t)nx   * (size_t)nzvx * sizeof(double));
    if (rv) memset(rv, 0, (size_t)nxvz * (size_t)nz   * sizeof(double));
    if (rp) memset(rp, 0, (size_t)ncx  * (size_t)ncz  * sizeof(double));

    // ---------------- U-momentum ----------------
    for (int l = 0; l < nzvx; ++l) {
        for (int k = 0; k < nx; ++k) {
            const int c1 = k + l * nx;
            const int c2 = (k - 1) + (l - 1) * ncx;
            const int c3 = k + l * nxvz;
            if (!(mesh->BCu.type[c1] == -1 || mesh->BCu.type[c1] == 2 || mesh->BCu.type[c1] == -2))
                continue;
            if (!(k >= 0 && k <= nx - 1 && l > 0 && l < nzvx - 1))
                continue;
            if (model.Newton == 1) {
                ru[c1] = xmom_row_matvec_newton(mesh, model, k, l, stab, comp, theta,
                                                one_dx, one_dz, celvol,
                                                c1, c2, c3, nx, nz, u, v, p);
            } else {
                ru[c1] = xmom_row_matvec(mesh, model, k, l, stab, comp, theta,
                                         one_dx, one_dz, celvol,
                                         c1, c2, c3, nx, nz, u, v, p);
            }
        }
    }

    // ---------------- V-momentum ----------------
    for (int l = 0; l < nz; ++l) {
        for (int k = 0; k < nxvz; ++k) {
            const int c3 = k + l * nxvz;
            const int c1 = k + l * nx;
            const int c2 = (k - 1) + (l - 1) * ncx;
            if (!(mesh->BCv.type[c3] == -1 || mesh->BCv.type[c3] == 2))
                continue;
            if (!(k > 0 && k < nxvz - 1 && l >= 0 && l <= nz - 1))
                continue;
            if (model.Newton == 1) {
                rv[c3] = zmom_row_matvec_newton(mesh, model, k, l, stab, comp, theta,
                                                one_dx, one_dz, celvol,
                                                c1, c2, c3, nx, nz, u, v, p);
            } else {
                rv[c3] = zmom_row_matvec(mesh, model, k, l, stab, comp, theta,
                                         one_dx, one_dz, celvol,
                                         c1, c2, c3, nx, nz, u, v, p);
            }
        }
    }

    // ---------------- Continuity ----------------
    for (int l = 0; l < ncz; ++l) {
        for (int k = 0; k < ncx; ++k) {
            const int c2 = k + l * ncx;
            const int c1 = k + (l + 1) * nx;
            const int c3 = k + l * nxvz + 1;
            if (mesh->BCp.type[c2] != -1) continue;
            int cell_comp = (comp == 1) ? mesh->comp_cells[c2] : 0;
            rp[c2] = cont_row_matvec(mesh, model, cell_comp,
                                     one_dx, one_dz, celvol,
                                     c1, c2, c3, nx, nxvz, u, v, p);
        }
    }
}

// ---------------------------------------------------------------------------
// 5×5 Vanka block builder (task 15.6).
//
// Block DOF order: [Vx_L, Vx_R, Vz_B, Vz_T, P_C] for pressure cell (i, j)
// in MDOODZ-padded full-grid terms:
//   Vx_L = u[i   + (j+1)*nx]        Vx_R = u[i+1 + (j+1)*nx]
//   Vz_B = v[i+1 +  j   *nxvz]      Vz_T = v[i+1 + (j+1)*nxvz]
//   P_C  = p[i   +  j   *ncx]
//
// The five corresponding rows of the MDOODZ matrix are extracted via the
// same `*_fill_coeffs` helpers that `ApplyStokesOperatorMDOODZ` uses, so
// by construction the block's 5×5 sub-matrix plus the external `out_b[r]`
// contributions match the assembled-matrix pattern bit-for-bit (same
// correctness contract as D11: enforced by
// `TESTS/StokesMatvecEquivalence.cpp`'s new block-cross-check).
//
// The output is raw (no celvol) to match the `L->rhs_u` scaling produced
// by `BuildLevelRHSFromCompressed` in `MultigridStokes.c`. Callers that
// want celvol-scaled output can multiply after the fact.
//
// For a DOF that is NOT a real equation (Dirichlet/Neumann velocity
// tag 0/11/13/31, inactive tag 30, periodic ghost -12), the block row is
// the identity and `out_b[r]` is the BC value the caller should pin the
// DOF to. For a DOF that IS a real equation, `out_b[r]` is Σ coeff ·
// state for neighbours outside the block (matching the Vanka pattern
// `M·dx = rhs - M_ext·x_ext`). Dirichlet-neighbour contributions inside
// the block are zeroed in the column (they already contribute to
// `bbc`, which is baked into the caller's `rhs_u`).
//
// Returns 0 on success, non-zero if (i, j) is outside the dispatch range
// (i ∈ [0, Ncx), j ∈ [0, Ncz) is the valid domain; we require i ≤ Ncx-1
// and j ≤ Ncz-1).
// ---------------------------------------------------------------------------

// Whether a velocity BC tag corresponds to a real equation row in the
// assembled matrix (i.e., the outer dispatch loop in
// Xmomentum/Zmomentum_InnerNodesDecoupled would have kept this cell).
static inline bool gmg_vel_is_real_eqn(char bc) {
    return (bc == -1) || (bc == 2) || (bc == -2);
}
static inline bool gmg_prs_is_real_eqn(char bc) {
    return (bc == -1);
}

// Slot indices inside the 5×5 block.
enum { GMG_SLOT_VxL = 0, GMG_SLOT_VxR, GMG_SLOT_VzB, GMG_SLOT_VzT, GMG_SLOT_P };

// ---------------------------------------------------------------------------
// Design D7 — symmetric smoothing rule (task 15.21).
//
// The Vanka smoother is stability-critical: V-cycle convergence requires
// a symmetric, spectral-equivalent smoother on the fine level. The full
// Newton Jacobian is NOT symmetric — it has off-diagonal cross terms
// coming from the rotation coupling (D13/D14 and D31/D32/D34). Using
// the Jacobian directly in the 5×5 block would destroy the symmetric-
// positive-definite property the smoother relies on.
//
// Design decision D7 resolves this by splitting roles:
//   - The matvec (ApplyStokesOperatorMDOODZ): uses the full Newton
//     Jacobian when model.Newton == 1. This drives the FGMRES residual
//     and gives quadratic nonlinear convergence.
//   - The Vanka block (StokesCellBlockMDOODZ): ALWAYS returns the
//     symmetric Picard part (xmom_fill_coeffs / zmom_fill_coeffs), even
//     in Newton mode. The smoother thus preserves SPD structure.
//
// This is the standard "Newton-in-matvec, Picard-in-smoother" pattern
// for nonlinear Stokes GMG. For consistency with the matvec's Newton
// branch (which forces comp = 1), we mirror the same override here so
// the block's diagonal reflects the same compressibility setting.
// ---------------------------------------------------------------------------
int StokesCellBlockMDOODZ(grid *mesh, params model,
                          const double *u, const double *v, const double *p,
                          int i, int j,
                          double out_A[5][5], double out_b[5]) {
    if (mesh == NULL || u == NULL || v == NULL || p == NULL) return -1;
    if (out_A == NULL || out_b == NULL) return -1;

    const int nx = mesh->Nx, nz = mesh->Nz;
    const int nxvz = nx + 1, ncx = nx - 1, ncz = nz - 1;
    if (i < 0 || i >= ncx) return -1;
    if (j < 0 || j >= ncz) return -1;

    const double one_dx = 1.0 / mesh->dx;
    const double one_dz = 1.0 / mesh->dz;

    int comp = 0, stab = 0;
    if (model.free_surface_stab > 0) stab = 1;
    // Mirror ApplyStokesOperatorMDOODZ's Newton-mode override: forces
    // comp = 1 so the block's compressibility setting matches the matvec
    // the smoother is being fed residuals from. Keeps the smoother's
    // fixed-point operator spectrally compatible with the matvec.
    if (model.Newton == 1) {
        comp = 1;
    } else if (model.compressible == 1) {
        comp = 1;
    }
    const double theta = model.free_surface_stab;

    // Full-grid indices of the five block DOFs (MDOODZ padded layout).
    const int iVxL = i     + (j + 1) * nx;        // u-index for Vx_L
    const int iVxR = (i+1) + (j + 1) * nx;        // u-index for Vx_R
    const int iVzB = (i+1) + (j    ) * nxvz;      // v-index for Vz_B
    const int iVzT = (i+1) + (j + 1) * nxvz;      // v-index for Vz_T
    const int iPC  = i     + (j    ) * ncx;       // p-index for P_C

    const char tagVxL = mesh->BCu.type[iVxL];
    const char tagVxR = mesh->BCu.type[iVxR];
    const char tagVzB = mesh->BCv.type[iVzB];
    const char tagVzT = mesh->BCv.type[iVzT];
    const char tagPC  = mesh->BCp.type[iPC];

    const bool realVxL = gmg_vel_is_real_eqn(tagVxL);
    const bool realVxR = gmg_vel_is_real_eqn(tagVxR);
    const bool realVzB = gmg_vel_is_real_eqn(tagVzB);
    const bool realVzT = gmg_vel_is_real_eqn(tagVzT);
    const bool realPC  = gmg_prs_is_real_eqn(tagPC);
    const bool slot_real[5] = { realVxL, realVxR, realVzB, realVzT, realPC };
    const int  slot_dof  [5] = { iVxL, iVxR, iVzB, iVzT, iPC };
    const double slot_bc_val[5] = {
        mesh->BCu.val[iVxL], mesh->BCu.val[iVxR],
        mesh->BCv.val[iVzB], mesh->BCv.val[iVzT],
        mesh->BCp.val[iPC],
    };

    // Initialise block to zero; out_b to zero.
    for (int r = 0; r < 5; ++r) {
        for (int c = 0; c < 5; ++c) out_A[r][c] = 0.0;
        out_b[r] = 0.0;
    }

    // Convenience macro for dispatching an arbitrary (idx, kind, coeff,
    // Set, BCtag_array, state_array) tuple to the right block slot or
    // external bucket. `kind` is 0=Vx, 1=Vz, 2=P.
    //   - If Set && !gmg_skip_neighbour(tag): either put into block[r][c]
    //     for the matching slot, or accumulate coeff*state into out_b[r].
    //   - If !Set: skip entirely (tag-30 / outside-stencil).
    //   - If Set && gmg_skip_neighbour: skip entirely; the bbc goes into
    //     the caller's RHS already.
    //
    // For the inside-block mapping we need: is this (kind, idx) one of
    // the five local DOFs? We compare against the known slot_dof[]. For
    // a block-internal neighbour whose own block DOF is Dirichlet (not a
    // real equation), we STILL want the coefficient to appear in
    // block[r][c] — the identity-row decoupling below removes that
    // column consistently across all rows so the 5×5 solve treats the
    // Dirichlet DOF as "known = bc value". That matches the textbook
    // VankaBlockAssembleSolve pattern and keeps the block non-singular.

    // --- X-mom row at Vx_L (only if Vx_L carries an equation) ---
    if (realVxL) {
        XmomCoeffs XL;
        const int c1 = iVxL;                    // k = i, l = j+1
        const int c2 = (i - 1) + j * ncx;
        const int c3 = i + (j + 1) * nxvz;
        xmom_fill_coeffs(mesh, model, /*ix=*/i, /*iz=*/j + 1,
                         stab, comp, theta, one_dx, one_dz,
                         c1, c2, c3, nx, nz, &XL);
        // Inside-block: iVxC -> slot 0, iVxE -> slot 1, iVzSE -> slot 2,
        // iVzNE -> slot 3, iPrE -> slot 4.
        // Outside: iVxS / iVxW / iVxN / iVzSW / iVzNW / iPrW.
        const int r = GMG_SLOT_VxL;
        signed char *BCu = mesh->BCu.type, *BCv = mesh->BCv.type, *BCp = mesh->BCp.type;

        // Vx_C = Vx_L (slot 0)
        if (XL.SetVxC) out_A[r][GMG_SLOT_VxL] += XL.uC;  // always real (= this row's DOF)

        // Vx_E = Vx_R (slot 1)
        if (XL.SetVxE && !gmg_skip_neighbour(BCu[XL.iVxE])) out_A[r][GMG_SLOT_VxR] += XL.uE;

        // Vz_SE = Vz_B (slot 2)
        if (XL.SetVzSE && !gmg_skip_neighbour(BCv[XL.iVzSE])) out_A[r][GMG_SLOT_VzB] += XL.vSE;

        // Vz_NE = Vz_T (slot 3)
        if (XL.SetVzNE && !gmg_skip_neighbour(BCv[XL.iVzNE])) out_A[r][GMG_SLOT_VzT] += XL.vNE;

        // Pr_E = P_C (slot 4)
        if (XL.SetPrE && !gmg_skip_neighbour(BCp[XL.iPrE])) out_A[r][GMG_SLOT_P] += XL.pE;

        // Outside contributions.
        if (XL.SetVxS  && !gmg_skip_neighbour(BCu[XL.iVxS]))  out_b[r] += XL.uS  * u[XL.iVxS];
        if (XL.SetVxW  && !gmg_skip_neighbour(BCu[XL.iVxW]))  out_b[r] += XL.uW  * u[XL.iVxW];
        if (XL.SetVxN  && !gmg_skip_neighbour(BCu[XL.iVxN]))  out_b[r] += XL.uN  * u[XL.iVxN];
        if (XL.SetVzSW && !gmg_skip_neighbour(BCv[XL.iVzSW])) out_b[r] += XL.vSW * v[XL.iVzSW];
        if (XL.SetVzNW && !gmg_skip_neighbour(BCv[XL.iVzNW])) out_b[r] += XL.vNW * v[XL.iVzNW];
        if (XL.SetPrW  && !gmg_skip_neighbour(BCp[XL.iPrW]))  out_b[r] += XL.pW  * p[XL.iPrW];
    }

    // --- X-mom row at Vx_R ---
    if (realVxR) {
        XmomCoeffs XR;
        const int c1 = iVxR;                    // k = i+1, l = j+1
        const int c2 = i + j * ncx;
        const int c3 = (i + 1) + (j + 1) * nxvz;
        xmom_fill_coeffs(mesh, model, /*ix=*/i + 1, /*iz=*/j + 1,
                         stab, comp, theta, one_dx, one_dz,
                         c1, c2, c3, nx, nz, &XR);
        const int r = GMG_SLOT_VxR;
        signed char *BCu = mesh->BCu.type, *BCv = mesh->BCv.type, *BCp = mesh->BCp.type;

        // Vx_W = Vx_L (slot 0)
        if (XR.SetVxW && !gmg_skip_neighbour(BCu[XR.iVxW])) out_A[r][GMG_SLOT_VxL] += XR.uW;
        // Vx_C = Vx_R (slot 1)
        if (XR.SetVxC) out_A[r][GMG_SLOT_VxR] += XR.uC;
        // Vz_SW = Vz_B (slot 2)
        if (XR.SetVzSW && !gmg_skip_neighbour(BCv[XR.iVzSW])) out_A[r][GMG_SLOT_VzB] += XR.vSW;
        // Vz_NW = Vz_T (slot 3)
        if (XR.SetVzNW && !gmg_skip_neighbour(BCv[XR.iVzNW])) out_A[r][GMG_SLOT_VzT] += XR.vNW;
        // Pr_W = P_C (slot 4)
        if (XR.SetPrW && !gmg_skip_neighbour(BCp[XR.iPrW])) out_A[r][GMG_SLOT_P] += XR.pW;

        // Outside.
        if (XR.SetVxS  && !gmg_skip_neighbour(BCu[XR.iVxS]))  out_b[r] += XR.uS  * u[XR.iVxS];
        if (XR.SetVxE  && !gmg_skip_neighbour(BCu[XR.iVxE]))  out_b[r] += XR.uE  * u[XR.iVxE];
        if (XR.SetVxN  && !gmg_skip_neighbour(BCu[XR.iVxN]))  out_b[r] += XR.uN  * u[XR.iVxN];
        if (XR.SetVzSE && !gmg_skip_neighbour(BCv[XR.iVzSE])) out_b[r] += XR.vSE * v[XR.iVzSE];
        if (XR.SetVzNE && !gmg_skip_neighbour(BCv[XR.iVzNE])) out_b[r] += XR.vNE * v[XR.iVzNE];
        if (XR.SetPrE  && !gmg_skip_neighbour(BCp[XR.iPrE]))  out_b[r] += XR.pE  * p[XR.iPrE];
    }

    // --- Z-mom row at Vz_B ---
    if (realVzB) {
        ZmomCoeffs ZB;
        const int c3 = iVzB;                   // k = i+1, l = j
        const int c1 = (i + 1) + j * nx;
        const int c2 = i + (j - 1) * ncx;      // == (k-1) + (l-1)*ncx
        zmom_fill_coeffs(mesh, model, /*ix=*/i + 1, /*iz=*/j,
                         stab, comp, theta, one_dx, one_dz,
                         c1, c2, c3, nx, nz, &ZB);
        const int r = GMG_SLOT_VzB;
        signed char *BCu = mesh->BCu.type, *BCv = mesh->BCv.type, *BCp = mesh->BCp.type;

        // Vx_NW = Vx_L (slot 0)
        if (ZB.SetVxNW && !gmg_skip_neighbour(BCu[ZB.iVxNW])) out_A[r][GMG_SLOT_VxL] += ZB.uNW;
        // Vx_NE = Vx_R (slot 1)
        if (ZB.SetVxNE && !gmg_skip_neighbour(BCu[ZB.iVxNE])) out_A[r][GMG_SLOT_VxR] += ZB.uNE;
        // Vz_C = Vz_B (slot 2)
        if (ZB.SetVzC) out_A[r][GMG_SLOT_VzB] += ZB.vC;
        // Vz_N = Vz_T (slot 3)
        if (ZB.SetVzN && !gmg_skip_neighbour(BCv[ZB.iVzN])) out_A[r][GMG_SLOT_VzT] += ZB.vN;
        // Pr_N = P_C (slot 4)
        if (ZB.SetPrN && !gmg_skip_neighbour(BCp[ZB.iPrN])) out_A[r][GMG_SLOT_P] += ZB.pN;

        // Outside.
        if (ZB.SetVxSW && !gmg_skip_neighbour(BCu[ZB.iVxSW])) out_b[r] += ZB.uSW * u[ZB.iVxSW];
        if (ZB.SetVxSE && !gmg_skip_neighbour(BCu[ZB.iVxSE])) out_b[r] += ZB.uSE * u[ZB.iVxSE];
        if (ZB.SetVzS  && !gmg_skip_neighbour(BCv[ZB.iVzS]))  out_b[r] += ZB.vS  * v[ZB.iVzS];
        if (ZB.SetVzW  && !gmg_skip_neighbour(BCv[ZB.iVzW]))  out_b[r] += ZB.vW  * v[ZB.iVzW];
        if (ZB.SetVzE  && !gmg_skip_neighbour(BCv[ZB.iVzE]))  out_b[r] += ZB.vE  * v[ZB.iVzE];
        if (ZB.SetPrS  && !gmg_skip_neighbour(BCp[ZB.iPrS]))  out_b[r] += ZB.pS  * p[ZB.iPrS];
    }

    // --- Z-mom row at Vz_T ---
    if (realVzT) {
        ZmomCoeffs ZT;
        const int c3 = iVzT;                   // k = i+1, l = j+1
        const int c1 = (i + 1) + (j + 1) * nx;
        const int c2 = i + j * ncx;
        zmom_fill_coeffs(mesh, model, /*ix=*/i + 1, /*iz=*/j + 1,
                         stab, comp, theta, one_dx, one_dz,
                         c1, c2, c3, nx, nz, &ZT);
        const int r = GMG_SLOT_VzT;
        signed char *BCu = mesh->BCu.type, *BCv = mesh->BCv.type, *BCp = mesh->BCp.type;

        // Vx_SW = Vx_L (slot 0)
        if (ZT.SetVxSW && !gmg_skip_neighbour(BCu[ZT.iVxSW])) out_A[r][GMG_SLOT_VxL] += ZT.uSW;
        // Vx_SE = Vx_R (slot 1)
        if (ZT.SetVxSE && !gmg_skip_neighbour(BCu[ZT.iVxSE])) out_A[r][GMG_SLOT_VxR] += ZT.uSE;
        // Vz_S = Vz_B (slot 2)
        if (ZT.SetVzS && !gmg_skip_neighbour(BCv[ZT.iVzS])) out_A[r][GMG_SLOT_VzB] += ZT.vS;
        // Vz_C = Vz_T (slot 3)
        if (ZT.SetVzC) out_A[r][GMG_SLOT_VzT] += ZT.vC;
        // Pr_S = P_C (slot 4)
        if (ZT.SetPrS && !gmg_skip_neighbour(BCp[ZT.iPrS])) out_A[r][GMG_SLOT_P] += ZT.pS;

        // Outside.
        if (ZT.SetVxNW && !gmg_skip_neighbour(BCu[ZT.iVxNW])) out_b[r] += ZT.uNW * u[ZT.iVxNW];
        if (ZT.SetVxNE && !gmg_skip_neighbour(BCu[ZT.iVxNE])) out_b[r] += ZT.uNE * u[ZT.iVxNE];
        if (ZT.SetVzN  && !gmg_skip_neighbour(BCv[ZT.iVzN]))  out_b[r] += ZT.vN  * v[ZT.iVzN];
        if (ZT.SetVzW  && !gmg_skip_neighbour(BCv[ZT.iVzW]))  out_b[r] += ZT.vW  * v[ZT.iVzW];
        if (ZT.SetVzE  && !gmg_skip_neighbour(BCv[ZT.iVzE]))  out_b[r] += ZT.vE  * v[ZT.iVzE];
        if (ZT.SetPrN  && !gmg_skip_neighbour(BCp[ZT.iPrN]))  out_b[r] += ZT.pN  * p[ZT.iPrN];
    }

    // --- Continuity row at P_C ---
    if (realPC) {
        ContCoeffs CC;
        const int c1 = i + (j + 1) * nx;
        const int c3 = (i + 1) + j * nxvz;
        int cell_comp = (comp == 1) ? mesh->comp_cells[iPC] : 0;
        cont_fill_coeffs(mesh, model, cell_comp, one_dx, one_dz, c1, c3, nxvz, &CC);
        const int r = GMG_SLOT_P;
        signed char *BCu = mesh->BCu.type, *BCv = mesh->BCv.type;

        // Vx_W = Vx_L (slot 0)
        if (CC.SetVxW && !gmg_skip_neighbour(BCu[CC.iVxW])) out_A[r][GMG_SLOT_VxL] += CC.uW;
        // Vx_E = Vx_R (slot 1)
        if (CC.SetVxE && !gmg_skip_neighbour(BCu[CC.iVxE])) out_A[r][GMG_SLOT_VxR] += CC.uE;
        // Vz_S = Vz_B (slot 2)
        if (CC.SetVzS && !gmg_skip_neighbour(BCv[CC.iVzS])) out_A[r][GMG_SLOT_VzB] += CC.vS;
        // Vz_N = Vz_T (slot 3)
        if (CC.SetVzN && !gmg_skip_neighbour(BCv[CC.iVzN])) out_A[r][GMG_SLOT_VzT] += CC.vN;
        // P_C diagonal — Picard continuity has no pressure self-coupling
        // (out_A[r][GMG_SLOT_P] left at 0). Compressibility β/dt lives in
        // StokesC->F, not StokesD; see `cont_fill_coeffs`'s block comment.
        // No outside contributions exist for the continuity stencil.
    }

    // --- Identity-row handling for DOFs that are NOT real equations. ---
    //
    // MDOODZ's outer assembly loop (Xmomentum/Zmomentum/Continuity
    // dispatch in StokesAssemblyDecoupled.c) skips cells whose own BC
    // tag is not a real-equation tag (velocity: -1, 2, -2; pressure:
    // -1), meaning the global matrix has no row for those DOFs. At the
    // Vanka block level we install an identity row: the block solve
    // then returns the Dirichlet value exactly (= the current iterate,
    // which mesh->BC*.val pins during SetBCs). We also zero the column
    // of those Dirichlet DOFs in the OTHER rows — matching the
    // Assemble path's AddCoeff3 behaviour of routing Dirichlet-neighbour
    // contributions to `bbc` rather than into the matrix — but we do
    // this AFTER the coefficient-projection pass above, so the
    // non-gmg-skip Dirichlet case (unlikely but possible for tag 30)
    // stays consistent.
    for (int r = 0; r < 5; ++r) {
        if (slot_real[r]) continue;
        for (int c = 0; c < 5; ++c) out_A[r][c] = (c == r) ? 1.0 : 0.0;
        // Output b[r] = BC value; caller uses `rhs[r] = out_b[r]` for
        // identity rows rather than the rhs_global - out_b formula.
        // Slot 4 (pressure) is never a Dirichlet DOF in production mesh
        // setups (MDOODZ pins P via the null-space projection, not via
        // a BCp.val), but we still honour the convention.
        out_b[r] = slot_bc_val[r];
    }
    // Zero the column of each identity row in the OTHER rows. For
    // AddCoeff3-skipped neighbour tags (0/11/13/31) the coefficient
    // projection above already skipped the contribution; the only way
    // we would have a nonzero `out_A[rp][r]` for a non-real-eqn row r
    // is if the neighbour's tag was e.g. tag 30 (inactive) which the
    // Set* predicate already filters. So this loop is defensive — it
    // makes the block structurally clean and the identity-row
    // decoupling explicit, matching `VankaBlockAssembleSolve`'s textbook
    // pattern.
    for (int r = 0; r < 5; ++r) {
        if (slot_real[r]) continue;
        for (int rp = 0; rp < 5; ++rp) {
            if (rp == r) continue;
            out_A[rp][r] = 0.0;
        }
    }

    (void)slot_dof;
    return 0;
}
