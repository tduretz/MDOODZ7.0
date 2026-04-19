// Golden cross-check test for the GMG stencil bridge (design D11, task 15.9–15.13).
//
// Builds a small representative mesh, assembles the full sparse Stokes
// operator via BuildStokesOperatorDecoupled, and compares A·x against
// ApplyStokesOperatorMDOODZ row-by-row to within 1e-12 relative tolerance.
// Any drift between StokesAssemblyDecoupled.c and StokesAssemblyGMG.c
// coefficient formulas will fail this test loudly.

#include <cmath>
#include <complex>
#include <cstdint>
#include <cstring>
#include <gtest/gtest.h>
#include <random>
#include <string>
#include <vector>

extern "C" {
#include "mdoodz-private.h"
#include "StokesAssemblyGMG.h"
void BuildStokesOperatorDecoupled(grid *, params, int, double *, double *,
                                  double *, double *, SparseMat *, SparseMat *,
                                  SparseMat *, SparseMat *, SparseMat *, int);
void BuildJacobianOperatorDecoupled(grid *, params, int, double *, double *,
                                    double *, double *, SparseMat *, SparseMat *,
                                    SparseMat *, SparseMat *, SparseMat *, int);
void EvalNumberOfEquations(grid *, SparseMat *);
void AllocMat(SparseMat *, int);
void FreeMat(SparseMat *);
}

namespace {

// Holder for a self-contained 2-D mesh sufficient to drive
// BuildStokesOperatorDecoupled. All pointers inside `grid` that the
// assembler touches are backed by std::vectors in this struct so that
// memory lives exactly as long as the test body.
struct MeshFixture {
    int Nx, Nz;
    double dx, dz;

    std::vector<signed char> BCu_type, BCv_type, BCp_type, BCg_type;
    std::vector<double>      BCu_val,  BCv_val,  BCp_val;

    std::vector<double> D11_n, D22_n, D33_s;
    // Full D-tensor components used by the Newton Jacobian matvec
    // (Xjacobian_InnerNodesDecoupled3 / Zjacobian_InnerNodesDecoupled3).
    // The Picard path only touches D11_n / D22_n / D33_s; the extras
    // remain zero in that case and the read is harmless.
    std::vector<double> D12_n, D13_n, D14_n;
    std::vector<double> D21_n, D23_n, D24_n;
    std::vector<double> D31_s, D32_s, D34_s;
    std::vector<double> rho_n, gx, gz, bet_n;
    std::vector<int>    comp_cells;
    std::vector<int>    kvx, lvx, kvz, lvz, kp, lp;
    std::vector<double> roger_x, roger_z, rhs_p;
    std::vector<double> sxxd, szzd, sxz, p_corr;
    std::vector<double> u_in, v_in, p_in;

    grid mesh{};
    params model{};

    void build(int Nx_, int Nz_, double Lx = 1.0, double Lz = 1.0) {
        Nx = Nx_;
        Nz = Nz_;
        const int Ncx = Nx - 1, Ncz = Nz - 1, Nxvz = Nx + 1, Nzvx = Nz + 1;
        dx = Lx / Ncx;
        dz = Lz / Ncz;

        std::memset(&mesh, 0, sizeof(mesh));
        mesh.Nx = Nx;
        mesh.Nz = Nz;
        mesh.dx = dx;
        mesh.dz = dz;

        BCu_type.assign((size_t)Nx * Nzvx, -1);
        BCv_type.assign((size_t)Nxvz * Nz, -1);
        BCp_type.assign((size_t)Ncx * Ncz, -1);
        BCg_type.assign((size_t)Nx * Nz,   -1);
        BCu_val .assign((size_t)Nx * Nzvx, 0.0);
        BCv_val .assign((size_t)Nxvz * Nz, 0.0);
        BCp_val .assign((size_t)Ncx * Ncz, 0.0);

        // Closed-box Dirichlet default.
        for (int l = 0; l < Nzvx; ++l) {
            BCu_type[(size_t)0        + (size_t)l * Nx] = 0;
            BCu_type[(size_t)(Nx - 1) + (size_t)l * Nx] = 0;
        }
        for (int k = 0; k < Nx; ++k) {
            BCu_type[(size_t)k + (size_t)0           * Nx] = 11;
            BCu_type[(size_t)k + (size_t)(Nzvx - 1)  * Nx] = 11;
        }
        for (int l = 0; l < Nz; ++l) {
            BCv_type[(size_t)0           + (size_t)l * Nxvz] = 11;
            BCv_type[(size_t)(Nxvz - 1)  + (size_t)l * Nxvz] = 11;
        }
        for (int k = 0; k < Nxvz; ++k) {
            BCv_type[(size_t)k + (size_t)0          * Nxvz] = 0;
            BCv_type[(size_t)k + (size_t)(Nz - 1)   * Nxvz] = 0;
        }

        D11_n.assign((size_t)Ncx * Ncz, 2.0);
        D22_n.assign((size_t)Ncx * Ncz, 2.0);
        D33_s.assign((size_t)Nx  * Nz,  1.0);
        // Newton-only D-tensor entries. Default to zero (the cross-
        // coupling anisotropic / compressible-Jacobian contributions
        // reduce to the Picard stencil when these are zero), so the
        // Picard golden tests are unaffected. Newton-mode tests will
        // set these explicitly before invoking the Newton path.
        D12_n.assign((size_t)Ncx * Ncz, 0.0);
        D13_n.assign((size_t)Ncx * Ncz, 0.0);
        D14_n.assign((size_t)Ncx * Ncz, 0.0);
        D21_n.assign((size_t)Ncx * Ncz, 0.0);
        D23_n.assign((size_t)Ncx * Ncz, 0.0);
        D24_n.assign((size_t)Ncx * Ncz, 0.0);
        D31_s.assign((size_t)Nx  * Nz,  0.0);
        D32_s.assign((size_t)Nx  * Nz,  0.0);
        D34_s.assign((size_t)Nx  * Nz,  0.0);
        rho_n.assign((size_t)Ncx * Ncz, 0.0);
        gx   .assign((size_t)Nx  * Nzvx, 0.0);
        gz   .assign((size_t)Nxvz * Nz, 0.0);
        bet_n.assign((size_t)Ncx * Ncz, 0.0);
        comp_cells.assign((size_t)Ncx * Ncz, 0);

        roger_x.assign((size_t)Nx * Nzvx, 0.0);
        roger_z.assign((size_t)Nxvz * Nz, 0.0);
        rhs_p  .assign((size_t)Ncx * Ncz, 0.0);
        sxxd   .assign((size_t)Ncx * Ncz, 0.0);
        szzd   .assign((size_t)Ncx * Ncz, 0.0);
        sxz    .assign((size_t)Nx  * Nz,  0.0);
        p_corr .assign((size_t)Ncx * Ncz, 0.0);
        u_in   .assign((size_t)Nx  * Nzvx, 0.0);
        v_in   .assign((size_t)Nxvz * Nz, 0.0);
        p_in   .assign((size_t)Ncx * Ncz, 0.0);

        kvx.assign((size_t)Nx  * Nzvx, 0);
        lvx.assign((size_t)Nx  * Nzvx, 0);
        kvz.assign((size_t)Nxvz * Nz,  0);
        lvz.assign((size_t)Nxvz * Nz,  0);
        kp .assign((size_t)Ncx * Ncz,  0);
        lp .assign((size_t)Ncx * Ncz,  0);
        for (int l = 0; l < Nzvx; ++l)
            for (int k = 0; k < Nx; ++k) {
                kvx[(size_t)k + (size_t)l * Nx] = k;
                lvx[(size_t)k + (size_t)l * Nx] = l;
            }
        for (int l = 0; l < Nz; ++l)
            for (int k = 0; k < Nxvz; ++k) {
                kvz[(size_t)k + (size_t)l * Nxvz] = k;
                lvz[(size_t)k + (size_t)l * Nxvz] = l;
            }
        for (int l = 0; l < Ncz; ++l)
            for (int k = 0; k < Ncx; ++k) {
                kp[(size_t)k + (size_t)l * Ncx] = k;
                lp[(size_t)k + (size_t)l * Ncx] = l;
            }

        mesh.BCu.type = BCu_type.data();
        mesh.BCv.type = BCv_type.data();
        mesh.BCp.type = BCp_type.data();
        mesh.BCg.type = BCg_type.data();
        mesh.BCu.val  = BCu_val.data();
        mesh.BCv.val  = BCv_val.data();
        mesh.BCp.val  = BCp_val.data();
        mesh.D11_n    = D11_n.data();
        mesh.D22_n    = D22_n.data();
        mesh.D33_s    = D33_s.data();
        mesh.D12_n    = D12_n.data();
        mesh.D13_n    = D13_n.data();
        mesh.D14_n    = D14_n.data();
        mesh.D21_n    = D21_n.data();
        mesh.D23_n    = D23_n.data();
        mesh.D24_n    = D24_n.data();
        mesh.D31_s    = D31_s.data();
        mesh.D32_s    = D32_s.data();
        mesh.D34_s    = D34_s.data();
        mesh.rho_n    = rho_n.data();
        mesh.gx       = gx.data();
        mesh.gz       = gz.data();
        mesh.bet_n    = bet_n.data();
        mesh.comp_cells = comp_cells.data();
        mesh.kvx = kvx.data();
        mesh.lvx = lvx.data();
        mesh.kvz = kvz.data();
        mesh.lvz = lvz.data();
        mesh.kp  = kp .data();
        mesh.lp  = lp .data();
        mesh.roger_x = roger_x.data();
        mesh.roger_z = roger_z.data();
        mesh.rhs_p   = rhs_p  .data();
        mesh.sxxd    = sxxd   .data();
        mesh.szzd    = szzd   .data();
        mesh.sxz     = sxz    .data();
        mesh.p_corr  = p_corr .data();
        mesh.u_in    = u_in   .data();
        mesh.v_in    = v_in   .data();
        mesh.p_in    = p_in   .data();

        std::memset(&model, 0, sizeof(model));
        model.Nx = Nx;
        model.Nz = Nz;
        model.dx = dx;
        model.dz = dz;
        model.compressible      = 0;
        model.out_of_plane      = 0;
        model.free_surface_stab = 0.0;
        model.dt                = 1.0;
        model.Newton            = 0;
        model.periodic_x        = 0;
        model.writer_debug      = 0;
    }
};

// Copy full-grid arrays into the compressed Stokes `x` vector used by
// SparseMat multiplication (size Stokes->neq).
static void PackFullToCompressed(const grid *mesh, const SparseMat *Stokes,
                                 const double *u_full, const double *v_full,
                                 const double *p_full,
                                 std::vector<double> &x) {
    const int Nx = mesh->Nx, Nz = mesh->Nz;
    const int Ncx = Nx - 1, Ncz = Nz - 1, Nxvz = Nx + 1, Nzvx = Nz + 1;
    x.assign((size_t)Stokes->neq, 0.0);
    for (int l = 0; l < Nzvx; ++l)
        for (int k = 0; k < Nx; ++k) {
            int kk  = k + l * Nx;
            int eqn = Stokes->eqn_u[kk];
            if (eqn < 0) continue;
            if (mesh->BCu.type[kk] == -1 || mesh->BCu.type[kk] == 2 ||
                mesh->BCu.type[kk] == -2)
                x[eqn] = u_full[kk];
        }
    for (int l = 0; l < Nz; ++l)
        for (int k = 0; k < Nxvz; ++k) {
            int kk  = k + l * Nxvz;
            int eqn = Stokes->eqn_v[kk];
            if (eqn < 0) continue;
            if (mesh->BCv.type[kk] == -1 || mesh->BCv.type[kk] == 2)
                x[eqn] = v_full[kk];
        }
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            int kk  = k + l * Ncx;
            int eqn = Stokes->eqn_p[kk];
            if (eqn < 0) continue;
            if (mesh->BCp.type[kk] == -1)
                x[eqn] = p_full[kk];
        }
}

static void CsrMatvec(const SparseMat *M, const std::vector<double> &x,
                      std::vector<double> &y) {
    y.assign((size_t)M->neq, 0.0);
    for (int i = 0; i < M->neq; ++i) {
        int row_start = M->Ic[i];
        int row_end   = M->Ic[i + 1];
        double s = 0.0;
        for (int k = row_start; k < row_end; ++k) s += M->A[k] * x[M->J[k]];
        y[i] = s;
    }
}

// Core golden-check driver: assembles the sparse operator from the fixture,
// performs ten random matvec comparisons, and asserts relative error < tol.
// Returns the worst observed relative error (useful for diagnostics).
static double RunGoldenCheck(MeshFixture &F, const std::string &label,
                             double tol = 1e-12, int n_trials = 10) {
    SparseMat Stokes{}, StokesA{}, StokesB{}, StokesC{}, StokesD{};
    EvalNumberOfEquations(&F.mesh, &Stokes);
    StokesA.eqn_u = StokesB.eqn_u = StokesC.eqn_u = StokesD.eqn_u = Stokes.eqn_u;
    StokesA.eqn_v = StokesB.eqn_v = StokesC.eqn_v = StokesD.eqn_v = Stokes.eqn_v;
    StokesA.eqn_p = StokesB.eqn_p = StokesC.eqn_p = StokesD.eqn_p = Stokes.eqn_p;
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    std::vector<double> bbcA((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bbcB((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bbcC((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bbcD((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bA  ((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bB  ((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bC  ((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bD  ((size_t)Stokes.neq_cont, 0.0);
    StokesA.bbc = bbcA.data(); StokesB.bbc = bbcB.data();
    StokesC.bbc = bbcC.data(); StokesD.bbc = bbcD.data();
    StokesA.b = bA.data(); StokesB.b = bB.data();
    StokesC.b = bC.data(); StokesD.b = bD.data();

    BuildStokesOperatorDecoupled(&F.mesh, F.model, 0,
                                 F.p_corr.data(), F.p_in.data(),
                                 F.u_in.data(),   F.v_in.data(),
                                 &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1);

    std::mt19937 rng(42u + std::hash<std::string>{}(label));
    std::uniform_real_distribution<double> U(-1.0, 1.0);

    double worst_rel = 0.0;
    for (int trial = 0; trial < n_trials; ++trial) {
        std::vector<double> u_full(F.u_in.size(), 0.0);
        std::vector<double> v_full(F.v_in.size(), 0.0);
        std::vector<double> p_full(F.p_in.size(), 0.0);
        for (size_t k = 0; k < u_full.size(); ++k) u_full[k] = U(rng);
        for (size_t k = 0; k < v_full.size(); ++k) v_full[k] = U(rng);
        for (size_t k = 0; k < p_full.size(); ++k) p_full[k] = U(rng);

        for (size_t k = 0; k < u_full.size(); ++k)
            if (F.BCu_type[k] != -1 && F.BCu_type[k] != 2 && F.BCu_type[k] != -2)
                u_full[k] = 0.0;
        for (size_t k = 0; k < v_full.size(); ++k)
            if (F.BCv_type[k] != -1 && F.BCv_type[k] != 2)
                v_full[k] = 0.0;
        for (size_t k = 0; k < p_full.size(); ++k)
            if (F.BCp_type[k] != -1) p_full[k] = 0.0;

        std::vector<double> x;
        PackFullToCompressed(&F.mesh, &Stokes, u_full.data(), v_full.data(),
                             p_full.data(), x);
        std::vector<double> x_u(Stokes.neq_mom), x_p(Stokes.neq_cont);
        for (int i = 0; i < Stokes.neq_mom;  ++i) x_u[i] = x[i];
        for (int i = 0; i < Stokes.neq_cont; ++i) x_p[i] = x[Stokes.neq_mom + i];

        std::vector<double> Ax_u, Bp_u, Cx_p, Dp_p;
        CsrMatvec(&StokesA, x_u, Ax_u);
        CsrMatvec(&StokesB, x_p, Bp_u);
        CsrMatvec(&StokesC, x_u, Cx_p);
        CsrMatvec(&StokesD, x_p, Dp_p);

        std::vector<double> ru(F.u_in.size(), 0.0);
        std::vector<double> rv(F.v_in.size(), 0.0);
        std::vector<double> rp(F.p_in.size(), 0.0);
        ApplyStokesOperatorMDOODZ(&F.mesh, F.model,
                                  u_full.data(), v_full.data(), p_full.data(),
                                  ru.data(),     rv.data(),     rp.data());

        double max_abs = 0.0, max_err = 0.0;
        for (int l = 0; l < F.Nz + 1; ++l)
            for (int k = 0; k < F.Nx; ++k) {
                int kk  = k + l * F.Nx;
                int eqn = Stokes.eqn_u[kk];
                if (eqn < 0) continue;
                double expected = Ax_u[eqn] + Bp_u[eqn];
                double got      = ru[kk];
                double err      = std::fabs(expected - got);
                max_abs = std::max(max_abs, std::fabs(expected));
                max_err = std::max(max_err, err);
            }
        for (int l = 0; l < F.Nz; ++l)
            for (int k = 0; k < F.Nx + 1; ++k) {
                int kk  = k + l * (F.Nx + 1);
                int eqn = Stokes.eqn_v[kk];
                if (eqn < 0) continue;
                double expected = Ax_u[eqn] + Bp_u[eqn];
                double got      = rv[kk];
                double err      = std::fabs(expected - got);
                max_abs = std::max(max_abs, std::fabs(expected));
                max_err = std::max(max_err, err);
            }
        for (int l = 0; l < Ncz; ++l)
            for (int k = 0; k < Ncx; ++k) {
                int kk  = k + l * Ncx;
                int eqn = Stokes.eqn_p[kk];
                if (eqn < 0) continue;
                int p_row = eqn - Stokes.neq_mom;
                double expected = Cx_p[p_row] + Dp_p[p_row];
                double got      = rp[kk];
                double err      = std::fabs(expected - got);
                max_abs = std::max(max_abs, std::fabs(expected));
                max_err = std::max(max_err, err);
            }

        const double rel = (max_abs > 0.0) ? max_err / max_abs : max_err;
        worst_rel = std::max(worst_rel, rel);
        EXPECT_LT(rel, tol) << label << " trial " << trial
                            << " max_err=" << max_err
                            << " max_abs=" << max_abs;
    }

    FreeMat(&StokesA); FreeMat(&StokesB);
    FreeMat(&StokesC); FreeMat(&StokesD);
    DoodzFree(Stokes.eqn_u);
    DoodzFree(Stokes.eqn_v);
    DoodzFree(Stokes.eqn_p);
    return worst_rel;
}

}  // namespace

// ---------------------------------------------------------------------------
// Smoke test — constant-viscosity, closed-box Dirichlet (tags 0, 11).
// ---------------------------------------------------------------------------

TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledConstantViscosity) {
    MeshFixture F;
    F.build(21, 21);
    RunGoldenCheck(F, "constant_visc");
}

// ---------------------------------------------------------------------------
// BC-tag parameterisation (task 15.11).
// MDOODZ assigns BC types through SetBCs (user code); we simulate the
// production combinations by overriding BCu/BCv/BCp on the default closed
// box fixture and re-running the golden check.
// ---------------------------------------------------------------------------

// Tag-13 (Neumann velocity ghost) on top/bottom Vx ghost rows — exercises
// the RegN/RegS / NeuN/NeuS branches in xmom_row_matvec.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledNeumannVxGhost13) {
    MeshFixture F;
    F.build(21, 21);
    const int Nx = F.Nx, Nz = F.Nz, Nzvx = Nz + 1;
    for (int k = 1; k < Nx - 1; ++k) {
        F.BCu_type[(size_t)k + (size_t)0          * Nx] = 13;
        F.BCu_type[(size_t)k + (size_t)(Nzvx - 1) * Nx] = 13;
    }
    RunGoldenCheck(F, "neumann_vx_ghost_13");
}

// Tag-13 on left/right Vz ghost columns.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledNeumannVzGhost13) {
    MeshFixture F;
    F.build(21, 21);
    const int Nz = F.Nz, Nxvz = F.Nx + 1;
    for (int l = 1; l < Nz - 1; ++l) {
        F.BCv_type[(size_t)0          + (size_t)l * Nxvz] = 13;
        F.BCv_type[(size_t)(Nxvz - 1) + (size_t)l * Nxvz] = 13;
    }
    RunGoldenCheck(F, "neumann_vz_ghost_13");
}

// Tag-30 patch (deactivated/air cells) — exercises inW/inE/inS/inN gates in
// the momentum rows.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledAirPatch30) {
    MeshFixture F;
    F.build(21, 21);
    const int Nx = F.Nx, Nz = F.Nz, Ncx = Nx - 1, Ncz = Nz - 1;
    // Inner 3x3 patch of tag-30 cells in the upper half.
    for (int l = Ncz - 4; l < Ncz - 1; ++l)
        for (int k = Ncx / 2 - 1; k < Ncx / 2 + 2; ++k) {
            F.BCp_type[(size_t)k + (size_t)l * Ncx] = 30;
            F.BCg_type[(size_t)k + (size_t)l * Nx]  = 30;
        }
    RunGoldenCheck(F, "air_patch_30");
}

// ---------------------------------------------------------------------------
// Rheology parameterisation (task 15.12).
// ---------------------------------------------------------------------------

// Anisotropic D33 — vertex shear viscosity varies across the grid.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledAnisotropicD33) {
    MeshFixture F;
    F.build(21, 21);
    const int Nx = F.Nx, Nz = F.Nz;
    for (int l = 0; l < Nz; ++l)
        for (int k = 0; k < Nx; ++k) {
            double r = 1.0 + 0.5 * std::sin(3.0 * k / Nx) * std::cos(3.0 * l / Nz);
            F.D33_s[(size_t)k + (size_t)l * Nx] = r;
        }
    RunGoldenCheck(F, "anisotropic_D33");
}

// Heterogeneous D11/D22 — cell-centred normal viscosity varies.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledHeterogeneousD11D22) {
    MeshFixture F;
    F.build(21, 21);
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            double r = 1.0 + 0.8 * std::cos(2.0 * k / Ncx) * std::sin(2.0 * l / Ncz);
            F.D11_n[(size_t)k + (size_t)l * Ncx] = 2.0 * r;
            F.D22_n[(size_t)k + (size_t)l * Ncx] = 2.0 * r;
        }
    RunGoldenCheck(F, "heterogeneous_D11_D22");
}

// Compressible branch — comp_cells = 1, bet_n = 1, dt finite. Exercises the
// cont_row_matvec `pc = beta / dt` diagonal term and D11/D22 OOP*comp terms.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledCompressible) {
    MeshFixture F;
    F.build(21, 21);
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    F.model.compressible = 1;
    F.model.dt = 0.1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            F.comp_cells[(size_t)k + (size_t)l * Ncx] = 1;
            F.bet_n     [(size_t)k + (size_t)l * Ncx] = 0.3;
        }
    RunGoldenCheck(F, "compressible");
}

// Free-surface stabilisation — model.free_surface_stab > 0 activates the
// stab correction on the diagonal uC/vC coefficients.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledFreeSurfaceStab) {
    MeshFixture F;
    F.build(21, 21);
    F.model.free_surface_stab = 1.0;
    F.model.dt = 1.0;
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            F.rho_n[(size_t)k + (size_t)l * Ncx] = 1.0 + 0.5 * (double)l / Ncz;
        }
    const int Nx = F.Nx, Nzvx = F.Nz + 1, Nxvz = F.Nx + 1, Nz = F.Nz;
    for (size_t k = 0; k < (size_t)Nx * Nzvx; ++k) F.gx[k] = 0.1;
    for (size_t k = 0; k < (size_t)Nxvz * Nz; ++k) F.gz[k] = -9.81;
    RunGoldenCheck(F, "free_surface_stab");
}

// Out-of-plane (axisymmetric-ish) — OOP = 3/2 multiplier on the normal-stress
// coefficients.
TEST(StokesMatvecEquivalence, PicardMatvecMatchesAssembledOutOfPlane) {
    MeshFixture F;
    F.build(21, 21);
    F.model.out_of_plane = 1;
    RunGoldenCheck(F, "out_of_plane");
}

// ---------------------------------------------------------------------------
// Golden cell-block check (task 15.6 / D11).
//
// For every strict-interior pressure cell where the 5 block DOFs are all
// real equations, `StokesCellBlockMDOODZ` must produce:
//   out_A[r][c] * celvol == A-or-B-or-C-or-D[eqn_r][eqn_c]   (inside block)
//   out_b[r]   * celvol == (A·x_u + B·x_p)[eqn_r]
//                          - Σ_c out_A[r][c]*celvol * x_full_at_dof_c
// This binds the Vanka block builder to the same correctness contract as
// `ApplyStokesOperatorMDOODZ`: both consume the same coefficient fillers,
// both match the assembled SparseMat bit-for-bit.
// ---------------------------------------------------------------------------
namespace {

static double CsrAt(const SparseMat *M, int row, int col) {
    if (row < 0 || row >= M->neq) return 0.0;
    int rs = M->Ic[row], re = M->Ic[row + 1];
    for (int k = rs; k < re; ++k) if (M->J[k] == col) return M->A[k];
    return 0.0;
}

static double RunBlockGoldenCheck(MeshFixture &F, const std::string &label,
                                  double tol = 1e-11, int n_trials = 3) {
    SparseMat Stokes{}, StokesA{}, StokesB{}, StokesC{}, StokesD{};
    EvalNumberOfEquations(&F.mesh, &Stokes);
    StokesA.eqn_u = StokesB.eqn_u = StokesC.eqn_u = StokesD.eqn_u = Stokes.eqn_u;
    StokesA.eqn_v = StokesB.eqn_v = StokesC.eqn_v = StokesD.eqn_v = Stokes.eqn_v;
    StokesA.eqn_p = StokesB.eqn_p = StokesC.eqn_p = StokesD.eqn_p = Stokes.eqn_p;
    std::vector<double> bbcA((size_t)Stokes.neq_mom, 0.0), bbcB((size_t)Stokes.neq_mom, 0.0);
    std::vector<double> bbcC((size_t)Stokes.neq_cont, 0.0), bbcD((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bA  ((size_t)Stokes.neq_mom, 0.0), bB  ((size_t)Stokes.neq_mom, 0.0);
    std::vector<double> bC  ((size_t)Stokes.neq_cont, 0.0), bD  ((size_t)Stokes.neq_cont, 0.0);
    StokesA.bbc = bbcA.data(); StokesB.bbc = bbcB.data();
    StokesC.bbc = bbcC.data(); StokesD.bbc = bbcD.data();
    StokesA.b = bA.data(); StokesB.b = bB.data();
    StokesC.b = bC.data(); StokesD.b = bD.data();
    BuildStokesOperatorDecoupled(&F.mesh, F.model, 0,
                                 F.p_corr.data(), F.p_in.data(),
                                 F.u_in.data(),   F.v_in.data(),
                                 &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1);

    const int Nx = F.Nx, Nz = F.Nz;
    const int Ncx = Nx - 1, Ncz = Nz - 1, Nxvz = Nx + 1;
    const double celvol = F.dx * F.dz;

    std::mt19937 rng(12345u + std::hash<std::string>{}(label));
    std::uniform_real_distribution<double> U(-1.0, 1.0);

    double worst_err = 0.0;
    int total_checks = 0;

    for (int trial = 0; trial < n_trials; ++trial) {
        std::vector<double> u_full(F.u_in.size(), 0.0);
        std::vector<double> v_full(F.v_in.size(), 0.0);
        std::vector<double> p_full(F.p_in.size(), 0.0);
        for (size_t k = 0; k < u_full.size(); ++k) u_full[k] = U(rng);
        for (size_t k = 0; k < v_full.size(); ++k) v_full[k] = U(rng);
        for (size_t k = 0; k < p_full.size(); ++k) p_full[k] = U(rng);
        for (size_t k = 0; k < u_full.size(); ++k)
            if (F.BCu_type[k] != -1 && F.BCu_type[k] != 2 && F.BCu_type[k] != -2) u_full[k] = 0.0;
        for (size_t k = 0; k < v_full.size(); ++k)
            if (F.BCv_type[k] != -1 && F.BCv_type[k] != 2) v_full[k] = 0.0;
        for (size_t k = 0; k < p_full.size(); ++k)
            if (F.BCp_type[k] != -1) p_full[k] = 0.0;

        // Compute the assembled matvec in compressed space so we can
        // check `out_b` via the aggregate identity below.
        std::vector<double> x;
        PackFullToCompressed(&F.mesh, &Stokes, u_full.data(), v_full.data(), p_full.data(), x);
        std::vector<double> x_u(Stokes.neq_mom), x_p(Stokes.neq_cont);
        for (int i = 0; i < Stokes.neq_mom;  ++i) x_u[i] = x[i];
        for (int i = 0; i < Stokes.neq_cont; ++i) x_p[i] = x[Stokes.neq_mom + i];
        std::vector<double> Ax_u, Bp_u, Cx_p, Dp_p;
        CsrMatvec(&StokesA, x_u, Ax_u);
        CsrMatvec(&StokesB, x_p, Bp_u);
        CsrMatvec(&StokesC, x_u, Cx_p);
        CsrMatvec(&StokesD, x_p, Dp_p);

        for (int j = 1; j < Ncz - 1; ++j) {
            for (int i = 1; i < Ncx - 1; ++i) {
                const int iVxL = i     + (j + 1) * Nx;
                const int iVxR = (i+1) + (j + 1) * Nx;
                const int iVzB = (i+1) + (j    ) * Nxvz;
                const int iVzT = (i+1) + (j + 1) * Nxvz;
                const int iPC  = i     + (j    ) * Ncx;
                if (F.BCu_type[iVxL] != -1) continue;
                if (F.BCu_type[iVxR] != -1) continue;
                if (F.BCv_type[iVzB] != -1) continue;
                if (F.BCv_type[iVzT] != -1) continue;
                if (F.BCp_type[iPC]  != -1) continue;

                double block[5][5]; double out_b[5];
                int rc = StokesCellBlockMDOODZ(&F.mesh, F.model,
                                               u_full.data(), v_full.data(), p_full.data(),
                                               i, j, block, out_b);
                if (rc != 0) {
                    ADD_FAILURE() << label << " StokesCellBlockMDOODZ rc=" << rc
                                  << " cell=(" << i << "," << j << ")";
                    continue;
                }

                const int eqn_vel[4] = {
                    Stokes.eqn_u[iVxL], Stokes.eqn_u[iVxR],
                    Stokes.eqn_v[iVzB], Stokes.eqn_v[iVzT],
                };
                const int eqn_p  = Stokes.eqn_p[iPC];
                const int pc_row = eqn_p - Stokes.neq_mom;

                // --- Inside-block entry checks. ---
                // Velocity rows 0..3.
                for (int r = 0; r < 4; ++r) {
                    for (int c = 0; c < 4; ++c) {
                        double ref = CsrAt(&StokesA, eqn_vel[r], eqn_vel[c]) / celvol;
                        double err = std::fabs(block[r][c] - ref);
                        worst_err = std::max(worst_err, err);
                        EXPECT_LT(err, tol) << label
                            << " cell=(" << i << "," << j << ") r=" << r << " c=" << c
                            << " block=" << block[r][c] << " ref=" << ref;
                    }
                    // Pressure column (stored in B).
                    double ref_p = CsrAt(&StokesB, eqn_vel[r], pc_row) / celvol;
                    double err_p = std::fabs(block[r][4] - ref_p);
                    worst_err = std::max(worst_err, err_p);
                    EXPECT_LT(err_p, tol) << label
                        << " cell=(" << i << "," << j << ") r=" << r << " c=4"
                        << " block=" << block[r][4] << " ref=" << ref_p;
                }
                // Continuity row (row 4). C holds velocity coupling, D pressure self.
                for (int c = 0; c < 4; ++c) {
                    double ref = CsrAt(&StokesC, pc_row, eqn_vel[c]) / celvol;
                    double err = std::fabs(block[4][c] - ref);
                    worst_err = std::max(worst_err, err);
                    EXPECT_LT(err, tol) << label
                        << " cell=(" << i << "," << j << ") r=4 c=" << c
                        << " block=" << block[4][c] << " ref=" << ref;
                }
                double ref_pp = CsrAt(&StokesD, pc_row, pc_row) / celvol;
                double err_pp = std::fabs(block[4][4] - ref_pp);
                worst_err = std::max(worst_err, err_pp);
                EXPECT_LT(err_pp, tol) << label
                    << " cell=(" << i << "," << j << ") r=4 c=4"
                    << " block=" << block[4][4] << " ref=" << ref_pp;

                // --- Aggregate out_b identity: for each row r the block
                // builder promises out_b[r] = (assembled_matvec[eqn_r] /
                // celvol) - Σ_c out_A[r][c]*x_full_at_dof_c. ---
                const double x_dof[5] = {
                    u_full[iVxL], u_full[iVxR], v_full[iVzB], v_full[iVzT], p_full[iPC]
                };
                for (int r = 0; r < 4; ++r) {
                    double inside = 0.0;
                    for (int c = 0; c < 5; ++c) inside += block[r][c] * x_dof[c];
                    double total  = (Ax_u[eqn_vel[r]] + Bp_u[eqn_vel[r]]) / celvol;
                    double ref    = total - inside;
                    double err    = std::fabs(out_b[r] - ref);
                    worst_err = std::max(worst_err, err);
                    EXPECT_LT(err, tol) << label
                        << " cell=(" << i << "," << j << ") r=" << r
                        << " out_b=" << out_b[r] << " ref=" << ref;
                }
                {
                    double inside = 0.0;
                    for (int c = 0; c < 5; ++c) inside += block[4][c] * x_dof[c];
                    double total  = (Cx_p[pc_row] + Dp_p[pc_row]) / celvol;
                    double ref    = total - inside;
                    double err    = std::fabs(out_b[4] - ref);
                    worst_err = std::max(worst_err, err);
                    EXPECT_LT(err, tol) << label
                        << " cell=(" << i << "," << j << ") r=4"
                        << " out_b=" << out_b[4] << " ref=" << ref;
                }

                ++total_checks;
            }
        }
    }

    EXPECT_GT(total_checks, 0) << label << ": no strict-interior cells checked";

    FreeMat(&StokesA); FreeMat(&StokesB);
    FreeMat(&StokesC); FreeMat(&StokesD);
    DoodzFree(Stokes.eqn_u);
    DoodzFree(Stokes.eqn_v);
    DoodzFree(Stokes.eqn_p);
    return worst_err;
}

}  // namespace

TEST(StokesMatvecEquivalence, StokesCellBlockMatchesAssembledConstantViscosity) {
    MeshFixture F;
    F.build(21, 21);
    RunBlockGoldenCheck(F, "block_constant_visc");
}

TEST(StokesMatvecEquivalence, StokesCellBlockMatchesAssembledHeterogeneousD11D22) {
    MeshFixture F;
    F.build(21, 21);
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            double r = 1.0 + 0.8 * std::cos(2.0 * k / Ncx) * std::sin(2.0 * l / Ncz);
            F.D11_n[(size_t)k + (size_t)l * Ncx] = 2.0 * r;
            F.D22_n[(size_t)k + (size_t)l * Ncx] = 2.0 * r;
        }
    RunBlockGoldenCheck(F, "block_heterogeneous_D11_D22");
}

TEST(StokesMatvecEquivalence, StokesCellBlockMatchesAssembledAnisotropicD33) {
    MeshFixture F;
    F.build(21, 21);
    const int Nx = F.Nx, Nz = F.Nz;
    for (int l = 0; l < Nz; ++l)
        for (int k = 0; k < Nx; ++k) {
            double r = 1.0 + 0.5 * std::sin(3.0 * k / Nx) * std::cos(3.0 * l / Nz);
            F.D33_s[(size_t)k + (size_t)l * Nx] = r;
        }
    RunBlockGoldenCheck(F, "block_anisotropic_D33");
}

TEST(StokesMatvecEquivalence, StokesCellBlockMatchesAssembledCompressible) {
    MeshFixture F;
    F.build(21, 21);
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    F.model.compressible = 1;
    F.model.dt = 0.1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            F.comp_cells[(size_t)k + (size_t)l * Ncx] = 1;
            F.bet_n     [(size_t)k + (size_t)l * Ncx] = 0.3;
        }
    RunBlockGoldenCheck(F, "block_compressible");
}

// ---------------------------------------------------------------------------
// Newton-mode golden check (tasks 15.19 / 15.20).
//
// The Newton Jacobian matvec (`ApplyStokesOperatorMDOODZ` with
// `model.Newton == 1`) replaces the Picard stencils with:
//   - x-momentum: 9-Vx / 12-Vz / 6-P row from `Xjacobian_InnerNodesDecoupled3`
//   - z-momentum: 9-Vz / 12-Vx / 6-P row from `Zjacobian_InnerNodesDecoupled3`
// The reference operator is `BuildJacobianOperatorDecoupled` (not
// `BuildStokesOperatorDecoupled`).
//
// The assertion loop covers Vx AND Vz rows; continuity (P) rows are
// identical to the Picard path (continuity does not depend on D-tensor
// or `Newton`), so they are exercised by the existing
// `StokesMatvecMatchesAssembled*` suite and are not re-checked here.
// ---------------------------------------------------------------------------
namespace {

static double RunNewtonMomGoldenCheck(MeshFixture &F, const std::string &label,
                                       double tol = 1e-12, int n_trials = 10) {
    F.model.Newton = 1;
    // BuildJacobianOperatorDecoupled forces comp=1 internally regardless
    // of model.compressible; the ApplyStokesOperatorMDOODZ Newton branch
    // mirrors that. We leave model.compressible on the caller to set if
    // they additionally want mesh->comp_cells[] to drive anything (the
    // xmom Jacobian formulas don't read comp_cells).

    SparseMat Stokes{}, StokesA{}, StokesB{}, StokesC{}, StokesD{};
    EvalNumberOfEquations(&F.mesh, &Stokes);
    StokesA.eqn_u = StokesB.eqn_u = StokesC.eqn_u = StokesD.eqn_u = Stokes.eqn_u;
    StokesA.eqn_v = StokesB.eqn_v = StokesC.eqn_v = StokesD.eqn_v = Stokes.eqn_v;
    StokesA.eqn_p = StokesB.eqn_p = StokesC.eqn_p = StokesD.eqn_p = Stokes.eqn_p;
    std::vector<double> bbcA((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bbcB((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bbcC((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bbcD((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bA  ((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bB  ((size_t)Stokes.neq_mom,  0.0);
    std::vector<double> bC  ((size_t)Stokes.neq_cont, 0.0);
    std::vector<double> bD  ((size_t)Stokes.neq_cont, 0.0);
    StokesA.bbc = bbcA.data(); StokesB.bbc = bbcB.data();
    StokesC.bbc = bbcC.data(); StokesD.bbc = bbcD.data();
    StokesA.b = bA.data(); StokesB.b = bB.data();
    StokesC.b = bC.data(); StokesD.b = bD.data();

    BuildJacobianOperatorDecoupled(&F.mesh, F.model, 0,
                                   F.p_corr.data(), F.p_in.data(),
                                   F.u_in.data(),   F.v_in.data(),
                                   &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1);

    std::mt19937 rng(99u + std::hash<std::string>{}(label));
    std::uniform_real_distribution<double> U(-1.0, 1.0);

    double worst_rel = 0.0;
    for (int trial = 0; trial < n_trials; ++trial) {
        std::vector<double> u_full(F.u_in.size(), 0.0);
        std::vector<double> v_full(F.v_in.size(), 0.0);
        std::vector<double> p_full(F.p_in.size(), 0.0);
        for (size_t k = 0; k < u_full.size(); ++k) u_full[k] = U(rng);
        for (size_t k = 0; k < v_full.size(); ++k) v_full[k] = U(rng);
        for (size_t k = 0; k < p_full.size(); ++k) p_full[k] = U(rng);
        for (size_t k = 0; k < u_full.size(); ++k)
            if (F.BCu_type[k] != -1 && F.BCu_type[k] != 2 && F.BCu_type[k] != -2)
                u_full[k] = 0.0;
        for (size_t k = 0; k < v_full.size(); ++k)
            if (F.BCv_type[k] != -1 && F.BCv_type[k] != 2)
                v_full[k] = 0.0;
        for (size_t k = 0; k < p_full.size(); ++k)
            if (F.BCp_type[k] != -1) p_full[k] = 0.0;

        std::vector<double> x;
        PackFullToCompressed(&F.mesh, &Stokes, u_full.data(), v_full.data(),
                             p_full.data(), x);
        std::vector<double> x_u(Stokes.neq_mom), x_p(Stokes.neq_cont);
        for (int i = 0; i < Stokes.neq_mom;  ++i) x_u[i] = x[i];
        for (int i = 0; i < Stokes.neq_cont; ++i) x_p[i] = x[Stokes.neq_mom + i];

        std::vector<double> Ax_u, Bp_u;
        CsrMatvec(&StokesA, x_u, Ax_u);
        CsrMatvec(&StokesB, x_p, Bp_u);

        std::vector<double> ru(F.u_in.size(), 0.0);
        std::vector<double> rv(F.v_in.size(), 0.0);
        std::vector<double> rp(F.p_in.size(), 0.0);
        ApplyStokesOperatorMDOODZ(&F.mesh, F.model,
                                  u_full.data(), v_full.data(), p_full.data(),
                                  ru.data(),     rv.data(),     rp.data());

        double max_abs = 0.0, max_err = 0.0;
        // Vx rows.
        for (int l = 0; l < F.Nz + 1; ++l)
            for (int k = 0; k < F.Nx; ++k) {
                int kk  = k + l * F.Nx;
                int eqn = Stokes.eqn_u[kk];
                if (eqn < 0) continue;
                double expected = Ax_u[eqn] + Bp_u[eqn];
                double got      = ru[kk];
                double err      = std::fabs(expected - got);
                max_abs = std::max(max_abs, std::fabs(expected));
                max_err = std::max(max_err, err);
            }
        // Vz rows.
        for (int l = 0; l < F.Nz; ++l)
            for (int k = 0; k < F.Nx + 1; ++k) {
                int kk  = k + l * (F.Nx + 1);
                int eqn = Stokes.eqn_v[kk];
                if (eqn < 0) continue;
                double expected = Ax_u[eqn] + Bp_u[eqn];
                double got      = rv[kk];
                double err      = std::fabs(expected - got);
                max_abs = std::max(max_abs, std::fabs(expected));
                max_err = std::max(max_err, err);
            }

        const double rel = (max_abs > 0.0) ? max_err / max_abs : max_err;
        worst_rel = std::max(worst_rel, rel);
        EXPECT_LT(rel, tol) << label << " trial " << trial
                            << " max_err=" << max_err
                            << " max_abs=" << max_abs;
    }

    FreeMat(&StokesA); FreeMat(&StokesB);
    FreeMat(&StokesC); FreeMat(&StokesD);
    DoodzFree(Stokes.eqn_u);
    DoodzFree(Stokes.eqn_v);
    DoodzFree(Stokes.eqn_p);
    return worst_rel;
}

}  // namespace

// Baseline Newton test: closed-box Dirichlet, unit-isotropic D-tensor.
// Exercises only the Picard-equivalent part of Xjacobian (no anisotropy,
// no compressibility cross-term).
TEST(StokesMatvecEquivalence, JacobianMatvecMatchesAssembledConstantViscosity) {
    MeshFixture F;
    F.build(21, 21);
    RunNewtonMomGoldenCheck(F, "newton_mom_constant_visc");
}

// Heterogeneous D11 + D12 normal-stress coupling. D12 appears only in
// Newton's Jacobian (the Picard path ignores it), so this directly
// exercises the Newton-specific coefficient cross-terms in
// xmom_fill_coeffs_newton.
TEST(StokesMatvecEquivalence, JacobianMatvecMatchesAssembledNormalCoupling) {
    MeshFixture F;
    F.build(21, 21);
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            double s = 1.0 + 0.4 * std::sin(2.5 * k / Ncx) * std::cos(2.5 * l / Ncz);
            F.D11_n[(size_t)k + (size_t)l * Ncx] = 2.0 * s;
            F.D22_n[(size_t)k + (size_t)l * Ncx] = 2.0 * s;
            F.D12_n[(size_t)k + (size_t)l * Ncx] = 0.3 * s;
            F.D21_n[(size_t)k + (size_t)l * Ncx] = 0.3 * s;
        }
    RunNewtonMomGoldenCheck(F, "newton_mom_normal_coupling");
}

// Full anisotropic Jacobian: all D3* shear-stress vertex entries
// non-trivial. D13/D14 appear as rotation couplings that Picard drops.
// ---------------------------------------------------------------------------
// D7 symmetric-smoothing invariant (task 15.21).
//
// StokesCellBlockMDOODZ must return the symmetric Picard part in every
// mode. The test pins this property directly: the 5×5 block produced
// with `model.Newton == 1` (and any non-trivial anisotropic D-tensor
// loaded to exercise the Jacobian coefficients) must bit-identically
// match the 5×5 block produced with `model.Newton == 0 &
// model.compressible == 1` at the same physical state. Any leakage of
// Newton-Jacobian cross-terms into the smoother block would fail
// this check loudly.
// ---------------------------------------------------------------------------
TEST(StokesMatvecEquivalence, CellBlockNewtonReducesToSymmetricPicardPart) {
    MeshFixture F;
    F.build(21, 21);
    // Load a non-trivial anisotropic D-tensor — if the Newton-mode
    // block accidentally called the Jacobian branch, D13/D14/D31-D34
    // contributions would leak into out_A and the invariant would fail.
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1, Nx = F.Nx, Nz = F.Nz;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            F.D11_n[(size_t)k + (size_t)l * Ncx] = 2.0;
            F.D22_n[(size_t)k + (size_t)l * Ncx] = 2.0;
            F.D12_n[(size_t)k + (size_t)l * Ncx] = 0.3;
            F.D13_n[(size_t)k + (size_t)l * Ncx] = 0.2;
            F.D14_n[(size_t)k + (size_t)l * Ncx] = 0.1;
            F.D21_n[(size_t)k + (size_t)l * Ncx] = 0.3;
        }
    for (int l = 0; l < Nz; ++l)
        for (int k = 0; k < Nx; ++k) {
            F.D33_s[(size_t)k + (size_t)l * Nx] = 1.0;
            F.D31_s[(size_t)k + (size_t)l * Nx] = 0.2;
            F.D32_s[(size_t)k + (size_t)l * Nx] = 0.15;
            F.D34_s[(size_t)k + (size_t)l * Nx] = 0.1;
        }

    // Populate an arbitrary admissible state.
    std::mt19937 rng(777u);
    std::uniform_real_distribution<double> U(-1.0, 1.0);
    std::vector<double> u_full(F.u_in.size(), 0.0);
    std::vector<double> v_full(F.v_in.size(), 0.0);
    std::vector<double> p_full(F.p_in.size(), 0.0);
    for (auto &x : u_full) x = U(rng);
    for (auto &x : v_full) x = U(rng);
    for (auto &x : p_full) x = U(rng);

    params model_newton = F.model;
    model_newton.Newton       = 1;
    model_newton.compressible = 0;   // block forces comp=1 when Newton
    params model_picard = F.model;
    model_picard.Newton       = 0;
    model_picard.compressible = 1;   // match the block's forced comp

    double max_err_A = 0.0, max_err_b = 0.0;
    for (int j = 0; j < Ncz; ++j)
        for (int i = 0; i < Ncx; ++i) {
            double A_n[5][5], b_n[5];
            double A_p[5][5], b_p[5];
            int rn = StokesCellBlockMDOODZ(&F.mesh, model_newton,
                                           u_full.data(), v_full.data(), p_full.data(),
                                           i, j, A_n, b_n);
            int rp = StokesCellBlockMDOODZ(&F.mesh, model_picard,
                                           u_full.data(), v_full.data(), p_full.data(),
                                           i, j, A_p, b_p);
            ASSERT_EQ(rn, 0) << "i=" << i << " j=" << j;
            ASSERT_EQ(rp, 0) << "i=" << i << " j=" << j;
            for (int r = 0; r < 5; ++r) {
                for (int c = 0; c < 5; ++c)
                    max_err_A = std::max(max_err_A, std::fabs(A_n[r][c] - A_p[r][c]));
                max_err_b = std::max(max_err_b, std::fabs(b_n[r] - b_p[r]));
            }
        }
    EXPECT_LT(max_err_A, 1e-14) << "max_err_A=" << max_err_A;
    EXPECT_LT(max_err_b, 1e-14) << "max_err_b=" << max_err_b;
}

TEST(StokesMatvecEquivalence, JacobianMatvecMatchesAssembledFullAnisotropic) {
    MeshFixture F;
    F.build(21, 21);
    const int Ncx = F.Nx - 1, Ncz = F.Nz - 1;
    for (int l = 0; l < Ncz; ++l)
        for (int k = 0; k < Ncx; ++k) {
            double s = 1.0 + 0.3 * std::cos(3.0 * k / Ncx);
            F.D11_n[(size_t)k + (size_t)l * Ncx] = 2.0 * s;
            F.D22_n[(size_t)k + (size_t)l * Ncx] = 2.0 * s;
            F.D12_n[(size_t)k + (size_t)l * Ncx] = 0.2 * s;
            F.D13_n[(size_t)k + (size_t)l * Ncx] = 0.15 * s;
            F.D14_n[(size_t)k + (size_t)l * Ncx] = 0.1 * s;
        }
    const int Nx = F.Nx, Nz = F.Nz;
    for (int l = 0; l < Nz; ++l)
        for (int k = 0; k < Nx; ++k) {
            double s = 1.0 + 0.25 * std::sin(2.0 * k / Nx) * std::cos(2.0 * l / Nz);
            F.D33_s[(size_t)k + (size_t)l * Nx] = s;
            F.D31_s[(size_t)k + (size_t)l * Nx] = 0.2 * s;
            F.D32_s[(size_t)k + (size_t)l * Nx] = 0.15 * s;
            F.D34_s[(size_t)k + (size_t)l * Nx] = 0.1 * s;
        }
    RunNewtonMomGoldenCheck(F, "newton_mom_full_anisotropic");
}
