// add-gmg-stokes-defence §4.2 — SolKz (high-contrast Stokes) GMG↔CHOLMOD
// equivalence.
//
// Zhong & Gurnis (1996) "SolKz" benchmark:
//   * Prescribed viscosity  eta(z) = exp(2·B·z)           (B = user0 = 8)
//   * Prescribed density    rho(x,z) = -sin(pi z) cos(K pi x)  (K = user1)
//   * Box  [0,1] x [0,1]    free-slip on all four walls, gravity (0, -1).
//
// MDOODZ has no marker-level density/viscosity callback — the mesh body
// force and fine-level operator come from `materials->rho[p]` and
// `materials->eta0[p]` (phase tables). We discretise the continuous
// eta(x,z) × rho(x,z) fields on a N_BANDS_X × N_BANDS_Z grid of
// constant-viscosity phases (one phase per x_band × z_band). `SetPhase`
// assigns each marker the phase index `ix + N_BANDS_X * iz`;
// `MutateInput` overrides the placeholder eta0, rho declared in the .txt
// with the band-midpoint SolKz values.
//
// Scope of this test
//   §4.2 is a DUAL-SOLVER equivalence test: GMG must match CHOLMOD to
//   `1e-8` on the SAME discrete Stokes system. It is NOT an analytic
//   SolKz fidelity test — comparison against the analytic solution
//   (for defence figures + convergence study) is §5.
//
// Pressure bound
//   GMG projects the Stokes null space every V-cycle, CHOLMOD pins P
//   through the penalty coupling. The two pressures differ by a
//   constant after a single Picard iteration; we mean-subtract before
//   comparing, and hold the pressure to `1e-6` rather than `1e-8`
//   (same reasoning as GmgStokesEquivalence header).
//
// Solver choice — why `lin_solver = -1` (not 0) on the CHOLMOD twin
//   InputOutput.c auto-remaps lin_solver = 0 → 2 (KillerSolver, a
//   Powell–Hestenes / KSP hybrid) on problems that need a saddle-point
//   solve. On 10^7-contrast SolKz, the KSP inner stagnates the PH
//   outer loop at ~1e-6 relative momentum residual, missing the 1e-8
//   dual-solver bound. `lin_solver = -1` routes to
//   `DirectStokesDecoupledComp` (PH outer + direct CHOLMOD inner);
//   the direct inner factorises to machine precision and PH then
//   converges in ≤ 3 iterations. This is the right reference for a
//   solver-equivalence test: the CHOLMOD twin should be as close to
//   "the exact discrete solution" as possible.

extern "C" {
#include "mdoodz.h"
}

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <vector>

#include <gtest/gtest.h>

#include "TestHelpers.h"

class SolKzGmgEquivalence : public ::testing::Test {
protected:
  static constexpr int N_BANDS_X = 4;
  static constexpr int N_BANDS_Z = 4;
  static constexpr int N_PHASES  = N_BANDS_X * N_BANDS_Z;  // = 16, ≤ 20

  MdoodzSetup setup;

  void SetUp() override {
    setup = {
        .SetParticles = new SetParticles_ff{
            .SetPhase   = SetPhase,
            .SetDensity = SetDensity,
        },
        .SetBCs = new SetBCs_ff{
            .SetBCVx    = SetBCVx,
            .SetBCVz    = SetBCVz,
            .SetBCPType = SetBCPType,
        },
        .MutateInput = MutateInputCB,
    };
  }

  // Assign phase index from (x,z) band coordinates on the unit square.
  //   ix = floor((x - xmin) / Lx * N_BANDS_X)   (clamped to [0, N_BANDS_X-1])
  //   iz = floor((z - zmin) / Lz * N_BANDS_Z)   (clamped to [0, N_BANDS_Z-1])
  //   phase = ix + N_BANDS_X * iz
  static int SetPhase(MdoodzInput *input, Coordinates c) {
    const double xLen = input->model.xmax - input->model.xmin;
    const double zLen = input->model.zmax - input->model.zmin;
    int ix = (int)std::floor((c.x - input->model.xmin) / xLen *
                             (double)N_BANDS_X);
    int iz = (int)std::floor((c.z - input->model.zmin) / zLen *
                             (double)N_BANDS_Z);
    if (ix < 0)             ix = 0;
    if (ix > N_BANDS_X - 1) ix = N_BANDS_X - 1;
    if (iz < 0)             iz = 0;
    if (iz > N_BANDS_Z - 1) iz = N_BANDS_Z - 1;
    return ix + N_BANDS_X * iz;
  }

  // Kept for completeness so the mesh's per-marker rho matches the
  // phase-averaged rho at band boundaries. The mesh body force is
  // rebuilt every step from `materials->rho[phase]` via UpdateDensity
  // (RheologyDensity.c), so the return value of SetDensity is
  // effectively unused once the phase assignment is correct — we
  // just echo the phase table value for consistency.
  static double SetDensity(MdoodzInput *input, Coordinates, int phase) {
    return input->materials.rho[phase];
  }

  // Free-slip walls: Vx=0 on W/E, Vz=0 on S/N, tangential component free.
  // Same pattern as BlankenBenchTests.cpp and TopoBenchTests.cpp.
  static SetBC SetBCVx(MdoodzInput *, POSITION position, Coordinates) {
    SetBC bc;
    bc.value = 0.0;
    if      (position == S || position == SW || position == SE) bc.type = 11;
    else if (position == N || position == NW || position == NE) bc.type = 13;
    else if (position == W || position == E)                    bc.type = 0;
    else                                                         bc.type = -1;
    return bc;
  }

  static SetBC SetBCVz(MdoodzInput *, POSITION position, Coordinates) {
    SetBC bc;
    bc.value = 0.0;
    if (position == W || position == E || position == SW || position == SE ||
        position == NW || position == NE) bc.type = 13;
    else if (position == S || position == N) bc.type = 0;
    else                                      bc.type = -1;
    return bc;
  }

  // Pin pressure at NE / NW corner cells to anchor the constant-P gauge.
  static char SetBCPType(MdoodzInput *, POSITION position) {
    return (position == NE || position == NW) ? 0 : -1;
  }

  // MutateInput runs after the .txt has been parsed and materials table
  // has been scaled to code units. We overwrite eta0[p] and rho[p] for
  // each of the N_PHASES phases using the band-midpoint SolKz fields.
  // In the SolKz fixture scaling.{eta,rho} = 1 so the scaling divisions
  // are a no-op, but we include them for self-consistency.
  static void MutateInputCB(MdoodzInput *input) {
    if (input->model.Nb_phases < N_PHASES) {
      fprintf(stderr,
              "SolKzGmgEquivalence: Nb_phases=%d < N_PHASES=%d — the .txt "
              "must declare at least %d phases\n",
              input->model.Nb_phases, N_PHASES, N_PHASES);
      std::abort();
    }
    const double B = input->model.user0;
    const double K = input->model.user1;
    for (int iz = 0; iz < N_BANDS_Z; ++iz) {
      for (int ix = 0; ix < N_BANDS_X; ++ix) {
        const int    p     = ix + N_BANDS_X * iz;
        const double x_mid = ((double)ix + 0.5) / (double)N_BANDS_X;
        const double z_mid = ((double)iz + 0.5) / (double)N_BANDS_Z;
        const double eta   = exp(2.0 * B * z_mid);
        const double rho   = -sin(M_PI * z_mid) * cos(K * M_PI * x_mid);
        input->materials.eta0[p] = eta / input->scaling.eta;
        input->materials.rho [p] = rho / input->scaling.rho;
      }
    }
  }
};

TEST_F(SolKzGmgEquivalence, HighContrastSolKzMatchesCholmodWithin1e8) {
  char gmgInput [] = "AnalyticalBench/SolKz_gmg.txt";
  char cholInput[] = "AnalyticalBench/SolKz_chol.txt";
  RunMDOODZ(gmgInput , &setup);
  RunMDOODZ(cholInput, &setup);

  const char *gmgFile  = "SolKz_gmg/Output00001.gzip.h5";
  const char *cholFile = "SolKz_chol/Output00001.gzip.h5";

  auto vx_gmg  = readFieldAsArray(gmgFile , "VxNodes", "Vx");
  auto vx_chol = readFieldAsArray(cholFile, "VxNodes", "Vx");
  auto vz_gmg  = readFieldAsArray(gmgFile , "VzNodes", "Vz");
  auto vz_chol = readFieldAsArray(cholFile, "VzNodes", "Vz");
  auto p_gmg   = readFieldAsArray(gmgFile , "Centers", "P");
  auto p_chol  = readFieldAsArray(cholFile, "Centers", "P");

  ASSERT_EQ(vx_gmg.size(), vx_chol.size());
  ASSERT_EQ(vz_gmg.size(), vz_chol.size());
  ASSERT_EQ(p_gmg .size(), p_chol .size());

  for (double v : vx_gmg) ASSERT_TRUE(std::isfinite(v)) << "Vx_gmg NaN/Inf";
  for (double v : vz_gmg) ASSERT_TRUE(std::isfinite(v)) << "Vz_gmg NaN/Inf";
  for (double v : p_gmg ) ASSERT_TRUE(std::isfinite(v)) << "P_gmg  NaN/Inf";

  // Sanity: the solve must have produced non-trivial velocities (body
  // force is non-zero, rho contrast is O(1), so |V| should be O(eta^-1).
  // If both solvers return identically-zero velocities the equivalence
  // L2 is artificially 0 — guard against that with a CHOLMOD magnitude
  // check.
  double max_abs_vx = 0.0;
  for (double v : vx_chol) max_abs_vx = std::max(max_abs_vx, std::abs(v));
  ASSERT_GT(max_abs_vx, 1e-6)
      << "CHOLMOD twin produced trivially-zero velocities — fixture "
         "configuration failure, not a solver equivalence result";

  // Mean-subtract the pressures for gauge-independence.
  double m_g = 0.0, m_c = 0.0;
  for (size_t k = 0; k < p_gmg.size(); ++k) { m_g += p_gmg[k]; m_c += p_chol[k]; }
  m_g /= (double)p_gmg.size();
  m_c /= (double)p_chol.size();
  std::vector<double> p_g_shift(p_gmg.size()), p_c_shift(p_chol.size());
  for (size_t k = 0; k < p_gmg.size(); ++k) {
    p_g_shift[k] = p_gmg [k] - m_g;
    p_c_shift[k] = p_chol[k] - m_c;
  }

  const double dVx = computeL2Error(vx_gmg, vx_chol);
  const double dVz = computeL2Error(vz_gmg, vz_chol);
  const double dP  = computeL2Error(p_g_shift, p_c_shift);

  printf("SolKz 51x51 dual-solver consistency: "
         "|Vx_gmg-Vx_chol|/|Vx_chol| = %.3e, "
         "|Vz_gmg-Vz_chol|/|Vz_chol| = %.3e, "
         "|P_gmg -P_chol |/|P_chol | = %.3e (mean-subtracted)\n",
         dVx, dVz, dP);

  EXPECT_LT(dVx, 1e-8) << "§4.2 SolKz: Vx GMG↔CHOLMOD disagreement";
  EXPECT_LT(dVz, 1e-8) << "§4.2 SolKz: Vz GMG↔CHOLMOD disagreement";
  EXPECT_LT(dP,  1e-6) << "§4.2 SolKz: P (mean-shifted) GMG↔CHOLMOD "
                         "disagreement (gauge residual; see header)";
}
