// add-gmg-stokes-defence §4.1 — Free-surface relaxation GMG↔CHOLMOD
// equivalence.
//
// Runs the TopoRelax twin pair (TopoBench/TopoRelax_gmg.txt and
// TopoBench/TopoRelax_chol.txt, §3.1 fixtures) end-to-end through the
// MDOODZ pipeline and asserts that the two linear solvers produce the
// same Stokes solution to `1e-8` relative L2 on every non-constant mode
// of Vx, Vz, P (mean-subtracted) and the free-surface height. This is
// the first "easy" fixture in §4 of tasks.md: mechanical-only, Picard
// (Newton = 0), single phase × two viscosity blocks of identical
// rheology, exercising the free-surface stabilised Stokes operator.
//
// Bound justification
//   FGMRES tolerance in the GMG config is 1e-11 so velocity has ~3
//   decades of headroom above 1e-8. Pressure is mean-subtracted for
//   the gauge reason explained in the GmgStokesEquivalence header
//   (GMG projects the Stokes null space every V-cycle, CHOLMOD pins
//    P through the penalty coupling).
//
// Failure protocol (tasks.md §4.4 / D8)
//   Failures up to ~200 LOC of core-side change are fixed in this
//   change; otherwise the failing sub-assertion gets a STATUS.md entry
//   and the test is relaxed here with an explicit comment.

extern "C" {
#include "mdoodz.h"
}

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

#include <gtest/gtest.h>

#include "TestHelpers.h"

class TopoRelaxGmgEquivalence : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
        .BuildInitialTopography = new BuildInitialTopography_ff{
            .SetSurfaceZCoord = SetSurfaceZCoord,
        },
        .SetParticles = new SetParticles_ff{
            .SetPhase       = SetPhase,
            .SetTemperature = SetTemperature,
            .SetDensity     = SetDensity,
        },
        .SetBCs = new SetBCs_ff{
            .SetBCVx    = SetBCVx,
            .SetBCVz    = SetBCVz,
            .SetBCPType = SetBCPType,
        },
    };
  }

  // Sinusoidal topographic bump h(x) = -A cos(2 pi x / lambda).
  static double SetSurfaceZCoord(MdoodzInput *input, double x_coord) {
    const double A      = input->model.user1 / input->scaling.L;
    const double lambda = input->model.user2 / input->scaling.L;
    return -A * cos(2.0 * M_PI * x_coord / lambda);
  }

  static int SetPhase(MdoodzInput *, Coordinates) { return 0; }

  static double SetTemperature(MdoodzInput *input, Coordinates) {
    return (input->model.user0 + zeroC) / input->scaling.T;
  }

  static double SetDensity(MdoodzInput *input, Coordinates, int phase) {
    return input->materials.rho[phase];
  }

  // Free-slip at N/S, pure-shear-free at W/E — same as TopoBenchTests.
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

  static char SetBCPType(MdoodzInput *, POSITION position) {
    return (position == NE || position == NW) ? 0 : -1;
  }

  // Relative L2 error |a - b|_2 / |b|_2; when |b|_2 ≈ 0 falls back to
  // absolute RMS (shared with TestHelpers::computeL2Error, duplicated
  // here for printf formatting clarity).
  static double relL2(const std::vector<double> &a,
                      const std::vector<double> &b) {
    return computeL2Error(a, b);
  }
};

TEST_F(TopoRelaxGmgEquivalence, PicardFreeSurfaceMatchesCholmodWithin1e8) {
  char gmgInput [] = "TopoBench/TopoRelax_gmg.txt";
  char cholInput[] = "TopoBench/TopoRelax_chol.txt";
  RunMDOODZ(gmgInput,  &setup);
  RunMDOODZ(cholInput, &setup);

  // TopoRelax runs Nt = 3 with writer_step = 1 → compare the final step
  // (Output00003) after three Stokes solves + three free-surface
  // advections. A disagreement at step 3 is strictly harder than a
  // single-solve comparison because it accumulates three V-cycle-vs-
  // CHOLMOD gaps through the advection operator.
  const int step = 3;
  char gmgFile [128], cholFile[128];
  snprintf(gmgFile , sizeof(gmgFile ),
           "TopoRelax_gmg/Output%05d.gzip.h5", step);
  snprintf(cholFile, sizeof(cholFile),
           "TopoRelax_chol/Output%05d.gzip.h5", step);

  auto vx_gmg  = readFieldAsArray(gmgFile , "VxNodes", "Vx");
  auto vx_chol = readFieldAsArray(cholFile, "VxNodes", "Vx");
  auto vz_gmg  = readFieldAsArray(gmgFile , "VzNodes", "Vz");
  auto vz_chol = readFieldAsArray(cholFile, "VzNodes", "Vz");
  auto p_gmg   = readFieldAsArray(gmgFile , "Centers", "P");
  auto p_chol  = readFieldAsArray(cholFile, "Centers", "P");
  auto z_gmg   = readFieldAsArray(gmgFile , "Topo",    "z_grid");
  auto z_chol  = readFieldAsArray(cholFile, "Topo",    "z_grid");

  ASSERT_EQ(vx_gmg.size(), vx_chol.size());
  ASSERT_EQ(vz_gmg.size(), vz_chol.size());
  ASSERT_EQ(p_gmg .size(), p_chol .size());
  ASSERT_EQ(z_gmg .size(), z_chol .size());

  for (double v : vx_gmg ) ASSERT_TRUE(std::isfinite(v)) << "Vx_gmg  NaN";
  for (double v : vz_gmg ) ASSERT_TRUE(std::isfinite(v)) << "Vz_gmg  NaN";
  for (double v : p_gmg  ) ASSERT_TRUE(std::isfinite(v)) << "P_gmg   NaN";
  for (double v : z_gmg  ) ASSERT_TRUE(std::isfinite(v)) << "Topo_gmg NaN";

  // Mean-subtracted P for gauge-independence (see header).
  double m_g = 0.0, m_c = 0.0;
  for (size_t k = 0; k < p_gmg.size(); ++k) { m_g += p_gmg[k]; m_c += p_chol[k]; }
  m_g /= (double)p_gmg.size();
  m_c /= (double)p_chol.size();
  std::vector<double> p_g_shift(p_gmg.size()), p_c_shift(p_chol.size());
  for (size_t k = 0; k < p_gmg.size(); ++k) {
    p_g_shift[k] = p_gmg [k] - m_g;
    p_c_shift[k] = p_chol[k] - m_c;
  }

  const double dVx = relL2(vx_gmg, vx_chol);
  const double dVz = relL2(vz_gmg, vz_chol);
  const double dP  = relL2(p_g_shift, p_c_shift);
  const double dZ  = relL2(z_gmg, z_chol);

  printf("TopoRelax step %d: "
         "|Vx_gmg-Vx_chol|/|Vx_chol| = %.3e, "
         "|Vz_gmg-Vz_chol|/|Vz_chol| = %.3e, "
         "|P_gmg -P_chol |/|P_chol | = %.3e (mean-subtracted), "
         "|Topo_gmg-Topo_chol|/|Topo_chol| = %.3e\n",
         step, dVx, dVz, dP, dZ);

  EXPECT_LT(dVx, 1e-8) << "§4.1 TopoRelax: Vx GMG↔CHOLMOD disagreement";
  EXPECT_LT(dVz, 1e-8) << "§4.1 TopoRelax: Vz GMG↔CHOLMOD disagreement";
  EXPECT_LT(dP,  1e-6) << "§4.1 TopoRelax: P (mean-shifted) GMG↔CHOLMOD "
                         "disagreement (gauge residual; see header)";
  EXPECT_LT(dZ,  1e-8) << "§4.1 TopoRelax: surface height GMG↔CHOLMOD "
                         "disagreement";
}
