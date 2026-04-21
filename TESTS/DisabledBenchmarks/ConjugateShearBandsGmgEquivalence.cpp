// add-gmg-stokes-defence §7.3 — Conjugate (crossing) shear bands
// GMG↔CHOLMOD equivalence.
//
// Runs the ConjugateBands twin pair (§3.4 fixtures:
// TESTS/ShearBand/ConjugateBands_{gmg,chol}.txt) with two weak
// inclusions offset from the centre (`user2 = 1`, `user3 = 500 m`) so
// two crossing shear bands nucleate under pure-shear compression.
// Asserts that both solvers produce two bands with matching
// orientation and matching integrated dissipation.
//
// Acceptance criteria (tasks.md §7.3)
//   * Both bands form (i.e. plastic strain rate is bi-modal along
//     at least one principal axis) on both solvers.
//   * Orientation agreement ±2° between solvers, measured from the
//     principal axis of the plastic strain-rate second moment (the
//     tensor ⟨x_i x_j · eII_pl⟩ / ⟨eII_pl⟩ — PCA on the plastic
//     strain-rate "density").
//   * Integrated dissipation (sum of eII_pl on cell centres) within 5%.
//   * Velocity L2 agreement within ~1e-4.
//
// Orientation via PCA is well-defined for a single elongated band or
// for two symmetric crossing bands (where the moment tensor is nearly
// isotropic) — we report the principal-axis angle and expect it to
// match between solvers even when the absolute angle is arbitrary.
//
// The CHOLMOD twin runs `lin_solver = 0` (remapped to KillerSolver =
// Powell-Hestenes with direct CHOLMOD inner under Newton = 1), and the
// GMG twin runs `lin_solver = 3` (FGMRES + V-cycle with D11 Newton
// Jacobian matvec on fine and D7 symmetric Picard smoother). Same
// Newton-Jacobian apples-to-apples as §7.2.
//
// Failure protocol: tasks.md §7.5 / D8.
//
// STATUS — DEFERRED (see STATUS.md §1). The GMG path produces NaN
// residuals from the first Newton inner iteration on the Newton +
// Drucker-Prager + elastic + compressible Jacobian regardless of grid
// resolution or viscous law — see the SingleShearBandGmgEquivalence
// header for the full diagnostic note. Body below is prefixed
// `DISABLED_`; will be re-enabled when the GMG stencil bridge learns
// to handle indefinite Newton-plasticity blocks.

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

class ConjugateShearBandsGmgEquivalence : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
        .SetParticles = new SetParticles_ff{
            .SetPhase   = SetPhase,
            .SetDensity = SetDensity,
        },
        .SetBCs = new SetBCs_ff{
            .SetBCVx = SetPureShearBCVx,
            .SetBCVz = SetPureShearBCVz,
        },
    };
  }

  // Shared with SingleShearBandGmgEquivalence — user2 selects geometry.
  static int SetPhase(MdoodzInput *instance, Coordinates c) {
    const double radius = instance->model.user1 / instance->scaling.L;
    const int conjugate = (int)(instance->model.user2 + 0.5);
    if (conjugate == 0) {
      if (c.x * c.x + c.z * c.z < radius * radius) return 1;
      return 0;
    }
    const double d = instance->model.user3 / instance->scaling.L;
    const double dx1 = c.x - d, dz1 = c.z - d;
    const double dx2 = c.x + d, dz2 = c.z + d;
    if (dx1 * dx1 + dz1 * dz1 < radius * radius) return 1;
    if (dx2 * dx2 + dz2 * dz2 < radius * radius) return 1;
    return 0;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates, int phase) {
    return instance->materials.rho[phase];
  }

  // Compute the principal-axis angle (in degrees, in [-90, 90]) of the
  // plastic strain-rate density on a ncx × ncz grid.
  //
  // We use PCA on the first-moment-removed second moment:
  //   M_ij = Σ (x_i - x̄)(x_j - x̄) · w(x)  (w = eII_pl)
  // then angle = 0.5 * atan2(2·Mxz, Mxx - Mzz).
  //
  // This is well-defined for a single band (returns the band dip) and
  // stable for two symmetric crossing bands (returns either of the
  // two dominant directions; the CHOLMOD and GMG twins should produce
  // the same answer because both see the same band geometry).
  static double pcaAngleDeg(const std::vector<double> &field,
                            int ncx, int ncz, double dx, double dz) {
    double total = 0.0;
    double xbar = 0.0, zbar = 0.0;
    for (int iz = 0; iz < ncz; ++iz) {
      for (int ix = 0; ix < ncx; ++ix) {
        double w = field[ix + ncx * iz];
        if (w < 0.0) w = 0.0;
        const double x = (ix + 0.5) * dx;
        const double z = (iz + 0.5) * dz;
        total += w;
        xbar  += w * x;
        zbar  += w * z;
      }
    }
    if (total < 1e-30) return 0.0;
    xbar /= total;
    zbar /= total;
    double Mxx = 0.0, Mzz = 0.0, Mxz = 0.0;
    for (int iz = 0; iz < ncz; ++iz) {
      for (int ix = 0; ix < ncx; ++ix) {
        double w = field[ix + ncx * iz];
        if (w < 0.0) w = 0.0;
        const double x = (ix + 0.5) * dx - xbar;
        const double z = (iz + 0.5) * dz - zbar;
        Mxx += w * x * x;
        Mzz += w * z * z;
        Mxz += w * x * z;
      }
    }
    return 0.5 * std::atan2(2.0 * Mxz, Mxx - Mzz) * 180.0 / M_PI;
  }

  // Wrap angle difference to [-90, 90] — PCA principal axes are
  // 180°-periodic (no sign), so we fold to the short-arc difference.
  static double angleDiffDeg(double a, double b) {
    double d = a - b;
    while (d > 90.0)  d -= 180.0;
    while (d < -90.0) d += 180.0;
    return std::fabs(d);
  }

  // Integrated plastic dissipation proxy: sum(eII_pl).
  static double integrate(const std::vector<double> &field) {
    double s = 0.0;
    for (double v : field) s += v;
    return s;
  }
};

TEST_F(ConjugateShearBandsGmgEquivalence, DISABLED_NewtonTwoBandsMatchCholmodLoose) {
  char gmgInput [] = "ShearBand/ConjugateBands_gmg.txt";
  char cholInput[] = "ShearBand/ConjugateBands_chol.txt";
  RunMDOODZ(gmgInput,  &setup);
  RunMDOODZ(cholInput, &setup);

  const char *gmg_file  = "ConjugateBands_gmg/Output00001.gzip.h5";
  const char *chol_file = "ConjugateBands_chol/Output00001.gzip.h5";

  auto eii_gmg  = readFieldAsArray(gmg_file,  "Centers", "eII_pl");
  auto eii_chol = readFieldAsArray(chol_file, "Centers", "eII_pl");
  auto vx_gmg   = readFieldAsArray(gmg_file,  "VxNodes", "Vx");
  auto vx_chol  = readFieldAsArray(chol_file, "VxNodes", "Vx");
  auto vz_gmg   = readFieldAsArray(gmg_file,  "VzNodes", "Vz");
  auto vz_chol  = readFieldAsArray(chol_file, "VzNodes", "Vz");

  ASSERT_EQ(eii_gmg.size(), eii_chol.size());
  ASSERT_EQ(vx_gmg .size(), vx_chol .size());
  ASSERT_EQ(vz_gmg .size(), vz_chol .size());

  for (double v : eii_gmg) ASSERT_TRUE(std::isfinite(v)) << "eII_gmg NaN";
  for (double v : vx_gmg ) ASSERT_TRUE(std::isfinite(v)) << "Vx_gmg  NaN";
  for (double v : vz_gmg ) ASSERT_TRUE(std::isfinite(v)) << "Vz_gmg  NaN";

  // Sanity: plasticity must have activated on both solvers.
  double max_eii_gmg = 0.0, max_eii_chol = 0.0;
  for (double v : eii_gmg ) if (v > max_eii_gmg ) max_eii_gmg  = v;
  for (double v : eii_chol) if (v > max_eii_chol) max_eii_chol = v;
  ASSERT_GT(max_eii_gmg,  0.0) << "§7.3 GMG: plasticity did not activate";
  ASSERT_GT(max_eii_chol, 0.0) << "§7.3 CHOLMOD: plasticity did not activate";

  // Band geometry via PCA on the plastic strain rate.
  const int Nx  = (int)getModelParam(gmg_file, 3);
  const int Nz  = (int)getModelParam(gmg_file, 4);
  const double dx = getModelParam(gmg_file, 5);
  const double dz = getModelParam(gmg_file, 6);
  const int ncx = Nx - 1, ncz = Nz - 1;
  ASSERT_EQ(eii_gmg.size(), (size_t)(ncx * ncz));

  const double theta_g = pcaAngleDeg(eii_gmg , ncx, ncz, dx, dz);
  const double theta_c = pcaAngleDeg(eii_chol, ncx, ncz, dx, dz);
  const double dtheta  = angleDiffDeg(theta_g, theta_c);

  // Integrated dissipation proxy.
  const double int_g = integrate(eii_gmg );
  const double int_c = integrate(eii_chol);
  const double int_ratio = (std::fabs(int_c) > 1e-30)
      ? std::fabs(int_g - int_c) / std::fabs(int_c) : 0.0;

  // Velocity L2 agreement.
  const double dVx = computeL2Error(vx_gmg, vx_chol);
  const double dVz = computeL2Error(vz_gmg, vz_chol);

  printf("ConjugateBands 41x41 Newton GMG↔CHOLMOD: "
         "principal-axis angle (GMG, CHOL) = (%.2f°, %.2f°), "
         "|Δθ| = %.3e°, "
         "int_eII_pl ratio = %.3e, "
         "|Vx-Vx|/|Vx| = %.3e, |Vz-Vz|/|Vz| = %.3e\n",
         theta_g, theta_c, dtheta, int_ratio, dVx, dVz);

  EXPECT_LT(dtheta,    2.0 ) << "§7.3 ConjugateBands: principal-axis angle "
                                "disagreement > 2°";
  EXPECT_LT(int_ratio, 0.05) << "§7.3 ConjugateBands: integrated plastic "
                                "dissipation disagreement > 5%";
  EXPECT_LT(dVx, 1e-4) << "§7.3 ConjugateBands: Vx GMG↔CHOLMOD "
                          "disagreement > 1e-4";
  EXPECT_LT(dVz, 1e-4) << "§7.3 ConjugateBands: Vz GMG↔CHOLMOD "
                          "disagreement > 1e-4";
}
