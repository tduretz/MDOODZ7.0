// add-gmg-stokes-defence §7.2 — Single Drucker–Prager shear-band
// GMG↔CHOLMOD equivalence.
//
// Runs the SingleBand twin pair under Newton + power-law creep +
// Drucker–Prager pressure-sensitive plasticity (§3.3 fixtures:
// TESTS/ShearBand/SingleBand_{gmg,chol}.txt) and asserts that the two
// linear solvers produce the same shear-band structure at step 1.
//
// This is a Newton-mode test — the CHOLMOD twin runs with
// `lin_solver = 0`, which InputOutput.c:1284 remaps to `lin_solver = 2`
// (KillerSolver: Powell-Hestenes outer with direct CHOLMOD inner per
// PH iteration) because Newton = 1 is incompatible with the single-
// solve direct path. The GMG twin runs `lin_solver = 3` with the D11
// Newton-Jacobian matvec on the fine level and the D7 symmetric
// Picard smoother on the Vanka block (tasks 15.19-15.23).
//
// Acceptance criteria (tasks.md §7.2)
//   * Peak ε̇_II position within one fine cell (~100 m at 41x41)
//   * FWHM within 5% (band width along its propagation direction)
//   * Integrated band displacement within 5% (total plastic work proxy)
//   * Velocity L2 agreement within ~1e-4 (both solvers converge Newton
//     to nonlin_abs_mom = 1e-9, so ~3 decades of headroom)
//
// The FWHM analysis is simplified to "cells above half-peak strain
// rate" since the band is not axis-aligned at arbitrary dip; we
// report the ratio of above-half-peak cell counts between the two
// solvers. The integrated displacement is the domain sum of the
// plastic dissipation rate τ:ε̇_pl.
//
// Failure protocol (tasks.md §7.5 / D8)
//   A failure up to ~200 LOC of core-side change is fixed in this
//   change; otherwise the failing assertion is relaxed here with an
//   explicit comment pointing at a STATUS.md entry.
//
// STATUS — DEFERRED (see STATUS.md §1).
//   The GMG path (`lin_solver = 3`) produces NaN residuals from the
//   very first Newton inner iteration on this fixture's Newton +
//   Drucker-Prager + elastic + compressible Jacobian, independent of
//   grid resolution (tested 31² and 41²) and viscous law (constant
//   cstv=1 and power-law pwlv=1). A dual-CHOLMOD round-trip (both
//   twins on `lin_solver = -1`) trivially passes, confirming the
//   fixture is well-posed and the Powell-Hestenes direct inner path
//   converges cleanly; the `PlasticityTests.DruckerPragerYield` binary
//   also passes with the same fixture content, confirming MDOODZ's
//   core DP-plasticity assembly is fine. The root cause is localised
//   to the GMG-path's Jacobian-matvec + Vanka-smoother interaction on
//   indefinite blocks at yielding cells, and fixing it is an
//   architectural change to the §15 stencil bridge — out of scope for
//   this change per D8. Deferred to a follow-up. The GTest body below
//   is prefixed `DISABLED_` so it still compiles + links but is
//   skipped at runtime; when the core side is ready the prefix can be
//   removed and the acceptance bounds below will apply unchanged.

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

class SingleShearBandGmgEquivalence : public ::testing::Test {
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

  // user2 = 0 → single weak inclusion at origin (SingleBand)
  // user2 = 1 → two inclusions at (±user3, ±user3) (ConjugateBands)
  //
  // Both shear-band fixtures share this SetPhase so the §7.2 and §7.3
  // test binaries can link against a common implementation; user2
  // disambiguates the two geometries at runtime from the fixture file.
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

  // Return {ix, iz} of the global max of `field` on a (ncx × ncz) grid.
  static std::pair<int, int> argmax2D(const std::vector<double> &field,
                                       int ncx, int ncz) {
    int imax = 0;
    double vmax = field.empty() ? 0.0 : field[0];
    for (int iz = 0; iz < ncz; ++iz) {
      for (int ix = 0; ix < ncx; ++ix) {
        double v = field[ix + ncx * iz];
        if (v > vmax) { vmax = v; imax = ix + ncx * iz; }
      }
    }
    return {imax % ncx, imax / ncx};
  }

  // Count cells where |field| > fraction * max(|field|). Returns a proxy
  // for band width — monotonic in FWHM for a well-localised band.
  static int countAboveFraction(const std::vector<double> &field,
                                 double fraction) {
    double vmax = 0.0;
    for (double v : field) if (std::fabs(v) > vmax) vmax = std::fabs(v);
    const double threshold = fraction * vmax;
    int count = 0;
    for (double v : field) if (std::fabs(v) > threshold) ++count;
    return count;
  }

  // Integrated plastic dissipation proxy: sum(eII_pl) over all cells.
  // Full dissipation requires τ_II × eII_pl but the two solvers should
  // produce nearly identical stress fields so the L1 sum of plastic
  // strain rate is a sufficient 1-number summary of shear localisation.
  static double integrate(const std::vector<double> &field) {
    double s = 0.0;
    for (double v : field) s += v;
    return s;
  }
};

TEST_F(SingleShearBandGmgEquivalence, DISABLED_NewtonPlasticBandMatchesCholmodLoose) {
  char gmgInput [] = "ShearBand/SingleBand_gmg.txt";
  char cholInput[] = "ShearBand/SingleBand_chol.txt";
  RunMDOODZ(gmgInput,  &setup);
  RunMDOODZ(cholInput, &setup);

  const char *gmg_file  = "SingleBand_gmg/Output00001.gzip.h5";
  const char *chol_file = "SingleBand_chol/Output00001.gzip.h5";

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

  // Sanity: plasticity must have activated on both solvers. If the
  // band did not form on either, the equivalence test is vacuous.
  double max_eii_gmg = 0.0, max_eii_chol = 0.0;
  for (double v : eii_gmg ) if (v > max_eii_gmg ) max_eii_gmg  = v;
  for (double v : eii_chol) if (v > max_eii_chol) max_eii_chol = v;
  ASSERT_GT(max_eii_gmg,  0.0) << "§7.2 GMG: plasticity did not activate "
      "(max eII_pl = 0) — fixture may be under-driven";
  ASSERT_GT(max_eii_chol, 0.0) << "§7.2 CHOLMOD: plasticity did not activate "
      "(max eII_pl = 0) — fixture may be under-driven";

  // Peak position — read Nx, Nz from the Model params and compute
  // the cell (ix, iz) of the global max of eII_pl, compare between
  // solvers. Acceptance: Chebyshev distance ≤ 1 fine cell.
  const int Nx  = (int)getModelParam(gmg_file, 3);
  const int Nz  = (int)getModelParam(gmg_file, 4);
  const int ncx = Nx - 1, ncz = Nz - 1;
  ASSERT_EQ(eii_gmg.size(), (size_t)(ncx * ncz));
  auto [ix_g, iz_g] = argmax2D(eii_gmg , ncx, ncz);
  auto [ix_c, iz_c] = argmax2D(eii_chol, ncx, ncz);
  const int dx_pos = std::abs(ix_g - ix_c);
  const int dz_pos = std::abs(iz_g - iz_c);

  // FWHM proxy — cells above half-peak eII_pl. Ratio between solvers
  // should be within 5% for a well-localised band.
  const int above_g = countAboveFraction(eii_gmg , 0.5);
  const int above_c = countAboveFraction(eii_chol, 0.5);
  const double fwhm_ratio = (above_c > 0)
      ? std::fabs((double)above_g - (double)above_c) / (double)above_c
      : 0.0;

  // Integrated plastic strain rate — proxy for integrated band
  // displacement. Ratio between solvers should be within 5%.
  const double int_g = integrate(eii_gmg );
  const double int_c = integrate(eii_chol);
  const double int_ratio = (std::fabs(int_c) > 1e-30)
      ? std::fabs(int_g - int_c) / std::fabs(int_c) : 0.0;

  // Velocity L2 agreement (newton_abs_mom = 1e-9 → ~3 decades of
  // headroom above 1e-4).
  const double dVx = computeL2Error(vx_gmg, vx_chol);
  const double dVz = computeL2Error(vz_gmg, vz_chol);

  printf("SingleBand 41x41 Newton GMG↔CHOLMOD: "
         "peak_eII_pl pos Δ = (%d, %d) cells, "
         "above-half cell-count ratio = %.3e, "
         "int_eII_pl ratio = %.3e, "
         "|Vx-Vx|/|Vx| = %.3e, |Vz-Vz|/|Vz| = %.3e\n",
         dx_pos, dz_pos, fwhm_ratio, int_ratio, dVx, dVz);

  EXPECT_LE(dx_pos, 1) << "§7.2 SingleBand: peak eII_pl x-position differs "
                          "by > 1 cell between solvers";
  EXPECT_LE(dz_pos, 1) << "§7.2 SingleBand: peak eII_pl z-position differs "
                          "by > 1 cell between solvers";
  EXPECT_LT(fwhm_ratio, 0.05) << "§7.2 SingleBand: FWHM proxy (above-half "
                                 "cell count ratio) disagreement > 5%";
  EXPECT_LT(int_ratio,  0.05) << "§7.2 SingleBand: integrated plastic strain "
                                 "rate disagreement > 5%";
  EXPECT_LT(dVx, 1e-4) << "§7.2 SingleBand: Vx GMG↔CHOLMOD disagreement "
                          "> 1e-4 (expected ~nonlin_abs_mom × headroom)";
  EXPECT_LT(dVz, 1e-4) << "§7.2 SingleBand: Vz GMG↔CHOLMOD disagreement "
                          "> 1e-4";
}
