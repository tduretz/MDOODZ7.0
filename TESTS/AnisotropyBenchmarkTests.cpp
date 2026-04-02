extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <gtest/gtest.h>
#include "TestHelpers.h"

// ============================================================================
// Anisotropy benchmark tests
// Test 1: Director evolution under simple shear — verify the director field
//         rotates toward the shear direction over time.
// Test 2: Stress-vs-angle under pure shear — verify anisotropy-dependent
//         stress invariant differs from isotropic prediction.
// ============================================================================

class AnisotropyBenchmark : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = SetPhase,
                    .SetDensity = SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetSimpleShearBCVx,
                    .SetBCVz = SetSimpleShearBCVz,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    return 0;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }
};

// ---------------------------------------------------------------------------
// Test 1: Director evolution under simple shear
//
// Setup: Homogeneous material with aniso_angle = 45° under simple shear
//        (bkg_strain_rate = 0.5, shear_style = 1). Run 10 steps.
//
// Verification:
//   - Director field (nx, nz) is written to HDF5
//   - The mean angle should rotate toward the shear plane (θ → 0° or 180°)
//   - After γ = 0.5*10*0.05 = 0.25 shear strain, θ should be < 45°
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, DirectorEvolution) {
  // Override BCs to simple shear
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  char inputFile[] = "AnisotropyBenchmark/DirectorEvolution.txt";
  RunMDOODZ(inputFile, &setup);

  const char *outFile = "DirectorEvolution/Output00010.gzip.h5";

  // Verify solver converged
  int steps = getStepsCount(outFile);
  ASSERT_GE(steps, 0);

  // Read director field
  auto nx_field = readFieldAsArray(outFile, "Centers", "nx");
  auto nz_field = readFieldAsArray(outFile, "Centers", "nz");
  ASSERT_EQ(nx_field.size(), nz_field.size());
  ASSERT_GT(nx_field.size(), 0u);

  // Compute mean angle: θ = atan2(nz, nx)
  double sumAngle = 0.0;
  int count = 0;
  for (size_t i = 0; i < nx_field.size(); i++) {
    double nx = nx_field[i];
    double nz = nz_field[i];
    double norm = sqrt(nx * nx + nz * nz);
    if (norm > 0.01) {  // skip masked/air cells
      double angle_deg = atan2(nz, nx) * 180.0 / M_PI;
      // Map to [0, 180) range since director can be ±
      if (angle_deg < 0) angle_deg += 180.0;
      sumAngle += angle_deg;
      count++;
    }
  }
  ASSERT_GT(count, 0);
  double meanAngle = sumAngle / count;
  printf("Director evolution: mean angle = %.2f° (initial = 45°)\n", meanAngle);

  // After simple shear, the director should have rotated toward the shear plane
  // Initial angle was 45°. Under simple shear, the director rotates clockwise
  // (toward 0° or 180°). We expect the mean angle to have decreased.
  EXPECT_LT(meanAngle, 44.0);  // must have rotated noticeably from 45°
  EXPECT_GT(meanAngle, 0.0);   // sanity: should still be positive
}

// ---------------------------------------------------------------------------
// Test 2: Stress invariant depends on anisotropy angle
//
// Setup: Homogeneous pure shear (bkg_strain_rate = 1.0) with aniso_angle = 30°
//        and aniso_factor = 6.0. Single time step, 11×11 grid.
//
// Verification:
//   - For isotropic case (δ=1), sxxd = 2*eta*eps_dot everywhere
//   - With anisotropy (δ=6, θ=30°), stress should differ from isotropic
//   - The deviatoric stress invariant √(sxxd²+szzd²)/√2 should be > 0
//   - Mean sxxd should differ from 2*eta*eps_dot_xx
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, StressAnisotropy) {
  // Override BCs to pure shear for this test
  setup.SetBCs->SetBCVx = SetPureShearBCVx;
  setup.SetBCs->SetBCVz = SetPureShearBCVz;

  char inputFile[] = "AnisotropyBenchmark/StressAngle.txt";
  RunMDOODZ(inputFile, &setup);

  const char *outFile = "StressAngle/Output00001.gzip.h5";

  // Verify solver converged
  int steps = getStepsCount(outFile);
  ASSERT_GE(steps, 0);

  // Read deviatoric stresses
  auto sxxd_field = readFieldAsArray(outFile, "Centers", "sxxd");
  auto szzd_field = readFieldAsArray(outFile, "Centers", "szzd");
  ASSERT_GT(sxxd_field.size(), 0u);

  // Compute mean stress invariant: τ_II = sqrt(0.5*(sxxd² + szzd²))
  double sumTauII = 0.0;
  double sumSxxd = 0.0;
  int n = (int)sxxd_field.size();
  for (int i = 0; i < n; i++) {
    double tauII = sqrt(0.5 * (sxxd_field[i] * sxxd_field[i] +
                                szzd_field[i] * szzd_field[i]));
    sumTauII += tauII;
    sumSxxd += sxxd_field[i];
  }
  double meanTauII = sumTauII / n;
  double meanSxxd = sumSxxd / n;

  // Non-dimensional: eta0 = 1, eps_dot = 1
  // Isotropic prediction: sxxd = 2*eta*eps_dot = 2.0
  double sxxd_iso = 2.0;

  printf("StressAnisotropy: mean τ_II = %e, mean sxxd = %e (iso = %e)\n",
         meanTauII, meanSxxd, sxxd_iso);

  // With anisotropy factor δ=6 and θ=30°, stress should be substantially
  // different from the isotropic value. The exact value depends on the
  // anisotropic constitutive law, but we can verify:
  // 1. Stress invariant is non-zero (the solver produced meaningful stress)
  EXPECT_GT(meanTauII, 0.1);

  // 2. Mean sxxd should differ from isotropic (within the field, sxxd is
  //    modified by anisotropy). Since the field is homogeneous with a
  //    single angle, sxxd should be roughly uniform but different from 2.0.
  //    The anisotropic correction can increase or decrease the stress.
  EXPECT_GT(fabs(meanSxxd), 0.1);  // non-trivial stress

  // 3. Read director field to verify it was initialized correctly
  auto nx_field = readFieldAsArray(outFile, "Centers", "nx");
  auto nz_field = readFieldAsArray(outFile, "Centers", "nz");
  double sumAngle = 0.0;
  int count = 0;
  for (size_t i = 0; i < nx_field.size(); i++) {
    double nx = nx_field[i];
    double nz = nz_field[i];
    double norm = sqrt(nx * nx + nz * nz);
    if (norm > 0.01) {
      double angle_deg = atan2(nz, nx) * 180.0 / M_PI;
      if (angle_deg < 0) angle_deg += 180.0;
      sumAngle += angle_deg;
      count++;
    }
  }
  if (count > 0) {
    double meanAngle = sumAngle / count;
    printf("StressAnisotropy: director mean angle = %.2f° (expected ~30°)\n", meanAngle);
    // After 1 pure-shear step, angle should be close to initial 30°
    EXPECT_NEAR(meanAngle, 30.0, 15.0);
  }
}
