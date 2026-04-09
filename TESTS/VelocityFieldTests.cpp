extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class VelocityField : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = SetPhase,
                    .SetDensity = SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetPureOrSimpleShearBCVx,
                    .SetBCVz = SetPureOrSimpleShearBCVz,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }
};

TEST_F(VelocityField, PureShearVelocity) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "VelocityField/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // In pure shear: Vx = eps_dot * x, Vz = -eps_dot * z
  // Maximum |Vx| should be at the East/West boundaries: |Vx_max| = eps_dot * Lx/2
  // Maximum |Vz| should be at the North/South boundaries: |Vz_max| = eps_dot * Lz/2
  double maxVx = getMaxFieldValue(fileName, "VxNodes", "Vx");
  double minVx = getMinFieldValue(fileName, "VxNodes", "Vx");
  double maxVz = getMaxFieldValue(fileName, "VzNodes", "Vz");
  double minVz = getMinFieldValue(fileName, "VzNodes", "Vz");
  printf("Vx range: [%e, %e], Vz range: [%e, %e]\n", minVx, maxVx, minVz, maxVz);

  // With bkg_strain_rate = 1 and domain [-0.5, 0.5]:
  // Expected Vx_max ~ +0.5, Vx_min ~ -0.5 (at boundaries)
  // Expected Vz_max ~ +0.5, Vz_min ~ -0.5 (at boundaries)

  // Vx should be anti-symmetric: maxVx > 0, minVx < 0
  EXPECT_GT(maxVx, 0.0);
  EXPECT_LT(minVx, 0.0);

  // Vz should be anti-symmetric: maxVz > 0 (bottom moves up? no, Vz = -eps_dot * z)
  // At z < 0 (bottom): Vz > 0 (upward), at z > 0 (top): Vz < 0
  EXPECT_GT(maxVz, 0.0);
  EXPECT_LT(minVz, 0.0);

  // Symmetry check: |maxVx| ~ |minVx| within 20% tolerance
  EXPECT_NEAR(fabs(maxVx), fabs(minVx), fabs(maxVx) * 0.2);
  EXPECT_NEAR(fabs(maxVz), fabs(minVz), fabs(maxVz) * 0.2);

  // L2 error against analytical pure shear: Vx = -eps_dot * x
  // (positive bkg_strain_rate = shortening in x, extension in z)
  auto Vx_field = readFieldAsArray(fileName, "VxNodes", "Vx");
  int Nx_grid = (int)getModelParam(fileName, 3);
  int Nz_grid = (int)getModelParam(fileName, 4);
  double dx_grid = getModelParam(fileName, 5);
  double dz_grid = getModelParam(fileName, 6);
  double Lx_grid = getModelParam(fileName, 1);
  double Lz_grid = getModelParam(fileName, 2);
  double xmin_grid = -Lx_grid / 2.0;
  double zmin_grid = -Lz_grid / 2.0;
  double eps_dot = 1.0;  // bkg_strain_rate (non-dimensional: scaling V=1, L=1)
  // Vx staggered grid: Nx × (Nz+1) entries
  int nVx_z = (int)(Vx_field.size() / Nx_grid);
  std::vector<double> Vx_ana(Vx_field.size());
  for (int iz = 0; iz < nVx_z; iz++) {
    for (int ix = 0; ix < Nx_grid; ix++) {
      double x = xmin_grid + ix * dx_grid;
      Vx_ana[ix + Nx_grid * iz] = -eps_dot * x;
    }
  }
  double L2_Vx = computeL2Error(Vx_field, Vx_ana);
  double L1_Vx = computeL1Error(Vx_field, Vx_ana);
  printf("PureShear Vx L2 error: %e, L1 error: %e\n", L2_Vx, L1_Vx);
  EXPECT_LT(L2_Vx, 3e-2);  // with 10:1 inclusion: L2 ≈ 0.024

  auto Vz_field = readFieldAsArray(fileName, "VzNodes", "Vz");
  // Vz staggered grid: (Nx+1) × Nz entries
  int nVz_x = (int)(Vz_field.size() / Nz_grid);
  std::vector<double> Vz_ana(Vz_field.size());
  for (int iz = 0; iz < Nz_grid; iz++) {
    double z = zmin_grid + iz * dz_grid;
    for (int ix = 0; ix < nVz_x; ix++) {
      Vz_ana[ix + nVz_x * iz] = eps_dot * z;
    }
  }
  double L2_Vz = computeL2Error(Vz_field, Vz_ana);
  double L1_Vz = computeL1Error(Vz_field, Vz_ana);
  printf("PureShear Vz L2 error: %e, L1 error: %e\n", L2_Vz, L1_Vz);
  EXPECT_LT(L2_Vz, 3e-2);  // with 10:1 inclusion: L2 ≈ 0.024

  free(inputName);
  free(fileName);
}

TEST_F(VelocityField, SimpleShearVelocity) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "VelocityField/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // Simple shear: Vx(z) = 2 * bkg_strain_rate * z, Vz = 0
  // BCs impose Vx = ±bkg_strain_rate*Lz at top/bottom,
  // so dVx/dz = 2*bkg_strain_rate (bkg_strain_rate = strain rate ε̇_xz = ½ dVx/dz)
  double maxVx = getMaxFieldValue(fileName, "VxNodes", "Vx");
  double minVx = getMinFieldValue(fileName, "VxNodes", "Vx");
  double maxVz = getMaxFieldValue(fileName, "VzNodes", "Vz");
  double minVz = getMinFieldValue(fileName, "VzNodes", "Vz");
  printf("Vx range: [%e, %e], Vz range: [%e, %e]\n", minVx, maxVx, minVz, maxVz);

  // Boundary Vx: bottom ≈ -1.0, top ≈ +1.0 (bkg_strain_rate * Lz = 1)
  EXPECT_NEAR(maxVx, 1.0, 0.1);
  EXPECT_NEAR(minVx, -1.0, 0.1);

  // Vz should be ~0 everywhere
  EXPECT_LT(fabs(maxVz), 1e-6);
  EXPECT_LT(fabs(minVz), 1e-6);

  // L2 error: Vx against analytical Vx(z) = bkg_strain_rate * z
  auto Vx_field = readFieldAsArray(fileName, "VxNodes", "Vx");
  int Nx_grid = (int)getModelParam(fileName, 3);
  int Nz_grid = (int)getModelParam(fileName, 4);
  double Lz_grid = getModelParam(fileName, 2);
  double dz_grid = getModelParam(fileName, 6);
  double zmin_grid = -Lz_grid / 2.0;
  double gamma_dot = 1.0;  // bkg_strain_rate; actual dVx/dz = 2*gamma_dot

  // Vx staggered grid: Nx × (Nz+1) entries
  // z-coordinates of Vx nodes: zmin - dz/2, zmin + dz/2, ..., zmax + dz/2
  int nVx_z = (int)(Vx_field.size() / Nx_grid);
  std::vector<double> Vx_ana(Vx_field.size());
  for (int iz = 0; iz < nVx_z; iz++) {
    double z = zmin_grid + iz * dz_grid - dz_grid / 2.0;
    for (int ix = 0; ix < Nx_grid; ix++) {
      Vx_ana[ix + Nx_grid * iz] = 2.0 * gamma_dot * z;
    }
  }
  double L2_Vx = computeL2Error(Vx_field, Vx_ana);
  printf("SimpleShear Vx L2 error: %e\n", L2_Vx);
  EXPECT_LT(L2_Vx, 1e-6);

  // Vz should be zero everywhere
  auto Vz_field = readFieldAsArray(fileName, "VzNodes", "Vz");
  std::vector<double> Vz_ana(Vz_field.size(), 0.0);
  double L2_Vz = computeL2Error(Vz_field, Vz_ana);
  printf("SimpleShear Vz L2 error: %e\n", L2_Vz);
  EXPECT_LT(L2_Vz, 1e-6);

  free(inputName);
  free(fileName);
}
