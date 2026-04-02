extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class Density : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase       = SetPhase,
                    .SetTemperature = SetTemperature,
                    .SetDensity     = SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetPureShearBCVx,
                    .SetBCVz = SetPureShearBCVz,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (radius > 0.0 && coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    } else {
      return 0;
    }
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    // user0 = background temperature, user2 = inclusion temperature perturbation
    const double Tbkg = (instance->model.user0 + zeroC) / instance->scaling.T;
    const double radius = instance->model.user1 / instance->scaling.L;
    if (radius > 0.0 && coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      // Inclusion gets different temperature: user0 + user2
      const double Tincl = (instance->model.user0 + instance->model.user2 + zeroC) / instance->scaling.T;
      return Tincl;
    }
    return Tbkg;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }

  static int SetPhaseUniform(MdoodzInput *instance, Coordinates coordinates) {
    return 0;
  }

  static SetBC SetNoSlipBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == N || position == NW || position == NE ||
        position == S || position == SW || position == SE ||
        position == W || position == E) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }

  static SetBC SetNoSlipBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == W || position == E ||
        position == SW || position == SE || position == NW || position == NE ||
        position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }
};

TEST_F(Density, ThermalExpansion) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Density/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");

  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // With thermal expansion alpha > 0, hotter material should have lower density
  // rho(T) = rho0 * (1 - alpha * (T - T0))
  double maxRho = getMaxFieldValue(fileName, "Centers", "rho_n");
  double minRho = getMinFieldValue(fileName, "Centers", "rho_n");
  printf("Density range: min=%e, max=%e\n", minRho, maxRho);

  // There should be a density contrast (min < max) due to thermal expansion
  EXPECT_LT(minRho, maxRho);
  // Both densities should be positive and physically reasonable
  EXPECT_GT(minRho, 0.0);

  // Verify density ratio is consistent with thermal expansion: rho_hot/rho_cold = (1-alpha*T_hot)/(1-alpha*T_cold)
  double alpha = 3e-5;
  double T_cold_SI = 873.15;   // (600+273.15) K
  double T_hot_SI  = 2473.15;  // (600+1600+273.15) K
  double rho_ratio_ana = (1.0 - alpha * T_hot_SI) / (1.0 - alpha * T_cold_SI);
  double rho_ratio_num = minRho / maxRho;
  printf("ThermalExpansion: rho_ratio num=%e, ana=%e\n", rho_ratio_num, rho_ratio_ana);
  EXPECT_NEAR(rho_ratio_num, rho_ratio_ana, 0.1);

  free(inputName);
  free(fileName);
}

TEST_F(Density, HydrostaticPressure) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Density/%s.txt", testName);

  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");

  int stepsCount = getStepsCount(fileName);
  ASSERT_GE(stepsCount, 0);

  // In hydrostatic equilibrium, P should increase with depth: P ~ rho * g * z
  double maxP = getMaxFieldValue(fileName, "Centers", "P");
  double minP = getMinFieldValue(fileName, "Centers", "P");
  printf("Pressure range: min=%e, max=%e\n", minP, maxP);

  // Pressure at the bottom should be larger than at the top
  EXPECT_GT(maxP, minP);
  // Max pressure should be positive (downward gravity, increasing with depth)
  EXPECT_GT(maxP, 0.0);

  // L2 error against analytical P(z) = rho * |g| * (z_top - z)
  auto P_field = readFieldAsArray(fileName, "Centers", "P");
  auto zc = readCoordArray(fileName, "zc_coord");
  int Nx_grid = (int)getModelParam(fileName, 3);
  int Nz_grid = (int)getModelParam(fileName, 4);
  int ncx = Nx_grid - 1;
  int ncz = Nz_grid - 1;
  double Lz = getModelParam(fileName, 2);
  double z_top_SI = Lz / 2.0;
  double rho_hyd = 2700.0;
  double gz_abs = 10.0;
  std::vector<double> P_ana(P_field.size());
  for (int iz = 0; iz < ncz; iz++) {
    double P_z = rho_hyd * gz_abs * (z_top_SI - zc[iz]);
    for (int ix = 0; ix < ncx; ix++) {
      P_ana[ix + ncx * iz] = P_z;
    }
  }
  double L2_P = computeL2Error(P_field, P_ana);
  printf("Hydrostatic P L2 error: %e\n", L2_P);
  EXPECT_LT(L2_P, 2e-1);

  free(inputName);
  free(fileName);
}
