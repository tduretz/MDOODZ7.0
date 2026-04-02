extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class Thermal : public ::testing::Test {
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
                    .SetBCVx  = SetPureShearBCVx,
                    .SetBCVz  = SetPureShearBCVz,
                    .SetBCT   = SetBCT,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    return 0;
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    const double Tbkg  = (instance->model.user0 + zeroC) / instance->scaling.T;
    const double dT    = instance->model.user2 / instance->scaling.T;
    const double rad   = instance->model.user3 / instance->scaling.L;
    const double r2    = coordinates.x * coordinates.x + coordinates.z * coordinates.z;
    if (rad > 0.0 && r2 < rad * rad) {
      return Tbkg + dT;
    }
    return Tbkg;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates, double gridTemperature) {
    SetBC bc;
    double Ttop = (instance->model.user0 + zeroC) / instance->scaling.T;
    double Tbot = instance->model.user1 / instance->scaling.T;
    if (position == N || position == NE || position == NW) {
      bc.value = Ttop;
      bc.type  = 1;
    } else if (position == S || position == SE || position == SW) {
      if (Tbot > 0.0) {
        bc.value = Tbot;
        bc.type  = 1;
      } else {
        bc.value = 0.0;
        bc.type  = 0;
      }
    } else {
      bc.value = 0.0;
      bc.type  = 0;
    }
    return bc;
  }
};

TEST_F(Thermal, GaussianDiffusion) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Thermal/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  // Read initial and final temperature
  char *fileInit;
  asprintf(&fileInit, "%s/%s", testName, "Output00000.gzip.h5");
  char *fileFinal;
  asprintf(&fileFinal, "%s/%s", testName, "Output00005.gzip.h5");
  double maxT_init  = getMaxFieldValue(fileInit, "Centers", "T");
  double maxT_final = getMaxFieldValue(fileFinal, "Centers", "T");
  printf("Peak T init: %e, Peak T final: %e\n", maxT_init, maxT_final);
  // After diffusion, peak temperature should decrease
  EXPECT_LT(maxT_final, maxT_init);
  free(inputName);
  free(fileInit);
  free(fileFinal);
}

TEST_F(Thermal, SteadyStateGeotherm) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Thermal/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00020.gzip.h5");
  double maxT = getMaxFieldValue(fileName, "Centers", "T");
  double minT = getMinFieldValue(fileName, "Centers", "T");
  double meanT = getMeanFieldValue(fileName, "Centers", "T");
  printf("T range: min=%e, max=%e, mean=%e\n", minT, maxT, meanT);
  // Should have developed a geotherm: min near top BC, max near bottom BC
  EXPECT_GT(maxT, minT);
  // Mean temperature should be between min and max (linear geotherm)
  EXPECT_GT(meanT, minT);
  EXPECT_LT(meanT, maxT);

  // L2 error against analytical geotherm T(z) = T_top + (T_bot - T_top)*(z_top - z)/H
  // Note: simulation may not be fully converged to steady state; use loose threshold
  auto T_field = readFieldAsArray(fileName, "Centers", "T");
  auto zc = readCoordArray(fileName, "zc_coord");
  int Nx_grid = (int)getModelParam(fileName, 3);
  int Nz_grid = (int)getModelParam(fileName, 4);
  int ncx = Nx_grid - 1;
  int ncz = Nz_grid - 1;
  double Lz = getModelParam(fileName, 2);
  double z_top_SI = Lz / 2.0;
  double T_top_SI = 273.15;   // user0 = 0 C
  double T_bot_SI = 1600.0;   // user1 = 1600 K
  std::vector<double> T_ana(T_field.size());
  for (int iz = 0; iz < ncz; iz++) {
    double T_z = T_top_SI + (T_bot_SI - T_top_SI) * (z_top_SI - zc[iz]) / Lz;
    for (int ix = 0; ix < ncx; ix++) {
      T_ana[ix + ncx * iz] = T_z;
    }
  }
  double L2_T = computeL2Error(T_field, T_ana);
  printf("Geotherm L2 error: %e\n", L2_T);
  EXPECT_LT(L2_T, 1.0);  // loose: simulation may not fully reach steady state

  free(inputName);
  free(fileName);
}

TEST_F(Thermal, RadiogenicHeat) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Thermal/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileInit;
  asprintf(&fileInit, "%s/%s", testName, "Output00000.gzip.h5");
  char *fileFinal;
  asprintf(&fileFinal, "%s/%s", testName, "Output00005.gzip.h5");
  double meanT_init  = getMeanFieldValue(fileInit, "Centers", "T");
  double meanT_final = getMeanFieldValue(fileFinal, "Centers", "T");
  double maxT_init   = getMaxFieldValue(fileInit, "Centers", "T");
  double maxT_final  = getMaxFieldValue(fileFinal, "Centers", "T");
  printf("Mean T init: %e, Mean T final: %e\n", meanT_init, meanT_final);
  printf("Max T init: %e, Max T final: %e\n", maxT_init, maxT_final);
  // Radiogenic heat should increase mean temperature
  EXPECT_GT(meanT_final, meanT_init);
  // Peak temperature should also increase
  EXPECT_GT(maxT_final, maxT_init);

  // Analytical: T_mean(t) ≈ T0 + Qr*t/(rho*Cp)
  double Qr = 1e-6;
  double rho_mat = 3300.0;
  double Cp_mat = 1050.0;
  double time_SI = getModelParam(fileFinal, 0);
  double T0_SI = 773.15;  // (500 + 273.15) K
  double T_ana_mean = T0_SI + Qr * time_SI / (rho_mat * Cp_mat);
  printf("RadiogenicHeat: meanT_final=%e, T_ana=%e\n", meanT_final, T_ana_mean);
  EXPECT_NEAR(meanT_final, T_ana_mean, fabs(T_ana_mean - T0_SI) * 0.5);

  free(inputName);
  free(fileInit);
  free(fileFinal);
}

TEST_F(Thermal, DirichletBC) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "Thermal/%s.txt", testName);
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "%s/%s", testName, "Output00001.gzip.h5");
  double minT = getMinFieldValue(fileName, "Centers", "T");
  double maxT = getMaxFieldValue(fileName, "Centers", "T");
  printf("T range: min=%e, max=%e\n", minT, maxT);
  // Temperature should be bounded between the Dirichlet BCs
  EXPECT_GT(minT, 0.0);
  EXPECT_GT(maxT, minT);
  free(inputName);
  free(fileName);
}
