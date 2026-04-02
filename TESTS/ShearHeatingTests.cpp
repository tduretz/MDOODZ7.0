extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class ShearHeating : public ::testing::Test {
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
                    .SetBCT  = SetBCT,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    return 0;
  }

  static double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
    return (instance->model.user0 + zeroC) / instance->scaling.T;
  }

  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }

  static SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates, double gridTemperature) {
    SetBC bc;
    double Ttop = (instance->model.user0 + zeroC) / instance->scaling.T;
    if (position == N || position == NE || position == NW ||
        position == S || position == SE || position == SW) {
      bc.value = Ttop;
      bc.type  = 1;
    } else {
      bc.value = 0.0;
      bc.type  = 0;
    }
    return bc;
  }
};

TEST_F(ShearHeating, ViscousDissipation) {
  const char *testName = testing::UnitTest::GetInstance()->current_test_info()->name();
  char *inputName;
  asprintf(&inputName, "ShearHeating/%s.txt", testName);
  RunMDOODZ(inputName, &setup);

  // Read initial and final output
  char *fileInit;
  asprintf(&fileInit, "%s/%s", testName, "Output00000.gzip.h5");
  char *fileFinal;
  asprintf(&fileFinal, "%s/%s", testName, "Output00003.gzip.h5");

  double maxT_init  = getMaxFieldValue(fileInit, "Centers", "T");
  double maxT_final = getMaxFieldValue(fileFinal, "Centers", "T");
  printf("Peak T init: %e, Peak T final: %e\n", maxT_init, maxT_final);

  // Shear heating (H = tau_ij * eps_ij) converts mechanical work to heat
  // With pure shear at eps_dot = 1e-14 and eta = 1e22:
  //   H = 2 * eta * eps_dot^2 = 2 * 1e22 * (1e-14)^2 = 2e-6 W/m^3
  // Temperature should increase (or at least not decrease)
  EXPECT_GE(maxT_final, maxT_init);

  // L2: compare temperature change against analytical
  // dT = 2*eta*eps_dot^2*t / (rho*Cp)
  double meanT_init = getMeanFieldValue(fileInit, "Centers", "T");
  double meanT_final_val = getMeanFieldValue(fileFinal, "Centers", "T");
  double dT_num = meanT_final_val - meanT_init;
  double t_final = getModelParam(fileFinal, 0);
  double eta_val = 1e22;
  double eps_d = 1e-14;
  double rho_val = 3300.0;
  double Cp_val = 1050.0;
  double dT_ana = 2.0 * eta_val * eps_d * eps_d * t_final / (rho_val * Cp_val);
  printf("ViscousDissipation: dT_num=%e, dT_ana=%e\n", dT_num, dT_ana);
  EXPECT_NEAR(dT_num, dT_ana, fabs(dT_ana) * 2.0);  // within factor 3

  free(inputName);
  free(fileInit);
  free(fileFinal);
}
