extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <gtest/gtest.h>
#include "TestHelpers.h"

// ===================== Oceanic Cooling Callbacks ===================== //

namespace OceanicCoolingCallbacks {
  int SetPhase(MdoodzInput *input, Coordinates coordinates) { return 0; }

  double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
    return (input->model.user0 + zeroC) / input->scaling.T;
  }

  double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
  }

  SetBC SetBCT(MdoodzInput *input, POSITION position, Coordinates coordinates, double gridT) {
    SetBC bc;
    double T_m = (input->model.user0 + zeroC) / input->scaling.T;
    if (position == N || position == NE || position == NW) {
      bc.value = zeroC / input->scaling.T;  // Top = 0 °C
      bc.type  = 1;
    } else if (position == S || position == SE || position == SW) {
      bc.value = T_m;  // Bottom = T_m
      bc.type  = 1;
    } else {
      bc.value = 0.0;
      bc.type  = 0;  // Neumann sides
    }
    return bc;
  }
}

// ===================== Thermo-Elastic Callbacks ===================== //

namespace ThermoElasticCallbacks {
  int SetPhase(MdoodzInput *input, Coordinates coordinates) { return 0; }

  double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
  }

  double FixTemperature(MdoodzInput *input, double pressure) {
    int thermal_evolution = (int)input->model.user0;
    if (thermal_evolution) {
      return 273.15 / input->scaling.T + 0.1 * input->model.time;
    }
    return input->model.bkg_temperature;
  }
}

// ===================== Test Fixture ===================== //

class ThermalVerification : public ::testing::Test {
protected:
  MdoodzSetup setup;
};

// ===================== Oceanic Cooling Tests ===================== //

TEST_F(ThermalVerification, OceanicCooling) {
  setup = {
    .SetParticles = new SetParticles_ff{
      .SetPhase       = OceanicCoolingCallbacks::SetPhase,
      .SetTemperature = OceanicCoolingCallbacks::SetTemperature,
      .SetDensity     = OceanicCoolingCallbacks::SetDensity,
    },
    .SetBCs = new SetBCs_ff{
      .SetBCVx = SetPureShearBCVx,
      .SetBCVz = SetPureShearBCVz,
      .SetBCT  = OceanicCoolingCallbacks::SetBCT,
    },
  };
  RunMDOODZ("Thermal/OceanicCooling.txt", &setup);

  char *fileFinal;
  asprintf(&fileFinal, "OceanicCooling/Output%05d.gzip.h5", 10);

  int Nx_grid = (int)getModelParam(fileFinal, 3);
  int Nz_grid = (int)getModelParam(fileFinal, 4);
  int ncx = Nx_grid - 1;
  int ncz = Nz_grid - 1;
  int ix_mid = ncx / 2;

  auto T_col = readColumnProfile(fileFinal, "Centers", "T", ix_mid, ncx, ncz);
  auto zc = readCoordArray(fileFinal, "zc_coord");
  double time_SI = getModelParam(fileFinal, 0);

  double k_SI   = 3.0;
  double rho_SI = 3300.0;
  double Cp_SI  = 1050.0;
  double kappa  = k_SI / (rho_SI * Cp_SI);
  double T_m    = 1400.0 + 273.15;
  double T_s    = 273.15;

  std::vector<double> T_ana(ncz);
  for (int iz = 0; iz < ncz; iz++) {
    double depth = -zc[iz];
    T_ana[iz] = T_s + (T_m - T_s) * erf(depth / sqrt(4.0 * kappa * time_SI));
  }

  double sum_sq = 0.0, sum_ref = 0.0;
  for (int iz = 0; iz < ncz; iz++) {
    double diff = T_col[iz] - T_ana[iz];
    sum_sq  += diff * diff;
    sum_ref += T_ana[iz] * T_ana[iz];
  }
  double L2 = sqrt(sum_sq / sum_ref);
  printf("OceanicCooling (CHOLMOD): L2 vs analytical = %e\n", L2);
  EXPECT_LT(L2, 5e-2);

  free(fileFinal);
}

TEST_F(ThermalVerification, OceanicCoolingPCG) {
  setup = {
    .SetParticles = new SetParticles_ff{
      .SetPhase       = OceanicCoolingCallbacks::SetPhase,
      .SetTemperature = OceanicCoolingCallbacks::SetTemperature,
      .SetDensity     = OceanicCoolingCallbacks::SetDensity,
    },
    .SetBCs = new SetBCs_ff{
      .SetBCVx = SetPureShearBCVx,
      .SetBCVz = SetPureShearBCVz,
      .SetBCT  = OceanicCoolingCallbacks::SetBCT,
    },
  };
  // Shared base fixture (Thermal/OceanicCooling.txt) + MutateInput override of
  // thermal_solver and writer_subfolder. Eliminates the OceanicCoolingPCG.txt twin.
  pcgOverride::set(1, "OceanicCoolingPCG");
  setup.MutateInput = pcgOverride::mutate;
  RunMDOODZ((char *)"Thermal/OceanicCooling.txt", &setup);
  setup.MutateInput = nullptr;
  pcgOverride::clear();

  char *fileFinal;
  asprintf(&fileFinal, "OceanicCoolingPCG/Output%05d.gzip.h5", 10);

  int Nx_grid = (int)getModelParam(fileFinal, 3);
  int Nz_grid = (int)getModelParam(fileFinal, 4);
  int ncx = Nx_grid - 1;
  int ncz = Nz_grid - 1;
  int ix_mid = ncx / 2;

  auto T_col = readColumnProfile(fileFinal, "Centers", "T", ix_mid, ncx, ncz);
  auto zc = readCoordArray(fileFinal, "zc_coord");
  double time_SI = getModelParam(fileFinal, 0);

  double k_SI   = 3.0;
  double rho_SI = 3300.0;
  double Cp_SI  = 1050.0;
  double kappa  = k_SI / (rho_SI * Cp_SI);
  double T_m    = 1400.0 + 273.15;
  double T_s    = 273.15;

  std::vector<double> T_ana(ncz);
  for (int iz = 0; iz < ncz; iz++) {
    double depth = -zc[iz];
    T_ana[iz] = T_s + (T_m - T_s) * erf(depth / sqrt(4.0 * kappa * time_SI));
  }

  double sum_sq = 0.0, sum_ref = 0.0;
  for (int iz = 0; iz < ncz; iz++) {
    double diff = T_col[iz] - T_ana[iz];
    sum_sq  += diff * diff;
    sum_ref += T_ana[iz] * T_ana[iz];
  }
  double L2_pcg = sqrt(sum_sq / sum_ref);
  printf("OceanicCoolingPCG: L2 vs analytical = %e\n", L2_pcg);
  EXPECT_LT(L2_pcg, 5e-2);

  // Cross-compare with CHOLMOD result
  char *fileCholmod;
  asprintf(&fileCholmod, "OceanicCooling/Output%05d.gzip.h5", 10);
  auto T_cholmod = readColumnProfile(fileCholmod, "Centers", "T", ix_mid, ncx, ncz);
  double sum_sq2 = 0.0, sum_ref2 = 0.0;
  for (int iz = 0; iz < ncz; iz++) {
    double diff = T_col[iz] - T_cholmod[iz];
    sum_sq2  += diff * diff;
    sum_ref2 += T_cholmod[iz] * T_cholmod[iz];
  }
  double L2_diff = sqrt(sum_sq2 / sum_ref2);
  printf("OceanicCooling CHOLMOD vs PCG L2 = %e\n", L2_diff);
  EXPECT_LT(L2_diff, 1e-6);

  free(fileFinal);
  free(fileCholmod);
}

// ===================== Thermo-Elastic Coupling Tests ===================== //

TEST_F(ThermalVerification, ThermoElastic) {
  setup = {
    .SetParticles = new SetParticles_ff{
      .SetPhase   = ThermoElasticCallbacks::SetPhase,
      .SetDensity = ThermoElasticCallbacks::SetDensity,
    },
    .SetBCs = new SetBCs_ff{
      .SetBCVx        = SetPureShearBCVx,
      .SetBCVz        = SetPureShearBCVz,
      .FixTemperature = ThermoElasticCallbacks::FixTemperature,
    },
  };
  RunMDOODZ("Thermal/ThermoElastic.txt", &setup);

  char *fileFinal;
  asprintf(&fileFinal, "ThermoElastic/Output%05d.gzip.h5", 5);

  auto T_field = readFieldAsArray(fileFinal, "Centers", "T");
  auto P_field = readFieldAsArray(fileFinal, "Centers", "P");
  int Nx_grid = (int)getModelParam(fileFinal, 3);
  int Nz_grid = (int)getModelParam(fileFinal, 4);
  int ncx = Nx_grid - 1;
  int ncz = Nz_grid - 1;

  double sumP = 0.0, sumT = 0.0;
  int n = ncx * ncz;
  for (int i = 0; i < n; i++) {
    sumP += P_field[i];
    sumT += T_field[i];
  }
  double meanP = sumP / n;
  double meanT = sumT / n;
  printf("ThermoElastic: meanT = %e K, meanP = %e Pa\n", meanT, meanP);

  // Analytical: dP = +K * alpha * dT (confined thermal expansion -> compression)
  double T_init = 273.15;
  double K_bulk = 1.0 / 1e-10;  // 10 GPa
  double alpha  = 1e-5;
  double deltaT = meanT - T_init;
  double P_ana  = K_bulk * alpha * deltaT;
  printf("ThermoElastic: dT = %e K, P_ana = %e Pa, meanP = %e Pa\n", deltaT, P_ana, meanP);

  double rel_err = fabs(meanP - P_ana) / fabs(P_ana);
  printf("ThermoElastic: relative error = %e\n", rel_err);
  EXPECT_GT(meanP, 0.0);  // Thermal expansion -> positive pressure (compression)
  EXPECT_LT(rel_err, 0.25);

  free(fileFinal);
}

TEST_F(ThermalVerification, ThermoElasticPCG) {
  setup = {
    .SetParticles = new SetParticles_ff{
      .SetPhase   = ThermoElasticCallbacks::SetPhase,
      .SetDensity = ThermoElasticCallbacks::SetDensity,
    },
    .SetBCs = new SetBCs_ff{
      .SetBCVx        = SetPureShearBCVx,
      .SetBCVz        = SetPureShearBCVz,
      .FixTemperature = ThermoElasticCallbacks::FixTemperature,
    },
  };
  pcgOverride::set(1, "ThermoElasticPCG");
  setup.MutateInput = pcgOverride::mutate;
  RunMDOODZ((char *)"Thermal/ThermoElastic.txt", &setup);
  setup.MutateInput = nullptr;
  pcgOverride::clear();

  char *fileFinal;
  asprintf(&fileFinal, "ThermoElasticPCG/Output%05d.gzip.h5", 5);

  auto T_field = readFieldAsArray(fileFinal, "Centers", "T");
  auto P_field = readFieldAsArray(fileFinal, "Centers", "P");
  int Nx_grid = (int)getModelParam(fileFinal, 3);
  int Nz_grid = (int)getModelParam(fileFinal, 4);
  int ncx = Nx_grid - 1;
  int ncz = Nz_grid - 1;

  double sumP = 0.0, sumT = 0.0;
  int n = ncx * ncz;
  for (int i = 0; i < n; i++) {
    sumP += P_field[i];
    sumT += T_field[i];
  }
  double meanP = sumP / n;
  double meanT = sumT / n;

  double T_init = 273.15;
  double K_bulk = 1.0 / 1e-10;
  double alpha  = 1e-5;
  double deltaT = meanT - T_init;
  double P_ana  = K_bulk * alpha * deltaT;
  double rel_err = fabs(meanP - P_ana) / fabs(P_ana);
  printf("ThermoElasticPCG: rel_err = %e\n", rel_err);
  EXPECT_GT(meanP, 0.0);
  EXPECT_LT(rel_err, 0.25);

  // Cross-compare with CHOLMOD
  char *fileCholmod;
  asprintf(&fileCholmod, "ThermoElastic/Output%05d.gzip.h5", 5);
  auto T_cholmod = readFieldAsArray(fileCholmod, "Centers", "T");
  auto P_cholmod = readFieldAsArray(fileCholmod, "Centers", "P");
  double T_L2 = computeL2Error(T_field, T_cholmod);
  double P_L2 = computeL2Error(P_field, P_cholmod);
  printf("ThermoElastic T L2 (CHOLMOD vs PCG) = %e\n", T_L2);
  printf("ThermoElastic P L2 (CHOLMOD vs PCG) = %e\n", P_L2);
  EXPECT_LT(T_L2, 1e-6);
  EXPECT_LT(P_L2, 1e-4);

  free(fileFinal);
  free(fileCholmod);
}
