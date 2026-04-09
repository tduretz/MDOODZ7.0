extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include "TestHelpers.h"

class TopoBench : public ::testing::Test {
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

  // Sinusoidal topography: h(x) = -A * cos(2*pi*x / lambda)
  static double SetSurfaceZCoord(MdoodzInput *input, double x_coord) {
    double Amplitude  = input->model.user1 / input->scaling.L;
    double Wavelength = input->model.user2 / input->scaling.L;
    return -Amplitude * cos(2.0 * M_PI * x_coord / Wavelength);
  }

  static int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    return 0;  // uniform single phase
  }

  static double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
    return (input->model.user0 + zeroC) / input->scaling.T;
  }

  static double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
  }

  // Free-slip BCs (same as TopoBenchCase1)
  static SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == S || position == SW || position == SE) {
      bc.value = 0.0;
      bc.type  = 11;
    } else if (position == N || position == NW || position == NE) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == W || position == E) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }

  static SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
    SetBC bc;
    if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
      bc.value = 0.0;
      bc.type  = 13;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type  = 0;
    } else {
      bc.value = 0.0;
      bc.type  = -1;
    }
    return bc;
  }

  static char SetBCPType(MdoodzInput *input, POSITION position) {
    if (position == NE || position == NW) {
      return 0;
    } else {
      return -1;
    }
  }

  // Analytical relaxation time for an internal density interface
  // between two equal-viscosity fluids (air has same eta as mantle in MDOODZ)
  // tau_r = (4*eta*k) / (rho*g) * coth(k*H)
  // This is 2x the free-surface formula because the viscous air doubles
  // the effective resistance to surface deflection.
  static double relaxationTime(double eta, double rho, double g, double lambda, double H) {
    double k = 2.0 * M_PI / lambda;
    return (4.0 * eta * k) / (rho * fabs(g)) * (1.0 / tanh(k * H));
  }

  // Run a TopoBench simulation and return amplitude at each output step
  // Returns pairs of (time_s, amplitude_m)
  struct TimeAmplitude {
    double time_s;
    double amplitude_m;
  };

  std::vector<TimeAmplitude> runAndExtractAmplitudes(const char *paramFile, const char *subfolder, int Nt) {
    char *input = strdup(paramFile);
    RunMDOODZ(input, &setup);
    free(input);
    std::vector<TimeAmplitude> results;
    for (int step = 1; step <= Nt; step++) {
      char fileName[256];
      snprintf(fileName, sizeof(fileName), "%s/Output%05d.gzip.h5", subfolder, step);
      double time_s = getModelParam(fileName, 0);  // Params[0] = time * scaling.t
      auto z_grid = readFieldAsArray(fileName, "Topo", "z_grid");
      // Amplitude = max(z_grid) in SI metres (already scaled in HDF5)
      double maxZ = *std::max_element(z_grid.begin(), z_grid.end());
      // Check for finite values
      for (size_t i = 0; i < z_grid.size(); i++) {
        EXPECT_TRUE(std::isfinite(z_grid[i])) << "NaN/Inf in z_grid at step " << step << " index " << i;
      }
      results.push_back({time_s, maxZ});
    }
    return results;
  }

  // Compute mean relative error of h(t)/h0 vs exp(-t/tau_r)
  static double meanRelativeError(const std::vector<TimeAmplitude> &data, double h0, double tau_r) {
    double sumErr = 0.0;
    int count = 0;
    for (size_t i = 0; i < data.size(); i++) {
      double t = data[i].time_s;
      double h_ana = h0 * exp(-t / tau_r);
      double h_num = data[i].amplitude_m;
      if (fabs(h_ana) > 1e-10) {
        sumErr += fabs(h_num - h_ana) / fabs(h_ana);
        count++;
      }
    }
    return (count > 0) ? sumErr / count : 0.0;
  }
};

TEST_F(TopoBench, TopoBenchRelaxation) {
  const char *paramFile = "TopoBench/TopoBenchRelaxation.txt";
  const char *subfolder = "TopoBenchRelaxation";
  const int Nt = 20;

  auto data = runAndExtractAmplitudes(paramFile, subfolder, Nt);
  ASSERT_EQ((int)data.size(), Nt);

  // Physical parameters (SI)
  double h0     = 7000.0;       // 7 km amplitude
  double eta    = 1e21;         // viscosity Pa·s
  double rho    = 3300.0;       // density kg/m³
  double g      = -10.0;        // gravity m/s²
  double lambda = 2800e3;       // wavelength m
  double H      = 700e3;        // domain depth m
  double tau_r  = relaxationTime(eta, rho, g, lambda, H);
  printf("Analytical tau_r = %e s (%.1f kyr)\n", tau_r, tau_r / (365.25 * 24 * 3600 * 1000));

  // Check monotonic decay
  for (size_t i = 1; i < data.size(); i++) {
    EXPECT_LT(data[i].amplitude_m, data[i - 1].amplitude_m)
            << "Amplitude not decreasing at step " << (i + 1)
            << ": " << data[i].amplitude_m << " >= " << data[i - 1].amplitude_m;
  }

  // Check relative error at each step (skip step 0 as it can have initialisation noise)
  for (size_t i = 1; i < data.size(); i++) {
    double t = data[i].time_s;
    double h_ana = h0 * exp(-t / tau_r);
    double h_num = data[i].amplitude_m;
    double relErr = fabs(h_num - h_ana) / fabs(h_ana);
    printf("Step %2zu: t=%10.3e s, h_num=%8.1f m, h_ana=%8.1f m, relErr=%.3f\n",
           i + 1, t, h_num, h_ana, relErr);
    EXPECT_LT(relErr, 0.15) << "Relative error > 15% at step " << (i + 1);
  }

  double meanErr = meanRelativeError(
          std::vector<TimeAmplitude>(data.begin() + 1, data.end()), h0, tau_r);
  printf("Mean relative error (steps 2-%d): %.4f\n", Nt, meanErr);
}

TEST_F(TopoBench, TopoBenchConvergence) {
  const char *paramFiles[] = {
          "TopoBench/TopoBenchConvergence31.txt",
          "TopoBench/TopoBenchRelaxation.txt",     // 51×26 (reuse base)
          "TopoBench/TopoBenchConvergence101.txt",
  };
  const char *subfolders[] = {
          "TopoBenchConvergence31",
          "TopoBenchRelaxation",
          "TopoBenchConvergence101",
  };
  const int Nx_vals[] = {31, 51, 101};
  const int Nt = 20;

  // Physical parameters
  double h0     = 7000.0;
  double eta    = 1e21;
  double rho    = 3300.0;
  double g      = -10.0;
  double lambda = 2800e3;
  double H      = 700e3;
  double tau_r  = relaxationTime(eta, rho, g, lambda, H);

  double errors[3];
  for (int r = 0; r < 3; r++) {
    auto data = runAndExtractAmplitudes(paramFiles[r], subfolders[r], Nt);
    ASSERT_EQ((int)data.size(), Nt);
    errors[r] = meanRelativeError(
            std::vector<TimeAmplitude>(data.begin() + 1, data.end()), h0, tau_r);
    printf("Nx=%3d: mean relative error = %.6f\n", Nx_vals[r], errors[r]);
  }

  // Monotonic error decrease
  for (int r = 1; r < 3; r++) {
    EXPECT_LT(errors[r], errors[r - 1])
            << "Error not decreasing: Nx=" << Nx_vals[r] << " (" << errors[r]
            << ") >= Nx=" << Nx_vals[r - 1] << " (" << errors[r - 1] << ")";
  }

  // Overall convergence: finest resolution should be notably better than coarsest
  double Lx = 2800e3;
  double h_coarse = Lx / Nx_vals[0];
  double h_fine   = Lx / Nx_vals[2];
  double order = log(errors[0] / errors[2]) / log(h_coarse / h_fine);
  printf("Overall order (Nx=%d->%d): %.2f\n", Nx_vals[0], Nx_vals[2], order);
  EXPECT_GE(order, 0.3) << "Overall convergence order < 0.3 from Nx=" << Nx_vals[0] << " to " << Nx_vals[2];
}
