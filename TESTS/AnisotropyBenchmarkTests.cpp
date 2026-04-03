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
// Test 1: Director evolution under simple shear — L2 error vs analytical ODE
// Test 2: Dt-convergence order for the forward Euler director integration
// Test 3: Stress-vs-angle under pure shear — analytical anisotropic stress
// ============================================================================

// ---------------------------------------------------------------------------
// Analytical helpers
// ---------------------------------------------------------------------------

// Director angle under constant simple shear rate γ̇:
//   dθ/dt = -γ̇ * cos²(θ)  =>  θ(t) = atan(tan(θ₀) - γ̇·t)
static double analyticalDirectorAngle(double theta0_rad, double gamma_dot, double t) {
  return atan(tan(theta0_rad) - gamma_dot * t);
}

// Analytical stress invariant τ_II for transverse isotropic material under
// pure shear (Exx, Ezz = -Exx, Exz = 0) with constant viscosity η,
// anisotropy factor δ, and director angle θ (angle of the director, not
// the foliation normal). The code stores the director angle and internally
// uses tet = angle - π/2 for the rotation.
static double analyticalStressTII(double theta_rad, double delta, double eta, double Exx) {
  double Ezz = -Exx;
  double Exz = 0.0;
  // Code convention: tet = angle - π/2
  double tet = theta_rad - M_PI / 2.0;
  double c = cos(tet), s = sin(tet);
  double lx2  = c * c;
  double nz2  = s * s;
  double lxlz = c * s;
  // Rotate strain rate to anisotropy frame
  double E_rot_xx =  lx2 * Exx + nz2 * Ezz + 2.0 * lxlz * Exz;
  double E_rot_zz =  nz2 * Exx + lx2 * Ezz - 2.0 * lxlz * Exz;
  double E_rot_xz = -lxlz * Exx + lxlz * Ezz + (lx2 - nz2) * Exz;
  // Apply anisotropic viscosity
  double T_rot_xx = 2.0 * eta * E_rot_xx;
  double T_rot_zz = 2.0 * eta * E_rot_zz;
  double T_rot_xz = 2.0 * eta * E_rot_xz / delta;
  // Back-rotate stress to lab frame
  double Txx =  lx2 * T_rot_xx + nz2 * T_rot_zz - 2.0 * lxlz * T_rot_xz;
  double Tzz =  nz2 * T_rot_xx + lx2 * T_rot_zz + 2.0 * lxlz * T_rot_xz;
  double Txz = lxlz * T_rot_xx - lxlz * T_rot_zz + (lx2 - nz2) * T_rot_xz;
  // τ_II = sqrt(0.5*(Txx² + Tzz²) + Txz²)
  return sqrt(0.5 * (Txx * Txx + Tzz * Tzz) + Txz * Txz);
}

// Same rotation but returns sxxd
static double analyticalStressSxxd(double theta_rad, double delta, double eta, double Exx) {
  double Ezz = -Exx;
  double Exz = 0.0;
  double tet = theta_rad - M_PI / 2.0;
  double c = cos(tet), s = sin(tet);
  double lx2  = c * c;
  double nz2  = s * s;
  double lxlz = c * s;
  double E_rot_xx =  lx2 * Exx + nz2 * Ezz + 2.0 * lxlz * Exz;
  double E_rot_zz =  nz2 * Exx + lx2 * Ezz - 2.0 * lxlz * Exz;
  double E_rot_xz = -lxlz * Exx + lxlz * Ezz + (lx2 - nz2) * Exz;
  double T_rot_xx = 2.0 * eta * E_rot_xx;
  double T_rot_zz = 2.0 * eta * E_rot_zz;
  double T_rot_xz = 2.0 * eta * E_rot_xz / delta;
  return lx2 * T_rot_xx + nz2 * T_rot_zz - 2.0 * lxlz * T_rot_xz;
}

// Same rotation but returns sxz
static double analyticalStressSxz(double theta_rad, double delta, double eta, double Exx) {
  double Ezz = -Exx;
  double Exz = 0.0;
  double tet = theta_rad - M_PI / 2.0;
  double c = cos(tet), s = sin(tet);
  double lx2  = c * c;
  double nz2  = s * s;
  double lxlz = c * s;
  double E_rot_xx =  lx2 * Exx + nz2 * Ezz + 2.0 * lxlz * Exz;
  double E_rot_zz =  nz2 * Exx + lx2 * Ezz - 2.0 * lxlz * Exz;
  double E_rot_xz = -lxlz * Exx + lxlz * Ezz + (lx2 - nz2) * Exz;
  double T_rot_xx = 2.0 * eta * E_rot_xx;
  double T_rot_zz = 2.0 * eta * E_rot_zz;
  double T_rot_xz = 2.0 * eta * E_rot_xz / delta;
  return lxlz * T_rot_xx - lxlz * T_rot_zz + (lx2 - nz2) * T_rot_xz;
}

// Compute L2 error of director angle field vs a single analytical value
// Returns the RMS angular error in radians
static double directorAngleL2(const std::vector<double> &nx_field,
                               const std::vector<double> &nz_field,
                               double theta_analytical_rad) {
  double sumErr2 = 0.0;
  int count = 0;
  for (size_t i = 0; i < nx_field.size(); i++) {
    double nx = nx_field[i];
    double nz = nz_field[i];
    double norm = sqrt(nx * nx + nz * nz);
    if (norm > 0.01) {
      double theta_num = atan2(nz, nx);
      // Handle director ambiguity: ±n gives angles differing by π
      double diff = theta_num - theta_analytical_rad;
      // Wrap to [-π, π]
      while (diff >  M_PI) diff -= 2.0 * M_PI;
      while (diff < -M_PI) diff += 2.0 * M_PI;
      // Also check the antipodal director
      double diff2 = theta_num + M_PI - theta_analytical_rad;
      while (diff2 >  M_PI) diff2 -= 2.0 * M_PI;
      while (diff2 < -M_PI) diff2 += 2.0 * M_PI;
      double minDiff = (fabs(diff) < fabs(diff2)) ? diff : diff2;
      sumErr2 += minDiff * minDiff;
      count++;
    }
  }
  if (count == 0) return 1e30;
  return sqrt(sumErr2 / count);
}

// Extract mean director angle (degrees) from an HDF5 output file
static double meanDirectorAngleDeg(const char *hdf5File) {
  auto nx_field = readFieldAsArray(hdf5File, "Centers", "nx");
  auto nz_field = readFieldAsArray(hdf5File, "Centers", "nz");
  double sumTheta = 0.0;
  int count = 0;
  for (size_t i = 0; i < nx_field.size(); i++) {
    double nx = nx_field[i], nz = nz_field[i];
    double norm = sqrt(nx * nx + nz * nz);
    if (norm > 0.01) {
      sumTheta += atan2(nz, nx) * 180.0 / M_PI;
      count++;
    }
  }
  return (count > 0) ? sumTheta / count : 0.0;
}

// Run a director evolution simulation and return the L2 angular error
static double runDirectorAndGetL2(MdoodzSetup *setup, const char *inputFile,
                                   const char *outFile,
                                   double theta0_rad, double gamma_dot,
                                   double total_time) {
  char buf[256];
  snprintf(buf, sizeof(buf), "%s", inputFile);
  RunMDOODZ(buf, setup);

  auto nx_field = readFieldAsArray(outFile, "Centers", "nx");
  auto nz_field = readFieldAsArray(outFile, "Centers", "nz");

  double theta_ana = analyticalDirectorAngle(theta0_rad, gamma_dot, total_time);
  double L2 = directorAngleL2(nx_field, nz_field, theta_ana);

  printf("  θ₀=%.1f° γ̇=%.2f t=%.3f → θ_ana=%.4f° L2=%.6e rad\n",
         theta0_rad * 180.0 / M_PI, gamma_dot, total_time,
         theta_ana * 180.0 / M_PI, L2);
  return L2;
}

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
// Test 1: Director evolution under simple shear — L2 error vs analytical
//
// Setup: Homogeneous material with aniso_angle = 45° under simple shear
//        (bkg_strain_rate = 0.5, shear_style = 1 → γ̇ = 1.0). Run 40 steps.
//
// Analytical: θ(t) = arctan(tan(45°) - 1.0·0.5) = arctan(0.5) ≈ 26.57°
// Golden standard: dt=0.0125 gives L2 ≈ 2.34e-3 rad (0.13° error on 18.4° rotation)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, DirectorEvolution) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double theta0     = 45.0 * M_PI / 180.0;  // initial angle
  const double gamma_dot  = 1.0;                    // dudz = 2 * bkg_strain_rate
  const double total_time = 40 * 0.0125;            // Nt * dt = 0.5

  double L2 = runDirectorAndGetL2(&setup,
    "AnisotropyBenchmark/DirectorEvolution_dt0125.txt",
    "DirectorEvolution_dt0125/Output00040.gzip.h5",
    theta0, gamma_dot, total_time);

  EXPECT_LT(L2, 5e-3) << "Director L2 error (radians) exceeds threshold";
}

// ---------------------------------------------------------------------------
// Test 2: Director dt-convergence order (5-point study)
//
// Run the same total time (t=0.5) at 5 dt values: 0.25, 0.1, 0.05, 0.025,
// 0.0125. The forward Euler scheme should converge at order ~1.0.
// Writes director_convergence.dat for gnuplot visualization.
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, DirectorDtConvergence) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double theta0     = 45.0 * M_PI / 180.0;
  const double gamma_dot  = 1.0;
  const double total_time = 0.5;

  const int N_DT = 5;
  const double dts[N_DT] = {0.25, 0.1, 0.05, 0.025, 0.0125};
  const int    Nts[N_DT] = {2, 5, 10, 20, 40};
  const char *inputFiles[N_DT] = {
    "AnisotropyBenchmark/DirectorEvolution_dt25.txt",
    "AnisotropyBenchmark/DirectorEvolution_dt01.txt",
    "AnisotropyBenchmark/DirectorEvolution.txt",
    "AnisotropyBenchmark/DirectorEvolution_dt025.txt",
    "AnisotropyBenchmark/DirectorEvolution_dt0125.txt"
  };
  const char *outFiles[N_DT] = {
    "DirectorEvolution_dt25/Output00002.gzip.h5",
    "DirectorEvolution_dt01/Output00005.gzip.h5",
    "DirectorEvolution/Output00010.gzip.h5",
    "DirectorEvolution_dt025/Output00020.gzip.h5",
    "DirectorEvolution_dt0125/Output00040.gzip.h5"
  };
  const char *outDirs[N_DT] = {
    "DirectorEvolution_dt25",
    "DirectorEvolution_dt01",
    "DirectorEvolution",
    "DirectorEvolution_dt025",
    "DirectorEvolution_dt0125"
  };
  const char *trajFiles[N_DT] = {
    "director_trajectory_dt250.dat",
    "director_trajectory_dt100.dat",
    "director_trajectory_dt050.dat",
    "director_trajectory_dt025.dat",
    "director_trajectory_dt0125.dat"
  };

  double L2s[N_DT];
  for (int k = 0; k < N_DT; k++) {
    L2s[k] = runDirectorAndGetL2(&setup, inputFiles[k], outFiles[k],
                                  theta0, gamma_dot, total_time);
  }

  // Monotonic error decrease
  for (int k = 1; k < N_DT; k++) {
    EXPECT_LT(L2s[k], L2s[k-1]) << "L2 should decrease: dt=" << dts[k]
                                  << " vs dt=" << dts[k-1];
  }

  // Print all pair-wise convergence orders
  printf("Dt-convergence (5-point):\n");
  for (int k = 0; k < N_DT; k++) {
    printf("  dt=%.4f  L2=%.6e rad\n", dts[k], L2s[k]);
  }
  for (int k = 1; k < N_DT; k++) {
    double o = log(L2s[k-1] / L2s[k]) / log(dts[k-1] / dts[k]);
    printf("  Order(%.4f → %.4f): %.2f\n", dts[k-1], dts[k], o);
  }

  // Assert convergence order from finest pair
  double order = log(L2s[N_DT-2] / L2s[N_DT-1]) / log(dts[N_DT-2] / dts[N_DT-1]);
  EXPECT_GE(order, 0.8) << "Forward Euler should converge at ~1st order";

  // Write convergence data for gnuplot (optional visualization)
  FILE *fp = fopen("director_convergence.dat", "w");
  if (fp) {
    fprintf(fp, "# dt          L2_rad\n");
    for (int k = 0; k < N_DT; k++) {
      fprintf(fp, "%.6e  %.6e\n", dts[k], L2s[k]);
    }
    fclose(fp);
    printf("Wrote: director_convergence.dat\n");
  }

  // Write trajectory data files from real MDOODZ HDF5 outputs for gnuplot
  const double theta0_deg = theta0 * 180.0 / M_PI;
  for (int k = 0; k < N_DT; k++) {
    FILE *tf = fopen(trajFiles[k], "w");
    if (!tf) continue;
    fprintf(tf, "# t  theta_deg  (MDOODZ dt=%.4f, Nt=%d)\n", dts[k], Nts[k]);
    fprintf(tf, "%.6e  %.6f\n", 0.0, theta0_deg);  // initial condition
    for (int step = 1; step <= Nts[k]; step++) {
      char hdf5Path[256];
      snprintf(hdf5Path, sizeof(hdf5Path), "%s/Output%05d.gzip.h5", outDirs[k], step);
      double theta_deg = meanDirectorAngleDeg(hdf5Path);
      fprintf(tf, "%.6e  %.6f\n", step * dts[k], theta_deg);
    }
    fclose(tf);
    printf("Wrote: %s (%d points)\n", trajFiles[k], Nts[k] + 1);
  }
}

// ---------------------------------------------------------------------------
// Test 3: Stress invariant under anisotropic pure shear — analytical check
//
// Setup: Homogeneous pure shear (bkg_strain_rate = 1.0) with aniso_angle = 30°
//        and aniso_factor = 6.0. eta0 = 1. Single time step.
//
// Analytical: rotate strain rate to anisotropy frame, apply σ_ns/δ, rotate back
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, StressAnisotropy) {
  setup.SetBCs->SetBCVx = SetPureShearBCVx;
  setup.SetBCs->SetBCVz = SetPureShearBCVz;

  char inputFile[] = "AnisotropyBenchmark/StressAngle.txt";
  RunMDOODZ(inputFile, &setup);

  const char *outFile = "StressAngle/Output00001.gzip.h5";

  int steps = getStepsCount(outFile);
  ASSERT_GE(steps, 0);

  // Analytical values for pure shear: Exx = -1.0, η = 1.0, θ = 30°, δ = 6
  const double theta_deg = 30.0;
  const double theta_rad = theta_deg * M_PI / 180.0;
  const double delta     = 6.0;
  const double eta       = 1.0;
  const double Exx       = -1.0;  // shortening in x

  double tauII_ana = analyticalStressTII(theta_rad, delta, eta, Exx);
  double sxxd_ana  = analyticalStressSxxd(theta_rad, delta, eta, Exx);

  // Read numerical stresses
  auto sxxd_field = readFieldAsArray(outFile, "Centers", "sxxd");
  auto szzd_field = readFieldAsArray(outFile, "Centers", "szzd");
  auto sxz_field  = readFieldAsArray(outFile, "Vertices", "sxz");
  ASSERT_GT(sxxd_field.size(), 0u);
  ASSERT_GT(sxz_field.size(), 0u);

  double sumSxxd = 0.0, sumSzzd = 0.0;
  int nc = (int)sxxd_field.size();
  for (int i = 0; i < nc; i++) {
    sumSxxd += sxxd_field[i];
    sumSzzd += szzd_field[i];
  }
  double meanSxxd = sumSxxd / nc;
  double meanSzzd = sumSzzd / nc;

  double sumSxz = 0.0;
  int nv = (int)sxz_field.size();
  for (int i = 0; i < nv; i++) {
    sumSxz += sxz_field[i];
  }
  double meanSxz = sumSxz / nv;

  // Full τ_II = sqrt(0.5*(sxxd² + szzd²) + sxz²)
  double meanTauII = sqrt(0.5 * (meanSxxd * meanSxxd + meanSzzd * meanSzzd)
                          + meanSxz * meanSxz);

  // Analytical sxz
  double sxz_ana = analyticalStressSxz(theta_rad, delta, eta, Exx);

  printf("StressAnisotropy: τ_II num=%.6e ana=%.6e\n", meanTauII, tauII_ana);
  printf("StressAnisotropy: sxxd num=%.6e ana=%.6e\n", meanSxxd, sxxd_ana);
  printf("StressAnisotropy: sxz  num=%.6e ana=%.6e\n", meanSxz, sxz_ana);

  double relErr_tauII = fabs(meanTauII - tauII_ana) / fabs(tauII_ana);
  double relErr_sxxd  = fabs(meanSxxd - sxxd_ana) / fabs(sxxd_ana);
  printf("StressAnisotropy: relErr(τ_II)=%.4e relErr(sxxd)=%.4e\n",
         relErr_tauII, relErr_sxxd);

  EXPECT_LT(relErr_tauII, 1e-4) << "τ_II relative error exceeds 0.01%";
  EXPECT_LT(relErr_sxxd, 1e-4)  << "sxxd relative error exceeds 0.01%";
}
