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

// Closure-style state for `mutateDirectorInput`. MDLIB's MutateInput receives
// only `MdoodzInput*`, so a static struct is the simplest way to thread the
// per-iteration override values through the callback. Tests run sequentially
// within a binary → no race.
namespace {
  struct DirectorOverride {
    double      dt;        // 0.0 → no override
    int         Nt;        // 0   → no override
    const char *subfolder; // nullptr → no override
  };
  static DirectorOverride g_dirOverride{0.0, 0, nullptr};

  // writer_subfolder is heap-owned by ReadChar and freed at sim end —
  // free() the original before strdup-ing the override (string literal would
  // crash free() on shutdown).
  static void mutateDirectorInput(MdoodzInput *input) {
    if (g_dirOverride.dt != 0.0) {
      input->model.dt       = g_dirOverride.dt;
      input->model.dt0      = g_dirOverride.dt;  // dt0 set from dt at parse (InputOutput.c:1301)
      input->model.dt_start = g_dirOverride.dt;  // dt_start used by AdvectionRoutines.c:685 to reset model.dt — must mirror
    }
    if (g_dirOverride.Nt != 0)   input->model.Nt = g_dirOverride.Nt;
    if (g_dirOverride.subfolder != nullptr) {
      free(input->model.writer_subfolder);
      input->model.writer_subfolder = strdup(g_dirOverride.subfolder);
    }
  }

  // Overrides for the AniFstrainEvolution base fixture. Three test cases share
  // a single .txt and inject their differentiating parameters here.
  // bkg_strain_rate is in the .txt's input units (non-dim, since the fixture
  // uses scales = 1). All zero/-1/nullptr fields → no override.
  // ani_fac_max is per-phase (phase 0 only — single-phase fixture).
  struct AniFstrainOverride {
    double      bkg_strain_rate;  // 0.0     → no override
    int         Nt;               // 0       → no override
    double      dt;               // 0.0     → no override (sets dt, dt0, dt_start together)
    double      ani_fac_max;      // 0.0     → no override
    int         shear_style;      // -1      → no override (0 = pure shear, 1 = simple)
    int         periodic_x;       // -1      → no override
    int         pure_shear_ALE;   // -99     → no override (valid values include -1)
    const char *subfolder;        // nullptr → no override
  };
  static AniFstrainOverride g_aniFstrainOverride{0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  static void mutateAniFstrainInput(MdoodzInput *input) {
    if (g_aniFstrainOverride.bkg_strain_rate != 0.0)
      input->model.bkg_strain_rate = g_aniFstrainOverride.bkg_strain_rate;
    if (g_aniFstrainOverride.Nt != 0)
      input->model.Nt = g_aniFstrainOverride.Nt;
    if (g_aniFstrainOverride.dt != 0.0) {
      input->model.dt       = g_aniFstrainOverride.dt;
      input->model.dt0      = g_aniFstrainOverride.dt;
      input->model.dt_start = g_aniFstrainOverride.dt;
    }
    if (g_aniFstrainOverride.ani_fac_max != 0.0)
      input->materials.ani_fac_max[0] = g_aniFstrainOverride.ani_fac_max;
    if (g_aniFstrainOverride.shear_style >= 0)
      input->model.shear_style = g_aniFstrainOverride.shear_style;
    if (g_aniFstrainOverride.periodic_x >= 0)
      input->model.periodic_x = g_aniFstrainOverride.periodic_x;
    if (g_aniFstrainOverride.pure_shear_ALE != -99)
      input->model.pure_shear_ALE = g_aniFstrainOverride.pure_shear_ALE;
    if (g_aniFstrainOverride.subfolder != nullptr) {
      free(input->model.writer_subfolder);
      input->model.writer_subfolder = strdup(g_aniFstrainOverride.subfolder);
    }
  }
}

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

  // Inject dt=0.0125, Nt=40 via MutateInput on the shared base fixture
  // (DirectorEvolution.txt). Eliminates the need for a separate
  // DirectorEvolution_dt0125.txt fixture.
  setup.MutateInput = mutateDirectorInput;
  g_dirOverride = { 0.0125, 40, "DirectorEvolution_dt0125" };
  double L2 = runDirectorAndGetL2(&setup,
    "AnisotropyBenchmark/DirectorEvolution.txt",
    "DirectorEvolution_dt0125/Output00040.gzip.h5",
    theta0, gamma_dot, total_time);
  setup.MutateInput = nullptr;
  g_dirOverride = {0.0, 0, nullptr};

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

  // All 5 dt values share a single base fixture (`DirectorEvolution.txt`).
  // dt, Nt and writer_subfolder are injected per-iteration via MutateInput —
  // eliminates 4 redundant fixture files.
  setup.MutateInput = mutateDirectorInput;
  const char *baseTxt = "AnisotropyBenchmark/DirectorEvolution.txt";

  double L2s[N_DT];
  for (int k = 0; k < N_DT; k++) {
    char outFile[256];
    snprintf(outFile, sizeof(outFile), "%s/Output%05d.gzip.h5", outDirs[k], Nts[k]);
    g_dirOverride = { dts[k], Nts[k], outDirs[k] };
    L2s[k] = runDirectorAndGetL2(&setup, baseTxt, outFile,
                                  theta0, gamma_dot, total_time);
  }
  setup.MutateInput = nullptr;
  g_dirOverride = {0.0, 0, nullptr};

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

// ---------------------------------------------------------------------------
// Test 4: Spatial L2 error norm for anisotropic stress components
//
// Same setup as Test 3 (StressAngle.txt: 11×11, θ=30°, δ=6, η=1, Exx=−1).
// Computes L2 error of sxxd, szzd, sxz against the spatially constant
// analytical values (homogeneous problem → same value at every grid point).
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, StressAnisotropyL2) {
  setup.SetBCs->SetBCVx = SetPureShearBCVx;
  setup.SetBCs->SetBCVz = SetPureShearBCVz;

  char inputFile[] = "AnisotropyBenchmark/StressAngle.txt";
  RunMDOODZ(inputFile, &setup);

  const char *outFile = "StressAngle/Output00001.gzip.h5";

  // Analytical values
  const double theta_rad = 30.0 * M_PI / 180.0;
  const double delta     = 6.0;
  const double eta       = 1.0;
  const double Exx       = -1.0;

  double sxxd_ana = analyticalStressSxxd(theta_rad, delta, eta, Exx);
  double sxz_ana  = analyticalStressSxz(theta_rad, delta, eta, Exx);
  double szzd_ana = -sxxd_ana;  // incompressibility: szzd = -sxxd

  // Read numerical fields
  auto sxxd_field = readFieldAsArray(outFile, "Centers", "sxxd");
  auto szzd_field = readFieldAsArray(outFile, "Centers", "szzd");
  auto sxz_field  = readFieldAsArray(outFile, "Vertices", "sxz");
  ASSERT_GT(sxxd_field.size(), 0u);
  ASSERT_GT(szzd_field.size(), 0u);
  ASSERT_GT(sxz_field.size(), 0u);

  // Build constant analytical vectors at matching grid sizes
  std::vector<double> sxxd_ana_vec(sxxd_field.size(), sxxd_ana);
  std::vector<double> szzd_ana_vec(szzd_field.size(), szzd_ana);
  std::vector<double> sxz_ana_vec(sxz_field.size(), sxz_ana);

  double l2_sxxd = computeL2Error(sxxd_field, sxxd_ana_vec);
  double l2_szzd = computeL2Error(szzd_field, szzd_ana_vec);
  double l2_sxz  = computeL2Error(sxz_field, sxz_ana_vec);

  printf("StressAnisotropyL2: sxxd L2=%.4e (ana=%.6e, N=%zu)\n",
         l2_sxxd, sxxd_ana, sxxd_field.size());
  printf("StressAnisotropyL2: szzd L2=%.4e (ana=%.6e, N=%zu)\n",
         l2_szzd, szzd_ana, szzd_field.size());
  printf("StressAnisotropyL2: sxz  L2=%.4e (ana=%.6e, N=%zu)\n",
         l2_sxz, sxz_ana, sxz_field.size());

  EXPECT_LT(l2_sxxd, 1e-6) << "sxxd spatial L2 error exceeds threshold";
  EXPECT_LT(l2_szzd, 1e-6) << "szzd spatial L2 error exceeds threshold";
  EXPECT_LT(l2_sxz, 1e-6)  << "sxz spatial L2 error exceeds threshold";
}


// --- Interpolation mode tests ---

TEST_F(AnisotropyBenchmark, DirectorEvolutionInterpMode1) {
  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/DirectorEvolutionInterpMode1.txt");
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "DirectorEvolutionInterpMode1/Output%05d.gzip.h5", 5);
  int steps = getStepsCount(fileName);
  EXPECT_GE(steps, 0);
  free(inputName);
  free(fileName);
}

TEST_F(AnisotropyBenchmark, DirectorEvolutionInterpMode2) {
  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/DirectorEvolutionInterpMode2.txt");
  RunMDOODZ(inputName, &setup);
  char *fileName;
  asprintf(&fileName, "DirectorEvolutionInterpMode2/Output%05d.gzip.h5", 5);
  int steps = getStepsCount(fileName);
  EXPECT_GE(steps, 0);
  free(inputName);
  free(fileName);
}

// ---------------------------------------------------------------------------
// Finite-strain anisotropy benchmarks (ani_fstrain = 1)
//
// The per-phase switch ani_fstrain = 1 selects evolving anisotropy strength
// δ = min(FS_AR, ani_fac_max) where FS_AR = e1/e2 is the strain-ellipse aspect
// ratio computed from the polar decomposition of the deformation gradient F.
//
// Closed-form references (derived in TESTS/AnalyticalSolutions.md):
//   pure  shear: F(t) = diag(exp(-ε̇·t), exp(ε̇·t)) → FS_AR(t) = exp(2·ε̇·t)
//   simple shear: F(t) = [[1, γ], [0, 1]] (γ = γ̇·t)
//                 → FS_AR(γ) = 1 + γ²/2 + γ·sqrt(γ²/4 + 1)
//
// MDOODZ BC convention:
//   pure  shear (SetPureShearBC): Vx = -x·bkg_strain_rate, Vz = +z·bkg_strain_rate
//                                  ⇒ ε̇ = bkg_strain_rate
//   simple shear (SetSimpleShearBC): Vx(top/bot) = ±bkg_strain_rate·Lz
//                                     ⇒ γ̇ = 2·bkg_strain_rate
// ---------------------------------------------------------------------------

static double anaFsArPureShear(double Eii, double t) {
  return exp(2.0 * Eii * t);
}

static double anaFsArSimpleShear(double gamma_dot, double t) {
  const double g  = gamma_dot * t;
  return 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
}

// Read all interior `Centers/ani_fac` values, return (mean, min, max, L2 vs constant).
struct AniFacStats {
  double mean;
  double min_v;
  double max_v;
  double l2_rel;
  size_t n;
};

static AniFacStats readAniFacStats(const char *hdf5File, double analytical) {
  auto field = readFieldAsArray(hdf5File, "Centers", "ani_fac");
  AniFacStats s{0.0, DBL_MAX, -DBL_MAX, 0.0, field.size()};
  std::vector<double> ana_vec(field.size(), analytical);
  double sum = 0.0;
  for (double v : field) {
    sum += v;
    if (v < s.min_v) s.min_v = v;
    if (v > s.max_v) s.max_v = v;
  }
  s.mean   = sum / (double)field.size();
  s.l2_rel = computeL2Error(field, ana_vec);
  return s;
}

// ---------------------------------------------------------------------------
// Test: AniFstrainPureShearSaturation — δ pinned at ani_fac_max for t > t_sat
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainSaturation) {
  // Saturation test: simple shear with a small ani_fac_max so the clamp kicks
  // in within a few steps. Simple shear is used (not pure shear) because pure
  // shear with ALE is unstable at the per-step strain rates needed to reach
  // saturation in a short test. The min(FS_AR, ani_fac_max) clamp itself is
  // deformation-mode-agnostic, so simple shear is sufficient to gate it.
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;             // → γ̇ = 2·bkg_sr = 2
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 5;
  const double dt          = 0.5;
  const double ani_fac_max = 4.0;
  // Saturation strain (γ_sat at FS_AR = ani_fac_max, simple-shear formula
  // FS_AR = 1 + γ²/2 + γ·sqrt(γ²/4+1) inverted): ani_fac_max=4 → γ_sat ≈ 1.5.
  const double gamma_sat   = sqrt(pow(ani_fac_max, 2.0) - 1.0)
                             / sqrt(ani_fac_max);

  setup.MutateInput  = mutateAniFstrainInput;
  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainSaturation" };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  printf("Saturation regime: γ̇=%.2f ani_fac_max=%.2f γ_sat≈%.4f\n",
         gamma_dot, ani_fac_max, gamma_sat);

  // .dat file for the CI plot — the clamp is the interesting visual signature
  // of the min() form, so we emit measurements from the saturation test.
  FILE *fp = fopen("aniso_factor_evolution.dat", "w");
  if (fp) {
    fprintf(fp,
            "# step  t          gamma       FS_AR_ana    delta_clamped_ana    "
            "delta_mdoodz_mean   ani_fac_max=%.4e gamma_dot=%.4f gamma_sat=%.4f\n",
            ani_fac_max, gamma_dot, gamma_sat);
    fprintf(fp, "%d  %.6e  %.6e  %.6e  %.6e  %.6e\n",
            0, 0.0, 0.0, 1.0, 1.0, 1.0);
  }

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainSaturation/Output%05d.gzip.h5", step);
    // Output<k> reflects F at end of step k-1 → analytical at γ = γ̇·(k-1)·dt
    const double t        = (step - 1) * dt;
    const double gamma    = gamma_dot * t;
    const double fs_ar    = anaFsArSimpleShear(gamma_dot, t);
    const double ana      = (fs_ar < ani_fac_max) ? fs_ar : ani_fac_max;
    AniFacStats  s        = readAniFacStats(hdf5Path, ana);
    const bool   saturated = (gamma > gamma_sat);

    printf("  step=%d γ=%.3f FS_AR=%.4f δ_ana=%.4f δ_mean=%.6e %s\n",
           step, gamma, fs_ar, ana, s.mean,
           saturated ? "(saturated)" : "(unsaturated)");

    if (saturated) {
      EXPECT_NEAR(s.mean, ani_fac_max, 1e-9)
          << "Saturated step " << step << ": δ must equal ani_fac_max";
    } else {
      EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
          << "Unsaturated step " << step << ": δ vs FS_AR(γ)";
    }
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "Saturation step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));

    if (fp) {
      fprintf(fp, "%d  %.6e  %.6e  %.6e  %.6e  %.6e\n",
              step, t, gamma, fs_ar, ana, s.mean);
    }
  }

  if (fp) {
    fclose(fp);
    printf("Wrote: aniso_factor_evolution.dat (Nt+1 points, ani_fac_max=%.2f)\n",
           ani_fac_max);
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainSimpleShear — δ(γ) = 1 + γ²/2 + γ·√(γ²/4+1), unsaturated.
// Also writes aniso_factor_evolution.dat for the gnuplot CI visualisation.
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainSimpleShear) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr     = 1.0;             // → γ̇ = 2·bkg_sr = 2
  const double gamma_dot  = 2.0 * bkg_sr;
  const int    Nt         = 5;
  const double dt         = 0.5;
  const double ani_fac_max = 1.0e6;

  setup.MutateInput  = mutateAniFstrainInput;
  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainSimpleShear" };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  // Output<k> reflects F at end of step k-1 → reference at γ = γ̇·(k-1)·dt.
  // The .dat file for the CI plot is written by AniFstrainSaturation (which
  // exhibits the clamp — more interesting visual than the unbounded curve).
  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainSimpleShear/Output%05d.gzip.h5", step);
    const double t      = (step - 1) * dt;
    const double gamma  = gamma_dot * t;
    const double fs_ar  = anaFsArSimpleShear(gamma_dot, t);
    const double ana    = (fs_ar < ani_fac_max) ? fs_ar : ani_fac_max;
    AniFacStats  s      = readAniFacStats(hdf5Path, ana);

    printf("SimpleShear step=%d t=%.3f γ=%.3f FS_AR=%.4f δ_mean=%.6e "
           "L2=%.4e\n",
           step, t, gamma, fs_ar, s.mean, s.l2_rel);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Simple shear step " << step
        << ": mean δ vs 1 + γ²/2 + γ·√(γ²/4+1)";
    EXPECT_LT(s.l2_rel, 1e-4) << "Simple shear step " << step << ": L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "Simple shear step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainPureShear — δ(t) = exp(2·ε̇·t), unsaturated.
//
// Pure-shear F(t) = diag(exp(-ε̇·t), exp(ε̇·t)) → FS_AR = exp(2·ε̇·t).
// With ani_fac_max = 1e6 the test stays unsaturated across all Nt steps.
//
// Output<k>.gzip.h5 reflects FS_AR computed at the START of step k from
// particles.F at the END of step k-1, so the analytical reference for
// Output<k> is exp(2·ε̇·(k-1)·dt).
//
// Eii = 0.1 keeps per-step strain at 5 % so the post-ALE-shrink mesh stays
// reasonably sized; at this strain rate, FS_AR_final = exp(2·0.1·2.5) ≈ 1.65,
// a small range but enough to verify the formula at FP precision.
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainPureShear) {
  setup.SetBCs->SetBCVx = SetPureShearBCVx;
  setup.SetBCs->SetBCVz = SetPureShearBCVz;

  const double Eii         = 0.1;
  const int    Nt          = 5;
  const double dt          = 0.5;
  const double ani_fac_max = 1.0e6;

  setup.MutateInput  = mutateAniFstrainInput;
  // shear_style = 0 (pure), periodic_x = 0, pure_shear_ALE = 1
  g_aniFstrainOverride = { Eii, Nt, dt, ani_fac_max, 0, 0, 1,
                           "AniFstrainPureShear" };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  // Threshold note: pure-shear F is integrated through MDOODZ's Stokes solve
  // + particle advection + DeformationGradient chain (in contrast to the
  // simple-shear test, where F = [[1, γ̇·t], [0, 1]] is exact-by-construction
  // because Vx is linear in z under periodic BCs and the Stokes residual is
  // ~1e-9). Pure shear with ALE accumulates O(Eii·dt) per-step integration
  // error (forward-Euler-class), so residuals grow ~linearly with step
  // count. Observed at Eii = 0.1, dt = 0.5: ~6e-4 (step 2) → ~3.5e-3
  // (step 5). Thresholds are set with 2× safety margin above observed values.
  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainPureShear/Output%05d.gzip.h5", step);
    const double t   = (step - 1) * dt;
    const double ana = anaFsArPureShear(Eii, t);
    AniFacStats s    = readAniFacStats(hdf5Path, ana);

    printf("PureShear step=%d t=%.3f δ_ana=%.6f δ_mean=%.6f L2=%.4e "
           "min=%.4f max=%.4f\n",
           step, t, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-2)
        << "Pure shear step " << step
        << ": mean δ vs exp(2·Eii·(step-1)·dt)";
    EXPECT_LT(s.l2_rel, 1e-2)
        << "Pure shear step " << step << ": spatial L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "Pure shear step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainHansenOlivine — Hansen+12 / +16 olivine-calibrated δ-form
//   δ(γ) = 24.5 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.536, γ_e = 3.96
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//
// Constants come from a least-squares fit to the combined 38-point Hansen-
// group olivine torsion dataset (Hansen+14 EPSL 387,157 Table 1 + Hansen+16
// EPSL 445,92 Part 1 Fig.6 + Table 1). See
// misc/aniso_fstrain/calibrate/calibrate_olivine_hansen.py.
//
// Asymptote validation: 24.5·M_inf+1 = 14.13 agrees with Hansen+16 Part 1
// Eq.3-4 stress-aware direct calc (1.39/0.73)^4.1 = 14.02 to 100.8%.
// ---------------------------------------------------------------------------
namespace {
  // Mirror MDLIB's calibrated constants. Kept private to this test file —
  // any drift between the values here and the values in
  // MDLIB/AnisotropyRoutines.c will surface as a 1e-6 assertion failure.
  static constexpr double kHansenMInf       = 0.536;
  static constexpr double kHansenGammaE     = 3.96;
  static constexpr double kHansenDeltaSlope = 24.5;

  static double anaDeltaHansenOlivine(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kHansenMInf * (1.0 - exp(-g_eff / kHansenGammaE));
    return kHansenDeltaSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainHansenOlivine) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  // Span the full Hansen olivine strain range — γ ∈ {0, 1, 2, ..., 20}
  // covers initial growth, the lab-data cluster around γ ≈ 5-10, and the
  // saturating tail. Simple-shear F is exact-by-construction under periodic
  // BCs, so dt = 0.5 stays FP-precise at any γ.
  const double bkg_sr      = 1.0;     // → γ̇ = 2·bkg_sr = 2
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 21;
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;   // uncapped relative to δ_∞ ≈ 14

  // Override per-phase ani_fstrain to 2 (Hansen-olivine form). The base
  // fixture uses ani_fstrain = 1 by default; mutateAniFstrainInput does not
  // touch ani_fstrain, so we use a wrapper that calls the existing override
  // and then patches ani_fstrain directly. Non-capturing lambda decays to a
  // function pointer compatible with setup.MutateInput.
  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainHansenOlivine" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  // Output<k> reflects F at end of step k-1 → reference at γ = γ̇·(k-1)·dt.
  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainHansenOlivine/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaHansenOlivine(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("HansenOlivine step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Hansen-olivine step " << step
        << ": mean δ vs 24.5·M(γ_eff)+1 (M_inf=0.536, γ_e=3.96)";
    EXPECT_LT(s.l2_rel, 1e-4)
        << "Hansen-olivine step " << step << ": spatial L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "Hansen-olivine step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainHansenOlivineDamped — damped olivine δ-form via aniso_db = 7
//   δ(γ) = 24.5 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.16, γ_e = 0.76, slope = 24.5
//   asymptote δ_∞ = 4.84
//
// Constants come from a free-LSQ fit to the combined 93-sample non-Hansen
// olivine dataset (Tasaka 2016 + Boneh & Skemer 2014 + Kumamoto+19 +
// Bernard+19). Same δ-vs-M slope (24.5) as case 1 because Hansen+12 Fig 3b's
// linear relation is universal across olivine fabric strength. The lower
// M_inf reflects damping by water content, pre-existing CPO, or complex
// natural strain history.
//
// 2026-05-11 REFIT: Bernard+19's M re-derived from J via natural-sample
// formula M=(J−1)/24 (R²=0.825, Bernard+19 natural-sample regression), superseding
// Bernard's as-reported M (implicit Skemer-(J−1)/15 bias for natural
// xenoliths). Free LSQ on refit: (0.1567, 0.7635, RMS 0.0687); committed
// to 2 dp (0.16, 0.76). Previous committed (0.20, 1.30) is in git history.
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kOlivineDampedMInf       = 0.16;
  static constexpr double kOlivineDampedGammaE     = 0.76;
  static constexpr double kOlivineDampedDeltaSlope = 24.5;

  static double anaDeltaOlivineDamped(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kOlivineDampedMInf
                         * (1.0 - exp(-g_eff / kOlivineDampedGammaE));
    return kOlivineDampedDeltaSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainHansenOlivineDamped) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  // γ ∈ {0, 1, ..., 8} — damped olivine saturates at γ ≈ 4 (γ_e=1.30 means
  // by γ=5 we're at 99% of asymptote δ_∞ = 5.90).
  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 9;
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;   // uncapped relative to δ_∞ = 5.90

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainHansenOlivineDamped" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 7;   // damped olivine
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainHansenOlivineDamped/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaOlivineDamped(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("OlivineDamped step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "OlivineDamped step " << step
        << ": mean δ vs 24.5·M(γ_eff)+1 (M_inf=0.20, γ_e=1.30)";
    EXPECT_LT(s.l2_rel, 1e-4)
        << "OlivineDamped step " << step << ": spatial L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "OlivineDamped step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainCalcite — Pieri+01 / Barnhoorn+04 calibrated δ-form via aniso_db = 2
//   δ(γ) = 6.0 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.41, γ_e = 3.15, slope = 6.0
//   asymptote δ_∞ = 3.46
//
// Constants come from the combined Pieri+01 (5 pts) + Barnhoorn+04
// J ≤ 15 (10 pts) J-index torsion dataset (J→M-converted via Skemer+05).
// 15-pt fit, post-audit (Bruijn+11 retracted 2026-05-08).  The slope
// value of 6.0 is bracketed by Tommasi+09-style VPSC bounds for calcite
// saturated CPO.  See misc/aniso_fstrain/data/calcite/.
//
// The test runs to γ = 10 — calcite is well-saturated by then (δ at
// ~92% of asymptote 3.46).
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kCalciteMInf   = 0.41;
  static constexpr double kCalciteGammaE = 3.15;   // updated 2026-05-08 (Bruijn+11 retracted)
  static constexpr double kCalciteSlope  = 6.0;

  static double anaDeltaCalcite(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kCalciteMInf * (1.0 - exp(-g_eff / kCalciteGammaE));
    return kCalciteSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainCalcite) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;     // → γ̇ = 2·bkg_sr = 2
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 11;      // γ ∈ {0, 1, 2, ..., 10}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;   // uncapped relative to δ_∞ = 4

  // Override per-phase ani_fstrain to 2 (calibrated form) AND aniso_db to 2
  // (calcite). The auto-default routes ani_fstrain=2 + aniso_db=0 to olivine,
  // so calcite users MUST set aniso_db = 2 explicitly.
  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainCalcite" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 2;   // calcite (Pieri+01 + Barnhoorn+04, post-audit)
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainCalcite/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaCalcite(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("Calcite step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Calcite step " << step
        << ": mean δ vs 6.0·M(γ_eff)+1 (M_inf=0.41, γ_e=3.15)";
    EXPECT_LT(s.l2_rel, 1e-4)
        << "Calcite step " << step << ": spatial L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "Calcite step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainQuartz — Pennacchioni+10 calibrated δ-form via aniso_db = 3
//   δ(γ) = 3.0 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.75, γ_e = 3.50, slope = 3.0
//   asymptote δ_∞ = 3.25
//
// Constants come from a free-LSQ fit to the Pennacchioni+10 17-point
// natural-shear-zone dataset (ODF J-index, J→M-converted via Skemer+05).
// The slope of 3.0 is fixed by VPSC-style bounds for quartz aggregate
// viscous anisotropy at saturated CPO. T-regime: ~500°C, prism-a slip.
//
// Test runs to γ = 8 (Nt = 9, dt = 0.5) — by γ = 8 we're at 90% of
// asymptote (δ ≈ 3.03 of δ_∞ = 3.25).
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kQuartzMInf   = 0.75;
  static constexpr double kQuartzGammaE = 3.50;
  static constexpr double kQuartzSlope  = 3.0;

  static double anaDeltaQuartz(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kQuartzMInf * (1.0 - exp(-g_eff / kQuartzGammaE));
    return kQuartzSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainQuartz) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;     // → γ̇ = 2·bkg_sr = 2
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 9;       // γ ∈ {0, 1, 2, ..., 8}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;   // uncapped relative to δ_∞ = 3.25

  // Override per-phase ani_fstrain to 2 + aniso_db to 3 (quartz). Auto-default
  // routes ani_fstrain=2 + aniso_db=0 to olivine, so quartz users MUST set
  // aniso_db = 3 explicitly.
  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainQuartz" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 3;   // quartz (Pennacchioni+10)
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainQuartz/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaQuartz(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("Quartz step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Quartz step " << step
        << ": mean δ vs 3.0·M(γ_eff)+1 (M_inf=0.75, γ_e=3.50)";
    EXPECT_LT(s.l2_rel, 1e-4)
        << "Quartz step " << step << ": spatial L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "Quartz step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainQuartzBlackford — Blackford+24 calibrated δ-form via
// aniso_db = 6
//   δ(γ) = 11.0 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.20, γ_e = 4.00, slope = 11.0
//   asymptote δ_∞ = 3.20  (matches case 3's δ_∞ = 3.25 to within rounding)
//
// Constants come from a centroid LSQ fit to the Blackford+24 38-point
// Northern Snake Range dataset (pole-figure J of c-axis, J→M-converted
// via Skemer+05; γ_oct→γ_eq via 2·sinh(ε/√2) plane-strain inversion).
// Slope is rescaled to 11.0 to compensate for pfJ < ODF J for the same
// physical texture, yielding a comparable physical δ_∞ to case 3.
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kQuartzBlaMInf   = 0.20;
  static constexpr double kQuartzBlaGammaE = 4.00;
  static constexpr double kQuartzBlaSlope  = 11.0;

  static double anaDeltaQuartzBlackford(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kQuartzBlaMInf * (1.0 - exp(-g_eff / kQuartzBlaGammaE));
    return kQuartzBlaSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainQuartzBlackford) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 9;       // γ ∈ {0, 1, ..., 8}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;   // uncapped relative to δ_∞ = 3.20

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainQuartzBlackford" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 6;   // quartz (Blackford+24)
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainQuartzBlackford/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaQuartzBlackford(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("QuartzBlackford step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "QuartzBlackford step " << step
        << ": mean δ vs 11.0·M(γ_eff)+1 (M_inf=0.20, γ_e=4.00)";
    EXPECT_LT(s.l2_rel, 1e-4)
        << "QuartzBlackford step " << step << ": spatial L2 error";
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4)
        << "QuartzBlackford step " << step << ": spatial homogeneity";
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainCalciteLowT — Barnhoorn+04 mid-low-T regime via aniso_db = 4
//   M_inf = 0.22, γ_e = 1.37, slope = 6.0, asymptote δ_∞ = 2.33
//
// Subgrain-rotation recrystallization dominates at T ≤ 650°C; CPO is
// randomized at high γ, capping fabric at modest strength. Test runs to
// γ = 6 (well past the γ_e=1.37 saturation strain).
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kCalciteLowTMInf   = 0.22;
  static constexpr double kCalciteLowTGammaE = 1.37;
  static constexpr double kCalciteLowTSlope  = 6.0;

  static double anaDeltaCalciteLowT(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kCalciteLowTMInf * (1.0 - exp(-g_eff / kCalciteLowTGammaE));
    return kCalciteLowTSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainCalciteLowT) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 7;       // γ ∈ {0, 1, 2, ..., 6}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainCalciteLowT" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 4;   // calcite mid-low-T (Barnhoorn+04)
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainCalciteLowT/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaCalciteLowT(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("CalciteLowT step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Calcite mid-low-T step " << step
        << ": mean δ vs 6.0·M(γ_eff)+1 (M_inf=0.22, γ_e=1.37)";
    EXPECT_LT(s.l2_rel, 1e-4);
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4);
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainCalciteHighT — Pieri+01 + Barnhoorn+04(727°C, J ≤ 15)
//   high-T regime via aniso_db = 5
//   M_inf = 0.73, γ_e = 5.28, slope = 6.0, asymptote δ_∞ = 5.40
//
// Grain-boundary-migration recrystallization dominates at T > 650°C; CPO
// continues to strengthen at high γ.  Saturation is faster (γ_e = 5.28)
// than the previous 8-pt + Bruijn fit (γ_e = 8.24) — see
// MDLIB/FlowLaws.c case 5 docstring for the audit-driven update
// (Bruijn+11 retracted 2026-05-08).  Test runs to γ = 14 — by then
// δ is at ~96% of asymptote 5.40.
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kCalciteHighTMInf   = 0.73;   // updated 2026-05-08
  static constexpr double kCalciteHighTGammaE = 5.28;   // updated 2026-05-08
  static constexpr double kCalciteHighTSlope  = 6.0;

  static double anaDeltaCalciteHighT(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kCalciteHighTMInf * (1.0 - exp(-g_eff / kCalciteHighTGammaE));
    return kCalciteHighTSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainCalciteHighT) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 15;      // γ ∈ {0, 1, 2, ..., 14}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainCalciteHighT" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 5;   // calcite high-T
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainCalciteHighT/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaCalciteHighT(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("CalciteHighT step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Calcite high-T step " << step
        << ": mean δ vs 6.0·M(γ_eff)+1 (M_inf=0.73, γ_e=5.28)";
    EXPECT_LT(s.l2_rel, 1e-4);
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4);
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainQuartzMicaBracketed — case 8 polyphase bracketed estimate
//   δ(γ) = 15.0 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.40, γ_e = 2.50, slope = 15.0
//   asymptote δ_∞ = 7.0
//
// IMPORTANT: this is a BRACKETED estimate, NOT a direct LSQ fit.  The CI
// test gates the IMPLEMENTATION (MDOODZ matches the formula to 1e-6) —
// it does NOT validate the science.  Bracketing logic is documented in
// misc/aniso_fstrain/notes/quartz_calibration.md §11.
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kQzMicaMInf       = 0.40;
  static constexpr double kQzMicaGammaE     = 2.50;
  static constexpr double kQzMicaDeltaSlope = 15.0;

  static double anaDeltaQuartzMicaBracketed(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kQzMicaMInf
                         * (1.0 - exp(-g_eff / kQzMicaGammaE));
    return kQzMicaDeltaSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainQuartzMicaBracketed) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 9;       // γ ∈ {0, 1, ..., 8}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainQuartzMicaBracketed" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 8;   // bracketed quartz-mica polyphase
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainQuartzMicaBracketed/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaQuartzMicaBracketed(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("QzMicaBracketed step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "QzMicaBracketed step " << step
        << ": mean δ vs 15.0·M(γ_eff)+1 (M_inf=0.40, γ_e=2.50)";
    EXPECT_LT(s.l2_rel, 1e-4);
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4);
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainPlagioclaseBracketed — case 9 polyphase bracketed estimate
//   δ(γ) = 8.0 · M_inf · (1 − exp(−γ_eff / γ_e)) + 1
//   M_inf = 0.12, γ_e = 2.00, slope = 8.0
//   asymptote δ_∞ = 1.96 ≈ 2.0
//
// IMPORTANT: BRACKETED estimate, NOT direct LSQ.  See
// misc/aniso_fstrain/notes/quartz_calibration.md §12 for the
// literature bracket details.  Saturating-exponential form is itself a
// simplification — natural plagioclase CPO can WEAKEN at very high γ
// due to DisGBS / diffusion creep.
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kPlagMInf       = 0.12;
  static constexpr double kPlagGammaE     = 2.00;
  static constexpr double kPlagDeltaSlope = 8.0;

  static double anaDeltaPlagioclaseBracketed(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kPlagMInf
                         * (1.0 - exp(-g_eff / kPlagGammaE));
    return kPlagDeltaSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainPlagioclaseBracketed) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 9;       // γ ∈ {0, 1, ..., 8}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainPlagioclaseBracketed" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 9;   // bracketed plagioclase polyphase
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainPlagioclaseBracketed/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaPlagioclaseBracketed(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("PlagBracketed step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "PlagBracketed step " << step
        << ": mean δ vs 8.0·M(γ_eff)+1 (M_inf=0.12, γ_e=2.00)";
    EXPECT_LT(s.l2_rel, 1e-4);
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4);
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainCalciteSuperplastic — case 13 calcite regime-3 bracketed
//   M_inf = 0.08, γ_e = 1.00, slope = 1.5, asymptote δ_∞ = 1.12
//
// Solnhofen-limestone-regime-3 analog: fine-grained calcite mylonite in
// GBS-dominant low-σ high-T regime where CPO is actively destroyed.
// σ_c(T) is LSQ-fitted (R² = 0.86, n = 4 Arrhenius); M_∞ / γ_e / slope are
// literature-bracketed from Schmid 1976, Schmid+77, De Bresser 2002 (same
// discipline as cases 8 / 9). See misc/aniso_fstrain/notes/case13_calibration.md.
//
// Very weak anisotropy (δ_∞=1.12, only +12% over isotropic) reflects the
// regime-3 GBS-dominated rheology where CPO is destroyed faster than it
// can build up. Test runs to γ = 6 (γ_e=1.0 → saturation at γ ≈ 5).
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kCalciteSuperplasticMInf   = 0.08;
  static constexpr double kCalciteSuperplasticGammaE = 1.00;
  static constexpr double kCalciteSuperplasticSlope  = 1.5;

  static double anaDeltaCalciteSuperplastic(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kCalciteSuperplasticMInf * (1.0 - exp(-g_eff / kCalciteSuperplasticGammaE));
    return kCalciteSuperplasticSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainCalciteSuperplastic) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 7;       // γ ∈ {0, 1, 2, ..., 6}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainCalciteSuperplastic" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 13;  // calcite superplastic (Schmid 1976/77 + De Bresser 2002, Tier-2 HYBRID)
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainCalciteSuperplastic/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaCalciteSuperplastic(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("CalciteSuperplastic step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "Calcite superplastic step " << step
        << ": mean δ vs 1.5·M(γ_eff)+1 (M_inf=0.08, γ_e=1.00)";
    EXPECT_LT(s.l2_rel, 1e-4);
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4);
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainQuartzMicaPolyphase — case 15 polyphase bracketed estimate
//   δ(γ) = 15.0 · M_init · exp(−γ_eff / γ_decay) + 1
//   M_init = 0.359, γ_decay = 15.35, slope = 15.0
//   δ(γ → 0) = 6.39, δ(γ → ∞) = 1.0
//
// Companion to case 8: case 8 = mica CPO STRENGTHENING with γ
// (saturating-exp), case 15 = quartz CPO WEAKENING with γ (decay form
// per Tokle+23 §3.2.3 "shielding" mechanism).
//
// IMPORTANT: this is a BRACKETED estimate, NOT a Tier-1 direct LSQ.
// The CI test gates the IMPLEMENTATION (MDOODZ matches the formula to
// 1e-6) — it does NOT validate the science.  Pooled-decay LSQ details
// + m.u.d. → M-index conversion choice documented in
// misc/aniso_fstrain/notes/case15_calibration.md.
// ---------------------------------------------------------------------------
namespace {
  static constexpr double kQzMicaPolyMInit      = 0.359;
  static constexpr double kQzMicaPolyGammaDecay = 15.35;
  static constexpr double kQzMicaPolyDeltaSlope = 15.0;

  static double anaDeltaQuartzMicaPolyphase(double gamma_dot, double t) {
    const double g     = gamma_dot * t;
    const double fs_ar = 1.0 + 0.5 * g * g + g * sqrt(0.25 * g * g + 1.0);
    const double s     = sqrt(fs_ar);
    const double g_eff = s - 1.0 / s;
    const double M     = kQzMicaPolyMInit
                         * exp(-g_eff / kQzMicaPolyGammaDecay);
    return kQzMicaPolyDeltaSlope * M + 1.0;
  }
}

TEST_F(AnisotropyBenchmark, AniFstrainQuartzMicaPolyphase) {
  setup.SetBCs->SetBCVx = SetSimpleShearBCVx;
  setup.SetBCs->SetBCVz = SetSimpleShearBCVz;

  const double bkg_sr      = 1.0;
  const double gamma_dot   = 2.0 * bkg_sr;
  const int    Nt          = 9;       // γ ∈ {0, 1, ..., 8}
  const double dt          = 0.5;
  const double ani_fac_max = 100.0;

  g_aniFstrainOverride = { bkg_sr, Nt, dt, ani_fac_max, 1, 1, 1,
                           "AniFstrainQuartzMicaPolyphase" };
  setup.MutateInput = [](MdoodzInput *input) {
    mutateAniFstrainInput(input);
    input->materials.ani_fstrain[0] = 2;
    input->materials.aniso_db[0]    = 15;  // bracketed quartz-mica polyphase, quartz perspective
  };

  char *inputName;
  asprintf(&inputName, "AnisotropyBenchmark/AniFstrainEvolution.txt");
  RunMDOODZ(inputName, &setup);
  free(inputName);

  setup.MutateInput  = nullptr;
  g_aniFstrainOverride = {0.0, 0, 0.0, 0.0, -1, -1, -99, nullptr};

  for (int step = 1; step <= Nt; step++) {
    char hdf5Path[256];
    snprintf(hdf5Path, sizeof(hdf5Path),
             "AniFstrainQuartzMicaPolyphase/Output%05d.gzip.h5", step);
    const double t     = (step - 1) * dt;
    const double gamma = gamma_dot * t;
    const double ana   = anaDeltaQuartzMicaPolyphase(gamma_dot, t);
    AniFacStats  s     = readAniFacStats(hdf5Path, ana);

    printf("QzMicaPolyphase step=%d t=%.3f γ=%.3f δ_ana=%.6f δ_mean=%.6f "
           "L2=%.4e min=%.4f max=%.4f\n",
           step, t, gamma, ana, s.mean, s.l2_rel, s.min_v, s.max_v);

    EXPECT_NEAR(s.mean / ana, 1.0, 1e-6)
        << "QzMicaPolyphase step " << step
        << ": mean δ vs 15.0·M(γ_eff)+1 (M_init=0.359, γ_decay=15.35, decay form)";
    EXPECT_LT(s.l2_rel, 1e-4);
    EXPECT_LT(s.max_v / s.min_v - 1.0, 1e-4);
    EXPECT_GE(s.min_v, 0.999);
    EXPECT_TRUE(std::isfinite(s.min_v));
    EXPECT_TRUE(std::isfinite(s.max_v));
  }
}

// ===========================================================================
// ani_fstrain == 3 — δ-relaxation analytical unit tests (capability:
// aniso-delta-relaxation). The relaxation is verified against its closed-form
// limits by exercising the MDLIB routine DeltaRelaxationTau directly together
// with the analytic exponential update — mirroring how anaDeltaHansenOlivine
// above mirrors MDLIB's calibrated δ-form. This is a 1-marker / homogeneous
// "box" in the sense the spec intends: a deterministic per-marker calculation
// with no spatial coupling. Running the math directly (rather than through a
// full sim that cannot easily pin per-marker T / strain_pwl) keeps the
// closed-form-limit assertions exact and gives the Δt-convergence and
// Arrhenius-ratio sub-tests the control they need.
//
// All tests use an identity scale struct (every scale = 1) so scaled units ==
// SI units; with that, materials->R == Rg. The gs flow law is OFF throughout
// (L_relax is supplied directly), which covers the spec's "runs with
// grain-size evolution disabled" scenario.
// ===========================================================================

extern "C" {
  // MDLIB/AnisotropyRoutines.c — relaxation timescale τ_relax = L_relax / V,
  // V = M0·exp(-Q/RT)·μ·b²·Δρ (Boneh et al. 2021). Declared here because it
  // lives in mdoodz-private.h; the symbol has external linkage in libmdoodz.
  // The Boneh kinetics constants are now per-phase material-database fields
  // (mat_prop::ani_relax_*, set in ReadDataAnisotropy / FlowLaws.c), so the
  // routine takes them as explicit arguments — see DeltaRelaxationTau() below.
  double DeltaRelaxationTau(double T_scaled, double L_relax_scaled,
                            double strain_pwl, double R_scaled, scale scaling,
                            double Q, double M0, double mu, double b,
                            double drho_min, double drho_max, double eps_ref);
}

namespace {
  // Identity scale — scaled units == SI units.
  static scale identityScale() {
    scale s;
    s.eta = 1.0; s.L = 1.0; s.V = 1.0; s.T = 1.0; s.t = 1.0; s.a = 1.0;
    s.E = 1.0; s.S = 1.0; s.m = 1.0; s.rho = 1.0; s.F = 1.0; s.J = 1.0;
    s.W = 1.0; s.Cp = 1.0; s.rhoE = 1.0; s.k = 1.0;
    return s;
  }

  // The analytic exponential update used inside FiniteStrainAspectRatio:
  //   δ_new = 1 + (δ_prod - 1)·exp(-dt/τ), clamped to [1, ani_fac_max].
  static double relaxStep(double delta_prod, double dt, double tau,
                          double ani_fac_max) {
    double d = 1.0 + (delta_prod - 1.0) * exp(-dt / tau);
    if (d > ani_fac_max) d = ani_fac_max;
    if (d < 1.0)         d = 1.0;
    return d;
  }

  // Boneh constants mirrored from MDLIB/FlowLaws.c (ReadDataAnisotropy, the
  // per-phase olivine calibration of mat_prop::ani_relax_*) — any drift
  // surfaces as an assertion failure here.
  static constexpr double kBonehQ_test       = 133.0e3;   // J/mol
  static constexpr double kBonehMu_test      = 50.0e9;    // Pa
  static constexpr double kBonehB_test       = 0.6e-9;    // m
  static constexpr double kBonehM0_test      = 2.0e-11;   // m^4 J^-1 s^-1
  static constexpr double kBonehDrhoMin_test = 1.0e11;    // m^-2
  static constexpr double kBonehDrhoMax_test = 1.0e13;    // m^-2
  static constexpr double kBonehEpsRef_test  = 4.0;       // [-]

  // Thin wrapper passing the olivine Boneh constants — keeps the relaxation
  // unit tests reading the same way they did when the constants were module
  // globals in AnisotropyRoutines.c. Mirrors what FiniteStrainAspectRatio now
  // does with materials->ani_relax_*[phase].
  static double DeltaRelaxationTauOlivine(double T_scaled, double L_relax_scaled,
                                          double strain_pwl, double R_scaled,
                                          scale scaling) {
    return DeltaRelaxationTau(T_scaled, L_relax_scaled, strain_pwl, R_scaled,
                              scaling, kBonehQ_test, kBonehM0_test,
                              kBonehMu_test, kBonehB_test, kBonehDrhoMin_test,
                              kBonehDrhoMax_test, kBonehEpsRef_test);
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxColdFreeze — a cold marker has τ_relax ≫ run time, so
// δ is constant to floating-point precision over the run. (Spec scenario:
// "Cold marker is effectively frozen".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxColdFreeze) {
  scale sc = identityScale();
  // Cold: 300 K. L_relax = 1 mm, modest stored strain.
  const double T          = 300.0;
  const double L_relax    = 1.0e-3;
  const double strain_pwl = 1.0;
  const double tau = DeltaRelaxationTauOlivine(T, L_relax, strain_pwl, Rg, sc);

  // τ_relax must dwarf any plausible run time (here: 1e6 s ≈ 11.6 days).
  const double run_time = 1.0e6;
  printf("RelaxColdFreeze T=%.1fK tau=%.6e s  run_time=%.1e s\n",
         T, tau, run_time);
  EXPECT_GT(tau, 1.0e6 * run_time)
      << "cold marker: τ_relax should be astronomically larger than run time";

  // Integrate with no production over many steps — δ must not move.
  const double delta0 = 5.0;
  double delta = delta0;
  const double dt = run_time / 100.0;
  for (int i = 0; i < 100; i++) delta = relaxStep(delta, dt, tau, 1.0e6);
  EXPECT_NEAR(delta, delta0, 1e-12)
      << "cold marker: δ must be constant to FP precision";
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxHot — a hot marker with no production relaxes as
//   δ(t) = 1 + (δ₀ - 1)·exp(-t/τ_relax)
// on the Boneh timescale. Includes the Arrhenius-ratio sub-assertion:
//   τ_relax(T₂)/τ_relax(T₁) = exp(Q/R·(1/T₁ - 1/T₂)).
// (Spec scenarios: "Hot marker relaxes toward isotropic" + "Relaxation rate
// follows the Arrhenius law".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxHot) {
  scale sc = identityScale();
  const double T1         = 1473.15;   // ~1200 °C
  const double T2         = 1673.15;   // ~1400 °C — hotter
  const double L_relax    = 1.0e-3;
  const double strain_pwl = 5.0;

  const double tau1 = DeltaRelaxationTauOlivine(T1, L_relax, strain_pwl, Rg, sc);
  const double tau2 = DeltaRelaxationTauOlivine(T2, L_relax, strain_pwl, Rg, sc);
  printf("RelaxHot tau(%.1fK)=%.6e s  tau(%.1fK)=%.6e s\n", T1, tau1, T2, tau2);
  EXPECT_TRUE(std::isfinite(tau1) && tau1 > 0.0);
  EXPECT_TRUE(std::isfinite(tau2) && tau2 > 0.0);

  // Hotter ⇒ faster ⇒ shorter τ.
  EXPECT_LT(tau2, tau1) << "Arrhenius: hotter marker must relax faster";

  // Arrhenius ratio: only the exp(-Q/RT) term depends on T, so
  // τ(T1)/τ(T2) = exp(Q/R·(1/T2 - 1/T1))  (note: τ ∝ 1/V ∝ exp(+Q/RT)).
  const double ratio_expected = exp(kBonehQ_test / Rg * (1.0 / T1 - 1.0 / T2));
  EXPECT_NEAR((tau1 / tau2) / ratio_expected, 1.0, 1e-9)
      << "τ ratio must equal the pure Arrhenius factor";

  // δ(t) decay matches the analytic exponential, monotone, never below 1.
  const double delta0 = 4.0;
  double delta = delta0;
  const double dt = tau2 / 20.0;     // fine steps relative to the fast τ
  double prev = delta0;
  for (int i = 1; i <= 200; i++) {
    delta = relaxStep(delta, dt, tau2, 1.0e6);
    const double t   = i * dt;
    const double ana = 1.0 + (delta0 - 1.0) * exp(-t / tau2);
    EXPECT_NEAR(delta / ana, 1.0, 1e-9)
        << "hot relaxation step " << i << ": δ(t) vs 1+(δ₀-1)exp(-t/τ)";
    EXPECT_LE(delta, prev + 1e-12) << "δ must decay monotonically";
    EXPECT_GE(delta, 1.0)          << "δ must not cross below 1";
    prev = delta;
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxHalfLife — under constant T, Δρ, L_relax the no-
// production update reaches δ = 1 + (δ₀-1)/2 at t = τ_relax·ln2. Includes a
// Δt-convergence sub-assertion (the analytic update is exact for the linear
// ODE → identical to FP regardless of Δt) and a large-Δt no-overshoot
// sub-assertion (Δt ≫ τ ⇒ 1 ≤ δ_new ≤ δ_old).
// (Spec scenarios: "Exponential half-life" + "No overshoot at large timestep".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxHalfLife) {
  scale sc = identityScale();
  const double T          = 1573.15;
  const double L_relax    = 1.0e-3;
  const double strain_pwl = 5.0;
  const double tau = DeltaRelaxationTauOlivine(T, L_relax, strain_pwl, Rg, sc);
  ASSERT_TRUE(std::isfinite(tau) && tau > 0.0);

  const double delta0 = 6.0;
  const double half   = 1.0 + (delta0 - 1.0) / 2.0;

  // Single exact step of length τ·ln2 lands on the half-life value.
  const double d_half = relaxStep(delta0, tau * M_LN2, tau, 1.0e6);
  printf("RelaxHalfLife tau=%.6e s  δ₀=%.3f  δ(τln2)=%.6f  expected=%.6f\n",
         tau, delta0, d_half, half);
  EXPECT_NEAR(d_half, half, 1e-9) << "δ(τ·ln2) must equal 1+(δ₀-1)/2";

  // Δt-convergence: reaching t = τ·ln2 in N sub-steps must give the same
  // answer for any N (the analytic exponential update telescopes exactly).
  for (int N : {1, 2, 4, 10, 100, 1000}) {
    double d = delta0;
    const double dt = tau * M_LN2 / N;
    for (int i = 0; i < N; i++) d = relaxStep(d, dt, tau, 1.0e6);
    EXPECT_NEAR(d, half, 1e-9)
        << "Δt-convergence: N=" << N << " sub-steps must reach the same δ";
  }

  // Large-Δt no-overshoot: Δt ≫ τ ⇒ 1 ≤ δ_new ≤ δ_old (cannot undershoot 1).
  for (double mult : {10.0, 100.0, 1.0e6}) {
    const double d = relaxStep(delta0, mult * tau, tau, 1.0e6);
    EXPECT_GE(d, 1.0)    << "large Δt: δ_new must not drop below 1";
    EXPECT_LE(d, delta0) << "large Δt: δ_new must not exceed δ_old";
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxColdLimitEqualsFstrain2 — at cold T (τ_relax ≫ run
// time) the ani_fstrain == 3 two-field operator split telescopes to exactly
// aniso_delta_fn(FS_AR) — i.e. it tracks ani_fstrain == 2 to FP precision.
// (Spec scenario: "Cold marker under ani_fstrain == 3 tracks the
// ani_fstrain == 2 form".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxColdLimitEqualsFstrain2) {
  scale sc = identityScale();
  const double T          = 300.0;     // cold → τ ≫ run time
  const double L_relax    = 1.0e-3;
  // strain_pwl grows with the simulated shear; τ stays astronomically large.
  const double dt         = 0.5;
  const double gamma_dot  = 2.0;

  // Two-field operator split, replicated from FiniteStrainAspectRatio:
  //   delta_FS   = aniso_delta_fn(FS_AR)
  //   delta_prod = aniso_delta + (delta_FS - aniso_delta_fs_prev)
  //   aniso_delta = relaxStep(delta_prod, dt, tau, ani_fac_max)
  double aniso_delta = 1.0, aniso_delta_fs_prev = 1.0;
  for (int step = 1; step <= 20; step++) {
    const double t        = step * dt;
    const double delta_FS = anaDeltaHansenOlivine(gamma_dot, t);  // ani_fstrain==2 value
    const double strain_pwl = gamma_dot * t;   // monotone proxy input
    const double tau = DeltaRelaxationTauOlivine(T, L_relax, strain_pwl, Rg, sc);
    const double delta_prod = aniso_delta + (delta_FS - aniso_delta_fs_prev);
    aniso_delta = relaxStep(delta_prod, dt, tau, 1.0e6);
    aniso_delta_fs_prev = delta_FS;
    // Cold limit: aniso_delta must equal the pure ani_fstrain==2 form.
    EXPECT_NEAR(aniso_delta / delta_FS, 1.0, 1e-9)
        << "cold-limit step " << step
        << ": ani_fstrain==3 δ must track aniso_delta_fn(FS_AR)";
  }
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxLengthMonotonic — τ_relax ∝ L_relax, so a larger
// L_relax relaxes more slowly. (Spec scenario: "L_relax monotonicity".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxLengthMonotonic) {
  scale sc = identityScale();
  const double T          = 1573.15;
  const double strain_pwl = 5.0;
  const double L1 = 1.0e-4;   // 0.1 mm
  const double L2 = 1.0e-3;   // 1 mm  (larger)
  const double L3 = 1.0e-2;   // 1 cm  (largest)
  const double tau1 = DeltaRelaxationTauOlivine(T, L1, strain_pwl, Rg, sc);
  const double tau2 = DeltaRelaxationTauOlivine(T, L2, strain_pwl, Rg, sc);
  const double tau3 = DeltaRelaxationTauOlivine(T, L3, strain_pwl, Rg, sc);
  printf("RelaxLengthMonotonic tau(L=%.1e)=%.4e  tau(L=%.1e)=%.4e  tau(L=%.1e)=%.4e\n",
         L1, tau1, L2, tau2, L3, tau3);
  EXPECT_LT(tau1, tau2) << "larger L_relax ⇒ slower relaxation (longer τ)";
  EXPECT_LT(tau2, tau3) << "larger L_relax ⇒ slower relaxation (longer τ)";
  // τ ∝ L_relax exactly (V is L-independent).
  EXPECT_NEAR((tau2 / tau1) / (L2 / L1), 1.0, 1e-9)
      << "τ_relax must be linear in L_relax";
  EXPECT_NEAR((tau3 / tau1) / (L3 / L1), 1.0, 1e-9)
      << "τ_relax must be linear in L_relax";

  // A larger L_relax ⇒ less decay after the same elapsed time.
  const double delta0 = 5.0;
  const double t = tau2;   // one e-fold for L2
  const double d1 = 1.0 + (delta0 - 1.0) * exp(-t / tau1);
  const double d2 = 1.0 + (delta0 - 1.0) * exp(-t / tau2);
  const double d3 = 1.0 + (delta0 - 1.0) * exp(-t / tau3);
  EXPECT_LT(d1, d2) << "after fixed t, larger L_relax retains more anisotropy";
  EXPECT_LT(d2, d3) << "after fixed t, larger L_relax retains more anisotropy";
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxClampAndFixedPoint — the ani_fac_max ceiling is
// honoured and δ = 1 is an exact fixed point of the relaxation update.
// (Spec: "ani_fac_max clamp honoured" + "isotropic fixed point".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxClampAndFixedPoint) {
  scale sc = identityScale();
  const double T          = 1573.15;
  const double L_relax    = 1.0e-3;
  const double strain_pwl = 5.0;
  const double tau = DeltaRelaxationTauOlivine(T, L_relax, strain_pwl, Rg, sc);
  ASSERT_TRUE(std::isfinite(tau) && tau > 0.0);

  // Fixed point: δ = 1 stays at 1 for any dt (1 + 0·exp(...) = 1).
  for (double dt : {1.0e-3 * tau, tau, 1.0e3 * tau}) {
    const double d = relaxStep(1.0, dt, tau, 1.0e6);
    EXPECT_NEAR(d, 1.0, 1e-15) << "δ = 1 must be an exact fixed point";
  }

  // Clamp: a production δ above ani_fac_max is capped. With dt ≪ τ almost no
  // relaxation happens, so delta_prod ≈ huge ⇒ result must be ani_fac_max.
  const double ani_fac_max = 10.0;
  const double d_clamped = relaxStep(1.0e3, 1.0e-6 * tau, tau, ani_fac_max);
  EXPECT_NEAR(d_clamped, ani_fac_max, 1e-9)
      << "δ above the ceiling must be clamped to ani_fac_max";

  // And the clamp never lifts δ below the isotropic floor.
  const double d_floor = relaxStep(1.0e3, 1.0e9 * tau, tau, ani_fac_max);
  EXPECT_GE(d_floor, 1.0) << "clamp must not push δ below 1";
  EXPECT_LE(d_floor, ani_fac_max);
}

// ---------------------------------------------------------------------------
// Test: AniFstrainRelaxHistoryDependence — two markers that reach the SAME
// current FS_AR via different thermal histories (one hot, one cold) end with
// different stored aniso_delta. This is the defining property of δ being a
// per-marker integrated state variable rather than a memoryless function of
// FS_AR. (Spec scenario: "History-dependence".)
// ---------------------------------------------------------------------------
TEST_F(AnisotropyBenchmark, AniFstrainRelaxHistoryDependence) {
  scale sc = identityScale();
  const double L_relax   = 1.0e-3;
  const double gamma_inc = 1.0;        // accumulated-shear increment per step
  const int    Nsteps    = 20;
  const double T_hot     = 1673.15;    // relaxes appreciably each step
  const double T_cold    = 300.0;      // effectively frozen

  // The two markers see the IDENTICAL FS_AR / strain_pwl trajectory (driven by
  // the same accumulated shear γ) — only their temperature differs. The
  // physical timestep dt is chosen as a meaningful fraction of the hot
  // marker's τ_relax so the hot marker actually relaxes within the run; the
  // γ-increment per step is independent of dt (in a real run dt would be the
  // Courant step and γ would accumulate at the local strain rate — here we
  // decouple them so the test is not at the mercy of identity-scale unit
  // coincidences). The two-field operator split is replicated from
  // FiniteStrainAspectRatio.
  const double strain_pwl_ref = gamma_inc * Nsteps;   // representative Δρ input
  const double tau_hot_ref = DeltaRelaxationTauOlivine(T_hot, L_relax, strain_pwl_ref, Rg, sc);
  const double dt = tau_hot_ref / 4.0;                // ~one e-fold over the run

  double hot_delta = 1.0,  hot_prev = 1.0;
  double cold_delta = 1.0, cold_prev = 1.0;
  double delta_FS_final = 0.0;
  for (int step = 1; step <= Nsteps; step++) {
    const double gamma      = gamma_inc * step;
    // anaDeltaHansenOlivine(γ̇, t) depends only on the product γ̇·t = γ.
    const double delta_FS   = anaDeltaHansenOlivine(gamma, 1.0);
    const double strain_pwl = gamma;
    delta_FS_final = delta_FS;

    const double tau_hot  = DeltaRelaxationTauOlivine(T_hot,  L_relax, strain_pwl, Rg, sc);
    const double tau_cold = DeltaRelaxationTauOlivine(T_cold, L_relax, strain_pwl, Rg, sc);

    hot_delta  = relaxStep(hot_delta  + (delta_FS - hot_prev),  dt, tau_hot,  1.0e6);
    cold_delta = relaxStep(cold_delta + (delta_FS - cold_prev), dt, tau_cold, 1.0e6);
    hot_prev = delta_FS; cold_prev = delta_FS;
  }

  printf("RelaxHistoryDependence FS_AR-form δ=%.6f  hot δ=%.6f  cold δ=%.6f  dt=%.3e s\n",
         delta_FS_final, hot_delta, cold_delta, dt);
  // Cold marker tracks the FS_AR form; hot marker has relaxed away from it.
  EXPECT_NEAR(cold_delta / delta_FS_final, 1.0, 1e-6)
      << "cold-history marker tracks aniso_delta_fn(FS_AR)";
  EXPECT_LT(hot_delta, cold_delta - 1e-3)
      << "hot-history marker must have a lower stored δ — same FS_AR, "
         "different thermal history ⇒ different δ (history dependence)";
  EXPECT_GE(hot_delta, 1.0);
}
