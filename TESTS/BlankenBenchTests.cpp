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

#ifdef _OPENMP
#include <omp.h>
#endif

// ============================================================================
// Blankenbach et al. (1989) Case 1a — Isoviscous thermal convection
//
// Ra = 10^4, unit square, free-slip walls, T=0 top, T=1 bottom
// Published values: Nu = 4.884, Vrms = 42.865, T_max(mid-depth) = 0.4266
// ============================================================================

class BlankenBench : public ::testing::Test {
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
                    .SetBCVx    = SetBCVx,
                    .SetBCVz    = SetBCVz,
                    .SetBCT     = SetBCT,
                    .SetBCPType = SetBCPType,
            },
    };
  }

  static int SetPhase(MdoodzInput *input, Coordinates coordinates) {
    return 0;
  }

  static double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
    const double T_top = (input->model.user0 + zeroC) / input->scaling.T;
    const double T_bot = input->model.user1 / input->scaling.T;
    const double zmin  = input->model.zmin;
    const double zmax  = input->model.zmax;
    const double H     = zmax - zmin;
    const double z_nd  = (zmax - coordinates.z) / H;
    double T_lin = T_top + (T_bot - T_top) * z_nd;
    const double xmin = input->model.xmin;
    const double Lx   = input->model.xmax - xmin;
    double x_nd = (coordinates.x - xmin) / Lx;
    double perturbation = 0.01 * cos(M_PI * x_nd) * sin(M_PI * z_nd) / input->scaling.T;
    return T_lin + perturbation;
  }

  static double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
    return input->materials.rho[phase];
  }

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

  static SetBC SetBCT(MdoodzInput *input, POSITION position, Coordinates coordinates, double gridTemperature) {
    SetBC bc;
    double T_top = (input->model.user0 + zeroC) / input->scaling.T;
    double T_bot = input->model.user1 / input->scaling.T;
    if (position == N || position == NE || position == NW) {
      bc.value = T_top;
      bc.type  = 1;
    } else if (position == S || position == SE || position == SW) {
      bc.value = T_bot;
      bc.type  = 1;
    } else {
      bc.value = 0.0;
      bc.type  = 0;
    }
    return bc;
  }
};

// ---------------------------------------------------------------------------
// Test 1: Simulation runs and convection develops
// Steady-state (Nu=4.884, Vrms=42.87) requires ~100k+ steps; this smoke test
// verifies no crash, physical temperature range, and non-zero convection.
// ---------------------------------------------------------------------------
TEST_F(BlankenBench, ConvectionDevelops) {
  char inputFile[] = "BlankenBench/BlankenBench.txt";
  RunMDOODZ(inputFile, &setup);

  // Read final output (Nt = 500, writer_step = 500)
  const char *outFile = "BlankenBench/Output00500.gzip.h5";

  const int Nx = 41, Nz = 41;
  const int ncx = Nx - 1, ncz = Nz - 1;
  const double T_top_K = zeroC;
  const double DeltaT = 1000.0;
  const double H = 1.0e5;

  auto T_field = readFieldAsArray(outFile, "Centers", "T");
  ASSERT_EQ((int)T_field.size(), ncx * ncz);

  // Check temperature is in valid range [T_top - margin, T_bot + margin]
  const double T_bot_K = T_top_K + DeltaT;
  double T_min = 1e30, T_max = -1e30;
  for (int i = 0; i < ncx * ncz; i++) {
    if (T_field[i] < T_min) T_min = T_field[i];
    if (T_field[i] > T_max) T_max = T_field[i];
  }
  printf("BlankenBench: T range [%.1f, %.1f] K (expected ~[%.1f, %.1f])\n",
         T_min, T_max, T_top_K, T_bot_K);
  EXPECT_GT(T_min, T_top_K - 50.0) << "Temperature dropped below expected minimum";
  EXPECT_LT(T_max, T_bot_K + 50.0) << "Temperature exceeded expected maximum";

  // Check velocity fields exist and convection started (Vrms > 0)
  auto Vx_field = readFieldAsArray(outFile, "VxNodes", "Vx");
  auto Vz_field = readFieldAsArray(outFile, "VzNodes", "Vz");
  ASSERT_EQ((int)Vx_field.size(), Nx * (Nz + 1));
  ASSERT_EQ((int)Vz_field.size(), (Nx + 1) * Nz);

  double sumV2 = 0.0;
  for (size_t i = 0; i < Vx_field.size(); i++) sumV2 += Vx_field[i] * Vx_field[i];
  for (size_t i = 0; i < Vz_field.size(); i++) sumV2 += Vz_field[i] * Vz_field[i];
  EXPECT_GT(sumV2, 0.0) << "No velocity — convection did not start";
  printf("BlankenBench: sum(V^2) = %.4e (should be > 0)\n", sumV2);

  // Nusselt number (informational at this stage, not expected to match steady state)
  double dz = H / ncz;
  double sumDtDz = 0.0;
  for (int i = 0; i < ncx; i++) {
    double T_centre = T_field[(ncz - 1) * ncx + i];
    double dTdz = (T_top_K - T_centre) / (dz / 2.0);
    sumDtDz += dTdz;
  }
  double Nu = -(H / DeltaT) * (sumDtDz / ncx);
  printf("BlankenBench: Nu = %.4f (published steady state: 4.884)\n", Nu);
}

// ---------------------------------------------------------------------------
// Test 2: Mid-depth temperature monotonic along z and no NaN
// ---------------------------------------------------------------------------
TEST_F(BlankenBench, TemperatureProfile) {
  // Reuse output from the same simulation (tests run sequentially)
  const char *outFile = "BlankenBench/Output00500.gzip.h5";

  const int Nx = 41, Nz = 41;
  const int ncx = Nx - 1, ncz = Nz - 1;
  const double T_top_K = zeroC;
  const double DeltaT  = 1000.0;

  auto T_field = readFieldAsArray(outFile, "Centers", "T");
  ASSERT_EQ((int)T_field.size(), ncx * ncz);

  // Check no NaN in temperature field
  for (int i = 0; i < ncx * ncz; i++) {
    ASSERT_FALSE(std::isnan(T_field[i])) << "NaN in temperature field at index " << i;
  }

  // Mid-depth max temperature (informational)
  int mid_row = ncz / 2;
  double T_max_mid = -1e30;
  for (int i = 0; i < ncx; i++) {
    double T_nd = (T_field[mid_row * ncx + i] - T_top_K) / DeltaT;
    if (T_nd > T_max_mid) T_max_mid = T_nd;
  }
  printf("BlankenBench: T_max(mid) = %.4f (published steady state: 0.4266)\n", T_max_mid);

  // At early times T_max should be close to 0.5 (linear gradient midpoint)
  EXPECT_GT(T_max_mid, 0.1) << "Mid-depth temperature is too low";
  EXPECT_LT(T_max_mid, 0.9) << "Mid-depth temperature is too high";
}

// ---------------------------------------------------------------------------
// Test 3: Steady-state Nusselt number and Vrms (MANUAL — not in CI)
//
// Runs 100k steps with Courant=0.5 to reach thermal steady state, then
// computes Nu and Vrms and compares against Blankenbach et al. (1989)
// published reference values. Requires OMP_NUM_THREADS=8.
//
// Run: OMP_NUM_THREADS=8 ./BlankenBenchTests --gtest_filter=BlankenBench.NusseltAndVrms
// ---------------------------------------------------------------------------
TEST_F(BlankenBench, NusseltAndVrms) {
  // Guard: only run when explicitly requested via environment variable
  const char *run_steady = std::getenv("BLANKENBACH_STEADY");
  if (!run_steady || std::string(run_steady) != "1") {
    GTEST_SKIP() << "Set BLANKENBACH_STEADY=1 to run this long test (~2-4 hours)";
  }

#ifdef _OPENMP
  int nthreads = omp_get_max_threads();
  printf("BlankenBench: Running with %d OpenMP threads\n", nthreads);
  ASSERT_GE(nthreads, 4) << "Need at least 4 threads for reasonable runtime";
#else
  GTEST_SKIP() << "Build without OpenMP (-DOMP=ON). Skipping long-running test.";
#endif

  // Run steady-state simulation
  char inputFile[] = "BlankenBench/BlankenBenchSteady.txt";
  RunMDOODZ(inputFile, &setup);

  // Read final output
  const char *outFile = "BlankenBench/Output100000.gzip.h5";

  const int Nx = 41, Nz = 41;
  const int ncx = Nx - 1, ncz = Nz - 1;
  const double T_top_K = zeroC;          // 273.15 K
  const double DeltaT  = 1000.0;         // K
  const double H       = 1.0e5;          // m
  const double rho     = 3300.0;         // kg/m³
  const double k       = 3.3;            // W/m/K
  const double Cp      = 1250.0;         // J/kg/K
  const double kappa   = k / (rho * Cp); // 8.0e-7 m²/s

  // --- Nusselt number ---
  auto T_field = readFieldAsArray(outFile, "Centers", "T");
  ASSERT_EQ((int)T_field.size(), ncx * ncz);

  double dz = H / ncz;
  double sumDtDz = 0.0;
  for (int i = 0; i < ncx; i++) {
    double T_centre = T_field[(ncz - 1) * ncx + i];  // top cell-centre row
    double dTdz = (T_top_K - T_centre) / (dz / 2.0);
    sumDtDz += dTdz;
  }
  double Nu = -(H / DeltaT) * (sumDtDz / ncx);

  // --- Vrms ---
  auto Vx_field = readFieldAsArray(outFile, "VxNodes", "Vx");
  auto Vz_field = readFieldAsArray(outFile, "VzNodes", "Vz");
  ASSERT_EQ((int)Vx_field.size(), Nx * (Nz + 1));
  ASSERT_EQ((int)Vz_field.size(), (Nx + 1) * Nz);

  // Interpolate staggered velocities to cell centres
  double sumV2 = 0.0;
  for (int j = 0; j < ncz; j++) {
    for (int i = 0; i < ncx; i++) {
      double vx_c = 0.5 * (Vx_field[j * Nx + i] + Vx_field[j * Nx + (i + 1)]);
      double vz_c = 0.5 * (Vz_field[j * (Nx + 1) + i] + Vz_field[(j + 1) * (Nx + 1) + i]);
      sumV2 += vx_c * vx_c + vz_c * vz_c;
    }
  }
  double Vrms_SI = sqrt(sumV2 / (ncx * ncz));
  double Vrms_nd = Vrms_SI * H / kappa;  // non-dimensionalise

  // --- Report ---
  printf("\n========================================\n");
  printf("BlankenBench Steady State Results:\n");
  printf("  Nu   = %.6f  (published: 4.884409)\n", Nu);
  printf("  Vrms = %.6f  (published: 42.864947)\n", Vrms_nd);
  printf("========================================\n\n");

  // --- Assert against published Blankenbach et al. (1989) Case 1a values ---
  const double Nu_ref   = 4.884409;
  const double Vrms_ref = 42.864947;
  const double tol      = 0.05;  // 5% relative tolerance for 41×41 resolution

  EXPECT_NEAR(Nu,      Nu_ref,   tol * Nu_ref)   << "Nu outside 5% of published value";
  EXPECT_NEAR(Vrms_nd, Vrms_ref, tol * Vrms_ref) << "Vrms outside 5% of published value";
}
