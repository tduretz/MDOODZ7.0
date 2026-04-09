// SolViMarkerComparison.cpp — Experimental test: marker density vs L2 accuracy
// Compares 4x4 (baseline), 8x8, and 16x16 markers per cell on SolVi benchmark.
// NOT part of the regular CI suite — run separately for analysis.
//
// ===== FINDINGS (April 2026, WSL2 Ubuntu 22.04, gcc 11.4) =====
//
// NOTE: A coordinate-mapping bug was found and fixed in April 2026.
// The original code used wrong array dimensions (Nx*(Nz-1) and (Nx-1)*Nz)
// and z/x offsets for VxNodes/VzNodes. Correct dimensions are Nx*(Nz+1)
// and (Nx+1)*Nz, with half-cell offsets (iz-0.5)*dz and (ix-0.5)*dx.
// All numbers below reflect the CORRECTED mapping.
//
// 1. L2 errors at 51×51 (fixed grid, varying Nx_part × Nz_part):
//
//    markers/cell |   L2(Vx)   |   L2(Vz)   |   L2(P)    |  L2(sxxd)  |  time(ms)
//    -------------|------------|------------|------------|------------|----------
//    4×4  = 16    | 2.24e-02   | 2.24e-02   | 7.50e-01   | 4.33e-01   |      254
//    8×8  = 64    | 2.38e-02   | 2.38e-02   | 6.58e-01   | 4.27e-01   |      570
//    16×16= 256   | 2.31e-02   | 2.31e-02   | 6.59e-01   | 4.26e-01   |     1781
//
//    Improvement 4×4 → 8×8:    Vx=0.9×  P=1.1×  time=2.2× slower
//    Improvement 4×4 → 16×16:  Vx=1.0×  P=1.1×  time=7.0× slower
//
//    Note: L2(Vx) = L2(Vz) by symmetry of the pure-shear SolVi problem.
//    Before the coordinate fix, L2(Vz) was inflated ~20× due to the
//    wrong stride in the flat-array decomposition.
//
// 2. Convergence orders (41→81):
//
//    Field |  4×4 markers  |  8×8 markers
//    ------|---------------|-------------
//    Vx    |  1.02         |  0.98
//    P     |  0.75         |  0.59
//
// 3. Conclusion:
//    Increasing marker density gives NEGLIGIBLE accuracy improvement:
//    - Vx: no change (within noise)
//    - P:  ~12% lower L2 at fixed resolution, but convergence ORDER worsens
//          (0.75 → 0.59), meaning the error reduction does not carry over
//          to finer grids — it's a one-off constant offset, not a rate gain
//    - sxxd: ~1.5% better
//    - Cost: 4× memory and ~7× runtime for 16×16 vs 4×4
//
//    The L2 error is dominated by the FD truncation error at the viscosity
//    discontinuity (inclusion boundary), NOT by marker interpolation noise.
//    More markers only give a better average of the same wrong stencil.
//
//    To improve convergence order, one would need:
//    - Immersed-boundary / ghost-fluid corrections (expensive, up to 2nd order)
//    - Cell-face D-tensor harmonic averaging was attempted and FAILED (see below)
//
//    Bottom line: keep Nx_part=4, Nz_part=4 for CI tests. The default 16
//    markers/cell is sufficient for the accuracy this scheme can deliver.
//
// ===== eta_average FINDINGS =====
//
// The `eta_average` parameter (0=arithmetic, 1=harmonic, 2=geometric)
// controls how phase-mixture viscosity is averaged from markers to nodes
// WITHIN each cell.  This is NOT cell-face averaging in the FD stencil.
//
// 1. L2 errors at 51×51 (4×4 markers/cell, varying eta_average):
//
//    eta_average |   L2(Vx)   |   L2(Vz)   |   L2(P)    |  L2(sxxd)
//    ------------|------------|------------|------------|----------
//    0 arithmetic| 2.24e-02   | 2.24e-02   | 7.50e-01   | 4.33e-01
//    1 harmonic  | 1.02e-03   | 1.02e-03   | 9.19e-01   | 4.43e-01
//    2 geometric | 1.27e-02   | 1.27e-02   | 6.12e-01   | 4.11e-01
//
// 2. Convergence orders (41→81):
//
//    Field | arithmetic(0) | harmonic(1) | geometric(2)
//    ------|---------------|-------------|-------------
//    Vx    |  1.022        |  0.454      |  0.874
//    P     |  0.753        |  0.441      |  0.717
//
// 3. Conclusion:
//    - Geometric (2): best absolute P error (6.12e-1, 18% better than arith),
//      Vx convergence order slightly lower (0.87 vs 1.02)
//    - Harmonic (1): lowest absolute Vx error (1.0e-3!) but worst P and
//      worst convergence orders (0.45 for both Vx and P)
//    - Arithmetic (0): best convergence orders overall
//    - Key insight: eta_average controls marker-to-node phase mixing, NOT
//      cell-face interpolation in the stencil.  True cell-face harmonic
//      averaging would require modifying StokesAssemblyDecoupled.c.
//
// ===== CELL-FACE D-TENSOR HARMONIC AVERAGING: FAILED =====
//
// Attempted (April 2026): Replace D11E/D11W (and D33N/D33S, D22S/D22N,
// D33E/D33W) in StokesAssemblyDecoupled.c with their harmonic mean,
// gated on a new `cell_avg` parameter.
//
// Implementation: In Xmomentum_InnerNodesDecoupled, set
//   D11E = D11W = 2*D11E*D11W / (D11E + D11W)
// and similarly for D33N/D33S. Same pattern for Zmomentum.
//
// Result: CHOLMOD fails with "matrix not positive definite".
//
// Root cause: In this staggered-grid stencil, the Vx node sits ON the
// face between cells iPrW and iPrE. D11E and D11W are values at the
// cell CENTRES on each side. The stencil computes:
//   (D11E * eps_dot_E  -  D11W * eps_dot_W) / dx
// Setting D11E = D11W = harmonic_mean makes the operator locally
// constant across the two cells, which:
//   (a) loses the viscosity contrast information entirely
//   (b) breaks the SPD structure of the stiffness matrix (CHOLMOD fails)
//
// This is NOT the same as "harmonic face averaging" from Deubelbeiss &
// Kaus (2008), which averages viscosity during the marker-to-node
// interpolation (i.e., what eta_average=1 already does).
//
// Conclusion: In MDOODZ's stencil, D values at centres ARE the
// natural face representation. There is no separate face viscosity to
// average. Improving P convergence order beyond ~0.75 would require
// fundamentally different stencil modifications (e.g., immersed boundary,
// ghost-fluid, or embedded discontinuity methods).
//
// =================================================================

extern "C" {
#include "mdoodz.h"
}

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <chrono>
#include <gtest/gtest.h>
#include "TestHelpers.h"

// Same analytical solution as SolViBenchmarkTests.cpp
static void eval_anal_Dani(double *vx, double *vz, double *p, double *eta,
                           double *sxx, double *syy,
                           double x, double z, int ps,
                           double rc, double mm, double mc) {
  using Complex = std::complex<double>;
  const Complex II(0.0, 1.0);
  double gr, er, A;
  *eta = *vx = *vz = *p = 0;
  if (ps == 1) { gr = 0; er = -1; } else { gr = -2.0; er = 0; }
  A = mm * (mc - mm) / (mc + mm);
  Complex Z(x, z);
  if (sqrt(x * x + z * z) <= rc) {
    *p    = 0.0;
    Complex V_tot = (mm / (mc + mm)) * (II * gr + 2.0 * er) * conj(Z) - (II / 2.0) * gr * Z;
    *vx   = V_tot.real();
    *vz   = V_tot.imag();
    *eta  = mc;
    *sxx  = 4.0 * er * mc * mm / (mc + mm);
    *syy  = -4.0 * er * mc * mm / (mc + mm);
  } else {
    *p    = -2.0 * mm * (mc - mm) / (mc + mm) *
            (pow(rc, 2) / pow(Z, 2) * (II * gr + 2.0 * er)).real();
    Complex phi_z        = -(II / 2.0) * mm * gr * Z
                           - (II * gr + 2.0 * er) * A * pow(rc, 2) * pow(Z, -1.0);
    Complex d_phi_z      = -(II / 2.0) * mm * gr
                           + (II * gr + 2.0 * er) * A * pow(rc, 2) / pow(Z, 2);
    Complex psi_z        = (II * gr - 2.0 * er) * mm * Z
                           - (II * gr + 2.0 * er) * A * pow(rc, 4) * pow(Z, -3.0);
    Complex V_tot        = (phi_z - Z * conj(d_phi_z) - conj(psi_z)) / (2.0 * mm);
    *vx  = V_tot.real();
    *vz  = V_tot.imag();
    *eta = mm;
    Complex d_d_phi_z_z = -(II * gr + 2.0 * er) * A * pow(rc, 2) / pow(Z, 3);
    Complex d_psi_z     = (II * gr - 2.0 * er) * mm
                         + (II * gr + 2.0 * er) * A * pow(rc, 4) * pow(Z, -4.0);
    *syy = 2.0 * d_phi_z.real() + (Complex(x, -z) * d_d_phi_z_z + d_psi_z).real();
    *sxx = 4.0 * d_phi_z.real() - *syy;
  }
}

class SolViMarkerComparison : public ::testing::Test {
protected:
  static constexpr double RC = 0.2;
  static constexpr double MM = 1.0;
  static constexpr double MC = 1e3;

  MdoodzSetup setup;

  void SetUp() override {
    setup = {
            .SetParticles = new SetParticles_ff{
                    .SetPhase   = SetPhase,
                    .SetDensity = SetDensity,
            },
            .SetBCs = new SetBCs_ff{
                    .SetBCVx = SetBCVx,
                    .SetBCVz = SetBCVz,
            },
    };
  }

  static int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
    const double radius = instance->model.user1 / instance->scaling.L;
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius)
      return 1;
    return 0;
  }
  static double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
    return instance->materials.rho[phase];
  }

  static SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coord) {
    SetBC bc;
    const double radius = instance->model.user1 / instance->scaling.L;
    double Vx, Vz, P, eta, sxx, szz;
    if (position == W || position == E) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, MM, MC);
      bc.type  = 0; bc.value = Vx;
    } else if (position == S || position == SE || position == SW) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z + instance->model.dx / 2.0, 1, radius, MM, MC);
      bc.type  = 11; bc.value = Vx;
    } else if (position == N || position == NE || position == NW) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z - instance->model.dx / 2.0, 1, radius, MM, MC);
      bc.type  = 11; bc.value = Vx;
    } else { bc.type = -1; bc.value = 0.0; }
    return bc;
  }

  static SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
    SetBC bc;
    const double radius = instance->model.user1 / instance->scaling.L;
    double Vx, Vz, P, eta, sxx, szz;
    if (position == S || position == N) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, MM, MC);
      bc.type  = 0; bc.value = Vz;
    } else if (position == W || position == SW || position == NW) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x + instance->model.dx / 2.0, coord.z, 1, radius, MM, MC);
      bc.type  = 11; bc.value = Vz;
    } else if (position == E || position == SE || position == NE) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x - instance->model.dx / 2.0, coord.z, 1, radius, MM, MC);
      bc.type  = 11; bc.value = Vz;
    } else { bc.type = -1; bc.value = 0.0; }
    return bc;
  }

  static double computeSolViL2(const char *hdf5File, const char *group,
                                const char *dataset, int fieldType) {
    auto field = readFieldAsArray(hdf5File, group, dataset);
    int    Nx   = (int)getModelParam(hdf5File, 3);
    int    Nz   = (int)getModelParam(hdf5File, 4);
    double dx   = getModelParam(hdf5File, 5);
    double dz   = getModelParam(hdf5File, 6);
    double Lx   = getModelParam(hdf5File, 1);
    double Lz   = getModelParam(hdf5File, 2);
    double xmin = -Lx / 2.0;
    double zmin = -Lz / 2.0;
    std::vector<double> analytical(field.size());
    for (size_t k = 0; k < field.size(); k++) {
      double x, z, aVx, aVz, aP, aEta, aSxx, aSyy;
      int ix, iz;
      if (fieldType == 0) {
        int ncx = Nx - 1; ix = k % ncx; iz = k / ncx;
        x = xmin + (ix + 0.5) * dx; z = zmin + (iz + 0.5) * dz;
      } else if (fieldType == 1) { // VxNodes: Nx*(Nz+1)
        ix = k % Nx; iz = k / Nx;
        x = xmin + ix * dx; z = zmin + (iz - 0.5) * dz;
      } else { // VzNodes: (Nx+1)*Nz
        int nxVz = Nx + 1; ix = k % nxVz; iz = k / nxVz;
        x = xmin + (ix - 0.5) * dx; z = zmin + iz * dz;
      }
      eval_anal_Dani(&aVx, &aVz, &aP, &aEta, &aSxx, &aSyy, x, z, 1, RC, MM, MC);
      if (fieldType == 1) analytical[k] = aVx;
      else if (fieldType == 2) analytical[k] = aVz;
      else if (strcmp(dataset, "P") == 0) analytical[k] = aP;
      else analytical[k] = aSxx;
    }
    return computeL2Error(field, analytical);
  }

  struct RunResult {
    double L2_Vx, L2_Vz, L2_P, L2_sxxd;
    double elapsed_ms;
  };

  RunResult runAndMeasure(const char *txtFile, const char *subfolder) {
    RunResult res;
    char outputFile[256];
    snprintf(outputFile, sizeof(outputFile), "%s/Output00001.gzip.h5", subfolder);

    auto t0 = std::chrono::high_resolution_clock::now();
    RunMDOODZ(const_cast<char*>(txtFile), &setup);
    auto t1 = std::chrono::high_resolution_clock::now();

    res.elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
    res.L2_Vx   = computeSolViL2(outputFile, "VxNodes", "Vx",   1);
    res.L2_Vz   = computeSolViL2(outputFile, "VzNodes", "Vz",   2);
    res.L2_P    = computeSolViL2(outputFile, "Centers", "P",     0);
    res.L2_sxxd = computeSolViL2(outputFile, "Centers", "sxxd",  0);
    return res;
  }
};

// =========================================================================
// Test 1: Compare L2 at 51x51 with 4x4, 8x8, and 16x16 markers/cell
// =========================================================================
TEST_F(SolViMarkerComparison, MarkerDensity51x51) {
  printf("\n");
  printf("=================================================================\n");
  printf("  SolVi 51x51: Marker Density Comparison\n");
  printf("=================================================================\n");

  auto r4  = runAndMeasure("SolViBenchmark/SolViRes51.txt",     "SolViRes51");
  auto r8  = runAndMeasure("SolViBenchmark/SolViRes51_m8.txt",  "SolViRes51_m8");
  auto r16 = runAndMeasure("SolViBenchmark/SolViRes51_m16.txt", "SolViRes51_m16");

  printf("\n");
  printf("  markers/cell |   L2(Vx)   |   L2(Vz)   |   L2(P)    |  L2(sxxd)  |  time(ms)\n");
  printf("  -------------|------------|------------|------------|------------|----------\n");
  printf("  4x4  = 16    | %10.4e | %10.4e | %10.4e | %10.4e | %8.0f\n",
         r4.L2_Vx, r4.L2_Vz, r4.L2_P, r4.L2_sxxd, r4.elapsed_ms);
  printf("  8x8  = 64    | %10.4e | %10.4e | %10.4e | %10.4e | %8.0f\n",
         r8.L2_Vx, r8.L2_Vz, r8.L2_P, r8.L2_sxxd, r8.elapsed_ms);
  printf("  16x16= 256   | %10.4e | %10.4e | %10.4e | %10.4e | %8.0f\n",
         r16.L2_Vx, r16.L2_Vz, r16.L2_P, r16.L2_sxxd, r16.elapsed_ms);
  printf("\n");
  printf("  Improvement 4x4 -> 8x8:   Vx=%.1fx  P=%.1fx  time=%.1fx\n",
         r4.L2_Vx / r8.L2_Vx, r4.L2_P / r8.L2_P, r8.elapsed_ms / r4.elapsed_ms);
  printf("  Improvement 4x4 -> 16x16: Vx=%.1fx  P=%.1fx  time=%.1fx\n",
         r4.L2_Vx / r16.L2_Vx, r4.L2_P / r16.L2_P, r16.elapsed_ms / r4.elapsed_ms);
  printf("=================================================================\n\n");

  // All should produce valid results
  EXPECT_LT(r4.L2_Vx,  1.0);
  EXPECT_LT(r8.L2_Vx,  1.0);
  EXPECT_LT(r16.L2_Vx, 1.0);
}

// =========================================================================
// Test 2: Convergence order with 8x8 markers (compare to 4x4 baseline)
// =========================================================================
TEST_F(SolViMarkerComparison, ConvergenceOrder8x8) {
  printf("\n");
  printf("=================================================================\n");
  printf("  SolVi Grid Convergence: 8x8 markers/cell\n");
  printf("=================================================================\n");

  const char *txtFiles[] = {
    "SolViBenchmark/SolViRes41_m8.txt",
    "SolViBenchmark/SolViRes81_m8.txt",
  };
  const char *subfolders[] = {"SolViRes41_m8", "SolViRes81_m8"};

  double L2_Vx[2], L2_P[2], h[2], elapsed[2];

  for (int i = 0; i < 2; i++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    RunMDOODZ(const_cast<char*>(txtFiles[i]), &setup);
    auto t1 = std::chrono::high_resolution_clock::now();
    elapsed[i] = std::chrono::duration<double, std::milli>(t1 - t0).count();

    char outputFile[256];
    snprintf(outputFile, sizeof(outputFile), "%s/Output00001.gzip.h5", subfolders[i]);
    L2_Vx[i] = computeSolViL2(outputFile, "VxNodes", "Vx", 1);
    L2_P[i]  = computeSolViL2(outputFile, "Centers", "P",  0);
    h[i]     = getModelParam(outputFile, 5);
    printf("  Resolution h=%e: L2(Vx)=%e, L2(P)=%e, time=%.0fms\n",
           h[i], L2_Vx[i], L2_P[i], elapsed[i]);
  }

  double order_Vx_m8 = log(L2_Vx[0] / L2_Vx[1]) / log(h[0] / h[1]);
  double order_P_m8  = log(L2_P[0]  / L2_P[1])  / log(h[0] / h[1]);

  printf("\n  Convergence orders (8x8 markers, 41->81):\n");
  printf("    Vx order = %.3f\n", order_Vx_m8);
  printf("    P  order = %.3f\n", order_P_m8);
  printf("  (Baseline 4x4: Vx~0.97, P~0.75)\n");
  printf("=================================================================\n\n");

  // Report but don't fail — this is experimental
  EXPECT_GT(order_Vx_m8, 0.5);
  EXPECT_GT(order_P_m8,  0.3);
}

// =========================================================================
// Test 3: eta_average comparison (arithmetic=0, harmonic=1, geometric=2)
// =========================================================================
TEST_F(SolViMarkerComparison, EtaAveraging) {
  printf("\n");
  printf("=================================================================\n");
  printf("  SolVi: eta_average comparison (0=arith, 1=harm, 2=geom)\n");
  printf("=================================================================\n");

  // --- L2 at 51x51 for each averaging ---
  const char *txt51[] = {
    "SolViBenchmark/SolViRes51.txt",       // arithmetic (0)
    "SolViBenchmark/SolViRes51_avg1.txt",  // harmonic (1)
    "SolViBenchmark/SolViRes51_avg2.txt",  // geometric (2)
  };
  const char *sub51[] = {"SolViRes51", "SolViRes51_avg1", "SolViRes51_avg2"};
  const char *names[] = {"arithmetic", "harmonic  ", "geometric "};

  double L2_Vx_51[3], L2_P_51[3], L2_sxxd_51[3], elapsed_51[3];

  for (int i = 0; i < 3; i++) {
    auto t0 = std::chrono::high_resolution_clock::now();
    RunMDOODZ(const_cast<char*>(txt51[i]), &setup);
    auto t1 = std::chrono::high_resolution_clock::now();
    elapsed_51[i] = std::chrono::duration<double, std::milli>(t1 - t0).count();

    char out[256];
    snprintf(out, sizeof(out), "%s/Output00001.gzip.h5", sub51[i]);
    L2_Vx_51[i]   = computeSolViL2(out, "VxNodes", "Vx",    1);
    L2_P_51[i]    = computeSolViL2(out, "Centers", "P",     0);
    L2_sxxd_51[i] = computeSolViL2(out, "Centers", "sxxd",  0);
  }

  printf("\n  51x51 L2 errors by averaging method:\n");
  printf("  avg method   |   L2(Vx)   |   L2(P)    |  L2(sxxd)  |  time(ms)\n");
  printf("  -------------|------------|------------|------------|----------\n");
  for (int i = 0; i < 3; i++) {
    printf("  %s | %10.4e | %10.4e | %10.4e | %8.0f\n",
           names[i], L2_Vx_51[i], L2_P_51[i], L2_sxxd_51[i], elapsed_51[i]);
  }

  // --- Convergence order (41→81) for each averaging ---
  const char *txt41[] = {
    "SolViBenchmark/SolViRes41.txt",
    "SolViBenchmark/SolViRes41_avg1.txt",
    "SolViBenchmark/SolViRes41_avg2.txt",
  };
  const char *txt81[] = {
    "SolViBenchmark/SolViRes81.txt",
    "SolViBenchmark/SolViRes81_avg1.txt",
    "SolViBenchmark/SolViRes81_avg2.txt",
  };
  const char *sub41[] = {"SolViRes41", "SolViRes41_avg1", "SolViRes41_avg2"};
  const char *sub81[] = {"SolViRes81", "SolViRes81_avg1", "SolViRes81_avg2"};

  printf("\n  Convergence orders (41->81):\n");
  printf("  avg method   |  Vx order  |   P order\n");
  printf("  -------------|------------|----------\n");

  for (int i = 0; i < 3; i++) {
    RunMDOODZ(const_cast<char*>(txt41[i]), &setup);
    RunMDOODZ(const_cast<char*>(txt81[i]), &setup);

    char out41[256], out81[256];
    snprintf(out41, sizeof(out41), "%s/Output00001.gzip.h5", sub41[i]);
    snprintf(out81, sizeof(out81), "%s/Output00001.gzip.h5", sub81[i]);

    double L2_Vx_41 = computeSolViL2(out41, "VxNodes", "Vx", 1);
    double L2_Vx_81 = computeSolViL2(out81, "VxNodes", "Vx", 1);
    double L2_P_41  = computeSolViL2(out41, "Centers", "P",  0);
    double L2_P_81  = computeSolViL2(out81, "Centers", "P",  0);
    double h41      = getModelParam(out41, 5);
    double h81      = getModelParam(out81, 5);

    double order_Vx = log(L2_Vx_41 / L2_Vx_81) / log(h41 / h81);
    double order_P  = log(L2_P_41  / L2_P_81)  / log(h41 / h81);

    printf("  %s |   %6.3f   |   %6.3f\n", names[i], order_Vx, order_P);
  }

  printf("=================================================================\n\n");

  // All should produce valid results
  EXPECT_LT(L2_Vx_51[0], 1.0);
  EXPECT_LT(L2_Vx_51[1], 1.0);
  EXPECT_LT(L2_Vx_51[2], 1.0);
}
