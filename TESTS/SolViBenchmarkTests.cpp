extern "C" {
#include "mdoodz.h"
}

#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <vector>
#include <gtest/gtest.h>
#include "TestHelpers.h"

// Ported from SETS/ViscousInclusion.c — Schmid & Podladchikov (2003)
// Uses std::complex<double> for C++ compatibility
static void eval_anal_Dani(double *vx, double *vz, double *p, double *eta,
                           double *sxx, double *syy,
                           double x, double z, int ps,
                           double rc, double mm, double mc) {
  using Complex = std::complex<double>;
  const Complex II(0.0, 1.0);
  double gr, er, A;

  *eta = *vx = *vz = *p = 0;

  if (ps == 1) {
    gr = 0;    // Pure shear
    er = -1;   // Strain rate
  } else {
    gr = -2.0; // Simple shear
    er = 0;
  }
  A = mm * (mc - mm) / (mc + mm);

  Complex Z(x, z);

  if (sqrt(x * x + z * z) <= rc) {
    // INSIDE CLAST
    *p    = 0.0;
    Complex V_tot = (mm / (mc + mm)) * (II * gr + 2.0 * er) * conj(Z) - (II / 2.0) * gr * Z;
    *vx   = V_tot.real();
    *vz   = V_tot.imag();
    *eta  = mc;
    // Correct analytical deviatoric stress inside clast
    *sxx  = 4.0 * er * mc * mm / (mc + mm);
    *syy  = -4.0 * er * mc * mm / (mc + mm);
  } else {
    // OUTSIDE CLAST (MATRIX)
    *p    = -2.0 * mm * (mc - mm) / (mc + mm) *
            (pow(rc, 2) / pow(Z, 2) * (II * gr + 2.0 * er)).real();

    Complex phi_z        = -(II / 2.0) * mm * gr * Z
                           - (II * gr + 2.0 * er) * A * pow(rc, 2) * pow(Z, -1.0);
    Complex d_phi_z      = -(II / 2.0) * mm * gr
                           + (II * gr + 2.0 * er) * A * pow(rc, 2) / pow(Z, 2);
    Complex conj_d_phi_z = conj(d_phi_z);
    Complex psi_z        = (II * gr - 2.0 * er) * mm * Z
                           - (II * gr + 2.0 * er) * A * pow(rc, 4) * pow(Z, -3.0);
    Complex conj_psi_z   = conj(psi_z);
    Complex V_tot        = (phi_z - Z * conj_d_phi_z - conj_psi_z) / (2.0 * mm);

    *vx  = V_tot.real();
    *vz  = V_tot.imag();
    *eta = mm;

    // Evaluate stresses
    Complex d_d_phi_z_z = -(II * gr + 2.0 * er) * A * pow(rc, 2) / pow(Z, 3);
    Complex d_psi_z     = (II * gr - 2.0 * er) * mm
                         + (II * gr + 2.0 * er) * A * pow(rc, 4) * pow(Z, -4.0);
    *syy = 2.0 * d_phi_z.real() + (Complex(x, -z) * d_d_phi_z_z + d_psi_z).real();
    *sxx = 4.0 * d_phi_z.real() - *syy;
  }
}

// SolVi test fixture
class SolViBenchmark : public ::testing::Test {
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
    if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
      return 1;
    }
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
      bc.type  = 0;
      bc.value = Vx;
    } else if (position == S || position == SE || position == SW) {
      double x = coord.x;
      double z = coord.z + instance->model.dx / 2.0;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
      bc.type  = 11;
      bc.value = Vx;
    } else if (position == N || position == NE || position == NW) {
      double x = coord.x;
      double z = coord.z - instance->model.dx / 2.0;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
      bc.type  = 11;
      bc.value = Vx;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
    return bc;
  }

  static SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
    SetBC bc;
    const double radius = instance->model.user1 / instance->scaling.L;
    double Vx, Vz, P, eta, sxx, szz;
    if (position == S || position == N) {
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, MM, MC);
      bc.type  = 0;
      bc.value = Vz;
    } else if (position == W || position == SW || position == NW) {
      double x = coord.x + instance->model.dx / 2.0;
      double z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
      bc.type  = 11;
      bc.value = Vz;
    } else if (position == E || position == SE || position == NE) {
      double x = coord.x - instance->model.dx / 2.0;
      double z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
      bc.type  = 11;
      bc.value = Vz;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
    return bc;
  }

  // Compute L2 error for a SolVi field
  // fieldType: 0=Centers, 1=VxNodes, 2=VzNodes
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
      double x, z;
      double aVx, aVz, aP, aEta, aSxx, aSyy;
      int ix, iz;

      if (fieldType == 0) {        // Centers: (Nx-1)*(Nz-1)
        int ncx = Nx - 1;
        ix = k % ncx;
        iz = k / ncx;
        x  = xmin + (ix + 0.5) * dx;
        z  = zmin + (iz + 0.5) * dz;
      } else if (fieldType == 1) { // VxNodes: Nx*(Nz-1)
        ix = k % Nx;
        iz = k / Nx;
        x  = xmin + ix * dx;
        z  = zmin + (iz + 0.5) * dz;
      } else {                      // VzNodes: (Nx-1)*Nz
        int ncx = Nx - 1;
        ix = k % ncx;
        iz = k / ncx;
        x  = xmin + (ix + 0.5) * dx;
        z  = zmin + iz * dz;
      }

      eval_anal_Dani(&aVx, &aVz, &aP, &aEta, &aSxx, &aSyy, x, z, 1, RC, MM, MC);

      if (fieldType == 1) {
        analytical[k] = aVx;
      } else if (fieldType == 2) {
        analytical[k] = aVz;
      } else if (strcmp(dataset, "P") == 0) {
        analytical[k] = aP;
      } else {
        analytical[k] = aSxx; // sxxd
      }
    }
    return computeL2Error(field, analytical);
  }
};

TEST_F(SolViBenchmark, L2Error) {
  char inputName[] = "SolViBenchmark/SolViRes51.txt";
  RunMDOODZ(inputName, &setup);

  const char *file = "SolViRes51/Output00001.gzip.h5";

  double L2_Vx   = computeSolViL2(file, "VxNodes", "Vx",   1);
  double L2_Vz   = computeSolViL2(file, "VzNodes", "Vz",   2);
  double L2_P    = computeSolViL2(file, "Centers", "P",     0);
  double L2_sxxd = computeSolViL2(file, "Centers", "sxxd",  0);

  printf("SolVi 51x51 L2 errors: Vx=%e, Vz=%e, P=%e, sxxd=%e\n",
         L2_Vx, L2_Vz, L2_P, L2_sxxd);

  EXPECT_LT(L2_Vx,   5e-1);  // marker-in-cell gives ~1e-2 at 51x51
  EXPECT_LT(L2_Vz,   5e-1);
  EXPECT_LT(L2_P,    2.0);   // pressure less accurate near inclusion
  EXPECT_LT(L2_sxxd, 2.0);
}

TEST_F(SolViBenchmark, GridConvergence) {
  char txtFile21[] = "SolViBenchmark/SolViRes21.txt";
  char txtFile41[] = "SolViBenchmark/SolViRes41.txt";
  char txtFile81[] = "SolViBenchmark/SolViRes81.txt";
  char *txtFiles[] = {txtFile21, txtFile41, txtFile81};
  const char *subfolders[] = {"SolViRes21", "SolViRes41", "SolViRes81"};
  double L2_Vx[3], L2_P[3], h[3];

  for (int i = 0; i < 3; i++) {
    RunMDOODZ(txtFiles[i], &setup);
    char outputFile[256];
    snprintf(outputFile, sizeof(outputFile), "%s/Output00001.gzip.h5", subfolders[i]);
    L2_Vx[i] = computeSolViL2(outputFile, "VxNodes", "Vx", 1);
    L2_P[i]  = computeSolViL2(outputFile, "Centers", "P",  0);
    h[i]     = getModelParam(outputFile, 5);
    printf("Resolution h=%e: L2(Vx)=%e, L2(P)=%e\n", h[i], L2_Vx[i], L2_P[i]);
  }

  // Convergence order between finest pair (41→81)
  double order_Vx = log(L2_Vx[1] / L2_Vx[2]) / log(h[1] / h[2]);
  double order_P  = log(L2_P[1]  / L2_P[2])  / log(h[1] / h[2]);

  printf("Convergence orders (41->81): Vx=%.2f, P=%.2f\n", order_Vx, order_P);

  EXPECT_GE(order_Vx, 0.7);  // marker-in-cell: ~first order near inclusion
  EXPECT_GE(order_P,  0.4);   // pressure converges slower
}
