extern "C" {
#include "mdoodz.h"
}

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <gtest/gtest.h>
#include <hdf5.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Angular velocity: one full revolution per second
static const double OMEGA = 2.0 * M_PI;

// Vortex amplitude
static const double V0 = 2.0;

// Disk parameters (non-dimensional, same as scaling L=1)
static const double DISK_CX = 0.25;
static const double DISK_CZ = 0.0;
static const double DISK_R  = 0.1;

// ---- HDF5 helpers for char (phase) and float arrays ----

static std::vector<signed char> readCharField(const char *fname, const char *group, const char *dset_name) {
  hid_t file  = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t grp   = H5Gopen(file, group, H5P_DEFAULT);
  hid_t dset  = H5Dopen(grp, dset_name, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  std::vector<signed char> buf(dims[0]);
  H5Dread(dset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(grp);
  H5Fclose(file);
  return buf;
}

static std::vector<float> readFloatField(const char *fname, const char *group, const char *dset_name) {
  hid_t file  = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
  hid_t grp   = H5Gopen(file, group, H5P_DEFAULT);
  hid_t dset  = H5Dopen(grp, dset_name, H5P_DEFAULT);
  hid_t space = H5Dget_space(dset);
  hsize_t dims[1];
  H5Sget_simple_extent_dims(space, dims, NULL);
  std::vector<float> buf(dims[0]);
  H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
  H5Sclose(space);
  H5Dclose(dset);
  H5Gclose(grp);
  H5Fclose(file);
  return buf;
}

// Compute L2 error between two char fields interpreted as doubles
static double compoL2(const std::vector<signed char> &a, const std::vector<signed char> &b) {
  double sumDiff2 = 0.0;
  size_t n = a.size();
  for (size_t i = 0; i < n; i++) {
    double d = (double)a[i] - (double)b[i];
    sumDiff2 += d * d;
  }
  return sqrt(sumDiff2 / (double)n);
}

// Count mismatched cells between two char fields
static int compoMismatch(const std::vector<signed char> &a, const std::vector<signed char> &b) {
  int count = 0;
  size_t n = a.size();
  for (size_t i = 0; i < n; i++) {
    if (a[i] != b[i]) count++;
  }
  return count;
}

// ---- MDOODZ callbacks ----

static int SetPhase(MdoodzInput *input, Coordinates coords) {
  double x = coords.x - DISK_CX / input->scaling.L;
  double z = coords.z - DISK_CZ / input->scaling.L;
  double r = DISK_R / input->scaling.L;
  if (x * x + z * z < r * r) {
    return 1;
  }
  return 0;
}

static double SetDensity(MdoodzInput *input, Coordinates coords, int phase) {
  return input->materials.rho[phase];
}

// Rigid-body rotation: Vx = -omega * z
static double SetHorizontalVelocity(MdoodzInput *input, Coordinates coords) {
  return -OMEGA * coords.z * input->scaling.L / input->scaling.V;
}

// Rigid-body rotation: Vz = +omega * x
static double SetVerticalVelocity(MdoodzInput *input, Coordinates coords) {
  return OMEGA * coords.x * input->scaling.L / input->scaling.V;
}

// Single vortex (divergence-free, strong shear):
//   Vx = -V0 * sin(2*pi*x) * cos(2*pi*z)
//   Vz =  V0 * cos(2*pi*x) * sin(2*pi*z)
static double SetHorizontalVelocityVortex(MdoodzInput *input, Coordinates coords) {
  double x = coords.x * input->scaling.L;  // back to SI (which is = non-dim here)
  double z = coords.z * input->scaling.L;
  return -V0 * sin(2.0 * M_PI * x) * cos(2.0 * M_PI * z) / input->scaling.V;
}

static double SetVerticalVelocityVortex(MdoodzInput *input, Coordinates coords) {
  double x = coords.x * input->scaling.L;
  double z = coords.z * input->scaling.L;
  return V0 * cos(2.0 * M_PI * x) * sin(2.0 * M_PI * z) / input->scaling.V;
}

// ---- Test fixture ----

class RotationAdvection : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
      .SetParticles = new SetParticles_ff{
        .SetHorizontalVelocity = SetHorizontalVelocity,
        .SetVerticalVelocity   = SetVerticalVelocity,
        .SetPhase              = SetPhase,
        .SetDensity            = SetDensity,
      },
      .SetBCs = new SetBCs_ff{
        .SetBCVx = SetPureShearBCVx,
        .SetBCVz = SetPureShearBCVz,
      },
    };
  }

  void TearDown() override {
    delete setup.SetParticles;
    delete setup.SetBCs;
  }

  // Run one simulation, return L2 error between initial and final composition
  double RunAndMeasureL2(const char *configName) {
    char *inputName;
    asprintf(&inputName, "RotationAdvection/%s.txt", configName);
    RunMDOODZ(inputName, &setup);
    free(inputName);

    char *fileInit;
    asprintf(&fileInit, "%s/Output00000.gzip.h5", configName);
    char *fileFinal;
    asprintf(&fileFinal, "%s/Output01000.gzip.h5", configName);

    // Read high-resolution composition maps
    auto compo_init  = readCharField(fileInit,  "VizGrid", "compo_hr");
    auto compo_final = readCharField(fileFinal, "VizGrid", "compo_hr");

    double l2 = compoL2(compo_init, compo_final);
    int mismatch = compoMismatch(compo_init, compo_final);
    printf("[%s] L2 = %.6e, mismatched cells = %d / %zu (%.2f%%)\n",
           configName, l2, mismatch, compo_init.size(),
           100.0 * mismatch / (double)compo_init.size());

    free(fileInit);
    free(fileFinal);
    return l2;
  }
};

TEST_F(RotationAdvection, CompareReseedModes) {
  double l2_mode0 = RunAndMeasureL2("RotationMode0");
  double l2_mode1 = RunAndMeasureL2("RotationMode1");
  double l2_mode2 = RunAndMeasureL2("RotationMode2");

  printf("\n=== Rotation Advection L2 Summary ===\n");
  printf("  Mode 0 (legacy):   L2 = %.6e\n", l2_mode0);
  printf("  Mode 1 (current):  L2 = %.6e\n", l2_mode1);
  printf("  Mode 2 (improved): L2 = %.6e\n", l2_mode2);
  printf("=====================================\n");

  // All modes should have reasonably small L2 error (< 0.1 = 10% RMS)
  EXPECT_LT(l2_mode0, 0.1) << "Mode 0 L2 error too large";
  EXPECT_LT(l2_mode1, 0.1) << "Mode 1 L2 error too large";
  EXPECT_LT(l2_mode2, 0.1) << "Mode 2 L2 error too large";

  // Mode 2 should be no worse than mode 0
  EXPECT_LE(l2_mode2, l2_mode0 * 1.5) << "Mode 2 should not be much worse than mode 0";
}

// =====================================================================
// Vortex Advection Test — strong shear, sparse particles
// =====================================================================

class VortexAdvection : public ::testing::Test {
protected:
  MdoodzSetup setup;

  void SetUp() override {
    setup = {
      .SetParticles = new SetParticles_ff{
        .SetHorizontalVelocity = SetHorizontalVelocityVortex,
        .SetVerticalVelocity   = SetVerticalVelocityVortex,
        .SetPhase              = SetPhase,
        .SetDensity            = SetDensity,
      },
      .SetBCs = new SetBCs_ff{
        .SetBCVx = SetPureShearBCVx,
        .SetBCVz = SetPureShearBCVz,
      },
    };
  }

  void TearDown() override {
    delete setup.SetParticles;
    delete setup.SetBCs;
  }

  // Run one vortex simulation, read final composition
  std::vector<signed char> RunAndGetFinal(const char *configName) {
    char *inputName;
    asprintf(&inputName, "RotationAdvection/%s.txt", configName);
    RunMDOODZ(inputName, &setup);
    free(inputName);

    char *fileFinal;
    asprintf(&fileFinal, "%s/Output00500.gzip.h5", configName);
    auto compo = readCharField(fileFinal, "VizGrid", "compo_hr");
    free(fileFinal);
    return compo;
  }
};

TEST_F(VortexAdvection, CompareReseedModes) {
  auto compo0 = RunAndGetFinal("VortexMode0");
  auto compo1 = RunAndGetFinal("VortexMode1");
  auto compo2 = RunAndGetFinal("VortexMode2");

  // Read initial field for reference
  auto compo_init = readCharField("VortexMode0/Output00000.gzip.h5", "VizGrid", "compo_hr");

  // Pairwise L2 between modes (measures how much reseeding strategy diverges)
  double l2_01 = compoL2(compo0, compo1);
  double l2_02 = compoL2(compo0, compo2);
  double l2_12 = compoL2(compo1, compo2);

  // L2 vs initial (measures total deformation + reseeding artifact)
  double l2_0_init = compoL2(compo0, compo_init);
  double l2_1_init = compoL2(compo1, compo_init);
  double l2_2_init = compoL2(compo2, compo_init);

  int mis_01 = compoMismatch(compo0, compo1);
  int mis_02 = compoMismatch(compo0, compo2);
  int mis_12 = compoMismatch(compo1, compo2);

  printf("\n=== Vortex Advection — Inter-Mode Comparison ===\n");
  printf("  Mode 0 vs Mode 1:  L2 = %.6e, mismatched = %d / %zu (%.2f%%)\n",
         l2_01, mis_01, compo0.size(), 100.0 * mis_01 / (double)compo0.size());
  printf("  Mode 0 vs Mode 2:  L2 = %.6e, mismatched = %d / %zu (%.2f%%)\n",
         l2_02, mis_02, compo0.size(), 100.0 * mis_02 / (double)compo0.size());
  printf("  Mode 1 vs Mode 2:  L2 = %.6e, mismatched = %d / %zu (%.2f%%)\n",
         l2_12, mis_12, compo0.size(), 100.0 * mis_12 / (double)compo0.size());
  printf("\n  vs initial (deformation measure):\n");
  printf("  Mode 0: L2 = %.6e\n", l2_0_init);
  printf("  Mode 1: L2 = %.6e\n", l2_1_init);
  printf("  Mode 2: L2 = %.6e\n", l2_2_init);
  printf("=================================================\n");

  // Modes should produce similar results (< 10% RMS difference)
  EXPECT_LT(l2_12, 0.10) << "Mode 1 vs Mode 2 diverged too much";
}
