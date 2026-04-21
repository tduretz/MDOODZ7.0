// add-gmg-stokes-defence §6.4 — V-cycle animation dump driver.
//
// Runs `TESTS/SolViBenchmark/SolViRes41_gmg_vcycle_anim.txt` (a 41×41
// SolCx-analog: piecewise-constant η with a 1000× contrast inclusion)
// through the GMG-FGMRES path with `gmg_dump_vcycle = 1`, producing
//
//   SolViRes41_gmg_vcycle_anim/vcycle_dump/level_{N}_{step}_{phase}_{seq}.h5
//
// under the test working directory. The fixture carries the `Nb_phases = 2`
// inclusion geometry from `SolViRes51_gmg.txt`; we reuse that test's
// SetParticles + SetBCs callbacks (analytic Dani boundary values) via a
// header include, so this driver stays at ~60 LOC.
//
// Why a dedicated test instead of piggy-backing on `GmgStokesEquivalence`:
//   * The one-shot dumper latches (`g_vcycle_dumper.fired = 1`) after
//     finish, so within a single process only one Stokes solve can
//     emit snapshots. Sharing a binary with other GMG tests would
//     make dump emission depend on test ordering, which is brittle.
//   * This test has no assertions — it just produces data. Making
//     it a separate binary keeps ctest's semantics clean.
//
// The test is tagged `experimental` so `ctest -L experimental` runs it
// explicitly. Default `ctest` skips it (3-second run but zero assertion
// value, only useful when regenerating defence figures).

extern "C" {
#include "mdoodz.h"
}

#include <cstring>
#include <gtest/gtest.h>

// SolVi (Schmid-Podladchikov) inclusion setup — ported from
// GmgStokesEquivalence.cpp and SolViBenchmarkTests.cpp. A minimal
// inline copy so this driver has no test-file coupling.
#include <complex>

static void eval_anal_Dani(double *vx, double *vz, double *p, double *eta,
                           double *sxx, double *syy,
                           double x, double z, int ps,
                           double rc, double mm, double mc) {
  using Complex = std::complex<double>;
  const Complex II(0.0, 1.0);
  double gr, er, A;
  *eta = *vx = *vz = *p = 0;
  if (ps == 1) { gr = 0;  er = -1; }
  else         { gr = -2; er =  0; }
  A = mm * (mc - mm) / (mc + mm);
  Complex Z(x, z);
  if (sqrt(x * x + z * z) <= rc) {
    *p    = 0.0;
    Complex V_tot = (mm / (mc + mm)) * (II * gr + 2.0 * er) * conj(Z)
                    - (II / 2.0) * gr * Z;
    *vx  = V_tot.real();
    *vz  = V_tot.imag();
    *eta = mc;
    *sxx =  4.0 * er * mc * mm / (mc + mm);
    *syy = -4.0 * er * mc * mm / (mc + mm);
  } else {
    *p = -2.0 * mm * (mc - mm) / (mc + mm)
         * (pow(rc, 2) / pow(Z, 2) * (II * gr + 2.0 * er)).real();
    Complex phi_z   = -(II / 2.0) * mm * gr * Z
                      - (II * gr + 2.0 * er) * A * pow(rc, 2) * pow(Z, -1.0);
    Complex d_phi_z = -(II / 2.0) * mm * gr
                      + (II * gr + 2.0 * er) * A * pow(rc, 2) / pow(Z, 2);
    Complex conj_d_phi_z = conj(d_phi_z);
    Complex psi_z      = (II * gr - 2.0 * er) * mm * Z
                         - (II * gr + 2.0 * er) * A * pow(rc, 4) * pow(Z, -3.0);
    Complex conj_psi_z = conj(psi_z);
    Complex V_tot = (phi_z - Z * conj_d_phi_z - conj_psi_z) / (2.0 * mm);
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

static constexpr double RC = 0.2;
static constexpr double MM = 1.0;
static constexpr double MC = 1e3;

static int SetPhase(MdoodzInput *instance, Coordinates c) {
  const double radius = instance->model.user1 / instance->scaling.L;
  return (c.x * c.x + c.z * c.z < radius * radius) ? 1 : 0;
}

static double SetDensity(MdoodzInput *instance, Coordinates, int phase) {
  return instance->materials.rho[phase];
}

static SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates c) {
  SetBC bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  double Vx, Vz, P, eta, sxx, szz;
  if (position == W || position == E) {
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, c.x, c.z, 1, radius, MM, MC);
    bc.type = 0; bc.value = Vx;
  } else if (position == S || position == SE || position == SW) {
    double x = c.x; double z = c.z + instance->model.dx / 2.0;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vx;
  } else if (position == N || position == NE || position == NW) {
    double x = c.x; double z = c.z - instance->model.dx / 2.0;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vx;
  } else { bc.type = -1; bc.value = 0.0; }
  return bc;
}

static SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates c) {
  SetBC bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  double Vx, Vz, P, eta, sxx, szz;
  if (position == S || position == N) {
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, c.x, c.z, 1, radius, MM, MC);
    bc.type = 0; bc.value = Vz;
  } else if (position == W || position == SW || position == NW) {
    double x = c.x + instance->model.dx / 2.0; double z = c.z;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vz;
  } else if (position == E || position == SE || position == NE) {
    double x = c.x - instance->model.dx / 2.0; double z = c.z;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vz;
  } else { bc.type = -1; bc.value = 0.0; }
  return bc;
}

TEST(VcycleAnimGenerate, RunsSolCx41AndEmitsDumps) {
  MdoodzSetup setup = {
      .SetParticles = new SetParticles_ff{
          .SetPhase   = SetPhase,
          .SetDensity = SetDensity,
      },
      .SetBCs = new SetBCs_ff{
          .SetBCVx = SetBCVx,
          .SetBCVz = SetBCVz,
      },
  };
  char inputName[] = "SolViBenchmark/SolViRes41_gmg_vcycle_anim.txt";
  RunMDOODZ(inputName, &setup);
  // Zero assertions: success == the solve produced
  //   SolViRes41_gmg_vcycle_anim/vcycle_dump/*.h5
  // downstream (benchmarks/vcycle_vis/make_gif.py) verifies the tree.
  SUCCEED();
}
