// add-gmg-upleg-fix §4 D1 — standalone diagnostic probe for gmg_levels
// bisection (§4.1) and per-stage residual snapshots (§4.2) on 201²
// SolVi-Dani with 1000× viscosity contrast.
//
// This is NOT a regression test. It is a developer-driven experiment
// binary: it accepts an input .txt path as argv[1], installs the
// SolVi-Dani setup callbacks (same as `GmgConvergesAtDefaultLevels.cpp`),
// and calls `RunMDOODZ`. All MDOODZ logging goes straight to stdout,
// so the driver bash script in §4.1 redirects per-level logs into
// `openspec/changes/add-gmg-upleg-fix/d1_probe_a/<label>.log` for later
// analysis.
//
// To drive a constant-viscosity run (§4.3 Probe C discriminator), pass
// an input file whose `eta0` is 1 for both phases — the Dani BCs
// degenerate to pure shear, and the setup stays geometrically identical
// to the variable-viscosity runs.
//
// Build: `cmake --build cmake-build --target GmgLevelsSweepProbe`
// Run:   `./GmgLevelsSweepProbe SolViBenchmark/D1ProbeA_Lvl4.txt`

extern "C" {
#include "mdoodz.h"
}

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>

static constexpr double MM = 1.0;
static constexpr double MC = 1e3;

static void EvalDani(double *vx, double *vz, double *p, double *eta,
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
  if (std::sqrt(x * x + z * z) <= rc) {
    *p = 0.0;
    Complex V_tot = (mm / (mc + mm)) * (II * gr + 2.0 * er) * std::conj(Z)
                    - (II / 2.0) * gr * Z;
    *vx = V_tot.real(); *vz = V_tot.imag(); *eta = mc;
    *sxx =  4.0 * er * mc * mm / (mc + mm);
    *syy = -4.0 * er * mc * mm / (mc + mm);
  } else {
    *p = -2.0 * mm * (mc - mm) / (mc + mm)
         * (std::pow(rc, 2) / std::pow(Z, 2) * (II * gr + 2.0 * er)).real();
    Complex phi_z   = -(II / 2.0) * mm * gr * Z
                      - (II * gr + 2.0 * er) * A * std::pow(rc, 2) * std::pow(Z, -1.0);
    Complex d_phi_z = -(II / 2.0) * mm * gr
                      + (II * gr + 2.0 * er) * A * std::pow(rc, 2) / std::pow(Z, 2);
    Complex conj_d_phi_z = std::conj(d_phi_z);
    Complex psi_z      = (II * gr - 2.0 * er) * mm * Z
                         - (II * gr + 2.0 * er) * A * std::pow(rc, 4) * std::pow(Z, -3.0);
    Complex conj_psi_z = std::conj(psi_z);
    Complex V_tot = (phi_z - Z * conj_d_phi_z - conj_psi_z) / (2.0 * mm);
    *vx = V_tot.real(); *vz = V_tot.imag(); *eta = mm;
    Complex d_d_phi_z_z = -(II * gr + 2.0 * er) * A * std::pow(rc, 2) / std::pow(Z, 3);
    Complex d_psi_z     = (II * gr - 2.0 * er) * mm
                          + (II * gr + 2.0 * er) * A * std::pow(rc, 4) * std::pow(Z, -4.0);
    *syy = 2.0 * d_phi_z.real() + (Complex(x, -z) * d_d_phi_z_z + d_psi_z).real();
    *sxx = 4.0 * d_phi_z.real() - *syy;
  }
}

static int SetPhase(MdoodzInput *input, Coordinates c) {
  const double radius = input->model.user1 / input->scaling.L;
  return (c.x * c.x + c.z * c.z < radius * radius) ? 1 : 0;
}

static double SetDensity(MdoodzInput *input, Coordinates, int phase) {
  return input->materials.rho[phase];
}

static SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates c) {
  SetBC bc;
  const double radius = input->model.user1 / input->scaling.L;
  double Vx, Vz, P, eta, sxx, szz;
  if (position == W || position == E) {
    EvalDani(&Vx, &Vz, &P, &eta, &sxx, &szz, c.x, c.z, 1, radius, MM, MC);
    bc.type = 0; bc.value = Vx;
  } else if (position == S || position == SE || position == SW) {
    double x = c.x; double z = c.z + input->model.dx / 2.0;
    EvalDani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vx;
  } else if (position == N || position == NE || position == NW) {
    double x = c.x; double z = c.z - input->model.dx / 2.0;
    EvalDani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vx;
  } else { bc.type = -1; bc.value = 0.0; }
  return bc;
}

static SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates c) {
  SetBC bc;
  const double radius = input->model.user1 / input->scaling.L;
  double Vx, Vz, P, eta, sxx, szz;
  if (position == S || position == N) {
    EvalDani(&Vx, &Vz, &P, &eta, &sxx, &szz, c.x, c.z, 1, radius, MM, MC);
    bc.type = 0; bc.value = Vz;
  } else if (position == W || position == SW || position == NW) {
    double x = c.x + input->model.dx / 2.0; double z = c.z;
    EvalDani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vz;
  } else if (position == E || position == SE || position == NE) {
    double x = c.x - input->model.dx / 2.0; double z = c.z;
    EvalDani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, MM, MC);
    bc.type = 11; bc.value = Vz;
  } else { bc.type = -1; bc.value = 0.0; }
  return bc;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    std::fprintf(stderr,
                 "usage: %s <input.txt>\n"
                 "  Runs a single MDOODZ Stokes solve with SolVi-Dani BCs\n"
                 "  (MC=1000, MM=1). All MDOODZ logging goes to stdout; the\n"
                 "  caller is expected to redirect into a per-run log file.\n",
                 argv[0]);
    return 2;
  }

  MdoodzSetup setup = {};
  setup.SetParticles = new SetParticles_ff{
      .SetPhase   = SetPhase,
      .SetDensity = SetDensity,
  };
  setup.SetBCs = new SetBCs_ff{
      .SetBCVx = SetBCVx,
      .SetBCVz = SetBCVz,
  };

  char input_buf[512];
  std::strncpy(input_buf, argv[1], sizeof(input_buf) - 1);
  input_buf[sizeof(input_buf) - 1] = '\0';

  std::fprintf(stdout, "[GmgLevelsSweepProbe] input=%s\n", input_buf);
  std::fflush(stdout);

  RunMDOODZ(input_buf, &setup);
  return 0;
}
