/* SolViPerf — parameterised SolVi driver for the local perf harness.
 *
 * Structurally identical to ViscousInclusion.c (same eval_anal_Dani BCs,
 * SetPhase, SetDensity) but main() takes argv[1] as the .txt fixture path
 * so the harness can hot-swap Nx / Nz / lin_solver overrides without
 * rebuilding MDOODZ between runs. On exit prints a single line beginning
 * "PERF_JSON:" carrying wall time (s) and peak RSS (KiB); callers grep
 * for that sentinel and discard everything else in the MDOODZ log.
 *
 * Added under add-gmg-stokes-defence §5 (macOS-local sweep replacing the
 * original EC2 design; see design.md §D1 amendment).
 */

#include "complex.h"
#include "mdoodz.h"
#include "math.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/resource.h>

/*----- eval_anal_Dani: verbatim from ViscousInclusion.c -----*/
void eval_anal_Dani(double *vx, double *vz, double *p, double *eta, double *sxx, double *syy, double x, double z, int ps, double rc, double mm, double mc) {
  double gr, er, A;
  _Complex double V_tot, phi_z, d_phi_z, conj_d_phi_z, psi_z, conj_psi_z, Z, d_d_phi_z_z, d_psi_z;

  *eta = *vx = *vz = *p = 0;

  if (ps == 1) { gr = 0; er = -1; }
  else         { gr = -2.0; er = 0; }
  A = mm * (mc - mm) / (mc + mm);

  Z = x + I * z;
  if (sqrt(pow(x, 2) + pow(z, 2)) <= rc) {
    *p   = 0.0;
    V_tot = (mm / (mc + mm)) * (I * gr + 2 * er) * conj(Z) - (I / 2) * gr * Z;
    *vx   = creal(V_tot);
    *vz   = cimag(V_tot);
    *eta  = mc;
    *syy  = 3.5;
    *sxx  = -3.5;
  } else {
    phi_z        = -(I / 2) * mm * gr * Z - I * gr * A * pow(rc, 2) / Z;
    d_phi_z      = -(I / 2) * mm * gr + I * gr * A * pow(rc, 2) / pow(Z, 2);
    conj_d_phi_z = conj(d_phi_z);
    psi_z        = (I * gr - 2 * er) * mm * Z - (I * gr + 2 * er) * A * pow(rc, 4) * cpow(Z, -3);
    conj_psi_z   = conj(psi_z);
    V_tot        = (phi_z - (x - I * z) * conj_d_phi_z - conj_psi_z) / (2 * mm);
    *p           = -2.0 * mm * creal(d_phi_z);
    *vx          = creal(V_tot);
    *vz          = cimag(V_tot);
    *eta         = mm;
    d_d_phi_z_z  = -(I * gr + 2 * er) * A * pow(rc, 2) / cpow(Z, 3);
    d_psi_z      = (I * gr - 2 * er) * mm + (I * gr + 2 * er) * A * pow(rc, 4) * cpow(Z, -4);
    *syy         = 2 * creal(d_phi_z) + creal((x - I * z) * d_d_phi_z_z + d_psi_z);
    *sxx         = 4 * creal(d_phi_z) - *syy;
  }
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) return 1;
  return 0;
}

double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  return instance->materials.rho[phase];
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm = 1.0, mc = 1e3;
  double Vx, Vz, P, eta, sxx, szz, x, z;
  if (position == W) {
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, mm, mc);
    bc.type = 2; bc.value = sxx;
  } else if (position == E) {
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, mm, mc);
    bc.type = 0; bc.value = Vx;
  } else if (position == S || position == SE || position == SW) {
    x = coord.x; z = coord.z + instance->model.dx / 2.0;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
    bc.type = 11; bc.value = Vx;
  } else if (position == N || position == NE || position == NW) {
    x = coord.x; z = coord.z - instance->model.dx / 2.0;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
    bc.type = 11; bc.value = Vx;
  } else { bc.type = -1; bc.value = 0.0; }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm = 1.0, mc = 1e3;
  double Vx, Vz, P, eta, sxx, szz, x, z;
  if (position == S) {
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, mm, mc);
    bc.type = 0; bc.value = Vz;
  } else if (position == N) {
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, coord.x, coord.z, 1, radius, mm, mc);
    bc.type = 2; bc.value = szz;
  } else if (position == W || position == SW || position == NW) {
    x = coord.x + instance->model.dx / 2.0; z = coord.z;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
    bc.type = 11; bc.value = Vz;
  } else if (position == E || position == SE || position == NE) {
    x = coord.x - instance->model.dx / 2.0; z = coord.z;
    eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
    bc.type = 11; bc.value = Vz;
  } else { bc.type = -1; bc.value = 0.0; }
  return bc;
}

int main(int argc, char **argv) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <path-to-.txt>\n", argv[0]);
    return 2;
  }
  const char *txt = argv[1];

  MdoodzSetup instance = {
    .SetParticles = &(SetParticles_ff){ .SetPhase = SetPhase, .SetDensity = SetDensity },
    .SetBCs       = &(SetBCs_ff){       .SetBCVx  = SetBCVx,  .SetBCVz    = SetBCVz   },
  };

  struct timespec t0, t1;
  clock_gettime(CLOCK_MONOTONIC, &t0);
  RunMDOODZ((char*)txt, &instance);
  clock_gettime(CLOCK_MONOTONIC, &t1);

  double wall_s = (t1.tv_sec - t0.tv_sec) + 1e-9 * (t1.tv_nsec - t0.tv_nsec);

  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  /* On macOS (Darwin) ru_maxrss is in bytes; on Linux it's already in KiB.
   * Normalise to KiB for the output. */
  long peak_rss_kb;
#ifdef __APPLE__
  peak_rss_kb = (long)(ru.ru_maxrss / 1024L);
#else
  peak_rss_kb = (long)ru.ru_maxrss;
#endif

  /* Sentinel + machine-readable line; the harness greps for "PERF_JSON:". */
  printf("\nPERF_JSON:{\"wall_time_s\":%.6f,\"peak_rss_kb\":%ld,\"txt\":\"%s\"}\n",
         wall_s, peak_rss_kb, txt);
  fflush(stdout);
  return 0;
}
