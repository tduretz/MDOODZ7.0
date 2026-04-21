// add-gmg-upleg-fix §1.9 / §3.2 — 201² default-`gmg_levels` FGMRES
// convergence regression fixture.
//
// Spec contract (`gmg-stokes-solver`, scenario
// "201² default-levels convergence"):
//   - At 201² with `lin_solver = 3` and the default auto-selected
//     `gmg_levels` (≥ 4 on this grid), FGMRES SHALL reach
//     `gmg_fgmres_tol = 1e-6` in no more than 60 total inner iterations
//     (at most 2 restarts at default `gmg_fgmres_restart = 30`).
//   - The dispatch layer SHALL NOT fall back to CHOLMOD on this run
//     (no "falling back to CHOLMOD" LOG_WARN).
//
// How the test works
// ------------------
// MDOODZ logs via `fprintf(stdout, ...)` through `MdoodzLog.c`. We redirect
// stdout to a temporary file via freopen(3) before calling RunMDOODZ,
// then slurp the file back and grep for:
//   "GMG-FGMRES converged: iters=<N>, restarts=<R>, final_res_rel=<X>, tol=<T>"
// (the new log format landed in §1.6) and
//   "GMG lin_solver = 3.*falling back to CHOLMOD"
// (StokesRoutines.c fallback LOG_WARN).
//
// Pre-fix expectation: this test FAILS. The 201² stall reported in
// PERFORMANCE_REPORT §7.5 runs 600 inner iterations without reaching
// 1e-6 relative residual; the dispatch layer then falls back to CHOLMOD
// and emits the fallback LOG_WARN. All three assertions trip.
//
// Post-fix expectation (§5.2): FGMRES converges in ≤ 60 iterations
// (typical healthy count ≈ 20, with restart budget 60 giving comfortable
// margin), no fallback, clean convergence line.

extern "C" {
#include "mdoodz.h"
}

#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <regex>
#include <sstream>
#include <string>

#include <gtest/gtest.h>

// SolVi-Dani setup with 1000× inclusion, same shape as
// `VcycleAnimGenerate.cpp` / `SingleShearBandGmgEquivalence.cpp`. Inlined
// so this fixture has no cross-test coupling.

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

// -------- stdout capture via freopen ------------------------------------

static std::string RunAndCapture(const char *input_txt) {
  MdoodzSetup setup = {};
  setup.SetParticles = new SetParticles_ff{
      .SetPhase   = SetPhase,
      .SetDensity = SetDensity,
  };
  setup.SetBCs = new SetBCs_ff{
      .SetBCVx = SetBCVx,
      .SetBCVz = SetBCVz,
  };

  // freopen stdout → a temp file; restore via /dev/tty is brittle on
  // headless CI, so we dup the original fd and fdopen back after.
  const char *tmpdir_env = std::getenv("TMPDIR");
  std::string tmpdir = tmpdir_env ? tmpdir_env : "/tmp";
  std::string capture_path = tmpdir + "/gmg_converges_log.XXXXXX";
  std::vector<char> path_buf(capture_path.begin(), capture_path.end());
  path_buf.push_back('\0');
  int fd = mkstemp(path_buf.data());
  std::string real_path(path_buf.data());
  if (fd < 0) {
    ADD_FAILURE() << "mkstemp failed";
    return "";
  }
  close(fd);

  std::fflush(stdout);
  int saved_stdout = dup(fileno(stdout));
  FILE *redirected = std::freopen(real_path.c_str(), "w", stdout);
  EXPECT_NE(redirected, nullptr) << "freopen stdout redirect failed";

  char input_buf[256];
  std::strncpy(input_buf, input_txt, sizeof(input_buf) - 1);
  input_buf[sizeof(input_buf) - 1] = '\0';
  RunMDOODZ(input_buf, &setup);

  std::fflush(stdout);
  dup2(saved_stdout, fileno(stdout));
  close(saved_stdout);

  std::ifstream in(real_path);
  std::stringstream ss;
  ss << in.rdbuf();
  std::string captured = ss.str();
  std::remove(real_path.c_str());

  // Echo the captured log so test diagnostics still show it.
  std::fputs(captured.c_str(), stdout);
  return captured;
}

// -------- Test -----------------------------------------------------------

TEST(GmgConvergesAtDefaultLevels, Res201_DefaultLevels_WithinIterationBudget) {
  const std::string log = RunAndCapture(
      "SolViBenchmark/SolViGmgDefaultLevels201.txt");

  // (a) No fallback to CHOLMOD on this GMG run.
  EXPECT_EQ(log.find("falling back to CHOLMOD"), std::string::npos)
      << "Dispatch layer fell back to CHOLMOD — GMG did not converge at 201² "
      << "with default gmg_levels (add-gmg-upleg-fix spec: production-grid "
      << "convergence at default hierarchy depth)";

  // (b) A "GMG-FGMRES converged" line is present, with iters=<N> ≤ 60.
  std::regex rx(R"(GMG-FGMRES converged: iters=(\d+), restarts=(\d+), final_res_rel=([0-9.eE+\-]+), tol=([0-9.eE+\-]+))");
  std::smatch m;
  const bool have_converged = std::regex_search(log, m, rx);
  ASSERT_TRUE(have_converged)
      << "No \"GMG-FGMRES converged: ...\" line in the captured log — FGMRES "
      << "did not declare convergence (pre-fix stall expected here; this is "
      << "the §1.9 failing-fixture checkpoint).";

  const int iters    = std::stoi(m[1].str());
  const int restarts = std::stoi(m[2].str());
  const double rel   = std::stod(m[3].str());
  const double tol   = std::stod(m[4].str());

  std::printf("GmgConvergesAtDefaultLevels: iters=%d, restarts=%d, "
              "final_res_rel=%.3e, tol=%.3e\n",
              iters, restarts, rel, tol);

  // (c) final_res_rel ≤ tol must hold whenever "converged" is logged
  // (§1.6 log alignment contract).
  EXPECT_LE(rel, tol)
      << "\"converged\" log line reports final_res_rel=" << rel
      << " > tol=" << tol
      << " — predicate/log alignment violation (spec: convergence predicate "
      << "and log message SHALL report the same quantity).";

  // (d) Iteration-count budget: ≤ 60 inner iterations (≤ 2 restarts).
  EXPECT_LE(iters, 60)
      << "FGMRES converged but exceeded the 60-iteration budget (got " << iters
      << "); spec contract for 201² default-levels convergence requires "
      << "≤ 2 restarts at default gmg_fgmres_restart = 30.";
}
