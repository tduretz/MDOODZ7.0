// Reference 0D integrator for the combined tensile-cap / Drucker-Prager
// yield surface of Popov et al. (2025, GMD, doi:10.5194/gmd-18-7035-2025).
//
// This is a self-contained scalar implementation of the implicit local stress
// update described in Popov §3.6 and mirrored by MDLIB/RheologyDensity.c:685.
// It is NOT a link against MDLIB — it duplicates the algorithm on purpose so
// that the CI test can compare the full field solver's output against an
// independent implementation of the same equations.
//
// The analytic closed-form solution is only available for the trivial pure-
// loading cases (volumetric extension with psi=0, or pure shear with psi=0).
// For the dilatant and mixed cases we integrate the Perzyna-regularised
// constitutive equations with implicit backward-Euler steps matching the
// paper's time discretisation.
//
// Equation references in the comments below point to the paper.

#include "Popov0DAnalytical.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace popov2025 {

namespace {

constexpr int    kLocalNewtonIterMax = 60;
constexpr double kLocalNewtonAtol    = 1.0e-10;
constexpr double kLocalNewtonRtol    = 1.0e-8;

// Solve a 3x3 linear system by Cramer's rule.  Returns false if singular.
bool solve3x3(double a11, double a12, double a13,
              double a21, double a22, double a23,
              double a31, double a32, double a33,
              double b1,  double b2,  double b3,
              double &x1, double &x2, double &x3) {
  const double det = a11 * (a22 * a33 - a23 * a32)
                   - a12 * (a21 * a33 - a23 * a31)
                   + a13 * (a21 * a32 - a22 * a31);
  if (std::fabs(det) < 1e-300) return false;
  x1 = ((b1)   * (a22 * a33 - a23 * a32)
      - (a12)  * (b2  * a33 - a23 * b3)
      + (a13)  * (b2  * a32 - a22 * b3)) / det;
  x2 = ((a11)  * (b2  * a33 - a23 * b3)
      - (b1)   * (a21 * a33 - a23 * a31)
      + (a13)  * (a21 * b3  - b2  * a31)) / det;
  x3 = ((a11)  * (a22 * b3  - b2  * a32)
      - (a12)  * (a21 * b3  - b2  * a31)
      + (b1)   * (a21 * a32 - a22 * a31)) / det;
  return true;
}

}  // namespace

YieldGeom geometry(const Params0D &p) {
  YieldGeom g;
  g.k     = std::sin(p.phi_rad);
  g.kq    = std::sin(p.psi_rad);
  g.c     = p.C_mc * std::cos(p.phi_rad);
  g.a     = std::sqrt(1.0 + g.k * g.k);
  g.b     = std::sqrt(1.0 + g.kq * g.kq);
  g.py    = (p.T_st + g.c / g.a) / (1.0 - g.k / g.a);  // Eq. 15
  g.Ry    = g.py - p.T_st;                              // Eq. 15
  g.pd    = g.py - g.Ry * g.k / g.a;                    // Eq. 16
  g.tau_d = g.k * g.pd + g.c;                           // Eq. 16
  return g;
}

std::vector<double> time_samples(const Params0D &p) {
  const std::size_t N = static_cast<std::size_t>(std::ceil(p.t_end / p.dt)) + 1;
  std::vector<double> t(N);
  for (std::size_t i = 0; i < N; ++i) t[i] = static_cast<double>(i) * p.dt;
  return t;
}

// Integrate one implicit time step.  Input: previous (tauII_n, p_n).  Output:
// updated (tauII, p).  Follows the residual form in Popov §3.6 / Eqs. 42-47
// and MDLIB/RheologyDensity.c:727.
static void stepImplicit(const Params0D &pr, const YieldGeom &g,
                         double tau_n, double p_n,
                         double &tau_out, double &p_out) {
  // Trial visco-elastic stress (no viscous creep in our 0D test setup —
  // A_D = A_N = 0, see eta0 = 1e50 + cstv = 1 in the parameter files).
  double tau_trial = tau_n + 2.0 * pr.G * pr.edot_dev * pr.dt;
  double p_trial   = p_n   - pr.K *       pr.edot_vol * pr.dt;
  if (tau_trial < 0.0) tau_trial = 0.0;

  // Pick active segment using the delimiter condition (paper Eq. 18).
  auto onDPyield = [&](double tau, double p) {
    return tau * (g.py - g.pd) >= g.tau_d * (g.py - p);
  };
  auto onDPflow = [&](double tau, double p) {
    const double pq = g.pd + g.kq * g.tau_d;
    return tau * (pq - g.pd) >= g.tau_d * (pq - p);
  };

  double F_trial;
  {
    const double R_haty = std::sqrt(tau_trial * tau_trial
                                    + (p_trial - g.py) * (p_trial - g.py));
    F_trial = onDPyield(tau_trial, p_trial)
                ? (tau_trial - g.k * p_trial - g.c)
                : (g.a * (R_haty - g.Ry));
  }

  if (F_trial <= 0.0) {
    tau_out = tau_trial;
    p_out   = p_trial;
    return;
  }

  // Corrector — 3x3 implicit Newton on (tau, p, lam_dot).
  double tau = tau_trial, p = p_trial, lam = 0.0;
  const double eta_ve = pr.G * pr.dt;  // With A_D = A_N = 0, A_L = 1/(2 G dt),
                                       // hence eta_ve = G * dt (Eq. A5).

  double R_prev = 2.0 * kLocalNewtonAtol;
  for (int it = 0; it < kLocalNewtonIterMax; ++it) {
    if (tau < 0.0) tau = 0.0;
    const double pq      = g.pd + g.kq * g.tau_d;
    const double R_haty  = std::sqrt(tau * tau + (p - g.py) * (p - g.py));
    const double R_hatq  = std::sqrt(tau * tau + (p - pq) * (p - pq));

    const bool DP_yield = onDPyield(tau, p);
    const bool DP_flow  = onDPflow(tau, p);

    double F, dFdtau, dFdp;
    if (DP_yield) {
      F      = tau - g.k * p - g.c;
      dFdtau = 1.0;
      dFdp   = -g.k;
    } else {
      F      = g.a * (R_haty - g.Ry);
      dFdtau = g.a * tau          / R_haty;
      dFdp   = g.a * (p - g.py)   / R_haty;
    }

    double A_tau, A_p, A_tt, A_pp, A_tp, A_pt;
    if (DP_flow) {
      A_tau = 0.5;  A_p = g.kq;
      A_tt  = 0.0;  A_pp = 0.0;  A_tp = 0.0;  A_pt = 0.0;
    } else {
      A_tau = g.b * tau       / (2.0 * R_hatq);
      A_p   = -g.b * (p - pq) / R_hatq;
      A_tt  =  g.b * (R_hatq * R_hatq - tau * tau) / (2.0 * R_hatq * R_hatq * R_hatq);
      A_pp  = -g.b * (R_hatq * R_hatq - (p - pq) * (p - pq)) / (R_hatq * R_hatq * R_hatq);
      A_tp  = -g.b * tau * (p - pq) / (2.0 * R_hatq * R_hatq * R_hatq);
      A_pt  =  g.b * tau * (p - pq) / (R_hatq * R_hatq * R_hatq);
    }

    // Residuals (Eq. 42)
    const double R1 = (tau_trial - tau) / (2.0 * eta_ve) - lam * A_tau;
    const double R2 = (p - p_trial) / (pr.K * pr.dt)     - lam * A_p;
    const double R3 = F                                  - lam * pr.eta_vp;

    const double R_norm = std::sqrt(R1 * R1 + R2 * R2 + R3 * R3);
    if (R_norm < kLocalNewtonAtol || R_norm / R_prev < kLocalNewtonRtol) break;
    R_prev = std::max(R_norm, kLocalNewtonAtol);

    // Jacobian (Eq. 47)
    const double J11 = -1.0 / (2.0 * eta_ve) - lam * A_tt;
    const double J12 =                         -lam * A_tp;
    const double J13 = -A_tau;
    const double J21 =                         -lam * A_pt;
    const double J22 =  1.0 / (pr.K * pr.dt)  - lam * A_pp;
    const double J23 = -A_p;
    const double J31 =  dFdtau;
    const double J32 =  dFdp;
    const double J33 = -pr.eta_vp;

    double dtau, dp, dlam;
    // Newton update: J * delta = -R, so pass (-R1, -R2, -R3).
    if (!solve3x3(J11, J12, J13, J21, J22, J23, J31, J32, J33,
                  -R1, -R2, -R3, dtau, dp, dlam)) break;

    // Damped Newton step (Armijo-style line search)
    double alpha = 1.0;
    for (int ls = 0; ls < 20; ++ls) {
      const double tau_t = tau + alpha * dtau;
      const double p_t   = p   + alpha * dp;
      const double lam_t = lam + alpha * dlam;
      const double R_haty_t = std::sqrt(tau_t * tau_t + (p_t - g.py) * (p_t - g.py));
      const double F_t = onDPyield(tau_t, p_t)
                           ? (tau_t - g.k * p_t - g.c)
                           : (g.a * (R_haty_t - g.Ry));
      double A_tau_t, A_p_t;
      if (onDPflow(tau_t, p_t)) {
        A_tau_t = 0.5; A_p_t = g.kq;
      } else {
        const double R_hatq_t = std::sqrt(tau_t * tau_t + (p_t - pq) * (p_t - pq));
        A_tau_t = g.b * tau_t       / (2.0 * R_hatq_t);
        A_p_t   = -g.b * (p_t - pq) / R_hatq_t;
      }
      const double Rt1 = (tau_trial - tau_t) / (2.0 * eta_ve) - lam_t * A_tau_t;
      const double Rt2 = (p_t - p_trial) / (pr.K * pr.dt)     - lam_t * A_p_t;
      const double Rt3 = F_t                                  - lam_t * pr.eta_vp;
      const double Rt_norm = std::sqrt(Rt1 * Rt1 + Rt2 * Rt2 + Rt3 * Rt3);
      if (Rt_norm < R_norm * (1.0 - 1e-4 * alpha)) break;
      alpha *= 0.5;
      if (alpha < 1e-4) break;
    }
    tau += alpha * dtau;
    p   += alpha * dp;
    lam += alpha * dlam;
  }

  tau_out = tau;
  p_out   = p;
}

static void integrateTrajectory(const Params0D &p,
                                std::vector<double> &tauII,
                                std::vector<double> &P) {
  const YieldGeom g = geometry(p);
  const std::vector<double> t = time_samples(p);
  tauII.assign(t.size(), 0.0);
  P.assign(t.size(), 0.0);
  double tau_n = 0.0, p_n = 0.0;
  tauII[0] = tau_n; P[0] = p_n;
  for (std::size_t i = 1; i < t.size(); ++i) {
    double tau_out, p_out;
    stepImplicit(p, g, tau_n, p_n, tau_out, p_out);
    tau_n    = tau_out;
    p_n      = p_out;
    tauII[i] = tau_out;
    P[i]     = p_out;
  }
}

std::vector<double> tauII_trajectory(const Params0D &p) {
  std::vector<double> tau, pp;
  integrateTrajectory(p, tau, pp);
  return tau;
}

std::vector<double> P_trajectory(const Params0D &p) {
  std::vector<double> tau, pp;
  integrateTrajectory(p, tau, pp);
  return pp;
}

double tauII_at(double t, const Params0D &p_in) {
  Params0D p = p_in;
  p.t_end    = t;
  const auto traj = tauII_trajectory(p);
  return traj.empty() ? 0.0 : traj.back();
}

double P_at(double t, const Params0D &p_in) {
  Params0D p = p_in;
  p.t_end    = t;
  const auto traj = P_trajectory(p);
  return traj.empty() ? 0.0 : traj.back();
}

}  // namespace popov2025

// -----------------------------------------------------------------------------
// Standalone sanity unit test (compiled only with -DPOPOV0D_STANDALONE).
// Verifies the helper against two Popov Fig. 5 checkpoints.
// -----------------------------------------------------------------------------
#ifdef POPOV0D_STANDALONE
#include <cassert>
#include <cstdio>

int main() {
  using namespace popov2025;

  // Volumetric extension steady state (paper Fig. 5a): after ~15 yr, p
  // saturates at T_st and tauII stays at 0.
  {
    Params0D p{};
    p.G        = 1.0e10;
    p.K        = 2.0e11;
    p.phi_rad  = 30.0 * M_PI / 180.0;
    p.psi_rad  = 10.0 * M_PI / 180.0;
    p.C_mc     = 1.0e6;
    p.T_st     = -5.0e5;
    p.eta_vp   = 0.0;
    p.edot_dev = 0.0;
    p.edot_vol = 7.0e-15;
    p.dt       = 2.0 * 365.25 * 86400.0;
    p.t_end    = 20.0 * 365.25 * 86400.0;
    std::vector<double> tau, pp;
    integrateTrajectory(p, tau, pp);
    std::printf("VolumetricExtension steady: tauII = %.3e Pa, p = %.3e Pa (expect ~0, ~%.3e)\n",
                tau.back(), pp.back(), p.T_st);
    assert(std::fabs(tau.back()) < 1.0e4);
    assert(std::fabs(pp.back() - p.T_st) / std::fabs(p.T_st) < 0.05);
  }

  // Deviatoric shear steady state (paper Fig. 5b): after ~25 yr, tauII
  // saturates near C*cos(phi) (minus small dilation coupling).
  {
    Params0D p{};
    p.G        = 1.0e10;
    p.K        = 2.0e11;
    p.phi_rad  = 30.0 * M_PI / 180.0;
    p.psi_rad  = 10.0 * M_PI / 180.0;
    p.C_mc     = 1.0e6;
    p.T_st     = -5.0e5;
    p.eta_vp   = 0.0;
    p.edot_dev = 7.0e-14;
    p.edot_vol = 0.0;
    p.dt       = 2.0 * 365.25 * 86400.0;
    p.t_end    = 30.0 * 365.25 * 86400.0;
    std::vector<double> tau, pp;
    integrateTrajectory(p, tau, pp);
    const YieldGeom g = geometry(p);
    const double tau_dp_at_p = g.k * pp.back() + g.c;
    std::printf("DeviatoricShear steady: tauII = %.3e Pa, p = %.3e Pa, DP(p) = %.3e Pa\n",
                tau.back(), pp.back(), tau_dp_at_p);
    assert(std::fabs(tau.back() - tau_dp_at_p) / tau_dp_at_p < 0.05);
  }

  std::printf("Popov0DAnalytical standalone checks passed.\n");
  return 0;
}
#endif  // POPOV0D_STANDALONE
