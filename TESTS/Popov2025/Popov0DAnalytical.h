#ifndef POPOV_2025_0D_ANALYTICAL_H
#define POPOV_2025_0D_ANALYTICAL_H

// Reference solutions for the three 0D stress-integration tests of Popov et al.
// (2025, GMD, doi:10.5194/gmd-18-7035-2025), Fig. 5 a/b/c.  These are
// closed-form piecewise expressions that integrate the Maxwell visco-elastic
// trial stress until the combined tensile-cap / Drucker-Prager yield surface
// is reached, then the Perzyna ODE along the active segment.  All arguments
// and return values are in SI units (Pa, s, 1/s).  Equation references below
// point to the paper.

#include <cstddef>
#include <vector>

namespace popov2025 {

// Parameters of the combined yield surface and loading path.
struct Params0D {
  // Elastic moduli
  double G;        // shear modulus [Pa]
  double K;        // bulk modulus [Pa]

  // Drucker-Prager shear envelope (paper Eq. 13)
  double phi_rad;  // friction angle [rad]
  double psi_rad;  // dilation angle [rad]
  double C_mc;     // Mohr-Coulomb cohesion [Pa]
  double T_st;     // tensile strength [Pa] — negative in the paper's convention

  // Perzyna regularization (paper Eq. 11)
  double eta_vp;   // viscoplastic regularization viscosity [Pa s]

  // 0D loading: constant background strain rates applied from t = 0.
  double edot_dev; // deviatoric second invariant of the applied strain rate [1/s]
  double edot_vol; // volumetric strain rate (trace) [1/s]  — positive = extension

  // Time discretization used by MDOODZ.
  double dt;       // step size [s]
  double t_end;    // total simulated time [s]
};

// Derived yield-surface geometry (paper §2.4, Eqs. 13-17).
struct YieldGeom {
  double k;    // sin(phi)
  double kq;   // sin(psi)
  double c;    // C_mc * cos(phi)
  double a;    // sqrt(1 + k^2)
  double b;    // sqrt(1 + kq^2)
  double py;   // (T_st + c/a) / (1 - k/a) — centre of yield cap
  double Ry;   // py - T_st                — radius of yield cap
  double pd;   // py - Ry * k/a            — delimiter-point pressure
  double tau_d;// k * pd + c               — delimiter-point stress
};

YieldGeom geometry(const Params0D &p);

// Analytical trajectories.  Each returns a vector of length N =
// ceil(t_end / dt) + 1 sampled at t_i = i * dt.
std::vector<double> time_samples(const Params0D &p);
std::vector<double> tauII_trajectory(const Params0D &p);
std::vector<double> P_trajectory(const Params0D &p);

// Convenience: point evaluation at time t (uses the same piecewise formulas
// as the trajectory functions).
double tauII_at(double t, const Params0D &p);
double P_at(double t, const Params0D &p);

}  // namespace popov2025

#endif  // POPOV_2025_0D_ANALYTICAL_H
