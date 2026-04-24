#ifndef MDOODZ_TESTS_SOLCX_ANALYTICAL_H
#define MDOODZ_TESTS_SOLCX_ANALYTICAL_H

// SolCx analytical solution for 2D Stokes flow with a step-function viscosity jump.
// Reference: Duretz, May, Gerya, Tackley (2011) "Discretization errors and free surface
// stabilization in the finite difference and marker-in-cell method for applied
// geodynamics", Geochemistry Geophysics Geosystems 12, Q07004.
//
// The underlying series evaluator `_Velic_solCx` is ported verbatim from ASPECT
// (benchmarks/solcx/solcx.h, GPL v2+), which in turn took it from the Underworld
// project (GPL). The function computes velocity, pressure, deviatoric stress and
// strain rate at any point (x, z) in the unit square.
//
// Domain: [0, 1] x [0, 1]
// Viscosity: eta_A for x < xc, eta_B for x >= xc (step function)
// Density forcing: rho(x, z) = -sin(n*pi*z) * cos(pi*x)  (n = wavenumber in z)
// Gravity: (0, 1)
// Boundary conditions: zero-normal-velocity, zero-tangential-stress (free-slip)
//   on all four boundaries. For MDOODZ test use we additionally impose the
//   analytical velocity as Dirichlet BCs to measure internal solver accuracy.
//
// Canonical headline parameters: eta_A=1, eta_B=1e6, xc=0.5, n=1.

namespace MdoodzSolCx {

struct SolCxValue {
  double vx;
  double vz;
  double p;
  double sxx;   // total stress xx
  double szz;   // total stress zz
  double sxz;   // total stress xz
  double exx;   // strain rate xx
  double ezz;   // strain rate zz
  double exz;   // strain rate xz
};

// Evaluate the Velic analytical solution at position (x, z).
//   eta_A  : viscosity for x < xc
//   eta_B  : viscosity for x >= xc
//   xc     : jump location (canonical 0.5)
//   n      : vertical density wavenumber (canonical 1)
SolCxValue EvalSolCx(double x, double z, double eta_A, double eta_B, double xc, int n);

}  // namespace MdoodzSolCx

#endif  // MDOODZ_TESTS_SOLCX_ANALYTICAL_H
