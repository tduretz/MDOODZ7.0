#include "math.h"
#include "mdoodz.h"

// ============================================================================
// BlankenBench — Blankenbach et al. (1989) Case 1a thermal convection
//
// Isoviscous Rayleigh-Bénard convection in a unit square box.
// Ra = α·ρ·g·ΔT·L³ / (κ·η) = 10⁴
//
// BCs:  Free-slip on all walls
//       T = T0 at top (Dirichlet), T = T0 + 1 at bottom (Dirichlet)
//       dT/dx = 0 on left/right (Neumann)
//
// Parameters:  η=1, ρ=1, α=1, g=Ra=1e4, ΔT=1, L=1, k=1, Cp=1 → κ=1
// T0 = zeroC = 273.15 K (MDOODZ convention for reference temperature)
// ============================================================================

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  return 0;
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  // Linear conductive profile: T = T_bot at bottom, T_top at top
  // T(z) = T_top + (T_bot - T_top) * (z_top - z) / H
  // Plus cosine perturbation to seed single convection cell
  const double T_top = (input->model.user0 + zeroC) / input->scaling.T;
  const double T_bot = input->model.user1 / input->scaling.T;
  const double zmin  = input->model.zmin;
  const double zmax  = input->model.zmax;
  const double H     = zmax - zmin;
  const double z_nd  = (zmax - coordinates.z) / H;            // 0 at top, 1 at bottom

  double T_lin = T_top + (T_bot - T_top) * z_nd;

  // Cosine perturbation to break symmetry and seed convection
  const double xmin = input->model.xmin;
  const double xmax = input->model.xmax;
  const double Lx   = xmax - xmin;
  double x_nd = (coordinates.x - xmin) / Lx;
  double perturbation = 0.01 * cos(M_PI * x_nd) * sin(M_PI * z_nd) / input->scaling.T;

  return T_lin + perturbation;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  return input->materials.rho[phase];
}

SetBC SetBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == S || position == SW || position == SE) {
    bc.value = 0.0;
    bc.type  = 11;     // free slip bottom
  } else if (position == N || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;     // free slip top
  } else if (position == W || position == E) {
    bc.value = 0.0;
    bc.type  = 0;      // free slip sides (Vx = 0)
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0.0;
    bc.type  = 13;     // free slip sides
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type  = 0;      // no-penetration top/bottom
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

char SetBCPType(MdoodzInput *input, POSITION position) {
  if (position == NE || position == NW) {
    return 0;          // pressure pin
  } else {
    return -1;
  }
}

SetBC SetBCT(MdoodzInput *input, POSITION position, Coordinates coordinates, double gridTemperature) {
  SetBC bc;
  double T_top = (input->model.user0 + zeroC) / input->scaling.T;
  double T_bot = input->model.user1 / input->scaling.T;
  if (position == N || position == NE || position == NW) {
    bc.value = T_top;
    bc.type  = 1;      // Dirichlet
  } else if (position == S || position == SE || position == SW) {
    bc.value = T_bot;
    bc.type  = 1;      // Dirichlet
  } else {
    bc.value = 0.0;
    bc.type  = 0;      // Neumann (zero heat flux)
  }
  return bc;
}

int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase       = SetPhase,
                  .SetTemperature = SetTemperature,
                  .SetDensity     = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx,
                  .SetBCVz    = SetBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
          },
  };
  RunMDOODZ("BlankenBench.txt", &setup);
}
