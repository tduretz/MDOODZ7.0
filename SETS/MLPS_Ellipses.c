#include "math.h"
#include "mdoodz.h"
#include "stdlib.h"


int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double HLit  = 70e3 / input->scaling.L;
  const double X0    = (75e3 * (double) rand() / RAND_MAX - 75e3 / 2) / input->scaling.L;
  const double Z0    = -HLit * (double) rand() / RAND_MAX;
  const double X     = coordinates.x;
  const double Z     = coordinates.z;
  const double a     = 1 / 0.005 / input->scaling.L;
  const double b     = 1 / 0.000016 / input->scaling.L;
  const double theta = 0 * M_PI / 180;
  if ((a * (X - X0) * (X - X0) + b * (Z - Z0) * (Z - Z0)) * pow(cos(theta), 2) + (b * (X - X0) * (X - X0) + a * (Z - Z0) * (Z - Z0)) * pow(sin(theta), 2) + (a - b) * (X - X0) * (Z - Z0) * sin(2 * theta) - 1 < 0) {
    if (Z0 > -HLit / 2) {
      return 1;
    }
    if (Z0 <= -HLit / 2) {
      return 3;
    }
  }
  if (coordinates.z < -HLit / 2) {
    return 2;
  }
  return 0;
}

double SetTemperature(MdoodzInput *input, Coordinates coordinates) {
  const double T0    = 273.15 / input->scaling.T;
  const double TmaxW = (1250 + 273.15) / input->scaling.T;
  const double HLit  = 70e3 / input->scaling.L;
  const double Htot  = 2 * HLit;
  return ((TmaxW - T0) / Htot) * (-coordinates.z - 0e3 / input->scaling.L) + T0;
}

double SetDensity(MdoodzInput *input, Coordinates coordinates, int phase) {
  const double TPart = SetTemperature(input, coordinates);
  if ( input->model.eqn_state ==  1 ) {
    return input->materials.rho[phase] * (1 -  input->materials.alp[phase] * (TPart - input->materials.T0[phase]) );
  } else {
    return input->materials.rho[phase];
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC  bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == FREE_SURFACE) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else {
    bc.value = 0.0;
    bc.type  = 0;
  }
  return bc;
}


SetBC SetBCTNew(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC  bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  double mantleTemperature  = (1330. + zeroC) / instance->scaling.T;
  if (position == S || position == SE || position == SW) {
    bc.value = particleTemperature;
    bc.type  = 1;
  } else if (position == N || position == NE || position == NW) {
    bc.value = surfaceTemperature;
    bc.type  = 1;
  } else if (position == W || position == E) {
    bc.value = mantleTemperature;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = 0;
  }
  return bc;
}


// B. Petri et al. / Earth and Planetary Science Letters 512 (2019) 147–162
int main() {
  MdoodzSetup setup = {
          .SetParticles = &(SetParticles_ff){
                  // The continental crust is initially 30 km thick and is modeled with a visco-plastic layer of moderate strength (anorthite flow law; Rybacki and Dresen, 2004)
                  // and includes embedded elliptical bod- ies of relatively weaker (4 × 67 km ellipses; wet quartzite; Kirby, 1983)
                  // and stronger rheologies (4–7 × 42–77 km ellipses; Maryland diabase; Mackwell et al., 1998)
                  // The upper subcontinental mantle is initially 40 km thick and is modeled using a visco-plastic layer dominated by dry olivine (Carter and Tsenn, 1987)
                  // and incorporates weaker lens-shaped bodies dominated by wet olivine (3 × 61 km ellipses; Carter and Tsenn, 1987)
                  // The lower subcontinental mantle is modeled using a purely viscous layer dominated by a wet olivine rheology.
                  .SetPhase       = SetPhase,
                  // The initial thermal field is that of an equilibrium field that includes radiogenic heat production
                  // in the continental crust and exhibits 500–550 ◦ C at the Moho.
                  .SetTemperature = SetTemperature,
                  .SetDensity     = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  // Extension is applied by prescribing horizontal velocities at the lateral model boundaries (left, right, bottom),
                  // satisfying a constant bulk extension rate of 10−15 s−1.
                  .SetBCVx    = SetPureShearBCVx,
                  .SetBCVz    = SetPureShearBCVz,
                  // the top and bottom boundaries are characterized by a constant temperature (0 and 1330◦C, respectively).
                  // Zero heat flow boundary conditions are prescribed at the lateral sides of the model;
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
  };
  RunMDOODZ("RiftingPauline.txt", &setup);
}
