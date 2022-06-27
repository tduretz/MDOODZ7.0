#include "math.h"
#include "mdoodz.h"
#include "stdlib.h"


int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  // The continental crust is initially 30 km thick and is modeled with a visco-plastic layer of moderate strength (anorthite flow law; Rybacki and Dresen, 2004)
  // and includes embedded elliptical bodies of relatively weaker (4 × 67 km ellipses; wet quartzite; Kirby, 1983)
  // and stronger rheologies (4–7 × 42–77 km ellipses; Maryland diabase; Mackwell et al., 1998)
  // The upper subcontinental mantle is initially 40 km thick and is modeled using a visco-plastic layer dominated by dry olivine (Carter and Tsenn, 1987)
  // and incorporates weaker lens-shaped bodies dominated by wet olivine (3 × 61 km ellipses; Carter and Tsenn, 1987)
  // The lower subcontinental mantle is modeled using a purely viscous layer dominated by a wet olivine rheology.
  const double HCrust  = 30e3 / input->scaling.L;
  const double HMantle  = 40e3 / input->scaling.L;
  if (coordinates.z > -HCrust) {
    return 0;
  } else if (coordinates.z > -HCrust - HMantle) {
    return 1;
  } else {
    return 2;
  }
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
          .BuildInitialTopography = &(BuildInitialTopography_ff){

          },
          .SetParticles = &(SetParticles_ff){
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
  RunMDOODZ("MLPS_Ellipses.txt", &setup);
}
