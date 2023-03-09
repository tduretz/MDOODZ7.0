#include "mdoodz.h"

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  Rectangle easternContinent = {
          .sizeZ   = 200e3,
          .sizeX   = 1500e3,
          .centreZ = -100e3,
          .centreX = 1250e3,
          .angle   = 0,
  };
  Rectangle westernContinent = {
          .sizeZ   = 200e3,
          .sizeX   = 1500e3,
          .centreZ = -100e3,
          .centreX = -1250e3,
          .angle   = 0,
  };
  const double HOceanicPlate = 20e3 / input->scaling.L;
  if (IsRectangleCoordinates(coordinates, easternContinent, input->scaling.L) || IsRectangleCoordinates(coordinates, westernContinent, input->scaling.L)) {
    return 2;
  } else if (coordinates.z > -HOceanicPlate) {
    return 1;
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
  /* Hi Roman, there is an awkward 1==0 below. In fact the global switch `eqn_state` was never used.
  Density is computed interally and the type of equation of state depends on the phase.
  To update density on particles one may use the function:
  EvaluateDensity( phase_ID, T, P, X,  model, materials );  where X is likely 0.0 in most cases (it's a depletion amount)
  --> Consequence: SetDensity() should be aware `of materials`
  */
  if ( 1==0 ) {
    return input->materials.rho[phase] * (1 - input->materials.alp[phase] * (TPart - input->materials.T0[phase]));
  } else {
    return input->materials.rho[phase];
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, double particleTemperature) {
  SetBC  bc;
  double surfaceTemperature = zeroC / instance->scaling.T;
  if (position == free_surface) {
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
                  .SetBCVx   = SetPureShearBCVx,
                  .SetBCVz   = SetPureShearBCVz,
                  // the top and bottom boundaries are characterized by a constant temperature (0 and 1330◦C, respectively).
                  // Zero heat flow boundary conditions are prescribed at the lateral sides of the model;
                  .SetBCT    = SetBCT,
                  .SetBCTNew = SetBCTNew,
          },
  };
  RunMDOODZ("CollisionZone.txt", &setup);
}
