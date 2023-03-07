#include "mdoodz.h"

static Ellipse GetWetOlivineEllipse(double centreX, double centreZ) {
  return (Ellipse) {
          .radiusZ = 3e3,
          .radiusX = 61e3,
          .centreX = centreX,
          .centreZ = centreZ,
          .angle   = 0,
  };
}

static Ellipse GetWetQuartziteEllipse(double centreX, double centreZ) {
  return (Ellipse) {
          .radiusZ = 4e3,
          .radiusX = 67e3,
          .centreX = centreX,
          .centreZ = centreZ,
          .angle   = 0,
  };
}

int SetPhase(MdoodzInput *input, Coordinates coordinates) {
  const double HCrust  = 30e3 / input->scaling.L;
  const double HMantle = 40e3 / input->scaling.L;
  // The continental crust is initially 30 km thick and is modeled with a visco-plastic layer of moderate strength (anorthite flow law; Rybacki and Dresen, 2004)
  if (coordinates.z > -HCrust) {
    if (IsEllipseCoordinates(coordinates, GetWetQuartziteEllipse(-50e3, -8e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetQuartziteEllipse(25e3, -7e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetQuartziteEllipse(-25e3, -18e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetQuartziteEllipse(48e3, -17e3), input->scaling.L)) {
      // and includes embedded elliptical bodies of relatively weaker (4 × 67 km ellipses; wet quartzite; Kirby, 1983)
      return 4;
    }
    // and stronger rheologies (4–7 × 42–77 km ellipses; Maryland diabase
    Ellipse diabase1 = (Ellipse){
            .angle   = 0,
            .radiusX = 42e3,
            .radiusZ = 6e3,
            .centreX = -16e3,
            .centreZ = -13e3,
    };
    Ellipse diabase2 = (Ellipse){
            .angle   = 0,
            .radiusX = 50e3,
            .radiusZ = 4e3,
            .centreX = -40e3,
            .centreZ = -25e3,
    };
    Ellipse diabase3 = (Ellipse){
            .angle   = 0,
            .radiusX = 77e3,
            .radiusZ = 7e3,
            .centreX = 35e3,
            .centreZ = -24e3,
    };
    if (IsEllipseCoordinates(coordinates, diabase1, input->scaling.L)
        || IsEllipseCoordinates(coordinates, diabase2, input->scaling.L)
        || IsEllipseCoordinates(coordinates, diabase3, input->scaling.L)) {
      return 6;
    }
    return 0;
  } else if (coordinates.z > -HCrust - HMantle) {
    // The upper subcontinental mantle is initially 40 km thick and is modeled using a visco-plastic layer dominated by dry olivine (Carter and Tsenn, 1987)
    // and incorporates weaker lens-shaped bodies dominated by wet olivine (3 × 61 km ellipses; Carter and Tsenn, 1987)
    if (IsEllipseCoordinates(coordinates, GetWetOlivineEllipse(-10e3, -37e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetOlivineEllipse(-55e3, -42e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetOlivineEllipse(32e3, -42e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetOlivineEllipse(-5e3, -47e3), input->scaling.L)
        || IsEllipseCoordinates(coordinates, GetWetOlivineEllipse(-28e3, -55e3), input->scaling.L)) {
      return 5;
    }
    return 1;
  } else {
    // The lower subcontinental mantle is modeled using a purely viscous layer dominated by a wet olivine rheology.
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
  if ( 1) {
    return input->materials.rho[phase] * (1 -  input->materials.alp[phase] * (TPart - input->materials.T0[phase]) );
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
