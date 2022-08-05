#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


double SetSurfaceZCoord(MdoodzInput *instance, double x_coord) {
  const double TopoLevel = -0.0e3 / instance->scaling.L;
  const double h_pert    = instance->model.user3 / instance->scaling.L;
  return TopoLevel + h_pert * (3330.0 - 2800.0) / 2800.0 * cos(2 * M_PI * x_coord / (instance->model.xmax - instance->model.xmin));
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness  = instance->model.user1 / instance->scaling.L;
  const double crustThickness        = instance->model.user2 / instance->scaling.L;
  const double perturbationAmplitude = instance->model.user3 / instance->scaling.L;
  const double mohoLevel             = -crustThickness - perturbationAmplitude * cos(2 * M_PI * coordinates.x / (instance->model.xmax - instance->model.xmin));
  const bool   isBelowLithosphere    = coordinates.z < -lithosphereThickness;
  const bool   isAboveMoho           = coordinates.z > mohoLevel;

  if (instance->model.user4 && isAboveMoho) {
    const bool is2500MAboveMoho = coordinates.z > mohoLevel + 2500 / instance->scaling.L;
    const bool is4500MAboveMoho = coordinates.z > mohoLevel + 4500 / instance->scaling.L;
    const bool is7000MAboveMoho = coordinates.z > mohoLevel + 7000 / instance->scaling.L;
    const bool is9000MAboveMoho = coordinates.z > mohoLevel + 9000 / instance->scaling.L;
    if (is2500MAboveMoho && !is4500MAboveMoho) {
      return 7;
    } else if (is7000MAboveMoho && !is9000MAboveMoho) {
      return 8;
    } else {
      return 1;
    }
  } else if (isAboveMoho) {
    return 1;
  } else if (isBelowLithosphere) {
    return 3;
  } else {
    return 2;
  }
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double lithosphereThickness = instance->model.user1 / instance->scaling.L;
  const double surfaceTemperature   = 273.15 / instance->scaling.T;
  const double mantleTemperature    = (1330.0 + 273.15) / instance->scaling.T;
  const double particleTemperature  = ((mantleTemperature - surfaceTemperature) / lithosphereThickness) * (-coordinates.z) + surfaceTemperature;
  if (particleTemperature > mantleTemperature) {
    return mantleTemperature;
  } else {
    return particleTemperature;
  }
}

double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
  const int asthenospherePhase = 3;
  return instance->materials.gs_ref[asthenospherePhase];
}

double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return -coordinates.x * instance->model.EpsBG;
}

double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
  return coordinates.z * instance->model.EpsBG;
}

char SetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
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

void AddCrazyConductivity(MdoodzInput *input) {
  int               *asthenospherePhases = (int *) malloc(sizeof(int));
  CrazyConductivity *crazyConductivity   = (CrazyConductivity *) malloc(sizeof(CrazyConductivity));
  asthenospherePhases[0]                 = 3;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000;
  input->crazyConductivity               = crazyConductivity;
}

void AddAnisotropy(MdoodzInput *input, MutateInputParams *mutateInputParams) {
  AddCrazyConductivity(input);
  input->model.aniso                              = 1;
  input->model.fstrain                            = 1;
  const int crustalPhase                          = 1;
  input->materials.cstv[crustalPhase]             = mutateInputParams->param1;
  if (input->materials.cstv[crustalPhase]) {
    input->materials.pwlv[crustalPhase]           = 0;
  } else {
    input->materials.pwlv[crustalPhase]           = 11;
  }
  input->materials.aniso_factor[crustalPhase]     = mutateInputParams->param3;
  input->materials.aniso_angle[crustalPhase]      = mutateInputParams->param4;
  const int upperMantlePhase                      = 2;
  input->materials.cstv[upperMantlePhase]         = mutateInputParams->param5;
  if (input->materials.cstv[upperMantlePhase]) {
    input->materials.pwlv[upperMantlePhase]           = 0;
  } else {
    input->materials.pwlv[upperMantlePhase]           = 11;
  }
  input->materials.aniso_factor[upperMantlePhase] = mutateInputParams->param7;
  input->materials.aniso_angle[upperMantlePhase]  = mutateInputParams->param8;
  
  snprintf(input->model.description, sizeof(input->model.description), "Nt: %i CrustalPhase: {cstv: %i, pwlv: %i, aniso_factor: %f, aniso_angle: %f}, UpperMantlePhase: {cstv: %i, pwlv: %i, aniso_factor: %f, aniso_angle: %f}}", input->model.Nt, mutateInputParams->param1, mutateInputParams->param2, mutateInputParams->param3, mutateInputParams->param4, mutateInputParams->param5, mutateInputParams->param6, mutateInputParams->param7, mutateInputParams->param8);
}

int main() {
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
                  .SetGrainSize          = SetGrainSize,
                  .SetHorizontalVelocity = SetHorizontalVelocity,
                  .SetVerticalVelocity   = SetVerticalVelocity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetPureShearBCVx,
                  .SetBCVz    = SetPureShearBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
                  .SetBCTNew  = SetBCTNew,
          },
          .MutateInput = AddCrazyConductivity,
  };
  //RunMDOODZ("RiftingPaulineAniso.txt", &setup);
  //rename("Output00100.gzip.h5", "RiftingPaulineIsotropic.h5");

  MutateInputParams *mutateInputParams     = (MutateInputParams *) malloc(sizeof(MutateInputParams));

  int                crustCstvs[2]         = {1, 0};
  double             crustAnisoFactors[2]  = {2.0, 4.0};
  double             crustAnisoAngles[2]   = {25.0, 45.0};
  int                mantleCstvs[2]        = {1, 0};
  double             mantleAnisoFactors[2] = {2.0, 4.0};
  double             mantleAnisoAngles[2]  = {25.0, 45.0};

  int                i                     = 0;
  for (int crustCstv = 0; crustCstv < 2; crustCstv++) {
    for (int crustAnisoFactor = 0; crustAnisoFactor < 2; crustAnisoFactor++) {
      for (int crustAnisoAngle = 0; crustAnisoAngle < 2; crustAnisoAngle++) {
        for (int mantleCstv = 0; mantleCstv < 2; mantleCstv++) {
          for (int mantleAnisoFactor = 0; mantleAnisoFactor < 2; mantleAnisoFactor++) {
            for (int mantleAnisoAngle = 0; mantleAnisoAngle < 2; mantleAnisoAngle++) {
              mutateInputParams->param1 = crustCstvs[crustCstv];
              mutateInputParams->param3 = crustAnisoFactors[crustAnisoFactor];
              mutateInputParams->param4 = crustAnisoAngles[crustAnisoAngle];
              mutateInputParams->param5 = mantleCstvs[mantleCstv];
              mutateInputParams->param7 = mantleAnisoFactors[mantleAnisoFactor];
              mutateInputParams->param8 = mantleAnisoAngles[mantleAnisoAngle];
              setup.mutateInputParams   = mutateInputParams;
              setup.MutateInput         = AddAnisotropy;
              RunMDOODZ("RiftingPaulineAniso.txt", &setup);
              char outputName[256];
              snprintf(outputName, sizeof(outputName), "RiftingPauline_%i.h5", i);
              rename("Output00100.gzip.h5", outputName);
              i++;
            }
          }
        }
      }
    }
  }
  free(mutateInputParams);
}
