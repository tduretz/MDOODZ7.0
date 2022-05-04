#include "mdoodz-private.h"
#include "mdoodz.h"
#include "stdio.h"
#include "stdlib.h"


void BuildInitialTopography(MdoodzInstance *instance, markers *topo_chain) {
  BuildInitialTopography_ff buildInitialTopography = *instance->BuildInitialTopography;
  for (int k = 0; k < topo_chain->Nb_part; k++) {
    const double x_coord = topo_chain->x[k];
    if (buildInitialTopography.SetSurfaceZCoord) {
      topo_chain->z[k] = buildInitialTopography.SetSurfaceZCoord(instance, x_coord);
    } else {
      topo_chain->z[k] = 1.0e3 / instance->scaling.L;
    }
    if (buildInitialTopography.SetSurfacePhase) {
      topo_chain->phase[k] = buildInitialTopography.SetSurfacePhase(instance, x_coord);
    } else {
      topo_chain->phase[k] = 0;
    }
  }
  printf("Topographic chain initialised with %d markers\n",
         topo_chain->Nb_part);
}

void ValidatePhase(int phaseId, int phasesCount) {
  if (phaseId > phasesCount) {
    printf("Lazy bastard! Fix your particle phase ID! \n");
    exit(144);
  }
}

void SetParticles(MdoodzInstance *instance, markers *particles) {
  SetParticles_ff setParticles = *instance->SetParticles;
  for (int np = 0; np < particles->Nb_part; np++) {
    Coordinates coordinates = {
            .x = particles->x[np],
            .z = particles->z[np]};
    if (setParticles.SetPhase) {
      particles->phase[np] = setParticles.SetPhase(instance, coordinates);
      particles->dual[np]  = setParticles.SetPhase(instance, coordinates);
    } else {
      particles->phase[np] = 0;
      particles->dual[np]  = 0;
    }
    if (setParticles.SetHorizontalVelocity) {
      particles->Vx[np] = setParticles.SetHorizontalVelocity(instance, coordinates);
    } else {
      particles->Vx[np] = -coordinates.x * instance->model.EpsBG;
    }
    if (setParticles.SetVerticalVelocity) {
      particles->Vz[np] = setParticles.SetVerticalVelocity(instance, coordinates);
    } else {
      particles->Vz[np] = coordinates.z * instance->model.EpsBG;
    }
    if (setParticles.SetGrainSize) {
      particles->d[np] = setParticles.SetGrainSize(instance, coordinates, particles->phase[np]);
    } else {
      particles->d[np] = 0.0;
    }
    if (setParticles.SetPorosity) {
      particles->phi[np] = setParticles.SetPorosity(instance, coordinates, particles->phase[np]);
    } else {
      particles->phi[np] = 0.0;
    }
    if (setParticles.SetDensity) {
      particles->rho[np] = setParticles.SetDensity(instance, coordinates, particles->phase[np]);
    } else {
      particles->rho[np] = instance->materials.rho[particles->phase[np]];
    }
    if (setParticles.SetXComponent) {
      particles->X[np] = setParticles.SetXComponent(instance, coordinates, particles->phase[np]);
    } else {
      particles->X[np] = 0.0;
    }
    if (setParticles.SetTemperature) {
      particles->T[np] = setParticles.SetTemperature(instance, coordinates);
    } else {
      particles->T[np] = zeroC / instance->scaling.T;
    }
    ValidatePhase(particles->phase[np], instance->model.Nb_phases);
  }
}

void ValidateInternalPoint(POSITION position, char bcType, Coordinates coordinates, char *setupFunctionName) {
  if (position == INTERNAL && bcType != -1) {
    printf("Internal point MUST be set as -1 but attempted to be set as %d. Please double check your SetBCs.%s setup\n", bcType, setupFunctionName);
    printf("Particle coordinates: X: %f, Z: %f \n", coordinates.x, coordinates.z);
    printf("Make sure you haven't missed position NORTHEAST, NORTHWEST, SOUTHEAST, SOUTHWEST positions");
    exit(144);
  } else if (position != INTERNAL && bcType == -1) {
    printf("Point is not internal but has a -1 type. Please double check your SetBCs.%s setup\n", setupFunctionName);
    printf("Particle coordinates: X: %f, Z: %f \n", coordinates.x, coordinates.z);
    printf("Make sure you haven't missed NORTHEAST, NORTHWEST, SOUTHEAST, SOUTHWEST positions");
    exit(144);
  }
}

void SetBCs(MdoodzInstance *instance, grid *mesh) {
  SetBCs_ff setBCs = *instance->SetBCs;
  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for Vx on all grid levels */
  /* Type  0: Dirichlet point that matches the physical boundary (Vx:
   * left/right, Vz: bottom/top)            */
  /* Type 11: Dirichlet point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)       */
  /* Type  2: Neumann point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)         */
  /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top,
   * Vz: left/right)              */
  /* Type -2: periodic in the x direction (matches the physical boundary) */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* --------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < mesh->Nz + 1; l++) {
    for (int k = 0; k < mesh->Nx; k++) {
      const int c = k + l * (mesh->Nx);
      if (mesh->BCu.type[c] != 30) {
        POSITION position;
        if (k == 0) {
          if (l == mesh->Nz) {
            position = NORTHWEST;
          } else if (l == 0) {
            position = SOUTHWEST;
          } else {
            position = WEST;
          }
        } else if (k == mesh->Nx - 1) {
          if (l == mesh->Nz) {
            position = NORTHEAST;
          } else if (l == 0) {
            position = SOUTHEAST;
          } else {
            position = EAST;
          }
        } else if (l == 0) {
          position = SOUTH;
        } else if (l == mesh->Nz) {
          position = NORTH;
        } else {
          position = INTERNAL;
        }
        mesh->BCu.type[c]       = setBCs.SetBCVxType(instance, position);
        Coordinates coordinates = {
                .x = mesh->xg_coord[k],
                .z = mesh->zvx_coord[l]};
        mesh->BCu.val[c] = setBCs.SetBCVxValue(instance, position, coordinates);
        ValidateInternalPoint(position, mesh->BCu.type[c], coordinates, "SetBCVxType");
      }
    }
  }

  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for Vz on all grid levels */
  /* Type  0: Dirichlet point that matches the physical boundary (Vx:
   * left/right, Vz: bottom/top)            */
  /* Type 11: Dirichlet point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)       */
  /* Type  2: Neumann point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)         */
  /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top,
   * Vz: left/right)              */
  /* Type -2: periodic in the x direction (does not match the physical boundary)
   */
  /* Type-10: useless point (set to zero) */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* --------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < mesh->Nz; l++) {
    for (int k = 0; k < mesh->Nx + 1; k++) {
      const int c = k + l * (mesh->Nx + 1);
      if (mesh->BCv.type[c] != 30) {
        POSITION position;
        if (k == 0) {
          if (l == mesh->Nz - 1) {
            position = NORTHWEST;
          } else if (l == 0) {
            position = SOUTHWEST;
          } else {
            position = WEST;
          }
        } else if (k == mesh->Nx) {
          if (l == mesh->Nz - 1) {
            position = NORTHEAST;
          } else if (l == 0) {
            position = SOUTHEAST;
          } else {
            position = EAST;
          }
        } else if (l == 0) {
          position = SOUTH;
        } else if (l == mesh->Nz - 1) {
          position = NORTH;
        } else {
          position = INTERNAL;
        }
        mesh->BCv.type[c]       = setBCs.SetBCVzType(instance, position);
        Coordinates coordinates = {
                .x = mesh->xvz_coord[k],
                .z = mesh->zg_coord[l]};
        mesh->BCv.val[c] = setBCs.SetBCVzValue(instance, position, coordinates);
        ValidateInternalPoint(position, mesh->BCv.type[c], coordinates, "SetBCVzType");
      }
    }
  }

  const int NCX = mesh->Nx - 1;
  const int NCZ = mesh->Nz - 1;

  /* --------------------------------------------------------------------------------------------------------*/
  /* Set the BCs for P on all grid levels */
  /* Type  0: Dirichlet within the grid */
  /* Type -1: not a BC point (tag for inner points) */
  /* Type 30: not calculated (part of the "air") */
  /* Type 31: surface pressure (Dirichlet) */
  /* --------------------------------------------------------------------------------------------------------*/

  for (int l = 0; l < NCZ; l++) {
    for (int k = 0; k < NCX; k++) {
      const int c = k + l * (NCX);
      POSITION  position;
      if (k == 0) {
        if (l == NCZ - 1) {
          position = NORTHWEST;
        } else if (l == 0) {
          position = SOUTHWEST;
        } else {
          position = WEST;
        }
      } else if (k == NCX - 1) {
        if (l == NCZ - 1) {
          position = NORTHEAST;
        } else if (l == 0) {
          position = SOUTHEAST;
        } else {
          position = EAST;
        }
      } else if (l == 0) {
        position = SOUTH;
      } else if (l == NCZ - 1) {
        position = NORTH;
      } else if ((mesh->BCp.type[c] == -1 || mesh->BCp.type[c] == 1 || mesh->BCp.type[c] == 0) && mesh->BCp.type[c + NCX] == 30) {
        position = FREE_SURFACE;
      } else {
        position = INTERNAL;
      }
      if (mesh->BCt.type[c] != 30) {
        if (setBCs.SetBCPType) {
          mesh->BCp.type[c] = setBCs.SetBCPType(instance, position);
        } else {
          mesh->BCp.type[c] = -1;
        }
        mesh->BCp.val[c] = 0.0;
      }
    }
  }

  if (instance->model.isthermal) {
    /* -------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for T on all grid levels */
    /* Type  1: Dirichlet point that do not match the physical boundary (Vx:
   * bottom/top, Vz: left/right)      */
    /* Type  0: Neumann point that matches the physical boundary (Vx: bottom/top,
   * Vz: left/right)             */
    /* Type -2: periodic in the x direction (matches the physical boundary) */
    /* Type -1: not a BC point (tag for inner points) */
    /* Type 30: not calculated (part of the "air") */
    /* -------------------------------------------------------------------------------------------------------*/

    for (int l = 0; l < NCZ; l++) {
      for (int k = 0; k < NCX; k++) {
        const int c = k + l * (NCX);
        POSITION  position;
        if (k == 0) {
          if (l == NCZ - 1) {
            position = NORTHWEST;
          } else if (l == 0) {
            position = SOUTHWEST;
          } else {
            position = WEST;
          }
        } else if (k == NCX - 1) {
          if (l == NCZ - 1) {
            position = NORTHEAST;
          } else if (l == 0) {
            position = SOUTHEAST;
          } else {
            position = EAST;
          }
        } else if (l == 0) {
          position = SOUTH;
        } else if (l == NCZ - 1) {
          position = NORTH;
        } else if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c + NCX] == 30) {
          position = FREE_SURFACE;
        } else {
          position = INTERNAL;
        }
        if (mesh->BCt.type[c] != 30) {
          mesh->BCt.type[c] = setBCs.SetBCTType(instance, position);
          mesh->BCt.val[c]  = setBCs.SetBCTValue(instance, position, mesh->T[c]);
          // TODO change size of BCt array and phase out SetBCTTypeNew
          if (setBCs.SetBCTTypeNew) {
            if (k == 0) {
              mesh->BCt.typW[l] = setBCs.SetBCTTypeNew(instance, position);
              mesh->BCt.valW[l] = setBCs.SetBCTValueNew(instance, position, mesh->T[c]);
            } else if (k == NCX - 1) {
              mesh->BCt.typE[l] = setBCs.SetBCTTypeNew(instance, position);
              mesh->BCt.valE[l] = setBCs.SetBCTValueNew(instance, position, mesh->T[c]);
            } else if (l == 0) {
              mesh->BCt.typS[k] = setBCs.SetBCTTypeNew(instance, position);
              mesh->BCt.valS[k] = setBCs.SetBCTValueNew(instance, position, mesh->T[c]);
            } else if (l == NCZ - 1) {
              mesh->BCt.typN[k] = setBCs.SetBCTTypeNew(instance, position);
              mesh->BCt.valN[k] = setBCs.SetBCTValueNew(instance, position, mesh->T[c]);
            }
          }
        }
      }
    }
  }

  printf("Velocity and pressure were initialised\n");
  printf("Boundary conditions were set up\n");
}

void ValidateSetup(MdoodzInstance *instance) {
  int   errorsCount   = 0;
  char *errors[10]    = {};
  int   warningsCount = 0;
  char *warnings[10]  = {};
  if (instance->model.free_surf) {
    if (!instance->BuildInitialTopography) {
      errors[errorsCount] = "If Free Surface mode is ON and BuildInitialTopography MUST be specified. Please set free_surf = 0 or set BuildInitialTopography";
      errorsCount++;
    } else {
      if (!instance->BuildInitialTopography->SetSurfaceZCoord) {
        warnings[warningsCount] = "BuildInitialTopography.SetSurfaceZCoord is not specified. Flat surface will be generated";
        warningsCount++;
      }
      if (!instance->BuildInitialTopography->SetSurfacePhase) {
        warnings[warningsCount] = "BuildInitialTopography.SetSurfacePhase is not specified. Phase 0 will be used as surface material";
        warningsCount++;
      }
    }
  }

  if (!instance->SetParticles) {
    errors[errorsCount] = "SetParticles MUST be specified. Please set SetParticles";
    errorsCount++;
  } else {
    if (!instance->SetParticles->SetHorizontalVelocity) {
      warnings[warningsCount] = "SetParticles.SetHorizontalVelocity is not specified. Horizontal velocity will be (-x * EpsBG)";
      warningsCount++;
    }
    if (!instance->SetParticles->SetVerticalVelocity) {
      warnings[warningsCount] = "SetParticles.SetVerticalVelocity is not specified. Horizontal velocity will be (z * EpsBG)";
      warningsCount++;
    }
    if (!instance->SetParticles->SetPhase) {
      warnings[warningsCount] = "SetParticles.SetPhase is not specified. Model will be homogeneous with phase 0";
      warningsCount++;
    }
    if (!instance->SetParticles->SetTemperature) {
      warnings[warningsCount] = "SetParticles.SetTemperature is not specified. Temperature will be set to 0°C";
      warningsCount++;
    }
    if (!instance->SetParticles->SetGrainSize) {
      warnings[warningsCount] = "SetParticles.SetGrainSize is not specified. Grain size will be set to 0.0";
      warningsCount++;
    }
    if (!instance->SetParticles->SetPorosity) {
      warnings[warningsCount] = "SetParticles.SetPorosity is not specified. Porosity will be set to 0.0";
      warningsCount++;
    }
    if (!instance->SetParticles->SetDensity) {
      warnings[warningsCount] = "SetParticles.SetDensity is not specified. Density will be set according to the particle phase";
      warningsCount++;
    }
    if (!instance->SetParticles->SetXComponent) {
      warnings[warningsCount] = "SetParticles.SetXComponent is not specified. XComponent will be set to 0.0";
      warningsCount++;
    }
  }

  if (!instance->SetBCs) {
    errors[errorsCount] = "SetBCs MUST be specified. Please set SetBCs";
    errorsCount++;
  } else {
    if (!instance->SetBCs->SetBCVxType) {
      errors[errorsCount] = "SetBCs.SetBCVxType MUST be specified";
      errorsCount++;
    }
    if (!instance->SetBCs->SetBCVxValue) {
      errors[errorsCount] = "SetBCs.SetBCVxValue MUST be specified";
      errorsCount++;
    }
    if (!instance->SetBCs->SetBCVzType) {
      errors[errorsCount] = "SetBCs.SetBCVzType MUST be specified";
      errorsCount++;
    }
    if (!instance->SetBCs->SetBCVzValue) {
      errors[errorsCount] = "SetBCs.SetBCVzValue MUST be specified";
      errorsCount++;
    }

    if (!instance->SetBCs->SetBCPType) {
      warnings[warningsCount] = "SetBCs.SetBCPType is not specified. BCP type will be set to -1";
      warningsCount++;
    }

    if (instance->model.isthermal) {
      if (!instance->SetBCs->SetBCTTypeNew) {
        errors[errorsCount] = "SetBCs.SetBCTTypeNew MUST be specified for Thermal model. Please set isthermal = 0 or specify SetBCTTypeNew";
        errorsCount++;
      }
      if (!instance->SetBCs->SetBCTValueNew) {
        errors[errorsCount] = "SetBCs.SetBCTTypeNew MUST be specified for Thermal model. Please set isthermal = 0 or specify SetBCTValueNew";
        errorsCount++;
      }
      if (!instance->SetBCs->SetBCTType) {
        errors[errorsCount] = "SetBCs.SetBCTType MUST be specified for Thermal model. Please set isthermal = 0 or specify SetBCTType (will be deprecated)";
        errorsCount++;
      }
      if (!instance->SetBCs->SetBCTValue) {
        errors[errorsCount] = "SetBCs.SetBCTValue MUST be specified for Thermal model. Please set isthermal = 0 or specify SetBCTValue (will be deprecated)";
        errorsCount++;
      }
    }
  }

  if (instance->crazyConductivity) {
    if (!instance->crazyConductivity->multiplier) {
      errors[errorsCount] = "Parameter multiplier MUST be set for crazy conductivity";
      errorsCount++;
    }
    if (!instance->crazyConductivity->phases) {
      errors[errorsCount] = "Parameter phases MUST be set for crazy conductivity";
      errorsCount++;
    }
    if (!instance->crazyConductivity->nPhases) {
      errors[errorsCount] = "Parameter nPhases MUST be set for crazy conductivity";
      errorsCount++;
    } else if (instance->crazyConductivity->nPhases > instance->model.Nb_phases) {
      errors[errorsCount] = "Asthenosphere phases for crazy conductivity is more than phases in a model. Please double check it";
      errorsCount++;
    }
  }

  if (warningsCount) {
    printf("\n\n******************  YOU HAVE %d SETUP WARNINGS  ******************\n", warningsCount);
    for (int i = 0; i < warningsCount; i++) {
      printf("%d) %s\n", i + 1, warnings[i]);
    }
    printf("******************************************************************\n");
  }

  if (errorsCount) {
    printf("\n\n*******************  YOU HAVE %d SETUP ERRORS  ********************\n", errorsCount);
    for (int i = 0; i < errorsCount; i++) {
      printf("%d) %s\n", i + 1, errors[i]);
    }
    printf("******************************************************************\n");
    exit(144);
  }
}
