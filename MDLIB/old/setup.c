#include "mdoodz-private.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdio.h"
#include "stdlib.h"


void BuildInitialTopography(BuildInitialTopography_ff buildInitialTopography, MdoodzInput *instance, markers *topo_chain) {
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

void SetParticles(SetParticles_ff setParticles, MdoodzInput *instance, markers *particles) {
  for (int np = 0; np < particles->Nb_part; np++) {
    Coordinates coordinates = {
            .x = particles->x[np],
            .z = particles->z[np]};
    if (setParticles.SetPhase) {
      particles->phase[np] = setParticles.SetPhase(instance, coordinates);
    } else {
      particles->phase[np] = 0;
      particles->dual[np]  = 0;
    }
    if (setParticles.SetDual) {
      particles->dual[np] = setParticles.SetPhase(instance, coordinates);
    } else {
      particles->dual[np] = particles->phase[np];
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
    if (setParticles.SetPressure) {
      particles->P[np] = setParticles.SetPressure(instance, coordinates, particles->phase[np]);
    } else {
      particles->P[np] = 0.0;
    }
    if (setParticles.SetNoise) {
      particles->noise[np] = setParticles.SetNoise(instance, coordinates, particles->phase[np]);
    } else {
      particles->noise[np] = 0.0;
    }
    ValidatePhase(particles->phase[np], instance->model.Nb_phases);
  }
}

void ValidateInternalPoint(POSITION position, char bcType, Coordinates coordinates, char *setupFunctionName) {
  if (position == INTERNAL && bcType != -1) {
    printf("Internal point MUST be set as -1 but attempted to be set as %d. Please double check your SetBCs.%s setup\n", bcType, setupFunctionName);
    printf("Particle coordinates: X: %f, Z: %f \n", coordinates.x, coordinates.z);
    printf("Make sure you haven't missed position NE, NW, SE, SW positions");
    exit(144);
  } else if (position != INTERNAL && bcType == -1) {
    printf("Point is not internal but has a -1 type. Please double check your SetBCs.%s setup\n", setupFunctionName);
    printf("Particle coordinates: X: %f, Z: %f \n", coordinates.x, coordinates.z);
    printf("Make sure you haven't missed NE, NW, SE, SW positions");
    exit(144);
  }
}

void SetBCs(SetBCs_ff setBCs, MdoodzInput *instance, grid *mesh) {
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
            position = NW;
          } else if (l == 0) {
            position = SW;
          } else {
            position = W;
          }
        } else if (k == mesh->Nx - 1) {
          if (l == mesh->Nz) {
            position = NE;
          } else if (l == 0) {
            position = SE;
          } else {
            position = E;
          }
        } else if (l == 0) {
          position = S;
        } else if (l == mesh->Nz) {
          position = N;
        } else {
          position = INTERNAL;
        }
        Coordinates coordinates = {
                .x = mesh->xg_coord[k],
                .z = mesh->zvx_coord[l]};
        SetBC bc          = setBCs.SetBCVx(instance, position, coordinates);
        mesh->BCu.type[c] = bc.type;
        mesh->BCu.val[c]  = bc.value;
        if (bc.type==11) mesh->BCu.val[c] = 2.0*bc.value;
        ValidateInternalPoint(position, bc.type, coordinates, "SetBCVxType");
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
            position = NW;
          } else if (l == 0) {
            position = SW;
          } else {
            position = W;
          }
        } else if (k == mesh->Nx) {
          if (l == mesh->Nz - 1) {
            position = NE;
          } else if (l == 0) {
            position = SE;
          } else {
            position = E;
          }
        } else if (l == 0) {
          position = S;
        } else if (l == mesh->Nz - 1) {
          position = N;
        } else {
          position = INTERNAL;
        }
        Coordinates coordinates = {
                .x = mesh->xvz_coord[k],
                .z = mesh->zg_coord[l]};
        SetBC bc          = setBCs.SetBCVz(instance, position, coordinates);
        mesh->BCv.type[c] = bc.type;
        mesh->BCv.val[c]  = bc.value;
        if (bc.type==11) mesh->BCu.val[c] = 2.0*bc.value;
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
          position = NW;
        } else if (l == 0) {
          position = SW;
        } else {
          position = W;
        }
      } else if (k == NCX - 1) {
        if (l == NCZ - 1) {
          position = NE;
        } else if (l == 0) {
          position = SE;
        } else {
          position = E;
        }
      } else if (l == 0) {
        position = S;
      } else if (l == NCZ - 1) {
        position = N;
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
          position = NW;
        } else if (l == 0) {
          position = SW;
        } else {
          position = W;
        }
      } else if (k == NCX - 1) {
        if (l == NCZ - 1) {
          position = NE;
        } else if (l == 0) {
          position = SE;
        } else {
          position = E;
        }
      } else if (l == 0) {
        position = S;
      } else if (l == NCZ - 1) {
        position = N;
      } else if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c + NCX] == 30) {
        position = FREE_SURFACE;
      } else {
        position = INTERNAL;
      }
      if (mesh->BCt.type[c] != 30) {
        if (setBCs.SetBCT) {
          SetBC bc          = setBCs.SetBCT(instance, position, mesh->T[c]);
          mesh->BCt.type[c] = bc.type;
          mesh->BCt.val[c]  = bc.value;
        }
        // TODO change size of BCt array and phase out SetBCTNew
        if (setBCs.SetBCTNew) {
          SetBC bc = setBCs.SetBCTNew(instance, position, mesh->T[c]);
          if (k == 0) {
            mesh->BCt.typW[l] = bc.type;
            mesh->BCt.valW[l] = bc.value;
          } else if (k == NCX - 1) {
            mesh->BCt.typE[l] = bc.type;
            mesh->BCt.valE[l] = bc.value;
          } else if (l == 0) {
            mesh->BCt.typS[k] = bc.type;
            mesh->BCt.valS[k] = bc.value;
          } else if (l == NCZ - 1) {
            mesh->BCt.typN[k] = bc.type;
            mesh->BCt.valN[k] = bc.value;
          }
        }
      }
    }
  }

  printf("Velocity and pressure were initialised\n");
  printf("Boundary conditions were set up\n");
  // Print2DArrayChar( mesh->BCu.type, mesh->Nx, (mesh->Nz+1), 1.0 );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ValidateSetup(MdoodzSetup *setup, MdoodzInput *instance) {
  int   errorsCount   = 0;
  char *errors[10]    = {};
  int   warningsCount = 0;
  char *warnings[10]  = {};
  if (instance->model.free_surf) {
    if (!setup->BuildInitialTopography) {
      errors[errorsCount] = "If Free Surface mode is ON and BuildInitialTopography MUST be specified. Please set free_surf = 0 or set BuildInitialTopography";
      errorsCount++;
    } else {
      if (!setup->BuildInitialTopography->SetSurfaceZCoord) {
        warnings[warningsCount] = "BuildInitialTopography.SetSurfaceZCoord is not specified. Flat surface will be generated";
        warningsCount++;
      }
      if (!setup->BuildInitialTopography->SetSurfacePhase) {
        warnings[warningsCount] = "BuildInitialTopography.SetSurfacePhase is not specified. Phase 0 will be used as surface material";
        warningsCount++;
      }
    }
  }

  if (!setup->SetParticles) {
    errors[errorsCount] = "SetParticles MUST be specified. Please set SetParticles";
    errorsCount++;
  } else {
    if (!setup->SetParticles->SetHorizontalVelocity) {
      warnings[warningsCount] = "SetParticles.SetHorizontalVelocity is not specified. Horizontal velocity will be (-x * EpsBG)";
      warningsCount++;
    }
    if (!setup->SetParticles->SetVerticalVelocity) {
      warnings[warningsCount] = "SetParticles.SetVerticalVelocity is not specified. Horizontal velocity will be (z * EpsBG)";
      warningsCount++;
    }
    if (!setup->SetParticles->SetPhase) {
      warnings[warningsCount] = "SetParticles.SetPhase is not specified. Model will be homogeneous with phase 0";
      warningsCount++;
    }
    if (!setup->SetParticles->SetDual) {
      warnings[warningsCount] = "SetParticles.SetDual is not specified. Phase value will be set";
      warningsCount++;
    }
    if (!setup->SetParticles->SetTemperature) {
      warnings[warningsCount] = "SetParticles.SetTemperature is not specified. Temperature will be set to 0°C";
      warningsCount++;
    }
    if (!setup->SetParticles->SetGrainSize) {
      warnings[warningsCount] = "SetParticles.SetGrainSize is not specified. Grain size will be set to 0.0";
      warningsCount++;
    }
    if (!setup->SetParticles->SetPorosity) {
      warnings[warningsCount] = "SetParticles.SetPorosity is not specified. Porosity will be set to 0.0";
      warningsCount++;
    }
    if (!setup->SetParticles->SetDensity) {
      warnings[warningsCount] = "SetParticles.SetDensity is not specified. Density will be set according to the particle phase";
      warningsCount++;
    }
    if (!setup->SetParticles->SetXComponent) {
      warnings[warningsCount] = "SetParticles.SetXComponent is not specified. XComponent will be set to 0.0";
      warningsCount++;
    }
    if (!setup->SetParticles->SetNoise) {
      warnings[warningsCount] = "SetParticles.SetNoise is not specified. SetNoise will be set to 0.0";
      warningsCount++;
    }
    if (!setup->SetParticles->SetPressure) {
      warnings[warningsCount] = "SetParticles.SetPressure is not specified. SetPressure will be set to 0.0";
      warningsCount++;
    }
  }

  if (!setup->SetBCs) {
    errors[errorsCount] = "SetBCs MUST be specified. Please set SetBCs";
    errorsCount++;
  } else {
    if (!setup->SetBCs->SetBCVx) {
      errors[errorsCount] = "SetBCs.SetBCVx MUST be specified";
      errorsCount++;
    }
    if (!setup->SetBCs->SetBCVz) {
      errors[errorsCount] = "SetBCs.SetBCVzType MUST be specified";
      errorsCount++;
    }
    if (!setup->SetBCs->SetBCPType) {
      warnings[warningsCount] = "SetBCs.SetBCPType is not specified. BCP type will be set to -1";
      warningsCount++;
    }
    if (instance->model.isthermal) {
      if (!setup->SetBCs->SetBCTNew) {
        errors[errorsCount] = "SetBCs.SetBCTNew MUST be specified for Thermal model. Please set isthermal = 0 or specify SetBCTTypeNew";
        errorsCount++;
      }
      if (!setup->SetBCs->SetBCT) {
        errors[errorsCount] = "SetBCs.SetBCT MUST be specified for Thermal model. Please set isthermal = 0 or specify SetBCTType (will be deprecated)";
        errorsCount++;
      }
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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// TEMPLATES

SetBC SetPureShearBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W || position == E) {
    bc.value = -coordinates.x * input->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetSimpleShearBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC        bc;
  const double Lz = (double) (input->model.zmax - input->model.zmin);
  if (position == S || position == SE || position == SW) {
    bc.value = -input->model.EpsBG * Lz;
    bc.type  = 11;
  } else if (position == N || position == NE || position == NW) {
    bc.value = input->model.EpsBG * Lz;
    bc.type  = 11;
  } else if (position == E) {
    bc.value = 0.0;
    bc.type  = -12;
  } else if (position == W) {
    bc.value = 0.0;
    bc.type  = -2;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetPureOrSimpleShearBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  if (input->model.shear_style) {
    return SetSimpleShearBCVx(input, position, coordinates);
  } else {
    return SetPureShearBCVx(input, position, coordinates);
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetPureShearBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == W || position == E || position == SW || position == SE || position == NW || position == NE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == S || position == N) {
    bc.value = coordinates.z * input->model.EpsBG;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetSimpleShearBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  SetBC bc;
  if (position == E || position == W || position == NE || position == NW || position == SE || position == SW) {
    bc.value = 0.0;
    bc.type  = -12;
  } else if (position == S || position == N) {
    bc.value = 0.0;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

SetBC SetPureOrSimpleShearBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates) {
  if (input->model.shear_style) {
    return SetSimpleShearBCVz(input, position, coordinates);
  } else {
    return SetPureShearBCVz(input, position, coordinates);
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

bool IsEllipseCoordinates(Coordinates coordinates, Ellipse ellipse, double scalingL) {
  const double theta   = ellipse.angle * M_PI / 180.0;
  const double radiusX = (ellipse.radiusX / 2) / scalingL;
  const double radiusZ = (ellipse.radiusZ / 2) / scalingL;
  const double Xn      = (coordinates.x - ellipse.centreX / scalingL) * cos(theta) - (coordinates.z + ellipse.centreZ / scalingL) * sin(theta);
  const double Zn      = (coordinates.x + ellipse.centreX / scalingL) * sin(theta) + (coordinates.z - ellipse.centreZ / scalingL) * cos(theta);
  if (pow(Xn / radiusX, 2) + pow(Zn / radiusZ, 2) - 1 < 0) {
    return true;
  } else {
    return false;
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

bool IsRectangleCoordinates(Coordinates coordinates, Rectangle rectangle, double scalingL) {
  const double wX = coordinates.x * cos(rectangle.angle) - coordinates.z * sin(rectangle.angle);
  const double wZ = coordinates.x * sin(rectangle.angle) + coordinates.z * cos(rectangle.angle);
  if (wX >= (rectangle.centreX / scalingL - (rectangle.sizeX / 2) / scalingL) && wX <= (rectangle.centreX / scalingL + (rectangle.sizeX / 2) / scalingL) && wZ >= (rectangle.centreZ / scalingL - (rectangle.sizeZ / 2) / scalingL) && wZ <= (rectangle.centreZ / scalingL + (rectangle.sizeZ / 2) / scalingL)) {
    return true;
  } else {
    return false;
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/