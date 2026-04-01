#include "math.h"
#include "mdoodz.h"
#include "stdbool.h"
#include "stdlib.h"
#include "stdio.h"

double SetSurfaceZCoord(MdoodzInput *instance, double x) {
  double h = 0.0, dh = 4e3/instance->scaling.L;
  const double Mars_radius = instance->model.planet_radius;

  // Default flat topo 
  h     = Mars_radius;

  // if (x>0) h += dh;

  if (instance->model.polar==1) {
    // Curved topography (see PolarCoordinatesStuff.py)
    h     =  sqrt((h - x)*(h + x));
  }                   
  return h;
}

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  int phase = 0; // Default: crust
  const double Lc = instance->scaling.L;
  const double x = coordinates.x, z = coordinates.z;
  double X, Z;
  const double Mars_radius = instance->model.planet_radius;
  double zMoho = Mars_radius + instance->model.user0/Lc;
  double zLAB  = Mars_radius + instance->model.user1/Lc;
  
  double w = 10e3/Lc, x1, x2, z1, z2, m;
  // x1 = -w/2, x2 = w/2;
  // z1 = Mars_radius + 1*instance->model.user0/Lc, z2 = Mars_radius + 2*instance->model.user0/Lc;
  // m = (z2 - z1) / (x2 - x1);
  // if (x<x1) zMoho = z1;
  // if (x>x2) zMoho = z2;
  // if (x>x1 && x<x2) zMoho = m * (x-x1) + z1;
  zMoho = Mars_radius - 10e3/Lc;

  // LAB
  x1 = -w/2, x2 = w/2;
  z1 = Mars_radius + 1*instance->model.user1/Lc, z2 = Mars_radius + 1.5*instance->model.user1/Lc;
  m = (z2 - z1) / (x2 - x1);
  // if (x<x1) zLAB = z1;
  // if (x>x2) zLAB = z2;
  // if (x>x1 && x<x2) zLAB = m * (x-x1) + z1
  
  zLAB = z1 + 40e3/Lc;
  if (x>x2) zLAB = z1 - 40e3/Lc;

  double ztransform = Mars_radius;
  if (x>x1 && x<x2) ztransform =  z1;

  // Weak
  double zweak = Mars_radius - 4e3/Lc;
  z1 = zweak, z2 = Mars_radius + 10e3/Lc;
  m = (z2 - z1) / (x2 - x1);
  // if (x<x1) zweak = z1;
  // if (x>x2) zweak = z2;
  // if (x>x1 && x<x2) zweak = m * (x-x1) + z1;

  if ( instance->model.polar==0 ) {     
      if (z < zMoho) {
          phase = 1;                     // Lithospheric mantle
      }
      if (z < zLAB) {
          phase = 2;                     // Asthenospheric mantle
      }
      if (z > zweak) {
          phase = 3;                     // Weak
      }
      if (z > ztransform) {
          phase = 3;                     // Weak
      }
  }
  if ( instance->model.polar==1 ) {
      Z = sqrt((zMoho - x)*(zMoho + x)); // Lithospheric mantle
      if ( z <  Z  ) phase = 1;
      Z = sqrt((zLAB  - x)*(zLAB  + x)); // Asthenospheric mantle
      if ( z <  Z  ) phase = 2;
  }

  // // Draw ellipse
  // const double angle = -35.*M_PI/180;
  // const double a_ell = 2.0*instance->model.user5/Lc, b_ell = 0.5*instance->model.user5/Lc;
  // double x_ell, z_ell;
  // Z = Mars_radius - 55e3/instance->scaling.L;
  // X = 20e3/instance->scaling.L;
  // x_ell = (x-X)*cos(angle) + (z-Z)*sin(angle);
  // z_ell =-(x-X)*sin(angle) + (z-Z)*cos(angle);
  // if (pow(x_ell/a_ell,2.0) + pow(z_ell/b_ell,2.0) < 1.0) phase = 3;
  // // Weak surface
  // const double Z_weak = Mars_radius - 5e3/instance->scaling.L;
  // if (z>Z_weak) phase = 3;
  return phase;
}

double SetTemperature(MdoodzInput *instance, Coordinates coordinates) {
  const double x = coordinates.x, z = coordinates.z;
  const double Lz = (double) (instance->model.zmax - instance->model.zmin) ;
  const double Mars_radius = instance->model.planet_radius;
  const double Ttop  = 293.0/(instance->scaling.T);
  const double Tbot  = (instance->model.user3 + zeroC)/instance->scaling.T;
  double Tgrad = (Ttop-Tbot)/ (Lz - (instance->model.zmax-Mars_radius));
  return Ttop + Tgrad*(z - Mars_radius);
}

// double SetGrainSize(MdoodzInput *instance, Coordinates coordinates, int phase) {
//   const int asthenospherePhase = 3;
//   return instance->materials.gs_ref[asthenospherePhase];
// }

// double SetHorizontalVelocity(MdoodzInput *instance, Coordinates coordinates) {
//   return -coordinates.x * instance->model.bkg_strain_rate;
// }

// double SetVerticalVelocity(MdoodzInput *instance, Coordinates coordinates) {
//   return coordinates.z * instance->model.bkg_strain_rate;
// }

char SetBCPType(MdoodzInput *instance, POSITION position) {
  if (position == NE || position == NW) {
    return 0;
  } else {
    return -1;
  }
}

SetBC SetBCT(MdoodzInput *instance, POSITION position, Coordinates coordinates,  double particleTemperature) {
  SetBC     bc;
  double surface_temperature =                   (20. + zeroC) / instance->scaling.T;
  double mantle_temperature  = (instance->model.user3 + zeroC) / instance->scaling.T;
  if (position == S) {
    bc.type  = constant_temperature;
    bc.value = mantle_temperature;
  }
  if (position == free_surface || position == N) {
    bc.type  = constant_temperature;
    bc.value = surface_temperature;
  } 
  if (position == W || position == E) {
    bc.type  = constant_heatflux;
    bc.value = 0.;
  }
  return bc;
}

void AddCrazyConductivity(MdoodzInput *input) {
  int               *asthenospherePhases = (int *) malloc(sizeof(int));
  CrazyConductivity *crazyConductivity   = (CrazyConductivity *) malloc(sizeof(CrazyConductivity));
  asthenospherePhases[0]                 = 2;
  crazyConductivity->phases              = asthenospherePhases;
  crazyConductivity->nPhases             = 1;
  crazyConductivity->multiplier          = 1000;
  input->crazyConductivity               = crazyConductivity;
}

double BoundaryVelocityProfile(double V, double z, double z_LAB, double dz_smooth) {
  return -0.5*V*erfc( (z-z_LAB)/dz_smooth ) + V;
}

double BoundaryVelocityProfilePrimitive(double V, double z, double z_LAB, double dz_smooth) {
  return V*(z - 0.5*(z-z_LAB) * erfc( (z-z_LAB)/dz_smooth ) );
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Mars_radius = instance->model.planet_radius, dz_smooth = 10e3/instance->scaling.L;
  const double Lx = instance->model.xmax - instance->model.xmin;
  double V_tot    =  Lx * instance->model.bkg_strain_rate; // |VxW| + |VxE|
  double maxAngle = asin(Lx/2/Mars_radius);     // Aperture angle
  double tet_W, alp_W, tet_E, alp_E;
  double VxW, VxE;
  double z_LAB, z_top;
  const double x = coordinates.x, z = coordinates.z;

  // Compute surface height and LAB position
  if (instance->model.polar == 1) {
    z_top = sqrt((Mars_radius - x)*(Mars_radius + x));
    z_LAB = z_top + instance->model.user1/instance->scaling.L;
  }
  else {
    z_top = Mars_radius;
    z_LAB = z_top + instance->model.user1/instance->scaling.L;
  }

  // Evaluate velocity of W and E boundaries
  if (instance->model.polar == 1) {
        const double tet_W  = -maxAngle;
        const double alp_W  = M_PI/2.0 - tet_W;
        VxW = (0.5*V_tot) * sin(alp_W) ;

        const double tet_E  = maxAngle;
        const double alp_E  = M_PI/2.0 - tet_E;
        VxE =-(0.5*V_tot) * sin(alp_E);
  }
  else {
    VxW =  0.5*V_tot;
    VxE = -0.5*V_tot;
  }

  // Apply smooth transition with depth
  VxW = BoundaryVelocityProfile(VxW, z, z_LAB, dz_smooth);
  VxE = BoundaryVelocityProfile(VxE, z, z_LAB, dz_smooth);

  // Assign BC values
  if (position == N || position == S || position == NW || position == SW || position == NE || position == SE) {
    bc.value = 0;
    bc.type  = 13;
  } else if (position == W) {
    bc.value = VxW;
    bc.type  = 0;
  } else if (position == E) {
    bc.value = VxE;
    bc.type  = 0;
  } else {
    bc.value = 0.0;
    bc.type  = -1;
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coordinates) {
  SetBC bc;
  const double Mars_radius = instance->model.planet_radius, dz_smooth = 10e3/instance->scaling.L;
  const double Lx = instance->model.xmax - instance->model.xmin;
  double V_tot    =  Lx * instance->model.bkg_strain_rate; // |VxW| + |VxE|
  double maxAngle = asin(Lx/2/Mars_radius);     // Aperture angle
  double tet_W, alp_W, tet_E, alp_E;
  double VxW, VxE, VzW, VzE, VzS=0.0;
  double z_LAB, z_top;
  const double x = coordinates.x, z = coordinates.z;  

  // Compute surface height and LAB position
  if (instance->model.polar == 1) {
    z_top = sqrt((Mars_radius - x)*(Mars_radius + x));
    z_LAB = z_top + instance->model.user1/instance->scaling.L;
  }
  else {
    z_top = Mars_radius;
    z_LAB = z_top + instance->model.user1/instance->scaling.L;
  }

  // Evaluate velocity of W and E boundaries
  if (instance->model.polar == 1) {
        const double tet_W  = -maxAngle;
        const double alp_W  = M_PI/2.0 - tet_W;
        VxW = (0.5*V_tot) * sin(alp_W) ;
        VzW = (0.5*V_tot) * cos(alp_W) ;

        const double tet_E  = maxAngle;
        const double alp_E  = M_PI/2.0 - tet_E;
        VxE =-(0.5*V_tot) * sin(alp_W) ;
        VzE =-(0.5*V_tot) * cos(alp_E);
  }
  else {
    VxW =  0.5*V_tot;
    VxE = -0.5*V_tot;
    VzW =  0.0;
    VzE = -0.0;
  }

  // // Compute compensating VzS value
  const double z_min       = instance->model.zmin; 
  const double prim_W_zmax = BoundaryVelocityProfilePrimitive(VxW, z_top, z_LAB, dz_smooth);
  const double prim_W_zmin = BoundaryVelocityProfilePrimitive(VxW, z_min, z_LAB, dz_smooth);
  const double prim_E_zmax = BoundaryVelocityProfilePrimitive(VxE, z_top, z_LAB, dz_smooth);
  const double prim_E_zmin = BoundaryVelocityProfilePrimitive(VxE, z_min, z_LAB, dz_smooth);
  const double intW = prim_W_zmax - prim_W_zmin;
  const double intE = prim_E_zmax - prim_E_zmin;
  VzS = - ( fabs(intW) + fabs(intE) ) / Lx;

  // Apply smooth transition with depth
  VzW = -BoundaryVelocityProfile(VzW, z, z_LAB, dz_smooth);
  VzE = -BoundaryVelocityProfile(VzE, z, z_LAB, dz_smooth);

  // printf("Vx %2.2e %2.2e %2.2e\n",   0*instance->scaling.V, VxW*instance->scaling.V, VxE*instance->scaling.V);
  // printf("Vz %2.2e %2.2e %2.2e\n", VzS*instance->scaling.V, VzW*instance->scaling.V, VzE*instance->scaling.V);

  if (position == W || position == SW || position == NW ) {
    bc.value = VzW;
    bc.type  = 11;
  } else if ( position == E || position == SE || position == NE) {
    bc.value = VzE;
    bc.type  = 11;
  } else if (position == S || position == N) {
    bc.value = VzS;
    bc.type  = 0;
  } else {
    bc.value = 0;
    bc.type  = -1;
  }
  return bc;
}


int main(int nargs, char *args[]) {
  // Input file name
  char *input_file;
  if ( nargs < 2 ) {
    asprintf(&input_file, "MarsPolarCartesian_retreat.txt"); // Default
  }
  else {
    asprintf(&input_file, "%s", args[1]);     // Custom
  }
  printf("Running MDoodz7.0 using %s\n", input_file);
  MdoodzSetup setup = {
          .BuildInitialTopography = &(BuildInitialTopography_ff){
                  .SetSurfaceZCoord = SetSurfaceZCoord,
          },
          .SetParticles = &(SetParticles_ff){
                  .SetPhase              = SetPhase,
                  .SetTemperature        = SetTemperature,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx    = SetBCVx,
                  .SetBCVz    = SetBCVz,
                  .SetBCPType = SetBCPType,
                  .SetBCT     = SetBCT,
          },
          .MutateInput = AddCrazyConductivity,

  };
  RunMDOODZ(input_file, &setup);
  free(input_file);
}
