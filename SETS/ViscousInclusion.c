#include "complex.h"
#include "mdoodz.h"

//---------------------------------------------//
void eval_anal_Dani(double *vx, double *vz, double *p, double *eta, double *sxx, double *syy, double x, double z, int ps, double rc, double mm, double mc) {

  double          gr, er, A;
  _Complex double V_tot, phi_z, d_phi_z, conj_d_phi_z, psi_z, conj_psi_z, Z, d_d_phi_z_z, d_psi_z;

  /*
	 % ---------------------------------------------------------------------------
	 % ANALYTICAL SOLUTION - PRESSURE AND VELOCITY AROUND A CIRCULAR INCLUSION:
	 %
	 % BASED ON DANI SCHMID'S 2002 CYL_P_MATRIX.M
	 % FAR FIELD FLOW - VISCOSITIES - GEOMETRY
	 % ---------------------------------------------------------------------------
	 */

  *eta = *vx = *vz = *p = 0;

  // INPUT:
  if (ps == 1) {
    gr = 0; // Pure shear: gr=0 er=-1
    er = -1;// Strain rate
  } 
  else {
    gr = -1;// Simple shear: gr=1, er=0
    er = 0; // Strain rate
  }
  A = mm * (mc - mm) / (mc + mm);

  /*
	 % --------------------------------------------------------------
	 % PRESSURE CALCULATION OUTSIDE OF AN INCLUSION IN THE Z-PLANE
	 % --------------------------------------------------------------
	 */

  Z = x + I * z;

  if (sqrt(pow(x, 2) + pow(z, 2)) <= rc) {
    // INSIDE CLAST
    *p    = 0.0;
    V_tot = (mm / (mc + mm)) * (I * gr + 2 * er) * conj(Z) - (I / 2) * gr * Z;
    *vx   = creal(V_tot);
    *vz   = cimag(V_tot);
    *eta  = mc;
    *syy  = 3.5;
    *sxx  = -3.5;
  } 
  else {
    // OUTSIDE CLAST, RESP. MATRIX
    *p           = (-2 * mm * (mc - mm) / (mc + mm) * creal(cpow(rc, 2) / cpow(Z, 2) * (I * gr + 2 * er)));
    phi_z        = -(I / 2) * mm * gr * Z - (I * gr + 2 * er) * A * cpow(rc, 2) * cpow(Z, (-1));
    d_phi_z      = -(I / 2) * mm * gr + (I * gr + 2 * er) * A * cpow(rc, 2) / cpow(Z, 2);
    conj_d_phi_z = conj(d_phi_z);
    psi_z        = (I * gr - 2 * er) * mm * Z - (I * gr + 2 * er) * A * cpow(rc, 4) * cpow(Z, (-3));
    conj_psi_z   = conj(psi_z);
    V_tot        = (phi_z - Z * conj_d_phi_z - conj_psi_z) / (2 * mm);
    *vx          = creal(V_tot);
    *vz          = cimag(V_tot);
    *eta         = mm;
    // Evaluate stresses
    d_d_phi_z_z  = -(I * gr + 2 * er) * A * pow(rc, 2) / cpow(Z, 3);
    d_psi_z      = (I * gr - 2 * er) * mm + (I * gr + 2 * er) * A * pow(rc, 4) * cpow(Z, -4);
    *syy         = 2 * creal(d_phi_z) + creal((x - I * z) * d_d_phi_z_z + d_psi_z);
    *sxx         = 4 * creal(d_phi_z) - *syy;
  }
}
//---------------------------------------------//

int SetPhase(MdoodzInstance *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInstance *instance, Coordinates coordinates, int phase) {
  return instance->materials.rho[phase];
}

char SetBCVxType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == SOUTH || position == NORTH || position == NORTHWEST || position == SOUTHWEST || position == NORTHEAST || position == SOUTHEAST) {
      return 11;
    } else if (position == WEST || position == EAST) {
      return 0;
    }
    else {
      return -1;
    }
  } else {
    if (position == WEST || position == NORTHWEST || position == SOUTHWEST) {
      return -2;
    } else if (position == EAST || position == NORTHEAST || position == SOUTHEAST) {
      return -12;
    } else if (position == SOUTH || position == NORTH) {
      return 11;
    } else {
      return -1;
    }
  }
}

double SetBCVxValue(MdoodzInstance *instance, POSITION position, Coordinates coord) {
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm     = 1.0;
  const double mc     = 1e3;
  double       Vx, Vz, P, eta, sxx, szz, x, z;
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST ) {
      x = coord.x;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      return Vx;
    } else if (position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {;
      x = coord.x;
      z = coord.z + instance->model.dx / 2.0; // make sure it lives on the boundary
      double       Vx, Vz, P, eta, sxx, szz;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      return Vx;
    } else if (position == NORTH || position == NORTHEAST || position == NORTHWEST ) {
      x = coord.x;
      z = coord.z - instance->model.dx / 2.0; // make sure it lives on the boundary
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      return Vx;
    } else {
      return 0.0;
    }
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == SOUTH) {
      return -instance->model.EpsBG * Lz;
    } else if (position == NORTH) {
      return instance->model.EpsBG * Lz;
    } else {
      return 0;
    }
  }
}

char SetBCVzType(MdoodzInstance *instance, POSITION position) {
  if (instance->model.shear_style == 0) {
    if (position == WEST || position == EAST || position == NORTHEAST || position == NORTHWEST || position == SOUTHEAST || position == SOUTHWEST) {
      return 11;
    } else if (position == SOUTH || position == NORTH) {
      return 0;
    } else {
      return -1;
    }
  }
  else {
    if (position == WEST || position == EAST) {
      return 0;
    } else if (position == SOUTH || position == NORTH) {
      return -12;
    } else {
      return -1;
    }
  }
}

double SetBCVzValue(MdoodzInstance *instance, POSITION position, Coordinates coord) {
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm     = 1.0;
  const double mc     = 1e3;
  double       Vx, Vz, P, eta, sxx, szz, x, z;
  if (instance->model.shear_style == 0) {
    if (position == NORTH || position == NORTHEAST || position == NORTHWEST || position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {
      x = coord.x;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      return Vz;
    } 
    if (position == WEST) {
      x = coord.x + instance->model.dx/2.0;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      return Vz;
    }
    if (position == EAST) {
      x = coord.x + instance->model.dx/2.0;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      return Vz;
    }
    else {
      return 0;
    }
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == WEST || position == EAST || position == NORTHWEST || position == SOUTH || position == SOUTHEAST || position == SOUTHWEST) {
      return 0.0 * instance->model.EpsBG * Lz;
    } else {
      return 0;
    }
  }
}

int main(int nargs, char *args[]) {
  MdoodzInstance instance = {
          .inputFileName = GetSetupFileName(nargs, args),
          .SetParticles  = &(SetParticles_ff){
                   .SetPhase   = SetPhase,
                   .SetDensity = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVxType  = SetBCVxType,
                  .SetBCVxValue = SetBCVxValue,
                  .SetBCVzType  = SetBCVzType,
                  .SetBCVzValue = SetBCVzValue,
          },
  };
  RunMDOODZ(&instance);
}
