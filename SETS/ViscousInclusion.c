#include "complex.h"
#include "mdoodz.h"
#include "math.h"

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

int SetPhase(MdoodzInput *instance, Coordinates coordinates) {
  const double radius = instance->model.user1 / instance->scaling.L;
  if (coordinates.x * coordinates.x + coordinates.z * coordinates.z < radius * radius) {
    return 1;
  } else {
    return 0;
  }
}

double SetDensity(MdoodzInput *instance, Coordinates coordinates, int phase) {
  return instance->materials.rho[phase];
}

SetBC SetBCVx(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm     = 1.0;
  const double mc     = 1e3;
  double       Vx, Vz, P, eta, sxx, szz, x, z;
  if (instance->model.shear_style == 0) {
    if (position == W || position == E) {
      x = coord.x;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 0;
      bc.value = Vx;
    } else if (position == S || position == SE || position == SW) {
      x = coord.x;
      z = coord.z + instance->model.dx / 2.0;// make sure it lives on the boundary
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 11;
      bc.value = 2.0*Vx;
    } else if (position == N || position == NE || position == NW) {
      x = coord.x;
      z = coord.z - instance->model.dx / 2.0;// make sure it lives on the boundary
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 11;
      bc.value = 2.0*Vx;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  } else {
    const double Lz = (double) (instance->model.zmax - instance->model.zmin);
    if (position == S || position == SE || position == SW) {
      x = coord.x;
      z = coord.z + instance->model.dx / 2.0;// make sure it lives on the boundary
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 11;
      bc.value = 2.0*Vx;
    } else if (position == N || position == NE || position == NW) {
      x = coord.x;
      z = coord.z - instance->model.dx / 2.0;// make sure it lives on the boundary
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 11;
      bc.value = 2.0*Vx;
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
  }
  return bc;
}

SetBC SetBCVz(MdoodzInput *instance, POSITION position, Coordinates coord) {
  SetBC           bc;
  const double radius = instance->model.user1 / instance->scaling.L;
  const double mm     = 1.0;
  const double mc     = 1e3;
  double       Vx, Vz, P, eta, sxx, szz, x, z;
  if (instance->model.shear_style == 0) {
    if (position == N || position == NE || position == NW || position == S || position == SE || position == SW) {
      x = coord.x;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 0;
      bc.value = Vz;
    }
    else if (position == W) {
      x = coord.x + instance->model.dx / 2.0;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 11;
      bc.value = 2.0*Vz;
    }
    else if (position == E) {
      x = coord.x + instance->model.dx / 2.0;
      z = coord.z;
      eval_anal_Dani(&Vx, &Vz, &P, &eta, &sxx, &szz, x, z, 1, radius, mm, mc);
      bc.type  = 11;
      bc.value = 2.0*Vz;
    } else {
      bc.type  = -1;
      bc.value = 0.0;
    }
  } 
  else {
    if (position == E || position == W || position == NE || position == NW || position == SE || position == SW) {
      bc.value = 0.0;
      bc.type = -12;
    } else if (position == S || position == N) {
      bc.value = 0.0;
      bc.type = 0;
    } else {
      bc.value = 0.0;
      bc.type = -1;
    }
  }
  return bc;
}

int main() {
  MdoodzSetup instance = {
          .SetParticles = &(SetParticles_ff){
                  .SetPhase   = SetPhase,
                  .SetDensity = SetDensity,
          },
          .SetBCs = &(SetBCs_ff){
                  .SetBCVx = SetBCVx,
                  .SetBCVz = SetBCVz,
          },
  };
  RunMDOODZ("ViscousInclusion.txt", &instance);
}
