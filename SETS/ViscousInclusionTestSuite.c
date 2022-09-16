#include "complex.h"
#include "mdoodz.h"
#include "math.h"

#if 0
// Set physical properties on the grid and boundary conditions
void SetBCs_MD6( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {

    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta ;
    double Lx, Lz, T1, T2, rate=model->EpsBG,  z_comp=-140e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    double Vfix = (50.0/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [50.0 == 5 cm/yr]
    
    // Define dimensions;
    Lx = (double) (model->xmax - model->xmin) ;
    Lz = (double) (model->zmax - model->zmin) ;
    
    // ---- T-Dependent marker types
    // -------------------- SPECIFIC TO YOANN's SETUP -------------------- //
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    X  = malloc (NX*sizeof(double));
    Z  = malloc (NZ*sizeof(double));
    XC = malloc (NCX*sizeof(double));
    ZC = malloc (NCZ*sizeof(double));
    
    for (k=0; k<NX; k++) {
        X[k] = mesh->xg_coord[k];
    }
    for (k=0; k<NCX; k++) {
        XC[k] = mesh->xc_coord[k];
    }
    for (l=0; l<NZ; l++) {
        Z[l] = mesh->zg_coord[l];
    }
    for (l=0; l<NCZ; l++) {
        ZC[l] = mesh->zc_coord[l];
    }
    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for Vx on all grid levels                                                                   */
    /* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
    /* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
    /* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
    /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
    /* Type -2: periodic in the x direction (matches the physical boundary)                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* --------------------------------------------------------------------------------------------------------*/
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<mesh->Nz+1; l++) {
            for (k=0; k<mesh->Nx; k++) {
                
                c = k + l*(mesh->Nx);
                
                if ( mesh->BCu.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCu.type[c] = -1;
                    mesh->BCu.val[c]  =  0;
                    
                    if (model->shear_style== 0 ) {
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
                    }
                    
                    // Free slip SOUTH
                    if (l==0  ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                    
                    // Free slip NORTH
                    if ( l==mesh->Nz ) {
                        mesh->BCu.type[c] = 13;
                        mesh->BCu.val[c]  =  0;
                    }
                        
                    }
                    if (model->shear_style== 1 ) {
                        
                        // Matching BC nodes WEST
                        if (k==0 ) {
                            mesh->BCu.type[c] = -2;
                            mesh->BCu.val[c]  = 0.0*model->EpsBG*Lx;
                        }
                        
                        // Matching BC nodes EAST
                        if (k==mesh->Nx-1 ) {
                            mesh->BCu.type[c] =  -12;
                            mesh->BCu.val[c]  = -0.0*model->EpsBG*Lx;
                        }
                        
                        // Free slip S
                        if (l==0 ) { //&& (k>0 && k<NX-1) ) {
                            mesh->BCu.type[c] =  11;
                            mesh->BCu.val[c]  = -1*model->EpsBG*Lz;
                        }
                        
                        // Free slip N
                        if ( l==mesh->Nz) {// && (k>0 && k<NX-1)) {
                            mesh->BCu.type[c] =  11;
                            mesh->BCu.val[c]  =  1*model->EpsBG*Lz;
                        }
                        
                    }
                }
                
            }
        }
        
    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for Vz on all grid levels                                                                   */
    /* Type  0: Dirichlet point that matches the physical boundary (Vx: left/right, Vz: bottom/top)            */
    /* Type 11: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)       */
    /* Type  2: Neumann point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)         */
    /* Type 13: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)              */
    /* Type -2: periodic in the x direction (does not match the physical boundary)                             */
    /* Type-10: useless point (set to zero)                                                                    */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* --------------------------------------------------------------------------------------------------------*/
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<mesh->Nz; l++) {
            for (k=0; k<mesh->Nx+1; k++) {
                
                c  = k + l*(mesh->Nx+1);
                
                if ( mesh->BCv.type[c] != 30 ) {
                    
                    // Internal points:  -1
                    mesh->BCv.type[c] = -1;
                    mesh->BCv.val[c]  =  0;
                    
                    if (model->shear_style== 0 ) {
                    
                        // Matching BC nodes SOUTH
                        if (l==0 ) {
                            mesh->BCv.type[c] = 0;
                            mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                        }
                        
                        // Matching BC nodes NORTH
                        if (l==mesh->Nz-1 ) {
                            mesh->BCv.type[c] = 0;
                            mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                        }
                        
                        // Non-matching boundary WEST
                        if ( (k==0) ) {
                            mesh->BCv.type[c] =   13;
                            mesh->BCv.val[c]  =   0;
                        }
                        
                        // Non-matching boundary EAST
                        if ( (k==mesh->Nx) ) {
                            mesh->BCv.type[c] =   13;
                            mesh->BCv.val[c]  =   0;
                        }
                    }
                    
                    
                    if (model->shear_style== 1 ) {
                        // Matching BC nodes SOUTH
                        if (l==0 ) {
                            mesh->BCv.type[c] = 0;
                            mesh->BCv.val[c]  = -0.0*model->EpsBG*Lz;
                        }
                        
                        // Matching BC nodes NORTH
                        if (l==mesh->Nz-1 ) {
                            mesh->BCv.type[c] = 0;
                            mesh->BCv.val[c]  = 0.0*model->EpsBG*Lz;
                        }
                        
                        // Non-matching boundary points
                        if ( (k==0)   ) {    //&& (l>0 && l<NZ-1)
                            mesh->BCv.type[c] =  -12;
                            mesh->BCv.val[c]  =   0;
                        }
                        
                        // Non-matching boundary points
                        if ( (k==mesh->Nx)  ) { // && (l>0 && l<NZ-1)
                            mesh->BCv.type[c] =  -12;
                            mesh->BCv.val[c]  =   0;
                        }
                    }
                    
                }
                
            }
        }
        
    
    /* --------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for P on all grid levels                                                                    */
    /* Type  0: Dirichlet within the grid                                                                      */
    /* Type -1: not a BC point (tag for inner points)                                                          */
    /* Type 30: not calculated (part of the "air")                                                             */
    /* Type 31: surface pressure (Dirichlet)                                                                   */
    /* --------------------------------------------------------------------------------------------------------*/
    
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<NCZ; l++) {
            for (k=0; k<NCX; k++) {
                
                c  = k + l*(NCX);
                
                if (mesh->BCt.type[c] != 30) {
                            
                    // Internal points:  -1
                    mesh->BCp.type[c] = -1;
                    mesh->BCp.val[c]  =  0;
                }
            }
        }
    
    /* -------------------------------------------------------------------------------------------------------*/
    /* Set the BCs for T on all grid levels                                                                   */
    /* Type  1: Dirichlet point that do not match the physical boundary (Vx: bottom/top, Vz: left/right)      */
    /* Type  0: Neumann point that matches the physical boundary (Vx: bottom/top, Vz: left/right)             */
    /* Type -2: periodic in the x direction (matches the physical boundary)                                   */
    /* Type -1: not a BC point (tag for inner points)                                                         */
    /* Type 30: not calculated (part of the "air")                                                            */
    /* -------------------------------------------------------------------------------------------------------*/
    
    double Ttop = 273.15/scaling.T;
    double Tbot, Tleft, Tright;
    
    
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        NXVZ = NX+1;
        NZVX = NZ+1;
        
        for (l=0; l<mesh->Nz-1; l++) {
            for (k=0; k<mesh->Nx-1; k++) {
                
                c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {
                    
                    // LEFT
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // RIGHT
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // BOT
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    
                    // TOP
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.val[c]  = mesh->T[c];
                    }
                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = Ttop;
                        }
                    }
                    
                    
                }
                
            }
        }
        
    
    free(X);
    free(Z);
    free(XC);
    free(ZC);
    printf("Velocity and pressure were initialised\n");
    printf("Boundary conditions were set up\n");
    
}
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


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
    printf("SETTING PURE SHEAR BC\n");
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
    printf("SETTING SIMPLE SHEAR BC\n");
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
                  .SetBCVx = SetPureOrSimpleShearBCVx,
                  .SetBCVz = SetPureOrSimpleShearBCVz,
          },
  };
  RunMDOODZ("ViscousInclusionTestSuite.txt", &instance);
}