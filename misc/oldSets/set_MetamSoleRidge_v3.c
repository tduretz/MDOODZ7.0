// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2018  MDOODZ Developper team
//
// This file is part of MDOODZ.
//
// MDOODZ is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// MDOODZ is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with MDOODZ.  If not, see <http://www.gnu.org/licenses/>.
// =========================================================================

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "mdoodz-private.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


double myErfInv2(double x){
   double tt1, tt2, lnx, sgn;
   sgn = (x < 0) ? -1.0f : 1.0f;

   x = (1 - x)*(1 + x);        // x = 1 - x*x;
   lnx = logf(x);

   tt1 = 2/(PI*0.147) + 0.5f * lnx;
   tt2 = 1/(0.147) * lnx;

   return(sgn*sqrtf(-tt1 + sqrtf(tt1*tt1 - tt2)));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid Mmesh, scale scaling ) {
    
    int k;
    double TopoLevel = 0.0e3/scaling.L;
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel;
        topo_chain->phase[k] = 0;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
    int np;
    FILE *read;
    int s1, s2;
    
    s1 = sizeof(int);
    s2 = sizeof(double);
    
    double H,xl,L, norm;
    int    k,nb_layers = 20;
    
    H = (model.zmax - model.zmin);
    L = (model.xmax - model.xmin);
    
    double gsbg = 2e-3/scaling.L;
    double Temperature = (700.0+zeroC)/scaling.T;
    double red_left = (double)(model.user0);
    double z_ast = (double)(model.user3)/scaling.L;
    double z_sed = (double)(model.user4)/scaling.L;
    
    double AspectRatio = 5.0;
    double x_ell, z_ell, a_ell = model.user1/scaling.L, b_ell = model.user1/scaling.L * AspectRatio, angle = 70*M_PI/180;
    double x0 = 0.0*350e3/scaling.L;
    double z0 = -5e3/scaling.L;
    double z_uc   = (double)(model.user1)/scaling.L;
    double z_lc   = (double)(model.user2)/scaling.L;
    double Hserp = 1e3/scaling.L;
    double Hstep = 1e3/scaling.L;
    double VplateW = (double)(model.user5)/scaling.V;
    double VplateE = (double)(model.user6)/scaling.V;
    double Kappa   = 1e-6/( (scaling.L*scaling.L)/scaling.t );
    
    double coeff   = -1.9439; // From inverse_erf to design base of lithosphere
    double a_surf  = -2821/scaling.L; // ocean bottem depth for 0 age
    double b_surf  = -5.6074e-05/(scaling.L/pow(scaling.t,0.5));
    double z_ob    = 0.0; // ocean bottom depth
    double z_uc1, z_lc1, z_sed1, z_ast1;
    double Tmant   = (1350+zeroC)/scaling.T, age, R, infsum, betj, c, jd;
    double Tlith   = (1300+zeroC)/scaling.T;
    double Tsurf   = zeroC/scaling.T;
    double a       = 125e3/scaling.L; // for plate model
    int thermal_model = 0, j;
    
//    double coeff   = myErfInv2( coeff )
    
    for( np=0; np<particles->Nb_part; np++ ) {
        
        particles->d[np]     = gsbg;                                            // same grain size everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->T[np]     = Temperature;
        particles->phase[np] = 2;
        
        //----------------------------------
        // INCLUSIONS INCLUSIONS INCLUSIONS
        //----------------------------------
        
        if (particles->x[np] <= 0.0 ) {
            z_ast = 2.0*sqrt( Kappa*fabs(particles->x[np]) / VplateW ) * coeff;
            z_ob  = a_surf + b_surf*pow( fabs(particles->x[np]) / VplateW, 0.5);
            age   = fabs(particles->x[np]) / VplateW;
            R     = VplateW*a/2.0/Kappa;
        }
               
        if (particles->x[np] > 0.0 ) {
            z_ast = 2.0*sqrt( Kappa*fabs(particles->x[np]) / VplateE ) * coeff;
            z_ob  = a_surf + b_surf*pow( fabs(particles->x[np]) / VplateE, 0.5);
            age   = fabs(particles->x[np]) / VplateE;
            R     = VplateW*a/2.0/Kappa;
        }
        
        if (thermal_model==0) {
            particles->T[np] =  (Tmant-Tsurf)*erf(-(particles->z[np]-z_ob) / sqrt(4.0 * Kappa * age) ) + Tsurf;
            if (particles->T[np]<Tsurf) particles->T[np] = Tsurf;
        }
        else {
            infsum = 0.0;
            for( j=1; j<10; j++ ) {
                jd = (double)j;
                betj   = sqrt(R*R + jd*jd*M_PI*M_PI) - R;

                c      = 2.0/jd/M_PI;
                infsum = infsum + c*exp(-betj*particles->x[np]/a) * sin(jd * M_PI * -particles->z[np]/a);
            }
            particles->T[np]   = (Tmant-zeroC/scaling.T) *(-particles->z[np]/a + infsum);
        }
//        z_ob = -5000/scaling.L;
        
//        printf("%2.2e\n", z_ob*scaling.L);
//
        z_uc1  = z_uc + z_ob;
        z_ast1 = z_ast + z_ob;
        z_lc1  = z_lc + z_ob;
        z_sed1 = z_sed + z_ob;

        // water
        if ( particles->z[np] > z_ob  ) particles->phase[np] = 0;
        
        // Basalt layer
        if ( particles->z[np] > z_uc1 && particles->z[np] < z_ob && particles->x[np] > 0.0) {
            particles->phase[np] = 4;
        }
        // Gabbro layer
        if ( particles->z[np] > z_lc1 && particles->z[np] < z_uc1 && particles->x[np] > 0.0 ) {
            particles->phase[np] = 5;
        }
        
        // Basalt layer
        if ( particles->z[np] > z_uc1*red_left && particles->z[np] < z_ob && particles->x[np] <= 0.0) {
            particles->phase[np] = 4;
        }
        // Gabbro layer
        if ( particles->z[np] > z_lc1*red_left && particles->z[np] < z_uc1*red_left && particles->x[np] <= 0.0 ) {
            particles->phase[np] = 5;
        }
        
//        // Serpentine
//        if ( particles->z[np] > -(Huc + Hlc + Hserp) && particles->z[np] < -(Huc + Hlc) ) {
//            particles->phase[np] = 6;
//        }

        // Astenosphere
//        if ( particles->z[np] < z_ast1 ) {
//            particles->phase[np] = 8;
//        }
        
        // Astenosphere
        if ( particles->T[np] > Tlith ) particles->phase[np] = 8;
        
//        // Step
//        if ( particles->x[np] > 0.0 && particles->z[np] < (z_ast+Hstep) ) {
//            particles->phase[np] = 8;
//        }
        
        // Sediments
        if ( particles->z[np] > z_sed1 && particles->z[np] < z_ob ) {
            particles->phase[np] = 7;
        }
        
        //------------------------------------------//
//        if ( layer == 1 ) {
//            // Central layer
//            if ( pow(particles->z[np],2) < pow(rad,2) ) {
//                particles->phase[np] = 0;
//            }
//        }
//        else {
            
//            x_ell = (particles->x[np]-x0)*cos(angle) + (particles->z[np]-z0)*sin(angle);
//            z_ell =-(particles->x[np]-x0)*sin(angle) + (particles->z[np]-z0)*cos(angle);
//            if (pow(x_ell/a_ell,2) + pow(z_ell/b_ell,2) < 1) particles->phase[np] = 6;

//        }

        //----------------------------------
        // STRIPES STRIPES STRIPES STRIPES
        //----------------------------------
        // Horizontal Layering
        if (particles->phase[np]==2 && particles->z[np]> 1.0*H/6.0 && particles->z[np]< 2.0*H/6.0){
            //if (particles->x[np]> 1.0*H/6.0 && particles->x[np]< 2.0*H/6.0)
            particles->phase[np] = 3;
        }
        if (particles->phase[np]==2 && particles->z[np]< 0.0*H/6.0 && particles->z[np]> -1.0*H/6.0){
            particles->phase[np] = 3;
        }
        if (particles->phase[np]==2 && particles->z[np]< -2.0*H/6.0 && particles->z[np]> -3.0*H/6.0){
            particles->phase[np] = 3;
        }
        //        // Vertical alternance:
        if (particles->x[np]> 1.0*L/6.0 && particles->x[np]< 2.0*L/6.0){
            if (particles->phase[np] == 3) particles->phase[np] = 30;
            if (particles->phase[np] == 2) particles->phase[np] = 20;
        }
        if (particles->x[np]> -1.0*L/6.0 && particles->x[np]< 0.0*L/6.0){
            if (particles->phase[np] == 3) particles->phase[np] = 30;
            if (particles->phase[np] == 2) particles->phase[np] = 20;
        }
        
        if (particles->x[np]> -3.0*L/6.0 && particles->x[np]< -2.0*L/6.0){
            if (particles->phase[np] == 3) particles->phase[np] = 30;
            if (particles->phase[np] == 2) particles->phase[np] = 20;
        }
        
        if (particles->phase[np] == 20) particles->phase[np] = 3;
        if (particles->phase[np] == 30) particles->phase[np] = 2;
        
        
        // here we define the initial particules velocity that should comply with the background homogeneous shear conditions
        particles->Vx[np]    = -0.5*model.bkg_strain_rate*(particles->z[np] + (model.zmax-model.zmin)/2);
        particles->Vz[np]    = 0;
        particles->rho[np]   = materials->rho[particles->phase[np]];
        
        
        //--------------------------//
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }
        //--------------------------//
        
    }
    
    MinMaxArray(particles->Vx,  scaling.V, particles->Nb_part,   "Vxp init" );
    MinMaxArray(particles->Vz,  scaling.V, particles->Nb_part,   "Vzp init" );
    MinMaxArray(particles->rho, scaling.rho, particles->Nb_part, "Rho init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "T init" );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {
    
    int   kk, k, l, c, c1;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double dmin, VzBC, width = 1 / scaling.L, eta = 1e4 / scaling.eta ;
    double Lx, Lz, T1, T2, rate=model->bkg_strain_rate,  z_comp=-140e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    double Vfix = (50.0/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [50.0 == 5 cm/yr]
    double shift = model->user8;
    
    // Thermal crazyness
    int    phase_ast = 8;
    double k_crazy = 1000*materials->k[phase_ast];
    
    if (model->step == 0) {
        materials->k_eff[phase_ast] = k_crazy;
        printf("Running with crazy conductivity for the asthenosphere!!\n");
    }
    
    else {
        materials->k_eff[phase_ast] = materials->k[phase_ast];
        printf("Running with normal conductivity for the asthenosphere...\n");
    }
    
    // Define dimensions;
    Lx = (double) (model->xmax - model->xmin) ;
    Lz = (double) (model->zmax - model->zmin) ;
        
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
    
    // Change particle color if in field of interest
    double Pmin = 0.6e9/scaling.S;
    double Pmax = 1.25e9/scaling.S;
    double Tmin = (700.0+zeroC)/scaling.T;
    double Tmax = (900.0+zeroC)/scaling.T;
    double TWaterSeds = (50.0+zeroC)/scaling.T;
    
    if (model->step>0 && (int)model->user7==1 ) {

    #pragma omp parallel for shared ( particles ) private ( k ) firstprivate ( Pmin, Pmax, Tmin, Tmax )
        for (k=0;k<particles->Nb_part;k++) {
            
            if ( particles->T[k]>=TWaterSeds ) {
                if ( particles->phase[k] == 0 ) particles->phase[k] = 7;
            }

            if ( particles->T[k]>=Tmin && particles->T[k]<=Tmax && particles->P[k]>=Pmin && particles->P[k]<=Pmax ) {

                if ( particles->phase[k] == 4 ) particles->phase[k] = 9;
                if ( particles->phase[k] == 5 ) particles->phase[k] = 10;
                if ( particles->phase[k] == 7 ) particles->phase[k] = 11;

            }
        }
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
    
    int show         = 0;
    double inflow_EW = 0.0;
    
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
//                        printf("shift*Lx/2.0 = %2.2e\n", shift*Lx/2.0*scaling.L);
                        mesh->BCu.val[c]  = -(mesh->xg_coord[k]+shift*Lx/2.0) * model->bkg_strain_rate;
                        inflow_EW  += fabs(mesh->BCu.val[c])*model->dz;
                        if (show==0) printf("Vx = %2.2e\n", mesh->BCu.val[c]*scaling.V);
                        if (show==0) show=1;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = -(mesh->xg_coord[k]+shift*Lx/2.0) * model->bkg_strain_rate;
                        inflow_EW  += fabs(mesh->BCu.val[c])*model->dz;
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
                        mesh->BCu.val[c]  = 0.0*model->bkg_strain_rate*Lx;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] =  -12;
                        mesh->BCu.val[c]  = -0.0*model->bkg_strain_rate*Lx;
                    }
                    
                    // Free slip S
                    if (l==0 ) { //&& (k>0 && k<NX-1) ) {
                        mesh->BCu.type[c] =  11;
                        mesh->BCu.val[c]  = -1*model->bkg_strain_rate*Lz;
                    }
                    
                    // Free slip N
                    if ( l==mesh->Nz) {// && (k>0 && k<NX-1)) {
                        mesh->BCu.type[c] =  11;
                        mesh->BCu.val[c]  =  1*model->bkg_strain_rate*Lz;
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
    
    double Vz_outflow_south = -inflow_EW/Lx;
    show =0;
    
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
                        mesh->BCv.val[c]  = Vz_outflow_south;////mesh->zg_coord[l] * model->bkg_strain_rate;
                        
                        if (show==0) printf("%2.2e %2.2e\n",  mesh->zg_coord[l] * model->bkg_strain_rate * scaling.V, Vz_outflow_south* scaling.V);
                        if (show==0) show=1;

                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0;//mesh->zg_coord[l] * model->bkg_strain_rate;
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
                        mesh->BCv.val[c]  = -0.0*model->bkg_strain_rate*Lz;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0*model->bkg_strain_rate*Lz;
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
    double Tbot = (1300 + 273.15)/scaling.T;
    double Tleft, Tright;
    
    
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
                
                // WEST
                if ( k==0 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typW[l] = 0;
                    mesh->BCt.valW[l] = 0;
                }
                
                // EAST
                if ( k==NCX-1 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typE[l] = 0;
                    mesh->BCt.valE[l] = 0;
                }
                
                // SOUTH
                if ( l==0 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typS[k] = 1;
                    mesh->BCt.valS[k] = Tbot;
                }
                
                // NORTH
                if ( l==NCZ-1 ) {
                    mesh->BCt.type[c] = 0;
                    mesh->BCt.typN[k] = 1;
                    mesh->BCt.valN[k] = Ttop;
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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
