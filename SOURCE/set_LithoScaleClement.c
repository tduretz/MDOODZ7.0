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
#include "header_MDOODZ.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SIMPLE
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = 2.e-1/scaling.L; // sets zero initial topography
    double sig = 0.5/scaling.L;
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = 0;//TopoLevel * exp(-topo_chain->x[k]*topo_chain->x[k]/sig/sig) ;
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
    
    // Define dimensions;
    double Lz = (double) (model.zmax - model.zmin) ;
    double HLit   = model.user1/scaling.L;
    double Hcrust = -60e3/scaling.L;
    double Tsurf  = 273.15/scaling.T, Tpart;
    double Htot   = Lz, Tmant = (model.user0 + 273.15)/scaling.T;
    double WZHW   = 10e3/scaling.L;
    double angle  = model.user4*M_PI/180.0;
//    double zlc = -35e3/scaling.L;
//    double zmc = -25e3/scaling.L;
//    double zuc = -15e3/scaling.L;
    double x_ell, z_ell, a_ell = 10e3/scaling.L, b_ell = 2e3/scaling.L;
    double x0 = 0.0*350e3/scaling.L;
    double z0 = -35e3/scaling.L;
    int il; //Number of layers in the crust
    
    double spacing = 5.0e3/scaling.L; // espacement layering crust
    
    //angle = 00*M_PI/180;
    
    // Loop on particles
    for ( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
        particles->Vx[np]    = -1.0*particles->x[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->Vz[np]    =      particles->z[np]*model.EpsBG;               // set initial particle velocity (unused)
        particles->phase[np] = 2;                                               // same phase number everywhere
        particles->phi[np]   = 0.0;                                             // zero porosity everywhere
        particles->X[np]     = 0.0;                                             // zero X everywhere
        
        
        //--------------------------//
        // TEMPERATURE
        Tpart = ((Tmant-Tsurf)/Htot)*(-particles->z[np]-0e3/scaling.L) + Tsurf; // linear geotherm
        if (Tpart>Tmant) Tpart = Tmant;                                         // correct if geotherm above Tmant
        particles->T[np]       = Tpart;                                         // assign temperature to particles
        //--------------------------//
        
        // If particles are in the lithosphere: change phase
        if (particles->z[np]>-HLit) particles->phase[np] = 0;
        
        // If particles are in the crust: change phase
        if (particles->z[np] > Hcrust)  particles->phase[np] = 3;
        
        // layering crust:
        if (particles->phase[np]==3){
            for( il=0; il<= 50; il=il+2 ) {
                if (particles->z[np]>= Hcrust && particles->z[np]>=-((il+1)*spacing) && particles->z[np]<-(il*spacing)) particles->phase[np] = 4;
            }
        }
        
        
        // inverse layering for checking streching/shortening
        
        if (particles->x[np] < -150.0e3/scaling.L && particles->phase[np]==3) particles->phase[np] = 11;
        if (particles->x[np] < -150.0e3/scaling.L && particles->phase[np]==4) particles->phase[np] = 10;
        if (particles->x[np] > +150.0e3/scaling.L && particles->phase[np]==3) particles->phase[np] = 11;
        if (particles->x[np] > +150.0e3/scaling.L && particles->phase[np]==4) particles->phase[np] = 10;
        
        if (particles->x[np] < -50.0e3/scaling.L && particles->phase[np]==3) particles->phase[np] = 10;
        if (particles->x[np] < -50.0e3/scaling.L && particles->phase[np]==4) particles->phase[np] = 11;
        if (particles->x[np] > +50.0e3/scaling.L && particles->phase[np]==3) particles->phase[np] = 10;
        if (particles->x[np] > +50.0e3/scaling.L && particles->phase[np]==4) particles->phase[np] = 11;
        
        if (particles->phase[np]==10) particles->phase[np] = 4;
        if (particles->phase[np]==11) particles->phase[np] = 3;
        

//        if (particles->z[np]>zlc) particles->phase[np] = 5;
//        if (particles->z[np]>zmc) particles->phase[np] = 4;
//        if (particles->z[np]>zuc) particles->phase[np] = 3;
        
        
//        x_ell = (particles->x[np]-x0)*cos(angle) + (particles->z[np]-z0)*sin(angle);
//        z_ell =-(particles->x[np]-x0)*sin(angle) + (particles->z[np]-z0)*cos(angle);
//        if (pow(x_ell/a_ell,2) + pow(z_ell/b_ell,2) < 1) particles->phase[np] = 1;
        
        
        // Weak zone
        //        if (particles->z[np]>-HLit && particles->x[np]>-WZHW + particles->z[np]*angle && particles->x[np]<WZHW + particles->z[np]*angle) particles->phase[np] = 1;
        
        
        //        // Crust layer
        //        if (particles->z[np]>-Hcrust) {
        //            particles->phase[np] = 1;
        ////            particles->X[np]     = 1.0;
        //        }
        
        particles->d[np]     = materials->gs_ref[particles->phase[np]];                                            // same grain size everywhere
        
        
        //        //--------------------------//
        //        // DENSITY
        //        if ( model.eqn_state > 0 ) {
        //            particles->rho[np] = materials->rho[particles->phase[np]] * (1 -  materials->alp[particles->phase[np]] * (Tpart - materials->T0[particles->phase[np]]) );
        //        }
        //        else {
        //            particles->rho[np] = materials->rho[particles->phase[np]];
        //        }
        //        //--------------------------//
    }
    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "Tp init" );
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
    double Lx, Lz, T1, T2, rate=model->EpsBG,  z_comp=-140e3/scaling.L;
    double Vx_r, Vx_l, Vz_b, Vz_t, Vx_tot, Vz_tot;
    double Lxinit = 1400e3/scaling.L, ShortSwitchV0 = 0.40;
    double Vfix = (50.0/(1000.0*365.25*24.0*3600.0))/(scaling.L/scaling.t); // [50.0 == 5 cm/yr]
    double k_crazy = 1000*materials->k[2];
    int    cctot=0, ccbcIN=0, ccbcOUT;
    double HLit  = model->user1/scaling.L;
    double VxIn  = model->user2/scaling.V, VxOut=0.0;
    double zlit_fact = 1.0;
    double influx, rat=model->user3;
    double Ttop = 273.15/scaling.T;
    double Tbot = (model->user0 + 273.15)/scaling.T;

    // Fix temperature
    double Tgrad = model->user3/(scaling.T/scaling.L);
    double Tsurf = 293.0/(scaling.T);
    
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    
    if (model->step == 0) {
        materials->k_eff[2] = k_crazy;
        printf("Running with crazy conductivity for the asthenosphere!!\n");
    }
    
    else {
        materials->k_eff[2] = materials->k[2];
        printf("Running with normal conductivity for the asthenosphere...\n");
    }
    
    //    // Define the plate detachement at the right side (Annelore 24.02.2017)
    //    if (fabs(Vcut)<10e-20) {
    //
    //        // Loop on particles
    //        for( np=0; np<particles->Nb_part; np++ ) {
    //
    //            // Plate detachement
    //            if (particles->x[np]>Xcutmin && particles->phase[np]==0) particles->phase[np] = 1;
    //            if (np==0) printf("Change particle numbers for the free subduction is OK! \n");
    //        }
    //    }
    
    // Define depth at which the velocity magnitude changes
    for (l=0; l<mesh->Nz+1; l++) {
        for (k=0; k<mesh->Nx; k++) {
            
            c = k + l*(mesh->Nx);
            if (k==mesh->Nx-1 && (l>0 && l<mesh->Nz)) {
                if ( mesh->BCu.type[c] != 30 ) {
                    cctot++;
                    if (mesh->zvx_coord[l]>-HLit*zlit_fact) ccbcIN++;
                }
            }
        }
    }
    // Calculate influx/outflux
    ccbcOUT = cctot - ccbcIN;
    influx  = ccbcIN*model->dz*VxIn;
    VxOut   = -influx/model->dz/(ccbcOUT+1);
    printf ("Number of cells on East side = %d\nNumber of cells with inflow  = %d\nNumber of cells with outflow = %d\n", cctot, ccbcIN, ccbcOUT);
    //        VxOut = -(ccbc*VxIn) / (cctot-ccbc-1);
    printf("Vx Outflow         = %2.2e Vx Inflow = %2.2e\n", VxOut*scaling.V,  VxIn*scaling.V);
    printf("Total influx       = %2.2e\n", influx*scaling.V*scaling.L);
    
    
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
    
    for(kk=0; kk<1; kk++) {
        NX  = mesh->Nx;
        NZ  = mesh->Nz;
        NCX = NX-1;
        NCZ = NZ-1;
        
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
                    // FREE SURFACE - Dirchlet on the free surface
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = Ttop;
                        }
                    }
                    
                    
                }
                
            }
        }
        
    }

    printf("Velocity and pressure were initialised\n");
    printf("Boundary conditions were set up\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

