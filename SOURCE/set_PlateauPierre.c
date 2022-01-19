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
    double TopoLevel = 0e3/scaling.L; // sets zero initial topography
    double x1=-25e3/scaling.L, x2 = x1 + model.user4/scaling.L;
    double dTopo = 1e3/scaling.L;

    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel;
        topo_chain->phase[k] = 0;
        
        if (topo_chain->x[k] < x2 && topo_chain->x[k] > x1) topo_chain->z[k] = TopoLevel + dTopo/(x2-x1)*(topo_chain->x[k]-x1);
        
        if (topo_chain->x[k] > x2) topo_chain->z[k] = TopoLevel + dTopo;
        if (topo_chain->x[k] < x1) topo_chain->z[k] = TopoLevel ;
    }
    
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
    int np;
    double HLit   = model.user0/scaling.L;                           // lithosphere thickness
    double HCrust = 49e3/scaling.L;                           // crust thickness
    double LCrust = HCrust - model.user1/scaling.L;                           // crust thickness
//    double LCrust = 25e3/scaling.L;                           // crust thickness
    double dHCrust = -model.user2/scaling.L;
    //double dHCrust = -28e3/scaling.L;                           // crust thickness
    double Tsurf  = 273.15/scaling.T, Tpart;                  // surface temperature
    double Tmant  = (model.user3+273.15)/scaling.T;                  // adiabatic mantle temperature
    double rad    = model.user0/scaling.L;
    double x1=-25e3/scaling.L, temp, temp2;
    double x2 = x1 + model.user4/scaling.L;
    double x3 = -750e3/scaling.L;
    double x4 = 750e3/scaling.L;
    double dTopo = 5e3/scaling.L;
    double dLab = 40e3/scaling.L;

    double zW  = model.user6/scaling.L;
    double dip = atan(model.user7*M_PI/180.0);
    double Vel = model.user8/scaling.V;
    
    double ztop, zbot, Hslab = 100e3/scaling.L;
    
//    int kinematics = (int)model.user9;
    
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
        particles->Vx[np]    = -particles->x[np]*model.EpsBG; // set initial particle velocity (unused)
        particles->Vz[np]    =  particles->z[np]*model.EpsBG; // set initial particle velocity (unused)
        particles->phase[np] = 0;                             // same phase number everywhere
        particles->phi[np]   = 0.0;                           // zero porosity everywhere
        
        //--------------------------//
        // TEMPERATURE - linear geotherm
        Tpart = ((Tmant-Tsurf)/HLit)*(-particles->z[np]) + Tsurf;
        if (Tpart>Tmant) Tpart = Tmant;
        particles->T[np]       = Tpart;
        
        //--------------------------//
        // Phases - lithosphere-asthenosphere
        if (particles->x[np] < x1 && particles->x[np] > x3) {
            temp2 = -HLit - dLab/(x4-x3)*(particles->x[np]-x3);
            if (particles->z[np] <(-LCrust))  particles->phase[np] = 3;
            if (particles->z[np] <-(HCrust))  particles->phase[np] = 1;
            if (particles->z[np] <temp2  )  particles->phase[np] = 2;
        }
        
        if (particles->x[np] < x2 && particles->x[np] > x1) {
           temp = -HCrust + dHCrust/(x2-x1)*(particles->x[np]-x1);
//            printf("%2.2e %2.2e\n", temp*scaling.L/1000, particles->x[np]*scaling.L/1000);
            temp2 = -(HLit+dLab/(x4-x3)*(x1-x3)) - dLab/(x4-x3)*(particles->x[np]-x1);
            
            //if (particles->z[np] >temp     )  particles->phase[np] = 0;
            if (particles->z[np] >temp && particles->z[np] < temp + model.user1/scaling.L      )  particles->phase[np] = 3;
            if (particles->z[np] <temp      )  particles->phase[np] = 1;
            if (particles->z[np] <temp2  )  particles->phase[np] = 2;
        }
        
        if (particles->x[np] > x2) {
            //if (particles->z[np] >(-HCrust+dHCrust))  particles->phase[np] = 0;
            temp2 = -(HLit+dLab/(x4-x3)*(x2-x3)) - dLab/(x4-x3)*(particles->x[np]-x2);
            particles->phase[np] = 5;
            if (particles->z[np] <(-LCrust+dHCrust))  particles->phase[np] = 6;
            if (particles->z[np] <(-HCrust+dHCrust))  particles->phase[np] = 1;
            if (particles->z[np] <temp2  )  particles->phase[np] = 2;

        }
        
        // Draw fake slab
        ztop = zW - dip*(particles->x[np] -model.xmin);
        zbot = zW-Hslab - dip*(particles->x[np] -model.xmin);
        if (particles->z[np] < ztop && particles->z[np] > zbot) {
//        if (particles->z[np] < ztop ) {
            particles->phase[np] = 4;
        }
        
        // Set Indian upper crust
        //if (particles->x[np] > x1 && particles->phase[np] == 0) particles->phase[np] = 4;
        
        // Set Indian middle crust
//        if (particles->phase[np] == 4 && particles->z[np] <(-(HCrust-LCrust+dTopo) ) ) particles->phase[np] = 6;
        //if (particles->phase[np] == 4 && particles->z[np] <(-(LCrust) ) ) particles->phase[np] = 6;

        
        
        // Set Indian lower crust
        //if (particles->x[np] > (x2+100e3/scaling.L) && particles->phase[np] == 3) particles->phase[np] = 5;
        
        
        //--------------------------//
        // DENSITY
        if ( model.eqn_state > 0 ) {
            particles->rho[np] = materials->rho[particles->phase[np]] * (1 -  materials->alp[particles->phase[np]] * (Tpart - materials->T0[particles->phase[np]]) );
        }
        else {
            particles->rho[np] = materials->rho[particles->phase[np]];
        }
        
        //--------------------------//
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }
        //--------------------------//
    }

    MinMaxArray(particles->Vx, scaling.V, particles->Nb_part, "Vxp init" );
    MinMaxArray(particles->Vz, scaling.V, particles->Nb_part, "Vzp init" );
    MinMaxArray(particles->T, scaling.T, particles->Nb_part,  "Tp init" );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Set physical properties on the grid and boundary conditions
void SetBCs( grid *mesh, params *model, scale scaling, markers* particles, mat_prop *materials, surface* topo ) {
    
    int   kk, k, l, c, c1, np;
    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    double TN =  273.15/scaling.T, TS = (1330.+273.15)/scaling.T;
    double TW = (1330.+273.15)/scaling.T, TE = (1330.+273.15)/scaling.T;
    double Tbot, Tleft, Tright;

    double xW = -100e3/scaling.L;
    double xE =  100e3/scaling.L;
    double zS = -600e3/scaling.L;
    double zN = -500e3/scaling.L;
    double VxBC =  1.5e-9/scaling.V;
    double VzBC = -1.5e-9/scaling.V;
    
    double zW  = model->user6/scaling.L;
    double dip = atan(model->user7*M_PI/180.0);
    double Vel = model->user8/scaling.V;
    
    // Linear gradient of velocity along West boundary
    double zlabW = -model->user0/scaling.L;
    double dVdzW = Vel / (zW-zlabW);
    
    // Linear gradient of velocity along East boundary
    double dLab  = 40e3/scaling.L;
    double zE    = zW - dip*(model->xmax-model->xmin);
    double zlabE = zlabW - dLab;
    double dVdzE = Vel / (zE-zlabE);
    double V;
    
    double Hslab = 100e3/scaling.L;
    double ztop, zbot;
    
    // super idée de Ben: utiliser des demi gaussiennes au nord/sud du slab sur une epaisseur donnée Hsig
    double sigma   = 10e3/scaling.L;
    double Hsigma  = 4.0*sigma;
    
    
    // Set unresaonable conductivity in the mantle to generate and adiabatic mantle as initial condition
    int    phase_ast1 = 2;                     // This has to be the ASTHENOSPHERE phase number
    int    phase_ast2 = 4;
    double k_crazy = 1000*materials->k[phase_ast1];

    if (model->step == 0) {
        materials->k_eff[phase_ast1] = k_crazy;
        materials->k_eff[phase_ast2] = k_crazy;
        printf("Running with crazy conductivity for the asthenosphere!!\n");
    }
    else {
        materials->k_eff[phase_ast1] = materials->k[phase_ast1];
        materials->k_eff[phase_ast2] = materials->k[phase_ast2];
        printf("Running with normal conductivity for the asthenosphere...\n");
    }
 
    // Some stuff...
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
                    
                    // Matching BC nodes WEST
                    if (k==0 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = 0.0;
//                        mesh->BCu.val[c]  = -mesh->xg_coord[k] * model->EpsBG;
//                        if (mesh->zvx_coord[l] < zlabW && mesh->zvx_coord[l] > zW  ) {
//                            V                 = dVdzW * (mesh->zvx_coord[l] - zlabW);
//                            mesh->BCu.val[c]  = V*cos(dip);
//                        }
                        
                        // Dans la limite dessus
                        if (mesh->zvx_coord[l] > zW && mesh->zvx_coord[l] < zW+Hsigma  ) {
                            mesh->BCu.val[c]  = Vel*cos(dip) * exp(-pow( (mesh->zvx_coord[l] - zW), 2) / (2.0*sigma*sigma) );
                        }
                        
                        // On est dans la steppe
                        if (mesh->zvx_coord[l] < zW && mesh->zvx_coord[l] > zW-Hslab  ) {
                            mesh->BCu.val[c]  = Vel*cos(dip);
                        }
                        
                        // Dans la limite dessous
                        if (mesh->zvx_coord[l] > zW-Hslab-Hsigma  && mesh->zvx_coord[l] < zW-Hslab  ) {
                            mesh->BCu.val[c]  = Vel*cos(dip) * exp(-pow( (mesh->zvx_coord[l] - (zW-Hslab)), 2) / (2.0*sigma*sigma) );
                        }
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = 0.0;
                        
                        // Dans la limite dessus
                        if (mesh->zvx_coord[l] > zE && mesh->zvx_coord[l] < zE+Hsigma  ) {
                            mesh->BCu.val[c]  = Vel*cos(dip) * exp(-pow( (mesh->zvx_coord[l] - zE), 2) / (2.0*sigma*sigma) );
                        }
                        
                        // On est dans la steppe
                        if (mesh->zvx_coord[l] < zE && mesh->zvx_coord[l] > zE-Hslab  ) {
                            mesh->BCu.val[c]  = Vel*cos(dip);
                        }
                        
                        // Dans la limite dessous
                        if (mesh->zvx_coord[l] > zE-Hslab-Hsigma  && mesh->zvx_coord[l] < zE-Hslab  ) {
                            mesh->BCu.val[c]  = Vel*cos(dip) * exp(-pow( (mesh->zvx_coord[l] - (zE-Hslab)), 2) / (2.0*sigma*sigma) );
                        }
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
                    
                    // !!!!!! Hand-of-god !!!!!!

                    // Draw fake slab
                    ztop =  zW        - dip*(mesh->xg_coord[k] - model->xmin);
                    zbot = (zW-Hslab) - dip*(mesh->xg_coord[k] - model->xmin);
                    
                    // Dans la limite dessus
                    if (mesh->zvx_coord[l] > ztop && mesh->zvx_coord[l] < ztop+Hsigma  ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = Vel*cos(dip) * exp(-pow( (mesh->zvx_coord[l] - ztop), 2) / (2.0*sigma*sigma) );
                    }
                    
                    // On est dans la steppe
                    if ( mesh->zvx_coord[l] < ztop && mesh->zvx_coord[l] > ztop-Hslab ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = Vel*cos(dip);
                    }
                    
                    // Dans la limite dessous
                    if (mesh->zvx_coord[l] > ztop-Hslab-Hsigma  && mesh->zvx_coord[l] < ztop-Hslab  ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  = Vel*cos(dip) * exp(-pow( (mesh->zvx_coord[l] - (ztop-Hslab)), 2) / (2.0*sigma*sigma) );
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
                        mesh->BCv.val[c]  = 0.0;//mesh->zg_coord[l] * model->EpsBG;
                        mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0;//mesh->zg_coord[l] * model->EpsBG;
                                                mesh->BCv.val[c]  = mesh->zg_coord[l] * model->EpsBG;
                    }
                    
                    // Non-matching boundary WEST
                    if ( (k==0) ) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0.0;
//                        if (mesh->zg_coord[l] < zlabW && mesh->zg_coord[l] > zW  ) {
//                            mesh->BCv.type[c] =   11;
//                            V                 = dVdzW * (mesh->zg_coord[l] - zlabW);
//                            mesh->BCv.val[c]  = -V*sin(dip);
//                        }
                        
                        // Dans la limite dessus
                        if (mesh->zg_coord[l] > zW && mesh->zg_coord[l] < zW+Hsigma  ) {
                            mesh->BCv.type[c] =   11;
                            mesh->BCv.val[c]  = -Vel*sin(dip) * exp(-pow( (mesh->zg_coord[l] - zW), 2) / (2.0*sigma*sigma) );
                        }
                        
                        // On est dans la steppe
                        if (mesh->zg_coord[l] < zW && mesh->zg_coord[l] > zW-Hslab  ) {
                            mesh->BCv.type[c] =   11;
                            mesh->BCv.val[c]  = -Vel*sin(dip);
                        }
                        
                        // Dans la limite dessous
                        if (mesh->zg_coord[l] > zW-Hslab-Hsigma  && mesh->zg_coord[l] < zW-Hslab  ) {
                            mesh->BCv.type[c] =   11;
                            mesh->BCv.val[c]  = -Vel*sin(dip) * exp(-pow( (mesh->zg_coord[l] - (zW-Hslab)), 2) / (2.0*sigma*sigma) );
                        }
                    }
                    
                    // Non-matching boundary EAST
                    if ( (k==mesh->Nx) ) {
                        mesh->BCv.type[c] =   13;
                        mesh->BCv.val[c]  =   0.0;
//                        if (mesh->zg_coord[l] < zlabE && mesh->zg_coord[l] > zE  ) {
//                            mesh->BCv.type[c] =   11;
//                            V                 = dVdzE * (mesh->zg_coord[l] - zlabE);
//                            mesh->BCv.val[c]  = -V*sin(dip);
//                        }
                        // Dans la limite dessus
                        if (mesh->zg_coord[l] > zE && mesh->zg_coord[l] < zE+Hsigma  ) {
                            mesh->BCv.type[c] =   11;
                            mesh->BCv.val[c]  = -Vel*sin(dip) * exp(-pow( (mesh->zg_coord[l] - zE), 2) / (2.0*sigma*sigma) );
                        }
                        
                        // On est dans la steppe
                        if (mesh->zg_coord[l] < zE && mesh->zg_coord[l] > zE-Hslab  ) {
                            mesh->BCv.type[c] =   11;
                            mesh->BCv.val[c]  = -Vel*sin(dip);
                        }
                        
                        // Dans la limite dessous
                        if (mesh->zg_coord[l] > zE-Hslab-Hsigma  && mesh->zg_coord[l] < zE-Hslab  ) {
                            mesh->BCv.type[c] =   11;
                            mesh->BCv.val[c]  = -Vel*sin(dip) * exp(-pow( (mesh->zg_coord[l] - (zE-Hslab)), 2) / (2.0*sigma*sigma) );
                        }
                    }
                    
                    // !!!!!! Hand-of-god !!!!!!
                    
                    // Draw fake slab
                    ztop =  zW        - dip*(mesh->xvz_coord[k] - model->xmin);
                    zbot = (zW-Hslab) - dip*(mesh->xvz_coord[k] - model->xmin);
                    
//                    // Dans la limite dessus
//                    if (mesh->zg_coord[l] > ztop && mesh->zg_coord[l] < ztop+Hsigma  ) {
//                        mesh->BCv.type[c] =   0;
//                        mesh->BCv.val[c]  = -Vel*sin(dip) * exp(-pow( (mesh->zg_coord[l] - ztop), 2) / (2.0*sigma*sigma) );
//                    }
                    
                    // On est dans la steppe
                    if ( mesh->zg_coord[l] < ztop && mesh->zg_coord[l] > ztop-Hslab ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = -Vel*sin(dip);
                    }
                    
//                    // Dans la limite dessous
//                    if (mesh->zg_coord[l] > ztop-Hslab-Hsigma  && mesh->zg_coord[l] < ztop-Hslab  ) {
//                        mesh->BCv.type[c] =  0;
//                        mesh->BCv.val[c]  = -Vel*sin(dip) * exp(-pow( (mesh->zg_coord[l] - (ztop-Hslab)), 2) / (2.0*sigma*sigma) );
//                    }
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
                    
                    if ( (k==0 || k==NCX-1) && l==NCZ-1 ) {
                        mesh->BCp.type[c] =  0;
                        mesh->BCp.val[c]  =  0;
                    }
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
                        mesh->BCt.valW[l] = TW;
                    }
                    
                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typE[l] = 0;
                        mesh->BCt.valE[l] = TE;
                    }
                    
                    // SOUTH
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typS[k] = 1;
                        mesh->BCt.valS[k] = mesh->T[c];
                    }
                    
                    // NORTH
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typN[k] = 1;
                        mesh->BCt.valN[k] = TN;
                    }
                    
                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = TN;
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
