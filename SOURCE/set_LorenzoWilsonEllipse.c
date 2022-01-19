#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "header_MDOODZ.h"

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SET INITIAL TOPOGRPAHIC MARKER CHAIN
void BuildInitialTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling ) {
    
    int k;
    double TopoLevel = -0.0e3/scaling.L; // sets zero initial topography
    double h_pert = model.user3/scaling.L;   
    double xmin   = model.xmin;                                     // xmin
    double xmax   = model.xmax;                                     // xmax
    
    //for ( k=0; k<topo_chain->Nb_part; k++ ) {
    //    topo_chain->z[k]     = TopoLevel;
    //    topo_chain->phase[k] = 0;
    //}
    
    for ( k=0; k<topo_chain->Nb_part; k++ ) {
        topo_chain->z[k]     = TopoLevel + h_pert*(3330.0-2800.0)/2800.0*cos(2*M_PI*topo_chain->x[k]/(xmax-xmin)); //Hplateau = Hroot*(RHOmantle-RHOcrust)/RHOcrust
        topo_chain->phase[k] = 0;
    }
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// SET INITIAL MARKER FIELD
void SetParticles( markers *particles, scale scaling, params model, mat_prop *materials  ) {
    
    int np, iEllipse, s1, s2;                                               // Running variables
    double xmin            = model.xmin;                                    // xmin
    double xmax            = model.xmax;                                    // xmax
    double depthLAB        = model.user1/scaling.L;                         // Depth lithosphere
    double depthMoho       = model.user2/scaling.L;                         // Depth crust
    double thickLowerCrust = model.user5/scaling.L;                         // Thickness of lower crust
    double tempSurf        = (15.0+273.15)/scaling.T, Tpart;                // surface temperature
    double tempLAB         = (1350.0+273.15)/scaling.T;                     // adiabatic mantle temperature
    double tempBot         = 1886.0/scaling.T;                              // Temperature at the bottom of the domain - from phase diagram
    double adiabGradMantle = (tempBot-tempLAB)/(fabs(xmin - depthLAB));            // Adiabatic gradient through the upper mantle
    // Ellipse related
    int    ellipseCrust    = model.user3;                                   // Number of ellipses in the crust
    int    ellipseLith     = model.user4;                                   // Number of ellipses in the lithosphere
    double semiMinor       = 2.5e3/scaling.L;                               // Semi minor axis of ellipses
    double semiMajor       = 30.0e3/scaling.L;                              // Semi major axis of ellipses
    // Seed random numbers
    s1 = sizeof(int);
    s2 = sizeof(double);
    srand(197);
    // Loop on particles
    for( np=0; np<particles->Nb_part; np++ ) {
        
        // Standart initialisation of particles
//        particles->Vx[np]    = -particles->x[np]*model.EpsBG;         // set initial particle velocity (unused)
//        particles->Vz[np]    =  particles->z[np]*model.EpsBG;         // set initial particle velocity (unused)
        particles->phase[np] = 0;                                       // same phase number everywhere
        particles->d[np]     = materials->gs_ref[particles->phase[np]]; // Grain size
        particles->phi[np]   = 0.0;                                     // zero porosity everywhere
        
        //----------------- TEMPERATURE FIELD ----------------
        if (particles->z[np]>=-depthLAB) {
            Tpart = ((tempLAB-tempSurf)/depthLAB)*(-particles->z[np]) + tempSurf;
        } else {
            Tpart = adiabGradMantle*(-particles->z[np]) + tempLAB;
        }
        particles->T[np]  = Tpart;

        //----------------- LOWER CRUST ---------------------
        if (particles->z[np]< -(depthMoho-thickLowerCrust) && particles->z[np]>= -depthMoho) {
            particles->phase[np] = 12;
        }
        //----------------- LITHOSPHERIC MANTLE -------------
        if (particles->z[np] < -depthMoho) {
            particles->phase[np] = 3;
        }

        //----------------- UPPER MANTLE -----------------
        if (  (particles->T[np]>tempLAB) ) {
            particles->phase[np] = 4;
        }
        //--------------------------//
        // DENSITY --- Will be overriden
        particles->rho[np] = materials->rho[particles->phase[np]];

        //--------------------------//
        // SANITY CHECK
        if (particles->phase[np] > model.Nb_phases) {
            printf("Lazy bastard! Fix your particle phase ID! \n");
            exit(144);
        }
        //--------------------------//
    }
    
    // Build ellipses in the crust
    double d = 400.0e3, c = 350.0e3, f = 300.0e3, h = 250.0e3;
    double xCenter, xCenterOld, zCenter, X, Z;
    int phaseChoice[2] = {1,2};
    // ---------- Layer 1 ------------
    for (iEllipse=0; iEllipse<ellipseCrust; iEllipse++) {
        // Current ellipse properties
        double rDx = semiMajor*1.75 + (semiMajor)*(double)rand()/RAND_MAX;
        double rStart = (d-c)/scaling.L*(double)rand()/RAND_MAX-d/scaling.L;
        double rEnd   = (f-h)/scaling.L*(double)rand()/RAND_MAX-f/scaling.L;
        // Horizontal coordinate of ellipse center
        if (iEllipse==0) {
            xCenter = ((rEnd-rStart)*(double)rand()/RAND_MAX + rStart);
        } else{
            xCenter = xCenterOld + rDx;
        }
        xCenterOld = xCenter;
        // Vertical coordinate of ellipse center
        zCenter = (1.0*depthMoho/3.0-3.0*semiMinor)*(double)rand()/RAND_MAX - 1.0*depthMoho/3.0+semiMinor;
//        zCenter = (-2.0*semiMinor+1.0*depthMoho/3.0-semiMinor)*(double)rand()/RAND_MAX-1.0*depthMoho/3.0+semiMinor;
        // Assign random phase (either 1 = weak, or 2 = strong)
        int ellipsePhase = phaseChoice[rand()%2];
        // Check if we are inside the ellipse -> set the phase
        for (np=0; np<particles->Nb_part; np++) {
            X = particles->x[np];
            Z = particles->z[np];
            if ( ((X-xCenter)*(X-xCenter)/(semiMajor*semiMajor)+(Z-zCenter)*(Z-zCenter)/(semiMinor*semiMinor))<1.0 ) {
                particles->phase[np] = ellipsePhase;
            }
        }
    }
    // ---------- Layer 2 ------------
    for (iEllipse=0; iEllipse<ellipseCrust; iEllipse++) {
        // Current ellipse properties
        double rDx    = semiMajor*1.75 + (semiMajor)*(double)rand()/RAND_MAX;
        double rStart = (d-c)/scaling.L*(double)rand()/RAND_MAX-d/scaling.L;
        double rEnd   = (f-h)/scaling.L*(double)rand()/RAND_MAX-f/scaling.L;
        // Horizontal coordinate of ellipse center
        if (iEllipse==0) {
            xCenter = ((rEnd-rStart)*(double)rand()/RAND_MAX + rStart);
        } else{
            xCenter = xCenterOld + rDx;
        }
        xCenterOld = xCenter;
        // Vertical coordinate of ellipse center
        zCenter = (1.0*depthMoho/3.0-2.0*semiMinor)*(double)rand()/RAND_MAX - 2.0*depthMoho/3.0+semiMinor;
        // Assign random phase (either 1 = weak, or 2 = strong)
        int ellipsePhase = phaseChoice[rand()%2];
        // Check if we are inside the ellipse -> set the phase
        for (np=0; np<particles->Nb_part; np++) {
            X = particles->x[np];
            Z = particles->z[np];
            if ( ((X-xCenter)*(X-xCenter)/(semiMajor*semiMajor)+(Z-zCenter)*(Z-zCenter)/(semiMinor*semiMinor))<1.0 ) {
                particles->phase[np] = ellipsePhase;
            }
        }
    }
    
    // Build ellipses in the mantle
    double Z_Lith = 60.0e3/scaling.L;
    // ---------- Layer 1 ------------
    for (iEllipse=0; iEllipse<ellipseLith; iEllipse++) {
        // Current ellipse properties
        double rDx    = semiMajor*4.0 + (semiMajor)*(double)rand()/RAND_MAX;
        double rStart = (d-c)/scaling.L*(double)rand()/RAND_MAX-d/scaling.L;
        double rEnd   = (f-h)/scaling.L*(double)rand()/RAND_MAX-f/scaling.L;
        // Horizontal coordinate of ellipse center
        if (iEllipse==0) {
            xCenter = ((rEnd-rStart)*(double)rand()/RAND_MAX + rStart);
        } else{
            xCenter = xCenterOld + rDx;
        }
        xCenterOld = xCenter;
        // Vertical coordinate of ellipse center
        zCenter = 1.0/3.0*( Z_Lith-depthMoho-9.0*semiMinor )*(double)rand()/RAND_MAX - 1.0/3.0*( Z_Lith-depthMoho ) - depthMoho + semiMinor;
        // Check if we are inside the ellipse -> set the phase
        for (np=0; np<particles->Nb_part; np++) {
            X = particles->x[np];
            Z = particles->z[np];
            if ( ((X-xCenter)*(X-xCenter)/(semiMajor*semiMajor)+(Z-zCenter)*(Z-zCenter)/(semiMinor*semiMinor))<1.0 ) {
                particles->phase[np] = 5;
            }
        }
    }
    // ---------- Layer 2 ------------
    for (iEllipse=0; iEllipse<ellipseLith; iEllipse++) {
        // Current ellipse properties
        double rDx    = semiMajor*4.0 + (semiMajor)*(double)rand()/RAND_MAX;
        double rStart = (d-c)/scaling.L*(double)rand()/RAND_MAX-d/scaling.L - 5.0*semiMajor;
        double rEnd   = (f-h)/scaling.L*(double)rand()/RAND_MAX-f/scaling.L;
        // Horizontal coordinate of ellipse center
        if (iEllipse==0) {
            xCenter = ((rEnd-rStart)*(double)rand()/RAND_MAX + rStart);
        } else{
            xCenter = xCenterOld + rDx;
        }
        xCenterOld = xCenter;
        // Vertical coordinate of ellipse center
        zCenter = (1.0/3.0*( Z_Lith-depthMoho) )*(double)rand()/RAND_MAX - Z_Lith;
        // Check if we are inside the ellipse -> set the phase
        for (np=0; np<particles->Nb_part; np++) {
            X = particles->x[np];
            Z = particles->z[np];
            if ( ((X-xCenter)*(X-xCenter)/(semiMajor*semiMajor)+(Z-zCenter)*(Z-zCenter)/(semiMinor*semiMinor))<1.0 ) {
                particles->phase[np] = 5;
            }
        }
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
    // Numerical params
    int   kk, k, l, c, c1, np, it, nt = 120;
//    double *X, *Z, *XC, *ZC;
    int   NX, NZ, NCX, NCZ, NXVZ, NZVX;
    // Physical params
    double *velXBCWest, *velXBCEast;
    double tempNorth       = (273.15+15.0)/scaling.T, tempSouth = (1886.0)/scaling.T;     // Temperatures at boundaries (north and south) [K]
    double tempWest        = (1843.0)/scaling.T,      tempEast  = (1843.0)/scaling.T;     // Temperatures at boundaries (west and east) [K]
    double tempLAB         = (1350.0+273.15)/scaling.T;                                   // adiabatic mantle temperature
    double tempBot         = 1886.0/scaling.T;                                            // Temperature at the bottom of the domain - from phase diagram
    double secondsYear     = 3600.0*24.0*365.25;                                          // Seconds per year
    double timeExtension   = 50.0e6*secondsYear/(scaling.L/scaling.V);                    // Duration of extension [Myrs]
    double timeCooling     = 60.0e6*secondsYear/(scaling.L/scaling.V);                    // Duration of cooling   [Myrs]
    double timeRedConv     = 500.0e6*secondsYear/(scaling.L/scaling.V);                   // Point in time at which convergence rate is reduced [Myrs]
    double plateVelExt     = 1.0e-2/secondsYear/scaling.V;                                // Absolute plate velocity for extension [cm/yr]
    double plateVelCool    = 0.0e-2/secondsYear/scaling.V;                                // Absolute plate velocity for convergence [cm/yr]
    double plateVelConv    = -1.5e-2/secondsYear/scaling.V;                               // Absolute plate velocity for convergence [cm/yr]
    double plateVelRedConv = -0.3e-2/secondsYear/scaling.V;                               // Absolute plate velocity for reduced convergence [cm/yr]
    double depthLAB        = model->user1/scaling.L;                                      // Depth of the LAB [m]
    double depthVelChange  = model->zmin/2.0;                                             // Turnover point for inflow/outflow velocity [m]
    double areaBelowVChan  = fabs(model->zmin)-fabs(depthVelChange);                      // Area below turnover point [m2/m]
    double depthSerpent    = -9.0e3/scaling.L;                                            // Depth of serpentinisation [m]
    int    serpentinisationSwitch = 1;                                                    // 1 = serpentinisation, 0 = no serpentinisation
    
    /* ------------- SETTINGS -------------*/
    // Set unresaonable conductivity in the mantle to generate an adiabatic mantle as initial condition
    int    phase_ast = 4;                     // This has to be the ASTHENOSPHERE phase number
    double k_crazy = 13.6*materials->k[phase_ast];
    // Check for first time step only
    if (model->step == 0) {
        materials->k_eff[phase_ast] = k_crazy;
        printf("Running with crazy conductivity for the asthenosphere!!\n");
    }
    else {
        materials->k_eff[phase_ast] = materials->k[phase_ast];
        printf("Running with normal conductivity for the asthenosphere...\n");
    }

    // Set boundary velocity for deformation phase
    double plateVelocity;                                                               // Applied plate velocity
    if (model->time <= timeExtension){                                                  // ------------ EXTENSION
        plateVelocity = plateVelExt;
        printf("Def. Phase --> Extension ( %.2f cm/yr )\n",plateVelocity*scaling.V*100.0*secondsYear);
    }
    if(model->time > timeExtension && model->time <= timeExtension + timeCooling) {     // ------------ COOLING
        plateVelocity = plateVelCool;
        printf("Def. Phase --> Cooling ( %.2f cm/yr )\n",plateVelocity*scaling.V*100.0*secondsYear);
    }
    if (model->time > timeExtension + timeCooling) {                                    // ------------ CONVERGENCE
        plateVelocity = plateVelConv;
        printf("Def. Phase --> Convergence ( %.2f cm/yr )\n",plateVelocity*scaling.V*100.0*secondsYear);
    }
    if (model->time >= timeRedConv) {                                                   // ------------ REDUCED CONVERGENCE
        plateVelocity = plateVelRedConv;
        printf("Def. Phase --> Reduced convergence ( %.2f cm/yr )\n",plateVelocity*scaling.V*100.0*secondsYear);
    }

    // Parameterised serpentinisation front
    if (model->time >= timeExtension+timeCooling-1.0e6*secondsYear/(scaling.L/scaling.V) && serpentinisationSwitch == 1) {
        for (np=0; np<particles->Nb_part; np++) {
            if (particles->phase[np]==3 && particles->z[np]>=depthSerpent || particles->phase[np]==4 && particles->z[np]>=depthSerpent || particles->phase[np]==5 && particles->z[np]>=depthSerpent) {
                particles->phase[np] = 8;
            }
        }
        serpentinisationSwitch = 0;
    }
    
    // Distinguish Adria & Europe
    if (model->time > timeExtension + timeCooling && model->time <= timeExtension + timeCooling + model->dt) {
     for( np=0; np<particles->Nb_part; np++ ) {
       if (particles->phase[np] == 0 && particles->x[np]>0.0e3/scaling.L) {
          particles->phase[np] = 9;
        }
        if (particles->phase[np] == 1 && particles->x[np]>0.0e3/scaling.L) {
           particles->phase[np] = 10;
         }
         if (particles->phase[np] == 2 && particles->x[np]>0.0e3/scaling.L) {
             particles->phase[np] = 11;
         }
         if (particles->phase[np] == 12 && particles->x[np]>0.0e3/scaling.L) {
             particles->phase[np] = 13;
         }
       }
     }

    // Numerics (staggered grid related)
    NX  = mesh->Nx;
    NZ  = mesh->Nz;
    NCX = NX-1;
    NCZ = NZ-1;
    NXVZ = NX+1;
    NZVX = NZ+1;
    
    /* %%%%%%%%%%%%%%%% SET BOUNDARY INFLOW/OUTFLOW VELOCITY %%%%%%%%%%%%%%%% */
    // Allocate boundary velocity arrays
    velXBCWest = DoodzCalloc(NZVX,sizeof(double));
    velXBCEast = DoodzCalloc(NZVX,sizeof(double));
    // Calculate boundary velocity (smoothed step function )
    k = 0; // Computed for WEST
    double dz = mesh->zg_coord[2] - mesh->zg_coord[1];                                                              // Vertical grid spacing
    double Dv = 1.0;                                                                                                // Velocity conductivity (for smoothing)
    double dt = (dz*dz)/Dv/2.1;                                                                                     // CFL for diffusion smoother
    double velFluxNorth, velFluxSouth, velFluxUp_smooth=0.0, velFluxLow_smooth=0.0, velFluxUp=0.0, velFluxLow=0.0;  // Fluxes
    double balanceCheck = 0.0, balanceCheck_smooth=0.0;                                                             // Check variables for velocity balance
    double velXNorth, velXSouth;                                                                                    // Velocity
    // Action
    for (it=0; it<nt; it++) {
        if (it==0) {
            for (l=NZVX-1; l>=0; l--) {                                                 // Counting top-down
                c = k + l*(NX);
                if ( mesh->BCu.type[c] != 30 ) {                                        // Ignore free air include ghost node at bottom
                    if (mesh->zvx_coord[l]>depthVelChange) {
                        velXBCWest[l] = plateVelocity;
                        velFluxUp = velFluxUp + velXBCWest[l]*dz;                       // Calculate flux in lithosphere
                    }
                    else if (mesh->zvx_coord[l]<=depthVelChange) {
                        velXBCWest[l] = -velFluxUp/((round(areaBelowVChan/dz)+1)*dz);   // Set velocity in transition zone such that ...
                        velFluxLow = velFluxLow + velXBCWest[l]*dz;                     // ... Flux in transition zone ...
                    }
//                    printf("qmL[%d] = %.2e; qmA[%d] = %.2e, Z[%d] = %.2f; BCu.type[%d] = %d; vxBCW[%d] = %.2e\n",l,velFluxUp,l,velFluxLow,l,mesh->zvx_coord[l]*scaling.L/1e3,l,mesh->BCu.type[c],l,velXBCWest[l]*scaling.V*1e2*(3600*24*365.25));
                    
                }
            }
            balanceCheck = velFluxLow + velFluxUp;                                      // ... ~= flux in lithosphere (chkq ~ 0)
            printf("balanceCheck = %.2e\n",balanceCheck);
        }
        // Start smoothing
        for (l=0; l<NZVX; l++) {
            c = k + l*(NX);
            if (mesh->BCu.type[c] != 30) {
                if (l==0) { velXBCWest[l] = velXBCWest[l+1]; velXSouth = velXBCWest[l+1];} // Ghost node value
                else {velXSouth = velXBCWest[l-1];}
                if (l==NZ || mesh->BCu.type[c+NX] == 30) {
                    velXBCWest[l] = plateVelocity;
                    velXNorth     = plateVelocity;
                }
                else {
                    velXNorth = velXBCWest[l+1];
                }
                velFluxNorth  = -Dv*(velXNorth     - velXBCWest[l])/dz;
                velFluxSouth  = -Dv*(velXBCWest[l] - velXSouth    )/dz;
                velXBCWest[l] = velXBCWest[l]      - dt*(velFluxNorth-velFluxSouth)/dz;
                if (it==nt-1) {
                    if (mesh->zvx_coord[l]>depthVelChange) {
                        velFluxUp_smooth = velFluxUp_smooth + velXBCWest[l]*dz;
                    }
                    else if (mesh->zvx_coord[l]<=depthVelChange) {
                        velFluxLow_smooth = velFluxLow_smooth + velXBCWest[l]*dz;
                    }
                }
            }
        }
        if(it==nt-1){
            balanceCheck_smooth = velFluxUp_smooth + velFluxLow_smooth;
            printf("WEST: velFluxUp_smooth = %.2e; velFluxLow_smooth = %.2e; balanceCheck_smooth = %.2e\n",velFluxUp_smooth,velFluxLow_smooth,balanceCheck_smooth);
        }
    }
    // Calculate boundary velocity (smoothed step function )
    k = NX-1; // Computed for East
    velFluxUp_smooth=0.0, velFluxLow_smooth=0.0, velFluxUp=0.0, velFluxLow=0.0;         // Fluxes
    balanceCheck = 0.0, balanceCheck_smooth=0.0;                                        // Check variables for velocity balance
    // Action
    for (it=0; it<nt; it++) {
        if (it==0) {
            for (l=NZVX-1; l>=0; l--) {                                                 // Counting top-down
                c = k + l*(NX);
                if ( mesh->BCu.type[c] != 30 ) {                                        // Ignore free air include ghost node at bottom
                    if (mesh->zvx_coord[l]>depthVelChange) {
                        velXBCEast[l] = plateVelocity;
                        velFluxUp = velFluxUp + velXBCEast[l]*dz;                       // Calculate flux in lithosphere
                    }
                    else if (mesh->zvx_coord[l]<=depthVelChange) {
                        velXBCEast[l] = -velFluxUp/((round(areaBelowVChan/dz)+1)*dz);   // Set velocity in transition zone such that ...
                        velFluxLow = velFluxLow + velXBCEast[l]*dz;                     // ... Flux in transition zone ...
                    }
//                    printf("qmL[%d] = %.2e; qmA[%d] = %.2e, Z[%d] = %.2f; BCu.type[%d] = %d; vxBCW[%d] = %.2e\n",l,velFluxUp,l,velFluxLow,l,mesh->zvx_coord[l]*scaling.L/1e3,l,mesh->BCu.type[c],l,velXBCEast[l]*scaling.V*1e2*(3600*24*365.25));
                    
                }
            }
            balanceCheck = velFluxLow + velFluxUp;                                      // ... ~= flux in lithosphere (chkq ~ 0)
            printf("balanceCheck = %.2e\n",balanceCheck);
        }
        // Start smoothing
        for (l=0; l<NZVX; l++) {
            c = k + l*(NX);
            if (mesh->BCu.type[c] != 30) {
                if (l==0) { velXBCEast[l] = velXBCEast[l+1]; velXSouth = velXBCEast[l+1];} // Ghost node value
                else {velXSouth = velXBCEast[l-1];}
                if (l==NZ || mesh->BCu.type[c+NX] == 30) {
                    velXBCEast[l] = plateVelocity;
                    velXNorth     = plateVelocity;
                }
                else {
                    velXNorth = velXBCEast[l+1];
                }
                velFluxNorth  = -Dv*(velXNorth     - velXBCEast[l])/dz;
                velFluxSouth  = -Dv*(velXBCEast[l] - velXSouth    )/dz;
                velXBCEast[l] = velXBCEast[l]      - dt*(velFluxNorth-velFluxSouth)/dz;
                if (it==nt-1) {
                    if (mesh->zvx_coord[l]>depthVelChange) {
                        velFluxUp_smooth = velFluxUp_smooth + velXBCEast[l]*dz;
                    }
                    else if (mesh->zvx_coord[l]<=depthVelChange) {
                        velFluxLow_smooth = velFluxLow_smooth + velXBCEast[l]*dz;
                    }
                }
            }
        }
        if(it==nt-1){
            balanceCheck_smooth = velFluxUp_smooth + velFluxLow_smooth;
            printf("EAST: velFluxUp_smooth = %.2e; velFluxLow_smooth = %.2e; balanceCheck_smooth = %.2e\n",velFluxUp_smooth,velFluxLow_smooth,balanceCheck_smooth);
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
                        mesh->BCu.val[c]  = -velXBCWest[l]/2.0;
                    }
                    
                    // Matching BC nodes EAST
                    if (k==mesh->Nx-1 ) {
                        mesh->BCu.type[c] = 0;
                        mesh->BCu.val[c]  =  velXBCEast[l]/2.0;
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
                    }
                    
                    // Matching BC nodes NORTH
                    if (l==mesh->Nz-1 ) {
                        mesh->BCv.type[c] = 0;
                        mesh->BCv.val[c]  = 0.0;//mesh->zg_coord[l] * model->EpsBG;
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
		
		for (l=0; l<mesh->Nz-1; l++) {
			for (k=0; k<mesh->Nx-1; k++) {
				
				c = k + l*(NCX);
                
                if ( mesh->BCt.type[c] != 30 ) {
                    
                    // WEST
                    if ( k==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typW[l] = 0;
                        mesh->BCt.valW[l] = tempWest;
                    }
                    
                    // EAST
                    if ( k==NCX-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typE[l] = 0;
                        mesh->BCt.valE[l] = tempEast;
                    }
                    
                    // SOUTH
                    if ( l==0 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typS[k] = 1;
                        mesh->BCt.valS[k] = tempSouth;
                    }
                    
                    // NORTH
                    if ( l==NCZ-1 ) {
                        mesh->BCt.type[c] = 0;
                        mesh->BCt.typN[k] = 1;
                        mesh->BCt.valN[k] = tempNorth;
                    }
                    
                    // FREE SURFACE
                    else {
                        if ((mesh->BCt.type[c] == -1 || mesh->BCt.type[c] == 1 || mesh->BCt.type[c] == 0) && mesh->BCt.type[c+NCX] == 30) {
                            mesh->BCt.type[c] = 1;
                            mesh->BCt.val[c]  = tempNorth;
                        }
                    }
                    
                    
                }
                
			}
		}

    /* %%%%%%%%%%% FREE ALLOCATED MEMORY  %%%%%%%%%%% */
    free(velXBCWest);
    free(velXBCEast);
    
//	free(X);
//	free(Z);
//	free(XC);
//	free(ZC);
	printf("Velocity and pressure were initialised\n");
	printf("Boundary conditions were set up\n");
	
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
