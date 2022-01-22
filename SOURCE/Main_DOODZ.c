// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2022  MDOODZ Developper team
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
#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "cholmod.h"
#include "header_MDOODZ.h"


#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int main( int nargs, char *args[] ) {
    
    int          istep, irestart, writer = 0, writer_step;
    char         *fin_name, *PartFileName;
    params       model;
    grid         mesh;
    Nparams      Nmodel;
    markers      particles, topo_chain, topo_chain_ini;
    mat_prop     materials;
    clock_t      t_omp, t_omp_step;
    scale        scaling;
    SparseMat    Stokes, Jacob;
    DirectSolver CholmodSolver;
    surface      topo, topo_ini;
    int          nstag;
    SparseMat    StokesA, StokesB, StokesC, StokesD;
    SparseMat    JacobA,  JacobB,  JacobC,  JacobD;
    int          Nx, Nz, Ncx, Ncz;
    int          IsNewtonStep, *TrackNewtonSteps, IsFirstNewtonStep, IsJacobianUsed;
    int cent=1, vert=0, prop=1, interp=0, vxnodes=-1, vznodes=-2;
    FILE        *GNUplotPipe;
    
    double *rx_abs, *rz_abs, *rp_abs, *rx_rel, *rz_rel, *rp_rel;
    
    // Initialise integrated quantities
    mesh.W    = 0.0; // Work
    mesh.Ut   = 0.0; // heat
    mesh.Ue   = 0.0; // elastic energy
    
#ifdef _NEW_INPUT_
    // Input file name
    if ( nargs < 3 ) {
        printf( "NEW INPUT: You should enter the setup file and the initial particle file names as command line arguments.\nExiting...\n" );
        exit(1);
    }
    else {
        asprintf(&fin_name,"%s", args[1]);
        asprintf(&PartFileName,"%s", args[2]);
    }
#else
    // Input file name
    if ( nargs < 2 ) {
        printf( "OLD INPUT: You should (at least) enter the setup file name as a command line argument.\nExiting...\n" );
        exit(1);
    }
    else {
        asprintf(&fin_name,"%s", args[1]);
    }
#endif
    
    printf("\n********************************************************\n");
    printf("************ Starting MDOODZ 6.0 simulation ************\n");
    printf("********************************************************\n");
    
    // Read input data
    ReadInputFile( fin_name, &istep, &irestart, &writer, &writer_step, &model, &scaling, &materials, &particles, &Nmodel);
    model.L0 = model.xmax - model.xmin; // Save initial length
    
    printf("*************************************\n");
    printf("****** Allocate and initialise ******\n");
    printf("*************************************\n");
    
    // Multi resolution allocation pointers
    GridAlloc( &mesh, &model );
    
    // Initialise grid coordinates
    SetGridCoordinates( &mesh, &model, model.Nx, model.Nz );
    
    //    // Initial solution fields (Fine mesh)
    //    SetBCs( &mesh, &model, scaling , &particles, &materials );
    //    InitialiseSolutionFields( &mesh, &model );
    
    // Get grid indices
    GridIndices( &mesh );
    
    rx_abs  = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); rx_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    rz_abs  = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); rz_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    rp_abs  = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); rp_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    TrackNewtonSteps = DoodzCalloc(Nmodel.nit_max+1, sizeof(int));
    
    Nx = mesh.Nx; Nz = mesh.Nz; Ncx = Nx-1; Ncz = Nz-1;
    if ( model.aniso  == 1 ) model.Newton = 1;
    if ( model.Newton == 1 ) IsNewtonStep = 1;
    if ( model.Newton == 0 ) IsNewtonStep = 0;
    
    printf("*************************************\n");
    printf("******* Initialize particles ********\n");
    printf("*************************************\n");
    
    // Allocate particle fields
    PartAlloc( &particles, &model );
    
    // Allocate marker chain
    if ( model.free_surf == 1 ) AllocateMarkerChain( &topo,     &topo_chain,     model );
    if ( model.free_surf == 1 ) AllocateMarkerChain( &topo_ini, &topo_chain_ini, model );
    
    // Set new particle distribution
    if ( irestart == 0 ) {
        
        model.step = 0;
        model.time = 0.0;
        
        if ( model.no_markers == 0 ) {
            
            // Initialise particle fields
            PartInit( &particles, &model );
            
#ifdef _NEW_INPUT_
            // Initial grid tags
            model.BC_setup_type = 1; // eventually it should be set from the input file
            SetBCs_new( &mesh, &model, scaling, &particles, &materials);
            
            LoadIniParticles( PartFileName, &particles, &mesh, &topo_chain, &topo_chain_ini, &model, scaling );
            
            if ( model.free_surf == 1 ) {
                // Project topography on vertices
                ProjectTopography( &topo, &topo_chain, model, mesh, scaling, mesh.xg_coord, 0 );
                
                // Marker chain polynomial fit
                MarkerChainPolyFit( &topo, &topo_chain, model, mesh );
                
                // Call cell flagging routine for free surface calculations
                CellFlagging( &mesh, model, topo, scaling );
            }
            
#else
            // Initial grid tags
            SetBCs( &mesh, &model, scaling , &particles, &materials, &topo);
            if ( model.free_surf == 1 ) {
                
                // Define the horizontal position of the surface marker chain
                SetTopoChainHorizontalCoords( &topo,     &topo_chain,     model, mesh, scaling );
                SetTopoChainHorizontalCoords( &topo_ini, &topo_chain_ini, model, mesh, scaling );
                
                // Define the vertical position of the surface marker chain
                BuildInitialTopography( &topo,     &topo_chain,     model, mesh, scaling );
                BuildInitialTopography( &topo_ini, &topo_chain_ini, model, mesh, scaling );
                
                // Project topography on vertices
                ProjectTopography( &topo, &topo_chain, model, mesh, scaling, mesh.xg_coord, 0 );
                
                // Marker chain polynomial fit
                MarkerChainPolyFit( &topo, &topo_chain, model, mesh );
                
                // Call cell flagging routine for free surface calculations
                CellFlagging( &mesh, model, topo, scaling );
            }
            // Set particles coordinates
            PutPartInBox( &particles, &mesh, model, topo, scaling );
            
            // Set phases on particles
            SetParticles( &particles, scaling, model, &materials );
            
#endif
            if ( model.free_surf == 1 ) CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );

            // Create phase percentage arrays
            if (model.cpc==-1) CountPartCell_BEN( &particles, &mesh, model, topo, 0, scaling );
            if (model.cpc== 0) CountPartCell_Old( &particles, &mesh, model, topo, 0, scaling  );
            if (model.cpc== 1) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 0, scaling  );
            if (model.cpc== 2) CountPartCell2    ( &particles, &mesh, model, topo, topo_ini, 0, scaling  );
            
            P2Mastah( &model, particles, materials.eta0, &mesh, mesh.eta_s, mesh.BCg.type,  0, 0, interp, vert, model.itp_stencil);
            P2Mastah( &model, particles, materials.eta0, &mesh, mesh.eta_n, mesh.BCp.type,  0, 0, interp, cent, model.itp_stencil);
            
            P2Mastah( &model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
            P2Mastah( &model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            
            P2Mastah( &model, particles, particles.T,  &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
            P2Mastah( &model, particles, materials.Cv, &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);
            
            P2Mastah( &model, particles, materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
            P2Mastah( &model, particles, materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);
            
            // UpdateDensity( &mesh, &particles, &materials, &model, &scaling );

            // if ( model.eqn_state > 0) {
            //     P2Mastah( &model, particles, particles.rho, &mesh, mesh.rho_s, mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
            //     P2Mastah( &model, particles, particles.rho, &mesh, mesh.rho_n, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            // }
            // else {
            //     P2Mastah( &model, particles, materials.rho, &mesh, mesh.rho_s, mesh.BCg.type,  0, 0, interp, vert, model.itp_stencil);
            //     P2Mastah( &model, particles, materials.rho, &mesh, mesh.rho_n, mesh.BCp.type,  0, 0, interp, cent, model.itp_stencil);
            // }

            // MinMaxArrayTag( mesh.rho_s,      scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
            // MinMaxArrayTag( mesh.rho_n,      scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            // exit(1);
            if ( model.noisy == 1 ) {
                MinMaxArray( mesh.u_in,  scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
                MinMaxArray( mesh.v_in,  scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
                MinMaxArray( mesh.p_in,  scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
            }
            // Initial solution fields (Fine mesh)
#ifdef _NEW_INPUT_
            SetBCs_new( &mesh, &model, scaling , &particles, &materials );
#else
            SetBCs( &mesh, &model, scaling , &particles, &materials, &topo);
#endif
            InitialiseSolutionFields( &mesh, &model );
            
            MinMaxArray( mesh.u_in,  scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in,  scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in,  scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
            
            printf("*************************************\n");
            printf("****** Initialize temperature *******\n");
            printf("*************************************\n");
            
            // Get energy and related material parameters from particles
            P2Mastah( &model, particles, materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
            P2Mastah( &model, particles, materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);
            P2Mastah( &model, particles, particles.T,     &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
            P2Mastah( &model, particles, materials.Cv,    &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);
            P2Mastah( &model, particles, materials.Qr,    &mesh, mesh.Qr, mesh.BCp.type,  0, 0, interp, cent, 1);
            
#ifdef _NEW_INPUT_
            SetBCs_new( &mesh, &model, scaling , &particles, &materials );
#else
            SetBCs( &mesh, &model, scaling , &particles, &materials, &topo);
#endif
            if ( model.thermal_eq == 1 ) ThermalSteps( &mesh, model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, model.cooling_time, scaling );
            if ( model.therm_pert == 1 ) SetThermalPert( &mesh, model, scaling );
            //            Interp_Grid2P_centroids ( particles, particles.T,    &mesh, mesh.T, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &model );
            
            Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &model );
            ArrayEqualArray( mesh.T0_n, mesh.T, (mesh.Nx-1)*(mesh.Nz-1) );
            
            //--------------------------------------------------------------------------------------------------------
            
            printf("*************************************\n");
            printf("******** Initialize pressure ********\n");
            printf("*************************************\n");

            // Compute dev. strain rate tensor components    
            StrainRateComponents( &mesh, scaling, &model );
            
            for (int iter=0; iter<10; iter++) {
                
                // if ( model.eqn_state > 0 ) {
                    UpdateDensity( &mesh, &particles, &materials, &model, &scaling );
                // }
                // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
                // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
                
                // Free surface - subgrid density correction
                if ( model.free_surf == 1 ) {
                    SurfaceDensityCorrection( &mesh, model, topo, scaling  );
                }
                // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
                
                // Lithostatic pressure for initial visco-plastic viscosity field
                ComputeLithostaticPressure( &mesh, &model, materials.rho[0], scaling, 1 );
                Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_lith, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
                P2Mastah( &model, particles, particles.P,     &mesh, mesh.p_in , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                ArrayEqualArray( mesh.p_in, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
                ArrayEqualArray( mesh.p0_n, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
                
                MinMaxArrayTag( mesh.p_in,       scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P initial ", mesh.BCp.type );
            }
            
            InterpCentroidsToVerticesDouble( mesh.p0_n, mesh.p0_s, &mesh, &model );
            //            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &model );
            //            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &model );
            
            
            printf("*************************************\n");
            printf("******* Initialize grain size *******\n");
            printf("*************************************\n");
            
            // Grain size
            InitialiseGrainSizeParticles( &particles, &materials );
            P2Mastah( &model, particles, particles.d,     &mesh, mesh.d_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            ArrayEqualArray( mesh.d0_n, mesh.d_n,  (mesh.Nx-1)*(mesh.Nz-1) );
            
            printf("*************************************\n");
            printf("******** Initialize density *********\n");
            printf("*************************************\n");
            
            // if ( model.eqn_state > 0 ) {
            UpdateDensity( &mesh, &particles, &materials, &model, &scaling );
            // }
            // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
            // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
            
            printf("*************************************\n");
            printf("****** Initialize composition *******\n");
            printf("*************************************\n");
            
            //            if ( model.diffuse_X == 1 ) {
            //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
            //                Diffuse_X(&mesh, &model, &scaling);
            //                Interp_Grid2P( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
            //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
            //                Interp_P2N ( particles, particles.X, &mesh, mesh.Xreac_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &model );
            //            }
            
            P2Mastah( &model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            P2Mastah( &model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
            
            if ( model.aniso == 1 ) {
                InitialiseDirectorVector ( &mesh, &particles, &model, &materials );
                P2Mastah( &model, particles, particles.nx,     &mesh, mesh.nx0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                P2Mastah( &model, particles, particles.nz,     &mesh, mesh.nz0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                P2Mastah( &model, particles, particles.nx,     &mesh, mesh.nx0_s , mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
                P2Mastah( &model, particles, particles.nz,     &mesh, mesh.nz0_s , mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
                NormalizeDirector( &mesh, mesh.nx0_n, mesh.nz0_n, mesh.nx0_s, mesh.nz0_s, &model );
                FiniteStrainAspectRatio ( &mesh, scaling, model, &particles );
                P2Mastah( &model, particles, materials.aniso_factor,     &mesh, mesh.aniso_factor_n , mesh.BCp.type,  0, 0, interp, cent, model.itp_stencil);
                P2Mastah( &model, particles, materials.aniso_factor,     &mesh, mesh.aniso_factor_s , mesh.BCg.type,  0, 0, interp, vert, model.itp_stencil);
            }
            //
            printf("*************************************\n");
            printf("******* Initialize viscosity ********\n");
            printf("*************************************\n");
            
            if ( model.iselastic == 1 ) ShearModCompExpGrid( &mesh, materials, model, scaling );
            
            // Compute cohesion and friction angle on the grid
            CohesionFrictionDilationGrid( &mesh, &particles, materials, model, scaling );
            ShearModCompExpGrid( &mesh, materials, model, scaling );
            Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_in, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
            Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T,    mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &model );
            NonNewtonianViscosityGrid (     &mesh, &materials, &model, Nmodel, &scaling );
            // MinMaxArrayTag( mesh.rho_s,      scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
            // MinMaxArrayTag( mesh.rho_n,      scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
            // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );
            // exit(1);
                        // P2Mastah( &model, particles, particles.rho,   &mesh, mesh.rho0_n, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                        // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );

        } // end of no_markers --- debug
        else {
            
            Initialise1DArrayChar( mesh.BCu.type,  (mesh.Nx+0)*(mesh.Nz+1), -1 ); // make sure dofs are activated
            Initialise1DArrayChar( mesh.BCv.type,  (mesh.Nx+1)*(mesh.Nz+0), -1 );
            Initialise1DArrayChar( mesh.BCp.type,  (mesh.Nx-1)*(mesh.Nz-1), -1 );
            SetUpModel_NoMarkers ( &mesh, &model, &scaling );
            InitialiseSolutionFields( &mesh, &model );
            StrainRateComponents( &mesh, scaling, &model );
            if ( model.iselastic == 1 ) ShearModCompExpGrid( &mesh, materials, model, scaling );
            SetBCs( &mesh, &model, scaling , &particles, &materials, &topo);
            ComputeLithostaticPressure( &mesh, &model, materials.rho[0], scaling, 1 );
            NonNewtonianViscosityGrid (     &mesh, &materials, &model, Nmodel, &scaling );
            ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

        }
        
        printf("Number of phases : %d\n", model.Nb_phases);
        if ( model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n,       scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d0        ", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n,        scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_lith,     scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P litho   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p0_n,       scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_in,       scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T,          scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s,       scaling.S,   (mesh.Nx-0)*(mesh.Nz-0), "mu_s      ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n,       scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n      ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_s,      scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_phys_s, scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_phys_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_phys_n, scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho_s,      scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n,      scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            for (int p=0; p<model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }
        }
        
        printf("*************************************\n");
        printf("******** Initialize timestep ********\n");
        printf("*************************************\n");
        
        DefineInitialTimestep( &model, &mesh, particles, materials, scaling );
        
        if ( model.rec_T_P_x_z == 1 ) {
            ArrayEqualArray( particles.T0,   particles.T, particles.Nb_part );
            ArrayEqualArray( particles.P0,   particles.P, particles.Nb_part );
            ArrayEqualArray( particles.x0,   particles.x, particles.Nb_part );
            ArrayEqualArray( particles.z0,   particles.z, particles.Nb_part );
            ArrayEqualArray( particles.Tmax, particles.T, particles.Nb_part );
            ArrayEqualArray( particles.Pmax, particles.P, particles.Nb_part );
        }
        
        printf("*************************************\n");
        printf("*** Write initial file or restart ***\n");
        printf("*************************************\n");
        
        // Write initial output
        if ( writer == 1 ) {
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output",  materials, scaling );
            if ( model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles",  materials, scaling );
        }
        
        // Set initial stresses and pressure to zero
        //        Initialise1DArrayDouble( particles.X,      particles.Nb_part, 0.0 ); // set X to zero for the first time step
        //        Initialise1DArrayDouble( particles.P,      particles.Nb_part, 0.0 ); // now dynamic pressure...
        //        Initialise1DArrayDouble( mesh.p_in,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxxd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.szzd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxz,   (mesh.Nx)  *(mesh.Nz)  , 0.0 );
        // Generate deformation maps
        if ( model.def_maps == 1 ) GenerateDeformationMaps( &mesh, &materials, &model, Nmodel, &scaling );
        particles.Nb_part_ini = particles.Nb_part;
    }
    else {
        // Which step do we restart from (BreakpointXXXX.dat)
        model.step = istep;
        printf("Restarting from step number %05d...\n", model.step);
        LoadBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, &model, &topo, &topo_ini, scaling  );
        SetGridCoordinates( &mesh, &model, model.Nx, model.Nz ); // Overwrite previous grid
    }
    
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                             TIME LOOP : en avant les Doud'series !
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    
    model.step += 1;
    
    if ( model.GNUplot_residuals == 1 ) GNUplotPipe = popen ("gnuplot -persistent", "w");
    
    for (; model.step<=model.Nt; model.step++) {
        
        printf("*****************************************************\n");
        printf("****************** Time step %05d ******************\n", model.step);
        printf("*****************************************************\n");
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        t_omp_step = (double)omp_get_wtime();
        
        // Old time step
        model.dt0 = model.dt;
        
        // Define new time step
        EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &model, scaling, &mesh, 0 );
        printf("Selected dt = %2.2e\n", model.dt*scaling.t);
        
        // Save initial dt
        model.dt0 = model.dt;
        
        // Track particule generation index
        Initialise1DArrayInt( particles.generation,   particles.Nb_part  , 0 );
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        if (model.no_markers == 0 ) {
            
            // Remove particles that would be above the surface
            if ( model.free_surf == 1 ) {
                CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );
                CellFlagging( &mesh, model, topo, scaling );
                ArrayEqualArray(     topo.height0,     topo.height, mesh.Nx );
                ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );
            }
            
            // Interpolate material properties from particles to nodes
            t_omp = (double)omp_get_wtime();
            
            // Energy - interpolate thermal parameters and advected energy
            
            // Get energy and related material parameters from particles
            P2Mastah( &model, particles, materials.Cv,     &mesh, mesh.Cv,     mesh.BCp.type,  0, 0, interp, cent, 1);
            P2Mastah( &model, particles, materials.Qr,     &mesh, mesh.Qr,     mesh.BCp.type,  0, 0, interp, cent, 1);
            
            P2Mastah ( &model, particles, materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
            P2Mastah ( &model, particles, materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);
            
            // Get T and dTdt from previous step from particles
            P2Mastah( &model, particles, particles.T,     &mesh, mesh.T0_n,     mesh.BCp.type,  1, 0, interp, cent, 1);
            P2Mastah( &model, particles, particles.divth, &mesh, mesh.divth0_n, mesh.BCp.type,  1, 0, interp, cent, 1);
            
            // Make sure T is up to date for rheology evaluation
            ArrayEqualArray( mesh.T, mesh.T0_n, (mesh.Nx-1)*(mesh.Nz-1) );
            
            //-----------------------------------------------------------------------------------------------------------
            // Interp P --> p0_n , p0_s
            P2Mastah( &model, particles, particles.P,     &mesh, mesh.p0_n,   mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            // P2Mastah( &model, particles, particles.rho,   &mesh, mesh.rho0_n, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);

            // Get physical properties that are constant throughout each timestep
            // if ( model.eqn_state  > 0 ) {
                UpdateDensity( &mesh, &particles, &materials, &model, &scaling );
            // }
            // else {
            //     P2Mastah( &model, particles, materials.rho, &mesh, mesh.rho_s, mesh.BCg.type,  0, 0, interp, vert, model.itp_stencil);
            //     P2Mastah( &model, particles, materials.rho, &mesh, mesh.rho_n, mesh.BCp.type,  0, 0, interp, cent, model.itp_stencil);
            // }
            
            
            
            // Free surface - subgrid density correction
            if ( model.free_surf == 1 ) {
                SurfaceDensityCorrection( &mesh, model, topo, scaling  );
            }
            
            // Lithostatic pressure
            ArrayEqualArray(  mesh.p_lith0,   mesh.p_lith, Ncx*Ncz );
            ComputeLithostaticPressure( &mesh, &model, materials.rho[0], scaling, 1 );
            
            // Elasticity - interpolate advected/rotated stresses
            if  ( model.iselastic == 1 ) {
                
                // Get old stresses from particles
                if (model.StressUpdate==1)                     OldDeviatoricStressesPressure( &mesh, &particles, scaling, &model );
                
                if (model.StressUpdate==0){
                    P2Mastah( &model, particles, particles.sxxd,    &mesh, mesh.sxxd0, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                    P2Mastah( &model, particles, particles.szzd,    &mesh, mesh.szzd0, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                    P2Mastah( &model, particles, particles.sxz,     &mesh, mesh.sxz0,  mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
                }
                
                //                ArrayEqualArray(  mesh.sxxd0,  mesh.sxxd, Ncx*Ncz );
                //                ArrayEqualArray(  mesh.szzd0,  mesh.szzd, Ncx*Ncz );
                //                ArrayEqualArray(  mesh.sxz0,   mesh.sxz,   Nx*Nz );
                
                InterpCentroidsToVerticesDouble( mesh.sxxd0, mesh.sxxd0_s, &mesh, &model );
                InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &model );
                InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &model );
                
                // Interpolate shear modulus
                ShearModCompExpGrid( &mesh, materials, model, scaling );
            }
            
            // Director vector
            if (model.aniso == 1 ) {
                P2Mastah( &model, particles, particles.nx,     &mesh, mesh.nx0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                P2Mastah( &model, particles, particles.nz,     &mesh, mesh.nz0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
                P2Mastah( &model, particles, particles.nx,     &mesh, mesh.nx0_s , mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
                P2Mastah( &model, particles, particles.nz,     &mesh, mesh.nz0_s , mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
                NormalizeDirector( &mesh, mesh.nx0_n, mesh.nz0_n, mesh.nx0_s, mesh.nz0_s, &model );
                FiniteStrainAspectRatio ( &mesh, scaling, model, &particles );
                P2Mastah( &model, particles, materials.aniso_factor,     &mesh, mesh.aniso_factor_n , mesh.BCp.type,  0, 0, interp, cent, model.itp_stencil);
                P2Mastah( &model, particles, materials.aniso_factor,     &mesh, mesh.aniso_factor_s , mesh.BCg.type,  0, 0, interp, vert, model.itp_stencil);
            }
            
            P2Mastah( &model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            P2Mastah( &model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
            
            P2Mastah( &model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, model.itp_stencil);
            P2Mastah( &model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            
            // Diffuse rheological contrasts
            //            if (model.diffuse_X == 1) {
            //                Interp_P2C ( particles, particles.X, &mesh, mesh.X0_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
            //                Diffuse_X(&mesh, &model, &scaling);
            //                Interp_Grid2P( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
            //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
            //                Interp_P2N ( particles, particles.X, &mesh, mesh.Xreac_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &model );
            //            }
            
            // Interpolate Grain size
            P2Mastah( &model, particles, particles.d,     &mesh, mesh.d0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            ArrayEqualArray(  mesh.d_n,  mesh.d0_n, Ncx*Ncz );
            
            // Interpolate Melt fraction
            P2Mastah( &model, particles, particles.phi,   &mesh, mesh.phi0_n , mesh.BCp.type,  1, 0, interp, cent, model.itp_stencil);
            
            //-------------------------------------------------------------------------------------------------------------
            
            // Compute cohesion and friction angle on the grid
            CohesionFrictionDilationGrid( &mesh, &particles, materials, model, scaling );
            
            // Detect compressible cells
            if ( model.compressible == 1 ) DetectCompressibleCells ( &mesh, &model );
            
            // Compute cohesion and friction angle on the grid
            CohesionFrictionDilationGrid( &mesh, &particles, materials, model, scaling );
            
            // Detect compressible cells
            if (model.compressible == 1) DetectCompressibleCells ( &mesh, &model );
            
        }
        else {
            ArrayEqualArray(  mesh.p0_n,   mesh.p_in, Ncx*Ncz );
            ArrayEqualArray(  mesh.sxxd0,  mesh.sxxd, Ncx*Ncz );
            ArrayEqualArray(  mesh.szzd0,  mesh.szzd, Ncx*Ncz );
            ArrayEqualArray(  mesh.sxz0,   mesh.sxz,   Nx*Nz );
            
            InterpCentroidsToVerticesDouble( mesh.sxxd0, mesh.sxxd0_s, &mesh, &model );
            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &model );
            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &model );
            
            UpdateDensity( &mesh, &particles, &materials, &model, &scaling );
            ShearModCompExpGrid( &mesh, materials, model, scaling );
            CohesionFrictionDilationGrid( &mesh, &particles, materials, model, scaling );
            // Detect compressible cells
            if (model.compressible == 1) DetectCompressibleCells ( &mesh, &model );
        }
        
        
        // Min/Max interpolated fields
        if ( model.noisy == 1 ) {
            MinMaxArray(particles.sxxd, scaling.S, particles.Nb_part, "sxxd part  ");
            MinMaxArray(particles.T,   scaling.T,   particles.Nb_part, "T part    ");
            MinMaxArrayTag( mesh.p0_n,         scaling.S,    (mesh.Nx-1)*(mesh.Nz-1),   "p0_n",   mesh.BCp.type );
            MinMaxArrayTag( mesh.sxz0,     scaling.S,   (mesh.Nx)*(mesh.Nz),     "sxz0    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.sxxd0,    scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "sxx0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.szzd0,    scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "szz0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s,     scaling.S,   (mesh.Nx)*(mesh.Nz),     "mu_s    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n,     scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.C_s,      scaling.S,   (mesh.Nx)*(mesh.Nz),     "C_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.C_n,      scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "C_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.fric_s,   180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "fric_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.fric_n,   180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "fric_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.dil_s,    180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "dil_s   ", mesh.BCg.type );
            MinMaxArrayTag( mesh.dil_n,    180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "dil_n   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.strain_s,   1.0,       (mesh.Nx)*(mesh.Nz),     "strain_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.strain_n,   1.0,       (mesh.Nx-1)*(mesh.Nz-1), "strain_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.bet_s,    1.0/scaling.S,   (mesh.Nx)*(mesh.Nz),     "beta_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.bet_n,    1.0/scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "beta_n  ", mesh.BCp.type );
            IsInfArray2DFP(mesh.bet_n, (mesh.Nx-1)*(mesh.Nz-1));
            IsNanArray2DFP(mesh.bet_n, (mesh.Nx-1)*(mesh.Nz-1));
            MinMaxArrayTag( mesh.T0_n,     scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T       ", mesh.BCt.type );
            MinMaxArrayTag( mesh.p_in,     scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P       ", mesh.BCt.type );
            MinMaxArrayI  ( mesh.comp_cells, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "comp_cells" );
            MinMaxArrayTag( mesh.rho_s,      scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n,      scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho0_n,     scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );
           
            
            for (int p=0; p<model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }
            
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.nx0_n,    1.0,   (mesh.Nx-1)*(mesh.Nz-1), "nx0_n  ", mesh.BCp.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.nz0_n,    1.0,   (mesh.Nx-1)*(mesh.Nz-1), "nz0_n  ", mesh.BCp.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.nx0_s,    1.0,   (mesh.Nx)*(mesh.Nz),     "nx0_s  ", mesh.BCg.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.nz0_s,    1.0,   (mesh.Nx)*(mesh.Nz),     "n0z_s  ", mesh.BCg.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.FS_AR_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "FS_AR_n", mesh.BCp.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.FS_AR_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "FS_AR_s", mesh.BCg.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.aniso_factor_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "aniso_factor_n", mesh.BCp.type );
            if  ( model.aniso == 1 ) MinMaxArrayTag( mesh.aniso_factor_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "aniso_factor_s", mesh.BCg.type );
        }
        
        printf("** Time for particles interpolations I = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );
        
        
        
        
//        struct timespec begin, end;
//        clock_gettime(CLOCK_REALTIME, &begin);
//
//        ViscosityDerivatives( &mesh, &materials, &model, Nmodel, &scaling );
//
//        // Stop measuring time and calculate the elapsed time
//        clock_gettime(CLOCK_REALTIME, &end);
//        long seconds = end.tv_sec - begin.tv_sec;
//        long nanoseconds = end.tv_nsec - begin.tv_nsec;
//        double elapsed = seconds + nanoseconds*1e-9;
//
//        printf("** Time for ViscosityDerivatives = %3f sec\n",  elapsed );
//
//        exit(0);
    
        
        if ( model.ismechanical == 1 ) {
            
            // Allocate and initialise solution and RHS vectors
#ifdef _NEW_INPUT_
            SetBCs_new( &mesh, &model, scaling , &particles, &materials );
#else
            SetBCs( &mesh, &model, scaling , &particles, &materials, &topo);
#endif
            // Reset fields and BC values if needed
            //        if ( model.ispureshear_ale == 1 ) InitialiseSolutionFields( &mesh, &model );
            InitialiseSolutionFields( &mesh, &model );
            EvalNumberOfEquations( &mesh, &Stokes );
            if ( IsNewtonStep == 1 ) EvalNumberOfEquations( &mesh, &Jacob  );
            SAlloc( &Stokes,  Stokes.neq );
            if ( IsNewtonStep == 1 ) SAlloc(  &Jacob,   Jacob.neq );
            if ( model.decoupled_solve == 1 ) {
                SAlloc( &StokesA, Stokes.neq_mom );
                SAlloc( &StokesB, Stokes.neq_mom );
                SAlloc( &StokesC, Stokes.neq_cont);
                SAlloc( &StokesD, Stokes.neq_cont );
            }
            if ( IsNewtonStep == 1 ) {
                SAlloc( &JacobA, Stokes.neq_mom );
                SAlloc( &JacobB, Stokes.neq_mom );
                SAlloc( &JacobC, Stokes.neq_cont);
                SAlloc( &JacobD, Stokes.neq_cont );
            }
            printf( "Linear systems allocated\n");
            printf( "neq_tot = %d, neq_mom = %d, neq_cont = %d\n", Stokes.neq, Stokes.neq_mom, Stokes.neq_cont );
            
            // Set vector x = [u;p]
            InitialiseSolutionVector( &mesh, &Stokes, &model );
            
            //------------------------------------------------------------------------------------------------------------------------------//
            
            // Non-linear iteration cycle
            Nmodel.nit       = 0;
            model.nit        = 0;
            Nmodel.stagnated = 0;
            nstag            = 0;
            IsJacobianUsed    = 0;
            
            ArrayEqualArray( mesh.p_start,    mesh.p_in,      (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.u_start,    mesh.u_in,      (mesh.Nx)  *(mesh.Nz+1) );
            ArrayEqualArray( mesh.v_start,    mesh.v_in,      (mesh.Nx+1)*(mesh.Nz)   );
            //            ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );
            
            // Set up solver context
            if ( model.decoupled_solve == 1 ) {
                cholmod_start( &CholmodSolver.c );
                if ( Nmodel.nit==0 ) CholmodSolver.Analyze = 1;
                printf("Run CHOLMOD analysis yes/no: %d \n", CholmodSolver.Analyze);
            }
            
            int Nmax_picard = Nmodel.nit_max;
            
            // If Picard 2 Newton is activated: force initial Newton step
            if ( IsNewtonStep == 1 && Nmodel.Picard2Newton == 1 ) {
                model.Newton = 0;
                IsFirstNewtonStep  = 1;
            }
            
            while ( Nmodel.nit <= Nmax_picard && nstag<model.nstagmax) {
                
                if ( Nmodel.nit > 0 && Nmodel.Picard2Newton == 1 ) {
                    if ( rx_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || rz_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || rp_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.nit>= Nmodel.nit_Pic_max) {
                        if ( IsFirstNewtonStep == 0 ) CholmodSolver.Analyze = 0;
                        if ( IsFirstNewtonStep == 1 ) {
                            cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                            CholmodSolver.Analyze = 1;
                            IsFirstNewtonStep    = 0;
                        }
                        model.Newton = 1;
                    }
                }

                // Determine whether Jacobian matrix should be assembled
                if ( model.Newton == 1 || model.aniso == 1 ) {
                    IsJacobianUsed = 1;
                }
                
                printf("**********************************************\n");
                if ( model.Newton == 0 ) { printf("*** Picard it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, model.step); TrackNewtonSteps[Nmodel.nit] = 0;}
                if ( model.Newton == 1 ) { printf("*** Newton it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, model.step); TrackNewtonSteps[Nmodel.nit] = 1;}
                printf("**********************************************\n");
                
                // Update non-linear rheology
                UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, materials, &model, &Nmodel, scaling, 0, 0.0 );
                RheologicalOperators( &mesh, &model, &scaling, 0 );                               // ??????????? déjà fait dans UpdateNonLinearity
                NonNewtonianViscosityGrid (     &mesh, &materials, &model, Nmodel, &scaling );    // ??????????? déjà fait dans UpdateNonLinearity

                if ( model.noisy == 1 ) {
                    MinMaxArrayTag( mesh.p0_n,      scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.p_in,      scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.div_u,     scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "div       ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exxd,      scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "exxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.ezzd,      scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "ezzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exz,       scaling.E, (mesh.Nx-0)*(mesh.Nz-0), "exz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.sxxd,      scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "sxxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.szzd,      scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "szzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.sxz,       scaling.S, (mesh.Nx-0)*(mesh.Nz-0), "sxz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.eta_phys_s, scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_phys_s", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_phys_n, scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho_s,      scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.rho_n,      scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho0_n,     scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.d_n,        scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
                }
                
                if ( model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output_BeforeSolve", materials, scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles_BeforeSolve", materials, scaling );
                }
                
                // Build discrete system of equations - Linearised Picard
                if ( model.decoupled_solve == 0 ) BuildStokesOperator           ( &mesh, model, 0, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, 1 );
                if ( model.decoupled_solve == 1 ) BuildStokesOperatorDecoupled  ( &mesh, model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1 );
                
                // Build discrete system of equations - Jacobian
                ViscosityDerivatives( &mesh, &materials, &model, Nmodel, &scaling );
                if ( model.Newton == 1 && Nmodel.nit > 0 ) RheologicalOperators( &mesh, &model, &scaling, 1 );
                if ( IsJacobianUsed == 1 )                  BuildJacobianOperatorDecoupled( &mesh, model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in,  &Jacob,  &JacobA,  &JacobB,  &JacobC,   &JacobD, 1 );
                
                //                MinMaxArrayTag( mesh.detadexx_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadexx_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadezz_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadezz_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadexx_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadgxz_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadexx_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadexx_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadezz_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadezz_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadgxz_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadgxz_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D11_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D21_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D12_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D22_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D13_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D23_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D14_n,      scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D24_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D31_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "D31_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D32_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "D32_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D33_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "D33_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D34_s,      scaling.eta, (mesh.Nx)*(mesh.Nz),     "D34_s     ", mesh.BCg.type );
                
                // Diagonal scaling
                if ( model.diag_scaling ) {
                    if ( model.Newton          == 0 )  ExtractDiagonalScale( &StokesA, &StokesB, &StokesC, &StokesD );
                    if ( model.Newton          == 1 )  ExtractDiagonalScale( &JacobA,  &JacobB,  &JacobC,   &JacobD );
                    if ( model.Newton          == 1 )  ArrayEqualArray(StokesA.d, JacobA.d, StokesA.neq);
                    if ( model.Newton          == 1 )  ArrayEqualArray(StokesC.d, JacobC.d, StokesC.neq);
                    ScaleMatrix( &StokesA, &StokesB, &StokesC, &StokesD );
                }
                
                // Needs to be done after matrix assembly since diagonal scaling is used in there
                printf("---- Non-linear residual ----\n");
                RheologicalOperators( &mesh, &model, &scaling, 0 );
                if ( model.decoupled_solve == 0 ) EvaluateStokesResidual( &Stokes, &Nmodel, &mesh, model, scaling, 0 );
                if ( model.decoupled_solve == 1 ) EvaluateStokesResidualDecoupled( &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &Nmodel, &mesh, model, scaling, 0 );
                printf("---- Non-linear residual ----\n");
                
                if ( model.Newton == 1 && model.aniso == 0 ) {
                    ArrayEqualArray( JacobA.F, StokesA.F,  StokesA.neq );
                    ArrayEqualArray( JacobC.F, StokesC.F,  StokesC.neq );
                }
                if ( model.Newton == 1 && model.diag_scaling ) ScaleMatrix( &JacobA,  &JacobB,  &JacobC,  &JacobD  );
                
                // Store residuals
                Nmodel.resx_f = Nmodel.resx; rx_abs[Nmodel.nit] = Nmodel.resx; rx_rel[Nmodel.nit] = Nmodel.resx/Nmodel.resx0;
                Nmodel.resz_f = Nmodel.resz; rz_abs[Nmodel.nit] = Nmodel.resz; rz_rel[Nmodel.nit] = Nmodel.resz/Nmodel.resz0;
                Nmodel.resp_f = Nmodel.resp; rp_abs[Nmodel.nit] = Nmodel.resp; rp_rel[Nmodel.nit] = Nmodel.resp/Nmodel.resp0;
                
                if ( model.write_debug == 1 ) WriteResiduals( mesh, model, Nmodel, scaling );
                
                // if pass --> clear matrix break
                if ( (Nmodel.resx < Nmodel.abs_tol_u || Nmodel.resx/Nmodel.resx0 < Nmodel.rel_tol_u) && (Nmodel.resz < Nmodel.abs_tol_u || Nmodel.resz/Nmodel.resz0 < Nmodel.rel_tol_u) && (Nmodel.resp < Nmodel.abs_tol_p || Nmodel.resp/Nmodel.resp0 < Nmodel.rel_tol_p) ) {
                    printf( "Non-linear solver converged to abs_tol_u = %2.2e abs_tol_p = %2.2e rel_tol_u = %2.2e rel_tol_p = %2.2e\n", Nmodel.abs_tol_u, Nmodel.abs_tol_p, Nmodel.rel_tol_u, Nmodel.rel_tol_p );
                    FreeSparseSystems( IsJacobianUsed, model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                    break;
                }
                
                // Direct solve
                t_omp = (double)omp_get_wtime();
                if ( model.decoupled_solve == 0 ) SolveStokesDefect( &Stokes, &CholmodSolver, &Nmodel, &mesh, &model, &particles, &topo_chain, &topo, materials, scaling );
                if ( model.decoupled_solve == 1 ) {
                    
                    if ( IsJacobianUsed==0 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &model, &particles, &topo_chain, &topo, materials, scaling, &StokesA, &StokesB, &StokesC );
                    if ( IsJacobianUsed==1 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &model, &particles, &topo_chain, &topo, materials, scaling,  &JacobA,  &JacobB,  &JacobC );
                    
                    if ( Nmodel.stagnated == 1 && model.safe_mode <= 0 ) {
                        printf( "Non-linear solver stagnated to res_u = %2.2e res_z = %2.2e\n", Nmodel.resx_f, Nmodel.resz_f );
                        printf( "You may want to try setting line_search_min > 0.0\n Good luck good man!\n");
                        FreeSparseSystems( IsJacobianUsed, model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                        break;
                        
                    }
                    
                }

                if ( Nmodel.stagnated == 1 && model.iselastic == 1 && model.safe_mode == 1 ) {
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (abs_tol_u = %2.2e abs_tol_p = %2.2e)\e[m\n", Nmodel.abs_tol_u, Nmodel.abs_tol_p );
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (rel_tol_u = %2.2e rel_tol_p = %2.2e)\e[m\n", Nmodel.rel_tol_u, Nmodel.rel_tol_p );
                    printf( "\e[1;31mReducing the timestep, and restart the iterations cycle...\e[m\n");
                    printf( "Before reduction: model.dt =, %2.2e\n", model.dt*scaling.t);
                    // ----------------------
                    model.dt /= model.safe_dt_div;
                    printf( "Timestep divided by %2.2f => NEW CURRENT model.dt =, %2.2e\n", model.safe_dt_div, model.dt*scaling.t);
                    Nmodel.stagnated = 0;
                    // ----------------------
                    Nmodel.nit = -1;
                    printf( "Restart solutions\n");
                    ArrayEqualArray( mesh.p_in, mesh.p_start,  (mesh.Nx-1)*(mesh.Nz-1) );
                    ArrayEqualArray( mesh.u_in, mesh.u_start,  (mesh.Nx)  *(mesh.Nz+1) );
                    ArrayEqualArray( mesh.v_in, mesh.v_start,  (mesh.Nx+1)*(mesh.Nz)   );
                    nstag++;
                    //-----------------------
                    printf( "nstag value = %02d - nstagmax = %02d\n", nstag, model.nstagmax);
                    if (nstag==model.nstagmax) {
                        printf( "CheckDoudzOut!!\n");
                        exit(0);
                    //-----------------------
                    }
                }
                FreeSparseSystems( IsJacobianUsed, model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );

                Nmodel.nit++; model.nit = Nmodel.nit;
            }
            
            // Clean solver context
            if ( model.decoupled_solve == 1 && Nmodel.nit>0 ) {
                printf("Cleaning up Cholesky factors --- nit = %02d\n",  Nmodel.nit);
                cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                cholmod_finish( &CholmodSolver.c );
            }
            
            // Update rheology
            UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, materials, &model, &Nmodel, scaling, 0, 1.0 );
            ArrayEqualArray( mesh.p_in, mesh.p_corr, (mesh.Nx-1)*(mesh.Nz-1) ); // make sure corrected pressure is used
            
            // Min/Max velocities
            MinMaxArray( mesh.u_in,  scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in,  scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in,  scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
            MinMaxArray( mesh.div_u, scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "  div(V)" );
            
            printf("--------------------------------------------------------------\n");
            int i, nit;
            if (Nmodel.nit>=Nmodel.nit_max)  nit = Nmodel.nit_max;
            if (Nmodel.nit< Nmodel.nit_max)  nit = Nmodel.nit;
            if (Nmodel.Picard2Newton == 1 )  printf("Picard 2 Newton is activated with condition: %2.2e\n", Nmodel.Pic2NewtCond);
            for (i=0; i<=nit; i++) {
                if (TrackNewtonSteps[i] == 0) printf("Pic. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, rx_abs[i], rz_abs[i], rp_abs[i], rx_rel[i], rz_rel[i], rp_rel[i]);
                if (TrackNewtonSteps[i] == 1) printf("New. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, rx_abs[i], rz_abs[i], rp_abs[i], rx_rel[i], rz_rel[i], rp_rel[i]);
                if (i == Nmodel.nit_max && model.safe_mode == 1) {
                    printf("Exit: Max iteration reached: Nmodel.nit_max = %02d! Check what you wanna do now...\n",Nmodel.nit_max);
                    if ( (Nmodel.resx < Nmodel.abs_tol_u) && (Nmodel.resz < Nmodel.abs_tol_u) && (Nmodel.resp < Nmodel.abs_tol_p) ) {}
                    else {exit(1);}
                }
            }
            printf("--------------------------------------------------------------\n");
            
            // plot residuals
            if ( model.GNUplot_residuals == 1 ) { // && model.step % writer_step == 0 
                
                printf("DOING GNU PLOTTING!\n");

                int NumCommands = 3;
                char *GNUplotCommands[] = {"set title \"Non-linear residuals\"", "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5", "set pointintervalbox 3"};
                for (i=0; i<NumCommands; i++) fprintf(GNUplotPipe, "%s \n", GNUplotCommands[i]); //Send commands to gnuplot one by one.
                
                fprintf(GNUplotPipe, "plot '-' with linespoints linestyle 1\n");
                for (i=0; i< Nmodel.nit+1; i++) {
                    fprintf(GNUplotPipe, "%lf %lf \n", (double)i, log10(rx_rel[i])); //Write the data to a temporary file
                }
                fprintf(GNUplotPipe, "e\n");
                fflush(GNUplotPipe);
            }
        
        }
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE --- Free surface - interpolate velocity components on the free surface
        if ( model.free_surf == 1 ) {
            SurfaceVelocity( &mesh, model, &topo, &topo_chain, scaling );
            MinMaxArray( topo_chain.Vx,      scaling.V, topo_chain.Nb_part,       "Vx surf." );
            MinMaxArray( topo_chain.Vz,      scaling.V, topo_chain.Nb_part,       "Vz surf." );
            MinMaxArray( topo_chain_ini.Vx,  scaling.V, topo_chain_ini.Nb_part,   "Vx surf. ini." );
            MinMaxArray( topo_chain_ini.Vz,  scaling.V, topo_chain_ini.Nb_part,   "Vz surf. ini." );
        }
        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE
        
        if (model.StressUpdate==1) TotalStresses( &mesh, &particles, scaling, &model );
        
        // Update stresses on markers
        if (model.StressUpdate==0) UpdateParticleStress(  &mesh, &particles, &model, &materials, &scaling );
        
        // Update pressure on markers
        UpdateParticlePressure( &mesh, scaling, model, &particles, &materials );
        //        Interp_Grid2P_centroids( particles, particles.P, &mesh, mesh.p_in, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
        UpdateParticleX( &mesh, scaling, model, &particles, &materials );
        
        // Grain size evolution
        UpdateParticleGrainSize( &mesh, scaling, model, &particles, &materials );
        
        if ( model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n   , scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d0", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n    , scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d ", mesh.BCp.type );
            MinMaxArrayPart( particles.d, scaling.L, particles.Nb_part, "d on markers", particles.phase ) ;
        }
                
        // Update phi on the particles
        UpdateParticlePhi( &mesh, scaling, model, &particles, &materials );
        
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        if (model.isthermal == 1 ) {
            
            printf("*************************************\n");
            printf("*********** Thermal solver **********\n");
            printf("*************************************\n");
            
            t_omp = (double)omp_get_wtime();
            
            // Matrix assembly and direct solve
            EnergyDirectSolve( &mesh, model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, model.dt, model.shear_heat, model.adiab_heat, scaling, 1 );
            MinMaxArray(particles.T, scaling.T, particles.Nb_part, "T part. before UpdateParticleEnergy");
            
            // Update energy on particles
            UpdateParticleEnergy( &mesh, scaling, model, &particles, &materials );
            MinMaxArray(particles.T, scaling.T, particles.Nb_part, "T part. after UpdateParticleEnergy");
            
            // Calculate energies
            if ( model.write_energies == 1 ) Energies( &mesh, model, scaling );
            
            printf("** Time for Thermal solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        }
        
        //--------------------------------------------------------------------------------------------------------------------------------//
        
        if (model.ProgReac == 1 ) {
            
            printf("*************************************\n");
            printf("********** Chemical solver **********\n");
            printf("*************************************\n");
            
            P2Mastah ( &model, particles, materials.k_chem, &mesh, mesh.kc_x, mesh.BCu.type,  0, 0, interp, vxnodes, model.itp_stencil);
            P2Mastah ( &model, particles, materials.k_chem, &mesh, mesh.kc_z, mesh.BCv.type,  0, 0, interp, vznodes, model.itp_stencil);
            ChemicalDirectSolve( &mesh, model, &particles, &materials, model.dt, scaling );
            
        }
        //--------------------------------------------------------------------------------------------------------------------------------//
        
        // Update maximum pressure and temperature on markers
        if ( model.rec_T_P_x_z == 1 )  UpdateMaxPT( scaling, model, &particles );
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        
        if ( model.advection == 1 ) {
            
            printf("*************************************\n");
            printf("************** Advection ************\n");
            printf("*************************************\n");
            
            t_omp = (double)omp_get_wtime();
            double dt_solve = model.dt;
            int    nsub=1.0, isub;
            double dt_sub;
            EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &model, scaling, &mesh, 0 );
            
            if ( model.dt < 0.95*dt_solve ) { // if dt advection is lower than dt used for thermo-mechanical solve then split dt
                nsub   = ceil(dt_solve/model.dt);
                dt_sub = dt_solve / nsub;
                printf("dt advection = %2.9e --- dt_solve = %2.9e\n", model.dt*scaling.t, dt_solve*scaling.t );
                printf("dt is %lf larger than dt_solve: need %d substeps of %2.2e s\n", dt_solve/model.dt, nsub , dt_sub*scaling.t );
                printf("So: nsub*dt_sub = %2.2e for dt_solve = %2.2e\n", nsub*dt_sub*scaling.t, dt_solve*scaling.t);
                model.dt = dt_sub;
            }
            else {
                nsub     = 1; // if dt advection is larger than dt used for thermo-mechanical solve then, forget it, and set dt to dt_solve
                model.dt = dt_solve;
            }
            
            // Loop on substeps
            for (isub=0;isub<nsub;isub++) {
                
                printf("************** Advection step %03d of %03d: dtsub = %2.2e **************\n", isub, nsub, model.dt*scaling.t );
                
                // Advect domain boundaries
                if ( model.ispureshear_ale > 0 ) {
                    PureShearALE( &model, &mesh, &topo_chain, scaling );
                }
                
                if ( model.free_surf == 1 ) {
                    // Save old topo
                    ArrayEqualArray( topo_chain.z0, topo_chain.z, topo_chain.Nb_part ); // save old z
                    ArrayEqualArray( topo_chain_ini.z0, topo_chain_ini.z, topo_chain.Nb_part ); // save old z
                    
                    
                    // Advect free surface with RK4
                    RogerGuntherII( &topo_chain,     model, mesh, 1, scaling );
                    RogerGuntherII( &topo_chain_ini, model, mesh, 1, scaling );
                }
                
                //                // Advect free surface
                //                if ( model.free_surf == 1 ) {
                //                    AdvectFreeSurf( &topo_chain,     model, scaling );
                //                    AdvectFreeSurf( &topo_chain_ini, model, scaling );
                //                    MinMaxArray( topo_chain.z,      scaling.L, topo_chain.Nb_part,       "z surf.     " );
                //                    MinMaxArray( topo_chain_ini.z,  scaling.L, topo_chain_ini.Nb_part,   "z surf. ini." );
                //                }
                
                // Correction for particle inflow 0
                if (model.ispureshear_ale == -1 && model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, model, topo, 0 );
                
                // Advect fluid particles
                RogerGuntherII( &particles, model, mesh, 1, scaling );
                
                // Correction for particle inflow 1
                if (model.ispureshear_ale == -1 && model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, model, topo, 1 );
                
                // Update accumulated strain
                AccumulatedStrainII( &mesh, scaling, model, &particles,  mesh.xc_coord,  mesh.zc_coord, mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
                
                // Update deformation gradient tensor components
                if ( model.fstrain == 1 ) DeformationGradient( mesh, scaling, model, &particles );
                
                if ( model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output_BeforeSurfRemesh", materials, scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles_BeforeSurfRemesh", materials, scaling );
                }
                
                if ( model.free_surf == 1 ) {
                    
                    //                    for (int k=0;k<topo_chain.Nb_part;k++) {
                    //
                    //                        if (fabs(topo_chain.x[k]) < 0.5*model.surf_Winc ) {
                    //                            topo_chain.z[k] -= model.dt*model.surf_Vinc;
                    //                        }
                    //                    }
                    //
                    
                    if ( model.surf_remesh == 1 ) {
                        // Get current topography
                        ProjectTopography( &topo,     &topo_chain,     model, mesh, scaling, mesh.xg_coord, 0 );
                        ProjectTopography( &topo_ini, &topo_chain_ini, model, mesh, scaling, mesh.xg_coord, 0 );
                        MarkerChainPolyFit( &topo,     &topo_chain,     model, mesh );
                        MarkerChainPolyFit( &topo_ini, &topo_chain_ini, model, mesh );
                        
                        // Remesh free surface I
                        RemeshMarkerChain( &topo_chain,     &topo,     model, scaling, &mesh, 1 );
                        RemeshMarkerChain( &topo_chain_ini, &topo_ini, model, scaling, &mesh, 1 );
                    }
                    
                    
                    // Project topography on vertices
                    ProjectTopography( &topo,     &topo_chain,     model, mesh, scaling, mesh.xg_coord, 0 );
                    ProjectTopography( &topo_ini, &topo_chain_ini, model, mesh, scaling, mesh.xg_coord, 0 );
                    //                    ArrayEqualArray( topo.height0, topo.height, mesh.Nx );
                    //                    ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );
                    
                    if ( model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output_AfterSurfRemesh", materials, scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles_AfterSurfRemesh", materials, scaling );
                    }
                    
                    // Diffuse topography
                    if ( model.surf_processes >= 1 )  DiffuseAlongTopography( &mesh, model, scaling, topo.height, topo.height, mesh.Nx, 0.0, model.dt );
                    
                    // Marker chain polynomial fit
                    MarkerChainPolyFit( &topo,     &topo_chain,     model, mesh );
                    CorrectTopoIni( &particles, materials, &topo_chain_ini, &topo, model, scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, model, mesh );
                    
                    // Sedimentation
                    if ( model.surf_processes >= 2 ) {
                        AddPartSed( &particles, materials, &topo_chain, &topo, model, scaling, &mesh);
                        if (model.cpc==-1) CountPartCell_BEN( &particles, &mesh, model, topo, 0, scaling );
                        if (model.cpc== 0) CountPartCell_Old( &particles, &mesh, model, topo, 0, scaling );
                        if (model.cpc== 1) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 0, scaling );
                    }
                    
                    if ( model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Outputx", materials, scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particlesx", materials, scaling );
                    }
                    
                    // Remesh free surface II
                    RemeshMarkerChain( &topo_chain,     &topo,     model, scaling, &mesh, 2 );
                    RemeshMarkerChain( &topo_chain_ini, &topo_ini, model, scaling, &mesh, 2 );
                    CorrectTopoIni( &particles, materials, &topo_chain_ini, &topo, model, scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, model, mesh );
                    
                    // Remove particles that are above the surface
                    CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );
                    
                    // Call cell flagging routine for free surface calculations
                    CellFlagging( &mesh, model, topo, scaling );
                }
                
                printf("** Time for advection solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );
                
#ifdef _HDF5_
                if ( model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Outputxx", materials, scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particlesxx", materials, scaling );
                }
#endif
                
                //                printf("*************************************\n");
                //                printf("************** Reseeding ************\n");
                //                printf("*************************************\n");
                
                // Count the number of particle per cell
                t_omp = (double)omp_get_wtime();
                if (model.cpc ==-1)                       CountPartCell_BEN( &particles, &mesh, model, topo, 0, scaling );
                if (model.cpc == 0)                       CountPartCell_Old( &particles, &mesh, model, topo, 0, scaling );
                if (model.cpc == 1 && model.Reseed == 1 ) CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 1, scaling );
                if (model.cpc == 1)                       CountPartCell    ( &particles, &mesh, model, topo, topo_ini, 0, scaling );
                
                if (model.cpc == 2 && model.Reseed == 1 ) CountPartCell2   ( &particles, &mesh, model, topo, topo_ini, 1, scaling );
                if (model.cpc == 2)                       CountPartCell2   ( &particles, &mesh, model, topo, topo_ini, 0, scaling );
                
                printf("After re-seeding :\n");
                printf("Initial number of particles = %d\n", particles.Nb_part_ini);
                printf("New number of particles     = %d\n", particles.Nb_part    );
                
                printf("** Time for CountPartCell = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );
                
                // Remove particles that would be above the surface
                if ( model.free_surf == 1 ) {
                    CleanUpSurfaceParticles( &particles, &mesh, topo, scaling );
                    CellFlagging( &mesh, model, topo, scaling );
                }
                
                //                for (int iphase=0; iphase<model.Nb_phases; iphase++) {
                //                CheckSym( mesh.phase_perc_n[iphase], 1.0, mesh.Nx-1, mesh.Nz-1, "perc_n" );
                //                CheckSym( mesh.phase_perc_s[iphase], 1.0, mesh.Nx-0, mesh.Nz-0, "perc_s" );
                //                }
            }
            model.dt = dt_solve;
        }
        
        // Update time
        model.time += model.dt;
        
        // Free solution arrays
        SFree( &Stokes );
        if ( IsNewtonStep == 1 ) SFree( &Jacob );
        if ( model.decoupled_solve == 1 ) {
            SFree( &StokesA );
            SFree( &StokesB );
            SFree( &StokesC );
            SFree( &StokesD );
        }
        if ( IsNewtonStep == 1 ) {
            SFree( &JacobA );
            SFree( &JacobB );
            SFree( &JacobC );
            SFree( &JacobD );
        }
        DoodzFree( Stokes.eqn_u );
        DoodzFree( Stokes.eqn_v );
        DoodzFree( Stokes.eqn_p );
        if ( IsNewtonStep == 1 ) {
            DoodzFree( Jacob.eqn_u );
            DoodzFree( Jacob.eqn_v );
            DoodzFree( Jacob.eqn_p );
        }
        
        //------------------------------------------------------------------------------------------------------------------------------//
        
        // Write output data
        if ( writer == 1 && model.step % writer_step == 0 ) {
            
            printf("*************************************\n");
            printf("********* Write output files ********\n");
            printf("*************************************\n");
            
            // Breakpoint file
            t_omp = (double)omp_get_wtime();
            if ( model.delete_breakpoints == 1 ) DeletePreviousBreakpoint( model.step, writer_step  );
            MakeBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, model, &topo, &topo_ini, scaling );
            UpdateInputFile( fin_name, model.step);
            printf("** Time for Breakpoint file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
            
            // Visualisation file
#ifndef _VG_
            t_omp = (double)omp_get_wtime();
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, model, "Output", materials, scaling );
            if ( model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, model, "Particles", materials, scaling );
            printf("** Time for Output file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
#endif
        }
        
        printf("** Total timestep calculation time = %lf sec\n", (double)((double)omp_get_wtime() - t_omp_step) );
        printf("** Model time = %2.2e sec\n", model.time*scaling.t );
        printf("** Current dt = %2.2e sec, Old dt = %2.2e sec\n", model.dt*scaling.t, model.dt0*scaling.t );
        
        // if (Np0 - particles.Nb_part != 0) exit (1);
        
        
    }
    
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                           END TIME LOOP : c'est fini les Doud'series...
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    
    if ( model.GNUplot_residuals == 1) pclose(GNUplotPipe);
    
    //    // Free markers chains
    if ( model.free_surf == 1 )    FreeMarkerChain( &topo,     &topo_chain     );
    if ( model.free_surf == 1 )    FreeMarkerChain( &topo_ini, &topo_chain_ini );
    
    // Free Phase diagrams if activated
    if ( model.isPD == 1 ) FreePhaseDiagrams( &model );
    //
    // Free particle arrays
    PartFree( &particles, &model );
    
    // Free arrays
    GridFree( &mesh, &model );
    
    // Free char*'s
    free(model.input_file);
    free(fin_name);
#ifdef _NEW_INPUT_
    free(PartFileName);
#endif
    
    // GNU plot
    DoodzFree( rx_abs );
    DoodzFree( rz_abs );
    DoodzFree( rp_abs );
    DoodzFree( rx_rel );
    DoodzFree( rz_rel );
    DoodzFree( rp_rel );
    DoodzFree( TrackNewtonSteps );
    
    // just for the quartz-coesite case
    for ( int k=0; k<model.Nb_phases; k++) {
        
        if ( materials.density_model[k] == 4 ) {
            printf("Unloading 1D diagram for coesite quartz only baby...\n");
            DoodzFree(model.PD1Drho[k]);
        }
    }
    DoodzFree(model.PD1DnP);
    DoodzFree(model.PD1Drho);
    DoodzFree(model.PD1Dmin);
    DoodzFree(model.PD1Dmax);
    
    printf("\n********************************************************\n");
    printf("************* Ending MDOODZ 6.0 simulation *************\n");
    printf("********************************************************\n");
    
    // Exit success signal
    return(EXIT_SUCCESS);
}
