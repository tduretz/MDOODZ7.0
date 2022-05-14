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
#include "math.h"
#include "time.h"
#include "cholmod.h"
#include "mdoodz-private.h"
#include "mdoodz.h"

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

void RunMDOODZ(MdoodzInstance *instance) {
    int          istep, irestart, writer = 0, writer_step;
    char         *fin_name, *PartFileName;
    Nparams      Nmodel;
    clock_t      t_omp, t_omp_step;
    SparseMat    Stokes, Jacob;
    DirectSolver CholmodSolver;
    surface      topo, topo_ini;
    int          nstag;
    SparseMat    StokesA, StokesB, StokesC, StokesD;
    SparseMat    JacobA,  JacobB,  JacobC,  JacobD;
    int          Nx, Nz, Ncx, Ncz;
    int          IsNewtonStep, IsFirstNewtonStep, IsJacobianUsed;
    int cent=1, vert=0, prop=1, interp=0, vxnodes=-1, vznodes=-2;
    FILE        *GNUplotPipe;
    
    grid mesh;
    markers particles, topo_chain, topo_chain_ini;

    // Initialise integrated quantities
    mesh.Work       = 0.0; // Work
    mesh.Uthermal   = 0.0; // heat
    mesh.Uelastic   = 0.0; // elastic energy

    fin_name = instance->inputFileName;
    asprintf(&fin_name,"%s", instance->inputFileName);

    printf("\n********************************************************\n");
    printf("************ Starting MDOODZ 7.0 simulation ************\n");
    printf("********************************************************\n");

    // Read input data
    ReadInputFile( fin_name, &istep, &irestart, &writer, &writer_step, &instance->model, &instance->scaling, &instance->materials, &particles, &Nmodel);
    ValidateSetup(instance);
    instance->model.L0 = instance->model.xmax - instance->model.xmin; // Save initial length

    printf("*************************************\n");
    printf("****** Allocate and initialise ******\n");
    printf("*************************************\n");

    // Multi resolution allocation pointers
    GridAlloc( &mesh, &instance->model );

    // Initialise grid coordinates
    SetGridCoordinates( &mesh, &instance->model, instance->model.Nx, instance->model.Nz );

    // Get grid indices
    GridIndices( &mesh );

    // Initialise data for logs of iterations
    Nmodel.nit             = 0;
    Nmodel.rx_abs          = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); Nmodel.rx_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    Nmodel.rz_abs          = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); Nmodel.rz_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    Nmodel.rp_abs          = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); Nmodel.rp_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    Nmodel.LogIsNewtonStep = DoodzCalloc(Nmodel.nit_max+1, sizeof(int));

    Nx = mesh.Nx; Nz = mesh.Nz; Ncx = Nx-1; Ncz = Nz-1;
    if (instance->model.aniso  == 1 ) instance->model.Newton = 1;
    if (instance->model.Newton == 1 ) IsNewtonStep = 1;
    if (instance->model.Newton == 0 ) IsNewtonStep = 0;

    printf("*************************************\n");
    printf("******* Initialize particles ********\n");
    printf("*************************************\n");

    // Allocate particle fields
    PartAlloc( &particles, &instance->model );

    // Allocate marker chain
    if (instance->model.free_surf == 1 ) AllocateMarkerChain( &topo,     &topo_chain, instance->model );
    if (instance->model.free_surf == 1 ) AllocateMarkerChain( &topo_ini, &topo_chain_ini, instance->model );

    // Set new particle distribution
    if ( irestart == 0 ) {

      instance->model.step = 0;
      instance->model.time = 0.0;

        // Initialise particle fields
        PartInit( &particles, &instance->model );

        if (instance->crazyConductivity) {
          for (int n = 0; n < instance->crazyConductivity->nPhases; n++) {
            const int    phase           = instance->crazyConductivity->phases[n];
            const double k_crazy         = instance->crazyConductivity->multiplier * instance->materials.k[phase];
            instance->materials.k_eff[phase] = k_crazy;
          }
          printf("Running with crazy conductivity for the asthenosphere!!\n");
        }

        // Initial grid tags
        SetBCs(instance, &mesh);
        if (instance->model.free_surf == 1 ) {

            // Define the horizontal position of the surface marker chain
            SetTopoChainHorizontalCoords( &topo,     &topo_chain, instance->model, mesh, instance->scaling );
            SetTopoChainHorizontalCoords( &topo_ini, &topo_chain_ini, instance->model, mesh, instance->scaling );

            // Define the vertical position of the surface marker chain
            BuildInitialTopography(instance, &topo_chain );
            BuildInitialTopography(instance, &topo_chain_ini );

            // Project topography on vertices
            ProjectTopography( &topo, &topo_chain, instance->model, mesh, instance->scaling, mesh.xg_coord, 0 );

            // Marker chain polynomial fit
            MarkerChainPolyFit( &topo, &topo_chain, instance->model, mesh );

            // Call cell flagging routine for free surface calculations
            CellFlagging( &mesh, instance->model, topo, instance->scaling );
        }
        // Set particles coordinates
        PutPartInBox( &particles, &mesh, instance->model, topo, instance->scaling );

        // Set phases on particles
        SetParticles(instance, &particles);
        //PlotPhases(instance, &particles);

        MinMaxArray(particles.Vx, instance->scaling.V, particles.Nb_part, "Vxp init" );
        MinMaxArray(particles.Vz, instance->scaling.V, particles.Nb_part, "Vzp init" );
        MinMaxArray(particles.T, instance->scaling.T, particles.Nb_part,  "Tp init" );

        if (instance->model.free_surf == 1 ) CleanUpSurfaceParticles( &particles, &mesh, topo, instance->scaling );

        // Create phase percentage arrays
        if (instance->model.cpc==-1) CountPartCell_BEN( &particles, &mesh, instance->model, topo, 0, instance->scaling );
        if (instance->model.cpc== 0) CountPartCell_Old( &particles, &mesh, instance->model, topo, 0, instance->scaling  );
        if (instance->model.cpc== 1) CountPartCell    ( &particles, &mesh, instance->model, topo, topo_ini, 0, instance->scaling  );
        if (instance->model.cpc== 2) CountPartCell2    ( &particles, &mesh, instance->model, topo, topo_ini, 0, instance->scaling  );

        P2Mastah( &instance->model, particles, instance->materials.eta0, &mesh, mesh.eta_s, mesh.BCg.type,  0, 0, interp, vert, instance->model.itp_stencil);
        P2Mastah( &instance->model, particles, instance->materials.eta0, &mesh, mesh.eta_n, mesh.BCp.type,  0, 0, interp, cent, instance->model.itp_stencil);

        P2Mastah( &instance->model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, instance->model.itp_stencil);
        P2Mastah( &instance->model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);

        P2Mastah( &instance->model, particles, particles.T,  &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &instance->model, particles, instance->materials.Cv, &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);

        P2Mastah( &instance->model, particles, instance->materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah( &instance->model, particles, instance->materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);

        // UpdateDensity( &mesh, &particles, &instance->materials, &instance->model, &instance->scaling );

        // if ( instance->model.eqn_state > 0) {
        //     P2Mastah( &instance->model, particles, particles.rho, &mesh, mesh.rho_s, mesh.BCg.type,  1, 0, interp, vert, instance->model.itp_stencil);
        //     P2Mastah( &instance->model, particles, particles.rho, &mesh, mesh.rho_n, mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
        // }
        // else {
        //     P2Mastah( &instance->model, particles, instance->materials.rho, &mesh, mesh.rho_s, mesh.BCg.type,  0, 0, interp, vert, instance->model.itp_stencil);
        //     P2Mastah( &instance->model, particles, instance->materials.rho, &mesh, mesh.rho_n, mesh.BCp.type,  0, 0, interp, cent, instance->model.itp_stencil);
        // }

        // MinMaxArrayTag( mesh.rho_s,      instance->scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
        // MinMaxArrayTag( mesh.rho_n,      instance->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
        // exit(1);
        if (instance->model.noisy == 1 ) {
            MinMaxArray( mesh.u_in, instance->scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in, instance->scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
        }
        if (instance->crazyConductivity) {
          for (int n = 0; n < instance->crazyConductivity->nPhases; n++) {
            const int    phase           = instance->crazyConductivity->phases[n];
            const double k_crazy         = instance->crazyConductivity->multiplier * instance->materials.k[phase];
            instance->materials.k_eff[phase] = k_crazy;
          }
          printf("Running with crazy conductivity for the asthenosphere!!\n");
        }
        // Initial solution fields (Fine mesh)
        SetBCs(instance, &mesh);
        InitialiseSolutionFields( &mesh, &instance->model );

        MinMaxArray( mesh.u_in, instance->scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
        MinMaxArray( mesh.v_in, instance->scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
        MinMaxArray( mesh.p_in, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );

        printf("*************************************\n");
        printf("****** Initialize temperature *******\n");
        printf("*************************************\n");

        // Get energy and related material parameters from particles
        P2Mastah( &instance->model, particles, instance->materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah( &instance->model, particles, instance->materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);
        P2Mastah( &instance->model, particles, particles.T,     &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &instance->model, particles, instance->materials.Cv,    &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);
        P2Mastah( &instance->model, particles, instance->materials.Qr,    &mesh, mesh.Qr, mesh.BCp.type,  0, 0, interp, cent, 1);

        SetBCs(instance, &mesh);
        if (instance->model.thermal_eq == 1 ) ThermalSteps( &mesh, instance->model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, instance->model.cooling_time, instance->scaling );
        if (instance->model.therm_pert == 1 ) SetThermalPert( &mesh, instance->model, instance->scaling );
        //            Interp_Grid2P_centroids ( particles, particles.T,    &mesh, mesh.T, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &instance->model );

        Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &instance->model );
        ArrayEqualArray( mesh.T0_n, mesh.T, (mesh.Nx-1)*(mesh.Nz-1) );

        //--------------------------------------------------------------------------------------------------------

        printf("*************************************\n");
        printf("******** Initialize pressure ********\n");
        printf("*************************************\n");

        // Compute dev. strain rate tensor components
        StrainRateComponents( &mesh, instance->scaling, &instance->model );

        for (int iter=0; iter<10; iter++) {

            // if ( instance->model.eqn_state > 0 ) {
                UpdateDensity( &mesh, &particles, &instance->materials, &instance->model, &instance->scaling );
            // }
            // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &instance->model );
            // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

            // Free surface - subgrid density correction
            if (instance->model.free_surf == 1 ) {
                SurfaceDensityCorrection( &mesh, instance->model, topo, instance->scaling  );
            }
            // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

            // Lithostatic pressure for initial visco-plastic viscosity field
            ComputeLithostaticPressure( &mesh, &instance->model, instance->materials.rho[0], instance->scaling, 1 );
            Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_lith, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &instance->model );
            P2Mastah( &instance->model, particles, particles.P,     &mesh, mesh.p_in , mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
            ArrayEqualArray( mesh.p_in, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.p0_n, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );

            MinMaxArrayTag( mesh.p_in, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P initial ", mesh.BCp.type );
        }

        InterpCentroidsToVerticesDouble( mesh.p0_n, mesh.p0_s, &mesh, &instance->model );
        //            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &instance->model );
        //            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &instance->model );


        printf("*************************************\n");
        printf("******* Initialize grain size *******\n");
        printf("*************************************\n");

        // Grain size
        InitialiseGrainSizeParticles( &particles, &instance->materials );
        P2Mastah( &instance->model, particles, particles.d,     &mesh, mesh.d_n , mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
        ArrayEqualArray( mesh.d0_n, mesh.d_n,  (mesh.Nx-1)*(mesh.Nz-1) );
        MinMaxArrayTag( mesh.d_n, instance->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );

        printf("*************************************\n");
        printf("******** Initialize density *********\n");
        printf("*************************************\n");

        // if ( instance->model.eqn_state > 0 ) {
        UpdateDensity( &mesh, &particles, &instance->materials, &instance->model, &instance->scaling );
        // }
        // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &instance->model );
        // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

        printf("*************************************\n");
        printf("****** Initialize composition *******\n");
        printf("*************************************\n");

        //            if ( instance->model.diffuse_X == 1 ) {
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Diffuse_X(&mesh, &instance->model, &instance->scaling);
        //                Interp_Grid2P( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Interp_P2N ( particles, particles.X, &mesh, mesh.Xreac_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &instance->model );
        //            }

        P2Mastah( &instance->model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
        P2Mastah( &instance->model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, instance->model.itp_stencil);

        if (instance->model.aniso == 1 ) {
            InitialiseDirectorVector ( &mesh, &particles, &instance->model, &instance->materials );
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d1_n , mesh.BCp.type, -1, 0, interp, cent, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d2_n , mesh.BCp.type, -2, 0, interp, cent, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d1_s , mesh.BCg.type, -1, 0, interp, vert, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d2_s , mesh.BCg.type, -2, 0, interp, vert, instance->model.itp_stencil);
            FiniteStrainAspectRatio ( &mesh, instance->scaling, instance->model, &particles );
            P2Mastah( &instance->model, particles, instance->materials.aniso_factor,     &mesh, mesh.aniso_factor_n , mesh.BCp.type,  0, 0, interp, cent, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, instance->materials.aniso_factor,     &mesh, mesh.aniso_factor_s , mesh.BCg.type,  0, 0, interp, vert, instance->model.itp_stencil);
        }

        printf("*************************************\n");
        printf("******* Initialize viscosity ********\n");
        printf("*************************************\n");

        // Compute shear modulus, expansivity, compressibiulity, cohesion and friction angle on the grid
        if (instance->model.iselastic == 1 ) ShearModCompExpGrid( &mesh, instance->materials, instance->model, instance->scaling );
        CohesionFrictionDilationGrid( &mesh, &particles, instance->materials, instance->model, instance->scaling );
        ShearModCompExpGrid( &mesh, instance->materials, instance->model, instance->scaling );
        Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_in, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &instance->model );
        Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T,    mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &instance->model );
        NonNewtonianViscosityGrid (     &mesh, &instance->materials, &instance->model, Nmodel, &instance->scaling );

        // Print informations!
        printf("Number of phases : %d\n", instance->model.Nb_phases);
        if (instance->model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n, instance->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d0        ", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n, instance->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_lith, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P litho   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p0_n, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_in, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T, instance->scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s, instance->scaling.S,   (mesh.Nx-0)*(mesh.Nz-0), "mu_s      ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n      ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_s, instance->scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_n, instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_phys_s, instance->scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_phys_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_phys_n, instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho_s, instance->scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n, instance->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            for (int p=0; p< instance->model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }
        }

        printf("*************************************\n");
        printf("******** Initialize timestep ********\n");
        printf("*************************************\n");

        DefineInitialTimestep( &instance->model, &mesh, particles, instance->materials, instance->scaling );

        if (instance->model.rec_T_P_x_z == 1 ) {
            ArrayEqualArray( particles.T0,   particles.T, particles.Nb_part );
            ArrayEqualArray( particles.P0,   particles.P, particles.Nb_part );
            ArrayEqualArray( particles.x0,   particles.x, particles.Nb_part );
            ArrayEqualArray( particles.z0,   particles.z, particles.Nb_part );
            ArrayEqualArray( particles.Tmax, particles.T, particles.Nb_part );
            ArrayEqualArray( particles.Pmax, particles.P, particles.Nb_part );
        }
        ComputeMeanQuantitesForTimeSeries( &mesh );
        LogTimeSeries( &mesh, instance->model, instance->scaling );

        printf("*************************************\n");
        printf("*** Write initial file or restart ***\n");
        printf("*************************************\n");

        // Write initial output
        if ( writer == 1 ) {
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Output", instance->materials, instance->scaling );
            if (instance->model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particles", instance->materials, instance->scaling );
        }

        // Set initial stresses and pressure to zero
        Initialise1DArrayDouble( mesh.sxxd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.szzd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxz,   (mesh.Nx)  *(mesh.Nz)  , 0.0 );
        // Generate deformation maps
        if (instance->model.def_maps == 1 ) GenerateDeformationMaps( &mesh, &instance->materials, &instance->model, Nmodel, &instance->scaling );
        particles.Nb_part_ini = particles.Nb_part;
    }
    else {
        // Which step do we restart from (BreakpointXXXX.dat)
        instance->model.step = istep;
        printf("Restarting from step number %05d...\n", instance->model.step);
        LoadBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, &instance->model, &topo, &topo_ini, instance->scaling  );
        SetGridCoordinates( &mesh, &instance->model, instance->model.Nx, instance->model.Nz ); // Overwrite previous grid
    }

    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                             TIME LOOP : en avant les Doud'series !
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//

    instance->model.step += 1;

    if (instance->model.GNUplot_residuals == 1 ) GNUplotPipe = popen ("gnuplot -persistent", "w");

    for (; instance->model.step<= instance->model.Nt; instance->model.step++) {

        printf("*****************************************************\n");
        printf("****************** Time step %05d ******************\n", instance->model.step);
        printf("*****************************************************\n");

        //------------------------------------------------------------------------------------------------------------------------------//

        t_omp_step = (double)omp_get_wtime();

        // Old time step
        instance->model.dt0 = instance->model.dt;

        // Define new time step
        EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &instance->model, instance->scaling, &mesh, 0 );
        printf("Selected dt = %2.2e\n", instance->model.dt* instance->scaling.t);

        // Save initial dt
        instance->model.dt0 = instance->model.dt;

        // Track particule generation index
        Initialise1DArrayInt( particles.generation,   particles.Nb_part  , 0 );

        //------------------------------------------------------------------------------------------------------------------------------//

        // Remove particles that would be above the surface
        if (instance->model.free_surf == 1 ) {
            CleanUpSurfaceParticles( &particles, &mesh, topo, instance->scaling );
            CellFlagging( &mesh, instance->model, topo, instance->scaling );
            ArrayEqualArray(     topo.height0,     topo.height, mesh.Nx );
            ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );
        }

        // Interpolate material properties from particles to nodes
        t_omp = (double)omp_get_wtime();

        // Energy - interpolate thermal parameters and advected energy

        // Get energy and related material parameters from particles
        P2Mastah( &instance->model, particles, instance->materials.Cv,     &mesh, mesh.Cv,     mesh.BCp.type,  0, 0, interp, cent, 1);
        P2Mastah( &instance->model, particles, instance->materials.Qr,     &mesh, mesh.Qr,     mesh.BCp.type,  0, 0, interp, cent, 1);

        P2Mastah ( &instance->model, particles, instance->materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah ( &instance->model, particles, instance->materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);

        // Get T and dTdt from previous step from particles
        P2Mastah( &instance->model, particles, particles.T,     &mesh, mesh.T0_n,     mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &instance->model, particles, particles.divth, &mesh, mesh.divth0_n, mesh.BCp.type,  1, 0, interp, cent, 1);

        // Make sure T is up to date for rheology evaluation
        ArrayEqualArray( mesh.T, mesh.T0_n, (mesh.Nx-1)*(mesh.Nz-1) );

        //-----------------------------------------------------------------------------------------------------------
        // Interp P --> p0_n , p0_s
        P2Mastah( &instance->model, particles, particles.P,     &mesh, mesh.p0_n,   mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);

        // Get physical properties that are constant throughout each timestep
        UpdateDensity( &mesh, &particles, &instance->materials, &instance->model, &instance->scaling );

        // Free surface - subgrid density correction
        if (instance->model.free_surf == 1 ) {
            SurfaceDensityCorrection( &mesh, instance->model, topo, instance->scaling  );
        }

        // Lithostatic pressure
        ArrayEqualArray(  mesh.p_lith0,   mesh.p_lith, Ncx*Ncz );
        ComputeLithostaticPressure( &mesh, &instance->model, instance->materials.rho[0], instance->scaling, 1 );

        // Elasticity - interpolate advected/rotated stresses
        if  (instance->model.iselastic == 1 ) {

            // Get old stresses from particles
            if (instance->model.StressUpdate==1)                     OldDeviatoricStressesPressure( &mesh, &particles, instance->scaling, &instance->model );

            if (instance->model.StressUpdate==0){
                P2Mastah( &instance->model, particles, particles.sxxd,    &mesh, mesh.sxxd0, mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
                P2Mastah( &instance->model, particles, particles.szzd,    &mesh, mesh.szzd0, mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
                P2Mastah( &instance->model, particles, particles.sxz,     &mesh, mesh.sxz0,  mesh.BCg.type,  1, 0, interp, vert, instance->model.itp_stencil);
            }

            InterpCentroidsToVerticesDouble( mesh.sxxd0, mesh.sxxd0_s, &mesh, &instance->model );
            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &instance->model );
            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &instance->model );

            // Interpolate shear modulus
            ShearModCompExpGrid( &mesh, instance->materials, instance->model, instance->scaling );
        }

        // Director vector
        if (instance->model.aniso == 1 ) {
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d1_n , mesh.BCp.type, -1, 0, interp, cent, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d2_n , mesh.BCp.type, -2, 0, interp, cent, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d1_s , mesh.BCg.type, -1, 0, interp, vert, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, NULL, &mesh, mesh.d2_s , mesh.BCg.type, -2, 0, interp, vert, instance->model.itp_stencil);
            FiniteStrainAspectRatio ( &mesh, instance->scaling, instance->model, &particles );
            P2Mastah( &instance->model, particles, instance->materials.aniso_factor,     &mesh, mesh.aniso_factor_n , mesh.BCp.type,  0, 0, interp, cent, instance->model.itp_stencil);
            P2Mastah( &instance->model, particles, instance->materials.aniso_factor,     &mesh, mesh.aniso_factor_s , mesh.BCg.type,  0, 0, interp, vert, instance->model.itp_stencil);
        }

        P2Mastah( &instance->model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
        P2Mastah( &instance->model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, instance->model.itp_stencil);

        P2Mastah( &instance->model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, instance->model.itp_stencil);
        P2Mastah( &instance->model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);

        // Diffuse rheological contrasts
        //            if (instance->model.diffuse_X == 1) {
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.X0_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Diffuse_X(&mesh, &instance->model, &instance->scaling);
        //                Interp_Grid2P( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Interp_P2N ( particles, particles.X, &mesh, mesh.Xreac_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &instance->model );
        //            }

        // Interpolate Grain size
        P2Mastah( &instance->model, particles, particles.d,     &mesh, mesh.d0_n , mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);
        ArrayEqualArray(  mesh.d_n,  mesh.d0_n, Ncx*Ncz );

        // Interpolate Melt fraction
        P2Mastah( &instance->model, particles, particles.phi,   &mesh, mesh.phi0_n , mesh.BCp.type,  1, 0, interp, cent, instance->model.itp_stencil);

        //-------------------------------------------------------------------------------------------------------------

        // Compute cohesion and friction angle on the grid
        CohesionFrictionDilationGrid( &mesh, &particles, instance->materials, instance->model, instance->scaling );

        // Detect compressible cells
        if (instance->model.compressible == 1 ) DetectCompressibleCells ( &mesh, &instance->model );

        // Compute cohesion and friction angle on the grid
        CohesionFrictionDilationGrid( &mesh, &particles, instance->materials, instance->model, instance->scaling );

        // Detect compressible cells
        if (instance->model.compressible == 1) DetectCompressibleCells ( &mesh, &instance->model );

        // Min/Max interpolated fields
        if (instance->model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n, instance->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d0        ", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n, instance->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
            MinMaxArray(particles.sxxd, instance->scaling.S, particles.Nb_part, "sxxd part  ");
            MinMaxArray(particles.T, instance->scaling.T,   particles.Nb_part, "T part    ");
            MinMaxArrayTag( mesh.p0_n, instance->scaling.S,    (mesh.Nx-1)*(mesh.Nz-1),   "p0_n",   mesh.BCp.type );
            MinMaxArrayTag( mesh.sxz0, instance->scaling.S,   (mesh.Nx)*(mesh.Nz),     "sxz0    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.sxxd0, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "sxx0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.szzd0, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "szz0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s, instance->scaling.S,   (mesh.Nx)*(mesh.Nz),     "mu_s    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.C_s, instance->scaling.S,   (mesh.Nx)*(mesh.Nz),     "C_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.C_n, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "C_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.fric_s,   180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "fric_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.fric_n,   180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "fric_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.dil_s,    180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "dil_s   ", mesh.BCg.type );
            MinMaxArrayTag( mesh.dil_n,    180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "dil_n   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.strain_s,   1.0,       (mesh.Nx)*(mesh.Nz),     "strain_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.strain_n,   1.0,       (mesh.Nx-1)*(mesh.Nz-1), "strain_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.bet_s,    1.0/ instance->scaling.S,   (mesh.Nx)*(mesh.Nz),     "beta_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.bet_n,    1.0/ instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "beta_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T0_n, instance->scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T       ", mesh.BCt.type );
            MinMaxArrayTag( mesh.p_in, instance->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P       ", mesh.BCt.type );
            MinMaxArrayI  ( mesh.comp_cells, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "comp_cells" );
            MinMaxArrayTag( mesh.rho_s, instance->scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n, instance->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho0_n, instance->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );

            for (int p=0; p< instance->model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }

            if  (instance->model.aniso == 1 ) MinMaxArrayTag( mesh.FS_AR_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "FS_AR_n", mesh.BCp.type );
            if  (instance->model.aniso == 1 ) MinMaxArrayTag( mesh.FS_AR_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "FS_AR_s", mesh.BCg.type );
            if  (instance->model.aniso == 1 ) MinMaxArrayTag( mesh.aniso_factor_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "aniso_factor_n", mesh.BCp.type );
            if  (instance->model.aniso == 1 ) MinMaxArrayTag( mesh.aniso_factor_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "aniso_factor_s", mesh.BCg.type );
        }

        printf("** Time for particles interpolations I = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );

//        struct timespec begin, end;
//        clock_gettime(CLOCK_REALTIME, &begin);
//
//        ViscosityDerivatives( &mesh, &instance->materials, &instance->model, Nmodel, &instance->scaling );
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

        if (instance->model.ismechanical == 1 ) {

          if (instance->crazyConductivity) {
            for (int n = 0; n < instance->crazyConductivity->nPhases; n++) {
              const int    phase           = instance->crazyConductivity->phases[n];
              instance->materials.k_eff[phase] = instance->materials.k[phase];
            }
            printf("Running with normal conductivity for the asthenosphere...\n");
          }

            // Allocate and initialise solution and RHS vectors
            SetBCs(instance, &mesh);

            // Reset fields and BC values if needed
            //        if ( instance->model.ispureshear_ale == 1 ) InitialiseSolutionFields( &mesh, &instance->model );
            InitialiseSolutionFields( &mesh, &instance->model );
            EvalNumberOfEquations( &mesh, &Stokes );
            if ( IsNewtonStep == 1 ) EvalNumberOfEquations( &mesh, &Jacob  );
            SAlloc( &Stokes,  Stokes.neq );
            if ( IsNewtonStep == 1 ) SAlloc(  &Jacob,   Jacob.neq );
            if (instance->model.decoupled_solve == 1 ) {
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
            InitialiseSolutionVector( &mesh, &Stokes, &instance->model );

            //------------------------------------------------------------------------------------------------------------------------------//

            // Non-linear iteration cycle
            Nmodel.nit       = 0;
            instance->model.nit        = 0;
            Nmodel.stagnated = 0;
            nstag            = 0;
            IsJacobianUsed    = 0;

            ArrayEqualArray( mesh.p_start,    mesh.p_in,      (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.u_start,    mesh.u_in,      (mesh.Nx)  *(mesh.Nz+1) );
            ArrayEqualArray( mesh.v_start,    mesh.v_in,      (mesh.Nx+1)*(mesh.Nz)   );
            //            ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );

            // Set up solver context
            if (instance->model.decoupled_solve == 1 ) {
                cholmod_start( &CholmodSolver.c );
                if ( Nmodel.nit==0 ) CholmodSolver.Analyze = 1;
                printf("Run CHOLMOD analysis yes/no: %d \n", CholmodSolver.Analyze);
            }

            int Nmax_picard = Nmodel.nit_max;

            // If Picard 2 Newton is activated: force initial Newton step
            if ( IsNewtonStep == 1 && Nmodel.Picard2Newton == 1 ) {
              instance->model.Newton = 0;
                IsFirstNewtonStep  = 1;
            }

            while ( Nmodel.nit <= Nmax_picard && nstag< instance->model.nstagmax) {

                if ( Nmodel.nit > 0 && Nmodel.Picard2Newton == 1 ) {
                    if ( Nmodel.rx_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.rz_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.rp_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.nit>= Nmodel.nit_Pic_max) {
                        if ( IsFirstNewtonStep == 0 ) CholmodSolver.Analyze = 0;
                        if ( IsFirstNewtonStep == 1 ) {
                            cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                            CholmodSolver.Analyze = 1;
                            IsFirstNewtonStep    = 0;
                        }
                        instance->model.Newton = 1;
                    }
                }

                // Determine whether Jacobian matrix should be assembled
                if (instance->model.Newton == 1 || instance->model.aniso == 1 ) IsJacobianUsed = 1;

                printf("**********************************************\n");
                if (instance->model.Newton == 1 ) {
                    printf("*** Newton it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, instance->model.step);
                    Nmodel.LogIsNewtonStep[Nmodel.nit] = 1;
                }
                else {
                    printf("*** Picard it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, instance->model.step);
                    Nmodel.LogIsNewtonStep[Nmodel.nit] = 0;
                }
                printf("**********************************************\n");

                // Update non-linear rheology
                UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, instance->materials, &instance->model, &Nmodel, instance->scaling, 0, 0.0 );
                RheologicalOperators( &mesh, &instance->model, &instance->scaling, 0 );                               // ??????????? déjà fait dans UpdateNonLinearity
                NonNewtonianViscosityGrid (     &mesh, &instance->materials, &instance->model, Nmodel, &instance->scaling );    // ??????????? déjà fait dans UpdateNonLinearity

                if (instance->model.noisy == 1 ) {
                    MinMaxArrayTag( mesh.T, instance->scaling.T, (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCt.type );
                    MinMaxArrayTag( mesh.p0_n, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.p_in, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.div_u, instance->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "div       ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exxd, instance->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "exxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.ezzd, instance->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "ezzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exz, instance->scaling.E, (mesh.Nx-0)*(mesh.Nz-0), "exz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.sxxd, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "sxxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.szzd, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "szzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.sxz, instance->scaling.S, (mesh.Nx-0)*(mesh.Nz-0), "sxz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_s, instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_n, instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.eta_phys_s, instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_phys_s", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_phys_n, instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho_s, instance->scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.rho_n, instance->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho0_n, instance->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.d_n, instance->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
                }

                if (instance->model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Output_BeforeSolve", instance->materials, instance->scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particles_BeforeSolve", instance->materials, instance->scaling );
                }

                // Build discrete system of equations - Linearised Picard
                if (instance->model.decoupled_solve == 0 ) BuildStokesOperator           ( &mesh, instance->model, 0, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, 1 );
                if (instance->model.decoupled_solve == 1 ) BuildStokesOperatorDecoupled  ( &mesh, instance->model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1 );

                // Build discrete system of equations - Jacobian
                ViscosityDerivatives( &mesh, &instance->materials, &instance->model, Nmodel, &instance->scaling );
                if (instance->model.Newton == 1 && Nmodel.nit > 0 ) RheologicalOperators( &mesh, &instance->model, &instance->scaling, 1 );
                if ( IsJacobianUsed == 1 )                  BuildJacobianOperatorDecoupled( &mesh, instance->model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in,  &Jacob,  &JacobA,  &JacobB,  &JacobC,   &JacobD, 1 );

                // IsNanArray2DFP(mesh.eta_n, (mesh.Nx-1)*(mesh.Nz-1));
                // IsInfArray2DFP(mesh.eta_n, (mesh.Nx-1)*(mesh.Nz-1));

                // IsNanArray2DFP(mesh.eta_s, (mesh.Nx-0)*(mesh.Nz));
                // IsInfArray2DFP(mesh.eta_s, (mesh.Nx-0)*(mesh.Nz));
                //                MinMaxArrayTag( mesh.detadexx_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadexx_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadezz_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadezz_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadexx_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadgxz_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadexx_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadexx_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadezz_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadezz_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadgxz_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadgxz_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D11_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D21_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D12_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D22_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D13_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D23_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D14_n,      instance->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D24_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D31_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D31_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D32_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D32_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D33_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D33_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D34_s,      instance->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D34_s     ", mesh.BCg.type );

                // Diagonal instance->scaling
                if (instance->model.diag_scaling ) {
                    if (instance->model.Newton          == 0 )  ExtractDiagonalScale( &StokesA, &StokesB, &StokesC, &StokesD );
                    if (instance->model.Newton          == 1 )  ExtractDiagonalScale( &JacobA,  &JacobB,  &JacobC,   &JacobD );
                    if (instance->model.Newton          == 1 )  ArrayEqualArray(StokesA.d, JacobA.d, StokesA.neq);
                    if (instance->model.Newton          == 1 )  ArrayEqualArray(StokesC.d, JacobC.d, StokesC.neq);
                    ScaleMatrix( &StokesA, &StokesB, &StokesC, &StokesD );
                }

                // Needs to be done after matrix assembly since diagonal instance->scaling is used in there
                printf("---- Non-linear residual ----\n");
                RheologicalOperators( &mesh, &instance->model, &instance->scaling, 0 );
                if (instance->model.decoupled_solve == 0 ) EvaluateStokesResidual( &Stokes, &Nmodel, &mesh, instance->model, instance->scaling, 0 );
                if (instance->model.decoupled_solve == 1 ) EvaluateStokesResidualDecoupled( &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &Nmodel, &mesh, instance->model, instance->scaling, 0 );
                printf("---- Non-linear residual ----\n");

                if (instance->model.Newton == 1 && instance->model.aniso == 0 ) {
                    ArrayEqualArray( JacobA.F, StokesA.F,  StokesA.neq );
                    ArrayEqualArray( JacobC.F, StokesC.F,  StokesC.neq );
                }
                if (instance->model.Newton == 1 && instance->model.diag_scaling ) ScaleMatrix( &JacobA,  &JacobB,  &JacobC,  &JacobD  );

                // Store residuals
                Nmodel.resx_f = Nmodel.resx; Nmodel.rx_abs[Nmodel.nit] = Nmodel.resx; Nmodel.rx_rel[Nmodel.nit] = Nmodel.resx/Nmodel.resx0;
                Nmodel.resz_f = Nmodel.resz; Nmodel.rz_abs[Nmodel.nit] = Nmodel.resz; Nmodel.rz_rel[Nmodel.nit] = Nmodel.resz/Nmodel.resz0;
                Nmodel.resp_f = Nmodel.resp; Nmodel.rp_abs[Nmodel.nit] = Nmodel.resp; Nmodel.rp_rel[Nmodel.nit] = Nmodel.resp/Nmodel.resp0;

                if (instance->model.write_debug == 1 ) WriteResiduals( mesh, instance->model, Nmodel, instance->scaling );

                // if pass --> clear matrix break
                if ( (Nmodel.resx < Nmodel.abs_tol_u || Nmodel.resx/Nmodel.resx0 < Nmodel.rel_tol_u) && (Nmodel.resz < Nmodel.abs_tol_u || Nmodel.resz/Nmodel.resz0 < Nmodel.rel_tol_u) && (Nmodel.resp < Nmodel.abs_tol_p || Nmodel.resp/Nmodel.resp0 < Nmodel.rel_tol_p) ) {
                    printf( "Non-linear solver converged to abs_tol_u = %2.2e abs_tol_p = %2.2e rel_tol_u = %2.2e rel_tol_p = %2.2e\n", Nmodel.abs_tol_u, Nmodel.abs_tol_p, Nmodel.rel_tol_u, Nmodel.rel_tol_p );
                    FreeSparseSystems( IsJacobianUsed, instance->model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                    break;
                }

                // Direct solve
                t_omp = (double)omp_get_wtime();
                if (instance->model.decoupled_solve == 0 ) SolveStokesDefect( &Stokes, &CholmodSolver, &Nmodel, &mesh, &instance->model, &particles, &topo_chain, &topo, instance->materials, instance->scaling );
                if (instance->model.decoupled_solve == 1 ) {

                    if ( IsJacobianUsed==1 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &instance->model, &particles, &topo_chain, &topo, instance->materials, instance->scaling,  &JacobA,  &JacobB,  &JacobC );
                    else                     SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &instance->model, &particles, &topo_chain, &topo, instance->materials, instance->scaling, &StokesA, &StokesB, &StokesC );

                    if ( Nmodel.stagnated == 1 && instance->model.safe_mode <= 0 ) {
                        printf( "Non-linear solver stagnated to res_u = %2.2e res_z = %2.2e\n", Nmodel.resx_f, Nmodel.resz_f );
                        printf( "You may want to try setting line_search_min > 0.0\n Good luck good man!\n");
                        FreeSparseSystems( IsJacobianUsed, instance->model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                        break;
                    }
                }

                if ( Nmodel.stagnated == 1 && instance->model.iselastic == 1 && instance->model.safe_mode == 1 ) {
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (abs_tol_u = %2.2e abs_tol_p = %2.2e)\e[m\n", Nmodel.abs_tol_u, Nmodel.abs_tol_p );
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (rel_tol_u = %2.2e rel_tol_p = %2.2e)\e[m\n", Nmodel.rel_tol_u, Nmodel.rel_tol_p );
                    printf( "\e[1;31mReducing the timestep, and restart the iterations cycle...\e[m\n");
                    printf( "Before reduction: instance->model.dt =, %2.2e\n", instance->model.dt* instance->scaling.t);
                    // ----------------------
                    instance->model.dt /= instance->model.safe_dt_div;
                    printf( "Timestep divided by %2.2f => NEW CURRENT instance->model.dt =, %2.2e\n", instance->model.safe_dt_div, instance->model.dt* instance->scaling.t);
                    Nmodel.stagnated = 0;
                    // ----------------------
                    Nmodel.nit = -1;
                    printf( "Restart solutions\n");
                    ArrayEqualArray( mesh.p_in, mesh.p_start,  (mesh.Nx-1)*(mesh.Nz-1) );
                    ArrayEqualArray( mesh.u_in, mesh.u_start,  (mesh.Nx)  *(mesh.Nz+1) );
                    ArrayEqualArray( mesh.v_in, mesh.v_start,  (mesh.Nx+1)*(mesh.Nz)   );
                    nstag++;
                    //-----------------------
                    printf( "nstag value = %02d - nstagmax = %02d\n", nstag, instance->model.nstagmax);
                    if (nstag== instance->model.nstagmax) {
                        printf( "CheckDoudzOut!!\n");
                        exit(0);
                    //-----------------------
                    }
                }
                FreeSparseSystems( IsJacobianUsed, instance->model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                Nmodel.nit++;
                instance->model.nit = Nmodel.nit;
            }

            // Clean solver context
            if (instance->model.decoupled_solve == 1 && Nmodel.nit>0 ) {
                printf("Cleaning up Cholesky factors --- nit = %02d\n",  Nmodel.nit);
                cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                cholmod_finish( &CholmodSolver.c );
            }

            // Update rheology
            UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, instance->materials, &instance->model, &Nmodel, instance->scaling, 0, 1.0 );
            ArrayEqualArray( mesh.p_in, mesh.p_corr, (mesh.Nx-1)*(mesh.Nz-1) ); // make sure corrected pressure is used

            // Min/Max velocities
            MinMaxArray( mesh.u_in, instance->scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in, instance->scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in, instance->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
            MinMaxArray( mesh.div_u, instance->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "  div(V)" );

            printf("--------------------------------------------------------------\n");
            int i, nit;
            if (Nmodel.nit>=Nmodel.nit_max)  nit = Nmodel.nit_max;
            if (Nmodel.nit< Nmodel.nit_max)  nit = Nmodel.nit;
            if (Nmodel.Picard2Newton == 1 )  printf("Picard 2 Newton is activated with condition: %2.2e\n", Nmodel.Pic2NewtCond);
            for (i=0; i<=nit; i++) {
                if (Nmodel.LogIsNewtonStep[i] == 1) printf("New. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, Nmodel.rx_abs[i], Nmodel.rz_abs[i], Nmodel.rp_abs[i], Nmodel.rx_rel[i], Nmodel.rz_rel[i], Nmodel.rp_rel[i]);
                else                         printf("Pic. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, Nmodel.rx_abs[i], Nmodel.rz_abs[i], Nmodel.rp_abs[i], Nmodel.rx_rel[i], Nmodel.rz_rel[i], Nmodel.rp_rel[i]);
                if (i == Nmodel.nit_max && instance->model.safe_mode == 1) {
                    printf("Exit: Max iteration reached: Nmodel.nit_max = %02d! Check what you wanna do now...\n",Nmodel.nit_max);
                    if ( (Nmodel.resx < Nmodel.abs_tol_u) && (Nmodel.resz < Nmodel.abs_tol_u) && (Nmodel.resp < Nmodel.abs_tol_p) ) {}
                    else exit(1);
                }
            }
            printf("--------------------------------------------------------------\n");

            // plot residuals
            if (instance->model.GNUplot_residuals == 1 ) { // && instance->model.step % writer_step == 0
              PlotResiduals(Nmodel, GNUplotPipe);
            }

        }

        //------------------------------------------------------------------------------------------------------------------------------//

        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE --- Free surface - interpolate velocity components on the free surface
        if (instance->model.free_surf == 1 ) {
            SurfaceVelocity( &mesh, instance->model, &topo, &topo_chain, instance->scaling );
            MinMaxArray( topo_chain.Vx, instance->scaling.V, topo_chain.Nb_part,       "Vx surf." );
            MinMaxArray( topo_chain.Vz, instance->scaling.V, topo_chain.Nb_part,       "Vz surf." );
            MinMaxArray( topo_chain_ini.Vx, instance->scaling.V, topo_chain_ini.Nb_part,   "Vx surf. ini." );
            MinMaxArray( topo_chain_ini.Vz, instance->scaling.V, topo_chain_ini.Nb_part,   "Vz surf. ini." );
        }
        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE

        if (instance->model.StressUpdate==1) TotalStresses( &mesh, &particles, instance->scaling, &instance->model );

        // Update stresses on markers
        if (instance->model.StressUpdate==0) UpdateParticleStress(  &mesh, &particles, &instance->model, &instance->materials, &instance->scaling );

        // Update pressure on markers
        UpdateParticlePressure( &mesh, instance->scaling, instance->model, &particles, &instance->materials );
        //        Interp_Grid2P_centroids( particles, particles.P, &mesh, mesh.p_in, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &instance->model );
        UpdateParticleX( &mesh, instance->scaling, instance->model, &particles, &instance->materials );

        // Grain size evolution
        UpdateParticleGrainSize( &mesh, instance->scaling, instance->model, &particles, &instance->materials );

        if (instance->model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n   , instance->scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d0", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n    , instance->scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d ", mesh.BCp.type );
            MinMaxArrayPart( particles.d, instance->scaling.L, particles.Nb_part, "d on markers", particles.phase ) ;
        }

        // Update phi on the particles
        UpdateParticlePhi( &mesh, instance->scaling, instance->model, &particles, &instance->materials );


        //------------------------------------------------------------------------------------------------------------------------------//

        if (instance->model.isthermal == 1 ) {

            printf("*************************************\n");
            printf("*********** Thermal solver **********\n");
            printf("*************************************\n");

            t_omp = (double)omp_get_wtime();

            // Matrix assembly and direct solve
            EnergyDirectSolve( &mesh, instance->model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, instance->model.dt, instance->model.shear_heat, instance->model.adiab_heat, instance->scaling, 1 );
            MinMaxArray(particles.T, instance->scaling.T, particles.Nb_part, "T part. before UpdateParticleEnergy");

            // Update energy on particles
            UpdateParticleEnergy( &mesh, instance->scaling, instance->model, &particles, &instance->materials );
            MinMaxArray(particles.T, instance->scaling.T, particles.Nb_part, "T part. after UpdateParticleEnergy");

            printf("** Time for Thermal solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        }

        //--------------------------------------------------------------------------------------------------------------------------------//

        if (instance->model.ProgReac == 1 ) {

            printf("*************************************\n");
            printf("********** Chemical solver **********\n");
            printf("*************************************\n");

            P2Mastah ( &instance->model, particles, instance->materials.k_chem, &mesh, mesh.kc_x, mesh.BCu.type,  0, 0, interp, vxnodes, instance->model.itp_stencil);
            P2Mastah ( &instance->model, particles, instance->materials.k_chem, &mesh, mesh.kc_z, mesh.BCv.type,  0, 0, interp, vznodes, instance->model.itp_stencil);
            ChemicalDirectSolve( &mesh, instance->model, &particles, &instance->materials, instance->model.dt, instance->scaling );

        }
        //--------------------------------------------------------------------------------------------------------------------------------//

        // Update maximum pressure and temperature on markers
        if (instance->model.rec_T_P_x_z == 1 )  UpdateMaxPT(instance->scaling, instance->model, &particles );

        ComputeMeanQuantitesForTimeSeries( &mesh );
        LogTimeSeries( &mesh, instance->model, instance->scaling );

        //------------------------------------------------------------------------------------------------------------------------------//


        if (instance->model.advection == 1 ) {

            printf("*************************************\n");
            printf("************** Advection ************\n");
            printf("*************************************\n");

            t_omp = (double)omp_get_wtime();
            double dt_solve = instance->model.dt;
            int    nsub=1.0, isub;
            double dt_sub;
            EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &instance->model, instance->scaling, &mesh, 0 );

            if (instance->model.dt < 0.95*dt_solve ) { // if dt advection is lower than dt used for thermo-mechanical solve then split dt
                nsub   = ceil(dt_solve/ instance->model.dt);
                dt_sub = dt_solve / nsub;
                printf("dt advection = %2.9e --- dt_solve = %2.9e\n", instance->model.dt* instance->scaling.t, dt_solve* instance->scaling.t );
                printf("dt is %lf larger than dt_solve: need %d substeps of %2.2e s\n", dt_solve/ instance->model.dt, nsub , dt_sub* instance->scaling.t );
                printf("So: nsub*dt_sub = %2.2e for dt_solve = %2.2e\n", nsub*dt_sub* instance->scaling.t, dt_solve* instance->scaling.t);
                instance->model.dt = dt_sub;
            }
            else {
                nsub     = 1; // if dt advection is larger than dt used for thermo-mechanical solve then, forget it, and set dt to dt_solve
                instance->model.dt = dt_solve;
            }

            // Loop on substeps
            for (isub=0;isub<nsub;isub++) {

                printf("************** Advection step %03d of %03d: dtsub = %2.2e **************\n", isub, nsub, instance->model.dt* instance->scaling.t );

                // Advect domain boundaries
                if (instance->model.ispureshear_ale > 0 ) {
                    PureShearALE( &instance->model, &mesh, &topo_chain, instance->scaling );
                }

                if (instance->model.free_surf == 1 ) {
                    // Save old topo
                    ArrayEqualArray( topo_chain.z0, topo_chain.z, topo_chain.Nb_part ); // save old z
                    ArrayEqualArray( topo_chain_ini.z0, topo_chain_ini.z, topo_chain.Nb_part ); // save old z


                    // Advect free surface with RK4
                    RogerGuntherII( &topo_chain, instance->model, mesh, 1, instance->scaling );
                    RogerGuntherII( &topo_chain_ini, instance->model, mesh, 1, instance->scaling );
                }

                //                // Advect free surface
                //                if ( instance->model.free_surf == 1 ) {
                //                    AdvectFreeSurf( &topo_chain,     instance->model, instance->scaling );
                //                    AdvectFreeSurf( &topo_chain_ini, instance->model, instance->scaling );
                //                    MinMaxArray( topo_chain.z,      instance->scaling.L, topo_chain.Nb_part,       "z surf.     " );
                //                    MinMaxArray( topo_chain_ini.z,  instance->scaling.L, topo_chain_ini.Nb_part,   "z surf. ini." );
                //                }

                // Correction for particle inflow 0
                if (instance->model.ispureshear_ale == -1 && instance->model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, instance->model, topo, 0 );

                // Advect fluid particles
                RogerGuntherII( &particles, instance->model, mesh, 1, instance->scaling );

                // Correction for particle inflow 1
                if (instance->model.ispureshear_ale == -1 && instance->model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, instance->model, topo, 1 );

                // Update accumulated strain
                AccumulatedStrainII( &mesh, instance->scaling, instance->model, &particles,  mesh.xc_coord,  mesh.zc_coord, mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );

                // Update deformation gradient tensor components
                if (instance->model.fstrain == 1 ) DeformationGradient( mesh, instance->scaling, instance->model, &particles );

                if (instance->model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Output_BeforeSurfRemesh", instance->materials, instance->scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particles_BeforeSurfRemesh", instance->materials, instance->scaling );
                }

                if (instance->model.free_surf == 1 ) {

                    //                    for (int k=0;k<topo_chain.Nb_part;k++) {
                    //
                    //                        if (fabs(topo_chain.x[k]) < 0.5*instance->model.surf_Winc ) {
                    //                            topo_chain.z[k] -= instance->model.dt*instance->model.surf_Vinc;
                    //                        }
                    //                    }
                    //

                    if (instance->model.surf_remesh == 1 ) {
                        // Get current topography
                        ProjectTopography( &topo,     &topo_chain, instance->model, mesh, instance->scaling, mesh.xg_coord, 0 );
                        ProjectTopography( &topo_ini, &topo_chain_ini, instance->model, mesh, instance->scaling, mesh.xg_coord, 0 );
                        MarkerChainPolyFit( &topo,     &topo_chain, instance->model, mesh );
                        MarkerChainPolyFit( &topo_ini, &topo_chain_ini, instance->model, mesh );

                        // Remesh free surface I
                        RemeshMarkerChain( &topo_chain,     &topo, instance->model, instance->scaling, &mesh, 1 );
                        RemeshMarkerChain( &topo_chain_ini, &topo_ini, instance->model, instance->scaling, &mesh, 1 );
                    }


                    // Project topography on vertices
                    ProjectTopography( &topo,     &topo_chain, instance->model, mesh, instance->scaling, mesh.xg_coord, 0 );
                    ProjectTopography( &topo_ini, &topo_chain_ini, instance->model, mesh, instance->scaling, mesh.xg_coord, 0 );
                    //                    ArrayEqualArray( topo.height0, topo.height, mesh.Nx );
                    //                    ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );

                    if (instance->model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Output_AfterSurfRemesh", instance->materials, instance->scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particles_AfterSurfRemesh", instance->materials, instance->scaling );
                    }

                    // Diffuse topography
                    if (instance->model.surf_processes >= 1 )  DiffuseAlongTopography( &mesh, instance->model, instance->scaling, topo.height, topo.height, mesh.Nx, 0.0, instance->model.dt );

                    // Marker chain polynomial fit
                    MarkerChainPolyFit( &topo,     &topo_chain, instance->model, mesh );
                    CorrectTopoIni( &particles, instance->materials, &topo_chain_ini, &topo, instance->model, instance->scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, instance->model, mesh );

                    // Sedimentation
                    if (instance->model.surf_processes >= 2 ) {
                        AddPartSed( &particles, instance->materials, &topo_chain, &topo, instance->model, instance->scaling, &mesh);
                        if (instance->model.cpc==-1) CountPartCell_BEN( &particles, &mesh, instance->model, topo, 0, instance->scaling );
                        if (instance->model.cpc== 0) CountPartCell_Old( &particles, &mesh, instance->model, topo, 0, instance->scaling );
                        if (instance->model.cpc== 1) CountPartCell    ( &particles, &mesh, instance->model, topo, topo_ini, 0, instance->scaling );
                    }

                    if (instance->model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Outputx", instance->materials, instance->scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particlesx", instance->materials, instance->scaling );
                    }

                    // Remesh free surface II
                    RemeshMarkerChain( &topo_chain,     &topo, instance->model, instance->scaling, &mesh, 2 );
                    RemeshMarkerChain( &topo_chain_ini, &topo_ini, instance->model, instance->scaling, &mesh, 2 );
                    CorrectTopoIni( &particles, instance->materials, &topo_chain_ini, &topo, instance->model, instance->scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, instance->model, mesh );

                    // Remove particles that are above the surface
                    CleanUpSurfaceParticles( &particles, &mesh, topo, instance->scaling );

                    // Call cell flagging routine for free surface calculations
                    CellFlagging( &mesh, instance->model, topo, instance->scaling );
                }

                printf("** Time for advection solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

#ifdef _HDF5_
                if ( instance->model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Outputxx", instance->materials, instance->scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particlesxx", instance->materials, instance->scaling );
                }
#endif

                //                printf("*************************************\n");
                //                printf("************** Reseeding ************\n");
                //                printf("*************************************\n");

                // Count the number of particle per cell
                t_omp = (double)omp_get_wtime();
                if (instance->model.cpc ==-1)                       CountPartCell_BEN( &particles, &mesh, instance->model, topo, 0, instance->scaling );
                if (instance->model.cpc == 0)                       CountPartCell_Old( &particles, &mesh, instance->model, topo, 0, instance->scaling );
                if (instance->model.cpc == 1 && instance->model.Reseed == 1 ) CountPartCell    ( &particles, &mesh, instance->model, topo, topo_ini, 1, instance->scaling );
                if (instance->model.cpc == 1)                       CountPartCell    ( &particles, &mesh, instance->model, topo, topo_ini, 0, instance->scaling );

                if (instance->model.cpc == 2 && instance->model.Reseed == 1 ) CountPartCell2   ( &particles, &mesh, instance->model, topo, topo_ini, 1, instance->scaling );
                if (instance->model.cpc == 2)                       CountPartCell2   ( &particles, &mesh, instance->model, topo, topo_ini, 0, instance->scaling );

                printf("After re-seeding :\n");
                printf("Initial number of particles = %d\n", particles.Nb_part_ini);
                printf("New number of particles     = %d\n", particles.Nb_part    );

                printf("** Time for CountPartCell = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

                // Remove particles that would be above the surface
                if (instance->model.free_surf == 1 ) {
                    CleanUpSurfaceParticles( &particles, &mesh, topo, instance->scaling );
                    CellFlagging( &mesh, instance->model, topo, instance->scaling );
                }

                //                for (int iphase=0; iphase<instance->model.Nb_phases; iphase++) {
                //                CheckSym( mesh.phase_perc_n[iphase], 1.0, mesh.Nx-1, mesh.Nz-1, "perc_n" );
                //                CheckSym( mesh.phase_perc_s[iphase], 1.0, mesh.Nx-0, mesh.Nz-0, "perc_s" );
                //                }
            }
            instance->model.dt = dt_solve;
        }

        // Update time
        instance->model.time += instance->model.dt;

        // Free solution arrays
        SFree( &Stokes );
        if ( IsNewtonStep == 1 ) SFree( &Jacob );
        if (instance->model.decoupled_solve == 1 ) {
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
        if ( writer == 1 && instance->model.step % writer_step == 0 ) {

            printf("*************************************\n");
            printf("********* Write output files ********\n");
            printf("*************************************\n");

            // Breakpoint file
            t_omp = (double)omp_get_wtime();
            if (instance->model.delete_breakpoints == 1 ) DeletePreviousBreakpoint(instance->model.step, writer_step  );
            MakeBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, instance->model, &topo, &topo_ini, instance->scaling );
            UpdateInputFile( fin_name, instance->model.step);
            printf("** Time for Breakpoint file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

            // Visualisation file
#ifndef _VG_
            t_omp = (double)omp_get_wtime();
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, instance->model, Nmodel, "Output", instance->materials, instance->scaling );
            if (instance->model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, instance->model, "Particles", instance->materials, instance->scaling );
            printf("** Time for Output file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
#endif
        }

        printf("** Total timestep calculation time = %lf sec\n", (double)((double)omp_get_wtime() - t_omp_step) );
        printf("** Model time = %2.2e sec\n", instance->model.time* instance->scaling.t );
        printf("** Current dt = %2.2e sec, Old dt = %2.2e sec\n", instance->model.dt* instance->scaling.t, instance->model.dt0* instance->scaling.t );

        // if (Np0 - particles.Nb_part != 0) exit (1);


    }

    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                           END TIME LOOP : c'est fini les Doud'series...
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//

    if (instance->model.GNUplot_residuals == 1) pclose(GNUplotPipe);

    //    // Free markers chains
    if (instance->model.free_surf == 1 )    FreeMarkerChain( &topo,     &topo_chain     );
    if (instance->model.free_surf == 1 )    FreeMarkerChain( &topo_ini, &topo_chain_ini );

    // Free Phase diagrams if activated
    if (instance->model.isPD == 1 ) FreePhaseDiagrams( &instance->model );
    //
    // Free particle arrays
    PartFree( &particles, &instance->model );

    // Free arrays
    GridFree( &mesh, &instance->model );

    // Free char*'s
    free(instance->model.input_file);
    free(fin_name);

    // GNU plot
    DoodzFree( Nmodel.rx_abs );
    DoodzFree( Nmodel.rz_abs );
    DoodzFree( Nmodel.rp_abs );
    DoodzFree( Nmodel.rx_rel );
    DoodzFree( Nmodel.rz_rel );
    DoodzFree( Nmodel.rp_rel );
    DoodzFree( Nmodel.LogIsNewtonStep );

    // just for the quartz-coesite case
    for ( int k=0; k< instance->model.Nb_phases; k++) {

        if (instance->materials.density_model[k] == 4 ) {
            printf("Unloading 1D diagram for coesite quartz only baby...\n");
            DoodzFree(instance->model.PD1Drho[k]);
        }
    }
    DoodzFree(instance->model.PD1DnP);
    DoodzFree(instance->model.PD1Drho);
    DoodzFree(instance->model.PD1Dmin);
    DoodzFree(instance->model.PD1Dmax);
    if (instance->model.kinetics==1 ) DoodzFree(instance->model.kin_dG);

    printf("\n********************************************************\n");
    printf("************* Ending MDOODZ 7.0 simulation *************\n");
    printf("********************************************************\n");
}