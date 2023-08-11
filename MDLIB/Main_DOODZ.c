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

void RunMDOODZ(char *inputFileName, MdoodzSetup *setup) {

    printf("\n********************************************************\n");
    printf("************ Starting MDOODZ 7.0 simulation ************\n");
    printf("********************************************************\n");

    Input       inputFile   = ReadInputFile(inputFileName);

    char *BaseOutputFileName, *BaseParticleFileName;
    asprintf(&BaseOutputFileName, "Output");
    asprintf(&BaseParticleFileName, "Particles");

    MdoodzInput input = (MdoodzInput) {
      .model = inputFile.model,
      .scaling = inputFile.scaling,
      .materials = inputFile.materials,
      .crazyConductivity = NULL,
    };

    if (setup->MutateInput) {
      setup->MutateInput(&input, setup->mutateInputParams);
    }

    // ValidateSetup(setup, &input);
    input.model.L0 = input.model.xmax - input.model.xmin;

    printf("*************************************\n");
    printf("****** Allocate and initialise ******\n");
    printf("*************************************\n");

    // Multi resolution allocation pointers
    grid mesh = GridAlloc(&input.model );

    // Initialise grid coordinates
    SetGridCoordinates( &mesh, &input.model, input.model.Nx, input.model.Nz );

    // Get grid indices
    GridIndices( &mesh );

    // Initialise data for logs of iterations
    Nparams Nmodel = NmodelAlloc(inputFile.Nmodel);

    bool NeedJacobian = input.model.Newton == 1 || input.model.anisotropy == 1;

    int IsFullNewton = input.model.Newton == 1;
    if (input.model.anisotropy  == 1) input.model.Newton = 1; // activate Newton context is anisotropy is activated
    int IsNewtonStep = input.model.Newton == 1;

    int num_threads = 1;
    #pragma omp parallel
    {
        num_threads = omp_get_num_threads();
    }

    printf("Num threads = %03d\n", omp_get_num_threads());

    printf("*************************************\n");
    printf("******* Initialize particles ********\n");
    printf("*************************************\n");

    // Initialize and allocate particle fields
    markers particles = PartAlloc(inputFile.particles, &input.model );

    // Allocate marker chain
    markers       topo_chain, topo_chain_ini;
    surface      topo, topo_ini;
    if (input.model.free_surface == 1 ) AllocateMarkerChain( &topo,     &topo_chain, input.model );
    if (input.model.free_surface == 1 ) AllocateMarkerChain( &topo_ini, &topo_chain_ini, input.model );

    // Set new particle distribution
    int cent=1, vert=0, prop=1, interp=0, vxnodes=-1, vznodes=-2;
    if ( input.model.irestart == 0 ) {

      input.model.step = 0;
      input.model.time = 0.0;

        // Initialise particle fields
        PartInit( &particles, &input.model );

        if (input.crazyConductivity) {
          for (int n = 0; n < input.crazyConductivity->nPhases; n++) {
            const int    phase           = input.crazyConductivity->phases[n];
            const double k_crazy         = input.crazyConductivity->multiplier * input.materials.k[phase];
            input.materials.k_eff[phase] = k_crazy;
          }
          printf("Running with crazy conductivity for the asthenosphere!!\n");
        }

        // Initial grid tags
        SetBCs(*setup->SetBCs, &input, &mesh);
        if (input.model.free_surface == 1 ) {

            // Define the horizontal position of the surface marker chain
            SetTopoChainHorizontalCoords( &topo,     &topo_chain, input.model, mesh, input.scaling );
            SetTopoChainHorizontalCoords( &topo_ini, &topo_chain_ini, input.model, mesh, input.scaling );

            // Define the vertical position of the surface marker chain
            BuildInitialTopography(*setup->BuildInitialTopography, &input, &topo_chain );
            BuildInitialTopography(*setup->BuildInitialTopography, &input, &topo_chain_ini );

            // Project topography on vertices
            ProjectTopography( &topo, &topo_chain, input.model, mesh, input.scaling, mesh.xg_coord, 0 );

            // Marker chain polynomial fit
            MarkerChainPolyFit( &topo, &topo_chain, input.model, mesh );

            // Call cell flagging routine for free surface calculations
            CellFlagging( &mesh, input.model, topo, input.scaling );
        }

        // Set particles coordinates
        PutPartInBox( &particles, &mesh, input.model, topo, input.scaling );

        // Set phases on particles
        SetParticles(*setup->SetParticles, &input, &particles);

        MinMaxArray(particles.Vx, input.scaling.V, particles.Nb_part, "Vxp init" );
        MinMaxArray(particles.Vz, input.scaling.V, particles.Nb_part, "Vzp init" );
        MinMaxArray(particles.T, input.scaling.T, particles.Nb_part,  "Tp init" );

        if (input.model.free_surface == 1 ) CleanUpSurfaceParticles( &particles, &mesh, topo, input.scaling );

        // Create phase percentage arrays
        CountPartCell( &particles, &mesh, input.model, topo, topo_ini, 0, input.scaling  );

        P2Mastah( &input.model, particles, input.materials.eta0, &mesh, mesh.eta_s, mesh.BCg.type,  0, 0, interp, vert, input.model.interp_stencil);
        P2Mastah( &input.model, particles, input.materials.eta0, &mesh, mesh.eta_n, mesh.BCp.type,  0, 0, interp, cent, input.model.interp_stencil);

        P2Mastah( &input.model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, input.model.interp_stencil);
        P2Mastah( &input.model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);

        P2Mastah( &input.model, particles, particles.T,  &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &input.model, particles, input.materials.Cv, &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);

        P2Mastah( &input.model, particles, input.materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah( &input.model, particles, input.materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);

        if (input.model.noisy == 1 ) {
            MinMaxArray( mesh.u_in, input.scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in, input.scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
        }
        if (input.crazyConductivity) {
          for (int n = 0; n < input.crazyConductivity->nPhases; n++) {
            const int    phase           = input.crazyConductivity->phases[n];
            const double k_crazy         = input.crazyConductivity->multiplier * input.materials.k[phase];
            input.materials.k_eff[phase] = k_crazy;
          }
          printf("Running with crazy conductivity for the asthenosphere!!\n");
        }
        // Initial solution fields
        SetBCs(*setup->SetBCs, &input, &mesh);

        if (input.model.mechanical == 1) {
            InitialiseSolutionFields( &mesh, &input.model );
        }
        else {
            P2Mastah( &input.model, particles, particles.Vx, &mesh, mesh.u_in, mesh.BCu.type,  1, 0, interp, vxnodes, 1);
            P2Mastah( &input.model, particles, particles.Vz, &mesh, mesh.v_in, mesh.BCv.type,  1, 0, interp, vznodes, 1);
            ApplyBC( &mesh, &input.model );
        }

        MinMaxArray( mesh.u_in, input.scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
        MinMaxArray( mesh.v_in, input.scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
        MinMaxArray( mesh.p_in, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );

        printf("*************************************\n");
        printf("****** Initialize temperature *******\n");
        printf("*************************************\n");

        // Print2DArrayChar( mesh.BCt.type, mesh.Nx-1, mesh.Nz-1, 1.0 );
        // printf("West");
        // Print2DArrayChar( mesh.BCt.typE, mesh.Nz-1, 1, 1.0 );
        // printf("East");
        // Print2DArrayChar( mesh.BCt.typW, mesh.Nz-1, 1, 1.0 );
        // printf("South");
        // Print2DArrayChar( mesh.BCt.typS, mesh.Nx-1, 1, 1.0 );
        // printf("North");
        // Print2DArrayChar( mesh.BCt.typN, mesh.Nx-1, 1, 1.0 );
        // printf("BCT_exp");
        // Print2DArrayChar( mesh.BCC_exp.type, mesh.Nx+1, mesh.Nz+1, 1.0 );
        // Print2DArrayDouble( mesh.BCT_exp.val, mesh.Nx+1, mesh.Nz+1, input.scaling.T );

        // exit(33);

        // Get energy and related material parameters from particles
        P2Mastah( &input.model, particles, input.materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah( &input.model, particles, input.materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);
        P2Mastah( &input.model, particles, particles.T,     &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &input.model, particles, input.materials.Cv,    &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);
        P2Mastah( &input.model, particles, input.materials.Qr,    &mesh, mesh.Qr, mesh.BCp.type,  0, 0, interp, cent, 1);

        SetBCs(*setup->SetBCs, &input, &mesh);
        if (input.model.initial_cooling == 1 ) ThermalSteps( &mesh, input.model,  mesh.rhs_t, &particles, input.model.cooling_duration, input.scaling );
        if (input.model.therm_perturb == 1 ) SetThermalPert( &mesh, input.model, input.scaling );
        Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &input.model );
        ArrayEqualArray( mesh.T0_n, mesh.T, (mesh.Nx-1)*(mesh.Nz-1) );

        //--------------------------------------------------------------------------------------------------------

        printf("*************************************\n");
        printf("******** Initialize pressure ********\n");
        printf("*************************************\n");

        // Compute deviatoric strain rate tensor components
        ApplyBC( &mesh, &input.model );
        StrainRateComponents( &mesh, input.scaling, &input.model );

        // Iterate on pressure and density
        for (int iter=0; iter<10; iter++) {

            // Compute density
            UpdateDensity( &mesh, &particles, &input.materials, &input.model, &input.scaling );
            if (input.model.free_surface == 1 ) SurfaceDensityCorrection( &mesh, input.model, topo, input.scaling  );

            // Compute pressure
            ComputeLithostaticPressure( &mesh, &input.model, input.materials.rho[0], input.scaling, 1 );
            Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_lith, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &input.model );
            P2Mastah( &input.model, particles, particles.P,     &mesh, mesh.p_in , mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
            ArrayEqualArray( mesh.p_in, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.p0_n, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
            MinMaxArrayTag( mesh.p_in, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P initial ", mesh.BCp.type );
        }

        // Interpolate pressure from centroids to vertices
        InterpCentroidsToVerticesDouble( mesh.p0_n, mesh.p0_s, &mesh, &input.model );
  
        printf("*************************************\n");
        printf("******* Initialize grain size *******\n");
        printf("*************************************\n");
        InitialiseGrainSizeParticles( &particles, &input.materials );
        P2Mastah( &input.model, particles, particles.d,     &mesh, mesh.d_n , mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
        ArrayEqualArray( mesh.d0_n, mesh.d_n,  (mesh.Nx-1)*(mesh.Nz-1) );

        printf("*************************************\n");
        printf("******** Initialize density *********\n");
        printf("*************************************\n");
        UpdateDensity( &mesh, &particles, &input.materials, &input.model, &input.scaling );
        if (input.model.free_surface == 1 ) SurfaceDensityCorrection( &mesh, input.model, topo, input.scaling  );

        printf("*************************************\n");
        printf("****** Initialize composition *******\n");
        printf("*************************************\n");
        P2Mastah( &input.model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
        P2Mastah( &input.model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, input.model.interp_stencil);
        ArrayEqualArray( mesh.X_n, mesh.X0_n,  (mesh.Nx-1)*(mesh.Nz-1) );

        if (input.model.anisotropy == 1 ) {
            InitialiseDirectorVector ( &mesh, &particles, &input.model, &input.materials, input.scaling.L);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d1_n,    mesh.BCp.type, -1, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d2_n,    mesh.BCp.type, -2, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.angle_n, mesh.BCp.type, -3, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d1_s,    mesh.BCg.type, -1, 0, interp, vert, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d2_s,    mesh.BCg.type, -2, 0, interp, vert, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.angle_s, mesh.BCg.type, -3, 0, interp, vert, input.model.interp_stencil);
            FiniteStrainAspectRatio ( &mesh, input.scaling, input.model, &particles );
        }

        printf("*************************************\n");
        printf("******* Initialize viscosity ********\n");
        printf("*************************************\n");

        // Compute shear modulus, expansivity, compressibility, cohesion and friction angle on the grid
        if (input.model.elastic == 1 ) ShearModCompExpGrid( &mesh, &input.materials, &input.model, input.scaling );
        CohesionFrictionDilationGrid( &mesh, &particles, input.materials, input.model, input.scaling );
        ShearModCompExpGrid( &mesh, &input.materials, &input.model, input.scaling );
        Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_in, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &input.model );
        Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T,    mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &input.model );
        if (input.model.anisotropy==0) NonNewtonianViscosityGrid(      &mesh, &input.materials, &input.model, Nmodel, &input.scaling, 0 );
        if (input.model.anisotropy==1) NonNewtonianViscosityGridAniso( &mesh, &input.materials, &input.model, Nmodel, &input.scaling, 0 );

        // Print informations!
        printf("Number of phases : %d\n", input.model.Nb_phases);
        if (input.model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n, input.scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d0        ", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n, input.scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_lith, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P litho   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p0_n, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_in, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T, input.scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s, input.scaling.S,   (mesh.Nx-0)*(mesh.Nz-0), "mu_s      ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n      ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_s, input.scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_n, input.scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_phys_s, input.scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_phys_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_phys_n, input.scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho_s, input.scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n, input.scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.X_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1),              "X_n       ", mesh.BCp.type );
            MinMaxArrayTag( mesh.X0_n, 1.0,   (mesh.Nx-1)*(mesh.Nz-1),              "X0_n      ", mesh.BCp.type );
            MinMaxArrayTag( mesh.X0_s, 1.0,   (mesh.Nx-0)*(mesh.Nz-0),              "X0_s      ", mesh.BCg.type );
            MinMaxArray(particles.X, 1.0, particles.Nb_part, "X part" );
            if ( input.model.anisotropy == 1 ) {
                MinMaxArrayTag( mesh.FS_AR_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1),      "FS_AR_n   ", mesh.BCp.type );
                MinMaxArrayTag( mesh.FS_AR_s,  1.0,   (mesh.Nx)*(mesh.Nz),          "FS_AR_s   ", mesh.BCg.type );
                MinMaxArrayTag( mesh.angle_n,  180/M_PI,   (mesh.Nx-1)*(mesh.Nz-1), "angle_n   ", mesh.BCp.type );
                MinMaxArrayTag( mesh.angle_s,  180/M_PI,   (mesh.Nx)*(mesh.Nz),     "angle_s   ", mesh.BCg.type );
                MinMaxArray(particles.nx, 1.0, particles.Nb_part,                   "nx        " );
                MinMaxArray(particles.nz, 1.0, particles.Nb_part,                   "nz        " );
            }
            if (input.model.marker_noise == 1) MinMaxArrayTag( mesh.noise_s, 1.0, (mesh.Nx-0)*(mesh.Nz-0), "noise_s     ", mesh.BCg.type );
            if (input.model.marker_noise == 1) MinMaxArrayTag( mesh.noise_n, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "noise_n     ", mesh.BCp.type );
            for (int p=0; p< input.model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }
        }

        printf("*************************************\n");
        printf("******** Initialize timestep ********\n");
        printf("*************************************\n");

        DefineInitialTimestep( &input.model, &mesh, particles, input.materials, input.scaling );

        if (input.model.track_T_P_x_z == 1 ) {
            ArrayEqualArray( particles.T0,   particles.T, particles.Nb_part );
            ArrayEqualArray( particles.P0,   particles.P, particles.Nb_part );
            ArrayEqualArray( particles.x0,   particles.x, particles.Nb_part );
            ArrayEqualArray( particles.z0,   particles.z, particles.Nb_part );
            ArrayEqualArray( particles.Tmax, particles.T, particles.Nb_part );
            ArrayEqualArray( particles.Pmax, particles.P, particles.Nb_part );
        }
        ComputeMeanQuantitesForTimeSeries( &mesh );
        LogTimeSeries( &mesh, input.model, input.scaling );

        printf("*************************************\n");
        printf("*** Write initial file or restart ***\n");
        printf("*************************************\n");

        // Write initial output
        if ( input.model.writer == 1 ) {
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, BaseOutputFileName, input.materials, input.scaling );
            if (input.model.writer_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, BaseParticleFileName, input.materials, input.scaling );
        }

        // Set initial stresses and pressure to zero
        Initialise1DArrayDouble( mesh.sxxd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.szzd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxz,   (mesh.Nx)  *(mesh.Nz)  , 0.0 );
        // Generate deformation maps
        if (input.model.deformation_maps == 1 ) GenerateDeformationMaps( &mesh, &input.materials, &input.model, Nmodel, &input.scaling );
        particles.Nb_part_ini = particles.Nb_part;
    }
    else {
        // Which step do we restart from (BreakpointXXXX.dat)
        input.model.step = input.model.istep;
        printf("Restarting from step number %05d...\n", input.model.step);
        LoadBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, &input.model, &topo, &topo_ini, input.scaling  );
        SetGridCoordinates( &mesh, &input.model, input.model.Nx, input.model.Nz ); // Overwrite previous grid
    }

    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                             TIME LOOP : en avant les Doud'series !
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//

    input.model.step += 1;

    FILE        *GNUplotPipe;
    if (input.model.gnuplot_log_res == 1 ) GNUplotPipe = popen ("gnuplot -persistent", "w");

    SparseMat    Stokes, Jacob;
    SparseMat    StokesA, StokesB, StokesC, StokesD;
    SparseMat    JacobA,  JacobB,  JacobC,  JacobD;
    DirectSolver CholmodSolver;
    for (; input.model.step<= input.model.Nt; input.model.step++) {

        printf("*****************************************************\n");
        printf("****************** Time step %05d ******************\n", input.model.step);
        printf("*****************************************************\n");

        //------------------------------------------------------------------------------------------------------------------------------//

        clock_t t_omp_step = (double)omp_get_wtime();

        printf(GREEN "Number of particles     = %d\n" RESET, particles.Nb_part    );

        // Old time step
        input.model.dt0 = input.model.dt;

        // Define new time step
        EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &input.model, input.scaling, &mesh, 0 );
        input.model.time += input.model.dt;
        printf("Selected dt = %2.2e\n", input.model.dt* input.scaling.t);

        // Save initial dt
        input.model.dt0 = input.model.dt;

        // Track particule generation index
        Initialise1DArrayInt( particles.generation,   particles.Nb_part  , 0 );

        //------------------------------------------------------------------------------------------------------------------------------//

        // Remove particles that would be above the surface
        if (input.model.free_surface == 1 ) {
            CleanUpSurfaceParticles( &particles, &mesh, topo, input.scaling );
            CellFlagging( &mesh, input.model, topo, input.scaling );
            ArrayEqualArray(     topo.height0,     topo.height, mesh.Nx );
            ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );
            ArrayEqualArray( topo_chain.z0, topo_chain.z, topo_chain.Nb_part ); // save old z
            ArrayEqualArray( topo_chain_ini.z0, topo_chain_ini.z, topo_chain.Nb_part ); // save old z
        }

        // Interpolate material properties from particles to nodes
        clock_t t_omp = (double)omp_get_wtime();

        // Energy - interpolate thermal parameters and advected energy

        // Get energy and related material parameters from particles
        P2Mastah( &input.model, particles, input.materials.Cv,     &mesh, mesh.Cv,     mesh.BCp.type,  0, 0, interp, cent, 1);
        P2Mastah( &input.model, particles, input.materials.Qr,     &mesh, mesh.Qr,     mesh.BCp.type,  0, 0, interp, cent, 1);

        P2Mastah ( &input.model, particles, input.materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah ( &input.model, particles, input.materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);

        // Get T and dTdt from previous step from particles
        P2Mastah( &input.model, particles, particles.T,     &mesh, mesh.T0_n,     mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &input.model, particles, particles.divth, &mesh, mesh.divth0_n, mesh.BCp.type,  1, 0, interp, cent, 1);

        // Make sure T is up to date for rheology evaluation
        ArrayEqualArray( mesh.T, mesh.T0_n, (mesh.Nx-1)*(mesh.Nz-1) );

        //-----------------------------------------------------------------------------------------------------------
        // Interp P --> p0_n , p0_s
        P2Mastah( &input.model, particles, particles.P,     &mesh, mesh.p0_n,   mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);

        // Allocate and initialise solution and RHS vectors
        // SetBCs(*setup->SetBCs, &input, &mesh);
        UpdateDensity( &mesh, &particles, &input.materials, &input.model, &input.scaling );

        // Free surface - subgrid density correction
        if (input.model.free_surface == 1 ) { // TODO: include into UpdateDensity
            SurfaceDensityCorrection( &mesh, input.model, topo, input.scaling  );
        }

        // Lithostatic pressure
        int          Ncx = mesh.Nx - 1 , Ncz =  mesh.Nz - 1;
        ArrayEqualArray(  mesh.p_lith0,   mesh.p_lith, Ncx*Ncz );
        ComputeLithostaticPressure( &mesh, &input.model, input.materials.rho[0], input.scaling, 1 );

        // Elasticity - interpolate advected/rotated stresses
        if  (input.model.elastic == 1 ) {

            // Get old stresses from particles
            P2Mastah( &input.model, particles, particles.sxxd,    &mesh, mesh.sxxd0, mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, particles.szzd,    &mesh, mesh.szzd0, mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, particles.sxz,     &mesh, mesh.sxz0,  mesh.BCg.type,  1, 0, interp, vert, input.model.interp_stencil);

            InterpCentroidsToVerticesDouble( mesh.sxxd0, mesh.sxxd0_s, &mesh, &input.model );
            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &input.model );
            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &input.model );

            // Interpolate shear modulus
            ShearModCompExpGrid( &mesh, &input.materials, &input.model, input.scaling );
        }

        // Director vector
        if (input.model.anisotropy == 1 ) {
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d1_n,    mesh.BCp.type, -1, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d2_n,    mesh.BCp.type, -2, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.angle_n, mesh.BCp.type, -3, 0, interp, cent, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d1_s,    mesh.BCg.type, -1, 0, interp, vert, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.d2_s,    mesh.BCg.type, -2, 0, interp, vert, input.model.interp_stencil);
            P2Mastah( &input.model, particles, NULL, &mesh, mesh.angle_s, mesh.BCg.type, -3, 0, interp, vert, input.model.interp_stencil);
            FiniteStrainAspectRatio ( &mesh, input.scaling, input.model, &particles );
        }

        P2Mastah( &input.model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
        P2Mastah( &input.model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, input.model.interp_stencil);

        P2Mastah( &input.model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, input.model.interp_stencil);
        P2Mastah( &input.model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);

        // Interpolate Grain size
        P2Mastah( &input.model, particles, particles.d,     &mesh, mesh.d0_n , mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);
        ArrayEqualArray(  mesh.d_n,  mesh.d0_n, Ncx*Ncz );

        // Interpolate Melt fraction
        P2Mastah( &input.model, particles, particles.phi,   &mesh, mesh.phi0_n , mesh.BCp.type,  1, 0, interp, cent, input.model.interp_stencil);

        //-------------------------------------------------------------------------------------------------------------

        // Compute cohesion and friction angle on the grid
        CohesionFrictionDilationGrid( &mesh, &particles, input.materials, input.model, input.scaling );

        // Detect compressible cells
        if (input.model.compressible == 1 ) DetectCompressibleCells ( &mesh, &input.model );

        // Get physical properties that are constant throughout each timestep
        UpdateDensity( &mesh, &particles, &input.materials, &input.model, &input.scaling );

        // Free surface - subgrid density correction
        if (input.model.free_surface == 1 ) {
            SurfaceDensityCorrection( &mesh, input.model, topo, input.scaling  );
        }

        // Min/Max interpolated fields
        if (input.model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n, input.scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d0        ", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n, input.scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
            MinMaxArray(particles.sxxd, input.scaling.S, particles.Nb_part, "sxxd part  ");
            MinMaxArray(particles.T, input.scaling.T,   particles.Nb_part, "T part    ");
            MinMaxArrayTag( mesh.p0_n, input.scaling.S,    (mesh.Nx-1)*(mesh.Nz-1),   "p0_n",   mesh.BCp.type );
            MinMaxArrayTag( mesh.sxz0, input.scaling.S,   (mesh.Nx)*(mesh.Nz),     "sxz0    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.sxxd0, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "sxx0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.szzd0, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "szz0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s, input.scaling.S,   (mesh.Nx)*(mesh.Nz),     "mu_s    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.C_s, input.scaling.S,   (mesh.Nx)*(mesh.Nz),     "C_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.C_n, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "C_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.fric_s,   180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "fric_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.fric_n,   180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "fric_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.dil_s,    180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "dil_s   ", mesh.BCg.type );
            MinMaxArrayTag( mesh.dil_n,    180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "dil_n   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.strain_s,   1.0,       (mesh.Nx)*(mesh.Nz),     "strain_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.strain_n,   1.0,       (mesh.Nx-1)*(mesh.Nz-1), "strain_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.bet_s,    1.0/ input.scaling.S,   (mesh.Nx)*(mesh.Nz),     "beta_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.bet_n,    1.0/ input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "beta_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T0_n, input.scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T0      ", mesh.BCt.type );
            MinMaxArrayTag( mesh.T,    input.scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T       ", mesh.BCt.type );
            MinMaxArrayTag( mesh.p_in, input.scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P       ", mesh.BCt.type );
            MinMaxArrayI  ( mesh.comp_cells, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "comp_cells" );
            MinMaxArrayTag( mesh.rho_s, input.scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n, input.scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho0_n, input.scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.X0_s, 1.0, (mesh.Nx)*(mesh.Nz),     "X0_s      ", mesh.BCg.type );
            MinMaxArrayTag( mesh.X0_n, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "X0_n      ", mesh.BCp.type );

            for (int p=0; p< input.model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }

            if ( input.model.anisotropy == 1 ) {
                MinMaxArrayTag( mesh.FS_AR_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1),      "FS_AR_n   ", mesh.BCp.type );
                MinMaxArrayTag( mesh.FS_AR_s,  1.0,   (mesh.Nx)*(mesh.Nz),          "FS_AR_s   ", mesh.BCg.type );
                MinMaxArrayTag( mesh.angle_n,  180/M_PI,   (mesh.Nx-1)*(mesh.Nz-1), "angle_n   ", mesh.BCp.type );
                MinMaxArrayTag( mesh.angle_s,  180/M_PI,   (mesh.Nx)*(mesh.Nz),     "angle_s   ", mesh.BCg.type );
                MinMaxArray(particles.nx, 1.0, particles.Nb_part,                   "nx        " );
                MinMaxArray(particles.nz, 1.0, particles.Nb_part,                   "nz        " );
            }
        }

        printf("** Time for particles interpolations I = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );

        if (input.model.mechanical == 1 ) {

            if (input.crazyConductivity) {
                for (int n = 0; n < input.crazyConductivity->nPhases; n++) {
                const int    phase           = input.crazyConductivity->phases[n];
                input.materials.k_eff[phase] = input.materials.k[phase];
                }
                printf("Running with normal conductivity for the asthenosphere...\n");
            }
            // Allocate and initialise solution and RHS vectors
            SetBCs(*setup->SetBCs, &input, &mesh);

            // Reset fields and BC values if needed
            //        if ( input.model.pure_shear_ALE == 1 ) InitialiseSolutionFields( &mesh, &input.model );
            InitialiseSolutionFields( &mesh, &input.model );
            EvalNumberOfEquations( &mesh, &Stokes );
            if ( NeedJacobian ) EvalNumberOfEquations( &mesh, &Jacob  );
            SAlloc( &Stokes,  Stokes.neq );
            if ( NeedJacobian ) SAlloc(  &Jacob,   Jacob.neq );
            SAlloc( &StokesA, Stokes.neq_mom );
            SAlloc( &StokesB, Stokes.neq_mom );
            SAlloc( &StokesC, Stokes.neq_cont);
            SAlloc( &StokesD, Stokes.neq_cont );
            if ( NeedJacobian ) {
              SAlloc( &JacobA, Stokes.neq_mom );
              SAlloc( &JacobB, Stokes.neq_mom );
              SAlloc( &JacobC, Stokes.neq_cont);
              SAlloc( &JacobD, Stokes.neq_cont );
            }
            printf( "Linear systems allocated\n");
            printf( "neq_tot = %d, neq_mom = %d, neq_cont = %d\n", Stokes.neq, Stokes.neq_mom, Stokes.neq_cont );

            // Non-linear iteration cycle
            Nmodel.nit         = 0;
            input.model.nit    = 0;
            Nmodel.stagnated   = 0;
            int nstag          = 0;
            int IsJacobianUsed = 0;

            ArrayEqualArray( mesh.p_start,    mesh.p_in,      (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.u_start,    mesh.u_in,      (mesh.Nx)  *(mesh.Nz+1) );
            ArrayEqualArray( mesh.v_start,    mesh.v_in,      (mesh.Nx+1)*(mesh.Nz)   );

            // Set up solver context
            cholmod_start( &CholmodSolver.c );
            if ( Nmodel.nit==0 ) CholmodSolver.Analyze = 1;
            printf("Run CHOLMOD analysis yes/no: %d \n", CholmodSolver.Analyze);
            int Nmax_picard = Nmodel.nit_max;
            int IsFirstNewtonStep;
            // If Picard 2 Newton is activated: force initial Newton step
            if ( IsNewtonStep == 1 && Nmodel.Picard2Newton == 1 ) {
              input.model.Newton = 0;
              IsNewtonStep       = 0;
              IsFirstNewtonStep  = 1;
            }

            // Set vector x = [u;p]
            InitialiseSolutionVector( &mesh, &Stokes, &input.model );

            //------------------------------------------------------------------------------------------------------------------------------//

            while ( Nmodel.nit <= Nmax_picard && nstag< input.model.max_num_stag) {

                if ( Nmodel.nit > 0 && Nmodel.Picard2Newton == 1 ) {
                    printf("input.model.Newton = %d --- %d\n", input.model.Newton, Nmodel.rp_rel[Nmodel.nit-1] < Nmodel.Picard2Newton_tol);
                    if (Nmodel.rx_rel[Nmodel.nit-1] < Nmodel.Picard2Newton_tol || Nmodel.rz_rel[Nmodel.nit-1] < Nmodel.Picard2Newton_tol || Nmodel.rp_rel[Nmodel.nit-1] < Nmodel.Picard2Newton_tol || Nmodel.nit>= Nmodel.max_Pic_its) {
                        if ( IsFirstNewtonStep == 0 ) CholmodSolver.Analyze = 0;
                        if ( IsFirstNewtonStep == 1 ) {
                            cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                            CholmodSolver.Analyze = 1;
                            IsFirstNewtonStep    = 0;
                        }
                        input.model.Newton = 1;
                        IsNewtonStep = 1;
                    }
                }

                // Determine whether Jacobian matrix should be assembled
                if (input.model.Newton == 1 || input.model.anisotropy == 1) IsJacobianUsed = 1;

                printf("**********************************************\n");
                if ( IsNewtonStep == 1 ) {
                    printf("*** Newton it. %02d of %02d (step = %05d) on %03d threads ***\n", Nmodel.nit, Nmodel.nit_max, input.model.step, num_threads);
                    Nmodel.LogIsNewtonStep[Nmodel.nit] = 1;
                }
                else {
                    printf("*** Picard it. %02d of %02d (step = %05d) on %03d threads ***\n", Nmodel.nit, Nmodel.nit_max, input.model.step, num_threads);
                    Nmodel.LogIsNewtonStep[Nmodel.nit] = 0;
                }
                printf("**********************************************\n");

                // Update non-linear rheology
                UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, input.materials, &input.model, &Nmodel, input.scaling, 0, 0 );

                if (input.model.noisy == 1 ) {
                    MinMaxArrayTag( mesh.T, input.scaling.T, (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCt.type );
                    MinMaxArrayTag( mesh.p0_n, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.p_in, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.div_u, input.scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "div       ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exxd, input.scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "exxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.ezzd, input.scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "ezzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exz, input.scaling.E, (mesh.Nx-0)*(mesh.Nz-0), "exz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.sxxd, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "sxxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.szzd, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "szzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.sxz, input.scaling.S, (mesh.Nx-0)*(mesh.Nz-0), "sxz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_s, input.scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_n, input.scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.eta_phys_s, input.scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_phys_s", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_phys_n, input.scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho_s, input.scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.rho_n, input.scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho0_n, input.scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.d_n, input.scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.X_s, 1.0, (mesh.Nx)*(mesh.Nz),     "X_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.X_n, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "X_n     ", mesh.BCp.type );
                    if (input.model.anisotropy==1) MinMaxArrayTag( mesh.aniso_factor_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "ani_fac_n ",   mesh.BCp.type );
                    if (input.model.anisotropy==1) MinMaxArrayTag( mesh.aniso_factor_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "ani_fac_s ",   mesh.BCg.type );
                }

                if (input.model.writer_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, "Output_BeforeSolve", input.materials, input.scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, "Particles_BeforeSolve", input.materials, input.scaling );
                }

                // Build discrete system of equations - Linearised Picard
                BuildStokesOperatorDecoupled  ( &mesh, input.model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1 );

                // Build discrete system of equations - Jacobian (do it also for densification since drhodp is needed)
                if ( (IsFullNewton   == 1 && Nmodel.nit > 0) || input.model.density_variations == 1  ) RheologicalOperators( &mesh, &input.model, &input.materials, &input.scaling, 1 );
                if ( IsJacobianUsed == 1 ) BuildJacobianOperatorDecoupled( &mesh, input.model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in,  &Jacob,  &JacobA,  &JacobB,  &JacobC,   &JacobD, 1 );
                
                // Diagonal input.scaling
                if ( input.model.diag_scaling ) {
                    if (input.model.Newton          == 0 )  ExtractDiagonalScale( &StokesA, &StokesB, &StokesC, &StokesD );
                    if (input.model.Newton          == 1 )  ExtractDiagonalScale( &JacobA,  &JacobB,  &JacobC,   &JacobD );
                    if (input.model.Newton          == 1 )  ArrayEqualArray(StokesA.d, JacobA.d, StokesA.neq);
                    if (input.model.Newton          == 1 )  ArrayEqualArray(StokesC.d, JacobC.d, StokesC.neq);
                    ScaleMatrix( &StokesA, &StokesB, &StokesC, &StokesD );
                }

                // Needs to be done after matrix assembly since diagonal input.scaling is used in there
                printf("---- Non-linear residual ----\n");
                RheologicalOperators( &mesh, &input.model, &input.materials, &input.scaling, 0 );
                EvaluateStokesResidualDecoupled( &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &Nmodel, &mesh, input.model, input.scaling, 0 );
                printf("---- Non-linear residual ----\n");

                if (input.model.Newton == 1 && input.model.anisotropy == 0 ) {
                    ArrayEqualArray( JacobA.F, StokesA.F,  StokesA.neq );
                    ArrayEqualArray( JacobC.F, StokesC.F,  StokesC.neq );
                }
                if (input.model.Newton == 1 && input.model.diag_scaling ) ScaleMatrix( &JacobA,  &JacobB,  &JacobC,  &JacobD  );

                // Store residuals
                Nmodel.resx_f = Nmodel.resx; Nmodel.rx_abs[Nmodel.nit] = Nmodel.resx; Nmodel.rx_rel[Nmodel.nit] = Nmodel.resx/Nmodel.resx0;
                Nmodel.resz_f = Nmodel.resz; Nmodel.rz_abs[Nmodel.nit] = Nmodel.resz; Nmodel.rz_rel[Nmodel.nit] = Nmodel.resz/Nmodel.resz0;
                Nmodel.resp_f = Nmodel.resp; Nmodel.rp_abs[Nmodel.nit] = Nmodel.resp; Nmodel.rp_rel[Nmodel.nit] = Nmodel.resp/Nmodel.resp0;

                if (input.model.writer_debug == 1 ) WriteResiduals( mesh, input.model, Nmodel, input.scaling );

                // if pass --> clear matrix break
                if ( (Nmodel.resx < Nmodel.nonlin_abs_mom || Nmodel.resx/Nmodel.resx0 < Nmodel.nonlin_rel_mom) && (Nmodel.resz < Nmodel.nonlin_abs_mom || Nmodel.resz/Nmodel.resz0 < Nmodel.nonlin_rel_mom) && (Nmodel.resp < Nmodel.nonlin_abs_div || Nmodel.resp/Nmodel.resp0 < Nmodel.nonlin_rel_div) ) {
                    printf( "Non-linear solver converged to nonlin_abs_mom = %2.2e nonlin_abs_div = %2.2e nonlin_rel_mom = %2.2e nonlin_rel_div = %2.2e\n", Nmodel.nonlin_abs_mom, Nmodel.nonlin_abs_div, Nmodel.nonlin_rel_mom, Nmodel.nonlin_rel_div );
                    FreeSparseSystems( IsJacobianUsed, 1, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                    break;
                }

                // Direct solve
                t_omp = (double)omp_get_wtime();
                if ( IsJacobianUsed==1 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &input.model, &particles, &topo_chain, &topo, input.materials, input.scaling,  &JacobA,  &JacobB,  &JacobC );
                else                     SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &input.model, &particles, &topo_chain, &topo, input.materials, input.scaling, &StokesA, &StokesB, &StokesC );

                if ( Nmodel.stagnated == 1 && input.model.safe_mode <= 0 ) {
                    printf( "Non-linear solver stagnated to res_u = %2.2e res_z = %2.2e\n", Nmodel.resx_f, Nmodel.resz_f );
                    printf( "You may want to try setting line_search_min > 0.0\n Good luck good man!\n");
                    FreeSparseSystems( IsJacobianUsed, 1, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                    break;
                }
                
                // Check for stagnation and apply safety restrictions is requested and required
                if ( Nmodel.stagnated == 1 && input.model.elastic == 1 && input.model.safe_mode == 1 ) {
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (nonlin_abs_mom = %2.2e nonlin_abs_div = %2.2e)\e[m\n", Nmodel.nonlin_abs_mom, Nmodel.nonlin_abs_div );
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (nonlin_rel_mom = %2.2e nonlin_rel_div = %2.2e)\e[m\n", Nmodel.nonlin_rel_mom, Nmodel.nonlin_rel_div );
                    printf( "\e[1;31mReducing the timestep, and restart the iterations cycle...\e[m\n");
                    printf( "Before reduction: input.model.dt =, %2.2e\n", input.model.dt* input.scaling.t);
                    // ----------------------
                    input.model.dt /= input.model.safe_dt_div;
                    printf( "Timestep divided by %2.2f => NEW CURRENT input.model.dt =, %2.2e\n", input.model.safe_dt_div, input.model.dt* input.scaling.t);
                    Nmodel.stagnated = 0;
                    // ----------------------
                    Nmodel.nit = -1;
                    printf( "Restart solutions\n");
                    ArrayEqualArray( mesh.p_in, mesh.p_start,  (mesh.Nx-1)*(mesh.Nz-1) );
                    ArrayEqualArray( mesh.u_in, mesh.u_start,  (mesh.Nx)  *(mesh.Nz+1) );
                    ArrayEqualArray( mesh.v_in, mesh.v_start,  (mesh.Nx+1)*(mesh.Nz)   );
                    nstag++;
                    //-----------------------
                    printf( "nstag value = %02d - max_num_stag = %02d\n", nstag, input.model.max_num_stag);
                    if (nstag== input.model.max_num_stag) {
                        printf( "CheckDoudzOut!!\n");
                        exit(0);
                    //-----------------------
                    }
                }
                FreeSparseSystems( IsJacobianUsed, 1, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                Nmodel.nit++;
                input.model.nit = Nmodel.nit;
            }

            // Clean solver context
            if ( Nmodel.nit>0 ) {
                printf("Cleaning up Cholesky factors --- nit = %02d\n",  Nmodel.nit);
                cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                cholmod_finish( &CholmodSolver.c );
            }

            // Update rheology
            UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, input.materials, &input.model, &Nmodel, input.scaling, 0, 1 );
            ArrayEqualArray( mesh.p_in, mesh.p_corr, (mesh.Nx-1)*(mesh.Nz-1) ); // make sure corrected pressure is used

            // Min/Max velocities
            MinMaxArray( mesh.u_in, input.scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx grid" );
            MinMaxArray( mesh.v_in, input.scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz grid" );
            MinMaxArray( mesh.p_in, input.scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "      P" );
            MinMaxArray( mesh.div_u, input.scaling.E, (mesh.Nx-1)*(mesh.Nz-1), " div(V)" );

            printf("--------------------------------------------------------------\n");
            int i, nit;
            if (Nmodel.nit>=Nmodel.nit_max)  nit = Nmodel.nit_max;
            if (Nmodel.nit< Nmodel.nit_max)  nit = Nmodel.nit;
            if (Nmodel.Picard2Newton == 1 )  printf("Picard 2 Newton is activated with condition: %2.2e\n", Nmodel.Picard2Newton_tol);
            for (i=0; i<=nit; i++) {
                if (Nmodel.LogIsNewtonStep[i] == 1) printf("New. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, Nmodel.rx_abs[i], Nmodel.rz_abs[i], Nmodel.rp_abs[i], Nmodel.rx_rel[i], Nmodel.rz_rel[i], Nmodel.rp_rel[i]);
                else                                printf("Pic. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, Nmodel.rx_abs[i], Nmodel.rz_abs[i], Nmodel.rp_abs[i], Nmodel.rx_rel[i], Nmodel.rz_rel[i], Nmodel.rp_rel[i]);
                if (i == Nmodel.nit_max && input.model.safe_mode == 1) {
                    printf("Exit: Max iteration reached: Nmodel.nit_max = %02d! Check what you wanna do now...\n",Nmodel.nit_max);
                    if ( (Nmodel.resx < Nmodel.nonlin_abs_mom) && (Nmodel.resz < Nmodel.nonlin_abs_mom) && (Nmodel.resp < Nmodel.nonlin_abs_div) ) {}
                    else exit(1);
                }
            }
            printf("--------------------------------------------------------------\n");

            // plot residuals
            if (input.model.gnuplot_log_res == 1 ) { // && input.model.step % writer_step == 0

                printf("DOING GNU PLOTTING!\n");

                int NumCommands = 3;
                char *GNUplotCommands[] = {"set title \"Non-linear residuals\"", "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5", "set pointintervalbox 3"};
                for (i=0; i<NumCommands; i++) fprintf(GNUplotPipe, "%s \n", GNUplotCommands[i]); //Send commands to gnuplot one by one.

                fprintf(GNUplotPipe, "plot '-' with linespoints linestyle 1\n");
                for (i=0; i< Nmodel.nit+1; i++) {
                    fprintf(GNUplotPipe, "%lf %lf \n", (double)i, log10(Nmodel.rx_rel[i])); //Write the data to a temporary file
                }
                fprintf(GNUplotPipe, "e\n");
                fflush(GNUplotPipe);
            }

        }

        //------------------------------------------------------------------------------------------------------------------------------//

        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE --- Free surface - interpolate velocity components on the free surface
        if (input.model.free_surface == 1 ) {
            SurfaceVelocity( &mesh, input.model, &topo, &topo_chain, input.scaling );
            MinMaxArray( topo_chain.Vx, input.scaling.V, topo_chain.Nb_part,       "Vx surf." );
            MinMaxArray( topo_chain.Vz, input.scaling.V, topo_chain.Nb_part,       "Vz surf." );
            MinMaxArray( topo_chain_ini.Vx, input.scaling.V, topo_chain_ini.Nb_part,   "Vx surf. ini." );
            MinMaxArray( topo_chain_ini.Vz, input.scaling.V, topo_chain_ini.Nb_part,   "Vz surf. ini." );
        }
        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE

        // Update stresses on markers
        UpdateParticleStress(  &mesh, &particles, &input.model, &input.materials, &input.scaling );

        // Update pressure on markers
        UpdateParticlePressure( &mesh, input.scaling, input.model, &particles, &input.materials );

        // Grain size evolution
        UpdateParticleGrainSize( &mesh, input.scaling, input.model, &particles, &input.materials );

        if (input.model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n   , input.scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d0", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n    , input.scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d ", mesh.BCp.type );
            MinMaxArrayPart( particles.d, input.scaling.L, particles.Nb_part, "d on markers", particles.phase ) ;
        }

        // Update phi on the particles
        UpdateParticlePhi( &mesh, input.scaling, input.model, &particles, &input.materials );

        //------------------------------------------------------------------------------------------------------------------------------//

        // Apply BC's prior to thermal solve
        // SetBCs(*setup->SetBCs, &input, &mesh);

        if (input.model.thermal == 1 ) {

            printf("*************************************\n");
            printf("*********** Thermal solver **********\n");
            printf("*************************************\n");
            
            t_omp = (double)omp_get_wtime();

            // Matrix assembly and direct solve
            EnergyDirectSolve( &mesh, input.model,  mesh.rhs_t, &particles, input.model.dt, input.model.shear_heating, input.model.adiab_heating, input.scaling, 1 );
            MinMaxArray(particles.T, input.scaling.T, particles.Nb_part, "T part. before UpdateParticleEnergy");
            printf("** Time for Thermal solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        }

        // Update energy on particles
        UpdateParticleEnergy( &mesh, input.scaling, input.model, &particles, &input.materials );

        //--------------------------------------------------------------------------------------------------------------------------------//

        if (input.model.chemical_diffusion == 1)  {
            if (input.model.unsplit_diff_reac == 0 ) ArrayEqualArray( mesh.X_n, mesh.X0_n,  (mesh.Nx-1)*(mesh.Nz-1) );
            printf("*************************************\n");
            printf("********** Chemical solver **********\n");
            printf("*************************************\n");
            P2Mastah ( &input.model, particles, input.materials.k_chem, &mesh, mesh.kc_x, mesh.BCu.type,  0, 0, interp, vxnodes, input.model.interp_stencil);
            P2Mastah ( &input.model, particles, input.materials.k_chem, &mesh, mesh.kc_z, mesh.BCv.type,  0, 0, interp, vznodes, input.model.interp_stencil);
            ChemicalDirectSolve( &mesh, input.model, &particles, &input.materials, input.model.dt, input.scaling );
        }
        else {
            if (input.model.density_variations==0) ArrayEqualArray( mesh.X_n, mesh.X0_n,  (mesh.Nx-1)*(mesh.Nz-1) );
        }
        UpdateParticleX( &mesh, input.scaling, input.model, &particles, &input.materials );

        // ArrayEqualArray( mesh.X_n, mesh.X0_n,  (mesh.Nx-1)*(mesh.Nz-1) );
        // UpdateParticleXpips( &mesh, input.scaling, input.model, &particles, &input.materials );


        //--------------------------------------------------------------------------------------------------------------------------------//

        // Update maximum pressure and temperature on markers
        if (input.model.track_T_P_x_z == 1) UpdateMaxPT(input.scaling, input.model, &particles );

        ComputeMeanQuantitesForTimeSeries( &mesh );
        LogTimeSeries( &mesh, input.model, input.scaling );

        //------------------------------------------------------------------------------------------------------------------------------//

        if (input.model.advection == 1 ) {

            printf("*************************************\n");
            printf("************** Advection ************\n");
            printf("*************************************\n");

            t_omp = (double)omp_get_wtime();
            double dt_solve = input.model.dt;
            int    nsub=1.0, isub;
            double dt_sub;
            EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &input.model, input.scaling, &mesh, 0 );

            if (input.model.dt < 0.95*dt_solve ) { // if dt advection is lower than dt used for thermo-mechanical solve then split dt
                nsub   = ceil(dt_solve/ input.model.dt);
                dt_sub = dt_solve / nsub;
                printf("dt advection = %2.9e --- dt_solve = %2.9e\n", input.model.dt* input.scaling.t, dt_solve* input.scaling.t );
                printf("dt is %lf larger than dt_solve: need %d substeps of %2.2e s\n", dt_solve/ input.model.dt, nsub , dt_sub* input.scaling.t );
                printf("So: nsub*dt_sub = %2.2e for dt_solve = %2.2e\n", nsub*dt_sub* input.scaling.t, dt_solve* input.scaling.t);
                input.model.dt = dt_sub;
            }
            else {
                nsub     = 1; // if dt advection is larger than dt used for thermo-mechanical solve then, forget it, and set dt to dt_solve
                input.model.dt = dt_solve;
            }

            // Loop on substeps
            for (isub=0; isub<nsub; isub++) {

                printf("************** Advection step %03d of %03d: dtsub = %2.2e **************\n", isub, nsub, input.model.dt* input.scaling.t );

                // Advect domain boundaries
                if (input.model.pure_shear_ALE > 0 ) {
                    PureShearALE( &input.model, &mesh, &topo_chain, input.scaling );
                }

                if (input.model.free_surface == 1 ) {
                    // Advect free surface with RK2
                    RogerGuntherII( &topo_chain, input.model, mesh, 1, input.scaling );
                    RogerGuntherII( &topo_chain_ini, input.model, mesh, 1, input.scaling );
                }

                // Correction for particle inflow 0
                if (input.model.pure_shear_ALE == -1 && input.model.periodic_x == 0) ParticleInflowCheck( &particles, &mesh, input.model, topo, 0 );

                // Advect fluid particles
                RogerGuntherII( &particles, input.model, mesh, 1, input.scaling );

                // Correction for particle inflow 1
                if (input.model.pure_shear_ALE == -1 && input.model.periodic_x == 0) ParticleInflowCheck( &particles, &mesh, input.model, topo, 1 );

                // Update accumulated strain
                AccumulatedStrainII( &mesh, input.scaling, input.model, &particles,  mesh.xc_coord,  mesh.zc_coord, mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );

                // Update deformation gradient tensor components
                if (input.model.finite_strain == 1 ) DeformationGradient( mesh, input.scaling, input.model, &particles );

                if (input.model.writer_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, "Output_BeforeSurfRemesh", input.materials, input.scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, "Particles_BeforeSurfRemesh", input.materials, input.scaling );
                }

                if (input.model.free_surface == 1 ) {

                    // Get current topography
                    ProjectTopography(  &topo,     &topo_chain,     input.model, mesh, input.scaling, mesh.xg_coord, 0 );
                    ProjectTopography(  &topo_ini, &topo_chain_ini, input.model, mesh, input.scaling, mesh.xg_coord, 0 );
                    MarkerChainPolyFit( &topo,     &topo_chain,     input.model, mesh );
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, input.model, mesh );

                    // Remesh free surface I
                    RemeshMarkerChain( &topo_chain,     &topo, input.model, input.scaling, &mesh, 1 );
                    RemeshMarkerChain( &topo_chain_ini, &topo_ini, input.model, input.scaling, &mesh, 1 );

                    // Project topography on vertices
                    ProjectTopography( &topo,     &topo_chain, input.model, mesh, input.scaling, mesh.xg_coord, 0 );
                    ProjectTopography( &topo_ini, &topo_chain_ini, input.model, mesh, input.scaling, mesh.xg_coord, 0 );

                    if ( input.model.writer_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, "Output_AfterSurfRemesh", input.materials, input.scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, "Particles_AfterSurfRemesh", input.materials, input.scaling );
                    }

                    // Diffuse topography
                    if ( input.model.surface_processes >= 1 )  DiffuseAlongTopography( &mesh, input.model, input.scaling, topo.height, topo.height, mesh.Nx, 0.0, input.model.dt );

                    // Marker chain polynomial fit
                    MarkerChainPolyFit( &topo,     &topo_chain, input.model, mesh );
                    CorrectTopoIni( &particles, input.materials, &topo_chain_ini, &topo, input.model, input.scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, input.model, mesh );

                    // Sedimentation
                    if ( input.model.surface_processes == 2 ) {
                        AddPartSed( &particles, input.materials, &topo_chain, &topo, input.model, input.scaling, &mesh);
                        CountPartCell( &particles, &mesh, input.model, topo, topo_ini, 0, input.scaling );
                    }

                    if ( input.model.writer_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, "Outputx", input.materials, input.scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, "Particlesx", input.materials, input.scaling );
                    }

                    // Remesh free surface II
                    RemeshMarkerChain( &topo_chain,     &topo, input.model, input.scaling, &mesh, 2 );
                    RemeshMarkerChain( &topo_chain_ini, &topo_ini, input.model, input.scaling, &mesh, 2 );
                    CorrectTopoIni( &particles, input.materials, &topo_chain_ini, &topo, input.model, input.scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, input.model, mesh );

                    // Remove particles that are above the surface
                    CleanUpSurfaceParticles( &particles, &mesh, topo, input.scaling );

                    // Call cell flagging routine for free surface calculations
                    CellFlagging( &mesh, input.model, topo, input.scaling );
                }

                printf("** Time for advection solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

#ifdef _HDF5_
                if ( input.model.writer_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, "Outputxx", input.materials, input.scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, "Particlesxx", input.materials, input.scaling );
                }
#endif

                //                printf("*************************************\n");
                //                printf("************** Reseeding ************\n");
                //                printf("*************************************\n");

                // Count the number of particle per cell
                t_omp = (double)omp_get_wtime();
                int NumPartOld  = particles.Nb_part;
                if ( input.model.reseed_markers == 1 ) CountPartCell( &particles, &mesh, input.model, topo, topo_ini, 1, input.scaling );
                CountPartCell( &particles, &mesh, input.model, topo, topo_ini, 0, input.scaling );

                printf(GREEN "After re-seeding :\n" RESET);
                printf(GREEN "Initial number of particles = %d\n" RESET, particles.Nb_part_ini);
                printf(GREEN "Old number of particles     = %d\n" RESET, NumPartOld           );
                printf(GREEN "New number of particles     = %d\n" RESET, particles.Nb_part    );

                printf("** Time for CountPartCell = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

                // Remove particles that would be above the surface
                if (input.model.free_surface == 1 ) {
                    CleanUpSurfaceParticles( &particles, &mesh, topo, input.scaling );
                    CellFlagging( &mesh, input.model, topo, input.scaling );
                }
            }
            input.model.dt = dt_solve;
        }

        // Update time
        // input.model.time += input.model.dt;

        // Free solution arrays
        SFree( &Stokes );
        if ( NeedJacobian ) SFree( &Jacob );
        SFree( &StokesA );
        SFree( &StokesB );
        SFree( &StokesC );
        SFree( &StokesD );
        if ( NeedJacobian ) {
            SFree( &JacobA );
            SFree( &JacobB );
            SFree( &JacobC );
            SFree( &JacobD );
        }
        DoodzFree( Stokes.eqn_u );
        DoodzFree( Stokes.eqn_v );
        DoodzFree( Stokes.eqn_p );
        if ( NeedJacobian ) {
            DoodzFree( Jacob.eqn_u );
            DoodzFree( Jacob.eqn_v );
            DoodzFree( Jacob.eqn_p );
        }

        //------------------------------------------------------------------------------------------------------------------------------//

        // Write output data
        if ( input.model.writer == 1 && input.model.step % input.model.writer_step == 0 ) {

            printf("*************************************\n");
            printf("********* Write output files ********\n");
            printf("*************************************\n");

            // Breakpoint file
            t_omp = (double)omp_get_wtime();
            if (input.model.delete_breakpoints == 1 ) DeletePreviousBreakpoint(input.model.step, input.model.writer_step  );
            MakeBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, input.model, &topo, &topo_ini, input.scaling );
            // UpdateInputFile( inputFileName, input.model.step);
            printf("** Time for Breakpoint file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

            // Visualisation file
#ifndef _VG_
            t_omp = (double)omp_get_wtime();
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, input.model, Nmodel, BaseOutputFileName, input.materials, input.scaling );
            if (input.model.writer_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, input.model, BaseParticleFileName, input.materials, input.scaling );
            printf("** Time for Output file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
#endif
        }

        printf("** Total timestep calculation time = %lf sec\n", (double)((double)omp_get_wtime() - t_omp_step) );
        printf("** Model time = %2.2e sec\n", input.model.time* input.scaling.t );
        printf("** Current dt = %2.2e sec, Old dt = %2.2e sec\n", input.model.dt* input.scaling.t, input.model.dt0* input.scaling.t );

    }

    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                           END TIME LOOP : c'est fini les Doud'series...
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//

    if (input.model.gnuplot_log_res == 1) pclose(GNUplotPipe);

    //    // Free markers chains
    if (input.model.free_surface == 1 ) FreeMarkerChain( &topo,     &topo_chain     );
    if (input.model.free_surface == 1 ) FreeMarkerChain( &topo_ini, &topo_chain_ini );

    // Free Phase diagrams if activated
    if (input.model.isPD == 1 ) FreePhaseDiagrams( &input.model );
    //
    // Free particle arrays
    PartFree( &particles, &input.model );

    // Free arrays
    GridFree( &mesh, &input.model );

    if (input.crazyConductivity) {
        free(input.crazyConductivity->phases);
        free(input.crazyConductivity);
    }
    // Free char*'s
    free(BaseOutputFileName);
    free(BaseParticleFileName);
    free(input.model.import_file);
    free(input.model.import_files_dir);
    free(input.model.writer_subfolder);
    free(input.model.initial_markers_file);
    // free(input.model.description);

    NmodelFree(Nmodel);

    // just for the quartz-coesite case
    for ( int k=0; k< input.model.Nb_phases; k++) {

        if (input.materials.density_model[k] == 4 ) {
            printf("Unloading 1D diagram for coesite quartz only baby...\n");
            DoodzFree(input.model.PD1Drho[k]);
        }
    }
    DoodzFree(input.model.PD1DnP);
    DoodzFree(input.model.PD1Drho);
    DoodzFree(input.model.PD1Dmin);
    DoodzFree(input.model.PD1Dmax);
    if (input.model.kinetics==1 ) DoodzFree(input.model.kin_dG);

    printf("\n********************************************************\n");
    printf("************* Ending MDOODZ 7.0 simulation *************\n");
    printf("********************************************************\n");
}
