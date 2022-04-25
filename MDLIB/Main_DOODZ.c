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

void RunMDOODZ(MdoodzInstance *this) {
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

    fin_name = this->inputFileName;
    asprintf(&fin_name,"%s", this->inputFileName);

    printf("\n********************************************************\n");
    printf("************ Starting MDOODZ 7.0 simulation ************\n");
    printf("********************************************************\n");

    // Read input data
    ReadInputFile( fin_name, &istep, &irestart, &writer, &writer_step, &this->model, &this->scaling, &this->materials, &particles, &Nmodel);
    ValidateSetup(this);
    this->model.L0 = this->model.xmax - this->model.xmin; // Save initial length

    printf("*************************************\n");
    printf("****** Allocate and initialise ******\n");
    printf("*************************************\n");

    // Multi resolution allocation pointers
    GridAlloc( &mesh, &this->model );

    // Initialise grid coordinates
    SetGridCoordinates( &mesh, &this->model, this->model.Nx, this->model.Nz );

    // Get grid indices
    GridIndices( &mesh );

    // Initialise data for logs of iterations
    Nmodel.nit             = 0;
    Nmodel.rx_abs          = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); Nmodel.rx_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    Nmodel.rz_abs          = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); Nmodel.rz_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    Nmodel.rp_abs          = DoodzCalloc(Nmodel.nit_max+1, sizeof(double)); Nmodel.rp_rel = DoodzCalloc(Nmodel.nit_max+1, sizeof(double));
    Nmodel.LogIsNewtonStep = DoodzCalloc(Nmodel.nit_max+1, sizeof(int));

    Nx = mesh.Nx; Nz = mesh.Nz; Ncx = Nx-1; Ncz = Nz-1;
    if ( this->model.aniso  == 1 ) this->model.Newton = 1;
    if ( this->model.Newton == 1 ) IsNewtonStep = 1;
    if ( this->model.Newton == 0 ) IsNewtonStep = 0;

    printf("*************************************\n");
    printf("******* Initialize particles ********\n");
    printf("*************************************\n");

    // Allocate particle fields
    PartAlloc( &particles, &this->model );

    // Allocate marker chain
    if ( this->model.free_surf == 1 ) AllocateMarkerChain( &topo,     &topo_chain,     this->model );
    if ( this->model.free_surf == 1 ) AllocateMarkerChain( &topo_ini, &topo_chain_ini, this->model );

    // Set new particle distribution
    if ( irestart == 0 ) {

        this->model.step = 0;
        this->model.time = 0.0;

        // Initialise particle fields
        PartInit( &particles, &this->model );

        if (this->crazyConductivity) {
          for (int n = 0; n < this->crazyConductivity->nPhases; n++) {
            const int    phase           = this->crazyConductivity->phases[n];
            const double k_crazy         = this->crazyConductivity->multiplier * this->materials.k[phase];
            this->materials.k_eff[phase] = k_crazy;
          }
          printf("Running with crazy conductivity for the asthenosphere!!\n");
        }

        // Initial grid tags
        SetBCs( this, &mesh);
        if ( this->model.free_surf == 1 ) {

            // Define the horizontal position of the surface marker chain
            SetTopoChainHorizontalCoords( &topo,     &topo_chain,     this->model, mesh, this->scaling );
            SetTopoChainHorizontalCoords( &topo_ini, &topo_chain_ini, this->model, mesh, this->scaling );

            // Define the vertical position of the surface marker chain
            BuildInitialTopography(this, &topo_chain );
            BuildInitialTopography(this, &topo_chain_ini );

            // Project topography on vertices
            ProjectTopography( &topo, &topo_chain, this->model, mesh, this->scaling, mesh.xg_coord, 0 );

            // Marker chain polynomial fit
            MarkerChainPolyFit( &topo, &topo_chain, this->model, mesh );

            // Call cell flagging routine for free surface calculations
            CellFlagging( &mesh, this->model, topo, this->scaling );
        }
        // Set particles coordinates
        PutPartInBox( &particles, &mesh, this->model, topo, this->scaling );

        // Set phases on particles
        SetParticles( this, &particles);

        if ( this->model.free_surf == 1 ) CleanUpSurfaceParticles( &particles, &mesh, topo, this->scaling );

        // Create phase percentage arrays
        if (this->model.cpc==-1) CountPartCell_BEN( &particles, &mesh, this->model, topo, 0, this->scaling );
        if (this->model.cpc== 0) CountPartCell_Old( &particles, &mesh, this->model, topo, 0, this->scaling  );
        if (this->model.cpc== 1) CountPartCell    ( &particles, &mesh, this->model, topo, topo_ini, 0, this->scaling  );
        if (this->model.cpc== 2) CountPartCell2    ( &particles, &mesh, this->model, topo, topo_ini, 0, this->scaling  );

        P2Mastah( &this->model, particles, this->materials.eta0, &mesh, mesh.eta_s, mesh.BCg.type,  0, 0, interp, vert, this->model.itp_stencil);
        P2Mastah( &this->model, particles, this->materials.eta0, &mesh, mesh.eta_n, mesh.BCp.type,  0, 0, interp, cent, this->model.itp_stencil);

        P2Mastah( &this->model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, this->model.itp_stencil);
        P2Mastah( &this->model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);

        P2Mastah( &this->model, particles, particles.T,  &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &this->model, particles, this->materials.Cv, &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);

        P2Mastah( &this->model, particles, this->materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah( &this->model, particles, this->materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);

        // UpdateDensity( &mesh, &particles, &this->materials, &this->model, &this->scaling );

        // if ( this->model.eqn_state > 0) {
        //     P2Mastah( &this->model, particles, particles.rho, &mesh, mesh.rho_s, mesh.BCg.type,  1, 0, interp, vert, this->model.itp_stencil);
        //     P2Mastah( &this->model, particles, particles.rho, &mesh, mesh.rho_n, mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
        // }
        // else {
        //     P2Mastah( &this->model, particles, this->materials.rho, &mesh, mesh.rho_s, mesh.BCg.type,  0, 0, interp, vert, this->model.itp_stencil);
        //     P2Mastah( &this->model, particles, this->materials.rho, &mesh, mesh.rho_n, mesh.BCp.type,  0, 0, interp, cent, this->model.itp_stencil);
        // }

        // MinMaxArrayTag( mesh.rho_s,      this->scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
        // MinMaxArrayTag( mesh.rho_n,      this->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
        // exit(1);
        if ( this->model.noisy == 1 ) {
            MinMaxArray( mesh.u_in,  this->scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in,  this->scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in,  this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
        }
        // Initial solution fields (Fine mesh)
        SetBCs( this, &mesh);
        InitialiseSolutionFields( &mesh, &this->model );

        MinMaxArray( mesh.u_in,  this->scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
        MinMaxArray( mesh.v_in,  this->scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
        MinMaxArray( mesh.p_in,  this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );

        printf("*************************************\n");
        printf("****** Initialize temperature *******\n");
        printf("*************************************\n");

        // Get energy and related material parameters from particles
        P2Mastah( &this->model, particles, this->materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah( &this->model, particles, this->materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);
        P2Mastah( &this->model, particles, particles.T,     &mesh, mesh.T , mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &this->model, particles, this->materials.Cv,    &mesh, mesh.Cv, mesh.BCp.type,  0, 0, interp, cent, 1);
        P2Mastah( &this->model, particles, this->materials.Qr,    &mesh, mesh.Qr, mesh.BCp.type,  0, 0, interp, cent, 1);

        SetBCs( this, &mesh);
        if ( this->model.thermal_eq == 1 ) ThermalSteps( &mesh, this->model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, this->model.cooling_time, this->scaling );
        if ( this->model.therm_pert == 1 ) SetThermalPert( &mesh, this->model, this->scaling );
        //            Interp_Grid2P_centroids ( particles, particles.T,    &mesh, mesh.T, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &this->model );

        Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &this->model );
        ArrayEqualArray( mesh.T0_n, mesh.T, (mesh.Nx-1)*(mesh.Nz-1) );

        //--------------------------------------------------------------------------------------------------------

        printf("*************************************\n");
        printf("******** Initialize pressure ********\n");
        printf("*************************************\n");

        // Compute dev. strain rate tensor components
        StrainRateComponents( &mesh, this->scaling, &this->model );

        for (int iter=0; iter<10; iter++) {

            // if ( this->model.eqn_state > 0 ) {
                UpdateDensity( &mesh, &particles, &this->materials, &this->model, &this->scaling );
            // }
            // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &this->model );
            // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

            // Free surface - subgrid density correction
            if ( this->model.free_surf == 1 ) {
                SurfaceDensityCorrection( &mesh, this->model, topo, this->scaling  );
            }
            // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

            // Lithostatic pressure for initial visco-plastic viscosity field
            ComputeLithostaticPressure( &mesh, &this->model, this->materials.rho[0], this->scaling, 1 );
            Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_lith, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &this->model );
            P2Mastah( &this->model, particles, particles.P,     &mesh, mesh.p_in , mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
            ArrayEqualArray( mesh.p_in, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.p0_n, mesh.p_lith,  (mesh.Nx-1)*(mesh.Nz-1) );

            MinMaxArrayTag( mesh.p_in,       this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P initial ", mesh.BCp.type );
        }

        InterpCentroidsToVerticesDouble( mesh.p0_n, mesh.p0_s, &mesh, &this->model );
        //            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &this->model );
        //            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &this->model );


        printf("*************************************\n");
        printf("******* Initialize grain size *******\n");
        printf("*************************************\n");

        // Grain size
        InitialiseGrainSizeParticles( &particles, &this->materials );
        P2Mastah( &this->model, particles, particles.d,     &mesh, mesh.d_n , mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
        ArrayEqualArray( mesh.d0_n, mesh.d_n,  (mesh.Nx-1)*(mesh.Nz-1) );

        printf("*************************************\n");
        printf("******** Initialize density *********\n");
        printf("*************************************\n");

        // if ( this->model.eqn_state > 0 ) {
        UpdateDensity( &mesh, &particles, &this->materials, &this->model, &this->scaling );
        // }
        // Interp_Grid2P_centroids( particles, particles.rho, &mesh, mesh.rho_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &this->model );
        // ArrayEqualArray( mesh.rho0_n, mesh.rho_n, (mesh.Nx-1)*(mesh.Nz-1) );

        printf("*************************************\n");
        printf("****** Initialize composition *******\n");
        printf("*************************************\n");

        //            if ( this->model.diffuse_X == 1 ) {
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Diffuse_X(&mesh, &this->model, &this->scaling);
        //                Interp_Grid2P( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Interp_P2N ( particles, particles.X, &mesh, mesh.Xreac_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &this->model );
        //            }

        P2Mastah( &this->model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
        P2Mastah( &this->model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, this->model.itp_stencil);

        if ( this->model.aniso == 1 ) {
            InitialiseDirectorVector ( &mesh, &particles, &this->model, &this->materials );
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d1_n , mesh.BCp.type, -1, 0, interp, cent, this->model.itp_stencil);
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d2_n , mesh.BCp.type, -2, 0, interp, cent, this->model.itp_stencil);
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d1_s , mesh.BCg.type, -1, 0, interp, vert, this->model.itp_stencil);
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d2_s , mesh.BCg.type, -2, 0, interp, vert, this->model.itp_stencil);
            FiniteStrainAspectRatio ( &mesh, this->scaling, this->model, &particles );
            P2Mastah( &this->model, particles, this->materials.aniso_factor,     &mesh, mesh.aniso_factor_n , mesh.BCp.type,  0, 0, interp, cent, this->model.itp_stencil);
            P2Mastah( &this->model, particles, this->materials.aniso_factor,     &mesh, mesh.aniso_factor_s , mesh.BCg.type,  0, 0, interp, vert, this->model.itp_stencil);
        }

        printf("*************************************\n");
        printf("******* Initialize viscosity ********\n");
        printf("*************************************\n");

        // Compute shear modulus, expansivity, compressibiulity, cohesion and friction angle on the grid
        if ( this->model.iselastic == 1 ) ShearModCompExpGrid( &mesh, this->materials, this->model, this->scaling );
        CohesionFrictionDilationGrid( &mesh, &particles, this->materials, this->model, this->scaling );
        ShearModCompExpGrid( &mesh, this->materials, this->model, this->scaling );
        Interp_Grid2P_centroids2( particles, particles.P,    &mesh, mesh.p_in, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &this->model );
        Interp_Grid2P_centroids2( particles, particles.T,    &mesh, mesh.T,    mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCt.type, &this->model );
        NonNewtonianViscosityGrid (     &mesh, &this->materials, &this->model, Nmodel, &this->scaling );

        // Print informations!
        printf("Number of phases : %d\n", this->model.Nb_phases);
        if ( this->model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n,       this->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d0        ", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n,        this->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_lith,     this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P litho   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p0_n,       this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.p_in,       this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T,          this->scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s,       this->scaling.S,   (mesh.Nx-0)*(mesh.Nz-0), "mu_s      ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n,       this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n      ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_s,      this->scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.eta_phys_s, this->scaling.eta, (mesh.Nx-0)*(mesh.Nz-0), "eta_phys_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.eta_phys_n, this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho_s,      this->scaling.rho, (mesh.Nx-0)*(mesh.Nz-0), "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n,      this->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            for (int p=0; p<this->model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }
        }

        printf("*************************************\n");
        printf("******** Initialize timestep ********\n");
        printf("*************************************\n");

        DefineInitialTimestep( &this->model, &mesh, particles, this->materials, this->scaling );

        if ( this->model.rec_T_P_x_z == 1 ) {
            ArrayEqualArray( particles.T0,   particles.T, particles.Nb_part );
            ArrayEqualArray( particles.P0,   particles.P, particles.Nb_part );
            ArrayEqualArray( particles.x0,   particles.x, particles.Nb_part );
            ArrayEqualArray( particles.z0,   particles.z, particles.Nb_part );
            ArrayEqualArray( particles.Tmax, particles.T, particles.Nb_part );
            ArrayEqualArray( particles.Pmax, particles.P, particles.Nb_part );
        }
        ComputeMeanQuantitesForTimeSeries( &mesh );
        LogTimeSeries( &mesh, this->model, this->scaling );

        printf("*************************************\n");
        printf("*** Write initial file or restart ***\n");
        printf("*************************************\n");

        // Write initial output
        if ( writer == 1 ) {
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Output",  this->materials, this->scaling );
            if ( this->model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particles",  this->materials, this->scaling );
        }

        // Set initial stresses and pressure to zero
        Initialise1DArrayDouble( mesh.sxxd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.szzd,  (mesh.Nx-1)*(mesh.Nz-1), 0.0 );
        Initialise1DArrayDouble( mesh.sxz,   (mesh.Nx)  *(mesh.Nz)  , 0.0 );
        // Generate deformation maps
        if ( this->model.def_maps == 1 ) GenerateDeformationMaps( &mesh, &this->materials, &this->model, Nmodel, &this->scaling );
        particles.Nb_part_ini = particles.Nb_part;
    }
    else {
        // Which step do we restart from (BreakpointXXXX.dat)
        this->model.step = istep;
        printf("Restarting from step number %05d...\n", this->model.step);
        LoadBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, &this->model, &topo, &topo_ini, this->scaling  );
        SetGridCoordinates( &mesh, &this->model, this->model.Nx, this->model.Nz ); // Overwrite previous grid
    }

    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                             TIME LOOP : en avant les Doud'series !
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//

    this->model.step += 1;

    if ( this->model.GNUplot_residuals == 1 ) GNUplotPipe = popen ("gnuplot -persistent", "w");

    for (; this->model.step<=this->model.Nt; this->model.step++) {

        printf("*****************************************************\n");
        printf("****************** Time step %05d ******************\n", this->model.step);
        printf("*****************************************************\n");

        //------------------------------------------------------------------------------------------------------------------------------//

        t_omp_step = (double)omp_get_wtime();

        // Old time step
        this->model.dt0 = this->model.dt;

        // Define new time step
        EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &this->model, this->scaling, &mesh, 0 );
        printf("Selected dt = %2.2e\n", this->model.dt*this->scaling.t);

        // Save initial dt
        this->model.dt0 = this->model.dt;

        // Track particule generation index
        Initialise1DArrayInt( particles.generation,   particles.Nb_part  , 0 );

        //------------------------------------------------------------------------------------------------------------------------------//

        // Remove particles that would be above the surface
        if ( this->model.free_surf == 1 ) {
            CleanUpSurfaceParticles( &particles, &mesh, topo, this->scaling );
            CellFlagging( &mesh, this->model, topo, this->scaling );
            ArrayEqualArray(     topo.height0,     topo.height, mesh.Nx );
            ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );
        }

        // Interpolate material properties from particles to nodes
        t_omp = (double)omp_get_wtime();

        // Energy - interpolate thermal parameters and advected energy

        // Get energy and related material parameters from particles
        P2Mastah( &this->model, particles, this->materials.Cv,     &mesh, mesh.Cv,     mesh.BCp.type,  0, 0, interp, cent, 1);
        P2Mastah( &this->model, particles, this->materials.Qr,     &mesh, mesh.Qr,     mesh.BCp.type,  0, 0, interp, cent, 1);

        P2Mastah ( &this->model, particles, this->materials.k_eff, &mesh, mesh.kx, mesh.BCu.type,  0, 0, interp, vxnodes, 1);
        P2Mastah ( &this->model, particles, this->materials.k_eff, &mesh, mesh.kz, mesh.BCv.type,  0, 0, interp, vznodes, 1);

        // Get T and dTdt from previous step from particles
        P2Mastah( &this->model, particles, particles.T,     &mesh, mesh.T0_n,     mesh.BCp.type,  1, 0, interp, cent, 1);
        P2Mastah( &this->model, particles, particles.divth, &mesh, mesh.divth0_n, mesh.BCp.type,  1, 0, interp, cent, 1);

        // Make sure T is up to date for rheology evaluation
        ArrayEqualArray( mesh.T, mesh.T0_n, (mesh.Nx-1)*(mesh.Nz-1) );

        //-----------------------------------------------------------------------------------------------------------
        // Interp P --> p0_n , p0_s
        P2Mastah( &this->model, particles, particles.P,     &mesh, mesh.p0_n,   mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);

        // Get physical properties that are constant throughout each timestep
        UpdateDensity( &mesh, &particles, &this->materials, &this->model, &this->scaling );

        // Free surface - subgrid density correction
        if ( this->model.free_surf == 1 ) {
            SurfaceDensityCorrection( &mesh, this->model, topo, this->scaling  );
        }

        // Lithostatic pressure
        ArrayEqualArray(  mesh.p_lith0,   mesh.p_lith, Ncx*Ncz );
        ComputeLithostaticPressure( &mesh, &this->model, this->materials.rho[0], this->scaling, 1 );

        // Elasticity - interpolate advected/rotated stresses
        if  ( this->model.iselastic == 1 ) {

            // Get old stresses from particles
            if (this->model.StressUpdate==1)                     OldDeviatoricStressesPressure( &mesh, &particles, this->scaling, &this->model );

            if (this->model.StressUpdate==0){
                P2Mastah( &this->model, particles, particles.sxxd,    &mesh, mesh.sxxd0, mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
                P2Mastah( &this->model, particles, particles.szzd,    &mesh, mesh.szzd0, mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
                P2Mastah( &this->model, particles, particles.sxz,     &mesh, mesh.sxz0,  mesh.BCg.type,  1, 0, interp, vert, this->model.itp_stencil);
            }

            InterpCentroidsToVerticesDouble( mesh.sxxd0, mesh.sxxd0_s, &mesh, &this->model );
            InterpCentroidsToVerticesDouble( mesh.szzd0, mesh.szzd0_s, &mesh, &this->model );
            InterpVerticesToCentroidsDouble( mesh.sxz0_n,  mesh.sxz0,  &mesh, &this->model );

            // Interpolate shear modulus
            ShearModCompExpGrid( &mesh, this->materials, this->model, this->scaling );
        }

        // Director vector
        if (this->model.aniso == 1 ) {
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d1_n , mesh.BCp.type, -1, 0, interp, cent, this->model.itp_stencil);
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d2_n , mesh.BCp.type, -2, 0, interp, cent, this->model.itp_stencil);
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d1_s , mesh.BCg.type, -1, 0, interp, vert, this->model.itp_stencil);
            P2Mastah( &this->model, particles, NULL, &mesh, mesh.d2_s , mesh.BCg.type, -2, 0, interp, vert, this->model.itp_stencil);
            FiniteStrainAspectRatio ( &mesh, this->scaling, this->model, &particles );
            P2Mastah( &this->model, particles, this->materials.aniso_factor,     &mesh, mesh.aniso_factor_n , mesh.BCp.type,  0, 0, interp, cent, this->model.itp_stencil);
            P2Mastah( &this->model, particles, this->materials.aniso_factor,     &mesh, mesh.aniso_factor_s , mesh.BCg.type,  0, 0, interp, vert, this->model.itp_stencil);
        }

        P2Mastah( &this->model, particles, particles.X,     &mesh, mesh.X0_n , mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
        P2Mastah( &this->model, particles, particles.X,     &mesh, mesh.X0_s , mesh.BCg.type,  1, 0, interp, vert, this->model.itp_stencil);

        P2Mastah( &this->model, particles, particles.noise, &mesh, mesh.noise_s, mesh.BCg.type,  1, 0, interp, vert, this->model.itp_stencil);
        P2Mastah( &this->model, particles, particles.noise, &mesh, mesh.noise_n, mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);

        // Diffuse rheological contrasts
        //            if (this->model.diffuse_X == 1) {
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.X0_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Diffuse_X(&mesh, &this->model, &this->scaling);
        //                Interp_Grid2P( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );
        //                Interp_P2C ( particles, particles.X, &mesh, mesh.Xreac_n, mesh.xg_coord, mesh.zg_coord, 1, 0 );
        //                Interp_P2N ( particles, particles.X, &mesh, mesh.Xreac_s, mesh.xg_coord, mesh.zg_coord, 1, 0, &this->model );
        //            }

        // Interpolate Grain size
        P2Mastah( &this->model, particles, particles.d,     &mesh, mesh.d0_n , mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);
        ArrayEqualArray(  mesh.d_n,  mesh.d0_n, Ncx*Ncz );

        // Interpolate Melt fraction
        P2Mastah( &this->model, particles, particles.phi,   &mesh, mesh.phi0_n , mesh.BCp.type,  1, 0, interp, cent, this->model.itp_stencil);

        //-------------------------------------------------------------------------------------------------------------

        // Compute cohesion and friction angle on the grid
        CohesionFrictionDilationGrid( &mesh, &particles, this->materials, this->model, this->scaling );

        // Detect compressible cells
        if ( this->model.compressible == 1 ) DetectCompressibleCells ( &mesh, &this->model );

        // Compute cohesion and friction angle on the grid
        CohesionFrictionDilationGrid( &mesh, &particles, this->materials, this->model, this->scaling );

        // Detect compressible cells
        if (this->model.compressible == 1) DetectCompressibleCells ( &mesh, &this->model );

        // Min/Max interpolated fields
        if ( this->model.noisy == 1 ) {
            MinMaxArray(particles.sxxd, this->scaling.S, particles.Nb_part, "sxxd part  ");
            MinMaxArray(particles.T,   this->scaling.T,   particles.Nb_part, "T part    ");
            MinMaxArrayTag( mesh.p0_n,         this->scaling.S,    (mesh.Nx-1)*(mesh.Nz-1),   "p0_n",   mesh.BCp.type );
            MinMaxArrayTag( mesh.sxz0,     this->scaling.S,   (mesh.Nx)*(mesh.Nz),     "sxz0    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.sxxd0,    this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "sxx0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.szzd0,    this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "szz0    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.mu_s,     this->scaling.S,   (mesh.Nx)*(mesh.Nz),     "mu_s    ", mesh.BCg.type );
            MinMaxArrayTag( mesh.mu_n,     this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "mu_n    ", mesh.BCp.type );
            MinMaxArrayTag( mesh.C_s,      this->scaling.S,   (mesh.Nx)*(mesh.Nz),     "C_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.C_n,      this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "C_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.fric_s,   180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "fric_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.fric_n,   180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "fric_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.dil_s,    180.0/M_PI,  (mesh.Nx)*(mesh.Nz),     "dil_s   ", mesh.BCg.type );
            MinMaxArrayTag( mesh.dil_n,    180.0/M_PI,  (mesh.Nx-1)*(mesh.Nz-1), "dil_n   ", mesh.BCp.type );
            MinMaxArrayTag( mesh.strain_s,   1.0,       (mesh.Nx)*(mesh.Nz),     "strain_s", mesh.BCg.type );
            MinMaxArrayTag( mesh.strain_n,   1.0,       (mesh.Nx-1)*(mesh.Nz-1), "strain_n", mesh.BCp.type );
            MinMaxArrayTag( mesh.bet_s,    1.0/this->scaling.S,   (mesh.Nx)*(mesh.Nz),     "beta_s  ", mesh.BCg.type );
            MinMaxArrayTag( mesh.bet_n,    1.0/this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "beta_n  ", mesh.BCp.type );
            MinMaxArrayTag( mesh.T0_n,     this->scaling.T,   (mesh.Nx-1)*(mesh.Nz-1), "T       ", mesh.BCt.type );
            MinMaxArrayTag( mesh.p_in,     this->scaling.S,   (mesh.Nx-1)*(mesh.Nz-1), "P       ", mesh.BCt.type );
            MinMaxArrayI  ( mesh.comp_cells, 1.0, (mesh.Nx-1)*(mesh.Nz-1), "comp_cells" );
            MinMaxArrayTag( mesh.rho_s,      this->scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
            MinMaxArrayTag( mesh.rho_n,      this->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
            MinMaxArrayTag( mesh.rho0_n,     this->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );

            for (int p=0; p<this->model.Nb_phases; p++) {
                printf("Phase number %d:\n", p);
                MinMaxArrayTag( mesh.phase_perc_n[p],    1.0, (mesh.Nx-1)*(mesh.Nz-1), "ph_n      ", mesh.BCp.type );
                MinMaxArrayTag( mesh.phase_perc_s[p],    1.0, (mesh.Nx-0)*(mesh.Nz-0), "ph_s      ", mesh.BCg.type );
            }

            if  ( this->model.aniso == 1 ) MinMaxArrayTag( mesh.FS_AR_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "FS_AR_n", mesh.BCp.type );
            if  ( this->model.aniso == 1 ) MinMaxArrayTag( mesh.FS_AR_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "FS_AR_s", mesh.BCg.type );
            if  ( this->model.aniso == 1 ) MinMaxArrayTag( mesh.aniso_factor_n,  1.0,   (mesh.Nx-1)*(mesh.Nz-1), "aniso_factor_n", mesh.BCp.type );
            if  ( this->model.aniso == 1 ) MinMaxArrayTag( mesh.aniso_factor_s,  1.0,   (mesh.Nx)*(mesh.Nz),     "aniso_factor_s", mesh.BCg.type );
        }

        printf("** Time for particles interpolations I = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );

//        struct timespec begin, end;
//        clock_gettime(CLOCK_REALTIME, &begin);
//
//        ViscosityDerivatives( &mesh, &this->materials, &this->model, Nmodel, &this->scaling );
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

        if ( this->model.ismechanical == 1 ) {

          if (this->crazyConductivity) {
            for (int n = 0; n < this->crazyConductivity->nPhases; n++) {
              const int    phase           = this->crazyConductivity->phases[n];
              this->materials.k_eff[phase] = this->materials.k[phase];
            }
            printf("Running with normal conductivity for the asthenosphere...\n");
          }

            // Allocate and initialise solution and RHS vectors
            SetBCs( this, &mesh);

            // Reset fields and BC values if needed
            //        if ( this->model.ispureshear_ale == 1 ) InitialiseSolutionFields( &mesh, &this->model );
            InitialiseSolutionFields( &mesh, &this->model );
            EvalNumberOfEquations( &mesh, &Stokes );
            if ( IsNewtonStep == 1 ) EvalNumberOfEquations( &mesh, &Jacob  );
            SAlloc( &Stokes,  Stokes.neq );
            if ( IsNewtonStep == 1 ) SAlloc(  &Jacob,   Jacob.neq );
            if ( this->model.decoupled_solve == 1 ) {
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
            InitialiseSolutionVector( &mesh, &Stokes, &this->model );

            //------------------------------------------------------------------------------------------------------------------------------//

            // Non-linear iteration cycle
            Nmodel.nit       = 0;
            this->model.nit        = 0;
            Nmodel.stagnated = 0;
            nstag            = 0;
            IsJacobianUsed    = 0;

            ArrayEqualArray( mesh.p_start,    mesh.p_in,      (mesh.Nx-1)*(mesh.Nz-1) );
            ArrayEqualArray( mesh.u_start,    mesh.u_in,      (mesh.Nx)  *(mesh.Nz+1) );
            ArrayEqualArray( mesh.v_start,    mesh.v_in,      (mesh.Nx+1)*(mesh.Nz)   );
            //            ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );

            // Set up solver context
            if ( this->model.decoupled_solve == 1 ) {
                cholmod_start( &CholmodSolver.c );
                if ( Nmodel.nit==0 ) CholmodSolver.Analyze = 1;
                printf("Run CHOLMOD analysis yes/no: %d \n", CholmodSolver.Analyze);
            }

            int Nmax_picard = Nmodel.nit_max;

            // If Picard 2 Newton is activated: force initial Newton step
            if ( IsNewtonStep == 1 && Nmodel.Picard2Newton == 1 ) {
                this->model.Newton = 0;
                IsFirstNewtonStep  = 1;
            }

            while ( Nmodel.nit <= Nmax_picard && nstag<this->model.nstagmax) {

                if ( Nmodel.nit > 0 && Nmodel.Picard2Newton == 1 ) {
                    if ( Nmodel.rx_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.rz_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.rp_rel[Nmodel.nit-1] < Nmodel.Pic2NewtCond || Nmodel.nit>= Nmodel.nit_Pic_max) {
                        if ( IsFirstNewtonStep == 0 ) CholmodSolver.Analyze = 0;
                        if ( IsFirstNewtonStep == 1 ) {
                            cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                            CholmodSolver.Analyze = 1;
                            IsFirstNewtonStep    = 0;
                        }
                        this->model.Newton = 1;
                    }
                }

                // Determine whether Jacobian matrix should be assembled
                if ( this->model.Newton == 1 || this->model.aniso == 1 ) IsJacobianUsed = 1;

                printf("**********************************************\n");
                if ( this->model.Newton == 1 ) {
                    printf("*** Newton it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, this->model.step);
                    Nmodel.LogIsNewtonStep[Nmodel.nit] = 1;
                }
                else {
                    printf("*** Picard it. %02d of %02d (step = %05d) ***\n", Nmodel.nit, Nmodel.nit_max, this->model.step);
                    Nmodel.LogIsNewtonStep[Nmodel.nit] = 0;
                }
                printf("**********************************************\n");

                // Update non-linear rheology
                UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, this->materials, &this->model, &Nmodel, this->scaling, 0, 0.0 );
                RheologicalOperators( &mesh, &this->model, &this->scaling, 0 );                               // ??????????? déjà fait dans UpdateNonLinearity
                NonNewtonianViscosityGrid (     &mesh, &this->materials, &this->model, Nmodel, &this->scaling );    // ??????????? déjà fait dans UpdateNonLinearity

                if ( this->model.noisy == 1 ) {
                    MinMaxArrayTag( mesh.T,         this->scaling.T, (mesh.Nx-1)*(mesh.Nz-1), "T         ", mesh.BCt.type );
                    MinMaxArrayTag( mesh.p0_n,      this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P old     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.p_in,      this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "P         ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.div_u,     this->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "div       ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exxd,      this->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "exxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.ezzd,      this->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "ezzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.exz,       this->scaling.E, (mesh.Nx-0)*(mesh.Nz-0), "exz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.sxxd,      this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "sxxd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.szzd,      this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "szzd      ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.sxz,       this->scaling.S, (mesh.Nx-0)*(mesh.Nz-0), "sxz       ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.eta_phys_s, this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "eta_phys_s", mesh.BCg.type );
                    MinMaxArrayTag( mesh.eta_phys_n, this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1), "eta_phys_n", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho_s,      this->scaling.rho, (mesh.Nx)*(mesh.Nz),     "rho_s     ", mesh.BCg.type );
                    MinMaxArrayTag( mesh.rho_n,      this->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho_n     ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.rho0_n,     this->scaling.rho, (mesh.Nx-1)*(mesh.Nz-1), "rho0_n    ", mesh.BCp.type );
                    MinMaxArrayTag( mesh.d_n,        this->scaling.L,   (mesh.Nx-1)*(mesh.Nz-1), "d         ", mesh.BCp.type );
                }

                if ( this->model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Output_BeforeSolve", this->materials, this->scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particles_BeforeSolve", this->materials, this->scaling );
                }

                // Build discrete system of equations - Linearised Picard
                if ( this->model.decoupled_solve == 0 ) BuildStokesOperator           ( &mesh, this->model, 0, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, 1 );
                if ( this->model.decoupled_solve == 1 ) BuildStokesOperatorDecoupled  ( &mesh, this->model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, 1 );

                // Build discrete system of equations - Jacobian
                ViscosityDerivatives( &mesh, &this->materials, &this->model, Nmodel, &this->scaling );
                if ( this->model.Newton == 1 && Nmodel.nit > 0 ) RheologicalOperators( &mesh, &this->model, &this->scaling, 1 );
                if ( IsJacobianUsed == 1 )                  BuildJacobianOperatorDecoupled( &mesh, this->model, 0, mesh.p_corr, mesh.p_in, mesh.u_in, mesh.v_in,  &Jacob,  &JacobA,  &JacobB,  &JacobC,   &JacobD, 1 );

                // IsNanArray2DFP(mesh.eta_n, (mesh.Nx-1)*(mesh.Nz-1));
                // IsInfArray2DFP(mesh.eta_n, (mesh.Nx-1)*(mesh.Nz-1));

                // IsNanArray2DFP(mesh.eta_s, (mesh.Nx-0)*(mesh.Nz));
                // IsInfArray2DFP(mesh.eta_s, (mesh.Nx-0)*(mesh.Nz));
                //                MinMaxArrayTag( mesh.detadexx_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadexx_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadezz_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadezz_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadexx_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "detadgxz_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadexx_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadexx_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadezz_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadezz_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.detadgxz_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "detadgxz_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D11_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D21_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D12_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D22_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D13_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D23_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D14_n,      this->scaling.eta, (mesh.Nx-1)*(mesh.Nz-1),     "D24_n     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D31_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D31_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D32_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D32_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D33_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D33_s     ", mesh.BCg.type );
                //                MinMaxArrayTag( mesh.D34_s,      this->scaling.eta, (mesh.Nx)*(mesh.Nz),     "D34_s     ", mesh.BCg.type );

                // Diagonal this->scaling
                if ( this->model.diag_scaling ) {
                    if ( this->model.Newton          == 0 )  ExtractDiagonalScale( &StokesA, &StokesB, &StokesC, &StokesD );
                    if ( this->model.Newton          == 1 )  ExtractDiagonalScale( &JacobA,  &JacobB,  &JacobC,   &JacobD );
                    if ( this->model.Newton          == 1 )  ArrayEqualArray(StokesA.d, JacobA.d, StokesA.neq);
                    if ( this->model.Newton          == 1 )  ArrayEqualArray(StokesC.d, JacobC.d, StokesC.neq);
                    ScaleMatrix( &StokesA, &StokesB, &StokesC, &StokesD );
                }

                // Needs to be done after matrix assembly since diagonal this->scaling is used in there
                printf("---- Non-linear residual ----\n");
                RheologicalOperators( &mesh, &this->model, &this->scaling, 0 );
                if ( this->model.decoupled_solve == 0 ) EvaluateStokesResidual( &Stokes, &Nmodel, &mesh, this->model, this->scaling, 0 );
                if ( this->model.decoupled_solve == 1 ) EvaluateStokesResidualDecoupled( &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &Nmodel, &mesh, this->model, this->scaling, 0 );
                printf("---- Non-linear residual ----\n");

                if ( this->model.Newton == 1 && this->model.aniso == 0 ) {
                    ArrayEqualArray( JacobA.F, StokesA.F,  StokesA.neq );
                    ArrayEqualArray( JacobC.F, StokesC.F,  StokesC.neq );
                }
                if ( this->model.Newton == 1 && this->model.diag_scaling ) ScaleMatrix( &JacobA,  &JacobB,  &JacobC,  &JacobD  );

                // Store residuals
                Nmodel.resx_f = Nmodel.resx; Nmodel.rx_abs[Nmodel.nit] = Nmodel.resx; Nmodel.rx_rel[Nmodel.nit] = Nmodel.resx/Nmodel.resx0;
                Nmodel.resz_f = Nmodel.resz; Nmodel.rz_abs[Nmodel.nit] = Nmodel.resz; Nmodel.rz_rel[Nmodel.nit] = Nmodel.resz/Nmodel.resz0;
                Nmodel.resp_f = Nmodel.resp; Nmodel.rp_abs[Nmodel.nit] = Nmodel.resp; Nmodel.rp_rel[Nmodel.nit] = Nmodel.resp/Nmodel.resp0;

                if ( this->model.write_debug == 1 ) WriteResiduals( mesh, this->model, Nmodel, this->scaling );

                // if pass --> clear matrix break
                if ( (Nmodel.resx < Nmodel.abs_tol_u || Nmodel.resx/Nmodel.resx0 < Nmodel.rel_tol_u) && (Nmodel.resz < Nmodel.abs_tol_u || Nmodel.resz/Nmodel.resz0 < Nmodel.rel_tol_u) && (Nmodel.resp < Nmodel.abs_tol_p || Nmodel.resp/Nmodel.resp0 < Nmodel.rel_tol_p) ) {
                    printf( "Non-linear solver converged to abs_tol_u = %2.2e abs_tol_p = %2.2e rel_tol_u = %2.2e rel_tol_p = %2.2e\n", Nmodel.abs_tol_u, Nmodel.abs_tol_p, Nmodel.rel_tol_u, Nmodel.rel_tol_p );
                    FreeSparseSystems( IsJacobianUsed, this->model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                    break;
                }

                // Direct solve
                t_omp = (double)omp_get_wtime();
                if ( this->model.decoupled_solve == 0 ) SolveStokesDefect( &Stokes, &CholmodSolver, &Nmodel, &mesh, &this->model, &particles, &topo_chain, &topo, this->materials, this->scaling );
                if ( this->model.decoupled_solve == 1 ) {

                    if ( IsJacobianUsed==1 ) SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &this->model, &particles, &topo_chain, &topo, this->materials, this->scaling,  &JacobA,  &JacobB,  &JacobC );
                    else                     SolveStokesDefectDecoupled( &StokesA, &StokesB, &StokesC, &StokesD, &Stokes, &CholmodSolver, &Nmodel, &mesh, &this->model, &particles, &topo_chain, &topo, this->materials, this->scaling, &StokesA, &StokesB, &StokesC );

                    if ( Nmodel.stagnated == 1 && this->model.safe_mode <= 0 ) {
                        printf( "Non-linear solver stagnated to res_u = %2.2e res_z = %2.2e\n", Nmodel.resx_f, Nmodel.resz_f );
                        printf( "You may want to try setting line_search_min > 0.0\n Good luck good man!\n");
                        FreeSparseSystems( IsJacobianUsed, this->model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                        break;
                    }
                }

                if ( Nmodel.stagnated == 1 && this->model.iselastic == 1 && this->model.safe_mode == 1 ) {
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (abs_tol_u = %2.2e abs_tol_p = %2.2e)\e[m\n", Nmodel.abs_tol_u, Nmodel.abs_tol_p );
                    printf( "\e[1;31mWARNING : Non-linear solver stagnated (rel_tol_u = %2.2e rel_tol_p = %2.2e)\e[m\n", Nmodel.rel_tol_u, Nmodel.rel_tol_p );
                    printf( "\e[1;31mReducing the timestep, and restart the iterations cycle...\e[m\n");
                    printf( "Before reduction: this->model.dt =, %2.2e\n", this->model.dt*this->scaling.t);
                    // ----------------------
                    this->model.dt /= this->model.safe_dt_div;
                    printf( "Timestep divided by %2.2f => NEW CURRENT this->model.dt =, %2.2e\n", this->model.safe_dt_div, this->model.dt*this->scaling.t);
                    Nmodel.stagnated = 0;
                    // ----------------------
                    Nmodel.nit = -1;
                    printf( "Restart solutions\n");
                    ArrayEqualArray( mesh.p_in, mesh.p_start,  (mesh.Nx-1)*(mesh.Nz-1) );
                    ArrayEqualArray( mesh.u_in, mesh.u_start,  (mesh.Nx)  *(mesh.Nz+1) );
                    ArrayEqualArray( mesh.v_in, mesh.v_start,  (mesh.Nx+1)*(mesh.Nz)   );
                    nstag++;
                    //-----------------------
                    printf( "nstag value = %02d - nstagmax = %02d\n", nstag, this->model.nstagmax);
                    if (nstag==this->model.nstagmax) {
                        printf( "CheckDoudzOut!!\n");
                        exit(0);
                    //-----------------------
                    }
                }
                FreeSparseSystems( IsJacobianUsed, this->model.decoupled_solve, &Stokes, &StokesA, &StokesB, &StokesC, &StokesD, &JacobA, &JacobB, &JacobC, &JacobD );
                Nmodel.nit++; this->model.nit = Nmodel.nit;
            }

            // Clean solver context
            if ( this->model.decoupled_solve == 1 && Nmodel.nit>0 ) {
                printf("Cleaning up Cholesky factors --- nit = %02d\n",  Nmodel.nit);
                cholmod_free_factor ( &CholmodSolver.Lfact, &CholmodSolver.c);
                cholmod_finish( &CholmodSolver.c );
            }

            // Update rheology
            UpdateNonLinearity( &mesh, &particles, &topo_chain, &topo, this->materials, &this->model, &Nmodel, this->scaling, 0, 1.0 );
            ArrayEqualArray( mesh.p_in, mesh.p_corr, (mesh.Nx-1)*(mesh.Nz-1) ); // make sure corrected pressure is used

            // Min/Max velocities
            MinMaxArray( mesh.u_in,  this->scaling.V, (mesh.Nx)*(mesh.Nz+1),   "Vx. grid" );
            MinMaxArray( mesh.v_in,  this->scaling.V, (mesh.Nx+1)*(mesh.Nz),   "Vz. grid" );
            MinMaxArray( mesh.p_in,  this->scaling.S, (mesh.Nx-1)*(mesh.Nz-1), "       P" );
            MinMaxArray( mesh.div_u, this->scaling.E, (mesh.Nx-1)*(mesh.Nz-1), "  div(V)" );

            printf("--------------------------------------------------------------\n");
            int i, nit;
            if (Nmodel.nit>=Nmodel.nit_max)  nit = Nmodel.nit_max;
            if (Nmodel.nit< Nmodel.nit_max)  nit = Nmodel.nit;
            if (Nmodel.Picard2Newton == 1 )  printf("Picard 2 Newton is activated with condition: %2.2e\n", Nmodel.Pic2NewtCond);
            for (i=0; i<=nit; i++) {
                if (Nmodel.LogIsNewtonStep[i] == 1) printf("New. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, Nmodel.rx_abs[i], Nmodel.rz_abs[i], Nmodel.rp_abs[i], Nmodel.rx_rel[i], Nmodel.rz_rel[i], Nmodel.rp_rel[i]);
                else                         printf("Pic. it. %02d: abs: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e --- rel: |Fx| = %2.2e - |Fz| = %2.2e - |Fp| = %2.2e\n", i, Nmodel.rx_abs[i], Nmodel.rz_abs[i], Nmodel.rp_abs[i], Nmodel.rx_rel[i], Nmodel.rz_rel[i], Nmodel.rp_rel[i]);
                if (i == Nmodel.nit_max && this->model.safe_mode == 1) {
                    printf("Exit: Max iteration reached: Nmodel.nit_max = %02d! Check what you wanna do now...\n",Nmodel.nit_max);
                    if ( (Nmodel.resx < Nmodel.abs_tol_u) && (Nmodel.resz < Nmodel.abs_tol_u) && (Nmodel.resp < Nmodel.abs_tol_p) ) {}
                    else exit(1);
                }
            }
            printf("--------------------------------------------------------------\n");

            // plot residuals
            if ( this->model.GNUplot_residuals == 1 ) { // && this->model.step % writer_step == 0

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
        if ( this->model.free_surf == 1 ) {
            SurfaceVelocity( &mesh, this->model, &topo, &topo_chain, this->scaling );
            MinMaxArray( topo_chain.Vx,      this->scaling.V, topo_chain.Nb_part,       "Vx surf." );
            MinMaxArray( topo_chain.Vz,      this->scaling.V, topo_chain.Nb_part,       "Vz surf." );
            MinMaxArray( topo_chain_ini.Vx,  this->scaling.V, topo_chain_ini.Nb_part,   "Vx surf. ini." );
            MinMaxArray( topo_chain_ini.Vz,  this->scaling.V, topo_chain_ini.Nb_part,   "Vz surf. ini." );
        }
        // THIS IS ACTIVATED JUST FOR TOPO.VX IN OUTPUT FILE

        if (this->model.StressUpdate==1) TotalStresses( &mesh, &particles, this->scaling, &this->model );

        // Update stresses on markers
        if (this->model.StressUpdate==0) UpdateParticleStress(  &mesh, &particles, &this->model, &this->materials, &this->scaling );

        // Update pressure on markers
        UpdateParticlePressure( &mesh, this->scaling, this->model, &particles, &this->materials );
        //        Interp_Grid2P_centroids( particles, particles.P, &mesh, mesh.p_in, mesh.xc_coord,  mesh.zc_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &this->model );
        UpdateParticleX( &mesh, this->scaling, this->model, &particles, &this->materials );

        // Grain size evolution
        UpdateParticleGrainSize( &mesh, this->scaling, this->model, &particles, &this->materials );

        if ( this->model.noisy == 1 ) {
            MinMaxArrayTag( mesh.d0_n   , this->scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d0", mesh.BCp.type );
            MinMaxArrayTag( mesh.d_n    , this->scaling.L, (mesh.Nx-1)*(mesh.Nz-1), "d ", mesh.BCp.type );
            MinMaxArrayPart( particles.d, this->scaling.L, particles.Nb_part, "d on markers", particles.phase ) ;
        }

        // Update phi on the particles
        UpdateParticlePhi( &mesh, this->scaling, this->model, &particles, &this->materials );


        //------------------------------------------------------------------------------------------------------------------------------//

        if (this->model.isthermal == 1 ) {

            printf("*************************************\n");
            printf("*********** Thermal solver **********\n");
            printf("*************************************\n");

            t_omp = (double)omp_get_wtime();

            // Matrix assembly and direct solve
            EnergyDirectSolve( &mesh, this->model,  mesh.T,  mesh.dT,  mesh.rhs_t, mesh.T, &particles, this->model.dt, this->model.shear_heat, this->model.adiab_heat, this->scaling, 1 );
            MinMaxArray(particles.T, this->scaling.T, particles.Nb_part, "T part. before UpdateParticleEnergy");

            // Update energy on particles
            UpdateParticleEnergy( &mesh, this->scaling, this->model, &particles, &this->materials );
            MinMaxArray(particles.T, this->scaling.T, particles.Nb_part, "T part. after UpdateParticleEnergy");

            printf("** Time for Thermal solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
        }

        //--------------------------------------------------------------------------------------------------------------------------------//

        if (this->model.ProgReac == 1 ) {

            printf("*************************************\n");
            printf("********** Chemical solver **********\n");
            printf("*************************************\n");

            P2Mastah ( &this->model, particles, this->materials.k_chem, &mesh, mesh.kc_x, mesh.BCu.type,  0, 0, interp, vxnodes, this->model.itp_stencil);
            P2Mastah ( &this->model, particles, this->materials.k_chem, &mesh, mesh.kc_z, mesh.BCv.type,  0, 0, interp, vznodes, this->model.itp_stencil);
            ChemicalDirectSolve( &mesh, this->model, &particles, &this->materials, this->model.dt, this->scaling );

        }
        //--------------------------------------------------------------------------------------------------------------------------------//

        // Update maximum pressure and temperature on markers
        if ( this->model.rec_T_P_x_z == 1 )  UpdateMaxPT( this->scaling, this->model, &particles );

        ComputeMeanQuantitesForTimeSeries( &mesh );
        LogTimeSeries( &mesh, this->model, this->scaling );

        //------------------------------------------------------------------------------------------------------------------------------//


        if ( this->model.advection == 1 ) {

            printf("*************************************\n");
            printf("************** Advection ************\n");
            printf("*************************************\n");

            t_omp = (double)omp_get_wtime();
            double dt_solve = this->model.dt;
            int    nsub=1.0, isub;
            double dt_sub;
            EvaluateCourantCriterion( mesh.u_in, mesh.v_in, &this->model, this->scaling, &mesh, 0 );

            if ( this->model.dt < 0.95*dt_solve ) { // if dt advection is lower than dt used for thermo-mechanical solve then split dt
                nsub   = ceil(dt_solve/this->model.dt);
                dt_sub = dt_solve / nsub;
                printf("dt advection = %2.9e --- dt_solve = %2.9e\n", this->model.dt*this->scaling.t, dt_solve*this->scaling.t );
                printf("dt is %lf larger than dt_solve: need %d substeps of %2.2e s\n", dt_solve/this->model.dt, nsub , dt_sub*this->scaling.t );
                printf("So: nsub*dt_sub = %2.2e for dt_solve = %2.2e\n", nsub*dt_sub*this->scaling.t, dt_solve*this->scaling.t);
                this->model.dt = dt_sub;
            }
            else {
                nsub     = 1; // if dt advection is larger than dt used for thermo-mechanical solve then, forget it, and set dt to dt_solve
                this->model.dt = dt_solve;
            }

            // Loop on substeps
            for (isub=0;isub<nsub;isub++) {

                printf("************** Advection step %03d of %03d: dtsub = %2.2e **************\n", isub, nsub, this->model.dt*this->scaling.t );

                // Advect domain boundaries
                if ( this->model.ispureshear_ale > 0 ) {
                    PureShearALE( &this->model, &mesh, &topo_chain, this->scaling );
                }

                if ( this->model.free_surf == 1 ) {
                    // Save old topo
                    ArrayEqualArray( topo_chain.z0, topo_chain.z, topo_chain.Nb_part ); // save old z
                    ArrayEqualArray( topo_chain_ini.z0, topo_chain_ini.z, topo_chain.Nb_part ); // save old z


                    // Advect free surface with RK4
                    RogerGuntherII( &topo_chain,     this->model, mesh, 1, this->scaling );
                    RogerGuntherII( &topo_chain_ini, this->model, mesh, 1, this->scaling );
                }

                //                // Advect free surface
                //                if ( this->model.free_surf == 1 ) {
                //                    AdvectFreeSurf( &topo_chain,     this->model, this->scaling );
                //                    AdvectFreeSurf( &topo_chain_ini, this->model, this->scaling );
                //                    MinMaxArray( topo_chain.z,      this->scaling.L, topo_chain.Nb_part,       "z surf.     " );
                //                    MinMaxArray( topo_chain_ini.z,  this->scaling.L, topo_chain_ini.Nb_part,   "z surf. ini." );
                //                }

                // Correction for particle inflow 0
                if (this->model.ispureshear_ale == -1 && this->model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, this->model, topo, 0 );

                // Advect fluid particles
                RogerGuntherII( &particles, this->model, mesh, 1, this->scaling );

                // Correction for particle inflow 1
                if (this->model.ispureshear_ale == -1 && this->model.isperiodic_x == 0) ParticleInflowCheck( &particles, &mesh, this->model, topo, 1 );

                // Update accumulated strain
                AccumulatedStrainII( &mesh, this->scaling, this->model, &particles,  mesh.xc_coord,  mesh.zc_coord, mesh.Nx-1, mesh.Nz-1, mesh.BCp.type );

                // Update deformation gradient tensor components
                if ( this->model.fstrain == 1 ) DeformationGradient( mesh, this->scaling, this->model, &particles );

                if ( this->model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Output_BeforeSurfRemesh", this->materials, this->scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particles_BeforeSurfRemesh", this->materials, this->scaling );
                }

                if ( this->model.free_surf == 1 ) {

                    //                    for (int k=0;k<topo_chain.Nb_part;k++) {
                    //
                    //                        if (fabs(topo_chain.x[k]) < 0.5*this->model.surf_Winc ) {
                    //                            topo_chain.z[k] -= this->model.dt*this->model.surf_Vinc;
                    //                        }
                    //                    }
                    //

                    if ( this->model.surf_remesh == 1 ) {
                        // Get current topography
                        ProjectTopography( &topo,     &topo_chain,     this->model, mesh, this->scaling, mesh.xg_coord, 0 );
                        ProjectTopography( &topo_ini, &topo_chain_ini, this->model, mesh, this->scaling, mesh.xg_coord, 0 );
                        MarkerChainPolyFit( &topo,     &topo_chain,     this->model, mesh );
                        MarkerChainPolyFit( &topo_ini, &topo_chain_ini, this->model, mesh );

                        // Remesh free surface I
                        RemeshMarkerChain( &topo_chain,     &topo,     this->model, this->scaling, &mesh, 1 );
                        RemeshMarkerChain( &topo_chain_ini, &topo_ini, this->model, this->scaling, &mesh, 1 );
                    }


                    // Project topography on vertices
                    ProjectTopography( &topo,     &topo_chain,     this->model, mesh, this->scaling, mesh.xg_coord, 0 );
                    ProjectTopography( &topo_ini, &topo_chain_ini, this->model, mesh, this->scaling, mesh.xg_coord, 0 );
                    //                    ArrayEqualArray( topo.height0, topo.height, mesh.Nx );
                    //                    ArrayEqualArray( topo_ini.height0, topo_ini.height, mesh.Nx );

                    if ( this->model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Output_AfterSurfRemesh", this->materials, this->scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particles_AfterSurfRemesh", this->materials, this->scaling );
                    }

                    // Diffuse topography
                    if ( this->model.surf_processes >= 1 )  DiffuseAlongTopography( &mesh, this->model, this->scaling, topo.height, topo.height, mesh.Nx, 0.0, this->model.dt );

                    // Marker chain polynomial fit
                    MarkerChainPolyFit( &topo,     &topo_chain,     this->model, mesh );
                    CorrectTopoIni( &particles, this->materials, &topo_chain_ini, &topo, this->model, this->scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, this->model, mesh );

                    // Sedimentation
                    if ( this->model.surf_processes >= 2 ) {
                        AddPartSed( &particles, this->materials, &topo_chain, &topo, this->model, this->scaling, &mesh);
                        if (this->model.cpc==-1) CountPartCell_BEN( &particles, &mesh, this->model, topo, 0, this->scaling );
                        if (this->model.cpc== 0) CountPartCell_Old( &particles, &mesh, this->model, topo, 0, this->scaling );
                        if (this->model.cpc== 1) CountPartCell    ( &particles, &mesh, this->model, topo, topo_ini, 0, this->scaling );
                    }

                    if ( this->model.write_debug == 1 ) {
                        WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Outputx", this->materials, this->scaling );
                        WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particlesx", this->materials, this->scaling );
                    }

                    // Remesh free surface II
                    RemeshMarkerChain( &topo_chain,     &topo,     this->model, this->scaling, &mesh, 2 );
                    RemeshMarkerChain( &topo_chain_ini, &topo_ini, this->model, this->scaling, &mesh, 2 );
                    CorrectTopoIni( &particles, this->materials, &topo_chain_ini, &topo, this->model, this->scaling, &mesh);
                    MarkerChainPolyFit( &topo_ini, &topo_chain_ini, this->model, mesh );

                    // Remove particles that are above the surface
                    CleanUpSurfaceParticles( &particles, &mesh, topo, this->scaling );

                    // Call cell flagging routine for free surface calculations
                    CellFlagging( &mesh, this->model, topo, this->scaling );
                }

                printf("** Time for advection solver = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

#ifdef _HDF5_
                if ( this->model.write_debug == 1 ) {
                    WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Outputxx", this->materials, this->scaling );
                    WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particlesxx", this->materials, this->scaling );
                }
#endif

                //                printf("*************************************\n");
                //                printf("************** Reseeding ************\n");
                //                printf("*************************************\n");

                // Count the number of particle per cell
                t_omp = (double)omp_get_wtime();
                if (this->model.cpc ==-1)                       CountPartCell_BEN( &particles, &mesh, this->model, topo, 0, this->scaling );
                if (this->model.cpc == 0)                       CountPartCell_Old( &particles, &mesh, this->model, topo, 0, this->scaling );
                if (this->model.cpc == 1 && this->model.Reseed == 1 ) CountPartCell    ( &particles, &mesh, this->model, topo, topo_ini, 1, this->scaling );
                if (this->model.cpc == 1)                       CountPartCell    ( &particles, &mesh, this->model, topo, topo_ini, 0, this->scaling );

                if (this->model.cpc == 2 && this->model.Reseed == 1 ) CountPartCell2   ( &particles, &mesh, this->model, topo, topo_ini, 1, this->scaling );
                if (this->model.cpc == 2)                       CountPartCell2   ( &particles, &mesh, this->model, topo, topo_ini, 0, this->scaling );

                printf("After re-seeding :\n");
                printf("Initial number of particles = %d\n", particles.Nb_part_ini);
                printf("New number of particles     = %d\n", particles.Nb_part    );

                printf("** Time for CountPartCell = %lf sec\n", (double)((double)omp_get_wtime() - t_omp) );

                // Remove particles that would be above the surface
                if ( this->model.free_surf == 1 ) {
                    CleanUpSurfaceParticles( &particles, &mesh, topo, this->scaling );
                    CellFlagging( &mesh, this->model, topo, this->scaling );
                }

                //                for (int iphase=0; iphase<this->model.Nb_phases; iphase++) {
                //                CheckSym( mesh.phase_perc_n[iphase], 1.0, mesh.Nx-1, mesh.Nz-1, "perc_n" );
                //                CheckSym( mesh.phase_perc_s[iphase], 1.0, mesh.Nx-0, mesh.Nz-0, "perc_s" );
                //                }
            }
            this->model.dt = dt_solve;
        }

        // Update time
        this->model.time += this->model.dt;

        // Free solution arrays
        SFree( &Stokes );
        if ( IsNewtonStep == 1 ) SFree( &Jacob );
        if ( this->model.decoupled_solve == 1 ) {
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
        if ( writer == 1 && this->model.step % writer_step == 0 ) {

            printf("*************************************\n");
            printf("********* Write output files ********\n");
            printf("*************************************\n");

            // Breakpoint file
            t_omp = (double)omp_get_wtime();
            if ( this->model.delete_breakpoints == 1 ) DeletePreviousBreakpoint( this->model.step, writer_step  );
            MakeBreakpointParticles( &particles, &mesh, &topo_chain, &topo_chain_ini, this->model, &topo, &topo_ini, this->scaling );
            UpdateInputFile( fin_name, this->model.step);
            printf("** Time for Breakpoint file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));

            // Visualisation file
#ifndef _VG_
            t_omp = (double)omp_get_wtime();
            WriteOutputHDF5( &mesh, &particles, &topo, &topo_chain, this->model, Nmodel, "Output", this->materials, this->scaling );
            if ( this->model.write_markers == 1 ) WriteOutputHDF5Particles( &mesh, &particles, &topo, &topo_chain, &topo_ini, &topo_chain_ini, this->model, "Particles", this->materials, this->scaling );
            printf("** Time for Output file write = %lf sec\n", (double)((double)omp_get_wtime() - t_omp));
#endif
        }

        printf("** Total timestep calculation time = %lf sec\n", (double)((double)omp_get_wtime() - t_omp_step) );
        printf("** Model time = %2.2e sec\n", this->model.time*this->scaling.t );
        printf("** Current dt = %2.2e sec, Old dt = %2.2e sec\n", this->model.dt*this->scaling.t, this->model.dt0*this->scaling.t );

        // if (Np0 - particles.Nb_part != 0) exit (1);


    }

    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//
    //                                           END TIME LOOP : c'est fini les Doud'series...
    //------------------------------------------------------------------------------------------------------------------------------//
    //------------------------------------------------------------------------------------------------------------------------------//

    if ( this->model.GNUplot_residuals == 1) pclose(GNUplotPipe);

    //    // Free markers chains
    if ( this->model.free_surf == 1 )    FreeMarkerChain( &topo,     &topo_chain     );
    if ( this->model.free_surf == 1 )    FreeMarkerChain( &topo_ini, &topo_chain_ini );

    // Free Phase diagrams if activated
    if ( this->model.isPD == 1 ) FreePhaseDiagrams( &this->model );
    //
    // Free particle arrays
    PartFree( &particles, &this->model );

    // Free arrays
    GridFree( &mesh, &this->model );

    // Free char*'s
    free(this->model.input_file);
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
    for ( int k=0; k<this->model.Nb_phases; k++) {

        if ( this->materials.density_model[k] == 4 ) {
            printf("Unloading 1D diagram for coesite quartz only baby...\n");
            DoodzFree(this->model.PD1Drho[k]);
        }
    }
    DoodzFree(this->model.PD1DnP);
    DoodzFree(this->model.PD1Drho);
    DoodzFree(this->model.PD1Dmin);
    DoodzFree(this->model.PD1Dmax);
    if ( this->model.kinetics==1 ) DoodzFree(this->model.kin_dG);

    printf("\n********************************************************\n");
    printf("************* Ending MDOODZ 7.0 simulation *************\n");
    printf("********************************************************\n");
}

MdoodzInstance NewMdoodzInstance() {
  MdoodzInstance mdoodz = {
    .RunMDOODZ = RunMDOODZ,
    .crazyConductivity = NULL
  };
  return mdoodz;
}