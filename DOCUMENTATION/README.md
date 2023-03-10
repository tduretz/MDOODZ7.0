```C
// Simulation start/restart from Breakpoint
    model.istep              = ReadInt2( fin, "istep", 0 );
    model.irestart           = ReadInt2( fin, "irestart", 0 );
    if ( model.istep == 0 ) model.irestart = 0; // Override the restart step number written in the text file (istep)

    // Output
    model.writer             = ReadInt2( fin, "writer",              0 ); // Write files
    model.writer_step        = ReadInt2( fin, "writer_step",         1 ); // Frequency of output
    model.writer_markers     = ReadInt2( fin, "writer_markers",      0 ); // Writes marker files
    model.writer_debug       = ReadInt2( fin, "writer_debug",        0 ); // Writes debug files
    model.writer_subfolder   = ReadChar( fin, "writer_subfolder",  "./"); // Writes output in given subfolder
    // printf("%s\n", model.writer_subfolder); exit(1);
    model.noisy              = ReadInt2( fin, "noisy",               1 ); // Prints a lot of info to standard output
    model.track_T_P_x_z      = ReadInt2( fin, "track_T_P_x_z",       0 ); // Tracks initial T, P, x and z on particles 
    model.delete_breakpoints = ReadInt2( fin, "delete_breakpoints",  1 ); // Progressively deletes breakpoint files
    model.gnuplot_log_res    = ReadInt2( fin, "gnuplot_log_res",     0 ); // Activates GNU plot residuals visualisation

    // Input
    model.import_files_dir     = ReadChar( fin, "import_files_dir",    "../../IMPORT");
    model.import_file          = ReadChar( fin, "import_file",             "blah.bin");
    model.save_initial_markers = ReadInt2( fin, "save_initial_markers",            0 );
    model.load_initial_markers = ReadInt2( fin, "load_initial_markers",            0 );
    model.initial_markers_file = ReadChar( fin, "initial_markers_file", "markers.bin");

    // Read scales for non-dimensionalisation
    scale scaling            = (scale){
                        .eta = ReadDou2(fin, "eta", 1.0),
                        .L   = ReadDou2(fin, "L", 1.0),
                        .V   = ReadDou2(fin, "V", 1.0),
                        .T   = ReadDou2(fin, "T", 1.0),
    };
    ScaleMe( &scaling );
    double Ga = 1e9*365.25*3600*24/scaling.t;

    // Spatial domain
    model.Nx                 = ReadInt2( fin, "Nx",        10 );            // Number of vertices in x direction
    model.Nz                 = ReadInt2( fin, "Nz",        10 );            // Number of vertices in y direction
    model.xmin               = ReadDou2( fin, "xmin",     1.0 )/scaling.L;  // Spatial domain extent
    model.zmin               = ReadDou2( fin, "zmin",     1.0 )/scaling.L;  // Spatial domain extent 
    model.xmax               = ReadDou2( fin, "xmax",     1.0 )/scaling.L;  // Spatial domain extent 
    model.zmax               = ReadDou2( fin, "zmax",     1.0 )/scaling.L;  // Spatial domain extent
    // Time domain
    model.Nt                 = ReadInt2( fin, "Nt",         1 );            // Number of time steps    
    model.dt                 = ReadDou2( fin, "dt",       0.0 ) /scaling.t; // Time step
    model.Courant            = ReadDou2( fin, "Courant",             0.5 ); // Courant number
    model.RK                 = ReadInt2( fin, "RK",                    4 ); // Order of Runge-Kutta advection solver (1, 2 or 4)
    model.constant_dt        = ReadInt2( fin, "constant_dt",           0 ); // Activates constant time step
    model.stress_rotation    = ReadInt2( fin, "stress_rotation",       1 ); // 0: no stress rotation, 1: analytic rotation, 2: upper convected rate
    model.dt_max             = ReadDou2( fin, "dt_max", 1e20 ) /scaling.t;  // maximum allowed time step, the default value is set to ~infinite, it we become effective only if specificaly set in XXX.txt (see e.g. LithoScale.txt)
    model.dt_min             = ReadDou2( fin, "dt_min",-1e20 ) /scaling.t;  // minimum allowed time step, defaut is negative such that it will never be activated unless specifically set in XXX.txt file
    // Physics 
    model.mechanical         = ReadInt2( fin, "mechanical",            1 ); // Activates mechanical solver
    model.advection          = ReadInt2( fin, "advection",             1 ); // Activates advection
    model.elastic            = ReadInt2( fin, "elastic",               0 ); // Activates elasticity
    model.thermal            = ReadInt2( fin, "thermal",               0 ); // Activates thermal solver
    model.anisotropy         = ReadInt2( fin, "anisotropy",            0 ); // Turns on anisotropy
    model.polar              = ReadInt2( fin, "polar",                 0 ); // Activate polar-Cartesian coordinates
    model.finite_strain      = ReadInt2( fin, "finite_strain",         0 ); // Integrates finite stran and save deformation gradient tensor
    model.compressible       = ReadInt2( fin, "compressible",          0 ); // Turns on compressibility
    model.density_variations = ReadInt2( fin, "density_variations",    0 ); // Turns on volume change due to reaction if 1
    model.kinetics           = ReadInt2( fin, "kinetics",              0 ); // Activates reaction kinetics
    model.out_of_plane       = ReadInt2( fin, "out_of_plane",          0 ); // Out-of-plane strain
    // Numerics: linear solver
    model.penalty            = ReadDou2( fin, "penalty",          1.0e10 ); // Penalty factor
    model.auto_penalty       = ReadDou2( fin, "auto_penalty",        0.0 ); // Activates automatic penalty factor computation
    model.diag_scaling       = ReadInt2( fin, "diag_scaling",          1 ); // Activates diagonal scaling
    model.preconditioner     = ReadInt2( fin, "preconditioner",        0 ); // Preconditoner type for Newton ietrations, 0: Picard preconditionner
    model.lin_abs_div        = ReadDou2( fin, "lin_abs_div",      1.0e-9 ); // Tolerance for linear mechanical solver
    model.lin_rel_div        = ReadDou2( fin, "lin_rel_div",      1.0e-5 ); // Tolerance for linear mechanical solver
    model.lin_abs_mom        = ReadDou2( fin, "lin_abs_mom",      1.0e-9 ); // Tolerance for linear mechanical solver
    model.lin_rel_mom        = ReadDou2( fin, "lin_rel_mom",      1.0e-5 ); // Tolerance for linear mechanical solver
    model.lin_solver         = ReadInt2( fin, "lin_solver",            2 ); // 1: Powell-Hestenes, 2: Powell-Hestenes augmented (killer solver) 
    // Numerics: non-linear solver
    Nmodel.nit_max           = ReadInt2( fin, "nit_max",               1 ); // Maximum number of iterations
    model.Newton             = ReadInt2( fin, "Newton",                0 ); // Activates Newton iterations
    Nmodel.Picard2Newton     = ReadInt2( fin, "Picard2Newton",         0 ); // Switch from Picard to Newton iterations
    Nmodel.Picard2Newton_tol = ReadDou2( fin, "Picard2Newton_tol",  1e-1 ); // Condition for switching based on residual magnitude
    Nmodel.max_Pic_its       = ReadInt2( fin, "max_Pic_its",          10 ); // Condition for switching based on number of Picard iterations
    Nmodel.let_res_grow      = ReadInt2( fin, "let_res_grow",          0 ); // Allows residual to grow 
    model.rel_tol_KSP        = ReadDou2( fin, "rel_tol_KSP",        1e-4 ); // Relative tolerance for inner Krylov solver
    Nmodel.nonlin_abs_mom    = ReadDou2( fin, "nonlin_abs_mom",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    Nmodel.nonlin_abs_div    = ReadDou2( fin, "nonlin_abs_div",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    Nmodel.nonlin_rel_mom    = ReadDou2( fin, "nonlin_rel_mom",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    Nmodel.nonlin_rel_div    = ReadDou2( fin, "nonlin_rel_div",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    model.min_eta            = ReadDou2( fin, "min_eta", 1e18)/scaling.eta; // Minimum viscosity
    model.max_eta            = ReadDou2( fin, "max_eta", 1e24)/scaling.eta; // Maximum viscosity
    model.safe_mode          = ReadInt2( fin, "safe_mode",             0 ); // Activates safe mode: reduces time step if convergence fails
    model.safe_dt_div        = ReadDou2( fin, "safe_dt_div",         5.0 ); // Reduction factor for time step reduction
    model.max_num_stag       = ReadInt2( fin, "max_num_stag",          3 ); // maximum number of stagnation (safe mode)
    model.line_search        = ReadInt2( fin, "line_search",           0 ); // Activates line search
    model.line_search_min    = ReadDou2( fin, "line_search_min",     0.0 ); // Minimum alpha value for line search 
    model.residual_form      = ReadInt2( fin, "residual_form",         1 ); // Form of residual - TODO: delete if our models work with new default value (1)
    Nmodel.stagnated         = 0;
    // Numerics: marker-in-cell
    ParticlesInput particles;
    model.eta_average        = ReadInt2( fin, "eta_average",           0 ); // 0: arithmetic mean - 1: harmonic mean - 2: geometric mean
    model.interp_stencil     = ReadInt2( fin, "interp_stencil",        1 ); // 1: 1-Cell          - 9: 9-Cell
    model.subgrid_diffusion  = ReadInt2( fin, "subgrid_diffusion",     0 ); // 0: No subgrid diffusion, 1: temperature, 2: temperature + stress
    model.conserv_interp     = ReadInt2( fin, "conserv_interp",        0 ); // Activates Taras conservative interpolation
    model.direct_neighbour   = ReadInt2( fin, "direct_neighbour",      0 ); // Direct neighbour interpolation
    model.initial_noise      = ReadInt2( fin, "initial_noise",         0 ); // Add noise on initial marker locations
    model.marker_noise       = ReadInt2( fin, "marker_noise",          0 ); // Background noise field generated and tracked on the particles 
    model.reseed_markers     = ReadInt2( fin, "reseed_markers",        1 ); // Activates reseeding / particle injection
    particles.Nx_part        = ReadInt2( fin, "Nx_part",               4 ); // number of particle per cell in x
    particles.Nz_part        = ReadInt2( fin, "Nz_part",               4 ); // number of particle per cell in y
    particles.min_part_cell  = ReadInt2( fin, "min_part_cell",        16 ); // minimum number of particle per cell (if below: will trigger reseeding)
    particles.Nb_part        = (model.Nx-1)*(model.Nz-1) * particles.Nx_part * particles.Nz_part;
    particles.Nb_part_max    = 4.1*particles.Nb_part;
    // Boundary conditions
    model.shear_style        = ReadInt2( fin, "shear_style",           0 ); // BC type: 0: pure shear, 2: periodic simple shear
    model.periodic_x         = ReadInt2( fin, "periodic_x",            0 ); // Activates periodicity in x
    model.pure_shear_ALE     = ReadInt2( fin, "pure_shear_ALE",        0 ); // Activates Arbitrary Lagarangian Eulerian mode (pure shear box deformation)
    model.free_surface       = ReadInt2( fin, "free_surface",          0 ); // Activates free surface
    model.free_surface_stab  = ReadDou2( fin, "free_surface_stab",   0.0 ); // Activate free surface stabilisation: range 0.0-2.0
    // Model configurations
    model.initial_cooling    = ReadInt2( fin, "initial_cooling",       0 ); // Activates initial cooling
    model.cooling_duration   = ReadDou2( fin, "cooling_duration",     Ga ); // Initial cooling duration
    model.shear_heating      = ReadInt2( fin, "shear_heating",         1 ); // Activates shear heating
    model.adiab_heating      = ReadInt2( fin, "adiab_heating",         0 ); // 0: zero, 1: lithostatic P assumption, 2: full derivative
    model.surface_processes  = ReadInt2( fin, "surface_processes",     0 ); // 1: diffusion; 2: diffusion + sedimentation
    model.marker_aniso_angle = ReadInt2( fin, "marker_aniso_angle",    0 ); // Enables setting anisotropy angle per particles rather than phases
    // Transformations
    model.progress_transform = ReadInt2( fin, "progress_transform",    0 ); // Activate progressive reactions
    model.no_return          = ReadInt2( fin, "no_return",             0 ); // Turns off retrogression if 1.0
    model.unsplit_diff_reac  = ReadInt2( fin, "unsplit_diff_reac",     0 ); // Unsplits diffusion and reaction
    model.smooth_softening   = ReadInt2( fin, "smooth_softening",      1 ); // Activates smooth explicit kinematic softening function
    // Background ambient conditions
    model.bkg_strain_rate    = ReadDou2( fin, "bkg_strain_rate", 1e-30)/scaling.E; // Background tectonic rate, defaut is close to zero to avoid any Nans of Infs in rheology
    model.bkg_div_rate       = ReadDou2( fin, "bkg_div_rate",      0.0)/scaling.E; // Background divergence rate
    model.bkg_pressure       = ReadDou2( fin, "bkg_pressure",      0.0)/scaling.S; // Background pressure
    model.bkg_temperature    = ReadDou2( fin, "bkg_temperature",   0.0)/scaling.T; // Background temperature
    // Surface processes
    model.surf_diff          = ReadDou2( fin, "surf_diff",       0.0 ) / (pow(scaling.L,2.0)/scaling.t);
    model.surf_ised1         = ReadInt2( fin, "surf_ised1",      0.0 );
    model.surf_ised2         = ReadInt2( fin, "surf_ised2",      0.0 );
    model.surf_sedirate      = ReadDou2( fin, "surf_sedirate",   0.0 ) / scaling.V;
    model.surf_baselev       = ReadDou2( fin, "surf_baselev",    0.0 ) / scaling.L;
    model.surf_Winc          = ReadDou2( fin, "surf_Winc",       0.0 ) / scaling.L;
    model.surf_Vinc          = ReadDou2( fin, "surf_Vinc",       0.0 ) / scaling.V;
    // Initial thermal perturbation
    model.therm_perturb      = ReadInt2( fin, "therm_perturb",                 0 ); // Includes initial thermal perbation
    model.therm_perturb_x0   = ReadDou2( fin, "therm_perturb_x0",  0.0 )/scaling.L; // x position
    model.therm_perturb_z0   = ReadDou2( fin, "therm_perturb_z0",  0.0 )/scaling.L; // y position
    model.therm_perturb_rad  = ReadDou2( fin, "therm_perturb_rad", 0.0 )/scaling.L; // Radius
    model.therm_perturb_dT   = ReadDou2( fin, "therm_perturb_dT" , 0.0 )/scaling.T; // Temperature anomaly
    // For rheological database reasons...
    model.force_act_vol_ast  = ReadInt2( fin, "force_act_vol_ast",   0 ); // if 1 then:
    model.act_vol_dis_ast    = ReadDou2( fin, "act_vol_dis_ast" ,  0.0 ); // ... set dislocation creep to value
    model.act_vol_dif_ast    = ReadDou2( fin, "act_vol_dif_ast" ,  0.0 ); // ... set diffusion creep to value
    // Model user's delights
    model.user0              = ReadDou2( fin, "user0",           0.0 );
    model.user1              = ReadDou2( fin, "user1",           0.0 );
    model.user2              = ReadDou2( fin, "user2",           0.0 );
    model.user3              = ReadDou2( fin, "user3",           0.0 );
    model.user4              = ReadDou2( fin, "user4",           0.0 );
    model.user5              = ReadDou2( fin, "user5",           0.0 );
    model.user6              = ReadDou2( fin, "user6",           0.0 );
    model.user7              = ReadDou2( fin, "user7",           0.0 );
    model.user8              = ReadDou2( fin, "user8",           0.0 );

```
