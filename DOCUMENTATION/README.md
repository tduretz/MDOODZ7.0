# Parameters and parameter values read from input text file:

##  Simulation start/restart from breakpoint files
- `irestart`: 0 if simulation starts from beginning, 1 if simulation restarst from a step. Default: 00000
- `istep`: last step written, indicates step to restart from if irestart is 1. Defaut: 0

## Simulation start/restart 
- `writer`: Writes .hdf5 files for visualisation and .dat breakpoint files for restart. Default: 0
- `writer_step`: Frequency of output. Default 1.
- `writer_markers`: Writes hdf5 marker files (large files!). Default: 0
- `writer_debug`: Writes debug files. Default: 0
- `writer_subfolder`: Writes output in given subfolder. Default: ./
- `noisy`: Prints a lot of info to standard output. Default: 1
- `track_T_P_x_z`: Tracks initial T, P, x and z on particles 
- `delete_breakpoints`: Progressively deletes breakpoint files. Default: 0
- `gnuplot_log_res`: Activates GNU plot residuals visualisation (requires gnuplot). Default: 0

## External files
- `import_files_dir`: Default: `../../IMPORT`
- `import_file`: Default: `blah.bin`       
- `save_initial_markers`: saves initial marker configuration in a .bin file. Default: 0
- `load_initial_markers`: load the marker configurations if the .bin containing ,arker positions. Default: 0
- `initial_markers_file`: name of the file to be written/loaded. Default: `markers.bin`

## Physical scales
- `eta`: Viscosity [Pa.s]. Default: 1.0
- `L`: Length [m]. Default: 1.0
- `V`: Velocity [m/s]. Default: 1.0
- `T`: Temperature [K]. Default: 1.0
## Spatial domain
- `Nx`: Number of vertices in x direction. Default: 10
- `Nz`: Number of vertices in y direction. Default: 10
- `xmin`: Spatial domain extent. Default: -1.0
- `xmax`: Spatial domain extent
- `zmin`: Spatial domain extent. Default: -1.0 
- `zmax`: Spatial domain extent

## Time domain
- `Nt`:` Number of time steps. Default: 1  
- `dt`:` Time step [s]. Default: 0.0
- `Courant`: Courant number. Default: 0.5
- `RK`:` Order of Runge-Kutta advection solver (1, 2 or 4). Default 4.
- `constant_dt`:` Activates constant time step. Default: 0
- `stress_rotation`: 0: no stress rotation, 1: analytic rotation, 2: upper convected rate. Default: 1
- `dt_min`:` minimum allowed time step [s], defaut is negative such that it will never be activated unless specifically set in XXX.txt file. Default: -1e20
- `dt_max`:` maximum allowed time step [s], the default value is set to ~infinite, it we become effective only if specificaly set in XXX.txt (see e.g. LithoScale.txt). Default: 1e20

## Physics 
- `mechanical`: Activates mechanical solver. Default 1
- `advection`: Activates advection. Default 1
- `elastic`: Activates elasticity. Default 0
- `thermal`: Activates thermal solver. Default 0
- `anisotropy`: Turns on anisotropy. Default 0
- `polar`: Activate polar-Cartesian coordinates. Default 0
- `finite_strain`: Integrates finite stran and save deformation gradient tensor. Default 0
- `compressible`: Turns on compressible formulation. Default 0
- `density_variations`: Turns on volume change due to reaction if 1. Default 0
- `kinetics`: Activates reaction kinetics. Only for coesite --> quartz transformation so far. Default 0
- `out_of_plane: Out-of-plane strain. Default 0

## Numerics: linear solver
- `lin_solver`: 1: Powell-Hestenes, 2: Powell-Hestenes augmented (killer solver). Default: 2 
- `model.penalty`: Penalty factor. Default 1.0e3
- `auto_penalty`: Activates automatic penalty factor computation. Default: 0.0
- `diag_scaling`: Activates diagonal scaling. Default 1
- `preconditioner`: Preconditoner type for Newton.iterations, 0: Picard preconditionner. Default 0
- `lin_abs_div`: Tolerance for linear mechanical solver. Default 1e-9
- `lin_rel_div`: Tolerance for linear mechanical solver. Default 1e-5
- `lin_abs_mom`: Tolerance for linear mechanical solver. Default: 1e-9
- `lin_rel_mom`: Tolerance for linear mechanical solver. Default: 1e-5

## Numerics: non-linear solver
- `Newton`: Activates Newton iterations. 0: Picard, 1: Newton. Default: 0
- `nit_max`: Maximum number of iterations. Default 1
- `Picard2Newton`: Switch from Picard to Newton. iterations. Default: 0
- `Picard2Newton_tol`: Condition for switching based on residual magnitude. Default: 1e-1
- `max_Pic_its`: Condition for switching based on number of Picard iterations. Default: 10
- `line_search`: Activates line search. Default 0
- `line_search_min`: Minimum alpha value for line search. Default: 0.0
- `let_res_grow`: Allows residual to grow. Default 0
- `rel_tol_KSP`: Relative tolerance for inner Krylov solver. Default: 1e-4
- `nonlin_abs_mom`: Tolerance for non-linear mechanical solver. Default: 1e-6
- `nonlin_abs_div`: Tolerance for non-linear mechanical solver. Default: 1e-6
- `nonlin_rel_mom`: Tolerance for non-linear mechanical solver. Default: 1e-6
- `nonlin_rel_div`: Tolerance for non-linear mechanical solver. Default: 1e-6
- `min_eta`: Minimum viscosity. 1e18
- `max_eta`: Maximum viscosity. 1e24
- `safe_mode`: Activates safe mode: reduces time step if convergence fails. Default: 0
- `safe_dt_div`: Reduction factor for time step reduction. Default: 5.0
- `max_num_stag`: Maximum number of stagnation (safe mode). Default 3
- `residual_form`: Form of residual - TODO: delete if our models work with new default value (1). Default 1
```C    
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
