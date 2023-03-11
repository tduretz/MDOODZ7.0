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
- `out_of_plane`: Out-of-plane strain. Default 0

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
   
## Numerics: marker-in-cell
- `eta_average`: 0: arithmetic mean - 1: harmonic mean - 2: geometric mean. Default: 0
- `interp_stencil`: 1: 1-Cell          - 9: 9-Cell. Default: 1
- `subgrid_diffusion`: 0: No subgrid diffusion, 1: temperature, 2: temperature + stress. Default: 0
- `conserv_interp`: Activates Taras. Default: 0 conservative interpolation
- `direct_neighbour`: Direct neighbour. Default: 0 interpolation
- `initial_noise`: Add noise on initial. Default: 0 marker locations
- `marker_noise`: Background noise. Default: 0field generated and tracked on the particles. Default: 0
- `reseed_markers`: Activates reseeding / particle injection. Default: 1
- `Nx_part"`: Number of particle per cell in x. Default: 4
- `Nz_part"`: Number of particle per cell in y. Default: 4
- `min_part_cell`: Minimum number of particle per cell (if below: will trigger reseeding). Default: 16

## Boundary conditions
- `shear_style`: BC type: 0: pure shear, 2: periodic simple shear. Default: 0
- `periodic_x`: Activates periodicity in x. Default: 0
- `pure_shear_ALE`: Activates Arbitrary Lagarangian Eulerian mode (pure shear box deformation). Default: 0
- `free_surface`: Activates free surface. Default: 0
- `free_surface_stab`: Activate free surface stabilisation: range 0.0-2.0. Default: 0.0 

## Model configurations
- `initial_cooling`:  Activates initial cooling. Default 0
- `cooling_duration`:  Initial cooling duration. Default 1 Ga
- `shear_heating`:  Activates shear heating. Default 1
- `adiab_heating`:  0: zero, 1: lithostatic P assumption, 2: full derivative. Default 0
- `surface_processes`:  1: diffusion; 2: diffusion + sedimentation; 5: diffusion + localised incision. Default 0
- `marker_aniso_angle`:  Enables setting anisotropy angle per particles rather than phases. Default 0

## Transformations
- `density_variations`: Turns on volume change due to reaction if 1. Default 0
- `kinetics`: Activates reaction kinetics. Only for coesite --> quartz transformation so far. Default 0
- `progress_transform`: Activate progressive reactions. Default 0
- `no_return`: Turns off retrogression if 1 Default 0
- `unsplit_diff_reac`: Unsplits diffusion and reaction Default. 0
- `smooth_softening`: Activates smooth explicit kinematic softening function. Default 1

## Background ambient conditions
- `bkg_strain_rate`: Background tectonic rate, default is close to zero to avoid any Nans of Infs in rheology. Default 1e-30 [1/s]
- `bkg_div_rate`:   Background divergence rate. Default 0
- `bkg_pressure`:   Background pressure. Default 0
- `bkg_temperature`: Background temperature. Default 0

## Surface processes
- `surf_diff`: Surface diffusivity [m^2/s]. Default: 0
- `surf_ised1`: Marker sediment index 1. Default: 0
- `surf_ised2`: Marker sediment index 2. Default: 0
- `surf_sedirate`: Sedimentation rate [m/s]. Default: 0.0
- `surf_baselev`: Sedimentation base level [m]. Default: 0.0
- `surf_Winc`: Incision width [m] for `surface processes == 5`. Default: 0.0
- `surf_Vinc`: Incision rate [m] for `surface processes == 5`. Default: 0.0

## Initial thermal perturbation
- `therm_perturb`: Includes an initial circular thermal perturbation. Default: 0
- `therm_perturb_x0`: x center position [m]. Default: 0.0
- `therm_perturb_z0`: z center position [m]. Default: 0.0
- `therm_perturb_rad`: Radius [m]. Default: 0.0
- `therm_perturb_dT`: Temperature anomaly [K]. Default: 0.0
## For rheological database purposes...
- `force_act_vol_ast`: Default is 0. If 1 then:
- `act_vol_dis_ast`:  ... set dislocation creep to value [J/bar]
- `act_vol_dif_ast`:  ... set diffusion creep to value [J/bar]
## Model user's delights
- `user0`: Can be used to parse anything. Needs to be scaled internally 
- ...
- `user8`: Can be used to parse anything. Needs to be scaled internally 
