#ifndef MDOODZ_H
#define MDOODZ_H

#define zeroC     273.15
#define Rg        8.314510
#define PI        3.14159265359
#define Rad_Earth 6370000
#include "stdbool.h"

#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */

// address Williams comments
typedef enum { ARITHMETIC = 0,
               HARMONIC   = 1,
               GEOMETRIC  = 2 } average;

// Tensor 2D
typedef struct {
   double xx, zz, yy, xz, zx, ii, ii2;
} Tensor2D;

// Vector 2D
typedef struct {
   double x, z;
} Vector2D;

// params contains the model parameters
typedef struct {
  int    balance_boundaries, zero_mean_topo; 
  char   description[500];
  double xmin, zmin, xmax, zmax, time, dx, dz, dt, dt0, dt_start, dt_max, L0,
          dt_min;
  double  xmin0, zmin0, xmax0, zmax0;
  double  gx, gz;
  int     Nx, Nz, Nt, step, nit, Newton, noisy;
  average eta_average, ani_average;
  int     interp_stencil;
  double  nexp_radial_basis;
  int     mechanical, periodic_x, elastic, isnonnewtonian,
          thermal, pure_shear_ALE, free_surface, writer_markers, writer_debug, topo_update;
  double free_surface_stab;
  int    constant_dt, RK, line_search, initial_cooling, subgrid_diffusion, adiab_heating,
          shear_heating, advection, finite_strain, conserv_interp;
  int surface_processes, loc_iter, therm_perturb, surf_ised1,
          surf_ised2, MantleID, topografix, reseed_markers, smooth_softening, fix_temperature;
  double bkg_strain_rate, bkg_div_rate, user0, user1, user2, user3, user4, user5, user6, user7,
          user8;
  char  *import_file;
  char  *import_files_dir;
  int    Nb_phases;
  int    ncont;
  double Courant, min_eta, max_eta, eta_tol;
  // Particles
  int    initial_noise;
  // Linear solver
  int    lin_solver, diag_scaling, preconditioner;
  double penalty, lin_abs_div, lin_rel_div, lin_abs_mom, lin_rel_mom, auto_penalty, compressible,
          rel_tol_KSP;
  // Non-linear solver
  double line_search_min, safe_dt_div;
  int    safe_mode, max_num_stag;
  // Deformation maps
  int    nT, nE, nd, deformation_maps;
  double Pn, Tmin, Tmax, Emin, Emax, dmin, dmax, bkg_pressure, bkg_temperature;
  // Surface processes
  double surf_diff, surf_sedirate, surf_baselev, surf_Winc, surf_Vinc;
  // Initial thermal perturbation
  double therm_perturb_x0, therm_perturb_z0, therm_perturb_dT, therm_perturb_rad,
          cooling_duration;
  // For rheological database...
  int      force_act_vol_ast;
  double   act_vol_dis_ast, act_vol_dif_ast;
  // Phase diagrams
  int      isPD, num_PD, *PDMnT, *PDMnP, *PD1DnP;
  double **PDMrho, *PDMTmin, *PDMTmax, *PDMPmin, *PDMPmax;
  double **PD1Drho, *PD1Dmin, *PD1Dmax;
  // Kinetics
  int      kin_nP, kin_nT;
  double  *kin_dG, kin_Tmin, kin_Tmax, kin_Pmin, kin_Pmax;
  // Visualisation
  int      track_T_P_x_z, delete_breakpoints, gnuplot_log_res;
  // Boundary conditions type
  int      BC_setup_type, shear_style, polar;
  int      stress_rotation, direct_neighbour;
  // For diffused rheological constrasts
  int      diffuse_X, diffuse_avg;
  double   diffusion_length;
  // For Pips
  int      chemical_diffusion, no_return, density_variations, unsplit_diff_reac, kinetics;
  // Anisotropy
  int      anisotropy, out_of_plane, marker_noise; //aniso_fstrain
  int      residual_form;
  int      irestart, istep;
  int      writer, writer_step;
  const char     *writer_subfolder;
  int      save_initial_markers, load_initial_markers;
  char    *initial_markers_file;
  int     marker_aniso_angle;
} params;

// Stucture scale contains scaling parameters
typedef struct {
  double eta, L, V, T, t, a, E, S, m, rho, F, J, W, Cv, rhoE, k;
} scale;

// mat_prop contains information related to material phase properties
typedef struct {
  int    Nb_phases;
  double R;
  double  eps0[20], tau0[20], eta0[20], rho[20], G[20], Cv[20], k[20], Qr[20], C[20], phi[20],
          psi[20], Slim[20], n[20], A[20], Ea[20], Va[20], alp[20], bet[20], Qm[20],
          T0[20], P0[20], drho[20], k_eff[20];
  double tpwl[20], Qpwl[20], Vpwl[20], npwl[20], mpwl[20], Apwl[20], apwl[20],
          fpwl[20], rpwl[20], Fpwl[20], pref_pwl[20];
  double texp[20], Qexp[20], Vexp[20], Sexp[20], Eexp[20], Gexp[20], aexp[20],
          fexp[20], rexp[20], qexp[20], nexp[20];
  double tlin[20], Qlin[20], Vlin[20], nlin[20], mlin[20], Alin[20], alin[20],
          flin[20], rlin[20], Flin[20];
  double tgbs[20], Qgbs[20], Vgbs[20], ngbs[20], mgbs[20], Agbs[20], agbs[20],
          fgbs[20], rgbs[20], Fgbs[20];
  double ppzm[20], Kpzm[20], Qpzm[20], Vpzm[20], Gpzm[20], cpzm[20], Lpzm[20],
          gs_ref[20];
  double Skin[20], kkin[20], Qkin[20];
  int    kin[20], gs[20], cstv[20], pwlv[20], linv[20], expv[20], gbsv[20],
          phase_diagram[20], density_model[20];
  double C_end[20], phi_end[20], psi_end[20], pls_start[20], pls_end[20],
          eta_vp[20], n_vp[20], sig_tens[20], sig1[20], dsig1[20];
  int    plast[20], phi_soft[20], psi_soft[20], coh_soft[20], yield[20];
  double Pr[20], tau_kin[20], dPr[20], k_chem[20];
  int    reac_soft[20], reac_phase[20];
  int    phase_mix[20], phase_two[20];
  double aniso_angle[20], ani_fac_max[20], aniso_factor[20];   //ani_fac_v[20], ani_fac_e[20], ani_fac_p[20]
  double axx[20], azz[20], ayy[20];
  int    ani_fstrain[20];
  int    transmutation[20], transmutation_phase[20];
  double transmutation_temperature[20];
} mat_prop;

char                         *GetSetupFileName(int nargs, char *args[]);

typedef struct MdoodzSetup    MdoodzSetup;
typedef struct MdoodzInput    MdoodzInput;

typedef struct {
  double l;
  double k;
  double x;
  double z;
} Coordinates;

typedef double (*SetSurfaceZCoord_f)(MdoodzInput *input, double x_coord);
typedef int (*SetSurfacePhase_f)(MdoodzInput *input, double x_coord);

typedef struct {
  SetSurfaceZCoord_f SetSurfaceZCoord;
  SetSurfacePhase_f  SetSurfacePhase;
} BuildInitialTopography_ff;

typedef double (*SetHorizontalVelocity_f)(MdoodzInput *input, Coordinates coordinates);
typedef double (*SetVerticalVelocity_f)(MdoodzInput *input, Coordinates coordinates);
typedef double (*SetTemperature_f)(MdoodzInput *input, Coordinates coordinates);
typedef int    (*SetPhase_f)(MdoodzInput *input, Coordinates coordinates);
typedef int    (*SetDualPhase_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetGrainSize_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetPorosity_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetDensity_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetXComponent_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetPressure_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetNoise_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef double (*SetAnisoAngle_f)(MdoodzInput *input, Coordinates coordinates, int phase);
typedef Tensor2D (*SetDefGrad_f)(MdoodzInput *input, Coordinates coordinates, int phase);

typedef struct {
  SetHorizontalVelocity_f SetHorizontalVelocity;
  SetVerticalVelocity_f   SetVerticalVelocity;
  SetPhase_f              SetPhase;
  SetDualPhase_f          SetDualPhase;
  SetPressure_f           SetPressure;
  SetNoise_f              SetNoise;
  SetTemperature_f        SetTemperature;
  SetGrainSize_f          SetGrainSize;
  SetPorosity_f           SetPorosity;
  SetDensity_f            SetDensity;
  SetXComponent_f         SetXComponent;
  SetDefGrad_f            SetDefGrad;
  SetAnisoAngle_f         SetAnisoAngle;
} SetParticles_ff;

typedef enum {
  NE,
  NW,
  SE,
  SW,
  INTERNAL,
  N,
  S,
  W,
  E,
  free_surface
} POSITION;

// Thermal boundary condition types
typedef enum {
    constant_heatflux      = 0,
    constant_temperature   = 1,
    periodic_temperature   =-2,
} Thermal_BC;

// Chemical boundary condition types
typedef enum {
    constant_Xflux = 0,
    constant_X     = 1,
    periodic_X     =-2,
} Chemical_BC;

// Mechanical boundary condition types
typedef enum {
    inside                =  -1,
    constant_velocity     =   0,
    constant_shear_stress =  13,
    constant_velocity_NC  =  11,
    periodic_VxW          =  -2,
    periodic_VxE          = -12,
    periodic_VzW          = -12,
    periodic_VzE          = -12,
} Mechanical_BC;

typedef struct {
  double value;
  char   type;
} SetBC;

typedef SetBC (*SetBCVx_f)(MdoodzInput *input, POSITION position, Coordinates coordinates);
typedef SetBC (*SetBCVz_f)(MdoodzInput *input, POSITION position, Coordinates coordinates);
typedef SetBC (*SetBCT_f)(MdoodzInput *input, POSITION position, double gridTemperature);
typedef SetBC (*SetBCC_f)(MdoodzInput *input, POSITION position, double gridXvalue);
typedef double (*FixTemperature_f)(MdoodzInput *input, double pressure);
typedef char (*SetBCPType_f)(MdoodzInput *input, POSITION position);

typedef struct {
  SetBCVx_f    SetBCVx;
  SetBCVz_f    SetBCVz;
  SetBCT_f     SetBCT;
  SetBCC_f     SetBCC;
  SetBCPType_f SetBCPType;
  FixTemperature_f  FixTemperature;
} SetBCs_ff;

typedef struct {
  double multiplier;
  int   *phases;
  int    nPhases;
} CrazyConductivity;

typedef struct {
  int   nx;
  int   nz;
  int   nb_elems;
  char *ph_hr;
} Geometry;

typedef struct {
  double double1;
  double double2;
  double double3;
  double double4;
  int int1;
  int int2;
  int int3;
  int int4;
  const char *str1;
} MutateInputParams;

typedef void(MutateInput_f)(MdoodzInput *input, MutateInputParams *mutateInputParams);

struct MdoodzSetup {
  BuildInitialTopography_ff *BuildInitialTopography;
  SetParticles_ff           *SetParticles;
  SetBCs_ff                 *SetBCs;
  MutateInput_f             *MutateInput;
  MutateInputParams         *mutateInputParams;
};

typedef struct {
  double east;
  double west;
} LateralFlux;

struct MdoodzInput {
  char              *inputFileName;
  params             model;
  mat_prop           materials;
  scale              scaling;
  LateralFlux        *flux;
  CrazyConductivity *crazyConductivity;
  Geometry          *geometry;
};

void  RunMDOODZ(char *inputFileName, MdoodzSetup *setup);

// Setup templates

SetBC SetPureShearBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates);
SetBC SetPureShearBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates);
SetBC SetSimpleShearBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates);
SetBC SetSimpleShearBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates);
SetBC SetPureOrSimpleShearBCVx(MdoodzInput *input, POSITION position, Coordinates coordinates);
SetBC SetPureOrSimpleShearBCVz(MdoodzInput *input, POSITION position, Coordinates coordinates);

typedef struct {
  double centreX;
  double centreZ;
  double radiusX;
  double radiusZ;
  double angle;
} Ellipse;

typedef struct {
  double centreX;
  double centreZ;
  double sizeX;
  double sizeZ;
  double angle;
} Rectangle;

bool IsEllipseCoordinates(Coordinates coordinates, Ellipse ellipse, double scalingL);
bool IsRectangleCoordinates(Coordinates coordinates, Rectangle rectangle, double scalingL);

#endif
