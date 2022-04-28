#ifndef MDOODZ_H
#define MDOODZ_H

#define zeroC     273.15
#define Rg        8.314510
#define PI        3.14159265359
#define Rad_Earth 6370000

// address Williams comments
typedef enum { ARITHMETIC = 0,
               HARMONIC   = 1,
               GEOMETRIC  = 2 } ETA_AVG;

// params contains the model parameters
typedef struct {
  double xmin, zmin, xmax, zmax, time, dx, dz, dt, dt0, dt_start, dt_max, L0,
          dt_min;
  double  xmin0, zmin0, xmax0, zmax0;
  double  gx, gz;
  int     Nx, Nz, Nt, step, nit, Newton, noisy;
  ETA_AVG eta_avg;
  int     itp_stencil;
  double  nexp_radial_basis;
  int     ismechanical, isperiodic_x, isinertial, iselastic, isnonnewtonian,
          isthermal, ispureshear_ale, free_surf, write_markers, write_debug;
  double free_surf_stab;
  int    dt_constant, RK, line_search, thermal_eq, subgrid_diff, adiab_heat,
          shear_heat, advection, fstrain, ConservInterp;
  int surf_processes, cpc, surf_remesh, loc_iter, therm_pert, surf_ised1,
          surf_ised2, MantleID, topografix, Reseed, SmoothSoftening;
  double EpsBG, DivBG, user0, user1, user2, user3, user4, user5, user6, user7,
          user8;
  char  *input_file;
  int    Nb_phases;
  int    ncont;
  double Courant, mineta, maxeta;
  // Particles
  int initial_noise, initial_part;
  // Linear solver
  int    decoupled_solve, lsolver, diag_scaling, pc_type;
  double penalty, abs_tol_div, rel_tol_div, auto_penalty, compressible,
          rel_tol_KSP;
  // Non-linear solver
  double line_search_min, safe_dt_div;
  int    safe_mode, nstagmax;
  // Deformation maps
  int    nT, nE, nd, def_maps;
  double Pn, Tmin, Tmax, Emin, Emax, dmin, dmax, PrBG, TBG;
  // Surface processes
  double surf_diff, surf_sedirate, surf_baselev, surf_Winc, surf_Vinc;
  // Initial thermal perturbation
  double therm_pert_x0, therm_pert_z0, therm_pert_dT, therm_pert_rad,
          cooling_time;
  // For rheological database...
  int    force_act_vol_ast;
  double act_vol_dis_ast, act_vol_dif_ast;
  // Phase diagrams
  int      isPD, num_PD, *PDMnT, *PDMnP, *PD1DnP;
  double **PDMrho, *PDMTmin, *PDMTmax, *PDMPmin, *PDMPmax;
  double **PD1Drho, *PD1Dmin, *PD1Dmax;
  // Kinetics
  int     kin_nP, kin_nT;
  double *kin_dG, kin_Tmin, kin_Tmax, kin_Pmin, kin_Pmax;
  // Visualisation
  int rec_T_P_x_z, delete_breakpoints, GNUplot_residuals;
  // Boundary conditions type
  int BC_setup_type, shear_style, polar;
  int StressRotation, StressUpdate, DirectNeighbour;
  // For diffused rheological constrasts
  int    diffuse_X, diffuse_avg;
  double diffusion_length;
  // For Pips
  int ProgReac, NoReturn, VolChangeReac, Plith_trick, UnsplitDiffReac, kinetics;
  // Anisotropy
  int aniso, aniso_fstrain, oop, noise_bg;
  int eqn_state;
  int residual_form;
} params;

// Stucture scale contains scaling parameters
typedef struct {
  double eta, L, V, T, t, a, E, S, m, rho, F, J, W, Cv, rhoE, k;
} scale;

// mat_prop contains information related to material phase properties
typedef struct {
  int    Nb_phases;
  double R;
  double eta0[20], rho[20], mu[20], Cv[20], k[20], Qr[20], C[20], phi[20],
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
          eta_vp[20], n_vp[20];
  int    plast[20], phi_soft[20], psi_soft[20], coh_soft[20], is_tensile[20];
  double Pr[20], tau_kin[20], dPr[20], k_chem[20];
  int    reac_soft[20], reac_phase[20];
  int    phase_mix[20], phase_two[20];
  double aniso_factor[20], aniso_angle[20];
} mat_prop;

char *GetSetupFileName(int nargs, char *args[]);

typedef struct MdoodzInstance MdoodzInstance;

typedef struct {
  double x;
  double z;
} Coordinates;

typedef double (*SetSurfaceZCoord_f)(MdoodzInstance *instance, double x_coord);
typedef int (*SetSurfacePhase_f)(MdoodzInstance *instance, double x_coord);

typedef struct {
  SetSurfaceZCoord_f SetSurfaceZCoord;
  SetSurfacePhase_f  SetSurfacePhase;
} BuildInitialTopography_ff;

typedef double (*SetHorizontalVelocity_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef double (*SetVerticalVelocity_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef double (*SetTemperature_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef double (*SetGrainSize_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef double (*SetPorosity_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef double (*SetDensity_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef double (*SetXComponent_f)(MdoodzInstance *instance, Coordinates coordinates);
typedef int (*SetPhase_f)(MdoodzInstance *instance, Coordinates coordinates);

typedef struct {
  SetHorizontalVelocity_f SetHorizontalVelocity;
  SetVerticalVelocity_f   SetVerticalVelocity;
  SetPhase_f              SetPhase;
  SetTemperature_f        SetTemperature;
  SetGrainSize_f          SetGrainSize;
  SetPorosity_f           SetPorosity;
  SetDensity_f            SetDensity;
  SetXComponent_f         SetXComponent;
} SetParticles_ff;

typedef enum {
  NORTHEAST,
  NORTHWEST,
  SOUTHEAST,
  SOUTHWEST,
  INTERNAL,
  NORTH,
  SOUTH,
  WEST,
  EAST,
  FREE_SURFACE
} POSITION;

typedef char (*SetBCVxType_f)(MdoodzInstance *instance, POSITION position);
typedef char (*SetBCVzType_f)(MdoodzInstance *instance, POSITION position);
typedef char (*SetBCPType_f)(MdoodzInstance *instance, POSITION position);
typedef char (*SetBCTType_f)(MdoodzInstance *instance, POSITION position);
typedef char (*SetBCTTypeNew_f)(MdoodzInstance *instance, POSITION position);

typedef double (*SetBCVxValue_f)(MdoodzInstance *instance, POSITION position, Coordinates gridCoordinates);
typedef double (*SetBCVzValue_f)(MdoodzInstance *instance, POSITION position, Coordinates gridCoordinates);
typedef double (*SetBCTValue_f)(MdoodzInstance *instance, POSITION position, double particleTemperature);
typedef double (*SetBCTValueNew_f)(MdoodzInstance *instance, POSITION position, double particleTemperature);

typedef struct {
  SetBCVxType_f    SetBCVxType;
  SetBCVxValue_f   SetBCVxValue;
  SetBCVzType_f    SetBCVzType;
  SetBCVzValue_f   SetBCVzValue;
  SetBCPType_f     SetBCPType;
  SetBCTType_f     SetBCTType;
  SetBCTValue_f    SetBCTValue;

  SetBCTValueNew_f SetBCTValueNew;
  SetBCTTypeNew_f  SetBCTTypeNew;
} SetBCs_ff;

typedef void (*RunMDOODZ_f)(MdoodzInstance *this);

typedef struct {
  double multiplier;
  int   *phases;
  int    nPhases;
} CrazyConductivity;

struct MdoodzInstance {
  char                      *inputFileName;
  params                     model;
  mat_prop                   materials;
  scale                      scaling;

  BuildInitialTopography_ff *BuildInitialTopography;
  SetParticles_ff           *SetParticles;
  SetBCs_ff                 *SetBCs;

  RunMDOODZ_f               *RunMDOODZ;

  CrazyConductivity         *crazyConductivity;
};

void RunMDOODZ(MdoodzInstance *instance);

#endif
