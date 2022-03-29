#ifndef MDOODZ_H
#define MDOODZ_H

#define zeroC 273.15
#define Rg    8.314510
#define PI    3.14159265359
#define Rad_Earth 6370000

// BC is a boundary condition structure for the mechanical solver
typedef struct {
  char *type;
  double *val;
} BC;

// BCT is a boundary condition structure for the thermal solver
typedef struct {
  char *type, *typW, *typE, *typS, *typN;
  double *val, *valW, *valE, *valS, *valN;
} BCT;

// grid contains all the fine grid arrays (double *)
typedef struct {
  int Nx, Nz, NN, NC;
  double dx, dz;
  double *roger_x, *roger_z, *div_u, *div_u_s, *div_u_el, *div_u_pl, *div_u_r,
      *u_in, *v_in, *p_in, *p_corr, *sxxd, *szzd, *sxz, *exxd, *ezzd, *exz,
      *VE_s, *VE_n, *sxxd0, *szzd0, *sxz0, *mu_s, *mu_n, *u_adv, *v_adv,
      *eta_phys_n, *kx, *kz, *Cv, *Qr, *eta_phys_s, *u_start, *v_start,
      *p_start, *divth0_n, *T0_n;
  int *iter_smooth;
  int *nb_part_cell, *nb_part_vert;
  BC BCu, BCv, BCp, BCp_exp;
  // TODO rename BCu and rename BCv to BCVx and BCVz
  // TODO rename BCp to BCP (P as pressure is always capital P)
  // TODO find and and do the same for (T as temperature)
  BCT BCt, BCg, BCc;
  // TODO rename BCg to BCv (v is not capital and stands for vertex)
  double *xg_coord, *zg_coord, *xc_coord, *zc_coord, *xvz_coord, *zvx_coord,
      *xg_coord0, *zg_coord0, *xg_coord_ext, *zg_coord_ext;
  double *eta_s, *eta_n, *rho_s, *rho_n;
  // TODO eta_s -> eta_v (for vertex) eta_n -> eta_c (c for centroids)
  double *X_s, *X_n, *X0_s, *X0_n, *p0_n, *p0_s;
  double *OverS_n;
  double *strain_n, *strain_s;
  double *u, *v, *p;
  double *ru, *rv, *rp;
  double *rhs_u, *rhs_v, *rhs_p, *rhs_t, *gx, *gz;
  double p_scale;
  double *alp, *bet_n, *bet_s, *p_lith, *p_lith0, *dp, *Qrho;
  double *VxVz, *VzVx;
  int *P2N, *P2C;
  int *kvx, *lvx, *kvz, *lvz, *kp, *lp, *kn, *ln;
  double **phase_perc_n, **phase_perc_s, **phase_eta_n,
      **phase_eta_s; // TODO: refactor to _c and _v
  double *sxxd0_s, *szzd0_s, *sxz0_n, *exxd_s, *ezzd_s, *exz_n, *sxz_n;
  double *rho0_n;
  double Uthermal, Uelastic, Work, Tii_mean, Eii_mean, T_mean, P_mean;
  double *Work_time, *Uelastic_time, *Uthermal_time, *Time_time, *Short_time,
      *P_mean_time, *T_mean_time, *Tii_mean_time, *Eii_mean_time;
  double *T, *dT, *d_n, *d0_n, *phi_n, *phi0_n;
  double *eII_el, *eII_pl, *eII_pl_s, *eII_pwl, *eII_exp, *eII_lin, *eII_gbs,
      *eII_cst;
  double *eII_pwl_s;
  double *exx_el, *ezz_el, *exz_el, *exx_diss, *ezz_diss, *exz_diss;
  int *comp_cells;
  // For Newton iterations
  double *D11_n, *D12_n, *D13_n, *D14_n;
  double *D21_n, *D22_n, *D23_n, *D24_n;
  double *D31_s, *D32_s, *D33_s, *D34_s;
  double *detadexx_n, *detadezz_n, *detadgxz_n, *detadp_n;
  double *ddivpdexx_n, *ddivpdezz_n, *ddivpdgxz_n, *ddivpdp_n;
  double *detadexx_s, *detadezz_s, *detadgxz_s, *detadp_s;
  double *drhodp_n;
  double *phi0_s, *d0_s, *T_s, *P_s;
  // For anisotropy
  double *FS_AR_n, *FS_AR_s, *aniso_factor_n, *aniso_factor_s;
  double *d1_n, *d2_n, *d1_s, *d2_s;

  double *cell_min_z, *cell_max_z, *vert_min_z, *vert_max_z;
  double *dil_n, *dil_s, *fric_n, *fric_s, *C_n, *C_s;
  double *exz_n_el, *exz_n_diss, *exz_n_pl, *Wdiss, *Wel, *Wtot;
  double *kc_x, *kc_z;
  double *FreeSurfW_s, *FreeSurfW_n;
  double *noise_n, *noise_s;
} grid;

// address Williams comments
typedef enum { ARITHMETIC = 0, HARMONIC = 1, GEOMETRIC = 2 } ETA_AVG;

// params contains the model parameters
typedef struct {
  double xmin, zmin, xmax, zmax, time, dx, dz, dt, dt0, dt_start, dt_max, L0,
      dt_min;
  double xmin0, zmin0, xmax0, zmax0;
  double gx, gz;
  int Nx, Nz, Nt, step, nit, Newton, noisy;
  ETA_AVG eta_avg;
  int itp_stencil;
  double nexp_radial_basis;
  int ismechanical, isperiodic_x, isinertial, iselastic, isnonnewtonian,
      isthermal, ispureshear_ale, free_surf, write_markers, write_debug;
  double free_surf_stab;
  int dt_constant, RK, line_search, thermal_eq, subgrid_diff, adiab_heat,
      shear_heat, advection, fstrain, ConservInterp;
  int surf_processes, cpc, surf_remesh, loc_iter, therm_pert, surf_ised1,
      surf_ised2, MantleID, topografix, Reseed, SmoothSoftening;
  double EpsBG, DivBG, user0, user1, user2, user3, user4, user5, user6, user7,
      user8;
  char *input_file;
  int Nb_phases;
  int ncont;
  double Courant, mineta, maxeta;
  // Particles
  int initial_noise, initial_part;
  // Linear solver
  int decoupled_solve, lsolver, diag_scaling, pc_type;
  double penalty, abs_tol_div, rel_tol_div, auto_penalty, compressible,
      rel_tol_KSP;
  // Non-linear solver
  double line_search_min, safe_dt_div;
  int safe_mode, nstagmax;
  // Deformation maps
  int nT, nE, nd, def_maps;
  double Pn, Tmin, Tmax, Emin, Emax, dmin, dmax, PrBG, TBG;
  // Surface processes
  double surf_diff, surf_sedirate, surf_baselev, surf_Winc, surf_Vinc;
  // Initial thermal perturbation
  double therm_pert_x0, therm_pert_z0, therm_pert_dT, therm_pert_rad,
      cooling_time;
  // For rheological database...
  int force_act_vol_ast;
  double act_vol_dis_ast, act_vol_dif_ast;
  // Phase diagrams
  int isPD, num_PD, *PDMnT, *PDMnP, *PD1DnP;
  double **PDMrho, *PDMTmin, *PDMTmax, *PDMPmin, *PDMPmax;
  double **PD1Drho, *PD1Dmin, *PD1Dmax;
  // Kinetics
  int kin_nP, kin_nT;
  double *kin_dG, kin_Tmin, kin_Tmax, kin_Pmin, kin_Pmax;
  // Visualisation
  int rec_T_P_x_z, delete_breakpoints, GNUplot_residuals;
  // Boundary conditions type
  int BC_setup_type, shear_style, polar;
  int StressRotation, StressUpdate, DirectNeighbour;
  // For diffused rheological constrasts
  int diffuse_X, diffuse_avg;
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

// markers is the particles structure
typedef struct {
  int Nx_part, Nz_part, Nb_part, Nb_part_max, min_part_cell, Nb_part_ini;
  double *x, *z, *Vx, *Vz, *P, *sxxd, *szzd, *sxz, *progress, *T, *d, *phi, *X,
      *syy, *dsyy;
  double *strain, *strain_el, *strain_pl, *strain_pwl, *strain_exp, *strain_lin,
      *strain_gbs;
  int *phase, *generation, *dual;
  int *intag;
  double *Fxx, *Fxz, *Fzx, *Fzz, *nx, *nz;
  double *T0, *P0, *x0, *z0, *Tmax, *Pmax, *divth;
  double *dsxxd, *dszzd, *dsxz;
  double *noise, *rho;
} markers;

// mat_prop contains information related to material phase properties
typedef struct {
  int Nb_phases;
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
  int kin[20], gs[20], cstv[20], pwlv[20], linv[20], expv[20], gbsv[20],
      phase_diagram[20], density_model[20];
  double C_end[20], phi_end[20], psi_end[20], pls_start[20], pls_end[20],
      eta_vp[20], n_vp[20];
  int plast[20], phi_soft[20], psi_soft[20], coh_soft[20], is_tensile[20];
  double Pr[20], tau_kin[20], dPr[20], k_chem[20];
  int reac_soft[20], reac_phase[20];
  int phase_mix[20], phase_two[20];
  double aniso_factor[20], aniso_angle[20];
} mat_prop;

// Structure surface contains free surface data
typedef struct {
  double *a, *b, *height, *vx, *vz, *a0, *b0, *height0;
  int *VertInd;
} surface;

char *GetSetupFileName(int nargs, char *args[]);
void MinMaxArray( double * array, double scale, int size, char* text );

int RunMDOODZ(char *inputFileName,
              void BuildInitialTopography(markers *topo_chain, params model,
                                          scale scaling),
              void SetParticles(markers *particles, scale scaling, params model,
                                mat_prop *materials),
              void SetBCs(grid *mesh, params *model, scale scaling,
                          markers *particles, mat_prop *materials,
                          surface *topo));

#endif
