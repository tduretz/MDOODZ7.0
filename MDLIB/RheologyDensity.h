#ifndef MDOODZ_RHEOLOGYDENSITY_H
#define MDOODZ_RHEOLOGYDENSITY_H

typedef struct {
  double *eta;
  double *d1;
} LocalIterationMutables;

typedef struct {
  double Eii;
  double f_ani;
  double Ag;
  double gam;
  double lam;
  double cg;
  double pg;
  double C_pwl;
  double n_pwl;
  double C_gbs;
  double n_gbs;
  double C_lin;
  double n_lin;
  double m_lin;
  double C_exp;
  double ST;
  double n_exp;
  double eta_cst;
  double eta_el;
  int    elastic;
  int    peierls;
  int    dislocation;
  int    diffusion;
  int    constant;
  int    gbs;
  int    phase;
  double d;
} LocalIterationParams;

#endif