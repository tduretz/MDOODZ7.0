// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2018  MDOODZ Developper team
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
#include "stdbool.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "mdoodz-private.h"
#include "include/mdoodz.h"
#include "RheologyDensity.h"

#include "mdoodz-log.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

// =========================================================================
// Per-phase analytic inverses of aniso_delta_fn — used by the init-from-
// finite-strain helper (see RheologyParticles.c::aniso_init_finite_strain
// _for_marker). Each helper returns γ_eff such that aniso_delta_fn(FS_AR(
// γ_eff)) == δ, for the saturating-exponential family
//     δ(γ_eff) = c2 · (1 − exp(−γ_eff / γ_e)) + 1
// whose analytic inverse is
//     γ_eff(δ) = − γ_e · ln(1 − (δ − 1)/c2)            (δ ∈ [1, c2+1))
// Out-of-range (δ above the saturation ceiling) is clamped to (c2+1) − ε
// with a one-shot LOG_WARN; the caller proceeds with the saturated γ_eff.
// Below-1 input is clamped to 1.0 (returns 0, F_init = I).
// =========================================================================

// Hansen olivine (case 1): δ = 24.5·0.536·(1 − exp(−γ_eff/3.96)) + 1
//   c2 = 24.5·0.536 = 13.132,  γ_e = 3.96,  δ_max = 14.132
double aniso_delta_inv_hansen( double delta ) {
    const double c2      = 13.132;
    const double gamma_e = 3.96;
    const double eps     = 1.0e-9;
    const double delta_max = 1.0 + c2;        // 14.132
    double d = delta;
    if ( !isfinite(d) || d < 1.0 ) {
        // Below the physical floor — caller should have caught this; clamp.
        d = 1.0;
    }
    if ( d > delta_max - eps ) {
        // One-shot warn (suppressed on repeat) — clamp to just-below-saturation.
        static int warned = 0;
        if ( !warned ) {
            LOG_WARN("aniso_delta_inv_hansen: delta=%g exceeds Hansen olivine saturation %g; clamped to %g", delta, delta_max, delta_max - eps);
            warned = 1;
        }
        d = delta_max - eps;
    }
    // Analytic inverse:  γ_eff = − γ_e · ln(1 − (δ − 1)/c2)
    const double frac = (d - 1.0) / c2;
    return -gamma_e * log( 1.0 - frac );
}

// Compute the finite-strain anisotropy factor δ given the per-cell strain-
// ellipse aspect ratio FS_AR, the user-set ceiling aniso_fac_max, the
// per-phase form selector ani_fstrain, and (when ani_fstrain == 2) the
// per-phase mineral-specific function pointer aniso_delta_fn populated by
// ReadDataAnisotropy() in MDLIB/FlowLaws.c.
//
//   ani_fstrain == 1: δ = min(FS_AR, aniso_fac_max)        — MDOODZ default
//   ani_fstrain == 2: δ = aniso_delta_fn(FS_AR), capped    — mineral-specific,
//                     populated by FlowLaws.c (case 1 = Hansen olivine, case 2
//                     = future calcite, case 3 = future quartz, …)
//                     IF aniso_d_threshold > 0, an additional d-coupling
//                     modifier randomizes CPO toward δ = 1 when the local
//                     grain size drops below the threshold (DisGBS / GBM
//                     regime). Default threshold is 0 (modifier inactive)
//                     for cases 1–8; case 9 (plagioclase BRACKETED) sets
//                     it to 5e-6 m by default per Mehl & Hirth (2008).
//   ani_fstrain == 3: ani_fstrain==2 δ-dispatch + a temperature-dependent
//                     kinetic relaxation of δ toward the isotropic limit
//                     (δ → 1) on the Boneh et al. 2021 (G3) discontinuous-
//                     static-recrystallization timescale. The relaxed δ is
//                     a per-marker integrated state variable (markers->
//                     aniso_delta), updated inside FiniteStrainAspectRatio
//                     and P2G'd to mesh->aniso_delta_n/_s; the ==3 arm here
//                     just returns that grid value (passed in as the
//                     relaxed_delta argument) instead of recomputing
//                     aniso_delta_fn(FS_AR). See FiniteStrainAspectRatio in
//                     RheologyParticles.c and DeltaRelaxationTau below.
//
// All branches honour aniso_fac_max as a hard ceiling. The dispatch is
// per-phase, so heterogeneous mineral assemblages can mix forms. This file
// has no mineral-specific knowledge — every calibration lives in FlowLaws.c.
//
//   d-coupling modifier formula (when aniso_d_threshold > 0 and d > 0):
//     if d < d_threshold:
//         f      = (d / d_threshold)^d_decay     in [0, 1]
//         δ_eff  = 1 + f · (δ − 1)
//     else: δ_eff = δ unchanged.
//   Limit cases: as d → 0, δ → 1 (isotropic limit, complete CPO randomization).
//                as d → d_threshold from below, δ → δ continuously.
//                d_decay = 1.0 → linear blend; > 1.0 → modifier kicks in
//                more sharply only near d=0; < 1.0 → gentler decay.
//
//   relaxed_delta: the P2G'd, temperature-relaxed δ at this grid point
//                  (mesh->aniso_delta_n/_s) — consumed ONLY by the
//                  ani_fstrain == 3 arm. For ani_fstrain ∈ {0,1,2} it is
//                  ignored and those arms stay byte-identical.
double AnisoFactorEvolv( double FS_AR, double aniso_fac_max, int ani_fstrain,
                         double (*aniso_delta_fn)(double FS_AR),
                         double grain_size,
                         double aniso_d_threshold,
                         double aniso_d_decay,
                         double relaxed_delta ) {
  if (ani_fstrain == 3) {
    // ani_fstrain==3: the relaxed δ is computed elsewhere — it is a per-marker
    // integrated state variable (markers->aniso_delta), updated inside
    // FiniteStrainAspectRatio by the Boneh-2021 DiSRX relaxation law and then
    // P2G'd to mesh->aniso_delta_n/_s. UpdateAnisoFactor passes that grid
    // value in as relaxed_delta; here we just clamp it to the ceiling and
    // return it (the relaxation routine already floors it at 1.0). This keeps
    // the relaxed δ flowing through the same per-phase phase_perc weighting
    // and harmonic/geometric averaging as every other ani_fstrain arm.
    return MINV(relaxed_delta, aniso_fac_max);
  }
  if (ani_fstrain == 2) {
    if (aniso_delta_fn == NULL) {
      LOG_ERR("AnisoFactorEvolv: ani_fstrain=%d but aniso_delta_fn is NULL — "
              "did the input file set aniso_db?", ani_fstrain); exit(12);
    }
    double delta = aniso_delta_fn(FS_AR);
    // d-coupling modifier — only active when threshold is set AND d is
    // valid AND d is below threshold. Default behavior (threshold=0) is
    // byte-identical to pre-modifier code.
    if ( aniso_d_threshold > 0.0 && grain_size > 0.0
         && grain_size < aniso_d_threshold ) {
      const double f = pow(grain_size / aniso_d_threshold, aniso_d_decay);
      delta = 1.0 + f * (delta - 1.0);
    }
    return MINV(delta, aniso_fac_max);
  }
  // ani_fstrain == 1: default MDOODZ behaviour (raw FS_AR capped, no Hansen)
  return MINV(FS_AR, aniso_fac_max);
}

// ===========================================================================
//  ani_fstrain == 3 — temperature-dependent δ-relaxation (Boneh et al. 2021)
// ===========================================================================
//
// MODELLING LEAP (caveat, design D7c): Boneh et al. 2021 (G3, "The effect of
// discontinuous static recrystallization on olivine crystallographic preferred
// orientation") establish the *kinetics* of discontinuous static recrystalli-
// zation (DiSRX) — the grain-boundary-migration velocity V and the time for an
// annealed grain to grow — and they show DiSRX *qualitatively* weakens CPO.
// They do NOT give a closed-form δ(t) relaxation law. The construction below —
// "δ relaxes toward the isotropic limit (δ → 1) on the DiSRX grain-growth
// timescale τ_relax = L_relax / V" — is THIS MODEL'S OWN, not Boneh's. It is a
// defensible first-order coupling (when grain boundaries sweep a marker faster,
// its inherited CPO is overwritten faster), but it is a modelling choice; the
// quantitative δ(t) curve is not lab-validated. See
// misc/aniso_fstrain/notes/ for the derivation and the tiered validation plan.
//
// Boneh et al. 2021 grain-boundary-migration kinetics. These constants are
// mineral-specific and now live in the per-phase material database
// (mat_prop::ani_relax_*), set in ReadDataAnisotropy (FlowLaws.c) — mirroring
// how the ani_fstrain==2 δ-calibrations are dispatched per phase. The full
// citation, the M0-units note, and the Δρ-proxy caveats travel WITH the
// constants and now live next to their initialisation in FlowLaws.c. The
// values currently populated there are the OLIVINE calibration; this routine
// stays mineral-agnostic and just consumes whatever the caller passes in.

// Dislocation-density driving force Δρ, proxied from the per-marker accumulated
// dislocation-creep strain strain_pwl (design D3, O2). MDOODZ tracks no
// dislocation density, so Δρ is derived on the fly from a field it already
// carries. Saturating exponential mapping: Δρ → drho_min at zero stored strain,
// → drho_max for strain ≫ eps_ref; result is always finite and in
// [drho_min, drho_max].
//
// drho_min, drho_max, eps_ref are passed in per-phase from
// materials->ani_relax_drho_{min,max} / ani_relax_eps_ref[phase] (SI units).
// The Δρ-proxy CAVEATS (monotonicity / no-recovery; absolute stored strain vs.
// boundary contrast — design D7b) and the calibration rationale travel WITH
// the constants — see ReadDataAnisotropy in FlowLaws.c.
static double DeltaRhoProxy( double strain_pwl,
                             double drho_min, double drho_max, double eps_ref ) {
  double s = strain_pwl;
  if ( s < 0.0 || !isfinite(s) ) s = 0.0;          // guard: finite & non-negative
  const double frac = 1.0 - exp( -s / eps_ref );
  double drho = drho_min + (drho_max - drho_min) * frac;
  if ( drho < drho_min ) drho = drho_min;
  if ( drho > drho_max ) drho = drho_max;
  return drho;
}

// Relaxation timescale τ_relax = L_relax / V (design D1), with the grain-
// boundary-migration velocity V = M0 · exp(−Q/RT) · μ · b² · Δρ.
//
// V uses ONLY the strain-energy driving force F_s = μ·b²·Δρ; the surface-
// energy term F_b = 3γ/d from Boneh's full ΣF = F_s + F_b is dropped (caveat,
// design D7d): ΣF ≈ F_s, justified because F_s ≫ F_b for mantle grain sizes
// ≥ 1 mm (Boneh 2021 §4.2, Fig. 8).
//
// The Boneh-2021 kinetics constants (Q, M0, mu, b) and the Δρ-proxy range
// (drho_min, drho_max, eps_ref) are mineral-specific and are passed in
// per-phase from materials->ani_relax_*[phase] (all in SI units; see
// ReadDataAnisotropy in FlowLaws.c for the calibration and caveats).
//
// T_scaled, L_relax_scaled and R_scaled are in MDOODZ scaled units; the
// routine converts to SI, computes τ_relax in seconds, then converts back to
// scaled time. Returns a (scaled) τ_relax; the caller does the analytic
// exponential update. A non-positive or non-finite L_relax / T / Δρ yields a
// huge τ (→ effectively frozen), never a crash or a NaN.
double DeltaRelaxationTau( double T_scaled, double L_relax_scaled,
                           double strain_pwl, double R_scaled, scale scaling,
                           double Q, double M0, double mu, double b,
                           double drho_min, double drho_max, double eps_ref ) {
  const double T_SI       = T_scaled       * scaling.T;   // K
  const double L_relax_SI = L_relax_scaled * scaling.L;   // m
  // R_scaled is materials->R (scaled); R_SI = R_scaled * (scaling.J/scaling.T).
  const double R_SI       = R_scaled * (scaling.J / scaling.T);   // J/(mol·K)
  if ( !(T_SI > 0.0) || !(L_relax_SI > 0.0) || !isfinite(T_SI) || !isfinite(L_relax_SI) ) {
    return 1.0e300;  // degenerate input → effectively frozen
  }
  const double d_rho = DeltaRhoProxy( strain_pwl, drho_min, drho_max, eps_ref ); // m^-2, always > 0
  // V = M0 · exp(−Q/RT) · μ · b² · Δρ   [m/s]
  const double V_SI = M0 * exp( -Q / (R_SI * T_SI) )
                      * mu * b * b * d_rho;
  if ( !(V_SI > 0.0) || !isfinite(V_SI) ) {
    return 1.0e300;  // cold marker: exp(−Q/RT) underflowed → V → 0 → τ → ∞
  }
  const double tau_SI = L_relax_SI / V_SI;                 // s
  return tau_SI / scaling.t;                               // scaled time
}

// Squared second invariant of deviatoric stress
double Y2( Tensor2D *T, double ani_fac ) {
  // printf("Y2 %2.2e\n", ani_fac);
    return 0.5*( pow(T->xx,2) + pow(T->zz,2) + pow(T->yy,2)) + pow(ani_fac*T->xz,2);
}

// Squared second invariant of deviatoric strain rate
double I2( Tensor2D *E ) {
    return 0.5*( pow(E->xx,2) + pow(E->zz,2) + pow(E->yy,2)) + pow(E->xz,2);
}

// Solution of 2x2 system: TODO: move to MiscFunction.c, make x and f vectors and A matrix
void Solve2x2(double *x1, double *x2, double f1, double f2, double a11, double a12, double a21, double a22) {
  const double det = a11 * a22 - a12 * a21; // determinant
  const double ai  =  1.0 / det * a22;      // inverse matrix components
  const double bi  = -1.0 / det * a12;
  const double ci  = -1.0 / det * a21;
  const double di  =  1.0 / det * a11;
  *x1 = ai * f1 + bi * f2;                  // inverse times rhs
  *x2 = ci * f1 + di * f2;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ViscosityConciseAniso( int phase, double lxlz, double lx2, double angle, double ani_fac, double G, double T, double P, double d0, double phi, double X0, double Exx, double Ezz, double Exz, double Txx0, double Tzz0, double Txz0, mat_prop* materials, params *model, scale *scaling, double *Txx, double *Tzz, double *Txz, double* eta_vep, double* Eii_el, double* Eii_pl, double* Eii_pwl, double* Eii_exp, double* Eii_lin, double* Eii_gbs, double* Eii_cst, double *d, double strain_acc, double dil, double fric, double Coh, double P0, double T0,  double *X1, double *OverS, double *Pcorr, double *rho, double beta, double div, double *div_el, double *div_pl, double *div_r, double *Wtot, double *Wel, double *Wdiss, int post_process, int centroid, int final_update, int index ) {
    // !!!!!!!!!!!!!!!!!!!!!!!!
    // ACTUNG FOR GSE:: d is now d0 and d1 is now d
    // !!!!!!!!!!!!!!!!!!!!!!!!
    // General paramaters
    const double tol    = 1e-11, R = materials->R, dt = model->dt, min_eta = model->min_eta, max_eta = model->max_eta;
    // nitmax bumped 20 → 50 to match RheologyDensity.c:205 — bisection-then-
    // Newton occasionally needs >20 Newton iterations when bisection ends
    // ~5 dex from the visco-elastic root and Newton converges quadratically
    // in the final ~3 iters but linearly in the approach.
    const int    nitmax = 50;
    int noisy = 0;
    double eta = 0.0, eta_el = 1e20, eta_cst = 0.0;
    double TmaxPeierls = (1200.0 + zeroC) / scaling->T;// max. T for Peierls
    int    plastic = 0, constant = 0, dislocation = 0, peierls = 0, diffusion = 0, gbs = 0, elastic = model->elastic, kinetics = 0, is_pl = 0;
    double eta_pwl = 0.0, B_pwl = 0.0, C_pwl = 0.0, Q_pwl = materials->Qpwl[phase], V_pwl = materials->Vpwl[phase], n_pwl = materials->npwl[phase], m_pwl = materials->mpwl[phase], r_pwl = materials->rpwl[phase], A_pwl = materials->Apwl[phase], f_pwl = materials->fpwl[phase], a_pwl = materials->apwl[phase], F_pwl = materials->Fpwl[phase], pre_factor = materials->pref_pwl[phase], t_pwl = materials->tpwl[phase];
    const double E_lin = materials->Qlin[phase], V_lin = materials->Vlin[phase], n_lin = materials->nlin[phase], m_lin = materials->mlin[phase], r_lin = materials->rlin[phase], A_lin = materials->Alin[phase], f_lin = materials->flin[phase], a_lin = materials->alin[phase], F_lin = materials->Flin[phase];
    double eta_ve, Tii;
    Tensor2D E_rot, T_rot, T0_rot;

    double tet = angle - M_PI/2.0;

    // printf("%2.2e\n", ani_fac);

    // Clean up
    *eta_vep   = 0.0;
    // Initialise strain rate invariants to 0
    *Eii_exp = 0.0, *Eii_lin = 0.0, *Eii_pl = 0.0, *Eii_pwl = 0.0, *Eii_el = 0.0, *Eii_gbs = 0, *Eii_cst = 0.0;
    *Txx     = 0.0, *Tzz     = 0.0, *Txz    = 0.0;
    *Wtot    = 0.0, *Wdiss   = 0.0, *Wel    = 0.0, *d       = 0.0;
    *div_el  = 0.0, *div_pl  = 0.0, *div_r  = 0.0;
    *X1      = 0.0, *rho     = 0.0, *OverS  = 0.0, *Pcorr   = 0.0;
    // P corr will be corrected if plasticity feedbacks on pressure (dilation)
    *Pcorr = P;
    // Constant grain size
    *d = d0;
    // Activate deformation mechanisms
    if ( materials->cstv[phase] !=0                  ) constant    = 1;
    if ( materials->pwlv[phase] !=0                  ) dislocation = 1;
    if ( materials->expv[phase] !=0 && T<TmaxPeierls ) peierls     = 1;
    if ( materials->linv[phase] !=0                  ) diffusion   = 1;
    if ( materials->plast[phase]!=0                  ) plastic     = 1;
     // Turn of elasticity for the initialisation step  (viscous flow stress)
    if ( model->step    == 0                         ) elastic     = 0;
    // Visco-plastic incompressible limit
    if ( elastic==0                 ) { G = 1e1; dil = 0.0;}; //K = 1e1;
    // Zero C limit
    if ( T< zeroC/scaling->T        ) T = zeroC/scaling->T;
     // Precomputations
    if ( dislocation == 1 ) {
      B_pwl = pre_factor * F_pwl * pow(A_pwl,-1.0/n_pwl) * exp( (Q_pwl + P*V_pwl)/R/n_pwl/T ) * pow(d0, m_pwl/n_pwl) * pow(f_pwl, -r_pwl/n_pwl) * exp(-a_pwl*phi/n_pwl);
      C_pwl   = pow(2.0*B_pwl, -n_pwl);
    }

    double Ea_lin = materials->Qlin[phase], Va_lin = materials->Vlin[phase];
    double B_lin = 0.0, C_lin = 0.0;

    if ( diffusion == 1 ) {
      if (m_lin>0.0 && *d < 1e-13/scaling->L) {
        LOG_WARN("Cannot run with grain size dependent viscosity if grain size is set to 0 --> d = %2.2e!!!", *d*scaling->L);
        exit(1);
      }
      B_lin = F_lin * pow(A_lin,-1.0/n_lin) * exp( (Ea_lin + P*Va_lin)/R/n_lin/T ) * pow(f_lin, -r_lin/n_lin) * exp(-a_lin*phi/n_lin); // * pow(d, m_lin/n_lin) !!!!!!!!!!!!!!!!!!!!!!!!
      C_lin = pow(2.0*B_lin, -n_lin);
    }

    double Ea_exp = materials->Qexp[phase], S_exp = materials->Sexp[phase], E_exp = materials->Eexp[phase], F_exp = 0.0;
    double gamma = materials->Gexp[phase], ST=0.0, n_exp = materials->nexp[phase];
    double  B_exp = 0.0, C_exp = 0.0;

    if ( peierls   == 1 ) {
      ST                           = Ea_exp/R/T * 2.0*gamma*(1-gamma);
      // new
      double Arr_exp = exp(-Ea_exp/R/T*pow(1.0-gamma,2.0));
      F_exp = pow( pow(2.0,1.0-ST-n_exp) / pow(sqrt(3.0), ST+n_exp+1.0), 1.0/(ST+n_exp));
      B_exp = F_exp * ( pow(gamma*S_exp, ST/(ST+n_exp)) / pow( E_exp*Arr_exp, 1.0/(ST+n_exp)) );
      C_exp = pow(2.0*B_exp, -(ST+n_exp));
    }

    // double cos_fric = cos(fric);
    // double sin_fric = sin(fric);
    double Tyield = 0.0;
    double  K = 1.0 / beta;
    // double soon = 0.0;


    if ( plastic==1 ) {
      Tyield     =  Coh*cos(fric) + P*sin(fric); // -  soon*sin(fric)/3.0*( a1*T_rot.xx + a2*T_rot.zz + a3*T_rot.yy);

      // Von-Mises cut-off
      if (materials->Slim[phase] < Tyield) {
        Coh        = materials->Slim[phase];
        Tyield     = materials->Slim[phase];
      }
    }

    // Transform strain rate: Q*E*Qt ---> 13/11/24 checked against rotation formula
    const double nz2 = 1.0-lx2;
    E_rot.xx  =   lx2*Exx +  nz2*Ezz +   2.*lxlz*Exz;
    E_rot.zz  =   nz2*Exx +  lx2*Ezz -   2.*lxlz*Exz;
    E_rot.xz  = -lxlz*Exx + lxlz*Ezz + (lx2-nz2)*Exz;
    E_rot.yy  = -E_rot.xx - E_rot.zz;                 // !!! out-of-plane component !!!
    // Local iterations to determine eta_vep and ani_vep
    double f0_n, f0_s, f_n, f_s, dfnden, dfndes, dfsden, dfsdes;
    double I2_v, Y2_v, Eii_vis, W_n, W_s, ieta_pwl, deta_n, deta_s;
    // Construct initial guess for viscosity: Isolated viscosities
    double Eii = sqrt(I2( &E_rot ));
    if (constant)    eta_cst  = materials->eta0[phase];
    if (dislocation) eta_pwl  = B_pwl * pow( Eii, (1-n_pwl)/n_pwl );
    if (elastic)     eta_el   = G*dt;

    double eta_lin = 0.0, eta_exp = 0.0;
    if ( diffusion   == 1 )              eta_lin  = B_lin * pow( Eii, 1.0/n_lin - 1.0 ) * pow(*d, m_lin/n_lin);
    if ( peierls     == 1 )              eta_exp  = B_exp * pow( Eii, 1.0 / (ST + n_exp) - 1.0);

    // Cases where no local iterations are needed
    // if (dislocation) eta_ve = eta_pwl;
    // if (constant)    eta_ve = eta_cst;
    // T_rot.xx = 2.0 * eta_ve * E_rot.xx;
    // T_rot.zz = 2.0 * eta_ve * E_rot.zz;
    // T_rot.xz = 2.0 * eta_ve * E_rot.xz/ani_fac;
    // T_rot.yy = -T_rot.xx - T_rot.xx;

    // Define viscosity bounds
    double eta_up = 1.0e100 / scaling->eta;
    double eta_lo = -1.0e100 / scaling->eta;
    if (constant)    { eta_up = MINV(eta_up, eta_cst); eta_lo = MAXV(eta_lo, eta_cst); }
    if (dislocation) { eta_up = MINV(eta_up, eta_pwl); eta_lo = MAXV(eta_lo, eta_pwl); }
    if (elastic)     { eta_up = MINV(eta_up, eta_el);  eta_lo = MAXV(eta_lo, eta_el);  }
    if (peierls)     { eta_up = MINV(eta_up, eta_exp); eta_lo = MAXV(eta_lo, eta_exp); }
    if (diffusion)   { eta_up = MINV(eta_up, eta_lin); eta_lo = MAXV(eta_lo, eta_lin); }

    // if (centroid==0 && index==65) noisy = 1;

    // Initial guess
    eta_ve = eta_up;
    // Iterations for visco-elastic trial
    if (noisy) LOG_INFO("Start local VE iterations for phase %d (cst: %d --- el: %d --- pwl: %d)", phase, constant, elastic, dislocation);

    // GSE active under anisotropy: route through the wattmeter-aware local
    // iteration shared with the non-aniso path. The Eii passed in is the
    // rotated-frame I2 invariant which equals the lab-frame I2 by rotation
    // invariance, so the iteration mathematics is identical to the non-aniso
    // call site at RheologyDensity.c:644. See merge-gse-anisotropy/design.md
    // D2-D3 for the rationale.
    const int gs = (materials->gs[phase] != 0) ? 1 : 0;
    double r = 0.0;
    if (gs) {
      const double Kg  = materials->Kpzm[phase], Qg = materials->Qpzm[phase];
      const double Ag  = Kg * exp(-Qg / (R * T));
      const double gam = materials->Gpzm[phase];
      const double lam = materials->Lpzm[phase];
      const double cg  = materials->cpzm[phase];
      const double pg  = materials->ppzm[phase];

      LocalIterationParams params = {
              .noisy       = 0,
              .phase       = phase,
              .Eii         = Eii,
              .gam         = gam,
              .lam         = lam,
              .cg          = cg,
              .pg          = pg,
              .Ag          = Ag,
              .C_pwl       = C_pwl,
              .n_pwl       = n_pwl,
              .C_lin       = C_lin,
              .n_lin       = n_lin,
              .m_lin       = m_lin,
              .ST          = ST,
              .n_exp       = n_exp,
              .elastic     = elastic,
              .eta_el      = eta_el,
              .peierls     = peierls,
              .C_exp       = C_exp,
              .dislocation = dislocation,
              .diffusion   = diffusion,
              .constant    = constant,
              .eta_cst     = eta_cst,
              .gbs         = 0,             // aniso path doesn't track GBS yet
              .C_gbs       = 0.0,
              .n_gbs       = 1.0,
              .d           = *d,
              .eta_up      = eta_up,
              .eta_lo      = eta_lo,
              .T           = T };
      LocalIterationViscoElasticGrainSize(
        (LocalIterationMutables){ .eta = &eta_ve, .Tii = &Tii,
                                   .Eii_cst = Eii_cst, .Eii_lin = Eii_lin,
                                   .Eii_gbs = Eii_gbs, .Eii_pwl = Eii_pwl,
                                   .Eii_exp = Eii_exp, .d1 = d },
        params);
    } else {
      // Existing inline iteration — Newton on the visco-elastic residual.
      // Runs whenever `gs == 0` (i.e. the phase has no GSE law), which
      // covers all current aniso scenarios in `SETS/`.
      //
      // Bisection-then-Newton (added to fix AnisoNewton overshoot at cold
      // lithospheric nodes where the local-iteration root is multiple decades
      // away from the initial guess). Mirrors RheologyDensity.c:259-310. The initial guess
      // eta_ve = eta_up (= MIN of mechanism viscosities = elastic for cold
      // HK03 wet olivine) lands far outside Newton's basin of convergence
      // when eta_pwl ≫ eta_el (≈11 orders of magnitude at T < 800 K), so
      // r/drdeta overshoots into eta_ve < 0 → pow(Tii,n_pwl) NaN.
      // 20 log-space bisection steps narrow the (eta_up, eta_lo) bracket
      // from up to ~30 dex down to ~3e-5 dex — well inside Newton's basin.
      // Bracket semantics in this codebase: eta_up = MIN of mechanism
      // viscosities (small bound); eta_lo = MAX (large bound). Residual
      // r(eta_ve) = Eii - elastic*Tii/(2 eta_el) - Eii_vis is monotonically
      // decreasing in eta_ve, so r>0 at small eta and r<0 at large eta.
      double r0 = 0.0, drdeta;
      const int bisection_steps = 30;   // 30 → 9e-9 dex resolution on a 30-dex bracket
      // Bracket the visco-elastic root in [eta_bi_lo, eta_bi_hi]:
      //   • lower end = min_eta (user-supplied solver floor — at any
      //     eta_ve below this, Tii is small, viscous and elastic
      //     contributions are negligible, so r ≈ Eii > 0).
      //   • upper end = MAX(eta_lo, max_eta)·10 — eta_lo here is the
      //     MAX-of-mechanism-viscosities (codebase semantics; see
      //     RheologyDensity.c:292-293). At any eta_ve above max viscosity,
      //     Tii blows up and Eii_pwl ∝ Tii^n_pwl makes r ≪ 0. Multiplied
      //     by 10 for safety against degenerate single-mechanism cells
      //     where eta_lo equals the root exactly.
      // For HK03 wet olivine cold (M7a): bracket spans ~13 dex; root sits
      // near eta_el. For warm asthenosphere (M0 phase 0): bracket spans
      // ~8 dex; root sits at eta_pwl. Both contain the root.
      double eta_bi_lo = min_eta * 1.0e-3;            // 3 dex below floor for safety
      double eta_bi_hi = 10.0 * MAXV(eta_lo, max_eta);
      // (the existing post-Newton clamp at min_eta below re-applies the floor.)
      for (int bit = 0; bit < bisection_steps; bit++) {
        // Geometric mean — log-space bisection, halves the dex each step
        eta_ve = sqrt(eta_bi_lo * eta_bi_hi);
        // Same residual evaluation as the Newton loop below; reuse the
        // local Eii_* slots so post-loop state is consistent.
        Tii           = 2.0*eta_ve*Eii;
        if (constant)    *Eii_cst = Tii  / 2.0 / eta_cst;
        if (dislocation) *Eii_pwl = C_pwl * pow(Tii, n_pwl);
        if (peierls)     *Eii_exp = C_exp * pow(Tii, ST + n_exp);
        if (diffusion)   *Eii_lin = C_lin * pow(Tii, n_lin) * pow(*d, -m_lin);
        Eii_vis       =  *Eii_pwl + *Eii_cst + *Eii_exp + *Eii_lin;
        const double r_bi = Eii - elastic*Tii/2/eta_el - Eii_vis;
        // Monotonicity of r in eta_ve: r decreases as eta_ve grows. So
        // r > 0 ⇒ root is at larger eta ⇒ tighten the small end up to eta_ve.
        if (r_bi > 0.0) eta_bi_lo = eta_ve;
        else            eta_bi_hi = eta_ve;
        // Early exit if bisection already meets the tight tolerance.
        if (fabs(r_bi/Eii) < tol/100) break;
      }
      // Start Newton from the bisection-narrowed midpoint (geometric mean
      // of the final bracket) — bisection drove r near zero, Newton now
      // takes the quadratic-convergence tail.
      eta_ve = sqrt(eta_bi_lo * eta_bi_hi);
      for (int it = 0; it < nitmax; it++) {
        Tii           = 2.0*eta_ve*Eii;
        if (constant)    *Eii_cst = Tii  / 2.0 / eta_cst;
        if (dislocation) *Eii_pwl = C_pwl * pow(Tii, n_pwl);
        if (peierls)     *Eii_exp = C_exp * pow(Tii, ST + n_exp);                     // Peierls - power law
        if (diffusion)   *Eii_lin = C_lin * pow(Tii, n_lin) * pow(*d, -m_lin); // !!! gs - dependence !!!
        Eii_vis       =  *Eii_pwl + *Eii_cst + *Eii_exp + *Eii_lin;

        r             = Eii - elastic*Tii/2/eta_el - Eii_vis;
        if (isnan(r))LOG_INFO("r= %2.2e Tii=%2.2e eta_ve=%2.2e Eii=%2.2e %2.2e %2.2e", r, Tii, eta_ve, Eii, elastic*Tii/2/eta_el, Eii_vis);
        if (isnan(r)) exit(1);
        if (it==0) r0 = r;
        drdeta = 0.0;

        if (elastic)     drdeta += - elastic*Eii/eta_el;
        if (dislocation) drdeta += - *Eii_pwl*n_pwl/eta_ve;
        if (constant)    drdeta += - constant*Eii/eta_cst;
        if (peierls)     drdeta += -(*Eii_exp) * (ST + n_exp) / eta_ve;
        if (diffusion)   drdeta += -(*Eii_lin) *n_lin / eta_ve;
        if ((centroid==0 && index==0 && dislocation && !constant) || it > 10) LOG_INFO("AnisoNewton centroid=%d index=%d phase=%d it=%02d Tii=%2.4e eta_ve=%2.4e Eii_pwl=%2.4e r=%2.4e drdeta=%2.4e Eii=%2.4e", centroid, index, phase, it, Tii, eta_ve, *Eii_pwl, r, drdeta, Eii);
        const double eta_ve_pre = eta_ve;
        eta_ve       -= r/drdeta;
        // Positivity guard on eta_ve after Newton update — prevents pow(Tii, n_pwl) NaN at next iter when n_pwl is non-integer; lets Newton settle anywhere in (0, ∞) without interfering with the natural [eta_up, eta_lo] bracket
        if (eta_ve <= 0.0) {
          if ((centroid==0 && index==0) || it > 10) LOG_INFO("AnisoNewton centroid=%d index=%d phase=%d it=%02d RESET eta_ve_pre=%2.4e step=%2.4e new=%2.4e<=0 → min_eta=%2.4e", centroid, index, phase, it, eta_ve_pre, r/drdeta, eta_ve, min_eta);
          eta_ve = min_eta;
        }
        // Make some noise!!!!!
        if (noisy) LOG_INFO("VE It. %02d: r = %2.2e", it, fabs(r/Eii));
        // Exit criteria
        if ( fabs(r/Eii) < tol / 100 ) {
          if (it > 10) LOG_INFO("V-E L.I. Warnung: more that 10 local iterations, there might be a problem... centroid=%d index=%d phase=%d it=%02d", centroid, index, phase, it);
          break;
        }
        // Floating-point cancellation floor — when Eii_vis ≈ Eii at the root,
        // the residual r = Eii - elastic_term - Eii_vis loses ~12 digits of
        // precision via catastrophic cancellation, so |r/Eii| plateaus around
        // 1e-11 even after the iterate has converged to the true root within
        // double precision. Accept the result if we exhaust iterations but
        // are within 10x of the nominal tol — the iterate is correct to the
        // double-precision noise floor, and the strict failure exit below
        // would be a false negative.
        if (it == nitmax - 1 && fabs(r/Eii) < 10.0 * tol) {
          if (centroid == 0 && index == 0)
            LOG_INFO("AnisoNewton centroid=%d index=%d phase=%d hit nitmax at noise floor (r/Eii=%2.2e); accepting", centroid, index, phase, fabs(r/Eii));
          break;
        }
        if (it == nitmax - 1 && (fabs(r/Eii) > tol) ) {
          LOG_INFO("Visco-Elastic iterations failed! centroid=%d index=%d phase=%d Tii=%2.4e eta_ve=%2.4e Eii_pwl=%2.4e r=%2.4e drdeta=%2.4e Eii=%2.4e eta_lo=%2.4e eta_up=%2.4e min_eta=%2.4e", centroid, index, phase, Tii, eta_ve, *Eii_pwl, r, drdeta, Eii, eta_lo, eta_up, min_eta);
          exit(0);
        }
      }
      if (isnan(r))exit(0);
    }

    // Post-Newton physical floor on eta_ve — prevents viscosity runaway (from shear-band localization producing tiny eta_pwl via power-law at high Eii) from breaking Stokes solver SPD condition
    if (eta_ve < min_eta) {
      if (centroid==0 && index==0) LOG_INFO("AnisoNewton centroid=%d index=%d phase=%d POST-FLOOR eta_ve=%2.4e<min_eta=%2.4e → clamped", centroid, index, phase, eta_ve, min_eta);
      eta_ve = min_eta;
    }
    // Post-Newton physical ceiling on eta_ve — mirrors RheologyDensity.c
    // lines 1058-1073 (non-aniso path clamps both ends). Cold or strain-
    // localized anisotropic cells can produce eta_ve ≫ max_eta (observed up
    // to 1e80 with HK03 wet olivine at high ani_fac), inside the bisection
    // bracket [min_eta·1e-3, 10·max_eta] but making the Stokes operator
    // numerically singular (condition number ~ max_eta/min_eta blowing up by
    // ~55 dex). The Powell-Hestenes solver then plateaus at rel.max.div
    // ~ 1e-1 instead of converging to lin_rel_div (~1e-5).
    if (eta_ve > max_eta) {
      eta_ve = max_eta;
    }

    // if (centroid==0 && index==65) printf("ani 1 - eta_ve = %2.4e Eii_pwl = %2.4e Eii=%2.4e C_pwl=%2.4e tol=%2.2e\n", eta_ve, *Eii_pwl, Eii, C_pwl, model->eta_tol);


    // Evaluate stress components
    T_rot.xx = 2.0 * eta_ve * E_rot.xx;
    T_rot.zz = 2.0 * eta_ve * E_rot.zz;
    T_rot.xz = 2.0 * eta_ve * E_rot.xz/ani_fac;
    T_rot.yy = -T_rot.xx - T_rot.xx;
    *eta_vep = eta_ve;

    // Iterations for visco-elasto-viscoplastic correction
    double eta_vp0 = materials->eta_vp[phase], n_vp = materials->n_vp[phase], eta_vp = materials->eta_vp[phase];
    double Ft   = Tii - Tyield;
    double Fc = Ft, Tiic = Tii;
    double gdot = 0.0;
    double eta_pl, divp;
    // double Y1_p_c, Y2_p_c, dFdgdot;
    //double txx = T_rot.xx, tzz = T_rot.zz, tyy = T_rot.yy, txz=T_rot.xz;
    double Pc = P;

    if (Ft>1e-17 && plastic == 1) {
      is_pl = 1;
      if (noisy) LOG_INFO("Start local VP iterations (cst: %d --- el: %d --- pwl: %d)", constant, elastic, dislocation);
      if (noisy) LOG_INFO("Start local VP iterations (cst: %d --- el: %d --- pwl: %d)", constant, elastic, dislocation);
      gdot = Ft / (eta_ve + eta_vp0 + K*dt*sin(dil)*sin(fric)); // soon*2.0*eta_ve*pow(sin(fric),2)/9.0 * (a1*a1 - a1*a3 + a2*a2 - a2*a3));
      Tiic = Tii - gdot*eta_ve;
      Pc   = P   + K*dt*sin(dil)*gdot;
      // txx  = 2*eta_ve*(Exx-gdot*(txx/2/Tii + a1*sin(fric)/3.0) );
      // tzz  = 2*eta_ve*(Ezz-gdot*(tzz/2/Tii + a2*sin(fric)/3.0) );
      // tyy  = -txx-tzz;
      Fc   = Tiic - Coh*cos(fric) - Pc*sin(fric) - eta_vp0*gdot; // + soon*sin(fric)/3.0*( a1*txx + a2*tzz + a3*tyy )
      if (noisy) LOG_INFO("Ft = %2.2e --- Fc = %2.2e", Ft, Fc);

      // Evaluate stress components
      *eta_vep = Tiic/2.0/Eii;
      T_rot.xx = 2.0 * (*eta_vep) * E_rot.xx;
      T_rot.zz = 2.0 * (*eta_vep) * E_rot.zz;
      T_rot.xz = 2.0 * (*eta_vep) * E_rot.xz / ani_fac;
      T_rot.yy = -T_rot.xx - T_rot.zz;
      *div_pl  = gdot*sin(dil);
      *Eii_pl  = gdot/2.0;
      Tii      = Tiic;
      //printf("Ft = %2.2e --- Fc = %2.2e %2.2e\n", Tiic, sqrt(Y2( &T_rot, ani_fac )), T_rot.yy );
    }
    // Back-transform stress: Qt*E*Q ---> 13/11/24 checked against rotation formula
    *Txx =   lx2*T_rot.xx +  nz2*T_rot.zz -   2.*lxlz*T_rot.xz;
    *Tzz =   nz2*T_rot.xx +  lx2*T_rot.zz +   2.*lxlz*T_rot.xz;
    *Txz =  lxlz*T_rot.xx - lxlz*T_rot.zz + (lx2-nz2)*T_rot.xz;

    // Update effective viscosity and anisotropy factor
    *Pcorr   = Pc;
    // Defensive positivity floor on *eta_vep — after the visco-plastic
    // correction at line ~491, Tiic = Tii − gdot·eta_ve can occasionally
    // overshoot to a negative value when the Drucker-Prager plastic
    // strain-rate update is asked to absorb a Tii that exceeds Tyield by
    // more than the elastic stiffness allows. The negative *eta_vep
    // would propagate into mesh->eta_s via the arithmetic average at
    // RheologyDensity.c-style averaging, making min. Maxwell = eta_s/G
    // negative and corrupting DefineInitialTimestep. Floor at min_eta
    // matches the post-Newton clamp applied earlier on eta_ve.
    if (*eta_vep < min_eta) {
      *eta_vep = min_eta;
    }
    // Defensive ceiling on *eta_vep — mirror the eta_ve ceiling above and
    // the non-aniso *etaVE clamp in RheologyDensity.c:1058. With the
    // post-Newton eta_ve already ≤ max_eta, the plastic branch can still
    // produce *eta_vep = Tiic/(2·Eii) > max_eta when Tiic ≈ Tii (low-gdot
    // cells just past the yield surface), so the ceiling needs to be
    // applied a second time here.
    if (*eta_vep > max_eta) {
      *eta_vep = max_eta;
    }
    eta      = *eta_vep;
    eta_pl   = Tii/gdot;
    divp     = gdot*sin(dil);

      //-------- Post-Processing
  if ( post_process == 1) {

    // Strain rates: VEP partitioning
    if ( dislocation == 1 )              eta_pwl  = B_pwl * pow( *Eii_pwl, 1.0/n_pwl - 1.0 );
    if ( diffusion   == 1 )              eta_lin  = B_lin * pow( *Eii_lin, 1.0/n_lin - 1.0 ) * pow(*d, m_lin/n_lin);
    if ( peierls     == 1 )              eta_exp  = B_exp * pow( *Eii_exp, 1.0 / (ST + n_exp) - 1.0);

    // Compute dissipative strain rate components
    if (elastic) *Eii_el   = Tii*Tii/(4*eta_el*eta_el);
    if (is_pl)   *Eii_pl   = gdot/2.0;

    // Partitioning of volumetric strain
    *div_pl  = divp;
    *div_el  = - (Pc - P0) / (K*dt);
    // *div_r   = divr;

    double inv_eta_diss = 0.0;
    if (diffusion  == 1)  inv_eta_diss += (1.0/eta_lin);
    if (peierls    == 1)  inv_eta_diss += (1.0/eta_exp);
    if (dislocation== 1)  inv_eta_diss += (1.0/eta_pwl);
    if (constant   == 1)  inv_eta_diss += (1.0/eta_cst);
    if (is_pl      == 1)  inv_eta_diss += (1.0/eta_pl );
    eta        = 1.0/(inv_eta_diss);

    // Viscoplastic overstress
    *OverS = eta_vp*gdot;

    // Work(s)
    if (centroid==1) {
      *Wtot  = Tii*Tii/eta_ve;
      *Wel   = Tii*Tii/eta_el;
      *Wdiss = Tii*Tii/eta;
      // printf("Wtot = %2.6e Wdiss = %2.6e Wel = %2.6e\n", *Wtot, *Wdiss, *Wel);
      // printf("Wtot - (Wdiss + Wel) = %2.6e\n", *Wtot - (*Wdiss + *Wel));
      // printf("eta = %2.6e --- eta_ve = %2.6e --- eta_ve1 = %2.6e\n", eta, eta_ve, 1.0/(1./eta+1./eta_el));
    }
  }
  // Final viscosity-limiter on the returned eta_phys — mirrors
  // RheologyDensity.c lines 1067-1073 (non-aniso path). Caller averages
  // this into mesh->eta_phys_s, which feeds the Stokes assembly; without
  // the cap the operator can develop O(1e60+) condition-number entries on
  // cold-mantle high-aniso cells.
  if (eta > max_eta) eta = max_eta;
  if (eta < min_eta) eta = min_eta;
  return eta;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateAnisoFactor( grid *mesh, mat_prop *materials, params *model, scale *scaling) {

  LOG_INFO("Update anisotropy factor");
    int p, k, l, Nx, Nz, Ncx, Ncz, c0, c1;
  int average = model->ani_average; // SHOULD NOT BE ALLOWED TO BE ELSE THAN 1 - but why??

  Nx = mesh->Nx;
  Nz = mesh->Nz;
  Ncx = Nx-1;
  Ncz = Nz-1;

  // Calculate cell centers shear modulus
  for ( l=0; l<Ncz; l++ ) {
    for ( k=0; k<Ncx; k++ ) {

      // Cell center index
      c0 = k  + l*(Ncx);

      // First - initialize to 0
      mesh->aniso_factor_n[c0] = 0.0;

      // Compute only if below free surface
      if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {

        // Loop on phases
        for ( p=0; p<model->Nb_phases; p++) {

          // Arithmetic
          if (average == 0) {
            if (materials->ani_fstrain[p]==0) mesh->aniso_factor_n[c0] += mesh->phase_perc_n[p][c0] * materials->aniso_factor[p];
            if (materials->ani_fstrain[p]!=0) mesh->aniso_factor_n[c0] += mesh->phase_perc_n[p][c0] * AnisoFactorEvolv( mesh->FS_AR_n[c0], materials->ani_fac_max[p], materials->ani_fstrain[p], materials->aniso_delta_fn[p], mesh->d_n[c0], materials->aniso_d_threshold[p], materials->aniso_d_decay[p], mesh->aniso_delta_n[c0] );
          }
          // Harmonic
          if (average == 1) {
            if (materials->ani_fstrain[p]==0) mesh->aniso_factor_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/materials->aniso_factor[p];
            if (materials->ani_fstrain[p]!=0) mesh->aniso_factor_n[c0] += mesh->phase_perc_n[p][c0] * 1.0/AnisoFactorEvolv( mesh->FS_AR_n[c0], materials->ani_fac_max[p], materials->ani_fstrain[p], materials->aniso_delta_fn[p], mesh->d_n[c0], materials->aniso_d_threshold[p], materials->aniso_d_decay[p], mesh->aniso_delta_n[c0] );
          }
          // Geometric
          if (average == 2) {
            if (materials->ani_fstrain[p]==0) mesh->aniso_factor_n[c0] += mesh->phase_perc_n[p][c0] * log(materials->aniso_factor[p]);
            if (materials->ani_fstrain[p]!=0) mesh->aniso_factor_n[c0] += mesh->phase_perc_n[p][c0] * log(AnisoFactorEvolv( mesh->FS_AR_n[c0], materials->ani_fac_max[p], materials->ani_fstrain[p], materials->aniso_delta_fn[p], mesh->d_n[c0], materials->aniso_d_threshold[p], materials->aniso_d_decay[p], mesh->aniso_delta_n[c0] ));
          }
        }


        if ( isinf(1.0/mesh->aniso_factor_n[c0]) ) {
          LOG_INFO("Cell %d has no neighbouring particles: average=%d, value = %2.2e, 1/value=%2.2e ", c0, average, mesh->aniso_factor_n[c0], 1.0/mesh->aniso_factor_n[c0]);
          for ( p=0; p<model->Nb_phases; p++) {
            LOG_INFO("Phase %d, proportion %2.2e ----> UpdateAnisoFactor()", p, mesh->phase_perc_n[p][c0]);
          }
          exit(1);
        }

        // Post-process for geometric/harmonic averages
        if ( average==1 ) mesh->aniso_factor_n[c0] = 1.0/mesh->aniso_factor_n[c0];
        if ( average==2 ) mesh->aniso_factor_n[c0] = exp(mesh->aniso_factor_n[c0]);
      }
    }
  }


  // Calculate vertices shear modulus
  for ( l=0; l<Nz; l++ ) {
    for ( k=0; k<Nx; k++ ) {

      // Vertex index
      c1 = k + l*Nx;

      // First - initialize to 0
      mesh->aniso_factor_s[c1] = 0.0;

      // Compute only if below free surface
      if ( mesh->BCg.type[c1] != 30 ) {

        // Average grain size at vertex (k, l) from 4 surrounding cell centers
        // (cell flat index = i + j*Ncx). Used by the d-coupling modifier in
        // AnisoFactorEvolv when aniso_d_threshold > 0; otherwise inert.
        double d_at_vertex = 0.0;
        int    d_count     = 0;
        for ( int dj = -1; dj <= 0; dj++ ) {
          int jj = l + dj;
          if ( jj < 0 || jj >= Ncz ) continue;
          for ( int di = -1; di <= 0; di++ ) {
            int ii = k + di;
            if ( ii < 0 || ii >= Ncx ) continue;
            int idx_n = ii + jj * Ncx;
            if ( mesh->d_n[idx_n] > 0.0 ) {
              d_at_vertex += mesh->d_n[idx_n];
              d_count++;
            }
          }
        }
        if ( d_count > 0 ) d_at_vertex /= (double)d_count;

        // Loop on phases
        for ( p=0; p<model->Nb_phases; p++) {

          // Arithmetic
          if (average == 0) {
            if (materials->ani_fstrain[p]==0) mesh->aniso_factor_s[c1] += mesh->phase_perc_s[p][c1] * materials->aniso_factor[p];
            if (materials->ani_fstrain[p]!=0) mesh->aniso_factor_s[c1] += mesh->phase_perc_s[p][c1] * AnisoFactorEvolv( mesh->FS_AR_s[c1], materials->ani_fac_max[p], materials->ani_fstrain[p], materials->aniso_delta_fn[p], d_at_vertex, materials->aniso_d_threshold[p], materials->aniso_d_decay[p], mesh->aniso_delta_s[c1] );
          }
          // Harmonic
          if (average == 1) {
            if (materials->ani_fstrain[p]==0) mesh->aniso_factor_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/materials->aniso_factor[p];
            if (materials->ani_fstrain[p]!=0) mesh->aniso_factor_s[c1] += mesh->phase_perc_s[p][c1] *  1.0/AnisoFactorEvolv( mesh->FS_AR_s[c1], materials->ani_fac_max[p], materials->ani_fstrain[p], materials->aniso_delta_fn[p], d_at_vertex, materials->aniso_d_threshold[p], materials->aniso_d_decay[p], mesh->aniso_delta_s[c1] );
          }
          // Geometric
          if (average == 2) {
            if (materials->ani_fstrain[p]==0) mesh->aniso_factor_s[c1] += mesh->phase_perc_s[p][c1] *  log(materials->aniso_factor[p]);
            if (materials->ani_fstrain[p]!=0) mesh->aniso_factor_s[c1] += mesh->phase_perc_s[p][c1] *  log(AnisoFactorEvolv( mesh->FS_AR_s[c1], materials->ani_fac_max[p], materials->ani_fstrain[p], materials->aniso_delta_fn[p], d_at_vertex, materials->aniso_d_threshold[p], materials->aniso_d_decay[p], mesh->aniso_delta_s[c1] ));
          }

        }

        if ( isinf(1.0/mesh->aniso_factor_s[c1]) ) {
          LOG_INFO("Node %d has no neighbouring particles", c1);
          for ( p=0; p<model->Nb_phases; p++) {
            LOG_INFO("Phase %d, proportion %2.2e ----> UpdateAnisoFactor()", p, mesh->phase_perc_s[p][c1]);
          }
          exit(1);
        }

        // Post-process for geometric/harmonic averages
        if ( average==1 ) mesh->aniso_factor_s[c1] = 1.0/mesh->aniso_factor_s[c1];
        if ( average==2 ) mesh->aniso_factor_s[c1] = exp(mesh->aniso_factor_s[c1]);
      }
    }
  }

  // // Smooth
  // InterpVerticesToCentroidsDouble( mesh->aniso_factor_n,  mesh->aniso_factor_s,  mesh, model );
  // InterpCentroidsToVerticesDouble( mesh->aniso_factor_n,  mesh->aniso_factor_s,  mesh, model );

  // Periodic
  double av;
  if (model->periodic_x==1) {
    // printf("average=%d  %d", average, model->eta_average);
    // exit(1);
    for( l=0; l<Nz; l++) {
      c1 = l*Nx + Nx-1;
      av = 0.5*(mesh->aniso_factor_s[c1] + mesh->aniso_factor_s[l*Nx]);
      mesh->aniso_factor_s[c1] = av; mesh->aniso_factor_s[l*Nx] = av;
    }
  }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NonNewtonianViscosityGridAniso( grid *mesh, mat_prop *materials, params *model, Nparams Nmodel, scale *scaling, int final_update ) {

  int    p, k, l, Nx, Nz, Ncx, Ncz, c0, c1, k1;
  double eta, txx1, tzz1, txz1, eta_vep, ani_vep, ani_e, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, div_el, div_pl, div_r;
  double el = 0.0, Wtot, Wel, Wdiss;
  int    average = model->eta_average;
  double Xreac;
  double OverS;
  double Pcorr, rho;
  double Exx, Ezz, Exz, eta_e, ani, d1, d2, angle, lxlz, lx2, lay_ang;
  double Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det;
  if (model->elastic == 1) el = 1.0;

  double etaVE, VEcoeff, d0; // to delete

  Nx  = mesh->Nx;  Ncx = Nx - 1;
  Nz  = mesh->Nz;  Ncz = Nz-1;

  // Stuff to be interpolated to vertices
  InterpCentroidsToVerticesDouble( mesh->T,       mesh->T_s,     mesh, model );
  InterpCentroidsToVerticesDouble( mesh->p_in,    mesh->P_s,     mesh, model );
  InterpCentroidsToVerticesDouble( mesh->d0_n,    mesh->d0_s,    mesh, model );
  InterpCentroidsToVerticesDouble( mesh->phi0_n,  mesh->phi0_s,  mesh, model ); // ACHTUNG NOT FRICTION ANGLE

  // Evaluate cell center viscosities
#pragma omp parallel for shared( mesh ) private( k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, ani_vep, ani_e, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, Wtot, Wel, Wdiss, Xreac,OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, eta_e, ani, d1, d2, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2, lay_ang ) firstprivate( el, materials, scaling, average, model, Ncx, Ncz )
  for ( k1=0; k1<Ncx*Ncz; k1++ ) {

    k      = mesh->kp[k1];
    l      = mesh->lp[k1];
    c0     = k  + l*(Ncx);

    mesh->eta_n[c0]       = 0.0;
    mesh->eta_phys_n[c0]  = 0.0;
    mesh->VE_n[c0]        = 0.0;
    mesh->sxxd[c0]        = 0.0;
    mesh->szzd[c0]        = 0.0;
    mesh->sxz_n[c0]       = 0.0;
    mesh->sxz_n[c0]       = 0.0;
    mesh->eII_el[c0]      = 0.0;
    mesh->eII_pl[c0]      = 0.0;
    mesh->eII_pwl[c0]     = 0.0;
    mesh->eII_exp[c0]     = 0.0;
    mesh->eII_lin[c0]     = 0.0;
    mesh->eII_gbs[c0]     = 0.0;
    mesh->eII_cst[c0]     = 0.0;
    mesh->d_n[c0]         = 0.0;
    mesh->p_corr[c0]      = 0.0;                // Achtung baby
    mesh->div_u_el[c0]    = 0.0;
    mesh->div_u_pl[c0]    = 0.0;
    mesh->div_u_r[c0]     = 0.0;
    mesh->Wtot[c0]        = 0.0;
    mesh->Wel[c0]         = 0.0;
    mesh->Wdiss[c0]       = 0.0;
    //        X                     =  mesh->Xreac_n[c0]; // Save X first
    //        if (model->chemical_diffusion==1) mesh->Xreac_n[c0]    = 0.0;
    mesh->X_n[c0]        = 0.0;
    mesh->OverS_n[c0]    = 0.0;

    if ( model->density_variations == 1 ) {
      mesh->rho_n[c0]  = 0.0;
    }

    // Loop on grid nodes
    if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31 ) {

      //----------------------------------------------------------//
      if ( model->elastic==1 ) eta_e      = model->dt*mesh->mu_n[c0];
      else                     eta_e      = 1.0; // set to arbitrary value to avoid division by 0.0
      //----------------------------------------------------------//
      // Anisotropy
      d1      = mesh->d1_n[c0];  // d1 = 2*lx^2*lz^2
      d2      = mesh->d2_n[c0];  // d2 = lx*lz*(pow(-lx, 2.0) + pow(lz, 2.0));
      angle   = mesh->angle_n[c0];
      lay_ang = angle - M_PI/2.0;
      lxlz    = cos(lay_ang)*sin(lay_ang);
      lx2     = cos(lay_ang)*cos(lay_ang);
      //----------------------------------------------------------//
      // if ( model->elastic==1 ) {
      //   Exx = mesh->exxd[c0]  + mesh->sxxd0[c0] /eta_e/2.0;
      //   Ezz = mesh->ezzd[c0]  + mesh->szzd0[c0] /eta_e/2.0;
      //   Exz = mesh->exz_n[c0] + mesh->sxz0_n[c0]/eta_e/2.0;
      // }
      // else {
      //   Exx = mesh->exxd[c0] ;
      //   Ezz = mesh->ezzd[c0] ;
      //   Exz = mesh->exz_n[c0];
      // }

      double aniS_e, aniS_vep;
      aniS_e = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
      EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], d1, d2, aniS_e, eta_e, model->elastic );

      // Loop on phases
      for ( p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(mesh->phase_perc_n[p][c0])>min_fraction ) is_phase_active = true;

        if ( is_phase_active ) {
          eta =  ViscosityConciseAniso( p, lxlz, lx2, angle, mesh->aniso_factor_n[c0], mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], Exx, Ezz, Exz, mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials    , model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], mesh->p0_n[c0], mesh->T0_n[c0], &Xreac, &OverS, &Pcorr, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 1, final_update, c0 );
          // eta =  ViscosityConcise( p, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi0_n[c0], mesh->X0_n[c0], Exx, Ezz, Exz, mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_n[c0], mesh->dil_n[c0], mesh->fric_n[c0], mesh->C_n[c0], mesh->p0_n[c0], mesh->T0_n[c0], &Xreac, &OverS, &Pcorr, &rho, mesh->bet_n[c0], mesh->div_u[c0], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 1, final_update, c0 );

          mesh->phase_eta_n[p][c0] = eta_vep;

          switch ( average ) {
            case 0 :
              // ARITHMETIC AVERAGE
              mesh->eta_n[c0]            += mesh->phase_perc_n[p][c0] * eta_vep;
              mesh->eta_phys_n[c0]       += mesh->phase_perc_n[p][c0] * eta;
              break;
            case 1 :
              // HARMONIC AVERAGE
              mesh->eta_n[c0]            += mesh->phase_perc_n[p][c0] * 1.0/eta_vep;
              mesh->eta_phys_n[c0]       += mesh->phase_perc_n[p][c0] * 1.0/eta;
              break;
            case 2 :
              // GEOMETRIC AVERAGE
              mesh->eta_n[c0]            += mesh->phase_perc_n[p][c0] * log(eta_vep);
              mesh->eta_phys_n[c0]       += mesh->phase_perc_n[p][c0] * log(eta);
              break;
          }

          mesh->VE_n[c0]       += mesh->phase_perc_n[p][c0] * eta_vep/eta_e;
          mesh->eII_el[c0]     += mesh->phase_perc_n[p][c0] * eII_el;
          mesh->eII_pl[c0]     += mesh->phase_perc_n[p][c0] * eII_pl;
          mesh->eII_pwl[c0]    += mesh->phase_perc_n[p][c0] * eII_pwl;
          mesh->eII_exp[c0]    += mesh->phase_perc_n[p][c0] * eII_exp;
          mesh->eII_lin[c0]    += mesh->phase_perc_n[p][c0] * eII_lin;
          mesh->eII_gbs[c0]    += mesh->phase_perc_n[p][c0] * eII_gbs;
          mesh->eII_cst[c0]    += mesh->phase_perc_n[p][c0] * eII_cst;
          mesh->d_n[c0]        += mesh->phase_perc_n[p][c0] * 1.0/dnew;

          mesh->Wtot[c0]       += mesh->phase_perc_n[p][c0] * Wtot;
          mesh->Wdiss[c0]      += mesh->phase_perc_n[p][c0] * Wdiss;
          mesh->Wel[c0]        += mesh->phase_perc_n[p][c0] * Wel;

          if (mesh->Wdiss[c0]<0.0) {LOG_INFO("negative dissipation: you crazy! --> Wdiss = %2.2e", mesh->Wdiss[c0]*scaling->S*scaling->E); }

          mesh->p_corr[c0]      += mesh->phase_perc_n[p][c0] * Pcorr;
          mesh->div_u_el[c0]    += mesh->phase_perc_n[p][c0] * div_el;
          mesh->div_u_pl[c0]    += mesh->phase_perc_n[p][c0] * div_pl;
          mesh->div_u_r[c0]     += mesh->phase_perc_n[p][c0] * div_r;

          mesh->X_n[c0]         += mesh->phase_perc_n[p][c0] * Xreac;
          mesh->OverS_n[c0]     += mesh->phase_perc_n[p][c0] * OverS;

          // Volume changes
          if ( model->density_variations == 1 ) {
            mesh->rho_n[c0]     += mesh->phase_perc_n[p][c0] * (rho);
          }
        }
      }

      mesh->d_n[c0]          = 1.0/mesh->d_n[c0];

      // HARMONIC AVERAGE
      if ( average == 1 ) {
        mesh->eta_n[c0]            = 1.0/mesh->eta_n[c0];
        mesh->eta_phys_n[c0]       = 1.0/mesh->eta_phys_n[c0];
        // mesh->aniso_factor_n[c0]   = 1.0/mesh->aniso_factor_n[c0];

        if (isinf (mesh->eta_phys_n[c0]) ) {
          LOG_INFO("Inf: Problem on cell centers:");
          for ( p=0; p<model->Nb_phases; p++) LOG_INFO("phase %d vol=%2.2e", p, mesh->phase_perc_n[p][c0]);
          LOG_INFO("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e %2.2e", eta, mesh->mu_n[c0], mesh->T[c0], mesh->p_in[c0], mesh->d0_n[c0], mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
          LOG_INFO("flag %d nb part cell = %d cell index = %d", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
          LOG_INFO("x=%2.2e z=%2.2e", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
          exit(1);
        }
        if (isnan (mesh->eta_phys_n[c0]) ) {
          LOG_INFO("NaN: Problem on cell centers:");
          LOG_INFO("chemical_diffusion %d", model->chemical_diffusion);
          for ( p=0; p<model->Nb_phases; p++) LOG_INFO("phase %d vol=%2.2e", p, mesh->phase_perc_n[p][c0]);
          LOG_INFO("eta=%2.2e G=%2.2e T=%2.2e P=%2.2e d=%2.2e phi=%2.2e %2.2e %2.2e %2.2e %2.2e", eta*scaling->eta, mesh->mu_n[c0]*scaling->S, mesh->T[c0]*scaling->T, mesh->p_in[c0]*scaling->S, mesh->d0_n[c0]*scaling->L, mesh->phi_n[c0], mesh->exxd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->sxz0_n[c0]);
          LOG_INFO("flag %d nb part cell = %d cell index = %d", mesh->BCp.type[c0],mesh->nb_part_cell[c0], c0);
          LOG_INFO("x=%2.2e z=%2.2e", mesh->xc_coord[k]*scaling->L/1000.0, mesh->zc_coord[l]*scaling->L/1000.0);
          exit(1);
        }
      }
      // GEOMETRIC AVERAGE
      if ( average == 2 ) {
        mesh->eta_n[c0]            = exp(mesh->eta_n[c0]);
        mesh->eta_phys_n[c0]       = exp(mesh->eta_phys_n[c0]);
      }

      // Effective strain rate
      // double aniS_e, aniS_vep;
      // if ( model->aniso_fstrain  == 0 ) aniS_e = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
      // if ( model->aniso_fstrain  == 1 ) aniS_e = 1.0 - 1.0 / mesh->FS_AR_n[c0];
      // EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd[c0], mesh->ezzd[c0], mesh->exz_n[c0], mesh->sxxd0[c0], mesh->szzd0[c0], mesh->sxz0_n[c0], d1, d2, aniS_e, eta_e, model->elastic );
      // Final stress update
      aniS_vep = 1.0 - 1.0 / mesh->aniso_factor_n[c0];
      // Normal stress
      const double Dxx = Exx*(1.0 - aniS_vep*d1) + Ezz*aniS_vep*d1 - 2.0*Exz*aniS_vep*d2;
      const double Dzz = Ezz*(1.0 - aniS_vep*d1) + Exx*aniS_vep*d1 + 2.0*Exz*aniS_vep*d2;
      const double Dxz = -Exx*aniS_vep*d2 + Ezz*aniS_vep*d2 + 2*Exz*(aniS_vep*(d1 - 0.5) + 0.5);
      mesh->sxxd[c0]  = 2.0 * mesh->eta_n[c0] * Dxx;
      mesh->szzd[c0]  = 2.0 * mesh->eta_n[c0] * Dzz;
      mesh->sxz_n[c0] = 2.0 * mesh->eta_n[c0] * Dxz;
    }
  }

// Calculate vertices viscosity
#pragma omp parallel for shared( mesh ) private( k, l, k1, p, eta, c1, c0, txx1, tzz1, txz1, eta_vep, ani_vep, ani_e, eII_el, eII_pl, eII_pwl, eII_exp, eII_lin, eII_gbs, eII_cst, dnew, Wtot, Wel, Wdiss, Xreac, OverS, Pcorr, rho, div_el, div_pl, div_r, Exx, Ezz, Exz, eta_e, ani, d1, d2, Da11, Da12, Da13, Da22, Da23, Da33, iDa11, iDa12, iDa13, iDa22, iDa23, iDa33, a11, a12, a13, a22, a23, a33, det, angle, lxlz, lx2, lay_ang ) firstprivate( el, materials, scaling, average, model, Nx, Nz )
  for ( k1=0; k1<Nx*Nz; k1++ ) {

    k  = mesh->kn[k1];
    l  = mesh->ln[k1];
    c1 = k + l*Nx;

    mesh->sxxd_s[c1]           = 0.0;
    mesh->szzd_s[c1]           = 0.0;
    mesh->sxz[c1]              = 0.0;
    mesh->eta_phys_s[c1]       = 0.0;
    mesh->eta_s[c1]            = 0.0;
    mesh->VE_s[c1]             = 0.0;

    mesh->X_s[c1]        = 0.0;

    if ( mesh->BCg.type[c1] != 30 ) {

      if ( model->elastic==1   ) eta_e      = model->dt*mesh->mu_s[c1];
      else           eta_e      = 1.0; // set to arbitrary value to avoid division by 0.0
      //----------------------------------------------------------//
      // Anisotropy
      d1      = mesh->d1_s[c1];
      d2      = mesh->d2_s[c1];
      angle   = mesh->angle_s[c1];
      lay_ang = angle - M_PI/2.0;
      lxlz    = cos(lay_ang)*sin(lay_ang);
      lx2     = cos(lay_ang)*cos(lay_ang);
      //----------------------------------------------------------//
      // Effective strain rate
      double aniS_e, aniS_vep;
      aniS_e = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
      EffectiveStrainRate( &Exx, &Ezz, &Exz, mesh->exxd_s[c1], mesh->ezzd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], d1, d2, aniS_e, eta_e, model->elastic );

      // if ( model->elastic==1   ) {
      //   Exx = mesh->exxd_s[c1] + mesh->sxxd0_s[c1]/eta_e/2.0;
      //   Ezz = mesh->ezzd_s[c1] + mesh->szzd0_s[c1]/eta_e/2.0;
      //   Exz = mesh->exz[c1]    + mesh->sxz0[c1]   /eta_e/2.0;
      // }
      // else {
      //   Exx = mesh->exxd_s[c1];
      //   Ezz = mesh->ezzd_s[c1];
      //   Exz = mesh->exz[c1]   ;
      // }

      // Loop on phases
      for ( p=0; p<model->Nb_phases; p++) {

        // Detect if there is a fraction of phase p in the cell c: compute only if there is a non-zero fraction
        bool is_phase_active = false;
        const double min_fraction=1e-13;
        if ( fabs(mesh->phase_perc_s[p][c1])>min_fraction ) is_phase_active = true;

        if ( is_phase_active ) {
          eta =  ViscosityConciseAniso( p, lxlz, lx2, angle,  mesh->aniso_factor_s[c1], mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], Exx, Ezz, Exz, mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], mesh->p0_s[c1], 0.0, &Xreac, &OverS, &Pcorr, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 0, final_update, c1 );

          // if (c1==65) {
          //   printf("final %1.4e\n", eta_vep);
          // }

          // eta =  ViscosityConcise( p, mesh->mu_s[c1], mesh->T_s[c1], mesh->P_s[c1], mesh->d0_s[c1], mesh->phi0_s[c1], mesh->X0_s[c1], Exx, Ezz, Exz, mesh->sxxd0_s[c1], mesh->szzd0_s[c1], mesh->sxz0[c1], materials, model, scaling, &txx1, &tzz1, &txz1, &eta_vep, &VEcoeff, &eII_el, &eII_pl, &eII_pwl, &eII_exp, &eII_lin, &eII_gbs, &eII_cst, &dnew, mesh->strain_s[c1], mesh->dil_s[c1], mesh->fric_s[c1], mesh->C_s[c1], mesh->p0_s[c1], 0.0, &Xreac, &OverS, &Pcorr, &rho, mesh->bet_s[c1], mesh->div_u_s[c1], &div_el, &div_pl, &div_r, &Wtot, &Wel, &Wdiss, 1, 0, final_update, c1 );


          // if (c1==65) {
          //   printf("final %1.4e\n", eta_vep);
          // }

          mesh->phase_eta_s[p][c1] = eta_vep;

          switch ( average ) {
            case 0 :
              mesh->eta_s[c1]            += mesh->phase_perc_s[p][c1] * eta_vep;
              mesh->eta_phys_s[c1]       += mesh->phase_perc_s[p][c1] * eta;
              break;
            case 1:
              mesh->eta_s[c1]            += mesh->phase_perc_s[p][c1] * 1.0/eta_vep;
              mesh->eta_phys_s[c1]       += mesh->phase_perc_s[p][c1] * 1.0/eta;
              break;
            case 2:
              mesh->eta_s[c1]            += mesh->phase_perc_s[p][c1] * log(eta_vep);
              mesh->eta_phys_s[c1]       += mesh->phase_perc_s[p][c1] * log(eta);
              break;
          }
          mesh->VE_s[c1]                            += mesh->phase_perc_s[p][c1] * eta_vep/eta_e;
          mesh->X_s[c1] += mesh->phase_perc_s[p][c1] * Xreac;
        }
      }
      // HARMONIC AVERAGE
      if (average == 1) {
        mesh->eta_s[c1]            = 1.0/mesh->eta_s[c1];
        mesh->eta_phys_s[c1]       = 1.0/mesh->eta_phys_s[c1];
        if (isinf (mesh->eta_phys_s[c1]) ) {
          LOG_INFO("Inf: Problem on cell vertices:");
          for ( p=0; p<model->Nb_phases; p++) LOG_INFO("phase %d vol=%2.2e", p, mesh->phase_perc_s[p][c1]);
          LOG_INFO("%2.2e %2.2e %2.2e %2.2e %2.2e ", mesh->mu_s[c1], mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
          LOG_INFO("x=%2.2e z=%2.2e", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
          exit(1);
        }
        if (isnan (mesh->eta_phys_s[c1]) ) {
          LOG_INFO("Nan: Problem on cell vertices:");
          for ( p=0; p<model->Nb_phases; p++) LOG_INFO("phase %d vol=%2.2e", p, mesh->phase_perc_s[p][c1]);
          LOG_INFO("%2.2e %2.2e %2.2e %2.2e %2.2e ", mesh->mu_s[c1],  mesh->exxd_s[c1], mesh->exz[c1], mesh->sxxd0_s[c1], mesh->sxz0[c1]);
          LOG_INFO("x=%2.2e z=%2.2e", mesh->xg_coord[k]*scaling->L/1000, mesh->zg_coord[l]*scaling->L/1000);
          exit(1);
        }
      }
      // GEOMETRIC AVERAGE
      if (average == 2) {
        mesh->eta_s[c1]            = exp(mesh->eta_s[c1]);
        mesh->eta_phys_s[c1]       = exp(mesh->eta_phys_s[c1]);
      }

      // Final stress update
      aniS_vep = 1.0 - 1.0 / mesh->aniso_factor_s[c1];
      const double Dxx = Exx*(1.0 - aniS_vep*d1) + Ezz*aniS_vep*d1 - 2.0*Exz*aniS_vep*d2;
      const double Dzz = Ezz*(1.0 - aniS_vep*d1) + Exx*aniS_vep*d1 + 2.0*Exz*aniS_vep*d2;
      const double Dxz  = -Exx*aniS_vep*d2 + Ezz*aniS_vep*d2 + 2*Exz*(aniS_vep*(d1 - 0.5) + 0.5);
      // Shear stress
      mesh->sxxd_s[c1] = 2.0*mesh->eta_s[c1]*Dxx;
      mesh->szzd_s[c1] = 2.0*mesh->eta_s[c1]*Dzz;
      mesh->sxz[c1]    = 2.0*mesh->eta_s[c1]*Dxz;
    }
  }

  // User-specified grid-density override (no-op if SetGridDensity callback is NULL)
  ApplySetGridDensityOverride(mesh, model);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseDirectorVector (grid* mesh, markers* particles, params* model, mat_prop* materials, double scaling_L) {

  int    cent=1, vert=0, prop=1, interp=0;
  int k;
  double angle, norm;
  const double Earth_radius = 6370e3/scaling_L;

// #pragma omp parallel for shared( particles ) private( angle, norm )
  for (k=0; k<particles->Nb_part; k++) {

    if ( particles->phase[k] != -1 ) {

      // Set up director vector
      if (model->marker_aniso_angle) {
        angle           = particles->aniso_angle[k];
      } else {
        angle           = materials->aniso_angle[particles->phase[k]];
      }
      if (model->polar == 1 ) {
        angle = angle - atan(particles->x[k] / (particles->z[k] ));
      }
      particles->nx[k]  = cos(angle);
      particles->nz[k]  = sin(angle);
      norm              = sqrt(particles->nx[k]*particles->nx[k] + particles->nz[k]*particles->nz[k]);
      particles->nx[k] /= norm;
      particles->nz[k] /= norm;
    }
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void NormalizeDirector( grid* mesh, DoodzFP* nx_n, DoodzFP* nz_n, DoodzFP* nx_s, DoodzFP* nz_s, params *model ) {

  int Nx, Nz, k;
  double nx, nz, norm;
  Nx = model->Nx;
  Nz = model->Nz;

#pragma omp parallel for shared ( mesh, nx_n, nz_n ) \
private ( k, nx, nz, norm )              \
firstprivate( model )
  for ( k=0; k<(Nx-1)*(Nz-1); k++ ) {
    if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) {
      nx            = nx_n[k];
      nz            = nz_n[k];
      norm          = sqrt(nx*nx + nz*nz);
      nx            = nx/norm;
      nz            = nz/norm;
      nx_n[k]       = nx;
      nz_n[k]       = nz;
    }
  }

#pragma omp parallel for shared ( mesh, nx_s, nz_s ) \
private ( k, nx, nz, norm )              \
firstprivate( model )
  for ( k=0; k<Nx*Nz; k++ ) {
    if (mesh->BCg.type[k] != 30) {
      nx            = nx_s[k];
      nz            = nz_s[k];
      norm          = sqrt(nx*nx + nz*nz);
      nx            = nx/norm;
      nz            = nz/norm;
      nx_s[k]       = nx;
      nz_s[k]       = nz;
    }
  }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
