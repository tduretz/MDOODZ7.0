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

#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"
#include "cholmod.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

//---- M-Doodz header file
#include "mdoodz-private.h"

void PartialMelting( double *phi, double* hlat, double P, double T, double phi0, double tk, double dt, int melt_model, scale *scaling ) {
    
    // ---- ACHTUNG: computations BELOW are made in dimensional form (MPa and K) to avoid conversion errors  
    *phi  = 0.0; *hlat = 0.0;
    double phi_eq;
    const double P_MPa = P/1e6*scaling->S; // convert P to MPa
    const double T_K   = T*scaling->T;     // convert T to K
    double Ts = T_K, Tl = T_K, Ql;
    double alpha = 1e-7/(1/scaling->T); // very small number -> linear melting

    if (melt_model>0) {

        // Crust or sediment melting
        if (melt_model==1) {
            if (P_MPa>1600) {
                Ts = 935. + 0.0035*P_MPa + 0.0000062*pow(P_MPa,2);
            }
            else {
                Ts = 973. - 70400/(P_MPa + 354) + 77800000/pow((P_MPa + 354),2);
            }
            Tl = 1423 + 0.105*P_MPa;
            Ql = 380000.;
        }
        
        // Mantle melting
        if (melt_model==2) {
            if (P_MPa > 10000) {
                Ts = 2212 + 0.030819*(P_MPa - 10000); //en K
            }
            else {
                Ts = 1394 + 0.132899*P_MPa - 0.000005014*pow(P_MPa,2.0); // en K
            }
            Tl = 2073 + 0.114*P_MPa; // en K
            Ql = 400000.;
        }

        // Simple T-dependent melting (Stuewe 1995)
        if (melt_model==3) {
            Ts = 600.  + 273.15;
            Tl = 1200. + 273.15;
            Ql = 320000.;
            alpha = -0.1; // Strong effect
        }
        
        Ql /= (scaling->J/scaling->m);

        // ---- ACHTUNG: computations ABOVE are made in dimensional form (MPa and K) to avoid conversion errors  

        // Computation of melt fraction and laten heat release/consumption
        *phi  = 0.0;
        *hlat = 0.0;
        if (Tl>0.0) {
            // Solidus and liquidus must not intersect in the extrapolation region
            if (Ts > Tl - 100.0) {
                Ts = Tl - 100.0;
            }
            // Melt fraction
            // phi_eq = (T_K - Ts) / (Tl - Ts);
            phi_eq = (exp(alpha*T_K) - exp(alpha*Ts)) / (exp(alpha*Tl) - exp(alpha*Ts));
            *phi   = dt/(tk + dt) * phi_eq + tk/(tk + dt) * phi0;
            if (*phi<0.0) *phi = 0.0;
            if (*phi>1.0) *phi = 1.0;
            // Latent heat
            *hlat = Ql;
        }
    }
}

void MeltFractionGrid( grid* mesh, markers* particles, mat_prop *materials, params *model, scale *scaling ) {
    
    int Nx, Nz, Ncx, Ncz;
    int c0, p;
    double Ql;
    double phi;
    
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        mesh->phi_n[c0] = 0.0;
        for ( p=0; p<model->Nb_phases; p++) {
            if ( fabs(mesh->phase_perc_n[p][c0])>1.0e-13 ) {
                const int melt_model = materials->melt[p];
                const double tk = materials->tau_kin[p];
                PartialMelting( &phi, &Ql, mesh->p_in[c0], mesh->T[c0], mesh->phi0_n[c0], tk, model->dt, melt_model, scaling );
                mesh->phi_n[c0] += mesh->phase_perc_n[p][c0]*phi;
            }
        }
    }
}


void UpdateAlphaCp( grid* mesh, markers* particles, mat_prop *materials, params *model, scale *scaling ) {
    
    // alpha effective =  alpha - rho/T*Ql*dphidP
    // Cp effective    =  Cp    +       Ql*dphidT
    int Nx, Nz, Ncx, Ncz;
    int c0, p;
    double dMdP, dMdT, Ql;
    double phi_plus, phi_minus, phi, pert;
    double Cp_melt = 0., alp_melt = 0.;
    double Ql_plus = 0., Ql_minus = 0.;

    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        
        mesh->Cp[c0]  = 0.0;
        mesh->alp[c0] = 0.0;

        for ( p=0; p<model->Nb_phases; p++) {
            
            if ( fabs(mesh->phase_perc_n[p][c0])>1.0e-13 ) {
                
                // Compute only if below free surface
                if ( mesh->BCp.type[0] != 30 && mesh->BCp.type[0] != 31) {
                    const int melt_model = materials->melt[p];
                    const double tk = materials->tau_kin[p];
                    PartialMelting( &phi, &Ql, mesh->p_in[c0], mesh->T[c0], mesh->phi0_n[c0], tk, model->dt, melt_model, scaling );  
                    Cp_melt  = 0.0;
                    alp_melt = 0.0;

                    if ((phi)>1.0e-13 ) {
                       pert = mesh->p_in[c0]*1e-6;
                       PartialMelting( &phi_minus, &Ql_minus, mesh->p_in[c0]-pert, mesh->T[c0], mesh->phi0_n[c0], tk, model->dt, melt_model, scaling );
                       PartialMelting( &phi_plus,  &Ql_plus,  mesh->p_in[c0]+pert, mesh->T[c0], mesh->phi0_n[c0], tk, model->dt, melt_model, scaling );
                       dMdP = (phi_plus - phi_minus)/(2.0*pert);
                       pert = mesh->T[c0]*1e-6;
                       PartialMelting( &phi_minus, &Ql_minus, mesh->p_in[c0], mesh->T[c0]-pert, mesh->phi0_n[c0], tk, model->dt, melt_model, scaling );
                       PartialMelting( &phi_plus,  &Ql_plus,  mesh->p_in[c0], mesh->T[c0]+pert, mesh->phi0_n[c0], tk, model->dt, melt_model, scaling );
                       dMdT = (phi_plus - phi_minus)/(2.0*pert);
    
                        // const doucpble dAdM = Ql*mesh->rho_n[c0]/mesh->T[c0]; 
                        // const double dAdT = dAdM*dMdT;
                        // Cp_melt  =  mesh->T[c0]/mesh->rho_n[c0]*dAdT; // equivalent to Taras
                        // alp_melt =  dAdT*dMdP;                        // equivalent to Taras
                        Cp_melt  = Ql*dMdT;
                        alp_melt = Ql*dMdP/mesh->T[c0]*mesh->rho_n[c0]; 
                    }

                    // Update centroid values of Cp and alpha
                    mesh->Cp[c0]  += (mesh->phase_perc_n[p][c0]*(materials->Cp[p]  + Cp_melt) );
                    mesh->alp[c0] += (mesh->phase_perc_n[p][c0]*(materials->alp[p] - alp_melt) );
                }
            }
        }
    }
}
