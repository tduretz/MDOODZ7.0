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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PartialMelting( double *phi, double* hlat, double P, double T, double phi0, double tk, double dt, int melt_model, scale *scaling ) {
    
    // ---- ACHTUNG: computations BELOW are made in dimensional form (MPa and K) to avoid conversion errors  
    *phi  = 0.0; *hlat = 0.0;
    double phi_eq;
    const double P_MPa = P/1e6*scaling->S; // convert P to MPa
    const double T_K   = T*scaling->T;     // convert T to K
    double Ts = T_K, Tl = T_K, Ql;
    double alpha = 1e-7/(1/scaling->T); // very small number -> linear melting

    if (melt_model>0 && T>1e-13 && P>1e-13) {

        // Simple T-dependent melting (Stuewe 1995)
        if (melt_model==1) {
            Ts = 600.  + 273.15;
            Tl = 1200. + 273.15;
            Ql = 320000.;
            alpha = -0.1; // Strong effect
        }

        // Crust or sediment melting
        if (melt_model==10) {
            if (P_MPa>1600) {
                Ts = 935. + 0.0035*P_MPa + 0.0000062*pow(P_MPa,2);
            }
            else {
                Ts = 973. - 70400/(P_MPa + 354) + 77800000/pow((P_MPa + 354),2);
            }
            Tl = 1423 + 0.105*P_MPa;
            Ql = 380000.;
        }

        // From Schenker et al, 2012
        if (melt_model==11) {
            if (P_MPa>1200) {
                Ts = 831. + 0.06*P_MPa;
            }
            else {
                Ts = 889. - 17900/(P_MPa + 54) + 20,200/pow((P_MPa + 54),2);
            }
            Tl = 1262 + 0.09*P_MPa;
            Ql = 300000.;
        }
        
        // Mantle melting
        if (melt_model==40) {
            if (P_MPa > 10000) {
                Ts = 2212 + 0.030819*(P_MPa - 10000); //en K
            }
            else {
                Ts = 1394 + 0.132899*P_MPa - 0.000005014*pow(P_MPa,2.0); // en K
            }
            Tl = 2073 + 0.114*P_MPa; // en K
            Ql = 400000.;
        }

        // From Schenker et al, 2012
        if (melt_model==41) {
            if (P_MPa>1e4) {
                Ts = 2212. + 0.030819*(P_MPa-10000.);
            }
            else {
                Ts = 1394 + 0.132899*P_MPa - 0.000005014/pow(P_MPa,2);
            }
            Tl = 2073 + 0.114*P_MPa;
            Ql = 300000.;
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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

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

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MassSourceTerm( grid* mesh, markers* particles, mat_prop *materials, params *model, scale *scaling ) {
    
    // alpha*rho0/rho*dT/dt
    int Nx, Nz, Ncx, Ncz;
    int c0, p;

    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        
        mesh->divth_n[c0]  = 0.0;

        for ( p=0; p<model->Nb_phases; p++) {
            
            if ( fabs(mesh->phase_perc_n[p][c0])>1.0e-13 ) {
                
                // Compute only if below free surface
                if ( mesh->BCp.type[0] != 30 && mesh->BCp.type[0] != 31) {

                    const double divth = materials->alp[p]*materials->rho[p]/mesh->rho_n[c0]*(mesh->T[c0] - mesh->T0_n[c0])/model->dt;
                    
                    // Update centroid values
                    mesh->divth_n[c0]  += mesh->phase_perc_n[p][c0]*divth;
                 }
            }
         }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
// This is the main function
void InjectDikesAndSills(markers *particles, MdoodzInput *instance){
    // Initialize and fix the semi-major and minor axis of the dikes, read in the rate and current time step to calcuted the magma budget
    double semi_maj = 2.5e3 / instance->scaling.L; // This times 2 is the length of the dike
    double semi_min = 1.0e3 / instance->scaling.L; // This times 2 is the width of the dike
    double dike_ext = 0.5e3 / instance->scaling.L; // Extend of the dike in the third dimension. Used to calculate the area magma budget below
    double dt       = instance->model.dt;          // Current time step - ACHTUNG this value is scaled
    double conv     = pow(1000.0, 3) / (3600.0*24.0*365.25); // Convert km3/yr to m3/s
    double Vdot     = instance->model.user4 * conv / (pow(instance->scaling.L, 3)) * instance->scaling.t;       // Fixed magma injection rate
    double Abud     = Vdot / dike_ext * instance->model.dt;             // This is the area magma budget we need to get provide at the current time step
    
    // Initialize central coords, angle of the dikes to be injected, and the current area budget - these are not fixed but change during dike injection
    double Lx         = instance->model.xmax - instance->model.xmin;
    double Lz         = instance->model.zmax - instance->model.zmin;
    double Ainj       = 0.0;
    double xinj_s     = instance->model.xmin + Lx / 5.0; // Vertex coordinates of injection area
    double xinj_e     = instance->model.xmax - Lx / 5.0;
    double zinj_s     = instance->model.zmin + Lz / 5.0;
    double zinj_e     = instance->model.zmax - Lz / 5.0;
    double xc         = 0.0;
    double zc         = 0.0;
    double dx_dike    = 5.0 * instance->model.dx;
    double dz_dike    = 5.0 * instance->model.dz;
    double angle_dike = 0.0;
    int    ph_dike    = (int)instance->model.user5;
    double T_dike     = (instance->model.user6 + zeroC) / instance->scaling.T;

    // Inject dikes until the area magma budget is met
    int iter = 0, max_iter = 1e3;
    double tol = 1e-8;
    while (Ainj < Abud)
    {
        // Update counter
        iter += 1;

        // Calculate coordinates of dike center
        xc = rand_between_float_asep(xinj_s + semi_min, xinj_e - semi_min, dx_dike, xc);
        zc = rand_between_float_asep(zinj_s + semi_maj, zinj_e - semi_maj, dz_dike, zc);

        // Depth dependent dike angle
        if (zc * instance->scaling.L < -10.0e3)
        {
            angle_dike = rand_between_float(1.0, 15.0) / 180.0 * M_PI;
        }
        else{
            angle_dike = rand_between_float(60.0, 120.0) / 180.0 * M_PI;
        }

        // Update injected dike area
        Ainj += M_PI * semi_maj * semi_min;

        // Assign new phase and temperature for particles within the injection zone
        int    ip;
        double xp, zp, xr, zr;
#pragma omp parallel for shared(particles, xc, zc, xinj_s, xinj_e, zinj_s, zinj_e, semi_maj, semi_min, angle_dike, ph_dike, T_dike) private(ip, xp, zp, xr, zr)
        for (ip = 0; ip < particles->Nb_part; ip++)
        {
            // Get local coords
            xp = particles->x[ip];
            zp = particles->z[ip];
            // Only compute if we are in the injection zone
            if ( 
                (xinj_s <= xp && xinj_e >= xp) && 
                (zinj_s <= zp && zinj_e >= zp)
            )
            {
                // Compute particle coordinate w.r.t. dike location
                xr =  (xp - xc) * cos(angle_dike) + (zp - zc) * sin(angle_dike);
                zr = -(xp - xc) * sin(angle_dike) + (zp - zc) * cos(angle_dike);

                // Reset phase and temperature for dike particle
                if ( (xr*xr)/(semi_min*semi_min) + (zr*zr)/(semi_maj*semi_maj) <= 1.0 )
                {
                    particles->phase[ip] = ph_dike;
                    particles->T[ip]     = T_dike;
                }
            }
            
        }
        
        // Break criterion
        if (iter > max_iter)
        {
            printf(RED "WARNING: Maximum iteration reached befor the entire magma budget has been used\n");
            printf(RED "Area budget: \t %.2f m2 \t Area injected: \t %.2f m2\n", Abud * pow(instance->scaling.L, 2), Ainj * pow(instance->scaling.L, 2));
            break;
        }
    }
}

// This calculates a random number between two floats
double rand_between_float(double min, double max){
    return (min + ((double)rand() + 0.5) / ((double)RAND_MAX + 1.0) * (max - min)); // Adding 0.5 and 1.0 avoids the first value to be exactly min value
}

// This calculates a random number between two floats using an absolute separator
double rand_between_float_asep(double min, double max, double sep, double prev){
    double curr  = prev;    // Set initial guess for current value
    int max_iter = 1e6;     // Maximum iteration limit
    int count    = 0;       // counter

    // Repeat computation until the desired distance is met
    while (fabs(prev - curr) < sep)
    {
        // Increment counter
        count += 1;

        // Compute new value
        curr = min + ((double)rand() + 0.5) / ((double)RAND_MAX + 1.0) * (max - min);

        // Break condition
        if (count > max_iter)
        {
            printf(RED "WARNING: Maximum iterations for rejection sampling reached. Previous value: \t %.4f | Current value: \t %.4f\n", prev, curr);
            break;
        }
    }
    
    return curr; // Adding 0.5 and 1.0 avoids the first value to be exactly min value
}
/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/