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
#include "header_MDOODZ.h"
#include "ctype.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataPowerLaw( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            printf ("should not be here\n");
            success      = 0;
            break;
            
        /****************** SPECIAL CASE: user-defined power law flow law ******************/
            
        case 1 :
            printf("'Homemade' power law flow of the form: eta = eta0 * exp(-Q/n/R/T) * Eii^(1/n - 1):\n" );
            mat->tpwl[k] = 0.0;
            mat->npwl[k] = mat->npwl[k];
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = mat->Qpwl[k];
            mat->Vpwl[k] = 0.0;
            mat->Apwl[k] = pow(mat->eta0[k], -mat->npwl[k]);
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

        case 2 :
            printf("Maryland Diabase - Mackwell et al. (1998) Homemade Stronger:\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 485.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 0.11111111*5.0477e-28;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

        case 3 :
            printf("Maryland Diabase - Mackwell et al. (1998) Homemade Weaker:\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 485.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 0.05*5.0477e-28;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        /******************************** Crust flow laws *********************************/
                
        case 10:
            printf("Westerly Granite (dry) - Hansen & Carter (1983):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.3;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 186.5e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.1623e-26;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 11:
            printf("Maryland Diabase - Mackwell et al. (1998):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 485.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 5.0477e-28;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

        case 12:
            printf("Maryland Diabase - Carter & Tsenn (1987):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 276.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.2e-20;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 13:
            printf("Wet Quartzite - Ranalli (1995):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 2.3;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 154.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 5.0717e-18;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 14:
            printf("Rocksalt - Ranalli (1995):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 5.3;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 102.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 9.9848e-32;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 15:
            printf("Calcite - Renner et al. (2002):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 297.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.5849e-25;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 16:
            printf("Felsic granulite - Ranalli (1995):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.1;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 243.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 2.0095e-21;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 17:
            printf("Mafic granulite - Ranalli (1995):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.2;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 445.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 8.8334e-22;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 18:
            printf("Calcite - Schmid et al. (1977):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 297.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.5849e-25;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 19:
            printf("Anorthite 100 - 0.007 H20 Wt%% - Rybacki & Dresen (2000, 2004):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 356.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.9811e-16;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
         
        case 20:
            printf("Garnet - Ji and Martignole (1994):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 2.22;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 485.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 2.7952e-07;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

        case 21:
            printf("Omphacite - Zhang et al. (2006):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.5;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 310.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.0e-23;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 22:
            printf("Qtz - Koch et al. (1989):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 2.72;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 134.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 7.3191e-24;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 23:
            printf("Fsp - Shelton and Tullis (1981):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.2;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 238.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 2.0632e-23;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 24:
            printf("Mica - Kronenberg et al. (1990):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 18.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 51.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.0000e-138;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 25:
            printf("Plagioclase - Shelton (1981) - Ph.D. Thesis:\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 431.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 7.9433e-24;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
        
        case 26:
            printf("Wet Quartzite - Burov (2011), Brace and Kohlstedt (1980), Kirby and Kronenberg (1987), Kohlstedt et al. (1995):\n" );
            mat->tpwl[k] = 1; //1           // NEED TO CHECKED!!!
            mat->npwl[k] = 2.4;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 160.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.981e-19;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 27:
            printf("Dry Quartz - Ranalli et al. (1995, 1997), Ranalli (2003):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 2.4;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 156.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 2.667e-20;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
            
        case 28:
            printf("Qtz - Hirth et al. (2001) - uniaxial:\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 135.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 6.3096e-36;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 29 :
            printf("Basalt - Hacker & Chritie (1990):\n" );
            mat->tpwl[k] = 0.0;
            mat->npwl[k] = 3.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 244e3;
            mat->Vpwl[k] = 0.0;
            mat->Apwl[k] = 6.9405e-27;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        // Granulite Eclogite Setups:
        case 30:
            printf("Anorthite dry - Rybacki and Dresen (2000):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 648.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 5.0119e-06;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 31:
            printf("DRY CPX - Bystricky & Mackwell, 2001:\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.7;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 760.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.9811e-19;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 32:
            printf("Anorthite 60 - Rybacki & Dresen (2004):\n" ); // added by Pauline
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 235.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.16228e-20;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        /******************************** Mantle flow laws ********************************/
            
        case 40:
            printf("Olivine Dry Dislocation creep - Hirth & Kohlstedt (2003):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.5;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 530.0e3;
            mat->Vpwl[k] = 11.0e-6;
            if ( model->force_act_vol_ast == 1 ) mat->Vpwl[k] = model->act_vol_dis_ast;
            mat->Apwl[k] = 1.1000e-16;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 41:
            printf("Olivine Wet Dislocation creep - Hirth & Kohlstedt (2003):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.5;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 1.2;
            mat->Qpwl[k] = 480.0e3;
            mat->Vpwl[k] = 11.0e-6;
            mat->Apwl[k] = 5.6786e-27;
            mat->fpwl[k] = 1.0e9;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 42:
            printf("Dry Olivine average from Ranalli (1995) - Kirby & Kronenberg (1987):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.5;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 532.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 2.5119e-17;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 43:
            printf("Wet Olivine average from Ranalli (1995) - Kirby & Kronenberg (1987):\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 4.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 471.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.9953e-21;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 44:
            printf("Olivine - Chopra and Paterson 1984 :\n" );
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.0;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 520.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.0e-14;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 45:
            printf("Olivine Dry Dislocation creep - Hirth & Kohlstedt (2003):\n" ); // added by Pauline
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.5;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 520.0e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 1.60e-18;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
            
        case 46:
            printf("Serpentine 1 GPa - Hilairet et al. (2007):\n" ); // added by Pauline
            mat->tpwl[k] = 1;
            mat->npwl[k] = 5.8;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 17.6e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 5.3443e-41;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

    }
    
    // Scaling
    mat->tpwl[k]   /= 1.0;
    mat->npwl[k]   /= 1.0;
    mat->mpwl[k]   /= 1.0;
    mat->rpwl[k]   /= 1.0;
    mat->Qpwl[k]   /= scaling->J;
    mat->Vpwl[k]   /= pow(scaling->L,3);
    mat->Apwl[k]   /= (pow(scaling->S, -mat->npwl[k]) * pow(scaling->S, -mat->rpwl[k]) * pow(scaling->L, mat->mpwl[k]) / scaling->t);
    mat->fpwl[k]   /= scaling->S;
    mat->apwl[k]   /= 1.0;
    
    // Calculate correction factor for invariant formulation 
    if (mat->tpwl[k]==0) mat->Fpwl[k]   = 1.0;                                                                              // no assumed experimental correction
    if (mat->tpwl[k]==1) mat->Fpwl[k]   = 1.0/6.0*pow(2.0,1.0/mat->npwl[k]) * pow(3.0,(mat->npwl[k]-1.0)/2.0/mat->npwl[k]); // axial compression
    if (mat->tpwl[k]==2) mat->Fpwl[k]   = 1.0/4*pow(2.0,1.0/mat->npwl[k]);                                                    // simple shear
    
    // Cancel correction factor
    if (number<0) mat->Fpwl[k]   = 1.0;
    
    if ( success==0 ) { printf("Error: Non existing Power Law flow law number\n"); exit(12);}
    
    // Print to screen 
    printf("t = %1.0lf  n = %1.1lf  m = %1.1lf  r = %1.1lf  Q = %2.2e J  V = %2.2e m^3  A = %2.2e Pa^-n/s  f = %2.2e Pa  a = %1.1lf  F = %1.2lf\n", mat->tpwl[k], mat->npwl[k], mat->mpwl[k], mat->rpwl[k], mat->Qpwl[k]*scaling->J, mat->Vpwl[k]*pow(scaling->L,3.0), mat->Apwl[k]*(pow(scaling->S, -mat->npwl[k]) * pow(scaling->L, mat->mpwl[k]) / scaling->t), mat->fpwl[k]*scaling->S, mat->apwl[k], mat->Fpwl[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataLinear( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            printf ("should not be here\n");
            success      = 0;
            break;
            
        case 15 :
            printf("Calcite Diffusion creep - Herwegh (2003):\n");
            mat->tlin[k] = 1.0;
            mat->nlin[k] = 1.1;
            mat->mlin[k] = 3.3;
            mat->rlin[k] = 0.0;
            mat->Qlin[k] = 200e3;
            mat->Vlin[k] = 0e-6;
            mat->Alin[k] = 1.7119e-19;
            mat->flin[k] = 0.0;
            mat->alin[k] = 0.0;
            success      = 1;
            break;
            
        case 19 :
            printf("Anorthite 100 - 0.007 H20 Wt%% - Rybacki & Dresen (2000, 2004):\n");
            mat->tlin[k] = 1.0;
            mat->nlin[k] = 1.0;
            mat->mlin[k] = 3.0;
            mat->rlin[k] = 0.0;
            mat->Qlin[k] = 170e3;
            mat->Vlin[k] = 0e-6;
            mat->Alin[k] = 5.0119e-23;
            mat->flin[k] = 0.0;
            mat->alin[k] = 0.0;
            success      = 1;
            break;
            
        case 25 :
            printf("Plagioclase Diffusion creep - Rybacki and Dresen (2000):\n");
            mat->tlin[k] = 1.0;
            mat->nlin[k] = 1.0;
            mat->mlin[k] = 3.0;
            mat->rlin[k] = 0.0;
            mat->Qlin[k] = 268e3;
            mat->Vlin[k] = 0e-6;
            mat->Alin[k] = 1.2589e-19;
            mat->flin[k] = 0.0;
            mat->alin[k] = 0.0;
            success      = 1;
            break;
            
        case 40 :
            printf("Olivine Dry Diffusion creep - Hirth & Kohlstedt (2003):\n");
            mat->tlin[k] = 1;
            mat->nlin[k] = 1.0;
            mat->mlin[k] = 3.0;
            mat->rlin[k] = 0.0;
            mat->Qlin[k] = 375e3;
            mat->Vlin[k] = 4.0e-6;            // 4e-6 Annelore 7e-6 Lorenzo
            if ( model->force_act_vol_ast == 1 ) mat->Vlin[k] = model->act_vol_dif_ast;
            mat->Alin[k] = 1.5000e-15;
            mat->flin[k] = 0.0;
            mat->alin[k] = 0.0;
            success      = 1;
            break;
            
        case 41 :
            printf("Olivine Wet Diffusion creep - Hirth & Kohlstedt (2003):\n");
            mat->tlin[k] = 1;
            mat->nlin[k] = 1.0;
            mat->mlin[k] = 3.0;
            mat->rlin[k] = 1.0;
            mat->Qlin[k] = 375e3;
            mat->Vlin[k] = 7e-6;
            mat->Alin[k] = 2.5000e-23;
            mat->flin[k] = 1e9;
            mat->alin[k] = 0.0;
            success      = 1;
            break;
    }
    
    // Scaling
    mat->tlin[k]   /= 1.0;
    mat->nlin[k]   /= 1.0;
    mat->mlin[k]   /= 1.0;
    mat->rlin[k]   /= 1.0;
    mat->Qlin[k]   /= scaling->J;
    mat->Vlin[k]   /= pow(scaling->L,3.0);
    mat->Alin[k]   /= (pow(scaling->S, -mat->nlin[k]) * pow(scaling->S, -mat->rlin[k]) * pow(scaling->L, mat->mlin[k]) / scaling->t);
    mat->flin[k]   /= scaling->S;
    mat->alin[k]   /= 1.0;

    // Calculate correction factor for invariant formulation
    if (mat->tlin[k]==0) mat->Flin[k]   = 1.0;
    if (mat->tlin[k]==1) mat->Flin[k]   = 1.0/6.0*pow(2.0,1.0/mat->nlin[k]) * pow(3.0,(mat->nlin[k]-1)/2.0/mat->nlin[k]);
    if (mat->tlin[k]==2) mat->Flin[k]   = 1.0/4.0*pow(2.0,1.0/mat->nlin[k]);
    
    // Cancel correction factor
    if (number<0) mat->Flin[k]   = 1.0;
    
    if ( success==0 ) { printf("Error: Non existing Linear flow law number\n"); exit(12);}
    
    // Print to screen
    printf("t = %1.0lf  n = %1.1lf  m = %1.1lf  r = %1.1lf  Q = %2.2e J  V = %2.2e m^3  A = %2.2e Pa^-n/s  f = %2.2e Pa  a = %1.1lf  F = %1.2lf\n", mat->tlin[k], mat->nlin[k], mat->mlin[k], mat->rlin[k], mat->Qlin[k]*scaling->J, mat->Vlin[k]*pow(scaling->L,3), mat->Alin[k]*(pow(scaling->S, -mat->nlin[k]) * pow(scaling->L, mat->mlin[k]) / scaling->t), mat->flin[k]*scaling->S, mat->alin[k], mat->Flin[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataGBS( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            printf ("should not be here\n");
            success      = 0;
            break;
            
        case 40 :
            printf("Olivine Low T GBS T<1250°C - Hirth & Kohlstedt (2003):\n");
            mat->tgbs[k] = 1;
            mat->ngbs[k] = 3.5;
            mat->mgbs[k] = 2.0;
            mat->rgbs[k] = 0.0;
            mat->Qgbs[k] = 400e3;
            mat->Vgbs[k] = 0.0;
            mat->Agbs[k] = 6.5000e-30;
            mat->fgbs[k] = 0.0;
            mat->agbs[k] = 0.0;
            success      = 1;
            break;
            
        case 41 :
            printf("Olivine Low T GBS T>1250°C - Hirth & Kohlstedt (2003):\n");
            mat->tgbs[k] = 1;
            mat->ngbs[k] = 3.5;
            mat->mgbs[k] = 2.0;
            mat->rgbs[k] = 0.0;
            mat->Qgbs[k] = 600e3;
            mat->Vgbs[k] = 0.0;
            mat->Agbs[k] = 4.7000e-23;
            mat->fgbs[k] = 0.0;
            mat->agbs[k] = 0.0;
            success      = 1;
            break;
    }
    
    // Scaling
    mat->tgbs[k]   /= 1;
    mat->ngbs[k]   /= 1.0;
    mat->mgbs[k]   /= 1.0;
    mat->rgbs[k]   /= 1.0;
    mat->Qgbs[k]   /= scaling->J;
    mat->Vgbs[k]   /= pow(scaling->L,3.0);
    mat->Agbs[k]   /= (pow(scaling->S, -mat->ngbs[k]) * pow(scaling->S, -mat->rgbs[k]) * pow(scaling->L, mat->mgbs[k]) / scaling->t);
    mat->fgbs[k]   /= scaling->S;
    mat->agbs[k]   /= 1.0;
    
    // Calculate correction factor for invariant formulation
    if (mat->tgbs[k]==0) mat->Fgbs[k]   = 1.0;
    if (mat->tgbs[k]==1) mat->Fgbs[k]   = 1.0/6.0*pow(2.0,1.0/mat->ngbs[k]) * pow(3.0,(mat->ngbs[k]-1.0)/2.0/mat->ngbs[k]);
    if (mat->tgbs[k]==2) mat->Fgbs[k]   = 1.0/4.0*pow(2.0,1.0/mat->ngbs[k]);
    
    // Cancel correction factor
    if (number<0) mat->Fgbs[k]   = 1.0;
    
    if ( success==0 ) { printf("Error: Non existing GBS flow law number\n"); exit(12);}
    
    // Print to screen
    printf("t = %1.0lf  n = %1.1lf  m = %1.1lf  r = %1.1lf  Q = %2.2e J  V = %2.2e m^3  A = %2.2e Pa^-n/s  f = %2.2e Pa  a = %1.1lf  F = %1.2f\n", mat->tgbs[k], mat->ngbs[k], mat->mgbs[k], mat->rgbs[k], mat->Qgbs[k]*scaling->J, mat->Vgbs[k]*pow(scaling->L,3.0), mat->Agbs[k]*(pow(scaling->S, -mat->ngbs[k]) * pow(scaling->L, mat->mgbs[k]) / scaling->t), mat->fgbs[k]*scaling->S, mat->agbs[k], mat->Fgbs[k] );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void ReadDataExponential( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            printf ("should not be here\n");
            success      = 0;
            break;
            
        case 25 :
            printf("Plagioclase Peierls creep - Azuma et al., 2014 / Regularized - Kameyama et al. (1999):\n");
            mat->texp[k] = 1;
            mat->qexp[k] = 2.0;
            mat->Gexp[k] = 0.2;
            mat->Sexp[k] = 9.831e9;
            mat->Qexp[k] = 431.0e3;
            mat->Vexp[k] = 0e-6;
            mat->Eexp[k] = 3.02e-9;
            mat->nexp[k] = 2.0;         // sigma^2
            success      = 1;
            break;

            
        case 40 :
            printf("Olivine Peierls creep - Evans & Goetze (1979) / Regularized - Kameyama et al. (1999):\n");
            mat->texp[k] = 1;
            mat->qexp[k] = 2.0;
            mat->Gexp[k] = 0.1;
            mat->Sexp[k] = 8.5e9;
            mat->Qexp[k] = 5.4e5;
            mat->Vexp[k] = 0.0e-6;
            mat->Eexp[k] = 5.7e11;
            mat->nexp[k] = 0.0;
            success      = 1;
            break;
    }
    
    // Scaling
    mat->texp[k]   /=  1.0;
    mat->qexp[k]   /=  1.0;
    mat->Gexp[k]   /=  1.0;
    mat->Sexp[k]   /= scaling->S;
    mat->Qexp[k]   /= scaling->J;
    mat->Vexp[k]   /= pow(scaling->L, 3.0);
    mat->Eexp[k]   /= scaling->E*pow(scaling->S,-mat->nexp[k]);
    
    // Cancel correction factor
    if (number<0) mat->texp[k]   = 0;
    
    if ( success==0 ) { printf("Error: Non existing Exponential flow law number\n"); exit(12);}
    
    // Print to screen
    printf("t = %1.0lf  q = %1.1lf  G = %1.1lf  S = %2.2e Pa  Q = %2.2e J  V = %2.2e m^3  E = %2.2e 1/s\n", mat->texp[k], mat->qexp[k], mat->Gexp[k], mat->Sexp[k]*scaling->S, mat->Qexp[k]*scaling->J, mat->Vexp[k]*pow(scaling->L,3.0), mat->Eexp[k]*scaling->E);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataGSE( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            printf ("should not be here\n");
            success      = 0;
            break;
            
        case 10 :
            printf("Calcite paleowattmeter - Austin & Evans (2002); Covey-Crump (1997):\n");
            mat->ppzm[k] = 3.0;
            mat->Kpzm[k] = 2.5e9* pow(10.0,-6.0*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 175.0e3                            / scaling->J;
            mat->Gpzm[k] = 1                                  / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = M_PI;
            mat->Lpzm[k] = 0.1;
            success      = 1;
            break;
            
        case 11 :
            printf("Calcite paleopiezometer - Rutter (1995):\n");
            mat->ppzm[k] = 1.14;
            mat->Kpzm[k] = pow(10.0,3.31)*1.0e-6*pow(1.0e6,mat->ppzm[k]);
            mat->Kpzm[k]/= (scaling->L*pow(scaling->S,mat->ppzm[k]));
            mat->Qpzm[k] = 0.0;
            mat->Gpzm[k] = 0.0;
            mat->cpzm[k] = 0.0;
            mat->Lpzm[k] = 0.0;
            success      = 1;
            break;
            
        case 12 :
            printf("Calcite paleopiezometer - Schmid (1977):\n");
            mat->ppzm[k] = 0.9767;
            mat->Kpzm[k] = pow(10.0,2.63)*1.0e-6*pow(1.0e6,mat->ppzm[k]);
            mat->Kpzm[k]/= (scaling->L*pow(scaling->S,mat->ppzm[k]));
            mat->Qpzm[k] = 0.0;
            mat->Gpzm[k] = 0.0;
            mat->cpzm[k] = 0.0;
            mat->Lpzm[k] = 0.0;
            success      = 1;
            break;
            
        case 40 :
            printf("Olivine  - Thielmann et al., (2015):\n");
            mat->ppzm[k] = 2.0;
            mat->Kpzm[k] = 607.0 * pow(10.0, -6.0*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 200.0e3                              / scaling->J;
            mat->Gpzm[k] = 0.1                                  / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = 1.2197; // = 1/Fr from Rozel's papers
            mat->Lpzm[k] = 0.1;
            success      = 1;
            break;
            
        case 41 :
            printf("Olivine - pyroxene mixture - Thielmann et al., (2015):\n");
            mat->ppzm[k] = 4.0;
            mat->Kpzm[k] = 497075.0 * pow(10, -6*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 300.0e3                             / scaling->J;
            mat->Gpzm[k] = 0.1                                 / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = 1.2197; // = 1/Fr from Rozel's papers
            mat->Lpzm[k] = 0.1;
            success      = 1;
            break;
            
        case 42 :
            printf("Olivine - Speciale et al., (2020):\n");
            mat->ppzm[k] = 3.2;
            mat->Kpzm[k] = 1800.0 * pow(10, -6*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 620.0e3                           / scaling->J;
            mat->Vpzm[k] = 5e-6                              / pow(scaling->L,3.0);
            mat->Gpzm[k] = 1.4                               / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = 3.1415; // = 1/Fr from Rozel's papers
            mat->Lpzm[k] = 0.01;
            success      = 1;
            break;
    }
    
    // Scaling done above

    if ( success==0 ) { printf("Error: Non existing grain size evolution law number\n"); exit(12);}
    
    // Print to screen
    printf("p = %1.0lf  K = %2.2e  Q = %2.2e   G = %2.2e J  c = %2.2e m^3  L = %2.2e 1/s\n", mat->ppzm[k], mat->Kpzm[k]*(pow(scaling->L, mat->ppzm[k]) / scaling->t), mat->Qpzm[k]*scaling->J, mat->Gpzm[k]*(scaling->J * pow(scaling->L, -2.0) ), mat->cpzm[k], mat->Lpzm[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
