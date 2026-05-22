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
#include "mdoodz-private.h"
#include "ctype.h"

#include "mdoodz-log.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#endif

#ifdef _VG_
#define LOG_INFO(...) LOG_INFO("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataPowerLaw( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            LOG_INFO("should not be here");
            success      = 0;
            break;
            
        /****************** SPECIAL CASE: user-defined power law flow law ******************/
        case 1 :
            LOG_INFO("'Homemade' power law flow of the form: eta = eta0 * exp(-Q/n/R/T) * Eii^(1/n - 1):");
            mat->tpwl[k] = 0.0;
            mat->npwl[k] = mat->npwl[k];
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = mat->Qpwl[k];
            mat->Vpwl[k] = 0.0;
            mat->Apwl[k] = pow(mat->eta0[k]*scaling->eta, -mat->npwl[k]); // rescale eta0 here, as A will be scaled later on
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

        case 2 :
            LOG_INFO("'Homemade' power law flow with constrained reference stress (tau0)");
            mat->tpwl[k] = 0.0;
            mat->npwl[k] = mat->npwl[k];
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = mat->Qpwl[k];
            mat->Vpwl[k] = 0.0;
            mat->Apwl[k] = pow( pow(2.0, mat->npwl[k]) * mat->eps0[k]/pow(mat->tau0[k], mat->npwl[k]), -mat->npwl[k]);//pow(mat->eta0[k], -mat->npwl[k]);
            LOG_INFO("pre exp = %2.2e", mat->Apwl[k]);
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;

        case 3 :
            LOG_INFO("'Homemade' power law flow of the form: eta = A * exp(-Q/n/R/T) * Eii^(1/n - 1):");
            mat->tpwl[k] = 0.0;
            mat->npwl[k] = mat->npwl[k];
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = mat->Qpwl[k];
            mat->Vpwl[k] = 0.0;
            mat->Apwl[k] = mat->Apwl[k]; // rescale eta0 here, as A will be scaled later on
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = mat->apwl[k];
            success      = 1;
            break;

        /******************************** Crust flow laws *********************************/
                
        case 10:
            LOG_INFO("Westerly Granite (dry) - Hansen & Carter (1983):");
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
            LOG_INFO("Maryland Diabase - Mackwell et al. (1998):");
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
            LOG_INFO("Maryland Diabase - Carter & Tsenn (1987):");
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
            LOG_INFO("Wet Quartzite - Ranalli (1995):");
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
            LOG_INFO("Rocksalt - Ranalli (1995):");
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
            LOG_INFO("Calcite - Renner et al. (2002):");
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
            LOG_INFO("Felsic granulite - Ranalli (1995):");
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
            LOG_INFO("Mafic granulite - Ranalli (1995):");
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
            LOG_INFO("Calcite - Schmid et al. (1977):");
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
            LOG_INFO("Anorthite 100 - 0.007 H20 Wt%% - Rybacki & Dresen (2000, 2004):");
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
            LOG_INFO("Garnet - Ji and Martignole (1994):");
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
            LOG_INFO("Omphacite - Zhang et al. (2006):");
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
            LOG_INFO("Qtz - Koch et al. (1989):");
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
            LOG_INFO("Fsp - Shelton and Tullis (1981):");
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
            LOG_INFO("Mica - Kronenberg et al. (1990):");
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
            LOG_INFO("Plagioclase - Shelton (1981) - Ph.D. Thesis:");
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
            LOG_INFO("Wet Quartzite - Burov (2011), Brace and Kohlstedt (1980), Kirby and Kronenberg (1987), Kohlstedt et al. (1995):");
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
            LOG_INFO("Dry Quartz - Ranalli et al. (1995, 1997), Ranalli (2003):");
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
            LOG_INFO("Qtz - Hirth et al. (2001) - uniaxial:");
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
            LOG_INFO("Basalt - Hacker & Chritie (1990):");
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
            LOG_INFO("Anorthite dry - Rybacki and Dresen (2000):");
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
            LOG_INFO("DRY CPX - Bystricky & Mackwell, 2001:");
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
            LOG_INFO("Anorthite 60 - Rybacki & Dresen (2004):"); // added by Pauline
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
        
        case 33:
            LOG_INFO("Garnet - Ji and Martignole (1994):"); // Activation volume set to 10 for incompressible (like in Mei 2010)
            mat->tpwl[k] = 1;
            mat->npwl[k] = 2.22;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 485.0e3;
            mat->Vpwl[k] = 10e-6;
            mat->Apwl[k] = 2.7952e-07;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 0.0;
            success      = 1;
            break;
        
	    case 34:
            LOG_INFO("Westerly Granite + melt weakening (dry) - Hansen & Carter (1983):");
            mat->tpwl[k] = 1;
            mat->npwl[k] = 3.3;
            mat->mpwl[k] = 0.0;
            mat->rpwl[k] = 0.0;
            mat->Qpwl[k] = 186.5e3;
            mat->Vpwl[k] = 0.0e-6;
            mat->Apwl[k] = 3.1623e-26;
            mat->fpwl[k] = 0.0;
            mat->apwl[k] = 50.0;
            success      = 1;
            break;
            
        /******************************** Mantle flow laws ********************************/
            
        case 40:
            LOG_INFO("Olivine Dry Dislocation creep - Hirth & Kohlstedt (2003):");
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
            LOG_INFO("Olivine Wet Dislocation creep - Hirth & Kohlstedt (2003):");
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
            LOG_INFO("Dry Olivine average from Ranalli (1995) - Kirby & Kronenberg (1987):");
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
            LOG_INFO("Wet Olivine average from Ranalli (1995) - Kirby & Kronenberg (1987):");
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
            LOG_INFO("Olivine - Chopra and Paterson 1984 :");
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
            LOG_INFO("Olivine Dry Dislocation creep - Hirth & Kohlstedt (2003):"); // added by Pauline
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
            LOG_INFO("Serpentine 1 GPa - Hilairet et al. (2007):"); // added by Pauline
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
    
    // Force melt weakening factor for all phases
    if ( model->force_melt_weak == 1 ) mat->apwl[k] = model->melt_weak;
    
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
    
    if ( success==0 ) { LOG_ERR("Error: Non existing Power Law flow law number"); exit(12);}
    
    // Print to screen 
    LOG_INFO("t = %1.0lf  n = %1.1lf  m = %1.1lf  r = %1.1lf  Q = %2.2e J  V = %2.2e m^3  A = %2.2e Pa^-n/s  f = %2.2e Pa  a = %1.1lf  F = %1.2lf", mat->tpwl[k], mat->npwl[k], mat->mpwl[k], mat->rpwl[k], mat->Qpwl[k]*scaling->J, mat->Vpwl[k]*pow(scaling->L,3.0), mat->Apwl[k]*(pow(scaling->S, -mat->npwl[k]) * pow(scaling->L, mat->mpwl[k]) / scaling->t), mat->fpwl[k]*scaling->S, mat->apwl[k], mat->Fpwl[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataLinear( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            LOG_INFO("should not be here");
            success      = 0;
            break;

        /****************** SPECIAL CASE: user-defined power law flow law ******************/
        case 1 :
            LOG_INFO("Please, set the parameters for this user-defined flow law in MDLIB/FlowLaws.c");
            success      = 0;
            break;
        /******************************** Crust flow laws *********************************/
               
        case 15 :
            LOG_INFO("Calcite Diffusion creep - Herwegh (2003):");
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
            LOG_INFO("Anorthite 100 - 0.007 H20 Wt%% - Rybacki & Dresen (2000, 2004):");
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
            LOG_INFO("Plagioclase Diffusion creep - Rybacki and Dresen (2000):");
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
            LOG_INFO("Olivine Dry Diffusion creep - Hirth & Kohlstedt (2003):");
            mat->tlin[k] = 1.0;
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
            LOG_INFO("Olivine Wet Diffusion creep - Hirth & Kohlstedt (2003):");
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
    
    if ( success==0 ) { LOG_ERR("Error: Non existing Linear flow law number"); exit(12);}
    
    // Print to screen
    LOG_INFO("t = %1.0lf  n = %1.1lf  m = %1.1lf  r = %1.1lf  Q = %2.2e J  V = %2.2e m^3  A = %2.2e Pa^-n/s  f = %2.2e Pa  a = %1.1lf  F = %1.2lf", mat->tlin[k], mat->nlin[k], mat->mlin[k], mat->rlin[k], mat->Qlin[k]*scaling->J, mat->Vlin[k]*pow(scaling->L,3), mat->Alin[k]*(pow(scaling->S, -mat->nlin[k]) * pow(scaling->L, mat->mlin[k]) / scaling->t), mat->flin[k]*scaling->S, mat->alin[k], mat->Flin[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataGBS( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            LOG_INFO("should not be here");
            success      = 0;
            break;
            
        case 40 :
            LOG_INFO("Olivine Low T GBS T<1250°C - Hirth & Kohlstedt (2003):");
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
            LOG_INFO("Olivine Low T GBS T>1250°C - Hirth & Kohlstedt (2003):");
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
    
    if ( success==0 ) { LOG_ERR("Error: Non existing GBS flow law number"); exit(12);}
    
    // Print to screen
    LOG_INFO("t = %1.0lf  n = %1.1lf  m = %1.1lf  r = %1.1lf  Q = %2.2e J  V = %2.2e m^3  A = %2.2e Pa^-n/s  f = %2.2e Pa  a = %1.1lf  F = %1.2f", mat->tgbs[k], mat->ngbs[k], mat->mgbs[k], mat->rgbs[k], mat->Qgbs[k]*scaling->J, mat->Vgbs[k]*pow(scaling->L,3.0), mat->Agbs[k]*(pow(scaling->S, -mat->ngbs[k]) * pow(scaling->L, mat->mgbs[k]) / scaling->t), mat->fgbs[k]*scaling->S, mat->agbs[k], mat->Fgbs[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void ReadDataExponential( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            LOG_INFO("should not be here");
            success      = 0;
            break;
            
        case 25 :
            LOG_INFO("Plagioclase Peierls creep - Azuma et al., 2014 / Regularized - Kameyama et al. (1999):");
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
            LOG_INFO("Olivine Peierls creep - Evans & Goetze (1979) / Regularized - Kameyama et al. (1999):");
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
    
    if ( success==0 ) { LOG_ERR("Error: Non existing Exponential flow law number"); exit(12);}
    
    // Print to screen
    LOG_INFO("t = %1.0lf  q = %1.1lf  G = %1.1lf  S = %2.2e Pa  Q = %2.2e J  V = %2.2e m^3  E = %2.2e 1/s", mat->texp[k], mat->qexp[k], mat->Gexp[k], mat->Sexp[k]*scaling->S, mat->Qexp[k]*scaling->J, mat->Vexp[k]*pow(scaling->L,3.0), mat->Eexp[k]*scaling->E);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataGSE( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success=0;
    
    switch ( abs(number) ) {
            
        case 0 :
            LOG_INFO("should not be here");
            success      = 0;
            break;
            
        case 9 :
            LOG_INFO("MODIFIED - Low lambda: Calcite paleowattmeter - Austin & Evans (2002); Covey-Crump (1997):");
            mat->ppzm[k] = 3.0;
            mat->Kpzm[k] = 2.5e9* pow(10.0,-6.0*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 175.0e3                            / scaling->J;
            mat->Gpzm[k] = 1.0                                / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = M_PI;
            mat->Lpzm[k] = 0.01;
            success      = 1;
            break;

        case 10 :
            LOG_INFO("Calcite paleowattmeter - Austin & Evans (2002); Covey-Crump (1997):");
            mat->ppzm[k] = 3.0;
            mat->Kpzm[k] = 2.5e9* pow(10.0,-6.0*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 175.0e3                            / scaling->J;
            mat->Gpzm[k] = 1.0                                / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = M_PI;
            mat->Lpzm[k] = 0.1;
            success      = 1;
            break;
            
        case 11 :
            LOG_INFO("Calcite paleopiezometer - Rutter (1995):");
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
            LOG_INFO("Calcite paleopiezometer - Schmid (1977):");
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
            LOG_INFO("Olivine  - Thielmann et al., (2015):");
            mat->ppzm[k] = 2.0;
            mat->Kpzm[k] = 607.0 * pow(10.0, -6.0*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 200.0e3                              / scaling->J;
            mat->Gpzm[k] = 0.1                                  / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = 1.2197; // = 1/Fr from Rozel's papers
            mat->Lpzm[k] = 0.1;
            success      = 1;
            break;
            
        case 41 :
            LOG_INFO("Olivine - pyroxene mixture - Thielmann et al., (2015):");
            mat->ppzm[k] = 4.0;
            mat->Kpzm[k] = 497075.0 * pow(10, -6*mat->ppzm[k]) / (pow(scaling->L, mat->ppzm[k]) / scaling->t);
            mat->Qpzm[k] = 300.0e3                             / scaling->J;
            mat->Gpzm[k] = 0.1                                 / (scaling->J * pow(scaling->L, -2.0) );
            mat->cpzm[k] = 1.2197; // = 1/Fr from Rozel's papers
            mat->Lpzm[k] = 0.1;
            success      = 1;
            break;
            
        case 42 :
            LOG_INFO("Olivine - Speciale et al., (2020):");
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

    if ( success==0 ) { LOG_ERR("Error: Non existing grain size evolution law number"); exit(12);}
    
    // Print to screen
    LOG_INFO("p = %1.0lf  K = %2.2e  Q = %2.2e   G = %2.2e J  c = %2.2e m^3  L = %2.2e 1/s", mat->ppzm[k], mat->Kpzm[k]*(pow(scaling->L, mat->ppzm[k]) / scaling->t), mat->Qpzm[k]*scaling->J, mat->Gpzm[k]*(scaling->J * pow(scaling->L, -2.0) ), mat->cpzm[k], mat->Lpzm[k]);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadDataKinetics( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    
    int success = 0;
    
    switch ( abs(number) ) {

        case 1 :
            LOG_INFO("Directly use tau_kin %d %d", model->kinetics, mat->kin[k]);
            success      = 1;
            break;
            
        case 9 :
            LOG_INFO("Quartz-Coesite - Mosenfelder & Bohlen (1997);");
            mat->Skin[k] = 3.35/1e-4 / (1.0/scaling->L);
            mat->kkin[k] = 0.185     / (scaling->L/scaling->t/scaling->T);
            mat->Qkin[k] = 243e3     / scaling->J;
            success      = 1;
            break;
    }

    if ( success==0 ) { LOG_ERR("Error: Non existing kinetic number %d", number); exit(12);}

}

/*--------------------------------------------------------------------------------------------------------------------*/

// Hansen-group olivine viscous anisotropy (case 1). Calibrated against the
// combined 38-point Hansen+14 (EPSL 387, 157) + Hansen+16 Part 1 (EPSL 445,
// 92) dataset by least-squares fit on per-sample (γ, M) data: M_∞ = 0.536,
// γ_e = 3.96. Asymptote δ_∞ = 24.5·M_∞ + 1 = 14.132 agrees with Hansen+16
// Part 1 Eq.3-4 stress-aware direct calculation (1.39/0.73)^4.1 = 14.02
// within 0.8 % (independent two-method validation).
//
// Slope 24.5 ± 5.5 (1σ; Hansen+12 Nature Fig 3b annotation, n = 3). The
// asymptote inherits a ±22 % (1σ) band: δ_∞ = 14.132 ± 2.95, 1σ interval
// [11.18, 17.08]. Case 7 inherits the same calibration-wide slope uncertainty.
//
// Reproducible from misc/aniso_fstrain/calibrate/calibrate_olivine_hansen.py
// (raw data: misc/aniso_fstrain/data/olivine/raw_Hansen_etal_{2014,2016}_EPSL.csv).
//
// Calibrated for FRESH-CPO olivine (initially random, dry Fo90, lab torsion).
// For naturally-deformed peridotite with complex deformation history (mantle
// xenoliths, shear-zone mylonites) where pre-existing CPO and repeated
// deformation damp fabric strength, use case 7 (anisoDelta_HansenOlivine_Damped).
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR    (simple-shear inversion of FS_AR)
//   M(γ)         = 0.536 · (1 − exp(−γ/3.96))   (Skemer M-index saturation)
//   δ            = 24.5 · M + 1                 (Hansen+12 Fig 3b linear fit)
static double anisoDelta_HansenOlivine( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.536;
    const double gamma_e   = 3.96;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return 24.5 * M + 1.0;
}

// Damped olivine viscous anisotropy (case 7). Calibrated against the combined
// 93-sample non-Hansen literature dataset:
//   Tasaka 2016 (JGR-SE 121, 92):              12 samples (wet Fo50 lab torsion)
//   Boneh & Skemer 2014 (EPSL 406, 213):       16 samples (Aheim dunite, pre-CPO)
//   Kumamoto+19 (JGR-SE 124, 12763):           30 samples (natural Josephine, γ ∈ [0, 20])
//   Bernard+19 (G3 20, 3469):                  35 samples (natural xenoliths)
//
// Free LSQ on the 93-point dataset: M_∞ = 0.157, γ_e = 0.764, RMS(M) = 0.069
// → δ_∞ = 4.84. Five subset refits (full corpus, without Bernard, without
// Boneh, Tasaka+Kumamoto only, Kumamoto only) all land at δ_∞ ∈ [4.80, 5.10]
// — robust result. Bernard+19 M values are re-derived from J via the natural-
// sample fit M = (J−1)/24 (R² = 0.825) which supersedes Bernard's as-reported
// M (implicitly biased toward Skemer (J−1)/15 lab-fit); details in
// misc/aniso_fstrain/calibrate/aniso_data.py and the notes companion.
//
// Same δ-vs-M slope (24.5 ± 5.5, 1σ) as case 1 because the linear δ–M
// relation is set by single-crystal olivine anisotropy, not by M_∞ magnitude.
// Case 7: δ_∞ committed = δ_∞ LSQ = 4.84 (1σ band [3.98, 5.70]).
//
// Three independent damping mechanisms:
//   (a) hydrous + Fe-rich olivine (Tasaka 2016, Demouchy+09 GRL, Karato+86):
//       water and Fe shift slip systems and slow CPO development.
//   (b) pre-existing CPO disruption (Boneh & Skemer 2014): unfavorable starting
//       fabric DECREASES M before any new CPO grows.
//   (c) complex natural strain history (Kumamoto+19, Bernard+19): multiple
//       non-coaxial deformation events with re-foliation; J/M does not grow
//       systematically with bulk strain. Case 7 is regime-mismatched for
//       subduction-zone serpentinization (Kumamoto+19's water source is melt,
//       not slab fluid).
//
// Reproducible from misc/aniso_fstrain/calibrate/calibrate_olivine_damped.py
// (raw data: misc/aniso_fstrain/data/olivine/raw_{Tasaka,Boneh_Skemer,Kumamoto,
//  Bernard}_*.csv).
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//   M(γ)         = 0.16 · (1 − exp(−γ/0.76))
//   δ            = 24.5 · M + 1            (same slope as case 1)
static double anisoDelta_HansenOlivine_Damped( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.16;
    const double gamma_e   = 0.76;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return 24.5 * M + 1.0;
}

// Quartz viscous anisotropy — Pennacchioni+10 case 3. Calibrated against
//   Pennacchioni+10 (JGR 115, B12405): 17 natural-shear-zone samples,
//     γ_fol ∈ [0.3, 15], T ~ 500°C, ODF J ∈ [2.2, 12].
//   Blackford+24 (Tectonics 43, e2023TC008166): 38 Northern Snake Range
//     samples, γ_oct ∈ [0.85, 2.35] (γ_simple_eq ∈ [2.5, 17.7] under
//     plane-strain inversion), T ~ 500–650°C, pole-figure J ∈ [1.3, 4.85].
//
// J→M conversion via Skemer+05 M = (J−1)/15. Applied to Pennacchioni's ODF J
// (the metric Skemer's relation was developed for); Blackford's pfJ is on a
// different scale (pfJ < ODF J for the same texture) and is used for context.
// Both datasets agree qualitatively that quartz CPO grows steadily and
// saturates slowly (γ_e ≥ 3.5).
//
// Free LSQ on Pennacchioni: M_∞ = 0.74, γ_e = 3.59, RMS(M) = 0.032; committed
// rounded to (0.75, 3.50). NOTE — the 17 (γ, J) values are audit-flagged as
// likely hand-fitted from Pennacchioni's qualitative WDV/MDV/SDV regime
// descriptions rather than per-sample MTEX measurements (residual smoothness
// σ/M_∞ = 0.045 is 3–4× below the real-data floor seen in Hansen lab and
// Blackford natural-quartz corpora). Case 3 is therefore best treated as a
// PHYSICAL-PRIOR commitment consistent with natural-quartz-mylonite
// expectations; for data-driven natural quartz use case 6 (Blackford).
// See misc/aniso_fstrain/audit/Pennacchioni_etal_2010_audit.md for the audit.
//
// The δ-vs-M slope = 3.0 is fixed by VPSC-style bounds for quartz aggregate
// viscous anisotropy at saturated CPO (δ ≈ 2-3 from Heilbronner-Tullis-style
// results; quartz is a much weaker viscous-anisotropy mineral than olivine
// or even calcite). Asymptote δ_∞ = 3·M_inf + 1 = 3.25.
//
// NOT applicable below ~500°C where basal-a / rhomb-a slip changes the
// fabric kinematics — see misc/aniso_fstrain/notes/quartz_calibration.md
// §6 for the T-range caveat. For T < 500°C scenarios, fall back to
// ani_fstrain = 1 with ani_fac_max ≈ 3 per phase.
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR    (simple-shear inversion of FS_AR)
//   M(γ)         = 0.75 · (1 − exp(−γ/3.50))    (Skemer M-index saturation)
//   δ            = 3.0 · M + 1                  (quartz-specific slope)
static double anisoDelta_Quartz( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.75;
    const double gamma_e   = 3.50;
    const double slope     = 3.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Quartz viscous anisotropy — Blackford+24 calibration (case 6).
// Calibrated against Blackford et al. (2024) Tectonics 43, e2023TC008166 —
// 38 quartzite samples from the Northern Snake Range metamorphic core
// complex, NV with 3D strain ellipsoids on detrital quartz clasts and
// pole-figure J-index of c-axis (density norm pfJ ∈ [1.3, 4.85]).  Strain
// reported as Nadai octahedral shear strain ε; converted to simple-shear-
// equivalent γ via γ_eq = 2·sinh(ε/√2) under the plane-strain
// approximation (samples described as plane-strain to slightly
// constrictional).  γ_eq spans [1.28, 5.08].
//
// pfJ scale differs from Pennacchioni's ODF J (case 3): pfJ < ODF J for
// the same texture, so M = (J−1)/15 gives ~5× smaller M values.  The
// δ-vs-M slope is rescaled to 11.0 (vs 3.0 for case 3) so the saturation
// asymptote δ_∞ = 3.20 matches case 3's δ_∞ = 3.25 to within rounding —
// both modes describe the same physical mineral, just measured on
// different texture metrics.
//
// Free LSQ on Blackford alone gives M_inf ≈ 0.17–0.29 across γ_e = 3–10
// (γ_e is poorly constrained by the data — saturation is not strongly
// expressed in the γ_eq ∈ [1.3, 5.1] range). Committed (M_inf = 0.20,
// γ_e = 4.0) is centered in the LSQ envelope.
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//   M(γ)         = 0.20 · (1 − exp(−γ/4.00))
//   δ            = 11.0 · M + 1
static double anisoDelta_Quartz_Blackford( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.20;
    const double gamma_e   = 4.00;
    const double slope     = 11.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Calcite viscous anisotropy — mid-low-T regime (T ≤ 650°C, case 4).
// Calibrated against Barnhoorn+04 J. Struct. Geol. 26, 885 — Carrara marble
// torsion at 500°C (4 pts) + 600°C (3 pts) — 7 datapoints total. Subgrain-
// rotation recrystallization dominates at this T regime; CPO is continually
// randomized at high γ, capping fabric at modest strength.
//
// Free LSQ fit yields M_inf = 0.22, γ_e = 1.37 (fast saturation by γ ≈ 4).
// RMS(M) = 0.073 — much tighter than the T-averaged fit because the
// physics in this regime is consistent across the 7 datapoints.
// Asymptote δ_∞ = 6 · 0.22 + 1 = 2.33.
//
// Recommended for T < 650°C calcite scenarios. Above ~650°C, grain-boundary
// migration recrystallization takes over and δ continues to strengthen —
// use case 5 (anisoDelta_Calcite_HighT) instead.
static double anisoDelta_Calcite_LowT( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.22;
    const double gamma_e   = 1.37;
    const double slope     = 6.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Calcite viscous anisotropy — high-T regime (T > 650°C, case 5).
// Calibrated against Pieri+01 (5 pts at 727°C) + Barnhoorn+04 (3 pts at
// 727°C, J ≤ 15 filter) — 8 datapoints total.  Grain-boundary-migration
// recrystallization dominates at this T regime, preferentially eliminating
// unfavorable orientations and STRENGTHENING CPO continuously even at
// very high strain.
//
// Free LSQ fit yields M_inf = 0.7334, γ_e = 5.2775.  Asymptote
// δ_∞ = 6 · 0.73 + 1 = 5.40.  RMS(M) = 0.17.
//
// Bruijn+11 audit: 3 (γ, J) points originally attributed to Bruijn+11
// (Tectonophysics 503, 75) cannot be sourced from that paper and were
// retracted from the calibration corpus. Constants above are the LSQ result
// on the 8 verifiable datapoints (Pieri+01 + Barnhoorn+04 727°C with J ≤ 15);
// removal dropped δ_∞ from 5.99 to 5.40 (~10 %) and tightened saturation
// (γ_e from 8.24 to 5.28). See misc/aniso_fstrain/audit/Bruijn_etal_2011_audit.md.
//
// Pieri-vs-Barnhoorn 727°C disagreement at γ=5 (Pieri J=9.76 vs Barnhoorn
// J=12.5) remains a real source of ~17% scatter in the residual.  For
// very-high-γ scenarios (γ > 10) this extrapolates beyond the kept fit
// data; for typical γ < 10 scenarios the curve closely tracks Pieri+01.
static double anisoDelta_Calcite_HighT( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.73;
    const double gamma_e   = 5.28;
    const double slope     = 6.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Calcite viscous anisotropy — calibrated against the combined Pieri+01
// (Tectonophysics 330, 119 — Carrara marble torsion at 1000 K, 300 MPa,
// γ ∈ {0,1,2,5,11}) + Barnhoorn+04 (J. Struct. Geol. 26, 885 — Carrara
// marble torsion at 500/600/727°C, γ ∈ [0.5, 46]) J-index torsion dataset,
// with J→M-index conversion via Skemer+05 (linear range J ≤ 15).  High-J
// Barnhoorn samples (γ > 15 at 727°C, where the linear J→M conversion
// breaks down) are filtered out.
//
// Free LSQ fit on the combined 15-point dataset (5 Pieri + 10 filtered
// Barnhoorn) yields M_inf = 0.4098, γ_e = 3.1504.  Asymptote
// δ_∞ = slope · M_inf + 1 = 3.46.
//
// Bruijn+11 audit: 3 (γ, J) points originally attributed to Bruijn+11 were
// retracted (cannot be sourced from the paper). Constants above reflect the
// LSQ refit on the 15 verifiable datapoints (Pieri+01 + Barnhoorn+04, J ≤ 15);
// M_∞ and δ_∞ are unchanged, γ_e tightens 3.71 → 3.15. See
// misc/aniso_fstrain/audit/Bruijn_etal_2011_audit.md.
//
// The δ-vs-M slope (slope = 6.0) is calcite-specific, distinct from
// olivine's Hansen+12 Fig 3b coefficient of 24.5. It reflects calcite's
// lower viscous-anisotropy-per-fabric-strength ratio (Barnhoorn+04: CPO
// contributes only ~1/3 of calcite weakening; the rest comes from
// grain-size reduction, already captured by MDOODZ's wattmeter). The
// slope value is bracketed by Tommasi+09-style VPSC bounds for calcite
// saturated CPO (δ_saturated ≈ 2-5). NOT a clean Hansen+12-equivalent
// direct lab fit — see misc/aniso_fstrain/notes/calcite_calibration.md §3.
//
// NOT applicable to high-pressure calcite (Schuster+18 HPT regime, > 1.6
// GPa) where the calcite → CaCO₃-II phase transition makes it a different
// material.
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR    (simple-shear inversion of FS_AR)
//   M(γ)         = 0.41 · (1 − exp(−γ/3.15))    (Skemer M-index saturation)
//   δ            = 6.0 · M + 1                  (calcite-specific slope)
static double anisoDelta_Calcite( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.41;
    const double gamma_e   = 3.15;
    const double slope     = 6.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// =====================================================================
// BRACKETED polyphase modes — case 8 and case 9.
// =====================================================================
// IMPORTANT — these two cases differ in PROVENANCE from cases 1–7:
//
//   Cases 1–7: direct least-squares fit of M(γ) = M_inf·(1−exp(−γ/γ_e))
//              to a coherent multi-point γ-vs-J(or M) lab dataset
//              (Hansen+12/+14/+16, Pennacchioni+10, Blackford+24,
//              Pieri+01/Barnhoorn+04 (Bruijn+11 retracted; see
//              misc/aniso_fstrain/audit/Bruijn_etal_2011_audit.md),
//              Tasaka+Boneh+Kumamoto+Bernard).  RMS reported per case
//              in the calibration notes.  These are the rigorous direct-fit modes.
//
//   Cases 8–9: BRACKETED point estimates from disparate sources.
//              No co-located (γ, J_bulk) dataset exists for polyphase
//              quartz-mica or plagioclase-bearing rocks.  Constants
//              are bracketed between literature bounds (saturation γ
//              from one source, M_inf at saturation from another,
//              slope from theory) and committed as a defensible point
//              estimate — analogous to Ranalli (1995)'s felsic
//              granulite (pwlv=16) and mafic granulite (pwlv=17) flow
//              laws, which are also bracketed compilations rather
//              than direct lab fits.  Quality tier matches pwlv 16/17
//              in the power-law database.
//
// Users seeking direct-fit rigor should use cases 1–7.  Users needing
// a polyphase mode for granitoid-mylonite or gabbro-mylonite scenarios
// should use cases 8 or 9 with the understanding that the constants
// are bracketed estimates, not regression outputs.  See
// misc/aniso_fstrain/notes/quartz_calibration.md §11 (mica)
// and §12 (plagioclase) for the literature bracket details.
// =====================================================================

// Quartz-mica polyphase BRACKETED — case 8.  Granitoid mylonite,
// quartzite mylonite with muscovite or biotite films, two-mica granite.
// Bracketed from:
//   Tokle+23 (JSG 169, 104835): γ_e ≈ 2–4 from quartz CPO transition
//                                between low-strain (γ ≈ 0.6) and
//                                high-strain (γ ≈ 4) samples
//   Dempsey+11 (Geol. Soc. Lond. Spec. Publ. 360, 33): mica J at
//                                                     saturated
//                                                     mylonite stage
//                                                     mean ≈ 9 →
//                                                     M_mica_sat ≈ 0.53
//   Holyoke & Tullis (2006): only ~10% interconnected mica required
//                            to make mica the dominant rheological
//                            phase
//   Bos & Spiers (2002): foliated quartz-mica fault rock 2-5x weaker
//                        than Byerlee due to easy basal slip
//
// Committed (M_inf=0.40, γ_e=2.5, slope=15) is a volume-weighted
// blend appropriate for typical 10–25% mica content.  Asymptote
// δ_∞ = 15·0.40 + 1 = 7.0.
//
// NOT a direct LSQ fit.  Slope is genuinely under-determined (3 ≤
// slope ≤ 50 across literature); 15.0 is a defensible middle value.
// For mica fractions outside 10–25%, scale δ_∞ proportionally
// (5–7 at 5–10% mica, 7–10 at 15–25% mica, 10–15 at 25–40% mica).
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//   M(γ)         = 0.40 · (1 − exp(−γ/2.50))
//   δ            = 15.0 · M + 1
static double anisoDelta_QuartzMicaBracketed( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.40;
    const double gamma_e   = 2.50;
    const double slope     = 15.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Plagioclase-bearing polyphase BRACKETED — case 9.  Gabbro mylonite,
// anorthosite mylonite, mafic granulite, basalt, diabase.  Modest
// asymptote because plagioclase has intrinsically low single-crystal
// viscous anisotropy AND is further damped in nature by DisGBS at
// high strain.
//
// Bracketed from:
//   Mehl & Hirth (2008) JGR-SE 113, B05202: ODP Hole 735B SWIR oceanic
//                                            gabbro mylonites; monophase
//                                            plagioclase layers M ∈
//                                            [0.05, 0.19] (mean ≈ 0.12)
//   Marti et al. (2018) Solid Earth 9, 985: lab plag-cpx mixtures, weak
//                                            distinct CPO via DPC + GBS
//                                            up to γ_a ≈ 6
//   Gómez Barreiro, Wenk & Vogel (2015) JSG 71, 100: Morin Shear Zone
//                                                    anorthosite mylonite
//                                                    plag pole-figure
//                                                    max only 2.9 m.r.d.,
//                                                    bulk AVp 2.4%
//   Stünitz et al. (2003) Tectonophys 372, 215: dominant slip systems
//                                                (010)[001] and {001}<110>
//                                                (verbatim "[001](010) and
//                                                <110>{001}" per abstract;
//                                                prior "(010)[100]" was a
//                                                transcription error — that
//                                                string is the olivine A-type
//                                                slip system, not plagioclase)
//
// Committed (M_inf=0.12, γ_e=2.0, slope=8.0) reflects monophase
// plagioclase mylonite at saturation.  Asymptote δ_∞ = 8·0.12 + 1 =
// 1.96 ≈ 2.0 — much weaker than olivine (14), comparable to weak
// quartz CPO regime.
//
// NOT a direct LSQ fit, AND the saturating-exponential form is itself
// a simplification — at the highest natural strains, plagioclase CPO
// can WEAKEN as DisGBS / diffusion creep activates.  The monotonic
// curve here over-predicts δ at very high γ (>10) compared to natural
// observations (Mehl & Hirth report polyphase plag layers with M ≈
// 0.04, much lower than the M_inf = 0.12 used here).  For plag-cpx
// polyphase lower-crust scenarios beyond γ ≈ 5, this over-predicts
// anisotropy by ~2x.
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//   M(γ)         = 0.12 · (1 − exp(−γ/2.00))
//   δ            = 8.0 · M + 1
static double anisoDelta_PlagioclaseBracketed( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.12;
    const double gamma_e   = 2.00;
    const double slope     = 8.0;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Calcite-superplastic BRACKETED — case 13.  Fine-grained calcite mylonite
// in low-stress / high-T regime where grain-boundary sliding (GBS) dominates
// and crystallographic preferred orientation is actively destroyed.
// Solnhofen limestone analog (Schmid 1976/77 three-regime calcite rheology,
// regime 3).  Selected by the user when modeling regime-3 conditions
// (σ < σ_c(T) per Schmid 1976 Arrhenius, fine-grained calcite <10 µm,
// T > 600°C); MDLIB does not check regime-3-applicability at runtime —
// physics-correctness is the user's responsibility (matches case-1/2/3/4
// pattern).  For dislocation-creep regime calcite use case 2/4/5 instead.
//
// Bracketed from:
//   Schmid (1976) Tectonophysics 31, T21-T28: σ_c(T) Arrhenius regression
//                                              (4-T-bin LSQ R²=0.86,
//                                               E_a=23.4 kJ/mol,
//                                               σ_c0=34.6 bar);
//                                               regime-2-to-3 transition
//                                               σ_c = 1050 bar at 600°C
//                                               down to 300 bar at 900°C
//   Schmid et al. (1977) Tectonophysics 43, 257-291: Solnhofen limestone
//                                                    superplastic flow
//                                                    paleopiezometer;
//                                                    direct CPO-weakening
//                                                    statement "CPO in
//                                                    superplastic regime
//                                                    is weaker and of
//                                                    different geometry
//                                                    when compared with
//                                                    that in flow regimes
//                                                    1 and 2"; n ∈ (1, 3)
//                                                    + GBS-dominant +
//                                                    Davies+1970
//                                                    "destruction of
//                                                    preferred orientation"
//   De Bresser (2002) JGR 107, 2337: Carrara marble + single-crystal
//                                     compression study; INDEPENDENT
//                                     three-regime corroboration via
//                                     Schmid+1980 citation + Renner&Evans
//                                     2002 confirmation of "failure of
//                                     conventional power law"; n=4 at
//                                     σ~15 MPa, T=1000°C brackets the
//                                     regime-3 onset stress.
//
// Committed (M_inf=0.08, γ_e=1.0, slope=1.5) reflects fine-grained
// calcite mylonite at regime-3 GBS saturation.  Asymptote δ_∞ = 1.5·0.08 + 1
// = 1.12 — very weak anisotropy, much weaker than case-1 olivine (δ_∞ ≈ 2.5)
// and case-2 calcite dislocation-creep (δ_∞ = 3.46).  Consistent with
// Schmid 1976 regime-3 σ ≈ ⅓ of regime-1 σ at same ε̇ (stress-drop ratio
// → weak-anisotropy character) and Schmid+77 fast-saturation in GBS regime.
//
// NOT a direct LSQ fit on (γ, M). σ_c(T) is LSQ-fitted (R² = 0.86, n = 4
// Arrhenius); M_∞ / γ_e / slope are literature-bracketed from Schmid 1976,
// Schmid+77 and De Bresser 2002 (3-paper substrate; same discipline as
// cases 8 / 9). See misc/aniso_fstrain/notes/case13_calibration.md.
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//   M(γ)         = 0.08 · (1 − exp(−γ/1.00))
//   δ            = 1.5 · M + 1
static double anisoDelta_Calcite_Superplastic( double FS_AR ) {
    const double s         = sqrt(FS_AR);
    const double gamma_eff = s - 1.0 / s;
    const double M_inf     = 0.08;
    const double gamma_e   = 1.00;
    const double slope     = 1.5;
    const double M         = M_inf * (1.0 - exp(-gamma_eff / gamma_e));
    return slope * M + 1.0;
}

// Quartz-mica polyphase BRACKETED — case 15.  Granitoid mylonite,
// quartzite mylonite with muscovite films, viewed from the QUARTZ
// phase perspective.  COMPLEMENTARY companion to case-8 (which views
// the same rock from the MICA phase perspective).
//
// case-8  = mica CPO STRENGTHENING with γ — saturating-exp form
//             M(γ) = M_∞ · (1 - exp(-γ/γ_e)),  M_∞ = 0.40, γ_e = 2.50.
// case-15 = quartz CPO WEAKENING with γ — DECAY form
//             M(γ) = M_init · exp(-γ/γ_decay), M_init = 0.359, γ_decay = 15.35.
//
// Bracketed from:
//   Tokle, Hirth & Stunitz (2023) JSG 169, 104835:
//       9-datapoint pooled LSQ on polyphase X_ms ∈ {5, 10, 25} % × γ ∈
//       {0, 0.6, 4}. Reproducible from
//       misc/aniso_fstrain/data/quartz_mica/case15_lsq_pooled.{py,txt,png}.
//       Pole-figure max in m.u.d.: M_init = 6.02 ± 0.14 (2 %),
//       γ_decay = 15.35 ± 2.86 (19 %), R² = 0.82, RMSE = 0.28 m.u.d.
//   Tokle+23 §3.2.3 "shielding" mechanism: muscovite-bordered quartz
//       grains deform homogeneously via GBS WITHOUT recrystallizing —
//       easy basal slip in Ms relaxes strain-rate gradients at phase
//       boundaries → less qz dislocation glide → less qz CPO development.
//   Holyoke & Tullis (2006) Tectonophys.: only ~10 % interconnected mica
//       is required to make mica the dominant rheological phase.
//
// m.u.d. → Skemer M-index conversion: M = (m.u.d. − 1) / 14, documented
// in misc/aniso_fstrain/notes/case15_calibration.md.
//   M_init_Skemer = (6.02 − 1) / 14 ≈ 0.359
// An alternative Watson-MC theoretical conversion gives M_init ≈ 0.26;
// both are bracketed estimates pending a dual-metric paper for quartz.
//
// Slope = 15.0 mirrors case-8 (companion case; same δ-multiplier class).
// Hansen+14 olivine slope = 24.5 is mineral-specific and over-predicts
// for quartz CPO; 15.0 is the defensible middle value for quartz-bearing
// polyphase (3 ≤ slope ≤ 50 across literature, see case-8 docstring).
//
// SCOPE — case-15 is BRACKETED valid only for X_ms ∈ [5%, 25%]:
//   X_ms ∈ [5%, 25%]:  case-15 valid; use this DECAY form.
//   X_ms <  5% :        NOT case-15; use case 3 (Pennacchioni+10
//                        audit-control) or case 6 (Blackford+24 lab)
//                        for near-monophase quartz, saturating-exp.
//   X_ms > 25% :        NOT case-15; FALL BACK to case 8 (mica-dominated
//                        regime; mica CPO is the relevant rheological
//                        signal, not quartz).
// SCOPE is enforced as a DOCUMENTATION contract (no runtime check —
// MDLIB has no per-marker X_ms field today; user picks aniso_db number
// based on their geological setup).
//
// Pooled fit collapses X_ms ∈ {5, 10, 25} % into a single γ-decay; γ_decay
// = 15.35 is the volume-weighted average across the three mica fractions
// sampled. Per-X_ms γ_decay would refine further (~1 h LSQ per fraction;
// deferred).
//
// Caveats:
//   - Asymptote δ(γ→∞) = 1.0 (full CPO erasure) is unphysical — quartz
//     CPO should saturate at a non-zero floor. Tokle+23 γ ∈ [0, 4] doesn't
//     reach the floor regime, so the calibration is bracketed only within
//     γ ∈ [0, ~5]. Flag as scope caveat for γ > 5.
//   - Bracketed pooled fit: γ_max / γ_decay = 4 / 15.35 = 0.26 — saturation
//     not reached in the substrate, hence γ_decay uncertainty is 19 %
//     while M_init is constrained to 2 %.
//
//   γ_eff(FS_AR) = √FS_AR − 1/√FS_AR
//   M(γ)         = 0.359 · exp(−γ/15.35)         // DECAY form
//   δ            = 15.0 · M + 1
static double anisoDelta_QuartzMicaPolyphase( double FS_AR ) {
    const double s           = sqrt(FS_AR);
    const double gamma_eff   = s - 1.0 / s;
    const double M_init      = 0.359;
    const double gamma_decay = 15.35;
    const double slope       = 15.0;
    const double M           = M_init * exp(-gamma_eff / gamma_decay);
    return slope * M + 1.0;
}

void ReadDataAnisotropy( mat_prop* mat, params* model, int k, int number, scale* scaling ) {
    // Strain-dependent viscous-anisotropy database (used when ani_fstrain == 2).
    // Populates the per-phase function pointer mat->aniso_delta_fn[k] with a
    // mineral-specific δ(FS_AR) closed form. Each mineral picks its own
    // functional form here; AnisotropyRoutines.c stays mineral-agnostic.
    // Matches the ReadDataPowerLaw/ReadDataGSE pattern: switch on `number`,
    // populate per-phase arrays, log the citation.

    int success = 0;

    // -----------------------------------------------------------------------
    //  ani_fstrain == 3 — Boneh et al. 2021 DiSRX grain-boundary-migration
    //  kinetics constants → per-phase material database (mat->ani_relax_*).
    // -----------------------------------------------------------------------
    // These set the temperature-dependent δ-relaxation timescale
    // τ_relax = L_relax / V consumed by DeltaRelaxationTau (AnisotropyRoutines.c).
    // They are MINERAL-SPECIFIC. The values below are the OLIVINE calibration:
    // Boneh's Q, M0 are the *olivine* DiSRX grain-boundary-mobility values, and
    // the Δρ range is the *olivine* dislocation-density range. We set the
    // olivine default for EVERY phase here so any aniso_db case inherits it;
    // the olivine cases (1 = Hansen, 7 = damped) are the calibrated ones. A
    // future calcite/quartz aniso_db case can OVERRIDE these per-case below —
    // non-olivine DiSRX kinetics are future work. Stored in SI units (the
    // consumer DeltaRelaxationTau converts the scaled T / L_relax / R to SI and
    // expects these constants already in SI).
    //
    // MODELLING LEAP (caveat, design D7c): Boneh et al. 2021 (G3, "The effect
    // of discontinuous static recrystallization on olivine crystallographic
    // preferred orientation") establish the *kinetics* of discontinuous static
    // recrystallization (DiSRX) — the grain-boundary-migration velocity V and
    // the time for an annealed grain to grow — and they show DiSRX
    // *qualitatively* weakens CPO. They do NOT give a closed-form δ(t)
    // relaxation law. The construction in RheologyParticles.c —
    // "δ relaxes toward the isotropic limit (δ → 1) on the DiSRX grain-growth
    // timescale τ_relax = L_relax / V" — is THIS MODEL'S OWN, not Boneh's. It
    // is a defensible first-order coupling (when grain boundaries sweep a
    // marker faster, its inherited CPO is overwritten faster), but it is a
    // modelling choice; the quantitative δ(t) curve is not lab-validated. See
    // misc/aniso_fstrain/notes/ for the derivation and the tiered validation
    // plan.
    mat->ani_relax_Q[k]  = 133.0e3;   // activation energy [J/mol] (Boneh 2021: 133 ± 18 kJ/mol)
    mat->ani_relax_mu[k] = 50.0e9;    // shear modulus [Pa]       (Boneh 2021: μ = 50 GPa)
    mat->ani_relax_b[k]  = 0.6e-9;    // Burgers vector [m]       (Boneh 2021: b = 0.6 nm)
    // Grain-boundary mobility prefactor.
    // UNITS NOTE (caveat, design D7d / O1): Boneh et al. 2021 prints the
    // prefactor units as "m³/s·J", which is dimensionally inconsistent with its
    // own Eq. 1 (V = M_b · ΣF would then give a *rate*, not a *velocity*). The
    // Boneh et al. 2017 companion paper's Appendix A prints the dimensionally-
    // correct m³/(N·s) ≡ m⁴·J⁻¹·s⁻¹, and with that unit
    // V = M0·exp(−Q/RT)·μ·b²·Δρ comes out in m/s and reproduces Boneh's
    // published "weeks-to-years at ~1300 °C" anchor (10 µm → 1 mm growth in
    // ~0.23–23 yr for Δρ = 10¹³–10¹¹ m⁻²). We therefore use the corrected unit.
    // This is committed, not open.
    mat->ani_relax_M0[k] = 2.0e-11;   // mobility prefactor [m^4 J^-1 s^-1] (= m^3 N^-1 s^-1)
    // Δρ proxy reference range, tied to Boneh's dislocation-density range
    // (~10¹¹–10¹³ m⁻²). The Δρ driving force is proxied from the per-marker
    // accumulated dislocation-creep strain strain_pwl (design D3, O2) inside
    // DeltaRhoProxy (AnisotropyRoutines.c); MDOODZ tracks no dislocation
    // density, so Δρ is derived on the fly from a field it already carries.
    // Saturating exponential mapping: Δρ → drho_min at zero stored strain,
    // → drho_max for strain ≫ eps_ref; result is always finite and in
    // [drho_min, drho_max].
    //
    // TWO CAVEATS, documented inline (design D7b):
    //  (a) MONOTONICITY / NO-RECOVERY: strain_pwl is accumulated dislocation-
    //      creep strain — it can only increase. A real dislocation density
    //      *recovers* (drops) once deformation slows; this proxy cannot. So Δρ
    //      here is a "sticky" upper-ish estimate. For the post-deformation
    //      relaxation use case this stickiness is acceptable (and arguably
    //      correct: the stored strain energy that drives DiSRX is itself slow
    //      to recover), but a marker that deformed long ago then sat quiescent
    //      will still report a high Δρ and therefore relax *faster* than a true
    //      ρ-evolution model.
    //  (b) ABSOLUTE STORED STRAIN vs. BOUNDARY CONTRAST: Boneh's Δρ is the
    //      dislocation-density *contrast across a migrating grain boundary* —
    //      the thermodynamic driving force for that specific boundary.
    //      strain_pwl is a marker's *absolute* accumulated stored strain
    //      energy. The two coincide only when the marker is significantly more
    //      strained than its neighbourhood (a strong gradient drives boundary
    //      migration into it); for a uniformly-strained region the true driving
    //      contrast is smaller than this proxy implies. A τ_II-piezometer proxy
    //      was considered and rejected as the primary choice (it collapses to
    //      zero exactly when deformation stops — wrong for the post-deformation
    //      use case); a piezometer-cap refinement is left open (design D3).
    mat->ani_relax_drho_min[k] = 1.0e10;  // [m^-2]  floor (truly-undeformed mantle, lowered from
                                          //   1e11 to better match annealed-mantle Δρ estimates;
                                          //   Karato 2008 §5.5, Bai et al. 1991 single-crystal floor)
    mat->ani_relax_drho_max[k] = 1.0e13;  // [m^-2]  saturation (highly-strained marker)
    mat->ani_relax_eps_ref[k]  = 4.0;     // [-]     reference saturation strain (Hansen γ_e ≈ 3.96 scale)

    // STRAIN-RATE GATE: DiSRX is a static/post-deformation mechanism (every Boneh+21
    // M0/Q data point comes from hot-press or static-anneal experiments — Karato 1989,
    // Cooper-Kohlstedt 1984, Speciale et al. 2020). Concurrent dislocation creep
    // suppresses DiSRX (Hansen+12/14/16 build CPO at eII ≈ 1e-5 with no observed
    // DiSRX). The gate factor f = 1/(1+(eII/eps_max)^2) divides τ_relax so τ_eff → ∞
    // when actively deforming. Default 1e-13 s^-1 splits the "transient event" regime
    // (mantle wedge / shear-zone activation) from quiescent mantle. Per-phase override
    // via the .txt file.
    mat->ani_relax_eps_max[k]  = 1.0e-13; // [s^-1]  strain-rate gate threshold

    switch ( abs(number) ) {

        case 0 :
            LOG_INFO("should not be here");
            success      = 0;
            break;

        case 1 :
            LOG_INFO("Olivine viscous anisotropy - Hansen+12/+14/+16 (combined 38-point fit):");
            mat->aniso_delta_fn[k]     = anisoDelta_HansenOlivine;
            mat->aniso_delta_fn_inv[k] = aniso_delta_inv_hansen;
            success                = 1;
            break;

        case 2 :
            LOG_INFO("Calcite viscous anisotropy - Pieri+01 / Barnhoorn+04 (15-point combined T-averaged fit, delta_inf = 3.46):");
            mat->aniso_delta_fn[k] = anisoDelta_Calcite;
            success                = 1;
            break;

        case 3 :
            LOG_INFO("Quartz viscous anisotropy - Pennacchioni+10 (17-point natural-shear-zone fit, ODF J, delta_inf = 3.25, T~500C prism-a):");
            mat->aniso_delta_fn[k] = anisoDelta_Quartz;
            success                = 1;
            break;

        case 4 :
            LOG_INFO("Calcite viscous anisotropy - mid-low-T (T <= 650 degC, Barnhoorn+04 7-point fit, delta_inf = 2.33, SR-recryst regime):");
            mat->aniso_delta_fn[k] = anisoDelta_Calcite_LowT;
            success                = 1;
            break;

        case 5 :
            LOG_INFO("Calcite viscous anisotropy - high-T (T > 650 degC, Pieri+01 + Barnhoorn+04(727C) 8-point fit, delta_inf = 5.40, GBM-recryst regime):");
            mat->aniso_delta_fn[k] = anisoDelta_Calcite_HighT;
            success                = 1;
            break;

        case 6 :
            LOG_INFO("Quartz viscous anisotropy - Blackford+24 (38-point Northern Snake Range fit, pfJ of c-axis, delta_inf = 3.20, plane-strain):");
            mat->aniso_delta_fn[k] = anisoDelta_Quartz_Blackford;
            success                = 1;
            break;

        case 7 :
            LOG_INFO("Damped olivine viscous anisotropy - Tasaka+Boneh+Kumamoto+Bernard (93-point non-Hansen fit, delta_inf = 5.90, natural/pre-CPO/wet regime):");
            mat->aniso_delta_fn[k] = anisoDelta_HansenOlivine_Damped;
            // δ-relaxation length default for ani_fstrain==3 (Boneh et al. 2021
            // DiSRX kinetics). Sentinel logic mirrors case-9 d-coupling: parser
            // stores -1 when user didn't write ani_relax_length; we leave it at
            // -1 so the use site inherits gs_ref[phase] (the per-phase annealed
            // grain-size proxy). A user-set positive value overrides gs_ref.
            if ( mat->ani_relax_length[k] < 0.0 ) {
                mat->ani_relax_length[k] = -1.0;  // inherit gs_ref at the use site
                LOG_INFO("  case 7 default δ-relaxation length: ani_relax_length = -1 (inherit gs_ref)");
            }
            success                = 1;
            break;

        case 8 :
            LOG_INFO("Quartz-mica polyphase BRACKETED estimate - granitoid mylonite (Tokle+23 + Dempsey+11 + Holyoke&Tullis06 + Bos&Spiers02, NOT direct LSQ; comparable rigor to pwlv=16 Ranalli granulite, delta_inf = 7.0):");
            mat->aniso_delta_fn[k] = anisoDelta_QuartzMicaBracketed;
            success                = 1;
            break;

        case 9 :
            LOG_INFO("Plagioclase polyphase BRACKETED - gabbro/anorthosite/mafic granulite (Mehl&Hirth08, Marti+18, Barreiro+15, Stunitz+03; delta_inf = 1.96):");
            mat->aniso_delta_fn[k] = anisoDelta_PlagioclaseBracketed;
            // d-coupling default: plagioclase CPO weakens via DisGBS at small
            // grain size. aniso_d_threshold = 100 µm calibrated against
            // Mehl & Hirth (2008) Table 1 (reproduces the M ≈ 0.04 polyphase-
            // analog at d = 30 µm and the M ≈ 0.12 plateau at d ≥ 100 µm).
            //
            // STRUCTURAL LIMIT: the d-only modifier cannot reproduce M&H's
            // bimodal distribution within the population-overlap region
            // [30, 135] µm — both monophase and polyphase samples coexist
            // there at similar d. A complete model would require
            // f = f(d, X_cpx, σ, T) — beyond case 9 today.
            //
            // User overrides: aniso_d_threshold and aniso_d_decay are
            // per-phase fields readable in the .txt config. Parser stores -1
            // when the user didn't write an explicit aniso_d_threshold line;
            // the default below applies only in that case.
            if ( mat->aniso_d_threshold[k] < 0.0 ) {
                mat->aniso_d_threshold[k] = 1.0e-4 / scaling->L;  // 100 µm → scaled
                mat->aniso_d_decay[k]     = 1.0;
                LOG_INFO("  case 9 default d-coupling: aniso_d_threshold = 1e-4 m (100 um), aniso_d_decay = 1.0 (Mehl & Hirth 2008 calibration)");
            }
            success                = 1;
            break;

        case 13 :
            LOG_INFO("Calcite-superplastic BRACKETED - Solnhofen GBS-dominant regime (Schmid 1976/77, De Bresser 2002; delta_inf=1.12; for dislocation-creep regime use case 2/4/5):");
            mat->aniso_delta_fn[k] = anisoDelta_Calcite_Superplastic;
            success                = 1;
            break;

        case 15 :
            LOG_INFO("Quartz-mica polyphase BRACKETED, quartz perspective - granitoid mylonite (Tokle+23, decay form; valid for X_ms in [5, 25] vol pct mica):");
            mat->aniso_delta_fn[k] = anisoDelta_QuartzMicaPolyphase;
            success                = 1;
            break;
    }

    if ( success==0 ) { LOG_ERR("Error: Non existing aniso_db number %d", number); exit(12);}
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
