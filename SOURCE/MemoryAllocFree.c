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

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void* DoodzMalloc( size_t size ) {

    void *pointer;

#ifdef _MKL_
    pointer = mkl_malloc( size, 64  );
#else
    pointer = malloc( size );
#endif

    return pointer;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void* DoodzCalloc( int length, size_t type_size ) {

    void *pointer;

#ifdef _MKL_
    pointer = mkl_calloc( length, type_size, 64  );
#else
    pointer = calloc( length, type_size );
#endif

    return pointer;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void* DoodzRealloc( void *pointer, size_t new_size ) {

    void *new;

#ifdef _MKL_
    new = mkl_realloc( pointer, new_size );
#else
    new = realloc( pointer, new_size );
#endif

    return new;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DoodzFree( void *pointer ) {

#ifdef _MKL_
    mkl_free( pointer );
#else
    free( pointer );
#endif

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocMat( SparseMat *Mat, int nnz ) {
    Mat->Ic  = DoodzCalloc((Mat->neq+1), sizeof(int));
    Mat->J   = DoodzCalloc(nnz, sizeof(int));
    Mat->A   = DoodzCalloc(nnz, sizeof(double));
    Mat->bbc = DoodzCalloc(Mat->neq, sizeof(double));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeMat( SparseMat *Mat ) {
    DoodzFree(Mat->Ic);
    DoodzFree(Mat->J);
    DoodzFree(Mat->A);
    DoodzFree(Mat->bbc);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocatePhaseDiagrams( params* model ) {
    model->PDMnT   = DoodzCalloc( model->num_PD, sizeof(int));
    model->PDMnP   = DoodzCalloc( model->num_PD, sizeof(int));
    model->PDMTmin = DoodzCalloc( model->num_PD, sizeof(double));
    model->PDMTmax = DoodzCalloc( model->num_PD, sizeof(double));
    model->PDMPmin = DoodzCalloc( model->num_PD, sizeof(double));
    model->PDMPmax = DoodzCalloc( model->num_PD, sizeof(double));
    model->PDMrho  = DoodzCalloc( model->num_PD, sizeof(double*));
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreePhaseDiagrams( params* model  ) {
    int k;
    
    for (k=0; k<model->num_PD; k++) {
        DoodzFree( model->PDMrho[k] );
    }
    DoodzFree( model->PDMrho );

    DoodzFree( model->PDMnT   );
    DoodzFree( model->PDMnP   );
    DoodzFree( model->PDMTmin );
    DoodzFree( model->PDMTmax );
    DoodzFree( model->PDMPmin );
    DoodzFree( model->PDMPmax );
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SAlloc( SparseMat *Mat, int neq ) {
    Mat->d   = DoodzCalloc(neq, sizeof(double)); // diagonal
    Initialise1DArrayDouble( Mat->d ,  neq, 1.0 );
    Mat->b   = DoodzCalloc(neq, sizeof(double)); // rhs
    Mat->x   = DoodzCalloc(neq, sizeof(double)); // solution
    Mat->F   = DoodzCalloc(neq, sizeof(double)); // residual
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SFree( SparseMat *Mat ) {
    DoodzFree(Mat->d);
    DoodzFree(Mat->x);
    DoodzFree(Mat->F);
    DoodzFree(Mat->b);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PartAlloc( markers *particles, params* model  ) {

    // Allocate particle arrays
    particles->x          = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->z          = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->Vx         = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->Vz         = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->P          = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->sxxd       = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->szzd       = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->sxz        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->dsxxd      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->dszzd      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->dsxz       = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->syy        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->dsyy        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));

    particles->phase      = DoodzCalloc( particles->Nb_part_max,sizeof(int));
    particles->generation = DoodzCalloc( particles->Nb_part_max,sizeof(int));
    particles->dual       = DoodzCalloc( particles->Nb_part_max,sizeof(int));
    particles->progress   = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));

    particles->divth      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->T          = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->d          = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->phi        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->X          = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->noise      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->rho      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));


    particles->strain     = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->strain_el  = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->strain_pl  = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->strain_pwl = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->strain_exp = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->strain_lin = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    particles->strain_gbs = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));

    particles->intag      = DoodzCalloc( particles->Nb_part_max,sizeof(int));

    if (model->fstrain == 1) {
        particles->Fxx        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->Fxz        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->Fzx        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->Fzz        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        Initialise1DArrayDouble (particles->Fxx,  particles->Nb_part_max, 1.0 );
        Initialise1DArrayDouble (particles->Fzz,  particles->Nb_part_max, 1.0 );
    }

    if (model->rec_T_P_x_z == 1) {
        particles->T0        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->P0        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->z0        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->x0        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->Tmax      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->Pmax      = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    }

    if (model->aniso == 1) {
        particles->dnx       = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->dnz       = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->nx        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
        particles->nz        = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
    }
    
//    particles->ddivth = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
//    particles->dT     = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
//    particles->dP     = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
//    particles->dd     = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
//    particles->dX     = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
//    particles->dphi   = DoodzCalloc( particles->Nb_part_max,sizeof(DoodzFP));
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PartFree( markers *particles, params* model ) {

    // Deallocate particle arrays
    DoodzFree(particles->x);
    DoodzFree(particles->z);
    DoodzFree(particles->Vx);
    DoodzFree(particles->Vz);
    DoodzFree(particles->P);
    DoodzFree(particles->sxxd);
    DoodzFree(particles->szzd);
    DoodzFree(particles->sxz);
    DoodzFree(particles->dsxxd);
    DoodzFree(particles->dszzd);
    DoodzFree(particles->dsxz);
    DoodzFree(particles->syy);
    DoodzFree(particles->dsyy);
    
    DoodzFree(particles->phase);
    DoodzFree(particles->progress);
    DoodzFree(particles->dual);
    DoodzFree(particles->generation);

    DoodzFree(particles->divth);
    DoodzFree(particles->T);
    DoodzFree(particles->d);
    DoodzFree(particles->phi);
    DoodzFree(particles->X);
    DoodzFree(particles->noise);
        DoodzFree(particles->rho);


    DoodzFree(particles->strain);
    DoodzFree(particles->strain_el);
    DoodzFree(particles->strain_pl);
    DoodzFree(particles->strain_pwl);
    DoodzFree(particles->strain_exp);
    DoodzFree(particles->strain_lin);
    DoodzFree(particles->strain_gbs);

    DoodzFree(particles->intag);

    if (model->fstrain == 1) {
        DoodzFree(particles->Fxx);
        DoodzFree(particles->Fxz);
        DoodzFree(particles->Fzx);
        DoodzFree(particles->Fzz);
    }

    if (model->rec_T_P_x_z == 1) {
        DoodzFree(particles->T0);
        DoodzFree(particles->P0);
        DoodzFree(particles->x0);
        DoodzFree(particles->z0);
        DoodzFree(particles->Tmax);
        DoodzFree(particles->Pmax);
    }

    if (model->aniso == 1) {
        DoodzFree(particles->dnx);
        DoodzFree(particles->dnz);
        DoodzFree(particles->nx);
        DoodzFree(particles->nz);
    }
    
//    DoodzFree(particles->ddivth);
//    DoodzFree(particles->dT);
//    DoodzFree(particles->dP);
//    DoodzFree(particles->dd);
//    DoodzFree(particles->dX);
//    DoodzFree(particles->dphi);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GridAlloc ( grid* mesh, params* model ) {

    // Allocate all grid-based arrays

    int k, Nx, Nz, NxVz, NzVx;

    Nx       = model->Nx;
    Nz       = model->Nz;
    mesh->Nx = model->Nx;
    mesh->Nz = model->Nz;
    NxVz     = Nx + 1;
    NzVx     = Nz + 1;

    printf("Allocation of grid arrays !\n");
    
    mesh->FreeSurfW_s = (double*) DoodzCalloc(  mesh->Nx   * mesh->Nz    ,sizeof(double));
    mesh->FreeSurfW_n = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));

    mesh->xg_coord  = (double*) DoodzCalloc( mesh->Nx ,sizeof(double));
    mesh->zg_coord  = (double*) DoodzCalloc( mesh->Nz ,sizeof(double));
    mesh->xg_coord_ext  = (double*) DoodzCalloc( (mesh->Nx+2) ,sizeof(double));
    mesh->zg_coord_ext  = (double*) DoodzCalloc( (mesh->Nz+2) ,sizeof(double));
    mesh->xg_coord0 = (double*) DoodzCalloc( mesh->Nx ,sizeof(double));
    mesh->zg_coord0 = (double*) DoodzCalloc( mesh->Nz ,sizeof(double));
    mesh->xc_coord  = (double*) DoodzCalloc( (mesh->Nx-1) ,sizeof(double));
    mesh->zc_coord  = (double*) DoodzCalloc( (mesh->Nz-1) ,sizeof(double));
    mesh->xvz_coord = (double*) DoodzCalloc( (mesh->Nx+1) ,sizeof(double));
    mesh->zvx_coord = (double*) DoodzCalloc( (mesh->Nz+1) ,sizeof(double));

    mesh->eta_s     = (double*) DoodzCalloc(  mesh->Nx   * mesh->Nz    ,sizeof(double));
    mesh->eta_n     = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));
    mesh->rho_s     = (double*) DoodzCalloc(  mesh->Nx   * mesh->Nz    ,sizeof(double));
    mesh->rho_n     = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));
    mesh->rho0_n    = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));

    mesh->u     = (double*) DoodzCalloc( (mesh->Nx)  *(mesh->Nz+1) ,sizeof(double));
    mesh->ru    = (double*) DoodzCalloc( (mesh->Nx)  *(mesh->Nz+1) ,sizeof(double));
    mesh->rhs_u = (double*) DoodzCalloc( (mesh->Nx)  *(mesh->Nz+1) ,sizeof(double));
    mesh->gx    = (double*) DoodzCalloc( (mesh->Nx)  *(mesh->Nz+1) ,sizeof(double));
    mesh->v     = (double*) DoodzCalloc( (mesh->Nx+1)*(mesh->Nz)   ,sizeof(double));
    mesh->rv    = (double*) DoodzCalloc( (mesh->Nx+1)*(mesh->Nz)   ,sizeof(double));
    mesh->rhs_v = (double*) DoodzCalloc( (mesh->Nx+1)*(mesh->Nz)   ,sizeof(double));
    mesh->gz    = (double*) DoodzCalloc( (mesh->Nx+1)*(mesh->Nz)   ,sizeof(double));
    mesh->p     = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));
    mesh->rp    = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));
    mesh->rhs_p = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));

    mesh->BCu.type     = (char*)   DoodzCalloc( (mesh->Nx)  *(mesh->Nz+1)   ,sizeof(char));
    mesh->BCv.type     = (char*)   DoodzCalloc( (mesh->Nx+1)*(mesh->Nz)     ,sizeof(char));
    mesh->BCp.type     = (char*)   DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1)   ,sizeof(char));
    mesh->BCp_exp.type = (char*)   DoodzCalloc( (mesh->Nx+1)*(mesh->Nz+1)   ,sizeof(char));
    Initialise1DArrayChar(   mesh->BCp_exp.type, (mesh->Nx+1)*(mesh->Nz+1), -1  );
    mesh->BCu.val      = (double*) DoodzCalloc( (mesh->Nx)  *(mesh->Nz+1) ,sizeof(double));
    mesh->BCv.val      = (double*) DoodzCalloc( (mesh->Nx+1)*(mesh->Nz)   ,sizeof(double));
    mesh->BCp.val      = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1) ,sizeof(double));
    mesh->BCp_exp.val  = (double*) DoodzCalloc( (mesh->Nx+1)*(mesh->Nz+1) ,sizeof(double));

    //-------------------------------------------------------------------------------------------------//

    // Allocate 2D array : phase percentage per cell center (fine mesh only)
    mesh->phase_perc_n = (double**) DoodzCalloc( model->Nb_phases,sizeof(double*));
    for ( k=0; k<model->Nb_phases; k++ ) {
        mesh->phase_perc_n[k] = (double*) DoodzCalloc( (mesh->Nx-1)*(mesh->Nz-1), sizeof(double));
    }

    // Allocate 2D array : phase percentage per cell vertices (fine mesh only)
    mesh->phase_perc_s = (double**) DoodzCalloc( model->Nb_phases,sizeof(double*));
    for ( k=0; k<model->Nb_phases; k++ ) {
        mesh->phase_perc_s[k] = (double*) DoodzCalloc( (mesh->Nx)*(mesh->Nz), sizeof(double));
    }

    //--------------------------------------------------//

    // Roger's
    mesh->roger_x    = DoodzCalloc (Nx*(Nz+1),sizeof(double));
    mesh->roger_z    = DoodzCalloc ((Nx+1)*Nz,sizeof(double));
    mesh->div_u      = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->div_u_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->div_u_el   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->div_u_pl   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->div_u_r    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->Qrho       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    // Solution arrays
    mesh->u_in       = DoodzCalloc (Nx*NzVx,sizeof(double));
    mesh->v_in       = DoodzCalloc (NxVz*Nz,sizeof(double));
    mesh->p_in       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->p_corr    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->u_start    = DoodzCalloc (Nx*NzVx,sizeof(double));
    mesh->v_start    = DoodzCalloc (NxVz*Nz,sizeof(double));
    mesh->p_start    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    // Derived quantities
    mesh->sxxd       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->szzd       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->sxz        = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->exxd       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->ezzd       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->exz        = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    // Right-hand side viscoelastic coefficient
    mesh->mu_s   = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->mu_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->VE_s   = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->VE_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->sxxd0  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->szzd0  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->sxz0   = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->sxxd0_s= DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->szzd0_s= DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->sxz0_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->exxd_s= DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->ezzd_s= DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->exz_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->sxz_n = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    mesh->u_adv  = DoodzCalloc (Nx*NzVx,sizeof(double));
    mesh->v_adv  = DoodzCalloc (NxVz*Nz,sizeof(double));

    mesh->eta_phys_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->eta_phys_s = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    mesh->strain_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->strain_s   = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    // Energy equation
    mesh->Cv       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->kx       = DoodzCalloc (Nx*NzVx,sizeof(double));
    mesh->kz       = DoodzCalloc (NxVz*Nz,sizeof(double));
    mesh->kc_x     = DoodzCalloc (Nx*NzVx,sizeof(double));
    mesh->kc_z     = DoodzCalloc (NxVz*Nz,sizeof(double));
    mesh->Qr       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->rhs_t    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->BCt.type = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(char));
    mesh->BCt.typW = DoodzCalloc ((Nz-1),sizeof(char));
    mesh->BCt.typE = DoodzCalloc ((Nz-1),sizeof(char));
    mesh->BCt.typS = DoodzCalloc ((Nx-1),sizeof(char));
    mesh->BCt.typN = DoodzCalloc ((Nx-1),sizeof(char));
    mesh->BCt.val  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->BCt.valW = DoodzCalloc ((Nz-1),sizeof(double));
    mesh->BCt.valE = DoodzCalloc ((Nz-1),sizeof(double));
    mesh->BCt.valS = DoodzCalloc ((Nx-1),sizeof(double));
    mesh->BCt.valN = DoodzCalloc ((Nx-1),sizeof(double));
    mesh->Wdiss    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->Wel      = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->Wtot     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    
    // Chemical diffusion
    mesh->BCc.type = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(char));
    mesh->BCc.typW = DoodzCalloc ((Nz-1),sizeof(char));
    mesh->BCc.typE = DoodzCalloc ((Nz-1),sizeof(char));
    mesh->BCc.typS = DoodzCalloc ((Nx-1),sizeof(char));
    mesh->BCc.typN = DoodzCalloc ((Nx-1),sizeof(char));
    mesh->BCc.val  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->BCc.valW = DoodzCalloc ((Nz-1),sizeof(double));
    mesh->BCc.valE = DoodzCalloc ((Nz-1),sizeof(double));
    mesh->BCc.valS = DoodzCalloc ((Nx-1),sizeof(double));
    mesh->BCc.valN = DoodzCalloc ((Nx-1),sizeof(double));

    // Grid tagging (free surf)
    mesh->BCg.type = DoodzCalloc ((Nx)*(Nz),sizeof(char));
    mesh->BCg.val  = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    // Volumetric deformation
    mesh->p_lith   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->p_lith0  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->dp       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->alp      = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->bet_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->bet_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    // Grid indices
    mesh->kvx       = DoodzCalloc (Nx*NzVx,sizeof(int));
    mesh->lvx       = DoodzCalloc (Nx*NzVx,sizeof(int));
    mesh->kvz       = DoodzCalloc (NxVz*Nz,sizeof(int));
    mesh->lvz       = DoodzCalloc (NxVz*Nz,sizeof(int));
    mesh->kp        = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(int));
    mesh->lp        = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(int));
    mesh->kn        = DoodzCalloc ((Nx)*(Nz),sizeof(int));
    mesh->ln        = DoodzCalloc ((Nx)*(Nz),sizeof(int));

    // Array containing the number of particles in each cell
    mesh->nb_part_cell = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(int));
    mesh->nb_part_vert = DoodzCalloc ((Nx)*(Nz),sizeof(int));

    // Inertia : keep density on velocity points and velocity increments
    //    mesh->rhoVx  = DoodzMalloc (Nx*NzVx*sizeof(double));
    //    mesh->rhoVz  = DoodzMalloc (NxVz*Nz*sizeof(double));
    mesh->VzVx   = DoodzCalloc (Nx*NzVx,sizeof(double));
    mesh->VxVz   = DoodzCalloc (NxVz*Nz,sizeof(double));

    // Time-series
    mesh->Uthermal_time   = DoodzCalloc (model->Nt+1,sizeof(double)); // +1 to include step 000
    mesh->Uelastic_time   = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->Work_time       = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->Time_time       = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->Short_time      = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->P_mean_time     = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->T_mean_time     = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->Tii_mean_time   = DoodzCalloc (model->Nt+1,sizeof(double));
    mesh->Eii_mean_time   = DoodzCalloc (model->Nt+1,sizeof(double));

    mesh->T          = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->T0_n       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->divth0_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->dT         = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->eII_el     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->exx_el     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->ezz_el     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->exz_el     = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->exx_diss   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->ezz_diss   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->exz_diss   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->exz_n_el   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->exz_n_diss = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    mesh->comp_cells = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    // To remove
//    mesh->exx_pl     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->exz_pl     = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
//    mesh->exx_pwl_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->exz_pwl_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->A2_pwl_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->exx_pwl_s  = DoodzCalloc ((Nx)*(Nz),sizeof(double));
//    mesh->exz_pwl_s  = DoodzCalloc ((Nx)*(Nz),sizeof(double));
//    mesh->eII_pwl_s  = DoodzCalloc ((Nx)*(Nz),sizeof(double));
//    mesh->A2_pwl_s   = DoodzCalloc ((Nx)*(Nz),sizeof(double));
//    mesh->eii_s      = DoodzCalloc ((Nx)*(Nz),sizeof(double));
//    mesh->eii_n      = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->tii0_s     = DoodzCalloc ((Nx)*(Nz),sizeof(double));
//    mesh->tii0_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->exz_n_pl   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    mesh->eII_pl     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
//    mesh->eII_pl_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->eII_pwl    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->eII_exp    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->eII_lin    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->eII_gbs    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->eII_cst    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->d_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->d0_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->phi_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->phi0_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->cell_min_z = DoodzCalloc( (Nz-1), sizeof(double));
    mesh->cell_max_z = DoodzCalloc( (Nz-1), sizeof(double));
    mesh->vert_min_z = DoodzCalloc( (Nz-0), sizeof(double));
    mesh->vert_max_z = DoodzCalloc( (Nz-0), sizeof(double));

    // New arrays for plastic strain softening
    mesh->C_s        = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->fric_s     = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->dil_s      = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->C_n        = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->fric_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->dil_n      = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    // For Newton iterations
    mesh->D11_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->D12_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->D13_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->D14_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->D21_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->D22_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->D23_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->D24_n     = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->D31_s     = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->D32_s     = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->D33_s     = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->D34_s     = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    mesh->drhodp_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    
    mesh->detadexx_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->detadezz_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->detadgxz_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->detadp_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    
    mesh->ddivpdexx_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->ddivpdezz_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->ddivpdgxz_n  = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->ddivpdp_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));

    mesh->detadexx_s  = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->detadezz_s  = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->detadgxz_s  = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->detadp_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    mesh->d0_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->phi0_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->T_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    mesh->P_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    if ( model->aniso == 1 ) mesh->nx_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    if ( model->aniso == 1 ) mesh->nx_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    if ( model->aniso == 1 ) mesh->nz_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    if ( model->aniso == 1 ) mesh->nz_s    = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    if ( model->aniso == 1 ) mesh->nx0_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    if ( model->aniso == 1 ) mesh->nx0_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    if ( model->aniso == 1 ) mesh->nz0_n   = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    if ( model->aniso == 1 ) mesh->nz0_s   = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    if ( model->aniso == 1 ) mesh->FS_AR_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    if ( model->aniso == 1 ) mesh->FS_AR_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));
    if ( model->aniso == 1 ) Initialise1DArrayDouble (mesh->FS_AR_n,  (Nx-1)*(Nz-1), 1.0 );
    if ( model->aniso == 1 ) Initialise1DArrayDouble (mesh->FS_AR_s,  (Nx-0)*(Nz-0), 1.0 );
    if ( model->aniso == 1 ) mesh->aniso_factor_n = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    if ( model->aniso == 1 ) mesh->aniso_factor_s = DoodzCalloc ((Nx-0)*(Nz-0),sizeof(double));

    mesh->p0_n       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->p0_s       = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    // reaction
    mesh->X0_n       = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->X_n        = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->noise_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->X0_s       = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->X_s        = DoodzCalloc ((Nx)*(Nz),sizeof(double));
    mesh->noise_s    = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    mesh->OverS_n    = DoodzCalloc ((Nx-1)*(Nz-1),sizeof(double));
    mesh->OverS_s    = DoodzCalloc ((Nx)*(Nz),sizeof(double));

    printf("Memory succesfully allocated : \nDiantre, que c'est bon!\n");

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void GridFree( grid* mesh, params* model ) {

    // Free all grid-based arrays

    int k;
    
    DoodzFree(mesh->FreeSurfW_s);
    DoodzFree(mesh->FreeSurfW_n);

    DoodzFree(mesh->roger_x);
    DoodzFree(mesh->roger_z);
    DoodzFree(mesh->div_u);
    DoodzFree(mesh->div_u_s);
    DoodzFree(mesh->div_u_el);
    DoodzFree(mesh->div_u_pl);
    DoodzFree(mesh->div_u_r);
    DoodzFree(mesh->Qrho);

	// Solution arrays
    DoodzFree(mesh->u_in);
    DoodzFree(mesh->v_in);
    DoodzFree(mesh->p_in);
    DoodzFree(mesh->p_corr);

    DoodzFree(mesh->u_start);
    DoodzFree(mesh->v_start);
    DoodzFree(mesh->p_start);

	// Stresses and strain rates
    DoodzFree(mesh->sxxd);
    DoodzFree(mesh->szzd);
    DoodzFree(mesh->sxz);
    DoodzFree(mesh->exxd);
    DoodzFree(mesh->ezzd);
    DoodzFree(mesh->exz);

    // Previously multigrid structures
    DoodzFree(mesh->xg_coord);
    DoodzFree(mesh->zg_coord);
    DoodzFree(mesh->xg_coord_ext);
    DoodzFree(mesh->zg_coord_ext);
    DoodzFree(mesh->xg_coord0);
    DoodzFree(mesh->zg_coord0);
    DoodzFree(mesh->xc_coord);
    DoodzFree(mesh->zc_coord);
    DoodzFree(mesh->xvz_coord);
    DoodzFree(mesh->zvx_coord);
    DoodzFree(mesh->eta_s);
    DoodzFree(mesh->eta_n);

    DoodzFree(mesh->rho_s);
    DoodzFree(mesh->rho_n);
    DoodzFree(mesh->rho0_n);
    DoodzFree(mesh->u);
    DoodzFree(mesh->ru);
    DoodzFree(mesh->v);
    DoodzFree(mesh->rv);
    DoodzFree(mesh->p);
    DoodzFree(mesh->rp);

    DoodzFree(mesh->rhs_u);
    DoodzFree(mesh->gx);
    DoodzFree(mesh->rhs_v);
    DoodzFree(mesh->gz);
    DoodzFree(mesh->rhs_p);

    DoodzFree(mesh->BCu.val);
    DoodzFree(mesh->BCv.val);
    DoodzFree(mesh->BCp.val);
    DoodzFree(mesh->BCp_exp.val);
    DoodzFree(mesh->BCu.type);
    DoodzFree(mesh->BCv.type);
    DoodzFree(mesh->BCp.type);
    DoodzFree(mesh->BCp_exp.type);

    // Viscoelasticity
    DoodzFree(mesh->mu_s);
    DoodzFree(mesh->mu_n);
    DoodzFree(mesh->VE_s);
    DoodzFree(mesh->VE_n);

    DoodzFree(mesh->sxxd0);
    DoodzFree(mesh->szzd0);
    DoodzFree(mesh->sxz0);
    DoodzFree(mesh->sxxd0_s);
    DoodzFree(mesh->szzd0_s);
    DoodzFree(mesh->sxz0_n);

    DoodzFree(mesh->exxd_s);
    DoodzFree(mesh->ezzd_s);
    DoodzFree(mesh->exz_n);
    DoodzFree(mesh->sxz_n);

    DoodzFree(mesh->u_adv);
    DoodzFree(mesh->v_adv);
    DoodzFree(mesh->eta_phys_s);
    DoodzFree(mesh->eta_phys_n);

    DoodzFree(mesh->strain_s);
    DoodzFree(mesh->strain_n);

    // Energy equation
    DoodzFree(mesh->Cv);
    DoodzFree(mesh->kx);
    DoodzFree(mesh->kz);
    DoodzFree(mesh->kc_x);
    DoodzFree(mesh->kc_z);
    DoodzFree(mesh->Qr);
    DoodzFree(mesh->rhs_t);
    DoodzFree(mesh->BCt.val);
    DoodzFree(mesh->BCt.valW);
    DoodzFree(mesh->BCt.valE);
    DoodzFree(mesh->BCt.valS);
    DoodzFree(mesh->BCt.valN);
    DoodzFree(mesh->BCt.type);
    DoodzFree(mesh->BCt.typW);
    DoodzFree(mesh->BCt.typE);
    DoodzFree(mesh->BCt.typS);
    DoodzFree(mesh->BCt.typN);
    DoodzFree(mesh->Wdiss);
    DoodzFree(mesh->Wel);
    DoodzFree(mesh->Wtot);
    
    // Chemical diffusion
    DoodzFree(mesh->BCc.val);
    DoodzFree(mesh->BCc.valW);
    DoodzFree(mesh->BCc.valE);
    DoodzFree(mesh->BCc.valS);
    DoodzFree(mesh->BCc.valN);
    DoodzFree(mesh->BCc.type);
    DoodzFree(mesh->BCc.typW);
    DoodzFree(mesh->BCc.typE);
    DoodzFree(mesh->BCc.typS);
    DoodzFree(mesh->BCc.typN);

    // Grid tagging
    DoodzFree(mesh->BCg.val);
    DoodzFree(mesh->BCg.type);

    // Volumetric deformation
    DoodzFree(mesh->p_lith);
    DoodzFree(mesh->p_lith0);
    DoodzFree(mesh->dp);
    DoodzFree(mesh->alp);
    DoodzFree(mesh->bet_n);
    DoodzFree(mesh->bet_s);

    // Grid indices
    DoodzFree( mesh->kvx );
    DoodzFree( mesh->lvx );
    DoodzFree( mesh->kvz );
    DoodzFree( mesh->lvz );
    DoodzFree( mesh->kp  );
    DoodzFree( mesh->lp  );
    DoodzFree( mesh->kn  );
    DoodzFree( mesh->ln  );

    // Phases cell centers
    for ( k=0; k<model->Nb_phases; k++ ) {
        DoodzFree( mesh->phase_perc_n[k] );
    }
    DoodzFree( mesh->phase_perc_n );

    // Phases cell vertices
    for ( k=0; k<model->Nb_phases; k++ ) {
        DoodzFree( mesh->phase_perc_s[k] );
    }
    DoodzFree( mesh->phase_perc_s );

    // Inertia
    DoodzFree(mesh->VxVz);
    DoodzFree(mesh->VzVx);

    // Array containing the number of particles in each cell
    DoodzFree(mesh->nb_part_cell);
    DoodzFree(mesh->nb_part_vert);

    // Time-series
    DoodzFree(mesh->Uthermal_time);
    DoodzFree(mesh->Uelastic_time);
    DoodzFree(mesh->Work_time);
    DoodzFree(mesh->Time_time);
    DoodzFree(mesh->Short_time);
    DoodzFree(mesh->P_mean_time);
    DoodzFree(mesh->T_mean_time);
    DoodzFree(mesh->Tii_mean_time);
    DoodzFree(mesh->Eii_mean_time);

    DoodzFree(mesh->T);
    DoodzFree(mesh->T0_n);
    DoodzFree(mesh->divth0_n);
    DoodzFree(mesh->dT);
    DoodzFree(mesh->eII_el);
    DoodzFree(mesh->exx_el);
    DoodzFree(mesh->ezz_el);
    DoodzFree(mesh->exz_el);
    DoodzFree(mesh->exx_diss);
    DoodzFree(mesh->ezz_diss);
    DoodzFree(mesh->exz_diss);

    DoodzFree(mesh->eII_pl);
    DoodzFree(mesh->eII_pwl);
    DoodzFree(mesh->eII_exp);
    DoodzFree(mesh->eII_lin);
    DoodzFree(mesh->eII_gbs);
    DoodzFree(mesh->eII_cst);

    DoodzFree( mesh->comp_cells );

    DoodzFree(mesh->exz_n_el);
    DoodzFree(mesh->exz_n_diss);

    DoodzFree(mesh->d_n);
    DoodzFree(mesh->d0_n);
    DoodzFree(mesh->phi_n);
    DoodzFree(mesh->phi0_n);

    DoodzFree(mesh->cell_min_z);
    DoodzFree(mesh->cell_max_z);
    DoodzFree(mesh->vert_min_z);
    DoodzFree(mesh->vert_max_z);

    // New arrays for plastic strain softening
    DoodzFree(mesh->C_s);
    DoodzFree(mesh->fric_s);
    DoodzFree(mesh->dil_s);
    DoodzFree(mesh->C_n);
    DoodzFree(mesh->fric_n);
    DoodzFree(mesh->dil_n);

    // For Newton iterations
    DoodzFree(mesh->D11_n);
    DoodzFree(mesh->D12_n);
    DoodzFree(mesh->D13_n);
    DoodzFree(mesh->D14_n);

    DoodzFree(mesh->D21_n);
    DoodzFree(mesh->D22_n);
    DoodzFree(mesh->D23_n);
    DoodzFree(mesh->D24_n);

    DoodzFree(mesh->D31_s);
    DoodzFree(mesh->D32_s);
    DoodzFree(mesh->D33_s);
    DoodzFree(mesh->D34_s);
    
    DoodzFree(mesh->drhodp_n);

    DoodzFree(mesh->detadexx_n);
    DoodzFree(mesh->detadezz_n);
    DoodzFree(mesh->detadgxz_n);
    DoodzFree(mesh->detadp_n);
    
    DoodzFree(mesh->ddivpdexx_n);
    DoodzFree(mesh->ddivpdezz_n);
    DoodzFree(mesh->ddivpdgxz_n);
    DoodzFree(mesh->ddivpdp_n);

    DoodzFree(mesh->detadexx_s);
    DoodzFree(mesh->detadezz_s);
    DoodzFree(mesh->detadgxz_s);
    DoodzFree(mesh->detadp_s);

    DoodzFree(mesh->d0_s);
    DoodzFree(mesh->phi0_s);
    DoodzFree(mesh->T_s);
    DoodzFree(mesh->P_s);

    if ( model->aniso == 1 ) DoodzFree(mesh->nx_n);
    if ( model->aniso == 1 ) DoodzFree(mesh->nx_s);
    if ( model->aniso == 1 ) DoodzFree(mesh->nz_n);
    if ( model->aniso == 1 ) DoodzFree(mesh->nz_s);
    if ( model->aniso == 1 ) DoodzFree(mesh->nx0_n);
    if ( model->aniso == 1 ) DoodzFree(mesh->nx0_s);
    if ( model->aniso == 1 ) DoodzFree(mesh->nz0_n);
    if ( model->aniso == 1 ) DoodzFree(mesh->nz0_s);
    if ( model->aniso == 1 ) DoodzFree(mesh->FS_AR_n);
    if ( model->aniso == 1 ) DoodzFree(mesh->FS_AR_s);
    if ( model->aniso == 1 ) DoodzFree(mesh->aniso_factor_n);
    if ( model->aniso == 1 ) DoodzFree(mesh->aniso_factor_s);

    // reac
    DoodzFree(mesh->X_n);
    DoodzFree(mesh->noise_n);
    DoodzFree(mesh->X_s);
    DoodzFree(mesh->noise_s);
    DoodzFree(mesh->X0_n);
    DoodzFree(mesh->X0_s);
    DoodzFree(mesh->p0_n);
    DoodzFree(mesh->p0_s);

    DoodzFree(mesh->OverS_n);
    DoodzFree(mesh->OverS_s);


}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeSparseSystems( int IsJacobianUsed, int IsDecoupledSolver, SparseMat* Stokes, SparseMat* StokesA, SparseMat* StokesB, SparseMat* StokesC, SparseMat* StokesD, SparseMat* JacobA, SparseMat* JacobB, SparseMat* JacobC, SparseMat* JacobD ) {

    // This free memory that was dynamically allocated for the currect sparse systems (Stokes and/or Jacobian)

    if ( IsDecoupledSolver == 0 ) FreeMat( Stokes );
    else {
        FreeMat( StokesA );
        FreeMat( StokesB );
        FreeMat( StokesC );
        FreeMat( StokesD );
    }
    if ( IsJacobianUsed == 1 ) {
        FreeMat( JacobA );
        FreeMat( JacobB );
        FreeMat( JacobC );
        FreeMat( JacobD );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
