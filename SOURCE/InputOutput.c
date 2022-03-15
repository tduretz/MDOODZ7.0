// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2019  MDOODZ Developper team
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
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "ctype.h"
#include "header_MDOODZ.h"

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

void ExplicitDiffusion2D( double* rho, int nT, int nP, double dT, double dP, scale *scaling ) {
    // Apply some diffusion...
    const double Kdiff   = 1e0;              // diffusivity
    const double dt_exp  =  MINV(dP*dP,dT*dT)/Kdiff/2.1; // explicit time step
    const int    n_steps = 200;              // number of steps
    double       qW, qE, qS, qN;   
    double *rho0   = DoodzCalloc( nT*nP, sizeof(double)); 

    for ( int it=0; it<n_steps; it++) {
        ArrayEqualArray( rho0, rho, nT*nP);
        for (int iz = 1; iz<nT-1; iz++) {
            for (int ix = 1; ix<nP-1; ix++) {
                int c = nT*iz + ix;
                qW      = - Kdiff*(rho0[c] - rho0[c- 1])/dT;
                qE      = - Kdiff*(rho0[c+1] - rho0[ix])/dT;
                qS      = - Kdiff*(rho0[c] - rho0[c-nT])/dP;
                qN      = - Kdiff*(rho0[c+nT] - rho0[c])/dP;
                rho[c]  = rho0[c] - dt_exp*(qE - qW)/dT - dt_exp*(qN - qS)/dP;
            }
        }
    }
    DoodzFree(rho0);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void LoadIniParticles( char* name, markers* particles, grid* mesh, markers *topo_chain, markers *topo_chain_ini, params *model, scale scaling ) {

    //char *name;
    int s1=0, s2=0, s3=0, s4=0;
    int ss1=sizeof(int), ss2=sizeof(double), ss3=sizeof(DoodzFP), ss4=sizeof(char);
    int k, l, Nx=model->Nx, Nz=model->Nz, c;
    int Ncx = Nx-1, Ncz = Nz-1;

    //---------------------------------------------------------------------------------------------------------//
    FILE *file;

    // Filename
    //    asprintf(&name, "Breakpoint%05d.dat", model->step);
    //asprintf(&name, "IniParticles.dat");

    file = fopen(name, "rb");
   if ( file == NULL ) {
       printf("Error opening %s, this file is probably corrupted or non-existant\nExiting...\n", name);
       fclose(file);
       free(name);
       exit(1);
   }

    fread( &s1, 1, 1, file);
    fread( &s2, 1, 1, file);
    fread( &s3, 1, 1, file);
    fread( &s4, 1, 1, file);

    if (model->free_surf == 1) {
        fread(&topo_chain->Nb_part,   s1,                   1, file );
        fread( topo_chain->x,         s3, topo_chain->Nb_part, file );
        fread( topo_chain->z,         s3, topo_chain->Nb_part, file );
        //fread( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
        //fread( topo_chain->Vz,        s3, topo_chain->Nb_part, file );
    }

    int Nb_part0 = particles->Nb_part;
    fread( &particles->Nb_part,  s1, 1, file);
    fread( particles->x,    s3, particles->Nb_part, file);
    fread( particles->z,    s3, particles->Nb_part, file);
    //particles->Nb_part = Nb_part0;
    //fread( particles->P,    s3, particles->Nb_part, file);
    //fread( particles->Vx,   s3, particles->Nb_part, file);
    //fread( particles->Vz,   s3, particles->Nb_part, file);
    //fread( particles->phi,  s3, particles->Nb_part, file);
    //fread( particles->X  ,  s3, particles->Nb_part, file);
    fread( particles->phase, s1, particles->Nb_part, file);
    fread( particles->dual, s1, particles->Nb_part, file);

    fclose(file);
    free(name);
    //---------------------------------------------------------------------------------------------------------//

    if (model->free_surf == 1) {
        topo_chain_ini->Nb_part = topo_chain->Nb_part;
        // note: the topo_chain should be scaled as well
#pragma omp parallel for shared( topo_chain, model, scaling )
        for ( k=0; k<topo_chain->Nb_part; k++ ) {
            topo_chain->x[k]     /=scaling.L;
            topo_chain->z[k]     /=scaling.L;


            topo_chain_ini->x[k] = topo_chain->x[k];
            topo_chain_ini->z[k] = topo_chain->z[k];
            //particles->P[k]     /=scaling.S;
            //particles->Vx[k]    /=scaling.V;
            //particles->Vz[k]    /=scaling.V;
            //particles->phi[k]   /=1.0;
            //particles->X[k]     /=1.0;
        }
    }
#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /=scaling.L;
        particles->z[k]     /=scaling.L;
        //particles->P[k]     /=scaling.S;
        //particles->Vx[k]    /=scaling.V;
        //particles->Vz[k]    /=scaling.V;
        //particles->phi[k]   /=1.0;
        //particles->X[k]     /=1.0;
    }


}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DeletePreviousBreakpoint( int step, int writer_step ) {
    char *name, *command, *new_name;
    int success;
    if ((step- 2*writer_step)>1) {
        asprintf(&new_name, "Breakpoint%05d.dat", step- 0*writer_step);
        asprintf(&name,     "Breakpoint%05d.dat", step- 2*writer_step);
        asprintf(&command, "mv %s %s", name, new_name );
        success = system( command );
        printf("File %s replaced by %s\n", name, new_name);
    //    success = system( command );
    //        success = remove( name );
    //    if ( success!=-1 ) {
    //        printf("File %s was successfully deleted\n", name);
    //    }
        if ( success!=-1 ) printf("File %s was successfully renamed\n", name);
        else printf("File %s was not successfully renamed\n", name);
        free(name);
        free(new_name);
        free(command);
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void LoadBreakpointParticles( markers *particles, grid* mesh, markers *topo_chain, markers *topo_chain_ini, params *model, surface* topo, surface* topo_ini, scale scaling ) {

    char *name;
    int s1=0, s2=0, s3=0, s4=0;
    int k, l, Nx=model->Nx, Nz=model->Nz, c;
    int Ncx = Nx-1, Ncz = Nz-1;

    //---------------------------------------------------------------------------------------------------------//

    FILE *file;

    // Filename
    if( model->delete_breakpoints == -1 ) asprintf(&name, "BreakpointXXXXX.dat");
    else asprintf(&name, "Breakpoint%05d.dat", model->step);
    file = fopen(name, "rb");
    if ( file == NULL ) {
       printf("Error opening %s, this file is probably corrupted or non-existant\nExiting...\n", name);
       fclose(file);
       free(name);
       exit(1);
    }
    fread( &s1, 1, 1, file);
    fread( &s2, 1, 1, file);
    fread( &s3, 1, 1, file);
    fread( &s4, 1, 1, file);
    fread( &model->dt,   s2, 1, file);
    fread( &model->dt0,  s2, 1, file);
    fread( &model->time, s2, 1, file);
    fread( &model->xmin, s2, 1, file);
    fread( &model->xmax, s2, 1, file);
    fread( &model->zmin, s2, 1, file);
    fread( &model->zmax, s2, 1, file);


    if (model->free_surf == 1) {
        fread(&topo_chain->Nb_part,   s1,                   1, file );
        fread( topo_chain->x,         s3, topo_chain->Nb_part, file );
        fread( topo_chain->z,         s3, topo_chain->Nb_part, file );
        fread( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
        fread( topo_chain->Vz,        s3, topo_chain->Nb_part, file );

        fread(&topo_chain_ini->Nb_part,   s1,                       1, file );
        fread( topo_chain_ini->x,         s3, topo_chain_ini->Nb_part, file );
        fread( topo_chain_ini->z,         s3, topo_chain_ini->Nb_part, file );
        fread( topo_chain_ini->Vx,        s3, topo_chain_ini->Nb_part, file );
        fread( topo_chain_ini->Vz,        s3, topo_chain_ini->Nb_part, file );
    }

    fread( &particles->Nb_part,  s1, 1, file);
    fread( particles->x,    s3, particles->Nb_part, file);
    fread( particles->z,    s3, particles->Nb_part, file);
    fread( particles->P,    s3, particles->Nb_part, file);
    fread( particles->Vx,   s3, particles->Nb_part, file);
    fread( particles->Vz,   s3, particles->Nb_part, file);
    fread( particles->phi,  s3, particles->Nb_part, file);
    fread( particles->X  ,  s3, particles->Nb_part, file);
    fread( particles->noise, s3, particles->Nb_part, file);
    fread( particles->phase, s1, particles->Nb_part, file);
    fread( particles->dual, s1, particles->Nb_part, file);

    if (model->iselastic == 1) {
        fread( particles->sxxd,   s3, particles->Nb_part, file );
        fread( particles->szzd,   s3, particles->Nb_part, file );
        fread( particles->sxz,    s3, particles->Nb_part, file );
        fread( particles->dsxxd,   s3, particles->Nb_part, file );
        fread( particles->dszzd,   s3, particles->Nb_part, file );
        fread( particles->dsxz,    s3, particles->Nb_part, file );

        fread( mesh->eta_n, s3, (Nx-1)*(Nz-1), file );
        fread( mesh->VE_n, s3, (Nx-1)*(Nz-1), file );
        fread( mesh->eta_s, s3, (Nx)*(Nz), file );
        fread( mesh->VE_s, s3, (Nx)*(Nz), file );
    }

    fread( particles->strain,     s3, particles->Nb_part, file);
    fread( particles->strain_el,  s3, particles->Nb_part, file);
    fread( particles->strain_pl,  s3, particles->Nb_part, file);
    fread( particles->strain_pwl, s3, particles->Nb_part, file);
    fread( particles->strain_exp, s3, particles->Nb_part, file);
    fread( particles->strain_lin, s3, particles->Nb_part, file);
    fread( particles->strain_gbs, s3, particles->Nb_part, file);
    fread( particles->d         , s3, particles->Nb_part, file);

    if (model->fstrain == 1) {
        fread( particles->Fxx         , s3, particles->Nb_part, file);
        fread( particles->Fxz         , s3, particles->Nb_part, file);
        fread( particles->Fzx         , s3, particles->Nb_part, file);
        fread( particles->Fzz         , s3, particles->Nb_part, file);
    }

    if (model->aniso == 1) {
        fread( particles->nx         , s3, particles->Nb_part, file);
        fread( particles->nz         , s3, particles->Nb_part, file);
    }

    if (model->rec_T_P_x_z == 1) {
        fread( particles->T0         , s3, particles->Nb_part, file);
        fread( particles->P0         , s3, particles->Nb_part, file);
        fread( particles->x0         , s3, particles->Nb_part, file);
        fread( particles->z0         , s3, particles->Nb_part, file);
        fread( particles->Tmax       , s3, particles->Nb_part, file);
        fread( particles->Pmax       , s3, particles->Nb_part, file);
    }

    fread( particles->divth, s3, particles->Nb_part, file);
    fread( particles->T,     s3, particles->Nb_part, file);
    
    fread( mesh->p0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->T0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->divth0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->X0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->d0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->phi0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->rho0_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->sxxd, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->szzd, s3, (Nx-1)*(Nz-1), file );

    fread( mesh->X0_s, s3, (Nx)*(Nz),     file );
    fread( mesh->sxz, s3, (Nx)*(Nz),     file );

    fread( mesh->eta_phys_n, s3, (Nx-1)*(Nz-1), file );
    fread( mesh->eta_phys_s, s3, (Nx)*(Nz),     file );

    fread( mesh->u_in, s3, (Nx*(Nz+1)), file );
    fread( mesh->v_in, s3, (Nz*(Nx+1)), file );
    fread( mesh->p_in, s3, ((Nz-1)*(Nx-1)), file );
    fread( mesh->p_lith,  s3, ((Nz-1)*(Nx-1)), file );
    fread( mesh->p_lith0, s3, ((Nz-1)*(Nx-1)), file );

    fread( &mesh->Uthermal,  s3, 1, file );
    fread( &mesh->Work,   s3, 1, file );
    fread( &model->L0, s3, 1, file );

    // This is to avoid any problem with restarting - more data needs to be stored

        // Topo related
        if (model->free_surf == 1) {
            fread( topo->height,  s3,  Nx, file );
            fread( topo->height0, s3,  Nx, file );
            fread( topo->a,       s3, Ncx, file );
            fread( topo->a0,      s3, Ncx, file );
            fread( topo->b,       s3, Ncx, file );
            fread( topo->b0,      s3, Ncx, file );
            fread( topo->vx,      s3,  Nx, file );
            fread( topo->vz,      s3,Nx+1, file );
            fread( topo_ini->height,  s3,  Nx, file );
            fread( topo_ini->height0, s3,  Nx, file );
            fread( topo_ini->a,       s3, Ncx, file );
            fread( topo_ini->a0,      s3, Ncx, file );
            fread( topo_ini->b,       s3, Ncx, file );
            fread( topo_ini->b0,      s3, Ncx, file );
            fread( topo_ini->vx,      s3,  Nx, file );
            fread( topo_ini->vz,      s3,Nx+1, file );
        }
        // Cell flags
        fread( mesh->BCt.type,     s4,  Ncx*Ncz,   file );
        fread( mesh->BCp.type,  s4,  Ncx*Ncz,   file );
        fread( mesh->BCu.type,  s4,  Nx*(Nz+1), file );
        fread( mesh->BCv.type,  s4,  (Nx+1)*Nz, file );
        fread( mesh->BCg.type,     s4,  Nx *Nz ,   file );
        fread( mesh->BCp.val,   s3,  Ncx*Ncz,   file );
        fread( mesh->BCu.val,   s3,  Nx*(Nz+1), file );
        fread( mesh->BCv.val,   s3,  (Nx+1)*Nz, file );
        fread( mesh->BCg.val,      s3,  Nx *Nz ,   file );
//    MinMaxArray(mesh->BCg.val, 1.0, Nx *Nz, "mesh->BCg.val");


     // Phase proportions
        printf("Loading phase proportions - Nb_phases = %d:\n", model->Nb_phases);
//        fread( &model->Nb_phases,  s1,        1,   file );
        for ( k=0; k<  model->Nb_phases; k++ ) { //model->Nb_phases
            fread( mesh->phase_perc_n[k], s3,  Ncx*Ncz,   file );
            fread( mesh->phase_perc_s[k], s3,  Nx *Nz ,   file );
//            printf("phase %02d:\n", k);
//            MinMaxArray(mesh->phase_perc_n[k], 1.0, Ncx*Ncz, "phase_perc_n");
//            MinMaxArray(mesh->phase_perc_s[k], 1.0, Nx *Nz , "phase_perc_s");
        }
        fread( mesh->nb_part_cell, s1,  Ncx*Ncz,   file );
        fread( mesh->nb_part_vert, s1,  Nx *Nz ,   file );


    fclose(file);
    free(name);

    //---------------------------------------------------------------------------------------------------------//

    // scale
    model->dt   /= scaling.t;
    model->dt0  /= scaling.t;
    model->time /= scaling.t;
    model->xmin /= scaling.L;
    model->xmax /= scaling.L;
    model->zmin /= scaling.L;
    model->zmax /= scaling.L;

#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /=scaling.L;
        particles->z[k]     /=scaling.L;
        particles->P[k]     /=scaling.S;
        particles->Vx[k]    /=scaling.V;
        particles->Vz[k]    /=scaling.V;
        particles->phi[k]   /=1.0;
        particles->X[k]     /=1.0;
        if ( model->iselastic == 1 ) {
            particles->sxxd[k]    /= scaling.S;
            particles->szzd[k]    /= scaling.S;
            particles->sxz[k]     /= scaling.S;
            particles->dsxxd[k]    /= scaling.S;
            particles->dszzd[k]    /= scaling.S;
            particles->dsxz[k]     /= scaling.S;
        }
   
        particles->divth[k]  /= scaling.E;
        particles->T[k]      /= scaling.T;

        if ( model->rec_T_P_x_z == 1) {
            particles->T0[k]  /= scaling.T;
            particles->P0[k]  /= scaling.S;
            particles->x0[k]  /= scaling.L;
            particles->z0[k]  /= scaling.L;
            particles->Tmax[k]/= scaling.T;
            particles->Pmax[k]/= scaling.S;
        }
    }

    // Grid data
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz+1; l++) {
            c = k +l*Nx;
            mesh->u_adv[c] /= scaling.V;
            mesh->u_in[c]  /= scaling.V;
        }
    }

    for (k=0; k<Nx+1; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx+1);
            mesh->v_adv[c] /= scaling.V;
            mesh->v_in[c]  /= scaling.V;
        }
    }

    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c = k +l*(Nx-1);
            mesh->eta_n[c] /= scaling.eta;
            mesh->VE_n[c] /= 1.0;
            mesh->eta_phys_n[c] /= scaling.eta;
            mesh->p_in[c]    /= scaling.S;
            mesh->p_lith[c]  /= scaling.S;
            mesh->p_lith0[c] /= scaling.S;
            
            mesh->p0_n[c]      /= scaling.S;
            mesh->T0_n[c]      /= scaling.T;
            mesh->divth0_n[c]  /= scaling.E;
            mesh->X0_n[c]      /= 1.0;
            mesh->d0_n[c]      /= scaling.L;
            mesh->phi0_n[c]    /= 1.0;
            mesh->rho0_n[c]    /= scaling.rho;
            mesh->sxxd[c]      /= scaling.S;
            mesh->szzd[c]      /= scaling.S;
        }
    }

    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx);
            mesh->eta_phys_s[c] /= scaling.eta;
            mesh->eta_s[c] /= scaling.eta;
            mesh->VE_s[c] /= 1.0;
            
            mesh->X0_s[c]      /= 1.0;
            mesh->sxz[c]       /= scaling.S;
        }
    }

    if (model->free_surf == 1) {
        for (k=0; k<topo_chain->Nb_part; k++) {
            topo_chain->x[k]  /= scaling.L;
            topo_chain->z[k]  /= scaling.L;
            topo_chain->Vx[k] /= scaling.V;
            topo_chain->Vz[k] /= scaling.V;
        }
        for (k=0; k<topo_chain_ini->Nb_part; k++) {
            topo_chain_ini->x[k] /= scaling.L;
            topo_chain_ini->z[k] /= scaling.L;
            topo_chain_ini->Vx[k] /= scaling.V;
            topo_chain_ini->Vz[k] /= scaling.V;
        }
        //    }

        // This is to avoid any problem with restarting - more data needs to be stored
        for (k=0;k<mesh->Nx;k++) {
            topo->height[k]      /= scaling.L;
            topo->height0[k]     /= scaling.L;
            topo->vx[k]          /= scaling.V;
            topo_ini->height[k]  /= scaling.L;
            topo_ini->height0[k] /= scaling.L;
            topo_ini->vx[k]      /= scaling.V;
        }

        for (k=0;k<mesh->Nx+1;k++) {
            topo->vz[k]          /= scaling.V;
            topo_ini->vz[k]      /= scaling.V;
        }

        for (k=0;k<mesh->Nx-1;k++) {
            topo->b0[k] /= scaling.L;
            topo->b[k]  /= scaling.L;
            topo_ini->b0[k] /= scaling.L;
            topo_ini->b[k]  /= scaling.L;
        }
    }

    mesh->Uthermal   /= (scaling.rhoE*scaling.L*scaling.L);
    mesh->Work    /= (scaling.rhoE*scaling.L*scaling.L);
    model->L0  /= (scaling.L);

//    MinMaxArray(particles->T, scaling.T, particles->Nb_part, "T part");

}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MakeBreakpointParticles( markers *particles,  grid* mesh, markers *topo_chain, markers *topo_chain_ini, params model, surface *topo, surface* topo_ini, scale scaling ) {

    char *name;
    int s1=sizeof(int), s2=sizeof(double), s3=sizeof(DoodzFP), s4=sizeof(char);
    int k, l, Nx=model.Nx, Nz=model.Nz, c;
    int Ncx = Nx-1, Ncz = Nz-1;

    // Scale such that dimensional data goes to the breakpoint file
    model.dt   *= scaling.t;
    model.dt0  *= scaling.t;
    model.time *= scaling.t;
    model.xmin *= scaling.L;
    model.xmax *= scaling.L;
    model.zmin *= scaling.L;
    model.zmax *= scaling.L;

#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]      *= scaling.L;
        particles->z[k]      *= scaling.L;
        particles->P[k]      *= scaling.S;
        particles->Vx[k]     *= scaling.V;
        particles->Vz[k]     *= scaling.V;
        particles->phi[k]    *= 1.0;
        particles->X[k]      *= 1.0;
        if (model.iselastic == 1) {
            particles->sxxd[k] *= scaling.S;
            particles->szzd[k] *= scaling.S;
            particles->sxz[k]  *= scaling.S;
            particles->dsxxd[k] *= scaling.S;
            particles->dszzd[k] *= scaling.S;
            particles->dsxz[k]  *= scaling.S;
        }
        particles->divth[k]  *= scaling.E;
        particles->T[k]         *= scaling.T;
        particles->d[k]        *= scaling.L;

        if ( model.rec_T_P_x_z == 1) {
            particles->T0[k]  *= scaling.T;
            particles->P0[k]  *= scaling.S;
            particles->x0[k]  *= scaling.L;
            particles->z0[k]  *= scaling.L;
            particles->Tmax[k]*= scaling.T;
            particles->Pmax[k]*= scaling.S;
        }
//        particles->ddivth[k] *= scaling.E;
//        particles->drho[k]   *= scaling.rho;
//        particles->dd[k]     *= scaling.L;
//        particles->dP[k]     *= scaling.S;
//        particles->dT[k]     *= scaling.T;
    }

    // grid data
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz+1; l++) {
            c = k +l*Nx;
            mesh->u_adv[c] *= scaling.V;
            mesh->u_in[c]  *= scaling.V;
        }
    }

    for (k=0; k<Nx+1; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx+1);
            mesh->v_adv[c] *= scaling.V;
            mesh->v_in[c]  *= scaling.V;
        }
    }

    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c = k +l*(Nx-1);
            mesh->eta_n[c] *= scaling.eta;
            mesh->VE_n[c] *= 1.0;
            mesh->eta_phys_n[c] *= scaling.eta;
            mesh->p_in[c]    *= scaling.S;
            mesh->p_lith[c]  *= scaling.S;
            mesh->p_lith0[c] *= scaling.S;
            
            mesh->p0_n[c]      *= scaling.S;
            mesh->T0_n[c]      *= scaling.T;
            mesh->divth0_n[c]  *= scaling.E;
            mesh->X0_n[c]      *= 1.0;
            mesh->d0_n[c]      *= scaling.L;
            mesh->phi0_n[c]    *= 1.0;
            mesh->rho0_n[c]    *= scaling.rho;
            mesh->sxxd0[c]     *= scaling.S;
            mesh->szzd0[c]     *= scaling.S;
        }
    }

    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx);
            mesh->eta_phys_s[c] *= scaling.eta;
            mesh->eta_s[c] *= scaling.eta;
            mesh->VE_s[c] *= 1.0;
            
            mesh->X0_s[c]      *= 1.0;
            mesh->sxz0[c]      *= scaling.S;
        }
    }

    if (model.free_surf == 1) {
        for (k=0; k<topo_chain->Nb_part; k++) {
            topo_chain->x[k] *= scaling.L;
            topo_chain->z[k] *= scaling.L;
            topo_chain->Vx[k] *= scaling.V;
            topo_chain->Vz[k] *= scaling.V;
        }
        for (k=0; k<topo_chain_ini->Nb_part; k++) {
            topo_chain_ini->x[k] *= scaling.L;
            topo_chain_ini->z[k] *= scaling.L;
            topo_chain_ini->Vx[k] *= scaling.V;
            topo_chain_ini->Vz[k] *= scaling.V;
        }


        // This is to avoid any problem with restarting - more data needs to be stored
        for (k=0;k<mesh->Nx;k++) {
            topo->height[k]      *= scaling.L;
            topo->height0[k]     *= scaling.L;
            topo->vx[k]          *= scaling.V;
            topo_ini->height[k]  *= scaling.L;
            topo_ini->height0[k] *= scaling.L;
            topo_ini->vx[k]      *= scaling.V;
        }

        for (k=0;k<mesh->Nx+1;k++) {
            topo->vz[k]          *= scaling.V;
            topo_ini->vz[k]      *= scaling.V;
        }

        for (k=0;k<mesh->Nx-1;k++) {
            topo->b0[k] *= scaling.L;
            topo->b[k]  *= scaling.L;
            topo_ini->b0[k] *= scaling.L;
            topo_ini->b[k]  *= scaling.L;
        }
    }

    mesh->Uthermal *= (scaling.rhoE*scaling.L*scaling.L);
    mesh->Work  *= (scaling.rhoE*scaling.L*scaling.L);
    model.L0 *= (scaling.L);


    //---------------------------------------------------------------------------------------------------------//

    FILE *file;

    // Filename
    if( model.delete_breakpoints == -1 ) asprintf(&name, "BreakpointXXXXX.dat");
    else asprintf(&name, "Breakpoint%05d.dat", model.step);

    file = fopen(name, "wb");
    if ( file == NULL ) {
        printf("Error opening %s, this file is probably corrupted or non-existant\nExiting...\n", name);
        free(name);
        fclose(file);
        exit(1);
    }
    fwrite( &s1, 1, 1, file);
    fwrite( &s2, 1, 1, file);
    fwrite( &s3, 1, 1, file);
    fwrite( &s4, 1, 1, file);
    fwrite( &model.dt,  s2, 1, file);
    fwrite( &model.dt0, s2, 1, file);
    fwrite( &model.time, s2, 1, file);
    fwrite( &model.xmin, s2, 1, file);
    fwrite( &model.xmax, s2, 1, file);
    fwrite( &model.zmin, s2, 1, file);
    fwrite( &model.zmax, s2, 1, file);

    if (model.free_surf == 1) {
        fwrite( &topo_chain->Nb_part,  s1,                   1, file );
        fwrite( topo_chain->x,         s3, topo_chain->Nb_part, file );
        fwrite( topo_chain->z,         s3, topo_chain->Nb_part, file );
        fwrite( topo_chain->Vx,        s3, topo_chain->Nb_part, file );
        fwrite( topo_chain->Vz,        s3, topo_chain->Nb_part, file );

        fwrite( &topo_chain_ini->Nb_part,  s1,                       1, file );
        fwrite( topo_chain_ini->x,         s3, topo_chain_ini->Nb_part, file );
        fwrite( topo_chain_ini->z,         s3, topo_chain_ini->Nb_part, file );
        fwrite( topo_chain_ini->Vx,        s3, topo_chain_ini->Nb_part, file );
        fwrite( topo_chain_ini->Vz,        s3, topo_chain_ini->Nb_part, file );
    }

    fwrite( &particles->Nb_part,  s1, 1, file);
    fwrite( particles->x,     s3, particles->Nb_part, file);
    fwrite( particles->z,     s3, particles->Nb_part, file);
    fwrite( particles->P,     s3, particles->Nb_part, file);
    fwrite( particles->Vx,    s3, particles->Nb_part, file);
    fwrite( particles->Vz,    s3, particles->Nb_part, file);
    fwrite( particles->phi,   s3, particles->Nb_part, file);
    fwrite( particles->X  ,   s3, particles->Nb_part, file);
    fwrite( particles->noise, s3, particles->Nb_part, file);
    fwrite( particles->phase, s1, particles->Nb_part, file);
    fwrite( particles->dual,  s1, particles->Nb_part, file);

    if (model.iselastic == 1) {
        fwrite( particles->sxxd,   s3, particles->Nb_part, file );
        fwrite( particles->szzd,   s3, particles->Nb_part, file );
        fwrite( particles->sxz,    s3, particles->Nb_part, file );
        fwrite( particles->dsxxd,   s3, particles->Nb_part, file );
        fwrite( particles->dszzd,   s3, particles->Nb_part, file );
        fwrite( particles->dsxz,    s3, particles->Nb_part, file );

        fwrite( mesh->eta_n, s3, (Nx-1)*(Nz-1), file );
        fwrite( mesh->VE_n, s3, (Nx-1)*(Nz-1), file );
        fwrite( mesh->eta_s, s3, (Nx)*(Nz), file );
        fwrite( mesh->VE_s, s3, (Nx)*(Nz), file );
    }

    fwrite( particles->strain,     s3, particles->Nb_part, file);
    fwrite( particles->strain_el,  s3, particles->Nb_part, file);
    fwrite( particles->strain_pl,  s3, particles->Nb_part, file);
    fwrite( particles->strain_pwl, s3, particles->Nb_part, file);
    fwrite( particles->strain_exp, s3, particles->Nb_part, file);
    fwrite( particles->strain_lin, s3, particles->Nb_part, file);
    fwrite( particles->strain_gbs, s3, particles->Nb_part, file);
    fwrite( particles->d         , s3, particles->Nb_part, file);

    if (model.fstrain == 1) {
        fwrite( particles->Fxx         , s3, particles->Nb_part, file);
        fwrite( particles->Fxz         , s3, particles->Nb_part, file);
        fwrite( particles->Fzx         , s3, particles->Nb_part, file);
        fwrite( particles->Fzz         , s3, particles->Nb_part, file);
    }

    if (model.aniso == 1) {
        fwrite( particles->nx         , s3, particles->Nb_part, file);
        fwrite( particles->nz         , s3, particles->Nb_part, file);
    }

    if (model.rec_T_P_x_z == 1) {
        fwrite( particles->T0         , s3, particles->Nb_part, file);
        fwrite( particles->P0         , s3, particles->Nb_part, file);
        fwrite( particles->x0         , s3, particles->Nb_part, file);
        fwrite( particles->z0         , s3, particles->Nb_part, file);
        fwrite( particles->Tmax       , s3, particles->Nb_part, file);
        fwrite( particles->Pmax       , s3, particles->Nb_part, file);
    }

    fwrite( particles->divth, s3, particles->Nb_part, file);
    fwrite( particles->T,     s3, particles->Nb_part, file);
    
//    fwrite( particles->ddivth, s3, particles->Nb_part, file);
//    fwrite( particles->drho,   s3, particles->Nb_part, file);
//    fwrite( particles->dX,     s3, particles->Nb_part, file);
//    fwrite( particles->dphi,   s3, particles->Nb_part, file);
//    fwrite( particles->dd,     s3, particles->Nb_part, file);
//    fwrite( particles->dP,     s3, particles->Nb_part, file);
//    fwrite( particles->dT,     s3, particles->Nb_part, file);
    
    fwrite( mesh->p0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->T0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->divth0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->X0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->d0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->phi0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->rho0_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->sxxd0, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->szzd0, s3, (Nx-1)*(Nz-1), file );

    fwrite( mesh->X0_s, s3, (Nx)*(Nz),     file );
    fwrite( mesh->sxz0, s3, (Nx)*(Nz),     file );

    fwrite( mesh->eta_phys_n, s3, (Nx-1)*(Nz-1), file );
    fwrite( mesh->eta_phys_s, s3, (Nx)*(Nz), file );

    fwrite( mesh->u_in, s3, (Nx*(Nz+1)), file );
    fwrite( mesh->v_in, s3, (Nz*(Nx+1)), file );
    fwrite( mesh->p_in, s3, ((Nx-1)*(Nz-1)), file );
    fwrite( mesh->p_lith,  s3, ((Nx-1)*(Nz-1)), file );
    fwrite( mesh->p_lith0, s3, ((Nx-1)*(Nz-1)), file );

    fwrite( &mesh->Uthermal, s3, 1, file );
    fwrite( &mesh->Work,  s3, 1, file );
    fwrite( &model.L0, s3, 1, file );

    // This is to avoid any problem with restarting - more data needs to be stored
    if (model.free_surf == 1) {
        // Topo related
        fwrite( topo->height,  s3,  Nx, file );
        fwrite( topo->height0, s3,  Nx, file );
        fwrite( topo->a,       s3, Ncx, file );
        fwrite( topo->a0,      s3, Ncx, file );
        fwrite( topo->b,       s3, Ncx, file );
        fwrite( topo->b0,      s3, Ncx, file );
        fwrite( topo->vx,      s3,  Nx, file );
        fwrite( topo->vz,      s3,Nx+1, file );
        fwrite( topo_ini->height,  s3,  Nx, file );
        fwrite( topo_ini->height0, s3,  Nx, file );
        fwrite( topo_ini->a,       s3, Ncx, file );
        fwrite( topo_ini->a0,      s3, Ncx, file );
        fwrite( topo_ini->b,       s3, Ncx, file );
        fwrite( topo_ini->b0,      s3, Ncx, file );
        fwrite( topo_ini->vx,      s3,  Nx, file );
        fwrite( topo_ini->vz,      s3,Nx+1, file );
    }
    // Cell flags
    fwrite( mesh->BCt.type,     s4,  Ncx*Ncz,   file );
    fwrite( mesh->BCp.type,  s4,  Ncx*Ncz,   file );
    fwrite( mesh->BCu.type,  s4,  Nx*(Nz+1), file );
    fwrite( mesh->BCv.type,  s4,  (Nx+1)*Nz, file );
    fwrite( mesh->BCg.type,     s4,  Nx *Nz ,   file );
    fwrite( mesh->BCp.val,   s3,  Ncx*Ncz,   file );
    fwrite( mesh->BCu.val,   s3,  Nx*(Nz+1), file );
    fwrite( mesh->BCv.val,   s3,  (Nx+1)*Nz, file );
    fwrite( mesh->BCg.val,      s3,  Nx *Nz ,   file );

//     MinMaxArray(mesh->BCg.val, 1.0, Nx *Nz, "mesh->BCg.val");

    // Phase proportions
    printf("Loading phase proportions - Nb_phases = %d:\n", model.Nb_phases);
//    fwrite( &model.Nb_phases, s1,        1,   file );
    for ( k=0; k< model.Nb_phases; k++ ) {
        fwrite( mesh->phase_perc_n[k], s3,  Ncx*Ncz,   file );
        fwrite( mesh->phase_perc_s[k], s3,  Nx *Nz ,   file );
//        printf("phase %02d:\n", k);
//        MinMaxArray(mesh->phase_perc_n[k], 1.0, Ncx*Ncz, "phase_perc_n");
//        MinMaxArray(mesh->phase_perc_s[k], 1.0, Nx *Nz, "phase_perc_s");
    }
    fwrite( mesh->nb_part_cell, s1,  Ncx*Ncz,   file );
    fwrite( mesh->nb_part_vert, s1,  Nx *Nz ,   file );
    //    }

    fclose(file);
    free(name);

    //---------------------------------------------------------------------------------------------------------//

    // scale such that dimensions vanish
    model.dt   /= scaling.t;
    model.dt0  /= scaling.t;
    model.time /= scaling.t;
    model.xmin /= scaling.L;
    model.xmax /= scaling.L;
    model.zmin /= scaling.L;
    model.zmax /= scaling.L;

#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /= scaling.L;
        particles->z[k]     /= scaling.L;
        particles->P[k]     /= scaling.S;
        particles->Vx[k]    /= scaling.V;
        particles->Vz[k]    /= scaling.V;
        particles->phi[k]   /= 1.0;
        particles->X[k]     /= 1.0;

        if (model.iselastic == 1) {
            particles->sxxd[k]   /= scaling.S;
            particles->szzd[k]   /= scaling.S;
            particles->sxz[k]    /= scaling.S;
            particles->dsxxd[k]   /= scaling.S;
            particles->dszzd[k]   /= scaling.S;
            particles->dsxz[k]    /= scaling.S;
        }
        particles->divth[k]  /= scaling.E;
        particles->T[k]      /= scaling.T;
        particles->d[k]  /= scaling.L;

        if (model.rec_T_P_x_z == 1) {
            particles->T0[k]  /= scaling.T;
            particles->P0[k]  /= scaling.S;
            particles->x0[k]  /= scaling.L;
            particles->z0[k]  /= scaling.L;
            particles->Tmax[k]/= scaling.T;
            particles->Pmax[k]/= scaling.S;
        }
    }

    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz+1; l++) {
            c = k +l*Nx;
            mesh->u_adv[c] /= scaling.V;
            mesh->u_in[c]  /= scaling.V;
        }
    }

    for (k=0; k<Nx+1; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx+1);
            mesh->v_adv[c] /= scaling.V;
            mesh->v_in[c]  /= scaling.V;
        }
    }

    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c = k +l*(Nx-1);
            mesh->eta_n[c] /= scaling.eta;
            mesh->VE_n[c] /= 1.0;
            mesh->eta_phys_n[c] /= scaling.eta;
            mesh->p_in[c]    /= scaling.S;
            mesh->p_lith[c]  /= scaling.S;
            mesh->p_lith0[c] /= scaling.S;
            
            mesh->p0_n[c]      /= scaling.S;
            mesh->T0_n[c]      /= scaling.T;
            mesh->divth0_n[c]  /= scaling.E;
            mesh->X0_n[c]      /= 1.0;
            mesh->d0_n[c]      /= scaling.L;
            mesh->phi0_n[c]    /= 1.0;
            mesh->rho0_n[c]    /= scaling.rho;
            mesh->sxxd0[c]     /= scaling.S;
            mesh->szzd0[c]     /= scaling.S;
        }
    }

    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            c = k +l*(Nx);
            mesh->eta_phys_s[c] /= scaling.eta;
            mesh->eta_s[c]      /= scaling.eta;
            mesh->VE_s[c]       /= 1.0;
            
            mesh->X0_s[c]      /= 1.0;
            mesh->sxz0[c]      /= scaling.S;
        }
    }

    if (model.free_surf == 1) {
        for (k=0; k<topo_chain->Nb_part; k++) {
            topo_chain->x[k]  /= scaling.L;
            topo_chain->z[k]  /= scaling.L;
            topo_chain->Vx[k] /= scaling.V;
            topo_chain->Vz[k] /= scaling.V;
        }
        for (k=0; k<topo_chain_ini->Nb_part; k++) {
            topo_chain_ini->x[k]  /= scaling.L;
            topo_chain_ini->z[k]  /= scaling.L;
            topo_chain_ini->Vx[k] /= scaling.V;
            topo_chain_ini->Vz[k] /= scaling.V;
        }
        //}

        // This is to avoid any problem with restarting - more data needs to be stored
        for (k=0;k<mesh->Nx;k++) {
            topo->height[k]      /= scaling.L;
            topo->height0[k]     /= scaling.L;
            topo->vx[k]          /= scaling.V;
            topo_ini->height[k]  /= scaling.L;
            topo_ini->height0[k] /= scaling.L;
            topo_ini->vx[k]      /= scaling.V;
        }

        for (k=0;k<mesh->Nx+1;k++) {
            topo->vz[k]          /= scaling.V;
            topo_ini->vz[k]      /= scaling.V;
        }

        for (k=0;k<mesh->Nx-1;k++) {
            topo->b0[k] /= scaling.L;
            topo->b[k]  /= scaling.L;
            topo_ini->b0[k] /= scaling.L;
            topo_ini->b[k]  /= scaling.L;
        }
    }

    mesh->Uthermal /= (scaling.rhoE*scaling.L*scaling.L);
    mesh->Work     /= (scaling.rhoE*scaling.L*scaling.L);
    model.L0       /= (scaling.L);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ReadInputFile( char* fin_name, int *istep, int *irestart, int *writer, int *writer_step, params *model, scale *scaling, mat_prop *materials, markers *particles, Nparams *Nmodel ) {

    //------------------------------------------------------------------------------------------------------------------------------//
    // MODEL PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//

    FILE *fin;
    int k, gsel;

    fin = fopen(fin_name,"rt");
    if (fin == NULL) {
        printf("Setup file '%s' does not exist\nExiting...\n", fin_name);
        fclose(fin);
        exit(1);
    }

    // Simulation start/restart from Breakpoint
    *istep                 = ReadInt2( fin, "istep", 0 );
    *irestart              = ReadInt2( fin, "irestart", 0 );

    if ( *istep == 0 ) {
        *irestart = 0; // Override the restart step number written in the text file (istep)
    }

    // Output
    *writer                = ReadInt2( fin, "writer",          0 );
    if (*writer<0 || *writer>1) {
        printf("'writer' should be set to 1 or 0\n Exiting...\n");
        exit(1);
    }
    *writer_step           = ReadInt2( fin, "writer_step",     1 );
    model->write_markers   = ReadInt2( fin, "writer_markers",  0 );
    model->write_debug     = ReadInt2( fin, "writer_debug",    0 );

    // Input
    model->input_file      = ReadChar( fin, "input_file", "blah.bin");

    // Read scales for non-dimensionalisation
    scaling->eta           = ReadDou2( fin, "eta", 1.0  );
    scaling->L             = ReadDou2( fin, "L",   1.0  );
    scaling->V             = ReadDou2( fin, "V",   1.0  );
    scaling->T             = ReadDou2( fin, "T",   1.0  );
    ScaleMe( scaling );

    // Domain size
    model->Nx              = ReadInt2( fin, "Nx",   10 );
    model->Nz              = ReadInt2( fin, "Nz",   10 );
    model->Nt              = ReadInt2( fin, "Nt",    1 );
    model->xmin            = ReadDou2( fin, "xmin",  1.0  ) / scaling->L; model->xmin0 = model->xmin;
    model->zmin            = ReadDou2( fin, "zmin",  1.0  ) / scaling->L; model->zmin0 = model->zmin;
    model->xmax            = ReadDou2( fin, "xmax",  1.0  ) / scaling->L; model->xmax0 = model->xmax;
    model->zmax            = ReadDou2( fin, "zmax",  1.0  ) / scaling->L; model->zmax0 = model->zmax;
    model->dt              = ReadDou2( fin, "dt",    0.0  ) / scaling->t;
    model->Courant         = ReadDou2( fin, "Courant",       0.5 );
    model->penalty         = ReadDou2( fin, "penalty",      1.0e10 );
    model->abs_tol_div     = ReadDou2( fin, "abs_tol_div", 1.0e-14 );
    model->rel_tol_div     = ReadDou2( fin, "rel_tol_div",  1.0e-5 );
    model->auto_penalty    = ReadDou2( fin, "auto_penalty",    0.0  );
    model->decoupled_solve = ReadInt2( fin, "decoupled_solve",    1 );
    model->diag_scaling    = ReadInt2( fin, "diag_scaling",       1 );
    model->pc_type         = ReadInt2( fin, "pc_type",       0 );
    model->safe_mode       = ReadInt2( fin, "safe_mode",     0 );
    model->safe_dt_div     = ReadDou2( fin, "safe_dt_div",  5.0 );
    model->nstagmax        = ReadInt2( fin, "nstagmax",      3 );
    model->noisy           = ReadInt2( fin, "noisy",         1 );  // Prints a lot of info to standard output
    model->eqn_state       = ReadInt2( fin, "eqn_state",       0 );

    // Switches
    model->initial_part    = ReadInt2( fin, "initial_part",    1 ); // Initial particule distribution, 0: MD4.5 style, 1: MD6.0 style
    model->initial_noise   = ReadInt2( fin, "initial_noise",   0 ); // Add noise on initial marker locations
    model->ismechanical    = ReadInt2( fin, "ismechanical",    1 ); // Activates mechanical solver
    model->advection       = ReadInt2( fin, "advection",       1 ); // Activates advection
    model->dt_constant     = ReadInt2( fin, "dt_constant",     0 ); // Fixed time step
    model->RK              = ReadInt2( fin, "RK",              4 ); // Order of Runge-Kutta advection solver
    model->isperiodic_x    = ReadInt2( fin, "isperiodic_x",    0 ); // Activates periodicity in x
    model->ispureshear_ale = ReadInt2( fin, "ispureshear_ALE", 0 );
    model->isinertial      = ReadInt2( fin, "isinertial",      0 );
    model->iselastic       = ReadInt2( fin, "iselastic",       0 );
    model->isthermal       = ReadInt2( fin, "isthermal",       0 );
    model->line_search     = ReadInt2( fin, "line_search",     0 );
    model->line_search_min = ReadDou2( fin, "line_search_min", 0.0 );
    model->free_surf       = ReadInt2( fin, "free_surf",       0 );
    model->free_surf_stab  = ReadDou2( fin, "free_surf_stab",  0 );
    model->thermal_eq      = ReadInt2( fin, "thermal_eq",      0 );
    model->cooling_time    = ReadDou2( fin, "cooling_time", 1000e6*365.25*3600*24/scaling->t);
    model->subgrid_diff    = ReadInt2( fin, "subgrid_diff",    0 );
    model->shear_heat      = ReadInt2( fin, "shear_heat",      1 );
    model->adiab_heat      = ReadInt2( fin, "adiab_heat",      0 );
    model->surf_processes  = ReadInt2( fin, "surf_processes",  0 ); // 1 = diffusion; 2 = diffusion + sedimentation
    model->surf_remesh     = ReadInt2( fin, "surf_remesh",     1 );
    model->cpc             = ReadInt2( fin, "cpc",             1 );
    model->loc_iter        = ReadInt2( fin, "loc_iter",        1 );
    model->therm_pert      = ReadInt2( fin, "therm_pert",      0 );
    model->fstrain         = ReadInt2( fin, "fstrain",         0 );
    model->rec_T_P_x_z     = ReadInt2( fin, "rec_T_P_x_z",     0 );
    model->delete_breakpoints = ReadInt2( fin, "delete_breakpoints",        1 );
    model->topografix      = ReadInt2( fin, "topografix",      0 );
    model->aniso           = ReadInt2( fin, "aniso",           0 ); // Turns on anisotropy
    model->aniso_fstrain   = ReadInt2( fin, "aniso_fstrain",   0 ); // Make anisotropy factor dependent on finite strain aspect ratio
    model->compressible    = ReadInt2( fin, "compressible",    0 ); // Turns on compressibility
    model->GNUplot_residuals = ReadInt2( fin, "GNUplot_residuals",    0 ); // Activate GNU plot residuals visualisation
    model->shear_style     = ReadInt2( fin, "shear_style",     0 ); // 0: pure shear, 2: periodic simple shear
    model->StressRotation  = ReadInt2( fin, "StressRotation",  1 ); // 0: no stress rotation, 1: analytic rotation, 2: upper convected rate
    model->StressUpdate    = ReadInt2( fin, "StressUpdate",    0 );
    model->polar           = ReadInt2( fin, "polar",           0 ); // Activate polar-Cartesian coordinates
    model->ProgReac        = ReadInt2( fin, "ProgReac",        0 ); // Activate progressive reactions
    model->NoReturn        = ReadInt2( fin, "NoReturn",        0 ); // Turns off retrogression if 1.0
    model->UnsplitDiffReac = ReadInt2( fin, "UnsplitDiffReac", 0 ); // Unsplit diffusion reaction
    model->kinetics        = ReadInt2( fin, "kinetics",        0 ); // Unsplit diffusion reaction
    model->VolChangeReac   = ReadInt2( fin, "VolChangeReac",   0 ); // Turns on volume change due to reaction if 1
    model->Plith_trick     = ReadInt2( fin, "Plith_trick",     0 );
    model->DirectNeighbour = ReadInt2( fin, "DirectNeighbour", 0);
    model->Reseed          = ReadInt2( fin, "Reseed",          1); // Activates reseeding / particle injection
    model->ConservInterp   = ReadInt2( fin, "ConservInterp",   0); // Activates Taras conservative interpolation
    model->SmoothSoftening = ReadInt2( fin, "SmoothSoftening", 1); // Activates smooth explicit kinematic softening function
    model->oop             = ReadInt2( fin, "oop",             0); // Out-of-plane strain 
    model->noise_bg        = ReadInt2( fin, "noise_bg",        0); // Background noise generated on the particles --> mesh (used for cohesion)
    model->residual_form   = ReadInt2( fin, "residual_form",   1); // form of residual - TODO: delete if our models work with new default value (1)
    if ( model->shear_style == 1 ) model->isperiodic_x  = 1;
    if ( model->shear_style == 0 ) model->isperiodic_x  = 0;
    if ( model->aniso       == 1 ) model->fstrain       = 1;
    // Setup dependant
    model->EpsBG           = ReadDou2( fin, "EpsBG",         1e-30 ) / scaling->E; // Background tectonic rate, defaut is close to zero to avoid any Nans of Infs in rheology
    model->DivBG           = ReadDou2( fin, "DivBG",           0.0 ) / scaling->E;
    model->PrBG            = ReadDou2( fin, "PrBG",            0.0 ) / scaling->S;
    model->TBG             = ReadDou2( fin, "TBG",             0.0 ) / scaling->T;
    // Surface processes
    model->surf_diff       = ReadDou2( fin, "surf_diff",       0.0 ) / (pow(scaling->L,2.0)/scaling->t);
    model->surf_ised1      = ReadInt2( fin, "surf_ised1",      0.0 );
    model->surf_ised2      = ReadInt2( fin, "surf_ised2",      0.0 );
    model->surf_sedirate   = ReadDou2( fin, "surf_sedirate",   0.0 ) / scaling->V;
    model->surf_baselev    = ReadDou2( fin, "surf_baselev",    0.0 ) / scaling->L;
    model->surf_Winc       = ReadDou2( fin, "surf_Winc",       0.0 ) / scaling->L;
    model->surf_Vinc       = ReadDou2( fin, "surf_Vinc",       0.0 ) / scaling->V;
    // Initial thermal perturbation
    model->therm_pert_x0   = ReadDou2( fin, "therm_pert_x0",   0.0 ) / scaling->L;
    model->therm_pert_z0   = ReadDou2( fin, "therm_pert_z0",   0.0 ) / scaling->L;
    model->therm_pert_rad  = ReadDou2( fin, "therm_pert_rad",  0.0 ) / scaling->L;
    model->therm_pert_dT   = ReadDou2( fin, "therm_pert_dT" ,  0.0 ) / scaling->T;
    // For rheological database reasons...
    model->force_act_vol_ast = ReadInt2( fin, "force_act_vol_ast",   0 );
    model->act_vol_dis_ast   = ReadDou2( fin, "act_vol_dis_ast" ,  0.0 );
    model->act_vol_dif_ast   = ReadDou2( fin, "act_vol_dif_ast" ,  0.0 );
    // Model user's delights
    model->user0           = ReadDou2( fin, "user0",           0.0 );
    model->user1           = ReadDou2( fin, "user1",           0.0 );
    model->user2           = ReadDou2( fin, "user2",           0.0 );
    model->user3           = ReadDou2( fin, "user3",           0.0 );
    model->user4           = ReadDou2( fin, "user4",           0.0 );
    model->user5           = ReadDou2( fin, "user5",           0.0 );
    model->user6           = ReadDou2( fin, "user6",           0.0 );
    model->user7           = ReadDou2( fin, "user7",           0.0 );
    model->user8           = ReadDou2( fin, "user8",           0.0 );
    // Derived quantities
    model->dx                = (model->xmax - model->xmin) / (model->Nx - 1);
    model->dz                = (model->zmax - model->zmin) / (model->Nz - 1);
    model->dt0               = model->dt;
    model->dt_start          = model->dt;
    model->dt_max            = ReadDou2( fin, "dt_max",     1e20 ) / scaling->t; // maximum allowed time step, the default value is set to ~infinite, it we become effective only if specificaly set in XXX.txt (see e.g. LithoScale.txt)
    model->dt_min            = ReadDou2( fin, "dt_min",    -1e20 ) / scaling->t; // minimum allowed time step, defaut is negative such that it will never be activated unless specifically set in XXX.txt file
    model->eta_avg           = ReadInt2( fin, "eta_avg",       0 );              // 0: arithmetic mean - 1: harmonic mean - 2: geometric mean
    model->itp_stencil       = ReadInt2( fin, "itp_stencil",       1   );        // 1: 1-Cell          - 9: 9-Cell
    if (model->itp_stencil!=1 && model->itp_stencil!=9) { printf("Wrong value of itp_stencil: shoulbd be 1 or 9.\n"); exit(1); }
    model->nexp_radial_basis = ReadDou2( fin, "nexp_radial_basis", 1.0 ); // exponent for radial basis function interp. TODO: remove this unused option

       // For Cindy's setup
    model->diffuse_X       = ReadInt2( fin, "diffuse_X",     0 );              // 0 or 1
    model->diffuse_avg     = ReadInt2( fin, "diffuse_avg",   0 );              // 0: arithmetic mean - 1: harmonic mean - 2: geometric mean
    model->diffusion_length= ReadDou2( fin, "diffusion_length",  0.0 ) / scaling->L;

    // Gravity
    model->gx              = ReadDou2( fin, "gx",  0.0 ) / scaling->a;
    model->gz              = ReadDou2( fin, "gz",  0.0 ) / scaling->a;

    // Material properties
    model->Nb_phases = materials->Nb_phases =  ReadInt2( fin, "Nb_phases", 0 );
    for ( k=0; k<materials->Nb_phases; k++) {
        // Read general parameters
        materials->rho[k]  = ReadMatProps( fin, "rho", k,   2700.0 )  / scaling->rho;
        materials->mu[k]   = ReadMatProps( fin, "mu",  k,   1.0e10 )  / scaling->S;
        materials->Cv[k]   = ReadMatProps( fin, "Cv",  k,   1.0e3  )  / scaling->Cv;
        materials->k[k]    = ReadMatProps( fin, "k",   k,   1.0e-6 )  / scaling->k;
        materials->k_eff[k] = materials->k[k];
        materials->Qr[k]   = ReadMatProps( fin, "Qr",  k,   1.0e-30)  / (scaling->W / pow(scaling->L,3.0));
        materials->alp[k]  = ReadMatProps( fin, "alp", k,      0.0)  / (1.0/scaling->T);
        materials->bet[k]  = ReadMatProps( fin, "bet", k,  1.0e-40 )  / (1.0/scaling->S);
        materials->drho[k] = ReadMatProps( fin, "drho",k,      0.0 )  / (scaling->rho);
        materials->T0[k]   = (zeroC) / (scaling->T); // +20
        materials->P0[k]   = 1e5 / (scaling->S);
        // Read flow law parameters
        materials->cstv[k]  = ReadMatProps( fin, "cstv",k,    1.0  );
        materials->pwlv[k]  = ReadMatProps( fin, "pwlv",k,    0.0  );
        materials->linv[k]  = ReadMatProps( fin, "linv",k,    0.0  );
        materials->gbsv[k]  = ReadMatProps( fin, "gbsv",k,    0.0  );
        materials->expv[k]  = ReadMatProps( fin, "expv",k,    0.0  );
        gsel                = ReadMatProps( fin, "gsel",k,    0.0  );
        materials->eta0[k]  = ReadMatProps( fin, "eta0",k, 1.0e20  );
        materials->npwl[k]  = ReadMatProps( fin, "npwl",k,    1.0  );
        materials->Qpwl[k]  = ReadMatProps( fin, "Qpwl",k,    0.0  );
        materials->pref_pwl[k] = ReadMatProps( fin, "pref_pwl",k,    1.0 );    // weakening prefactor for power law
        materials->gs[k]    = ReadMatProps( fin, "gs",    k,    0.0   );
        materials->gs_ref[k]= ReadMatProps( fin, "gs_ref" ,k,  2.0e-3  ) /scaling->L;
        materials->kin[k]   = ReadMatProps( fin, "kin",    k,    0.0   );
        // Read generals plasticity parameters
        materials->plast[k]= ReadMatProps( fin, "plast",k,     0.0 );               // This indicates if plasticity is activated for a given phase
        materials->C[k]    = ReadMatProps( fin, "C",    k,   1.0e7 )  / scaling->S;
        materials->phi[k]  = ReadMatProps( fin, "phi",  k,    30.0 )  * M_PI/ 180.0;
        materials->psi[k]  = ReadMatProps( fin, "psi",  k,     0.0 )  * M_PI/ 180.0;
        materials->Slim[k] = ReadMatProps( fin, "Slim" ,k,  1.0e90 )  / scaling->S;
        // Strain softening
        materials->coh_soft[k]   = (int)ReadMatProps( fin, "coh_soft",   k,    0.0   );
        materials->phi_soft[k]   = (int)ReadMatProps( fin, "phi_soft",   k,    0.0   );
        materials->psi_soft[k]   = (int)ReadMatProps( fin, "psi_soft",   k,    0.0   );
        materials->is_tensile[k] = (int)ReadMatProps( fin, "is_tensile", k,    0.0   );
        materials->C_end[k]     = ReadMatProps( fin, "Ce",     k,    materials->C[k]*scaling->S    ) / scaling->S;
        materials->phi_end[k]   = ReadMatProps( fin, "phie",   k,    materials->phi[k]*180.0/M_PI  ) * M_PI / 180.0;
        materials->psi_end[k]   = ReadMatProps( fin, "psie",   k,    materials->psi[k]*180.0/M_PI  ) * M_PI / 180.0;
        double eps_coh = 1.0 / scaling->S;
        double eps_phi = 0.1  * M_PI / 180.0;
        double eps_psi = 0.1  * M_PI / 180.0;
        if ( materials->coh_soft[k] == 1 && fabs( materials->C_end[k]   - materials->C[k]  ) < eps_coh ) { printf("Please set a difference in cohesion, if not set coh_soft of phase %d to 0.0\n", k); exit(122); };
        if ( materials->phi_soft[k] == 1 && fabs( materials->phi_end[k] - materials->phi[k]) < eps_phi ) { printf("Please set a difference in friction angle, if not set phi_soft of phase %d to 0.0\n", k); exit(122); };
        if ( materials->psi_soft[k] == 1 && fabs( materials->psi_end[k] - materials->psi[k]) < eps_psi ) { printf("Please set a difference in dilation angle, if not set psi_soft of phase %d to 0.0\n", k); exit(122); };
        materials->pls_start[k] = ReadMatProps( fin, "plss",   k,    1.0    );
        materials->pls_end[k]   = ReadMatProps( fin, "plse",   k,    2.0    );
        // Reaction stuff
        materials->reac_soft[k]  = (int)ReadMatProps( fin, "reac_soft",   k,    0.0  );
        materials->reac_phase[k] = (int)ReadMatProps( fin, "reac_phase",   k,    0.0  );
        materials->Pr[k]         = ReadMatProps( fin, "Pr",      k,    0.0  ) / scaling->S;
        materials->dPr[k]        = ReadMatProps( fin, "dPr",     k,    0.0  ) / scaling->S;
        materials->tau_kin[k]    = ReadMatProps( fin, "tau_kin", k,    0.0  ) / scaling->t;
        materials->k_chem[k]     = ReadMatProps( fin, "k_chem",  k,    0.0  ) / (scaling->L*scaling->L/scaling->t);
        // Density models
        materials->density_model[k]     = (int)ReadMatProps( fin, "density_model",     k,    3  );
        materials->phase_diagram[k]     = (int)ReadMatProps( fin, "phase_diagram",     k,   -1  );
        // Viscoplasticity
        materials->n_vp[k]      = ReadMatProps( fin, "n_vp",   k,       1.0 ) ;
        materials->eta_vp[k]    = ReadMatProps( fin, "eta_vp", k,       0.0 ) / scaling->S / pow(scaling->t, 1.0/materials->n_vp[k]);
        if ( materials->is_tensile[k]==1 && materials->n_vp[k]>1.0001 ) {
            printf("Power-law visco-plasticity not yet compatible with tensile yielding\n");
            exit(1);
        }
        // Diffused rheological contrasts
        materials->phase_mix[k]  = (int)ReadMatProps( fin, "phase_mix",k,          0.0  );
        materials->phase_two[k]  = (int)ReadMatProps( fin, "phase_mix",k,    (double)k  );
        // Anisotropy
        materials->aniso_factor[k] =      ReadMatProps( fin, "aniso_factor", k,    1.0  );
        materials->aniso_angle[k]  =      ReadMatProps( fin, "aniso_angle"  ,k,   90.0  )  * M_PI/ 180.0;
        // Check if any flow law is active
        int sum = abs(materials->cstv[k]) + abs(materials->pwlv[k]) + abs(materials->linv[k]) + abs(materials->gbsv[k]) + abs(materials->expv[k]);
        if ( sum == 0 ) {
            printf ("Phase %0d has no determined flow mechanism\n Simulation will end now!\n", k);
            exit(12);
        }
        // Print material parameters
        printf("----------------------------------------- MODEL DOMAIN ------------------------------------------\n");
        printf("Xmin   = %2.1lf  km         Xmax   = %2.1lf  km     Nx   = %3d    dx   = %.2lf m\n", (model->xmin*scaling->L)/1e3, (model->xmax*scaling->L)/1e3, model->Nx, model->dx*scaling->L);
        printf("Zmin   = %2.1lf  km         Zmax   = %2.1lf  km      Nz   = %3d    dz   = %.2lf m\n", (model->zmin*scaling->L)/1e3, (model->zmax*scaling->L)/1e3, model->Nz, model->dz*scaling->L );
        printf("-------------------------------------------- PHASE: %d -------------------------------------------\n", k);
        printf("rho    = %2.2e kg/m^3     mu = %2.2e Pa\n", materials->rho[k]*scaling->rho, materials->mu[k]*scaling->S );
        printf("Cv     = %2.2e J/kg/K      k = %2.2e W/m/K      Qr = %2.2e W/m3\n", materials->Cv[k]*scaling->Cv, materials->k[k]*scaling->k, materials->Qr[k]*(scaling->W / pow(scaling->L,3)) );
        printf("C      = %2.2e Pa        phi = %2.2e deg      Slim = %2.2e Pa\n",  materials->C[k]*scaling->S, materials->phi[k]*180/M_PI, materials->Slim[k]*scaling->S );
        printf("alp    = %2.2e 1/T        T0 = %2.2e K         bet = %2.2e 1/Pa       P0 = %2.2e Pa       drho = %2.2e kg/m^3 \n", materials->alp[k]*(1/scaling->T), materials->T0[k]*(scaling->T), materials->bet[k]*(1/scaling->S), materials->P0[k]*(scaling->S), materials->drho[k]*scaling->rho );
        printf("prefactor for power-law: %2.2e\n", materials->pref_pwl[k]);
                 printf("C_end    = %2.2e Pa        Phi_end = %2.2e deg         pls_start = %2.2e        pls_end = %2.2e \n", materials->C_end[k]*scaling->S, materials->phi_end[k]*180/M_PI, materials->pls_start[k],  materials->pls_end[k] );
        printf("eta0_vp   = %2.2e  Pa.s^(1/n)         n_vp   = %2.2e\n", materials->eta_vp[k]* (scaling->S*pow(scaling->t,1.0/materials->n_vp[k])) , materials->n_vp[k]);

        printf("Flow law settings:\n");
        if ( abs(materials->cstv[k])>0 ) printf("--->    Constant viscosity activated \n");
        if ( abs(materials->pwlv[k])>0 ) printf("--->   Power law viscosity activated \n");
        if ( abs(materials->linv[k])>0 ) printf("--->      Linear viscosity activated \n");
        if ( abs(materials->gbsv[k])>0 ) printf("--->         GBS viscosity activated \n");
        if ( abs(materials->expv[k])>0 ) printf("---> Exponential viscosity activated \n");
        if ( abs(gsel)              >0 ) printf("--->  Grain size evolution activated \n");

        // Call flow law data base
        if ( abs(materials->pwlv[k])>0 ) ReadDataPowerLaw   ( materials, model, k, materials->pwlv[k], scaling );
        if ( abs(materials->linv[k])>0 ) ReadDataLinear     ( materials, model, k, materials->linv[k], scaling );
        if ( abs(materials->gbsv[k])>0 ) ReadDataGBS        ( materials, model, k, materials->gbsv[k], scaling );
        if ( abs(materials->expv[k])>0 ) ReadDataExponential( materials, model, k, materials->expv[k], scaling );
        if ( abs(materials->gs[k])  >0 ) ReadDataGSE        ( materials, model, k, materials->gs[k],   scaling );
        if ( abs(materials->kin[k]) >0 ) ReadDataKinetics   ( materials, model, k, materials->kin[k],  scaling );

        if ( abs(materials->cstv[k])>0 ) {
            materials->eta0[k]  /= scaling->eta;
            printf("eta0 = %2.2e Pa.s\n", materials->eta0[k]*scaling->eta);
        }
    }

    materials->R = Rg / (scaling->J/scaling->T);
    
    //------------------------------------------------------------------------------------------------------------------------------//
    // PHASE DIAGRAM INFO - simple pressure-dependent density model == 4 --- for quartz/coesite study
    //------------------------------------------------------------------------------------------------------------------------------//
    
    model->PD1DnP  = DoodzCalloc( model->Nb_phases, sizeof(int));
    model->PD1Drho = DoodzCalloc( model->Nb_phases, sizeof(double*));
    model->PD1Dmin = DoodzCalloc( model->Nb_phases, sizeof(double));
    model->PD1Dmax = DoodzCalloc( model->Nb_phases, sizeof(double));
    
    for ( k=0; k<materials->Nb_phases; k++) {
    
        if ( materials->density_model[k] == 4 ) {
     
            printf("Phase #%d has density_model %d, so:\n", k, materials->density_model[k]);
            printf("Loading 1D diagram for coesite quartz only baby...\n");
        
            model->PD1DnP[k]         = 10000;
            model->PD1Dmin[k]        = 1e5/scaling->S;
            model->PD1Dmax[k]        = 1e10/scaling->S;
            model->PD1Drho[k] = DoodzCalloc( model->PD1DnP[k], sizeof(double));
            char *fname = "PHASE_DIAGRAMS/QuartzCoesite600C.bin";
            if ( fopen(fname, "rb") != NULL ) {
                FILE* read = fopen(fname, "rb");
                fread ( model->PD1Drho[k],     sizeof(double), model->PD1DnP[k] , read);
                fclose(read);
                
                double p[model->PD1DnP[k]];
                double dP = ( model->PD1Dmax[k] - model->PD1Dmin[k]) / (model->PD1DnP[k] -1);
                p[0] =  model->PD1Dmin[k];
                
                // Dimensionalize
                for (int i=0; i<model->PD1DnP[k]; i++) {
                    model->PD1Drho[k][i] /= scaling->rho;
                    if (i>0) p[i] = p[i-1] + dP;
                }

                // Apply some diffusion...
                const double Kdiff   = 1e0;              // diffusivity
                const double dt_exp  =  dP*dP/Kdiff/2.1; // explicit time step
                const int    n_steps = 200;              // number of steps
                double qxW, qxE;   
                double *rho0   = DoodzCalloc( model->PD1DnP[k], sizeof(double));
                
                for (int it=0; it<n_steps; it++) {
                    ArrayEqualArray( rho0, model->PD1Drho[k], model->PD1DnP[k]);
                    for (int ix = 1; ix<model->PD1DnP[k]-1; ix++) {
                        qxW = - Kdiff*(rho0[ix] - rho0[ix-1])/dP;
                        qxE = - Kdiff*(rho0[ix+1] - rho0[ix])/dP;
                        model->PD1Drho[k][ix] = rho0[ix] - dt_exp*(qxE - qxW)/dP;
                    }
                }
                DoodzFree(rho0);
                
//                printf("scaling->rho = %2.10e\n",scaling->rho);
//                FILE* write = fopen("PHASE_DIAGRAMS/QuartzCoesite600C_smoothed.bin", "wb");
//                fwrite ( model->PD1Drho[k],     sizeof(double), model->PD1DnP[k] , write);
//                fclose(write);
                
//                //-------- In situ VISU for debugging: do not delete ------------------------//
//                FILE        *GNUplotPipe;
//                GNUplotPipe = popen ("gnuplot -persistent", "Work");
//                int NumCommands = 2;
//                char *GNUplotCommands[] = { "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5", "set pointintervalbox 3"};
//                for (int i=0; i<NumCommands; i++) fprintf(GNUplotPipe, "%s \n", GNUplotCommands[i]); //Send commands to gnuplot one by one.
//
//                fprintf(GNUplotPipe, "plot '-' with lines linestyle 1\n");
//                for (int i=0; i<model->PD1DnP[k]; i++) {
//                    fprintf(GNUplotPipe, "%lf %lf \n", p[i]*scaling->S, model->PD1Drho[k][i]*scaling->rho); //Write the data to a temporary file
//                }
//                fprintf(GNUplotPipe, "e\n");
//                fflush(GNUplotPipe);
//                 //-------- In situ VISU for debugging: do not delete ------------------------//
            }
            else {
                printf("Cannot open file %s, check if the file exists in the current location !\n Exiting", fname);
                exit(1);
            }
        }
    }
        
    //------------------------------------------------------------------------------------------------------------------------------//
    // PHASE DIAGRAM INFO
    //------------------------------------------------------------------------------------------------------------------------------//
    model->isPD = 0;

    for ( k=0; k<materials->Nb_phases; k++) {
        printf("Phase %d ---  density model %d \n",k, materials->density_model[k]);
        if ( materials->density_model[k] == 2 ) {
            printf("The density model of phase %0d relies on phase diagrams\n", k);
            model->isPD = 1;
        }
        if ( materials->density_model[k] == 2 && materials->phase_diagram[k] == -1 ) {
            printf("However the phase diagram index was not set (default -1)\n");
            printf("Therefore the simulation will not run, go fix 'phase_diagram' of phase %02d\n", k);
            exit(5);
        }
    }

    // Phase diagrams
    if ( model->isPD == 1 ) {

        printf("Loading phase diagrams...\n");
        int pid;
        model->num_PD = 10;

        // Allocate
        AllocatePhaseDiagrams( model );

        /**** PHASE DIAGRAMS #00 - Mantle (Jenadi_stx.dat)  ****/
        pid                       = 0;         // Kaus & Connolly, 2005: Effect of mineral phase transitions on sedimentary basin subsidence and uplift
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1000;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1000;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 473.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 2273.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 100e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 15e9 /scaling->S;         // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Hawaiian_Pyrolite_rho_bin.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #01 - Mantle (Jenadi_stx_HR.dat)  ****/
        pid                       = 1;         // Kaus & Connolly, 2005: Effect of mineral phase transitions on sedimentary basin subsidence and uplift
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1500;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1500;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 273.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 2273.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 1e5/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 25e9 /scaling->S;         // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Hawaiian_Pyrolite_HR_rho_bin.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #02 - Basalt (MORB_L.dat)  ****/
        pid                       = 2;  // Water saturated MORB - Bulk composition taken from Schmidt & Poli 1998 EPSL (Table 1)
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1000;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1000;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 573.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1273.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 100e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.1e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/MORB_H2Osat_rho_bin.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #03 - Andesite (Andesite.dat)  ****/
        pid                       = 3;  // Andesite
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1249;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1249;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 373.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1373.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 10.13e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.5816e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Andesite.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #04 - Hydrated Peridotite (Hydrated_Pdt.dat)  ****/
        pid                       = 4;  // Hydrated peridotite
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 793;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 793;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 373.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1373.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 10.13e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.5816e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Hydrated_Pdt.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #05 - MORB (MORB.dat)  ****/
        pid                       = 5;  // MORB
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1249;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1249;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 373.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1373.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 10.13e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.5816e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/MORB.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #06 - Pelite (Pelite.dat)  ****/
        pid                       = 6;  // Pelite
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1249;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1249;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 373.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1373.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 10.13e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.5816e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Pelite.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #07 - Rhyolite (Rhyolite.dat)  ****/
        pid                       = 7;  // Rhyolite
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 1249;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 1249;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 373.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1373.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 10.13e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.5816e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Rhyolite.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #08 - Serpentinite (Serpentinite.dat)  ****/
        pid                       = 8;  // Serpentinite
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 793;                     // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 793;                     // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = 373.0/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = 1373.0/scaling->T;        // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = 10.13e6/scaling->S;         // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = 5.5816e9 /scaling->S;        // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/Serpentinite.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);

        /**** PHASE DIAGRAMS #09 - Si02 (Si02.dat)  ****/
        pid                       = 9;  // Si02
        if (pid > (model->num_PD-1) ) {
            printf ("One should increment 'model->num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model->PDMnT[pid]         = 675;                           // Resolution for temperature (MANTLE) []
        model->PDMnP[pid]         = 2500;                          // Resolution for pressure    (MANTLE) []
        model->PDMTmin[pid]       = (398+1e-3)/scaling->T;         // Minimum temperature        (MANTLE) [K]
        model->PDMTmax[pid]       = (800+273)/scaling->T;          // Maximum temperature        (MANTLE) [K]
        model->PDMPmin[pid]       = (0.0090/10*1e9)/scaling->S;    // Minimum pressure           (MANTLE) [Pa]
        model->PDMPmax[pid]       = (49.9890/10*1e9)/scaling->S;   // Maximum pressure           (MANTLE) [Pa]
        model->PDMrho[pid]        = ReadBin( "PHASE_DIAGRAMS/SiO2_nsm010.dat", model->PDMnT[pid], model->PDMnP[pid], scaling->rho);
    }

    //------------------------------------------------------------------------------------------------------------------------------//
    // KINETIC DATA
    //------------------------------------------------------------------------------------------------------------------------------//

    if ( model->kinetics == 1 ) {

        printf("Loading kinetic data...\n");

        /**** Quartz Coesite (dG_QuartzCoesite.dat)  ****/
        model->kin_nT        = 675;                          // Resolution for temperature (MANTLE) []
        model->kin_nP        = 2500;                         // Resolution for pressure    (MANTLE) []
        model->kin_Tmin      = (398+1e-3)/scaling->T;        // Minimum temperature        (MANTLE) [K]
        model->kin_Tmax      = (800+273)/scaling->T;         // Maximum temperature        (MANTLE) [K]
        model->kin_Pmin      = (0.0090/10*1e9)/scaling->S;   // Minimum pressure           (MANTLE) [Pa]
        model->kin_Pmax      = (49.9890/10*1e9)/scaling->S;  // Maximum pressure           (MANTLE) [Pa]
        model->kin_dG        = ReadBin( "PHASE_DIAGRAMS/dG_QuartzCoesite.dat", model->kin_nT, model->kin_nP, scaling->J);
    }

    //------------------------------------------------------------------------------------------------------------------------------//
    // DEFORMATION MAP PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//

    // Create deformation maps or not (default no)
    model->def_maps        = ReadInt2(fin, "def_maps", 0);

    // Resolution
    model->nT              = ReadInt2(fin, "nT", 11);
    model->nE              = ReadInt2(fin, "nE", 11);
    model->nd              = ReadInt2(fin, "nd", 11);

    // Temperature, strain rate, grain size MIN/MAX & pressure
    model->Tmin            = ReadDou2(fin, "Tmin", 100.0);  model->Tmin += zeroC; model->Tmin /= scaling->T;    // C -> K & non-dimensionalization
    model->Tmax            = ReadDou2(fin, "Tmax", 1000.0); model->Tmax += zeroC; model->Tmax /= scaling->T;    // C -> K & non-dimensionalization
    model->Emin            = ReadDou2(fin, "Emin", -30.0);  //model->Emin /= scaling->E;                           // Non-dimensionalization
    model->Emax            = ReadDou2(fin, "Emax", -4.0 );   //model->Emax /= scaling->E;                           // Non-dimensionalization
    model->dmin            = ReadDou2(fin, "dmin", -7.0 );   //model->dmin /= scaling->L;                           // Non-dimensionalization
    model->dmax            = ReadDou2(fin, "dmax", -2.0 );   //model->dmax /= scaling->L;                           // Non-dimensionalization
    model->Pn              = ReadDou2(fin, "Pn",  5.0e8 );      model->Pn   /= scaling->S;                           // Non-dimensionalization

    //------------------------------------------------------------------------------------------------------------------------------//
    // NUMERICAL PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//

    // Particles
    particles->Nx_part       = ReadInt2( fin, "Nx_part", 4 );
    particles->Nz_part       = ReadInt2( fin, "Nz_part", 4 );
    particles->min_part_cell = ReadInt2( fin, "min_part_cell", 16 );
    particles->Nb_part       = (model->Nx-1)*(model->Nz-1) * particles->Nx_part * particles->Nz_part;
    particles->Nb_part_max   = 4.1*particles->Nb_part;

    // Nonlinear iteration parameters
    model->Newton           = ReadInt2( fin, "Newton", 0 );
    Nmodel->Picard2Newton   = ReadInt2( fin, "Picard2Newton", 0 );
    Nmodel->let_res_grow    = ReadInt2( fin, "let_res_grow",  0 );
    Nmodel->Pic2NewtCond    = ReadDou2( fin, "Pic2NewtCond", 1e-1 );
    Nmodel->nit_Pic_max     = ReadInt2( fin, "nit_Pic_max", 10 );
    if (model->Newton==0) Nmodel->Picard2Newton = 0;
    
    model->rel_tol_KSP      = ReadDou2( fin, "rel_tol_KSP", 1e-4 );
    Nmodel->nit_max         = ReadInt2( fin, "nit_max", 1 );
    Nmodel->abs_tol_u       = ReadDou2( fin, "abs_tol_u", 1.0e-6 );// / (scaling->S * scaling->L);                  // Fx * cel_vol
    Nmodel->abs_tol_p       = ReadDou2( fin, "abs_tol_p", 1.0e-6 );// / (scaling->E * scaling->L * scaling->L );    // Fp * cel_vol
    Nmodel->rel_tol_u       = ReadDou2( fin, "rel_tol_u", 1.0e-6 );
    Nmodel->rel_tol_p       = ReadDou2( fin, "rel_tol_p", 1.0e-6 );
    model->mineta           = ReadDou2( fin, "mineta", 1.0e18 ) / scaling->eta;
    model->maxeta           = ReadDou2( fin, "maxeta", 1.0e24 ) / scaling->eta;
    Nmodel->stagnated       = 0;

    // Direct solver parameters
    model->lsolver          = ReadInt2( fin, "lsolver", 2 );
    if ( model->lsolver == 0) {
        printf("WARNING!! Changing from solver type 0 to solver type 2!!! That's the new standard in MDOODZ 6.0.\n");
        model->lsolver = 2;
    }
    if ( model->Newton == 1 ) model->lsolver         = 2;

    // Close input file
    fclose(fin);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ScaleMe( scale* scale) {

    scale->t    = scale->L / scale->V;
    scale->a    = scale->V / scale->t;
    scale->E    = 1.0 / scale->t;
    scale->S    = scale->eta / scale->t;
    scale->m    = scale->S * scale->L * pow(scale->t,2.0);
    scale->rho  = scale->m / pow(scale->L,3.0);
    scale->F    = scale->m * scale->L / pow(scale->t, 2.0);
    scale->J    = scale->m * pow(scale->L,2.0) / pow(scale->t,2.0);
    scale->W    = scale->J / scale->t;
    scale->Cv   = scale->J / scale->m / scale->T;
    scale->k    = scale->W / scale->L / scale->T;
    scale->rhoE = scale->J / scale->m * scale->rho;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double* ReadBin( char A_name[], int nx, int ny, double scale ){
    double *A;
    char* bname; size_t nb_elems = nx*ny; FILE* fid; asprintf(&bname, "%s", A_name);

    A = malloc((nb_elems)*sizeof(double));
    fid=fopen(bname, "rb"); // Open file
    if (!fid){
        fprintf(stderr, "\nUnable to open file %s. I will exit here ... \n", bname); exit(2);
    }
    fread(A, sizeof(double), nb_elems, fid); fclose(fid); free(bname);
    ArrayTimesScalar( A, 1.0/scale, nb_elems );
    return A;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

char* ReadChar( FILE *fin, char FieldName[], char Default[] ) {
    int h, find=0, bufmax=50, length;
    char *string1, *string2;
    char line[bufmax];

    string1 = malloc(sizeof(char)*bufmax);

    // Get the length of the demanded string
    length = strlen( FieldName );

    // Buffer array to contain the string to compare with
    char *param1, *param2;
    param2 = malloc( (length+1)*sizeof(char));
    asprintf(&param1, "%s", FieldName);

    // Initialise line
    for (h=0;h<bufmax;h++) {
        line[h]='\0';
    }

    // Seach for sign 'equal' in the lines
    while ( find == 0 ) {

        // Read new line
        fgets ( line, sizeof(line), fin );

        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file, running with default value %s\n", FieldName, Default);
            rewind (fin);
            int str_size = strlen(Default);
            int h1;
            string2 = malloc((str_size+1)*sizeof(char));
            for (h1=0; h1<str_size; h1++) {
                string2[h1] = Default[h1];
            }
            free(string1);
            free(param1);
            free(param2);
            return string2;
        }

        // Get the first 'length' characters of the line
        for (h=0;h<length;h++) {
            param2[h] = line[h];
        }
        param2[length] = '\0';

        // Check if we found the right parameter name
        if ( strcmp(param1, param2) == 0 ) {
            find = 1;
            int str_size = 0;
            int h1;
            // Search for equal sign
            for (h=0;h<bufmax;h++) {
                if(strlen(line)> 0 && line[h]=='=') {

                    for (h1=0; h1<30; h1++) {
                        if ( isspace(line[h+2+h1]) ) {
                            string1[h1] = '\0';
                            str_size++;
                            break;
                        }
                        else {
                            string1[h1] = (line[h+2+h1]);
                            str_size++;
                        }
                    }

                    string2 = malloc((str_size+1)*sizeof(char));
                    for (h1=0; h1<str_size+1; h1++) {
                        string2[h1] = string1[h1];
                    }

                    free(param1);
                    free(param2);
                    free(string1);
                    return string2;
                }
            }
        }
    }
    free(param1);
    free(param2);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

char* ReadPhaseDiagram( FILE *fin, char FieldName[] ) {
    int h, find=0, bufmax=50, length;
    char *string1, *string2;
    char line[bufmax];

    string1 = malloc(sizeof(char)*bufmax);

    // Get the length of the demanded string
    length = strlen( FieldName );

    // Buffer array to contain the string to compare with
    char *param1, *param2;
    param2 = malloc( (length+1)*sizeof(char));
    asprintf(&param1, "%s", FieldName);

    // Initialise line
    for (h=0;h<bufmax;h++) {
        line[h]='\0';
    }

    // Seach for sign 'equal' in the lines
    while ( find == 0 ) {

        // Read new line
        fgets ( line, sizeof(line), fin );

        if (feof(fin)) {
            printf("Error: The phase diagram '%s' could not be found in the setup file. I will exit here.\n", FieldName);
            rewind (fin    );
            free   (param1 );
            free   (param2 );
            free   (string1);
            exit(2);
        }

        // Get the first 'length' characters of the line
        for (h=0;h<length;h++) {
            param2[h] = line[h];
        }
        param2[length] = '\0';

        // Check if we found the right parameter name
        if ( strcmp(param1, param2) == 0 ) {
            find = 1;
            int str_size = 0;
            int h1;
            // Search for equal sign
            for (h=0;h<bufmax;h++) {
                if(strlen(line)> 0 && line[h]=='=') {

                    for (h1=0; h1<30; h1++) {
                        if ( isspace(line[h+2+h1]) ) {
                            //printf("found space last char is %c\n", line[h+2+h1-1]);
                            string1[h1] = '\0';
                            str_size++;
                            break;
                        }
                        else {
                            string1[h1] = (line[h+2+h1]);
                            str_size++;
                        }
                    }

                    string2 = malloc((str_size+1)*sizeof(char));
                    for (h1=0; h1<str_size+1; h1++) {
                        string2[h1] = string1[h1];
                    }
                    //printf("returned string is %s %s %d\n", string1, string2, str_size );
                    free(param1);
                    free(param2);
                    free(string1);
                    return string2;
                }
            }
        }
    }
    free(param1);
    free(param2);
    return NULL;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int ReadInt2( FILE *fin, char FieldName[], int Default )
{
    // Some declaration.
    int     bufmax=1000;
    int     h = 0;
    char    line[bufmax];
    int     value = 0;

    // Start from beginning of the file.
    rewind(fin);

    // Buffer array to contain the string to compare with
    char *param1;
    asprintf(&param1, "%s", FieldName);

    // Loop over all lines in input file
    while(value == 0)
    {
        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file, running with default value %d\n", FieldName, Default);
            rewind (fin);
            free(param1);
            return Default;
        }

        // Determine the length of the parameter string in the current line.
        int InPar_length = 0;

        while(line[InPar_length] != ' ')
        {
            InPar_length = InPar_length + 1;
        }

        // Allocate memory to save the current parameter string.
        char *param2;
        param2 = malloc( (InPar_length + 1)*sizeof(char));

        // Save the current parameter string.
        for (h = 0; h < InPar_length; h++)
        {
            param2[h] = line[h];
        }
        param2[InPar_length] = '\0';

        // Find match.
        if( (strcmp(param1, param2)) == 0 )
        {
            // Search for equal sign
            for (h = 0; h < bufmax ;h++)
            {
                if(strlen(line)> 0 && line[h]=='=')
                {
                    value = atoi(&line[h+1]);
                    free(param1);
                    free(param2);
                    rewind(fin);
                    return value;
                }
            }
        }
        free(param2);
    }
    rewind(fin);
    free(param1);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ReadDou2( FILE *fin, char FieldName[], double Default )
{
    // Some declaration.
    int     bufmax=1000;
    int     h = 0;
    char    line[bufmax];
    double  value = 0.0;

    // Start from beginning of the file.
    rewind(fin);

    // Buffer array to contain the string to compare with
    char *param1;
    asprintf(&param1, "%s", FieldName);

    // Loop over all lines in input file
    while(value == 0)
    {
        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file, running with default value %2.2e\n", FieldName, Default);
            rewind (fin);
            free(param1);
            return Default;
        }

        // Determine the length of the parameter string in the current line.
        int InPar_length = 0;
        while(line[InPar_length] != ' ')
        {
            InPar_length = InPar_length + 1;
        }

        // Allocate memory to save the current parameter string.
        char *param2;
        param2 = malloc( (InPar_length + 1)*sizeof(char));

        // Save the current parameter string.
        for (h = 0; h < InPar_length; h++)
        {
            param2[h] = line[h];
        }
        param2[InPar_length] = '\0';

        // Find match.
        if( (strcmp(param1, param2)) == 0 )
        {
            // Search for equal sign
            for (h = 0; h < bufmax ;h++)
            {
                if(strlen(line)> 0 && line[h]=='=')
                {
                    value = atof(&line[h+1]);
                    free(param1);
                    free(param2);
                    rewind(fin);
                    return value;
                }
            }
        }
        free(param2);
    }
    rewind(fin);
    free(param1);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double ReadMatProps( FILE *fin, char FieldName[], int PhaseID, double Default )
{
    // Some declarations.
    int     bufmax = 1000;
    int     h = 0, ID_value, find = 0, subfind = 0;
    char    line[bufmax], phase_line[bufmax];
    double  value = 0.0;

    // Find the current phase ID starting from the beginning of the input file.
    rewind(fin);

    // Buffer array to contain the string to compare with
    char *param1;
    asprintf(&param1, "%s", FieldName);

    // Loop over all lines in input file.
    while(find == 0)
    {
        // Read new line and check for existence of phases. If not exit.
        fgets ( line, sizeof(line), fin );
        if (feof(fin))
        {
            printf("Warning : No phase ID found! I will exit here.");
            free(param1);
            rewind (fin);
            exit(0);
        }

        // Determine the length of the parameter string in the current line.
        int InPar_length = 0;

        while(line[InPar_length] != ' ')
        {
            InPar_length = InPar_length + 1;
        }

        // Allocate memory to save the current parameter string.
        char *param2;
        param2 = malloc( (InPar_length + 1)*sizeof(char) );

        // Save the current parameter string.
        for (h = 0; h < InPar_length; h++)
        {
            param2[h] = line[h];
        }
        param2[InPar_length] = '\0';

        // Find match.
        if( (strcmp("ID", param2)) == 0 )
        {
            // Search for equal sign
            for (h = 0; h < bufmax ;h++)
            {
                if(strlen(line)> 0 && line[h] == '=')
                {
                    ID_value = atoi(&line[h+1]);
                    break;
                }
            }

            // Read parameters if current phase.
            if (ID_value == PhaseID)
            {
                // Loop over all lines in input file
                while(subfind == 0)
                {
                    // Read new line
                    fgets ( phase_line, sizeof(phase_line), fin );

                    // Determine the length of the parameter string in the current line.
                    InPar_length = 0;
                    while(phase_line[InPar_length] != ' ')
                    {
                        InPar_length = InPar_length + 1;
                    }

                    // Allocate memory to save the current parameter string.
                    char *param3;
                    param3 = malloc( (InPar_length + 1)*sizeof(char));

                    // Save the current parameter string.
                    for (h = 0; h < InPar_length; h++)
                    {
                        param3[h] = phase_line[h];
                    }
                    param3[InPar_length] = '\0';

                    // Break in case the parameter has not been defined for the current phase.
                    if ( strcmp(param3,"ID") == 0 || feof(fin) ) {
                        if ( fabs(Default) <  100 ) printf("Warning : Parameter '%s' not found in the setup file, running with default value %.2lf\n", FieldName, Default);
                        if ( fabs(Default) >= 100 ) printf("Warning : Parameter '%s' not found in the setup file, running with default value %2.2e\n", FieldName, Default);
                        rewind (fin);
                        free(param1);
                        free(param2);
                        free(param3);
                        return Default;
                    }

                    // Find match.
                    if( (strcmp(param1, param3)) == 0 )
                    {
                        // Search for equal sign
                        for (h = 0; h < bufmax ;h++)
                        {

                            if(strlen(line)> 0 && phase_line[h]=='=')
                            {
                                value = atof(&phase_line[h+1]);
                                free(param1);
                                free(param2);
                                free(param3);
                                return value;
                            }
                        }
                        subfind = 1;
                    }
                    free(param3);
                }
                free(param1);
            }
        }
        free(param2);
    }
    free(param1);
    return Default;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateInputFile( char fin_name[], int NewNumber ) {

    FILE* fin;
    char  FieldName[] = "istep";
    int h, find=0, bufmax=50, length;
    char line[bufmax];

    // Get the length of the demanded string
    length = strlen( FieldName );

    // Buffer array to contain the string to compare with
    char *param1, *param2;
    param2 = malloc( (length+1)*sizeof(char));
    asprintf(&param1, "%s", FieldName);

    // Open File
    fin = fopen(fin_name,"r+");
	if ( fin == NULL) {
		printf("Setup file '%s' does not exist\nExiting...\n", fin_name);
        fclose(fin);
		exit(1);
	}
    // Seach for sign 'equal' in the lines
    while ( find == 0 && feof(fin) !=1 ) {

        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            printf("Warning : Parameter '%s' not found in the setup file\n", FieldName);
        }

        // Get the first 'length' characters of the line
        for (h=0;h<length;h++) {
            param2[h] = line[h];
        }
        param2[length] = '\0';

        // printf("p1 %s p2 %s\n", param1, param2);

        // Check if we found the right parameter name
        if ( strcmp(param1, param2) == 0 ) {
            find = 1;
            // Search for equal sign and replace by the current step
            for (h=0;h<bufmax;h++) {

                //printf(" %c ", line[h]);

                if(strlen(line)> 0 && line[h]=='=') {

                    //printf("pos =%d SEEK_CUR=%d NewNumber=%d\n", pos, SEEK_CUR, NewNumber);
                    fseek(fin, -6, SEEK_CUR);
                    fprintf( fin, "%05d", NewNumber);
                    break;
                }
            }
        }
    }

    // Close file
    fclose(fin);
    free(param1);
    free(param2);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
