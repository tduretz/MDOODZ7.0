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
#include "mdoodz-private.h"

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

char* GetSetupFileName(int nargs, char *args[]) {
  char *fileExtension = ".txt";
  if (nargs == 2 && strstr(args[1], fileExtension)) {
    return args[1];
  } else {
    return strcat(args[0], fileExtension);
  }
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

    if (model->free_surface == 1) {
        fread(&topo_chain->Nb_part,   s1,                   1, file );
        fread( topo_chain->x,         s3, topo_chain->Nb_part, file );
        fread( topo_chain->z,         s3, topo_chain->Nb_part, file );
    }

    int Nb_part0 = particles->Nb_part;
    fread( &particles->Nb_part,  s1, 1, file);
    fread( particles->x,    s3, particles->Nb_part, file);
    fread( particles->z,    s3, particles->Nb_part, file);
    fread( particles->phase, s1, particles->Nb_part, file);
    fread( particles->dual, s1, particles->Nb_part, file);

    fclose(file);
    free(name);
    //---------------------------------------------------------------------------------------------------------//

    if (model->free_surface == 1) {
        topo_chain_ini->Nb_part = topo_chain->Nb_part;
        // note: the topo_chain should be scaled as well
#pragma omp parallel for shared( topo_chain, model, scaling )
        for ( k=0; k<topo_chain->Nb_part; k++ ) {
            topo_chain->x[k]     /=scaling.L;
            topo_chain->z[k]     /=scaling.L;
            topo_chain_ini->x[k] = topo_chain->x[k];
            topo_chain_ini->z[k] = topo_chain->z[k];
        }
    }
#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /=scaling.L;
        particles->z[k]     /=scaling.L;
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


    if (model->free_surface == 1) {
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
    fread( particles->rho,  s3, particles->Nb_part, file);
    fread( particles->Vx,   s3, particles->Nb_part, file);
    fread( particles->Vz,   s3, particles->Nb_part, file);
    fread( particles->phi,  s3, particles->Nb_part, file);
    fread( particles->X  ,  s3, particles->Nb_part, file);
    fread( particles->noise, s3, particles->Nb_part, file);
    fread( particles->phase, s1, particles->Nb_part, file);
    fread( particles->dual, s1, particles->Nb_part, file);

    if (model->elastic == 1) {
        fread( particles->sxxd,   s3, particles->Nb_part, file );
        fread( particles->szzd,   s3, particles->Nb_part, file );
        fread( particles->sxz,    s3, particles->Nb_part, file );
        // fread( particles->dsxxd,   s3, particles->Nb_part, file );
        // fread( particles->dszzd,   s3, particles->Nb_part, file );
        // fread( particles->dsxz,    s3, particles->Nb_part, file );

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

    if (model->finite_strain == 1) {
        fread( particles->Fxx         , s3, particles->Nb_part, file);
        fread( particles->Fxz         , s3, particles->Nb_part, file);
        fread( particles->Fzx         , s3, particles->Nb_part, file);
        fread( particles->Fzz         , s3, particles->Nb_part, file);
    }

    if (model->anisotropy == 1) {
        fread( particles->nx         , s3, particles->Nb_part, file);
        fread( particles->nz         , s3, particles->Nb_part, file);
    }

    if (model->track_T_P_x_z == 1) {
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
        if (model->free_surface == 1) {
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

        // Phase proportions
        printf("Loading phase proportions - Nb_phases = %d:\n", model->Nb_phases);
        for ( k=0; k<  model->Nb_phases; k++ ) { //model->Nb_phases
            fread( mesh->phase_perc_n[k], s3,  Ncx*Ncz,   file );
            fread( mesh->phase_perc_s[k], s3,  Nx *Nz ,   file );
        }
        fread( mesh->nb_part_cell, s1,  Ncx*Ncz,   file );
        fread( mesh->nb_part_vert, s1,  Nx *Nz ,   file );

        // Timeseries
        fread( mesh->Time_time      , s3,   model->Nt+1 ,   file );
        fread( mesh->Short_time     , s3,   model->Nt+1 ,   file );
        fread( mesh->Work_time      , s3,   model->Nt+1 ,   file );
        fread( mesh->Uthermal_time  , s3,   model->Nt+1 ,   file );
        fread( mesh->Uelastic_time  , s3,   model->Nt+1 ,   file );
        fread( mesh->T_mean_time    , s3,   model->Nt+1 ,   file );
        fread( mesh->P_mean_time    , s3,   model->Nt+1 ,   file );
        fread( mesh->sxxd_mean_time , s3,   model->Nt+1 ,   file );
        fread( mesh->szzd_mean_time , s3,   model->Nt+1 ,   file );
        fread( mesh->sxz_mean_time  , s3,   model->Nt+1 ,   file );
        fread( mesh->Tii_mean_time  , s3,   model->Nt+1 ,   file );
        fread( mesh->exxd_mean_time , s3,   model->Nt+1 ,   file );
        fread( mesh->ezzd_mean_time , s3,   model->Nt+1 ,   file );
        fread( mesh->exz_mean_time  , s3,   model->Nt+1 ,   file );
        fread( mesh->Eii_mean_time  , s3,   model->Nt+1 ,   file );


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

    for ( k=0; k<model->Nt+1; k++ ) {
        mesh->Time_time[k]      /= scaling.t;
        mesh->Short_time[k]     /= 1.0;
        mesh->Work_time[k]      /= (scaling.rhoE*scaling.L*scaling.L);
        mesh->Uthermal_time[k]  /= (scaling.rhoE*scaling.L*scaling.L);
        mesh->Uelastic_time[k]  /= (scaling.rhoE*scaling.L*scaling.L);
        mesh->T_mean_time[k]    /= scaling.T;
        mesh->P_mean_time[k]    /= scaling.S;
        mesh->sxxd_mean_time[k] /= scaling.S;
        mesh->szzd_mean_time[k] /= scaling.S;
        mesh->sxz_mean_time[k]  /= scaling.S;
        mesh->Tii_mean_time[k]  /= scaling.S;
        mesh->exxd_mean_time[k] /= scaling.E;
        mesh->ezzd_mean_time[k] /= scaling.E;
        mesh->exz_mean_time[k]  /= scaling.E;
        mesh->Eii_mean_time[k]  /= scaling.E;
    }

#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /=scaling.L;
        particles->z[k]     /=scaling.L;
        particles->P[k]     /=scaling.S;
        particles->rho[k]   /=scaling.rho;
        particles->Vx[k]    /=scaling.V;
        particles->Vz[k]    /=scaling.V;
        particles->phi[k]   /=1.0;
        particles->X[k]     /=1.0;
        if ( model->elastic == 1 ) {
            particles->sxxd[k]    /= scaling.S;
            particles->szzd[k]    /= scaling.S;
            particles->sxz[k]     /= scaling.S;
            // particles->dsxxd[k]    /= scaling.S;
            // particles->dszzd[k]    /= scaling.S;
            // particles->dsxz[k]     /= scaling.S;
        }

        particles->d[k]      /= scaling.L;
        particles->divth[k]  /= scaling.E;
        particles->T[k]      /= scaling.T;

        if ( model->track_T_P_x_z == 1) {
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

    if (model->free_surface == 1) {
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

    for ( k=0; k<model.Nt+1; k++ ) {
        mesh->Time_time[k]      *= scaling.t;
        mesh->Short_time[k]     *= 1.0;
        mesh->Work_time[k]      *= (scaling.rhoE*scaling.L*scaling.L);
        mesh->Uthermal_time[k]  *= (scaling.rhoE*scaling.L*scaling.L);
        mesh->Uelastic_time[k]  *= (scaling.rhoE*scaling.L*scaling.L);
        mesh->T_mean_time[k]    *= scaling.T;
        mesh->P_mean_time[k]    *= scaling.S;
        mesh->sxxd_mean_time[k] *= scaling.S;
        mesh->szzd_mean_time[k] *= scaling.S;
        mesh->sxz_mean_time[k]  *= scaling.S;
        mesh->Tii_mean_time[k]  *= scaling.S;
        mesh->exxd_mean_time[k] *= scaling.E;
        mesh->ezzd_mean_time[k] *= scaling.E;
        mesh->exz_mean_time[k]  *= scaling.E;
        mesh->Eii_mean_time[k]  *= scaling.E;
    }

#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]      *= scaling.L;
        particles->z[k]      *= scaling.L;
        particles->P[k]      *= scaling.S;
        particles->rho[k]    *= scaling.rho;
        particles->Vx[k]     *= scaling.V;
        particles->Vz[k]     *= scaling.V;
        particles->phi[k]    *= 1.0;
        particles->X[k]      *= 1.0;
        if (model.elastic == 1) {
            particles->sxxd[k] *= scaling.S;
            particles->szzd[k] *= scaling.S;
            particles->sxz[k]  *= scaling.S;
            // particles->dsxxd[k] *= scaling.S;
            // particles->dszzd[k] *= scaling.S;
            // particles->dsxz[k]  *= scaling.S;
        }
        particles->divth[k]  *= scaling.E;
        particles->T[k]         *= scaling.T;
        particles->d[k]        *= scaling.L;

        if ( model.track_T_P_x_z == 1) {
            particles->T0[k]  *= scaling.T;
            particles->P0[k]  *= scaling.S;
            particles->x0[k]  *= scaling.L;
            particles->z0[k]  *= scaling.L;
            particles->Tmax[k]*= scaling.T;
            particles->Pmax[k]*= scaling.S;
        }
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

    if (model.free_surface == 1) {
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

    if (model.free_surface == 1) {
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
    fwrite( particles->rho,   s3, particles->Nb_part, file);
    fwrite( particles->Vx,    s3, particles->Nb_part, file);
    fwrite( particles->Vz,    s3, particles->Nb_part, file);
    fwrite( particles->phi,   s3, particles->Nb_part, file);
    fwrite( particles->X  ,   s3, particles->Nb_part, file);
    fwrite( particles->noise, s3, particles->Nb_part, file);
    fwrite( particles->phase, s1, particles->Nb_part, file);
    fwrite( particles->dual,  s1, particles->Nb_part, file);

    if (model.elastic == 1) {
        fwrite( particles->sxxd,   s3, particles->Nb_part, file );
        fwrite( particles->szzd,   s3, particles->Nb_part, file );
        fwrite( particles->sxz,    s3, particles->Nb_part, file );
        // fwrite( particles->dsxxd,   s3, particles->Nb_part, file );
        // fwrite( particles->dszzd,   s3, particles->Nb_part, file );
        // fwrite( particles->dsxz,    s3, particles->Nb_part, file );

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

    if (model.finite_strain == 1) {
        fwrite( particles->Fxx         , s3, particles->Nb_part, file);
        fwrite( particles->Fxz         , s3, particles->Nb_part, file);
        fwrite( particles->Fzx         , s3, particles->Nb_part, file);
        fwrite( particles->Fzz         , s3, particles->Nb_part, file);
    }

    if (model.anisotropy == 1) {
        fwrite( particles->nx         , s3, particles->Nb_part, file);
        fwrite( particles->nz         , s3, particles->Nb_part, file);
    }

    if (model.track_T_P_x_z == 1) {
        fwrite( particles->T0         , s3, particles->Nb_part, file);
        fwrite( particles->P0         , s3, particles->Nb_part, file);
        fwrite( particles->x0         , s3, particles->Nb_part, file);
        fwrite( particles->z0         , s3, particles->Nb_part, file);
        fwrite( particles->Tmax       , s3, particles->Nb_part, file);
        fwrite( particles->Pmax       , s3, particles->Nb_part, file);
    }

    fwrite( particles->divth, s3, particles->Nb_part, file);
    fwrite( particles->T,     s3, particles->Nb_part, file);
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
    if (model.free_surface == 1) {
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

    // Phase proportions
    printf("Loading phase proportions - Nb_phases = %d:\n", model.Nb_phases);
//    fwrite( &model.Nb_phases, s1,        1,   file );
    for ( k=0; k< model.Nb_phases; k++ ) {
        fwrite( mesh->phase_perc_n[k], s3,  Ncx*Ncz,   file );
        fwrite( mesh->phase_perc_s[k], s3,  Nx *Nz ,   file );
    }
    fwrite( mesh->nb_part_cell, s1,  Ncx*Ncz,   file );
    fwrite( mesh->nb_part_vert, s1,  Nx *Nz ,   file );

    // Timeseries
    fwrite( mesh->Time_time      , s3,   model.Nt+1 ,   file );
    fwrite( mesh->Short_time     , s3,   model.Nt+1 ,   file );
    fwrite( mesh->Work_time      , s3,   model.Nt+1 ,   file );
    fwrite( mesh->Uthermal_time  , s3,   model.Nt+1 ,   file );
    fwrite( mesh->Uelastic_time  , s3,   model.Nt+1 ,   file );
    fwrite( mesh->T_mean_time    , s3,   model.Nt+1 ,   file );
    fwrite( mesh->P_mean_time    , s3,   model.Nt+1 ,   file );
    fwrite( mesh->sxxd_mean_time , s3,   model.Nt+1 ,   file );
    fwrite( mesh->szzd_mean_time , s3,   model.Nt+1 ,   file );
    fwrite( mesh->sxz_mean_time  , s3,   model.Nt+1 ,   file );
    fwrite( mesh->Tii_mean_time  , s3,   model.Nt+1 ,   file );
    fwrite( mesh->exxd_mean_time , s3,   model.Nt+1 ,   file );
    fwrite( mesh->ezzd_mean_time , s3,   model.Nt+1 ,   file );
    fwrite( mesh->exz_mean_time  , s3,   model.Nt+1 ,   file );
    fwrite( mesh->Eii_mean_time  , s3,   model.Nt+1 ,   file );

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

    for ( k=0; k<model.Nt+1; k++ ) {
        mesh->Time_time[k]      /= scaling.t;
        mesh->Short_time[k]     /= 1.0;
        mesh->Work_time[k]      /= (scaling.rhoE*scaling.L*scaling.L);
        mesh->Uthermal_time[k]  /= (scaling.rhoE*scaling.L*scaling.L);
        mesh->Uelastic_time[k]  /= (scaling.rhoE*scaling.L*scaling.L);
        mesh->T_mean_time[k]    /= scaling.T;
        mesh->P_mean_time[k]    /= scaling.S;
        mesh->sxxd_mean_time[k] /= scaling.S;
        mesh->szzd_mean_time[k] /= scaling.S;
        mesh->sxz_mean_time[k]  /= scaling.S;
        mesh->Tii_mean_time[k]  /= scaling.S;
        mesh->exxd_mean_time[k] /= scaling.E;
        mesh->ezzd_mean_time[k] /= scaling.E;
        mesh->exz_mean_time[k]  /= scaling.E;
        mesh->Eii_mean_time[k]  /= scaling.E;
    }

#pragma omp parallel for shared( particles, model, scaling )
    for ( k=0; k<particles->Nb_part; k++ ) {
        particles->x[k]     /= scaling.L;
        particles->z[k]     /= scaling.L;
        particles->P[k]     /= scaling.S;
        particles->rho[k]   /= scaling.rho;
        particles->Vx[k]    /= scaling.V;
        particles->Vz[k]    /= scaling.V;
        particles->phi[k]   /= 1.0;
        particles->X[k]     /= 1.0;

        if (model.elastic == 1) {
            particles->sxxd[k]   /= scaling.S;
            particles->szzd[k]   /= scaling.S;
            particles->sxz[k]    /= scaling.S;
            // particles->dsxxd[k]   /= scaling.S;
            // particles->dszzd[k]   /= scaling.S;
            // particles->dsxz[k]    /= scaling.S;
        }
        particles->divth[k]  /= scaling.E;
        particles->T[k]      /= scaling.T;
        particles->d[k]  /= scaling.L;

        if (model.track_T_P_x_z == 1) {
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

    if (model.free_surface == 1) {
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

Input ReadInputFile( char *fileName ) {

    //------------------------------------------------------------------------------------------------------------------------------//
    // MODEL PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//

    params model;
    mat_prop materials;
    Nparams Nmodel;
    FILE *fin;
    int k, gsel;

    fin = fopen(fileName,"rt");
    if (fin == NULL) {
        printf("Setup file '%s' does not exist\nExiting...\n", fileName);
        fclose(fin);
        exit(1);
    }

    // Simulation start/restart from Breakpoint
    model.istep              = ReadInt2( fin, "istep", 0 );
    model.irestart           = ReadInt2( fin, "irestart", 0 );
    if ( model.istep == 0 ) model.irestart = 0; // Override the restart step number written in the text file (istep)

    // Output
    model.writer             = ReadInt2( fin, "writer",              0 ); // Write files
    model.writer_step        = ReadInt2( fin, "writer_step",         1 ); // Frequency of output
    model.writer_markers     = ReadInt2( fin, "writer_markers",      0 ); // Writes marker files
    model.writer_debug       = ReadInt2( fin, "writer_debug",        0 ); // Writes debug files
    model.writer_subfolder   = ReadChar( fin, "writer_subfolder",  "./"); // Writes output in given subfolder
    // printf("%s\n", model.writer_subfolder); exit(1);
    model.noisy              = ReadInt2( fin, "noisy",               1 ); // Prints a lot of info to standard output
    model.track_T_P_x_z      = ReadInt2( fin, "track_T_P_x_z",       0 ); // Tracks initial T, P, x and z on particles 
    model.delete_breakpoints = ReadInt2( fin, "delete_breakpoints",  1 ); // Progressively deletes breakpoint files
    model.gnuplot_log_res    = ReadInt2( fin, "gnuplot_log_res",     0 ); // Activates GNU plot residuals visualisation

    // Input
    model.import_files_dir     = ReadChar( fin, "import_files_dir",    "../../IMPORT");
    model.import_file          = ReadChar( fin, "import_file",             "blah.bin");
    model.save_initial_markers = ReadInt2( fin, "save_initial_markers",            0 );
    model.load_initial_markers = ReadInt2( fin, "load_initial_markers",            0 );
    model.initial_markers_file = ReadChar( fin, "initial_markers_file", "markers.bin");

    // Read scales for non-dimensionalisation
    scale scaling            = (scale){
                        .eta = ReadDou2(fin, "eta", 1.0),
                        .L   = ReadDou2(fin, "L", 1.0),
                        .V   = ReadDou2(fin, "V", 1.0),
                        .T   = ReadDou2(fin, "T", 1.0),
    };
    ScaleMe( &scaling );
    double Ga = 1e9*365.25*3600*24;

    // Spatial domain
    model.Nx                 = ReadInt2( fin, "Nx",        10 );            // Number of vertices in x direction
    model.Nz                 = ReadInt2( fin, "Nz",        10 );            // Number of vertices in y direction
    model.xmin               = ReadDou2( fin, "xmin",    -1.0 )/scaling.L;  // Spatial domain extent
    model.xmax               = ReadDou2( fin, "xmax",     1.0 )/scaling.L;  // Spatial domain extent
    model.zmin               = ReadDou2( fin, "zmin",    -1.0 )/scaling.L;  // Spatial domain extent  
    model.zmax               = ReadDou2( fin, "zmax",     1.0 )/scaling.L;  // Spatial domain extent
    model.balance_boundaries = ReadInt2( fin, "balance_boundaries",   0 );  // Switch this to activate boundary velocity balancing
    model.zero_mean_topo     = ReadInt2( fin, "zero_mean_topo",       0 );  // Switch this to force zero mean topography
    // Time domain
    model.Nt                 = ReadInt2( fin, "Nt",         1 );            // Number of time steps    
    model.dt                 = ReadDou2( fin, "dt",       0.0 ) /scaling.t; // Time step
    model.Courant            = ReadDou2( fin, "Courant",             0.5 ); // Courant number
    model.RK                 = ReadInt2( fin, "RK",                    4 ); // Order of Runge-Kutta advection solver (1, 2 or 4)
    model.constant_dt        = ReadInt2( fin, "constant_dt",           0 ); // Activates constant time step
    model.stress_rotation    = ReadInt2( fin, "stress_rotation",       1 ); // 0: no stress rotation, 1: analytic rotation, 2: upper convected rate
    model.dt_max             = ReadDou2( fin, "dt_max", 1e20 ) /scaling.t;  // maximum allowed time step, the default value is set to ~infinite, it we become effective only if specificaly set in XXX.txt (see e.g. LithoScale.txt)
    model.dt_min             = ReadDou2( fin, "dt_min",-1e20 ) /scaling.t;  // minimum allowed time step, defaut is negative such that it will never be activated unless specifically set in XXX.txt file
    model.dt_reduction_factor= ReadDou2( fin, "dt_reduction_factor", 1);
    // Physics 
    model.mechanical         = ReadInt2( fin, "mechanical",            1 ); // Activates mechanical solver
    model.advection          = ReadInt2( fin, "advection",             1 ); // Activates advection
    model.elastic            = ReadInt2( fin, "elastic",               0 ); // Activates elasticity
    model.thermal            = ReadInt2( fin, "thermal",               0 ); // Activates thermal solver
    model.anisotropy         = ReadInt2( fin, "anisotropy",            0 ); // Turns on anisotropy
    model.polar              = ReadInt2( fin, "polar",                 0 ); // Activate polar-Cartesian coordinates
    model.finite_strain      = ReadInt2( fin, "finite_strain",         0 ); // Integrates finite stran and save deformation gradient tensor
    model.compressible       = ReadInt2( fin, "compressible",          0 ); // Turns on compressibility
    model.density_variations = ReadInt2( fin, "density_variations",    0 ); // Turns on volume change due to reaction if 1
    model.kinetics           = ReadInt2( fin, "kinetics",              0 ); // Activates reaction kinetics
    model.out_of_plane       = ReadInt2( fin, "out_of_plane",          0 ); // Out-of-plane strain
    model.melting            = ReadInt2( fin, "melting",               0 ); // Activates melting

    // Numerics: linear solver
    model.penalty            = ReadDou2( fin, "penalty",           1.0e3 ); // Penalty factor
    model.auto_penalty       = ReadDou2( fin, "auto_penalty",        0.0 ); // Activates automatic penalty factor computation
    model.diag_scaling       = ReadInt2( fin, "diag_scaling",          1 ); // Activates diagonal scaling
    model.preconditioner     = ReadInt2( fin, "preconditioner",        0 ); // Preconditoner type for Linear solver, 0: Picard preconditionner, -1: symmetrised Picard, 1: symmetrised Newton
    model.lin_abs_div        = ReadDou2( fin, "lin_abs_div",      1.0e-9 ); // Tolerance for linear mechanical solver
    model.lin_rel_div        = ReadDou2( fin, "lin_rel_div",      1.0e-5 ); // Tolerance for linear mechanical solver
    model.lin_abs_mom        = ReadDou2( fin, "lin_abs_mom",      1.0e-9 ); // Tolerance for linear mechanical solver
    model.lin_rel_mom        = ReadDou2( fin, "lin_rel_mom",      1.0e-5 ); // Tolerance for linear mechanical solver
    model.lin_solver         = ReadInt2( fin, "lin_solver",            2 ); // 1: Powell-Hestenes, 2: Powell-Hestenes augmented (killer solver) 
    model.max_its_PH         = ReadInt2( fin, "max_its_PH",           10 ); // Max number of Powell-Hestenes iterations
    // Numerics: non-linear solver
    Nmodel.nit_max           = ReadInt2( fin, "nit_max",               1 ); // Maximum number of iterations
    model.Newton             = ReadInt2( fin, "Newton",                0 ); // Activates Newton iterations
    Nmodel.Picard2Newton     = ReadInt2( fin, "Picard2Newton",         0 ); // Switch from Picard to Newton iterations
    Nmodel.Picard2Newton_tol = ReadDou2( fin, "Picard2Newton_tol",  1e-1 ); // Condition for switching based on residual magnitude
    Nmodel.max_its_Pic       = ReadInt2( fin, "max_its_Pic",          10 ); // Condition for switching based on number of Picard iterations
    Nmodel.let_res_grow      = ReadInt2( fin, "let_res_grow",          0 ); // Allows residual to grow 
    model.max_its_KSP        = ReadInt2( fin, "max_its_KSP",         100 ); // Condition for switching based on number of Picard iterations
    model.rel_tol_KSP        = ReadDou2( fin, "rel_tol_KSP",        1e-4 ); // Relative tolerance for inner Krylov solver
    Nmodel.nonlin_abs_mom    = ReadDou2( fin, "nonlin_abs_mom",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    Nmodel.nonlin_abs_div    = ReadDou2( fin, "nonlin_abs_div",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    Nmodel.nonlin_rel_mom    = ReadDou2( fin, "nonlin_rel_mom",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    Nmodel.nonlin_rel_div    = ReadDou2( fin, "nonlin_rel_div",   1.0e-6 ); // Tolerance for non-linear mechanical solver
    model.min_eta            = ReadDou2( fin, "min_eta", 1e18)/scaling.eta; // Minimum viscosity
    model.max_eta            = ReadDou2( fin, "max_eta", 1e24)/scaling.eta; // Maximum viscosity
    model.eta_tol            = ReadDou2( fin, "eta_tol", 1e11);             // Tolerance for Viscosity calculation
    model.safe_mode          = ReadInt2( fin, "safe_mode",             0 ); // Activates safe mode: reduces time step if convergence fails
    model.safe_dt_div        = ReadDou2( fin, "safe_dt_div",         5.0 ); // Reduction factor for time step reduction
    model.max_num_stag       = ReadInt2( fin, "max_num_stag",          3 ); // maximum number of stagnation (safe mode)
    model.line_search        = ReadInt2( fin, "line_search",           0 ); // Activates line search
    model.line_search_min    = ReadDou2( fin, "line_search_min",     0.0 ); // Minimum alpha value for line search 
    model.residual_form      = ReadInt2( fin, "residual_form",         1 ); // Form of residual - TODO: delete if our models work with new default value (1)
    Nmodel.stagnated         = 0;
    // Numerics: marker-in-cell
    ParticlesInput particles;
    model.ani_average        = ReadInt2( fin, "ani_average",           1 ); // 0: arithmetic mean - 1: harmonic mean - 2: geometric mean
    model.eta_average        = ReadInt2( fin, "eta_average",           0 ); // 0: arithmetic mean - 1: harmonic mean - 2: geometric mean
    model.interp_stencil     = ReadInt2( fin, "interp_stencil",        1 ); // 1: 1-Cell          - 9: 9-Cell
    model.subgrid_diffusion  = ReadInt2( fin, "subgrid_diffusion",     0 ); // 0: No subgrid diffusion, 1: temperature, 2: temperature + stress
    model.conserv_interp     = ReadInt2( fin, "conserv_interp",        0 ); // Activates Taras conservative interpolation
    model.direct_neighbour   = ReadInt2( fin, "direct_neighbour",      0 ); // Direct neighbour interpolation
    model.initial_noise      = ReadInt2( fin, "initial_noise",         0 ); // Add noise on initial marker locations
    model.marker_noise       = ReadInt2( fin, "marker_noise",          0 ); // Background noise field generated and tracked on the particles 
    model.reseed_markers     = ReadInt2( fin, "reseed_markers",        1 ); // Activates reseeding / particle injection
    model.reseed_mode        = ReadInt2( fin, "reseed_mode",           0 ); // Activates reseeding / particle injection 
    particles.Nx_part        = ReadInt2( fin, "Nx_part",               4 ); // number of particle per cell in x
    particles.Nz_part        = ReadInt2( fin, "Nz_part",               4 ); // number of particle per cell in y
    particles.min_part_cell  = ReadInt2( fin, "min_part_cell",        16 ); // minimum number of particle per cell (if below: will trigger reseeding)
    particles.Nb_part        = (model.Nx-1)*(model.Nz-1) * particles.Nx_part * particles.Nz_part;
    particles.Nb_part_max    = 4.1*particles.Nb_part;
    // Boundary conditions
    model.shear_style        = ReadInt2( fin, "shear_style",           0 ); // BC type: 0: pure shear, 2: periodic simple shear
    model.periodic_x         = ReadInt2( fin, "periodic_x",            0 ); // Activates periodicity in x
    model.pure_shear_ALE     = ReadInt2( fin, "pure_shear_ALE",        0 ); // Activates Arbitrary Lagarangian Eulerian mode (pure shear box deformation)
    model.free_surface       = ReadInt2( fin, "free_surface",          0 ); // Activates free surface
    model.free_surface_stab  = ReadDou2( fin, "free_surface_stab",   0.0 ); // Activate free surface stabilisation: range 0.0-2.0
    model.topo_update        = ReadInt2( fin, "topo_update",           1 ); // 0: total topography update (diffusive); 1: incremental
    // Model configurations
    model.initial_cooling    = ReadInt2( fin, "initial_cooling",       0 ); // Activates initial cooling
    model.cooling_duration   = ReadDou2( fin, "cooling_duration",     Ga ) /scaling.t; // Initial cooling duration
    model.shear_heating      = ReadInt2( fin, "shear_heating",         1 ); // Activates shear heating
    model.adiab_heating      = ReadInt2( fin, "adiab_heating",         0 ); // 0: zero, 1: lithostatic P assumption, 2: full derivative
    model.surface_processes  = ReadInt2( fin, "surface_processes",     0 ); // 1: diffusion; 2: diffusion + sedimentation
    model.marker_aniso_angle = ReadInt2( fin, "marker_aniso_angle",    0 ); // Enables setting anisotropy angle per particles rather than phases
    model.layering           = ReadInt2( fin, "layering",              0 ); // Activation of Layering (Anais setup)
    // Transformations
    model.chemical_diffusion  = ReadInt2( fin, "chemical_diffusion",              0 ); // Activate progressive reactions
    model.chemical_production = ReadInt2( fin, "chemical_production",              0 ); // Activate progressive reactions
    model.smooth_softening    = ReadInt2( fin, "smooth_softening",      1 ); // Activates smooth explicit kinematic softening function
    // Background ambient conditions
    model.bkg_strain_rate    = ReadDou2( fin, "bkg_strain_rate", 1e-30)/scaling.E; // Background tectonic rate, default is close to zero to avoid any Nans of Infs in rheology
    model.bkg_div_rate       = ReadDou2( fin, "bkg_div_rate",      0.0)/scaling.E; // Background divergence rate
    model.bkg_pressure       = ReadDou2( fin, "bkg_pressure",      0.0)/scaling.S; // Background pressure
    model.bkg_temperature    = ReadDou2( fin, "bkg_temperature",   0.0)/scaling.T; // Background temperature
    model.fix_temperature    = ReadInt2( fin, "fix_temperature",              0 ); // If 1: calls user defined function that sets temperature at each step  
    // Surface processes
    model.surf_diff          = ReadDou2( fin, "surf_diff",       0.0 ) / (pow(scaling.L,2.0)/scaling.t);
    model.surf_ised1         = ReadInt2( fin, "surf_ised1",      0.0 );
    model.surf_ised2         = ReadInt2( fin, "surf_ised2",      0.0 );
    model.surf_sedirate      = ReadDou2( fin, "surf_sedirate",   0.0 ) / scaling.V;
    model.surf_baselev       = ReadDou2( fin, "surf_baselev",    0.0 ) / scaling.L;
    model.surf_Winc          = ReadDou2( fin, "surf_Winc",       0.0 ) / scaling.L;
    model.surf_Vinc          = ReadDou2( fin, "surf_Vinc",       0.0 ) / scaling.V;
    // Initial stress field
    model.preload           = ReadInt2( fin, "preload",    0 ); // put 1 if you want initial stresses
    model.preload_sxxd      = ReadDou2( fin, "preload_sxxd",    0.0 ) / scaling.S;
    model.preload_szzd      = ReadDou2( fin, "preload_szzd",    0.0 ) / scaling.S;
    model.preload_sxz       = ReadDou2( fin, "preload_sxz",     0.0 ) / scaling.S;
    // Initial thermal perturbation
    model.therm_perturb      = ReadInt2( fin, "therm_perturb",                 0 ); // Includes initial thermal perbation
    model.therm_perturb_x0   = ReadDou2( fin, "therm_perturb_x0",  0.0 )/scaling.L; // x position
    model.therm_perturb_z0   = ReadDou2( fin, "therm_perturb_z0",  0.0 )/scaling.L; // y position
    model.therm_perturb_rad_x= ReadDou2( fin, "therm_perturb_rad_x", 0.0 )/scaling.L; // Radius
    model.therm_perturb_rad_z= ReadDou2( fin, "therm_perturb_rad_z", 0.0 )/scaling.L; // Radius
    model.therm_perturb_dT   = ReadDou2( fin, "therm_perturb_dT" , 0.0 )/scaling.T; // Temperature anomaly
    // For rheological database reasons...
    model.force_act_vol_ast  = ReadInt2( fin, "force_act_vol_ast",   0 ); // if 1 then:
    model.act_vol_dis_ast    = ReadDou2( fin, "act_vol_dis_ast" ,  0.0 ); // ... set dislocation creep to value
    model.act_vol_dif_ast    = ReadDou2( fin, "act_vol_dif_ast" ,  0.0 ); // ... set diffusion creep to value
    model.force_melt_weak    = ReadInt2( fin, "force_melt_weak",     0 ); // if 1 then:
    model.melt_weak          = ReadDou2( fin, "melt_weak",         0.0 ); // ... set melt waekening factor
    // Model user's delights
    model.user0              = ReadDou2( fin, "user0",           0.0 );
    model.user1              = ReadDou2( fin, "user1",           0.0 );
    model.user2              = ReadDou2( fin, "user2",           0.0 );
    model.user3              = ReadDou2( fin, "user3",           0.0 );
    model.user4              = ReadDou2( fin, "user4",           0.0 );
    model.user5              = ReadDou2( fin, "user5",           0.0 );
    model.user6              = ReadDou2( fin, "user6",           0.0 );
    model.user7              = ReadDou2( fin, "user7",           0.0 );
    model.user8              = ReadDou2( fin, "user8",           0.0 );
    model.user9              = ReadDou2( fin, "user9",           0.0 );
    // Static loading before ADVECTION and TRANSFORMATION?
    model.time_advection_start = ReadDou2( fin, "time_advection_start", 0.0 ) * (365.25*24.0*3600.0) / scaling.t; // time in yrs when advection = 1 (advection = 0 before)
    model.time_switch_phase    = ReadDou2( fin, "time_switch_phase",    0.0 ) * (365.25*24.0*3600.0) / scaling.t; // time in yrs when phase switch from phase_switch_init to phase_switch_final
    model.phase_switch_init    = ReadInt2( fin, "phase_switch_init",    0 ); // 2
    model.phase_switch_final   = ReadInt2( fin, "phase_switch_final",   0 ); // 3
    
    // Derived quantities
    model.dx                 = (model.xmax - model.xmin) / (model.Nx - 1);
    model.dz                 = (model.zmax - model.zmin) / (model.Nz - 1);
    model.dt0                = model.dt;
    model.dt_start           = model.dt;
    model.xmin0              = model.xmin;
    model.zmin0              = model.zmin;
    model.xmax0              = model.xmax;
    model.zmax0              = model.zmax;
    // For Cindy's setup
    model.diffuse_X          = ReadInt2( fin, "diffuse_X",     0 );              // 0 or 1
    model.diffuse_avg        = ReadInt2( fin, "diffuse_avg",   0 );              // 0: arithmetic mean - 1: harmonic mean - 2: geometric mean
    model.diffusion_length   = ReadDou2( fin, "diffusion_length",  0.0 ) / scaling.L;
    // Gravity
    model.gx                 = ReadDou2( fin, "gx",  0.0 ) / scaling.a;
    model.gz                 = ReadDou2( fin, "gz",  0.0 ) / scaling.a;
    // Consequential behaviour
    if (model.interp_stencil!=1 && model.interp_stencil!=9) { printf("Wrong value of interp_stencil: should be 1 or 9.\n"); exit(1); }
    if ( model.shear_style == 1 ) model.periodic_x     = 1; // If simple shear, it must  be periodic in x
    if ( model.shear_style == 0 ) model.periodic_x     = 0; // If simple shear, it can't be periodic in x
    if ( model.anisotropy  == 1 ) model.finite_strain  = 1; // If anisotropy, then also track finite strain
    // if ( model.anisotropy  == 1 ) model.preconditioner = -1; // If anisotropy, then use symmetrised preconditioner
    // printf("%d\n", model.preconditioner); exit(1);
    if (model.Newton==0) Nmodel.Picard2Newton = 0; // If Picard is activated, do not switch to Newton
    if (model.Newton==1) model.line_search    = 1; // If Newton is activated, switch to line search
    if ( model.lin_solver == 0 || model.Newton == 1 || model.anisotropy == 1) {
        printf("WARNING!! Changing from solver type 0 to solver type 2!!! That's the new standard in MDOODZ 6.0.\n");
        model.lin_solver = 2;
    }
    //------------------------------------------------------------------------------------------------------------------------------//
    // DEFORMATION MAP PARAMETERS
    //------------------------------------------------------------------------------------------------------------------------------//
    // Create deformation maps or not (default no)
    model.deformation_maps = ReadInt2(fin, "deformation_maps", 0);
    // Resolution
    model.nT               = ReadInt2(fin, "nT", 11);
    model.nE               = ReadInt2(fin, "nE", 11);
    model.nd               = ReadInt2(fin, "nd", 11);
    // Temperature, strain rate, grain size MIN/MAX & pressure
    model.Tmin             = ReadDou2(fin, "Tmin", 100.0);  model.Tmin += zeroC; model.Tmin /= scaling.T;    // C . K & non-dimensionalization
    model.Tmax             = ReadDou2(fin, "Tmax", 1000.0); model.Tmax += zeroC; model.Tmax /= scaling.T;    // C . K & non-dimensionalization
    model.Emin             = ReadDou2(fin, "Emin", -30.0);   //model.Emin /= scaling.E;                           // Non-dimensionalization
    model.Emax             = ReadDou2(fin, "Emax", -4.0 );   //model.Emax /= scaling.E;                           // Non-dimensionalization
    model.dmin             = ReadDou2(fin, "dmin", -7.0 );   //model.dmin /= scaling.L;                           // Non-dimensionalization
    model.dmax             = ReadDou2(fin, "dmax", -2.0 );   //model.dmax /= scaling.L;                           // Non-dimensionalization
    model.Pn               = ReadDou2(fin, "Pn",  5.0e8 );      model.Pn   /= scaling.S;                          // Non-dimensionalization
    //------------------------------------------------------------------------------------------------------------------------------//
    // MATERIAL PROPERTIES
    //------------------------------------------------------------------------------------------------------------------------------//
    model.Nb_phases = materials.Nb_phases =  ReadInt2( fin, "Nb_phases", 0 );
    for ( k=0; k<materials.Nb_phases; k++) {
        // Read general parameters
        materials.rho[k]  = ReadMatProps( fin, "rho", k,   2700.0 )  / scaling.rho;
        materials.G[k]    = ReadMatProps( fin, "G",  k,   1.0e10  )  / scaling.S;
        materials.Cp[k]   = ReadMatProps( fin, "Cp",  k,   1.0e3  )  / scaling.Cp;
        materials.k[k]    = ReadMatProps( fin, "k",   k,   1.0e-6 )  / scaling.k;
        materials.k_eff[k]= materials.k[k];
        materials.Qr[k]   = ReadMatProps( fin, "Qr",  k,   1.0e-30)  / (scaling.W / pow(scaling.L,3.0));
        materials.alp[k]  = ReadMatProps( fin, "alp", k,      0.0 )  / (1.0/scaling.T);
        materials.bet[k]  = ReadMatProps( fin, "bet", k,  1.0e-40 )  / (1.0/scaling.S);
        materials.drho[k] = ReadMatProps( fin, "drho",k,      0.0 )  / (scaling.rho);
        materials.T0[k]   = (zeroC) / (scaling.T); 
        materials.P0[k]   = 1e5 / (scaling.S);
        // Read melting model
        materials.melt[k]     = (int)ReadMatProps( fin, "melt",     k,    0.0   );     // 0: No plasticity --- >1: Yes, the type should be selected accordingly 
        // Read plasticity switches
        materials.plast[k]    = (int)ReadMatProps( fin, "plast",    k,    1.0   );     // 0: No plasticity --- 1: Yes 
        materials.yield[k]    = (int)ReadMatProps( fin, "yield",    k,    1.0   );     // 1: Drucker-Prager
        // Read plasticity parameters
        materials.sig_tens[k] = ReadMatProps( fin, "sig_tens", k,  1.0e7 )/ scaling.S; // Tension stress for hyperbolic Drucker-Prager
        materials.sig1[k]     = ReadMatProps( fin, "sig1",     k,  1.0e7 )/ scaling.S; // Transition stress to 0 dilatancy for hyperbolic Drucker-Prager
        materials.dsig1[k]    = ReadMatProps( fin, "dsig1",    k,  1.0e7 )/ scaling.S; // Width of transition for hyperbolic Drucker-Prager
        materials.C[k]        = ReadMatProps( fin, "C",        k,  1.0e7 )  / scaling.S;
        materials.phi[k]      = ReadMatProps( fin, "phi",      k,   30.0 )  * M_PI/ 180.0;
        materials.psi[k]      = ReadMatProps( fin, "psi",      k,    0.0 )  * M_PI/ 180.0;
        if (materials.psi[k]>0.0 && model.compressible==0) { printf("Set compressible=1 to activate dilation\n"); exit(1); }
        materials.Slim[k] = ReadMatProps( fin, "Slim" ,k,  1.0e90 )  / scaling.S;
        if (materials.Slim[k]<materials.C[k]) {
            printf("Upper stress limiter of phase %d is lower than cohesion: it makes no sense, please correct!\n", k);
            exit(1);
        }
        // Viscoplasticity
        materials.n_vp[k]      = ReadMatProps( fin, "n_vp",   k,       1.0 ) ;
        materials.eta_vp[k]    = ReadMatProps( fin, "eta_vp", k,       0.0 ) / scaling.S / pow(scaling.t, 1.0/materials.n_vp[k]);
        if ( materials.n_vp[k]>1.0001 ) {
            printf("Power-law viscoplasticity is not implemented in MD7, please complain!\n");
            exit(1);
        }
        // Strain softening
        materials.coh_soft[k]   = (int)ReadMatProps( fin, "coh_soft",   k,    0.0   );
        materials.phi_soft[k]   = (int)ReadMatProps( fin, "phi_soft",   k,    0.0   );
        materials.psi_soft[k]   = (int)ReadMatProps( fin, "psi_soft",   k,    0.0   );
        if (materials.psi_soft[k]>0 && model.compressible==0) { printf("Set compressible=1 to activate dilation softening\n"); exit(1); }
        materials.C_end[k]     = ReadMatProps( fin, "Ce",     k,    materials.C[k]*scaling.S    ) / scaling.S;      // TODO: rename Ce --> C_end
        materials.phi_end[k]   = ReadMatProps( fin, "phie",   k,    materials.phi[k]*180.0/M_PI  ) * M_PI / 180.0;  // TODO: rename phie --> phi_end
        materials.psi_end[k]   = ReadMatProps( fin, "psie",   k,    materials.psi[k]*180.0/M_PI  ) * M_PI / 180.0;  // TODO: rename psie --> psi_end
        double eps_coh = 1.0 / scaling.S;
        double eps_phi = 0.1  * M_PI / 180.0;
        double eps_psi = 0.1  * M_PI / 180.0;
        printf("%d materials.coh_soft[k]=%d  materials.C[k] = %2.2e  materials.C_end[k] = %2.2e\n", k, materials.coh_soft[k],  materials.C[k]*scaling.S,  materials.C_end[k]*scaling.S);
        if ( materials.coh_soft[k] == 1 && fabs( materials.C_end[k]   - materials.C[k]  ) < eps_coh ) { printf("Please set a difference in cohesion, if not set coh_soft of phase %d to 0.0\n", k); exit(122); };
        if ( materials.phi_soft[k] == 1 && fabs( materials.phi_end[k] - materials.phi[k]) < eps_phi ) { printf("Please set a difference in friction angle, if not set phi_soft of phase %d to 0.0\n", k); exit(122); };
        if ( materials.psi_soft[k] == 1 && fabs( materials.psi_end[k] - materials.psi[k]) < eps_psi ) { printf("Please set a difference in dilation angle, if not set psi_soft of phase %d to 0.0\n", k); exit(122); };
        materials.pls_start[k] = ReadMatProps( fin, "plss",   k,    1.0    );
        materials.pls_end[k]   = ReadMatProps( fin, "plse",   k,    2.0    );
        // Read flow law parameters
        materials.cstv[k]     = ReadMatProps( fin,     "cstv", k,    0.0  );
        materials.pwlv[k]     = ReadMatProps( fin,     "pwlv", k,    0.0  );
        materials.linv[k]     = ReadMatProps( fin,     "linv", k,    0.0  );
        materials.gbsv[k]     = ReadMatProps( fin,     "gbsv", k,    0.0  );
        materials.expv[k]     = ReadMatProps( fin,     "expv", k,    0.0  );
        gsel                  = ReadMatProps( fin,     "gsel", k,    0.0  );
        materials.eps0[k]     = ReadMatProps( fin,     "eps0", k,  1e-14  ) / scaling.E;
        materials.tau0[k]     = ReadMatProps( fin,     "tau0", k,    2e6  ) / scaling.S;
        materials.eta0[k]     = ReadMatProps( fin,     "eta0", k, 1.0e20  ) / scaling.eta;
        materials.npwl[k]     = ReadMatProps( fin,     "npwl", k,    1.0  );
        materials.Qpwl[k]     = ReadMatProps( fin,     "Qpwl", k,    0.0  );
        materials.pref_pwl[k] = ReadMatProps( fin, "pref_pwl", k,    1.0  );    // weakening prefactor for power law
        materials.gs[k]       = ReadMatProps( fin,       "gs", k,    0.0  );
        materials.gs_ref[k]   = ReadMatProps( fin,   "gs_ref", k, 2.0e-3  ) / scaling.L;
        materials.kin[k]      = ReadMatProps( fin,      "kin", k,    0.0  );
        // Reaction stuff
        materials.reac_soft[k]  = (int)ReadMatProps( fin, "reac_soft",   k,    0.0  );
        materials.reac_phase[k] = (int)ReadMatProps( fin, "reac_phase",   k,    0.0  );
        if (materials.reac_phase[k]>=model.Nb_phases && materials.reac_soft[k]==1) {printf("The target phase for reaction softening of phase %d is not set up. Edit the .txt file!\n", k  ); exit(13);}
        materials.Pr[k]         = ReadMatProps( fin, "Pr",      k,    0.0  ) / scaling.S;
        materials.dPr[k]        = ReadMatProps( fin, "dPr",     k,    0.0  ) / scaling.S;
        materials.tau_kin[k]    = ReadMatProps( fin, "tau_kin", k,    0.0  ) / scaling.t;
        materials.k_chem[k]     = ReadMatProps( fin, "k_chem",  k,    0.0  ) / (scaling.L*scaling.L/scaling.t);
        // Density models
        materials.density_model[k]     = (int)ReadMatProps( fin, "density_model",     k,    3  );
        materials.phase_diagram[k]     = (int)ReadMatProps( fin, "phase_diagram",     k,   -1  );
        // Diffused rheological contrasts
        materials.phase_mix[k]  = (int)ReadMatProps( fin, "phase_mix",k,          0.0  );
        materials.phase_two[k]  = (int)ReadMatProps( fin, "phase_mix",k,    (double)k  );
        // Anisotropy
        materials.aniso_angle[k]    =  ReadMatProps( fin, "aniso_angle",     k,   90.0  )  * M_PI/ 180.0;
        materials.aniso_factor[k]   =  ReadMatProps( fin, "aniso_factor",    k,    1.0  ); 
        // materials.ani_fac_v[k]      =  ReadMatProps( fin, "ani_fac_v",    k,    1.0  );        // viscous anisotropy strength
        // materials.ani_fac_e[k]      =  ReadMatProps( fin, "ani_fac_e",    k,    1.0  );        // elastic anisotropy strength
        // materials.ani_fac_p[k]      =  ReadMatProps( fin, "ani_fac_p",    k,    1.0  );        // plastic anisotropy strength
        materials.ani_fac_max[k]    =  ReadMatProps( fin, "ani_fac_max",  k, 1000.0  );        // maximum anisotropy strength
        materials.ani_fstrain[k]    =  (int)ReadMatProps( fin, "ani_fstrain",    k,    0.0  ); // strain dependent anisotropy per phase
        materials.axx[k]            =  ReadMatProps( fin, "axx",    k,    1.0  );              // plastic directional factor xx
        materials.azz[k]            =  ReadMatProps( fin, "azz",    k,    1.0  );              // plastic directional factor zz
        materials.ayy[k]            =  ReadMatProps( fin, "ayy",    k,    1.0  );              // plastic directional factor yy
        // temperature driven transition between phases
        materials.transmutation[k]             = (int)ReadMatProps( fin, "transmutation",    k,    0  ); // switcher: 0 - no transition; -1 transition happens when below transition_temperature; 1 transition happens when above transition_temperature;
        materials.transmutation_phase[k]       = (int)ReadMatProps( fin, "transmutation_phase",    k,    1  ); // phase to set
        materials.transmutation_temperature[k] = ReadMatProps( fin, "transmutation_temperature",    k,    1300.0  ) / scaling.T; // temperature boundary
        // Check if any flow law is active
        int sum = abs(materials.cstv[k]) + abs(materials.pwlv[k]) + abs(materials.linv[k]) + abs(materials.gbsv[k]) + abs(materials.expv[k]);
        if ( sum == 0 ) {
            printf ("Phase %0d has no determined flow mechanism\n Simulation will end now!\n", k);
            exit(12);
        }
        // Print material parameters
        printf("----------------------------------------- MODEL DOMAIN ------------------------------------------\n");
        printf("Xmin   = %2.1lf  km         Xmax   = %2.1lf  km     Nx   = %3d    dx   = %.2lf m\n", (model.xmin*scaling.L)/1e3, (model.xmax*scaling.L)/1e3, model.Nx, model.dx*scaling.L);
        printf("Zmin   = %2.1lf  km         Zmax   = %2.1lf  km      Nz   = %3d    dz   = %.2lf m\n", (model.zmin*scaling.L)/1e3, (model.zmax*scaling.L)/1e3, model.Nz, model.dz*scaling.L );
        printf("-------------------------------------------- PHASE: %d -------------------------------------------\n", k);
        printf("rho    = %2.2e kg/m^3     G = %2.2e Pa\n", materials.rho[k]*scaling.rho, materials.G[k]*scaling.S );
        printf("Cp     = %2.2e J/kg/K      k = %2.2e W/m/K      Qr = %2.2e W/m3\n", materials.Cp[k]*scaling.Cp, materials.k[k]*scaling.k, materials.Qr[k]*(scaling.W / pow(scaling.L,3)) );
        printf("C      = %2.2e Pa        phi = %2.2e deg      Slim = %2.2e Pa\n",  materials.C[k]*scaling.S, materials.phi[k]*180/M_PI, materials.Slim[k]*scaling.S );
        printf("alp    = %2.2e 1/T        T0 = %2.2e K         bet = %2.2e 1/Pa       P0 = %2.2e Pa       drho = %2.2e kg/m^3 \n", materials.alp[k]*(1/scaling.T), materials.T0[k]*(scaling.T), materials.bet[k]*(1/scaling.S), materials.P0[k]*(scaling.S), materials.drho[k]*scaling.rho );
        printf("prefactor for power-law: %2.2e\n", materials.pref_pwl[k]);
                 printf("C_end    = %2.2e Pa        Phi_end = %2.2e deg         pls_start = %2.2e        pls_end = %2.2e \n", materials.C_end[k]*scaling.S, materials.phi_end[k]*180/M_PI, materials.pls_start[k],  materials.pls_end[k] );
        printf("eta0_vp   = %2.2e  Pa.s^(1/n)         n_vp   = %2.2e\n", materials.eta_vp[k]* (scaling.S*pow(scaling.t,1.0/materials.n_vp[k])) , materials.n_vp[k]);
        printf("aniso_factor = %2.2e []        aniso_angle = %2.2e deg         ani_fac_max = %2.2e []        ani_fstrain = %d\n", materials.aniso_factor[k], materials.aniso_angle[k]/(M_PI/180.0), materials.ani_fac_max[k], materials.ani_fstrain[k]);

        printf("Flow law settings:\n");
        if ( abs(materials.cstv[k])>0 ) printf("--->    Constant viscosity activated \n");
        if ( abs(materials.pwlv[k])>0 ) printf("--->   Power law viscosity activated \n");
        if ( abs(materials.linv[k])>0 ) printf("--->      Linear viscosity activated \n");
        if ( abs(materials.gbsv[k])>0 ) printf("--->         GBS viscosity activated \n");
        if ( abs(materials.expv[k])>0 ) printf("---> Exponential viscosity activated \n");
        if ( abs(gsel)              >0 ) printf("--->  Grain size evolution activated \n");

        // Call flow law data base
        if ( abs(materials.pwlv[k])>0 ) ReadDataPowerLaw   ( &materials, &model, k, materials.pwlv[k], &scaling );
        if ( abs(materials.linv[k])>0 ) ReadDataLinear     ( &materials, &model, k, materials.linv[k], &scaling );
        if ( abs(materials.gbsv[k])>0 ) ReadDataGBS        ( &materials, &model, k, materials.gbsv[k], &scaling );
        if ( abs(materials.expv[k])>0 ) ReadDataExponential( &materials, &model, k, materials.expv[k], &scaling );
        if ( abs(materials.gs[k])  >0 ) ReadDataGSE        ( &materials, &model, k, materials.gs[k],   &scaling );
        if ( abs(materials.kin[k]) >0 ) ReadDataKinetics   ( &materials, &model, k, materials.kin[k],  &scaling );
    }

    materials.R = Rg / (scaling.J/scaling.T);

    //------------------------------------------------------------------------------------------------------------------------------//
    // PHASE DIAGRAM INFO - simple pressure-dependent density model == 4 --- for quartz/coesite study
    //------------------------------------------------------------------------------------------------------------------------------//

    model.PD1DnP  = DoodzCalloc( model.Nb_phases, sizeof(int));
    model.PD1Drho = DoodzCalloc( model.Nb_phases, sizeof(double*));
    model.PD1Dmin = DoodzCalloc( model.Nb_phases, sizeof(double));
    model.PD1Dmax = DoodzCalloc( model.Nb_phases, sizeof(double));

    for ( k=0; k<materials.Nb_phases; k++) {

        if ( materials.density_model[k] == 4 ) {

            printf("Phase #%d has density_model %d, so:\n", k, materials.density_model[k]);
            printf("Loading 1D diagram for coesite quartz only baby...\n");

            model.PD1DnP[k]         = 10000;
            model.PD1Dmin[k]        = 1e5/scaling.S;
            model.PD1Dmax[k]        = 1e10/scaling.S;
            model.PD1Drho[k] = DoodzCalloc( model.PD1DnP[k], sizeof(double));
            char *fname = "PHASE_DIAGRAMS/QuartzCoesite600C.bin";
            if ( fopen(fname, "rb") != NULL ) {
                FILE* read = fopen(fname, "rb");
                fread ( model.PD1Drho[k],     sizeof(double), model.PD1DnP[k] , read);
                fclose(read);

                double p[model.PD1DnP[k]];
                double dP = ( model.PD1Dmax[k] - model.PD1Dmin[k]) / (model.PD1DnP[k] -1);
                p[0] =  model.PD1Dmin[k];

                // Dimensionalize
                for (int i=0; i<model.PD1DnP[k]; i++) {
                    model.PD1Drho[k][i] /= scaling.rho;
                    if (i>0) p[i] = p[i-1] + dP;
                }

                // Apply some diffusion...
                const double Kdiff   = 1e0;              // diffusivity
                const double dt_exp  =  dP*dP/Kdiff/2.1; // explicit time step
                const int    n_steps = 200;              // number of steps
                double qxW, qxE;
                double *rho0   = DoodzCalloc( model.PD1DnP[k], sizeof(double));

                for (int it=0; it<n_steps; it++) {
                    ArrayEqualArray( rho0, model.PD1Drho[k], model.PD1DnP[k]);
                    for (int ix = 1; ix<model.PD1DnP[k]-1; ix++) {
                        qxW = - Kdiff*(rho0[ix] - rho0[ix-1])/dP;
                        qxE = - Kdiff*(rho0[ix+1] - rho0[ix])/dP;
                        model.PD1Drho[k][ix] = rho0[ix] - dt_exp*(qxE - qxW)/dP;
                    }
                }
                DoodzFree(rho0);

//                printf("scaling.rho = %2.10e\n",scaling.rho);
//                FILE* write = fopen("PHASE_DIAGRAMS/QuartzCoesite600C_smoothed.bin", "wb");
//                fwrite ( model.PD1Drho[k],     sizeof(double), model.PD1DnP[k] , write);
//                fclose(write);

//                //-------- In situ VISU for debugging: do not delete ------------------------//
//                FILE        *GNUplotPipe;
//                GNUplotPipe = popen ("gnuplot -persistent", "Work");
//                int NumCommands = 2;
//                char *GNUplotCommands[] = { "set style line 1 lc rgb '#0060ad' lt 1 lw 2 pt 7 pi -1 ps 1.5", "set pointintervalbox 3"};
//                for (int i=0; i<NumCommands; i++) fprintf(GNUplotPipe, "%s \n", GNUplotCommands[i]); //Send commands to gnuplot one by one.
//
//                fprintf(GNUplotPipe, "plot '-' with lines linestyle 1\n");
//                for (int i=0; i<model.PD1DnP[k]; i++) {
//                    fprintf(GNUplotPipe, "%lf %lf \n", p[i]*scaling.S, model.PD1Drho[k][i]*scaling.rho); //Write the data to a temporary file
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
    model.isPD = 0;

    for ( k=0; k<materials.Nb_phases; k++) {
        printf("Phase %d ---  density model %d \n",k, materials.density_model[k]);
        if ( materials.density_model[k] == 2 ) {
            printf("The density model of phase %0d relies on phase diagrams\n", k);
            model.isPD = 1;
        }
        if ( materials.density_model[k] == 2 && materials.phase_diagram[k] == -1 ) {
            printf("However the phase diagram index was not set (default -1)\n");
            printf("Therefore the simulation will not run, go fix 'phase_diagram' of phase %02d\n", k);
            exit(5);
        }
    }

    // Phase diagrams
    if ( model.isPD == 1 ) {

        printf("Loading phase diagrams...\n");
        int pid;
        model.num_PD = 10;

        // Allocate
        AllocatePhaseDiagrams( &model );

        /**** PHASE DIAGRAMS #00 - Mantle (Jenadi_stx.dat)  ****/
        pid                       = 0;         // Kaus & Connolly, 2005: Effect of mineral phase transitions on sedimentary basin subsidence and uplift
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1000;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1000;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 473.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 2273.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 100e6/scaling.S;         // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 15e9 /scaling.S;         // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin( model.import_files_dir, "Hawaiian_Pyrolite_rho_bin.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #01 - Mantle (Jenadi_stx_HR.dat)  ****/
        pid                       = 1;         // Kaus & Connolly, 2005: Effect of mineral phase transitions on sedimentary basin subsidence and uplift
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1500;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1500;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 273.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 2273.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 1e5/scaling.S;           // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 25e9 /scaling.S;         // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin( model.import_files_dir, "Hawaiian_Pyrolite_HR_rho_bin.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #02 - Basalt (MORB_L.dat)  ****/
        pid                       = 2;  // Water saturated MORB - Bulk composition taken from Schmidt & Poli 1998 EPSL (Table 1)
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1000;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1000;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 573.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1273.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 100e6/scaling.S;         // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.1e9 /scaling.S;        // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "MORB_H2Osat_rho_bin.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #03 - Andesite (Andesite.dat)  ****/
        pid                       = 3;  // Andesite
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1249;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1249;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "Andesite.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #04 - Hydrated Peridotite (Hydrated_Pdt.dat)  ****/
        pid                      = 4;  // Hydrated peridotite
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 793;                     // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 793;                     // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "Hydrated_Pdt.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #05 - MORB (MORB.dat)  ****/
        pid                       = 5;  // MORB
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1249;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1249;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "MORB.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #06 - Pelite (Pelite.dat)  ****/
        pid                      = 6;  // Pelite
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1249;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1249;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "Pelite.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #07 - Rhyolite (Rhyolite.dat)  ****/
        pid                       = 7;  // Rhyolite
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 1249;                    // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 1249;                    // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "Rhyolite.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #08 - Serpentinite (Serpentinite.dat)  ****/
        pid                       = 8;  // Serpentinite
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 793;                     // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 793;                     // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = 373.0/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = 1373.0/scaling.T;        // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = 10.13e6/scaling.S;       // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = 5.5816e9 /scaling.S;     // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "Serpentinite.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);

        /**** PHASE DIAGRAMS #09 - Si02 (Si02.dat)  ****/
        pid                       = 9;  // Si02
        if (pid > (model.num_PD-1) ) {
            printf ("One should increment 'model.num_PD' to allocate enough memory and store the database\n");
            exit(5);
        }
        model.PDMnT[pid]         = 675;                          // Resolution for temperature (MANTLE) []
        model.PDMnP[pid]         = 2500;                         // Resolution for pressure    (MANTLE) []
        model.PDMTmin[pid]       = (398+1e-3)/scaling.T;         // Minimum temperature        (MANTLE) [K]
        model.PDMTmax[pid]       = (800+273)/scaling.T;          // Maximum temperature        (MANTLE) [K]
        model.PDMPmin[pid]       = (0.0090/10*1e9)/scaling.S;    // Minimum pressure           (MANTLE) [Pa]
        model.PDMPmax[pid]       = (49.9890/10*1e9)/scaling.S;   // Maximum pressure           (MANTLE) [Pa]
        model.PDMrho[pid]        = ReadBin(model.import_files_dir, "SiO2_nsm010.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);
        // model.PDMrho[pid]        = ReadBin(model.import_files_dir, "SiO2_nsm000.dat", model.PDMnT[pid], model.PDMnP[pid], scaling.rho);
    }

    //------------------------------------------------------------------------------------------------------------------------------//
    // KINETIC DATA
    //------------------------------------------------------------------------------------------------------------------------------//

    if ( model.kinetics == 1 ) {

        printf("Loading kinetic data...\n");

        /**** Quartz Coesite (dG_QuartzCoesite.dat)  ****/
        model.kin_nT        = 675;                         // Resolution for temperature (MANTLE) []
        model.kin_nP        = 2500;                        // Resolution for pressure    (MANTLE) []
        model.kin_Tmin      = (398+1e-3)/scaling.T;        // Minimum temperature        (MANTLE) [K]
        model.kin_Tmax      = (800+273)/scaling.T;         // Maximum temperature        (MANTLE) [K]
        model.kin_Pmin      = (0.0090/10*1e9)/scaling.S;   // Minimum pressure           (MANTLE) [Pa]
        model.kin_Pmax      = (49.9890/10*1e9)/scaling.S;  // Maximum pressure           (MANTLE) [Pa]
        model.kin_dG        = ReadBin(model.import_files_dir, "dG_QuartzCoesite.dat", model.kin_nT, model.kin_nP, scaling.J);
    }

    // Close input file
    fclose(fin);
    return (Input) {
      .model     = model,
      .particles = particles,
      .scaling   = scaling,
      .Nmodel    = Nmodel,
      .materials = materials,
    };
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
    scale->Cp   = scale->J / scale->m / scale->T;
    scale->k    = scale->W / scale->L / scale->T;
    scale->rhoE = scale->J / scale->m * scale->rho;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double* ReadBin(char importDir[], char A_name[], int nx, int ny, double scale ){
    double *A;
    char* bname; size_t nb_elems = nx*ny; FILE* fid; asprintf(&bname, "%s/%s", importDir, A_name);

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
            rewind(fin);
            free(string1);
            free(param1);
            free(param2);
            // return Default;
            asprintf(&string2, "%s", Default);
            // int str_size = strlen(Default);
            // int h1;
            // string2 = malloc((str_size+1)*sizeof(char));
            // for (h1=0; h1<str_size; h1++) {
            //     string2[h1] = Default[h1];
            // }
            // free(string1);
            // free(param1);
            // free(param2);
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
    // Declarations
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
    while(value == 0.0)
    {
        // Read new line
        fgets ( line, sizeof(line), fin );
        if (feof(fin)) {
            // printf("Warning : Parameter '%s' not found in the setup file, running with default value %2.2e\n", FieldName, Default);
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

                    // Break in case the parameter has not been defined before the next phase starts.
                    if ( strcmp(param3,"ID") == 0) {
                        if ( fabs(Default) <  100.0 ) printf("Warning : Parameter '%s' not found in the setup file, running with default value %.2lf\n", FieldName, Default);
                        if ( fabs(Default) >= 100.0 ) printf("Warning : Parameter '%s' not found in the setup file, running with default value %2.2e\n", FieldName, Default);
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

                    // Break in case the parameter has not been defined until the end of the file.
                    if (feof(fin))
                    {
                        if ( fabs(Default) <  100.0 ) printf("Warning : Parameter '%s' not found in the setup file, running with default value %.2lf\n", FieldName, Default);
                        if ( fabs(Default) >= 100.0 ) printf("Warning : Parameter '%s' not found in the setup file, running with default value %2.2e\n", FieldName, Default);
                        rewind (fin);
                        free(param1);
                        free(param2);
                        free(param3);
                        return Default;
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
