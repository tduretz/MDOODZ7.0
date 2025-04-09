// =========================================================================
// MDOODZ - Visco-Elasto-Plastic Thermo-Mechanical solver
//
// Copyright (C) 2021  MDOODZ Developper team
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
#include "time.h"
#include "mdoodz-private.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void InitialiseGrainSizeParticles( markers* particles, mat_prop *materials ){
    // Set initial particle grain size
#pragma omp parallel for shared( particles, materials )
    for( int k=0; k<particles->Nb_part; k++ ) {
        const int phase = particles->phase[k];
        if (phase != -1) particles->d[k] = materials->gs_ref[phase];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// // Old deviatoric stress field and pressure --- DEPRECATED
// void OldDeviatoricStressesPressure( grid* mesh, markers* particles, scale scaling, params* model ) {
    
//     int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
//     double *sxx0,  *syy0, *szz0, *sxz0;
//     int    cent=1, vert=0, prop=1, interp=0;
    
//     Nx  = mesh->Nx;
//     Nz  = mesh->Nz;
//     Ncx = Nx-1;
//     Ncz = Nz-1;
    
//     sxx0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     syy0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     szz0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     sxz0 = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
//     P2Mastah( model, *particles, particles->sxxd,    mesh, sxx0,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil);
//     P2Mastah( model, *particles, particles->szzd,    mesh, szz0,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil);
//     P2Mastah( model, *particles, particles->syy,     mesh, syy0,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil);
//     P2Mastah( model, *particles, particles->sxz,     mesh, mesh->sxz0,   mesh->BCg.type,  1, 0, interp, vert, model->interp_stencil);
    
// #pragma omp parallel for shared( mesh, sxx0, syy0, szz0 ) private( k, k1, l, c0, c1, c2 ) firstprivate( Nx, Ncx, Ncz )
//     for ( k1=0; k1<Ncx*Ncz; k1++ ) {
//         k  = mesh->kp[k1];
//         l  = mesh->lp[k1];
//         c0 = k  + l*(Nx-1);
//         c1 = k  + l*(Nx);
//         c2 = k  + l*(Nx+1);
        
//         if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
//             mesh->p0_n[c0]       = -1.0/3.0*(sxx0[c0] + syy0[c0] + szz0[c0] );
//             mesh->sxxd0[c0]      =  mesh->p0_n[c0] + sxx0[c0];
//             mesh->szzd0[c0]      =  mesh->p0_n[c0] + szz0[c0];
//         }
//     }
    
//     InterpCentroidsToVerticesDouble( mesh->sxxd0,  mesh->sxxd0_s, mesh, model );
//     InterpCentroidsToVerticesDouble( mesh->szzd0,  mesh->szzd0_s, mesh, model );
//     InterpCentroidsToVerticesDouble( mesh->p0_n,   mesh->p0_s,    mesh, model );
//     InterpVerticesToCentroidsDouble( mesh->sxz0_n, mesh->sxz0,    mesh, model );
    
//     DoodzFree( sxx0 );
//     DoodzFree( szz0 );
//     DoodzFree( sxz0 );
// }

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// // Total stress field --- DEPRECATED
// void TotalStresses( grid* mesh, markers* particles, scale scaling, params* model ) {
    
//     int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, k1;
//     double *dsxx, *dsyy, *dszz, *dsxz, sxx, syy, szz, tyy, tyy0, sxx0, syy0, szz0;
//     double *mdsxx, *mdszz, *mdsyy, *mdsxz;
    
//     Nx  = mesh->Nx;
//     Nz  = mesh->Nz;
//     Ncx = Nx-1;
//     Ncz = Nz-1;
    
//     dsxx = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dsyy = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dszz = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dsxz = DoodzCalloc(Nx *Nz , sizeof(DoodzFP));
    
//     mdsxx = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     mdsyy = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     mdszz = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     mdsxz = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    
// #pragma omp parallel for shared( mesh, dsxx, dsyy, dszz ) private( k, k1, l, c0, c1, c2, sxx, syy, szz, tyy, tyy0, sxx0, syy0, szz0  ) firstprivate( Nx, Ncx, Ncz )
//     for ( k1=0; k1<Ncx*Ncz; k1++ ) {
//         k  = mesh->kp[k1];
//         l  = mesh->lp[k1];
//         c0 = k  + l*(Nx-1);
//         c1 = k  + l*(Nx);
//         c2 = k  + l*(Nx+1);
        
//         if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
//             tyy      = -(mesh->sxxd[c0]  + mesh->szzd[c0]);
//             tyy0     = -(mesh->sxxd0[c0] + mesh->szzd0[c0]);
//             sxx0     = -mesh->p0_n[c0] + mesh->sxxd0[c0];
//             szz0     = -mesh->p0_n[c0] + mesh->szzd0[c0];
//             syy0     = -mesh->p0_n[c0] + tyy0;
//             sxx      = -mesh->p_in[c0] + mesh->sxxd[c0];
//             szz      = -mesh->p_in[c0] + mesh->szzd[c0];
//             syy      = -mesh->p_in[c0] + tyy;
//             dsxx[c0] = sxx - sxx0;
//             dsyy[c0] = syy - syy0;
//             dszz[c0] = szz - szz0;
//         }
//     }
    
//     // Vertex: shear stress change
//     for (k=0; k<Nx; k++) {
//         for (l=0; l<Nz; l++) {
//             c1 = k  + l*(Nx);
//             if (mesh->BCg.type[c1] !=30 ) dsxz[c1] =  mesh->sxz[c1] - mesh->sxz0[c1];
//         }
//     }
    
//     // Interpolate stress changes to markers
//     Interp_Grid2P_centroids2( *particles, mdsxx,  mesh, dsxx, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
//     Interp_Grid2P_centroids2( *particles, mdszz,  mesh, dszz, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
//     Interp_Grid2P_centroids2( *particles, mdsyy,  mesh, dsyy, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
//     Interp_Grid2P( *particles, mdsxz,  mesh, dsxz, mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );
    
//     // Update marker stresses
//     ArrayPlusArray( particles->sxxd, mdsxx, particles->Nb_part );
//     ArrayPlusArray( particles->szzd, mdszz, particles->Nb_part );
//     ArrayPlusArray( particles->syy , mdsyy, particles->Nb_part );
//     ArrayPlusArray( particles->sxz,  mdsxz, particles->Nb_part );
    
//     // Freedom
//     DoodzFree( dsxx );
//     DoodzFree( dsyy );
//     DoodzFree( dszz );
//     DoodzFree( dsxz );
//     DoodzFree(mdsxx);
//     DoodzFree(mdszz);
//     DoodzFree(mdsyy);
//     DoodzFree(mdsxz);
// }

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AccumulatedStrainII( grid* mesh, scale scaling, params model, markers* particles, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag ) {
    
    double *strain_inc;
    int k, l, c1;
    DoodzFP *strain_inc_el, *strain_inc_pl, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;
    
    //    printf("Accumulating strain\n");
    
    // Allocate
    strain_inc      = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_el   = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pl   = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_pwl  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_exp  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_lin  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    strain_inc_gbs  = DoodzMalloc(sizeof(double)*(mesh->Nx-1)*(mesh->Nz-1));
    
    // Interpolate exz to cell centers
#pragma omp parallel for shared( mesh, model, strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs  ) private( c1 )
    for (c1=0; c1<(mesh->Nx-1)*(mesh->Nz-1); c1++) {
        
        strain_inc[c1]     = model.dt*(mesh->eII_pl[c1]+mesh->eII_pwl[c1]+mesh->eII_exp[c1]+mesh->eII_lin[c1]+mesh->eII_gbs[c1]+mesh->eII_el[c1]+mesh->eII_cst[c1]);

        if (strain_inc[c1]<0) {
            printf("negative strain increment\n");
            printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n" ,mesh->eII_pl[c1],mesh->eII_pwl[c1],mesh->eII_exp[c1],mesh->eII_lin[c1],mesh->eII_gbs[c1],mesh->eII_el[c1]);
            exit(0);
        }
        
        strain_inc_el[c1]  = model.dt*mesh->eII_el[c1];
        strain_inc_pl[c1]  = model.dt*mesh->eII_pl[c1];
        strain_inc_pwl[c1] = model.dt*mesh->eII_pwl[c1];
        strain_inc_exp[c1] = model.dt*mesh->eII_exp[c1];
        strain_inc_lin[c1] = model.dt*mesh->eII_lin[c1];
        strain_inc_gbs[c1] = model.dt*mesh->eII_gbs[c1];
    }
    
    //---------------------------------------------------------------//
    
    double dE_tot, dE_el, dE_pl, dE_pwl, dE_exp, dE_lin, dE_gbs;
    double dx, dz, dxm, dzm, dst, sumW;
    int    i_part, j_part, iSW, iNW, iSE, iNE;
    dx=mesh->dx;
    dz=mesh->dz;
    
    //#pragma omp parallel for shared( particles , strain_inc, strain_inc_el, strain_inc_pl, strain_inc_pwl, strain_inc_exp, strain_inc_lin, strain_inc_gbs, tag  ) private( dE_tot, dE_el, dE_pl, dE_pwl, dE_exp, dE_lin, dE_gbs, sumW, dst, dxm, dzm, iSW, iSE, iNW, iNE, i_part, j_part  ) firstprivate (dx, dz, Nx, Nz)
    for (k=0;k<particles->Nb_part;k++) {
        
        dE_tot = 0.0;
        dE_el  = 0.0;
        dE_pl  = 0.0;
        dE_pwl = 0.0;
        dE_exp = 0.0;
        dE_lin = 0.0;
        dE_gbs = 0.0;
        
        
        if (particles->phase[k] != -1) {
            
            dst = (particles->x[k]-X_vect[0]);
            j_part   = ceil((dst/dx)) - 1;
            if (j_part<0) {
                j_part = 0;
            }
            if (j_part>Nx-2) {
                j_part = Nx-2;
            }
            
            // Get the line:
            dst = (particles->z[k]-Z_vect[0]);
            i_part   = ceil((dst/dz)) - 1;
            if (i_part<0) {
                i_part = 0;
            }
            if (i_part>Nz-2) {
                i_part = Nz-2;
            }
            
            dxm = (particles->x[k] - X_vect[j_part]);
            dzm = (particles->z[k] - Z_vect[i_part]);
            
            iSW = j_part+i_part*Nx;
            iSE = j_part+i_part*Nx+1;
            iNW = j_part+(i_part+1)*Nx;
            iNE = j_part+(i_part+1)*Nx+1;
            
            dE_tot = 0.0;
            dE_el  = 0.0;
            dE_pl  = 0.0;
            dE_pwl = 0.0;
            dE_exp = 0.0;
            dE_lin = 0.0;
            dE_gbs = 0.0;
            
            sumW         = 0.0;
            
            //            if (j_part>0 && j_part<Nx-2 && i_part>0 && i_part<Nz-2) {
            
            if (tag[iSW]!=30 && tag[iSW]!=31) {
                dE_tot +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc[iSW];
                dE_el  +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_el[iSW];
                dE_pl  +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pl[iSW];
                dE_pwl +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_pwl[iSW];
                dE_exp +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_exp[iSW];
                dE_lin +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_lin[iSW];
                dE_gbs +=  (1.0-dxm/dx) * (1.0-dzm/dz) * strain_inc_gbs[iSW];
                sumW   +=  (1.0-dxm/dx) * (1.0-dzm/dz);
            }
            if (tag[iSE]!=30 && tag[iSE]!=31) {
                dE_tot +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc[iSE];
                dE_el  +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_el[iSE];
                dE_pl  +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pl[iSE];
                dE_pwl +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_pwl[iSE];
                dE_exp +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_exp[iSE];
                dE_lin +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_lin[iSE];
                dE_gbs +=  (dxm/dx) * (1.0-dzm/dz) * strain_inc_gbs[iSE];
                sumW   +=  (dxm/dx) * (1.0-dzm/dz);
            }
            if (tag[iNW]!=30 && tag[iNW]!=31) {
                dE_tot +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc[iNW];
                dE_el  +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_el[iNW];
                dE_pl  +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pl[iNW];
                dE_pwl +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_pwl[iNW];
                dE_exp +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_exp[iNW];
                dE_lin +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_lin[iNW];
                dE_gbs +=  (1.0-dxm/dx) * (dzm/dz) * strain_inc_gbs[iNW];
                sumW   +=  (1.0-dxm/dx) * (dzm/dz);
            }
            if (tag[iNE]!=30 && tag[iNE]!=31) {
                dE_tot +=  (dxm/dx) * (dzm/dz) * strain_inc[iNE];
                dE_el  +=  (dxm/dx) * (dzm/dz) * strain_inc_el[iNE];
                dE_pl  +=  (dxm/dx) * (dzm/dz) * strain_inc_pl[iNE];
                dE_pwl +=  (dxm/dx) * (dzm/dz) * strain_inc_pwl[iNE];
                dE_exp +=  (dxm/dx) * (dzm/dz) * strain_inc_exp[iNE];
                dE_lin +=  (dxm/dx) * (dzm/dz) * strain_inc_lin[iNE];
                dE_gbs +=  (dxm/dx) * (dzm/dz) * strain_inc_gbs[iNE];
                sumW += (dxm/dx)* (dzm/dz);
            }
            
            //            printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n", sumW, dE_tot, strain_inc[iSW], strain_inc[iSE], strain_inc[iNW], strain_inc[iNE]);
            //
            //             if (dE_tot<0.0)exit(1);
            //            }
            
            //            if (j_part==0    && i_part==0   ) {
            //                dE_tot = strain_inc[iNE];
            //                dE_el  = strain_inc[iNE];
            //                dE_pl  = strain_inc[iNE];
            //                dE_pwl = strain_inc[iNE];
            //                dE_exp = strain_inc[iNE];
            //                dE_lin = strain_inc[iNE];
            //                dE_gbs = strain_inc[iNE];
            //            }
            //            if (j_part==Nx-2 && i_part==0   ) {
            //                dE_tot = strain_inc[iNW];
            //                dE_el  = strain_inc[iNW];
            //                dE_pl  = strain_inc[iNW];
            //                dE_pwl = strain_inc[iNW];
            //                dE_exp = strain_inc[iNW];
            //                dE_lin = strain_inc[iNW];
            //                dE_gbs = strain_inc[iNW];
            //            }
            //            if (j_part==0    && i_part==Nz-2) {
            //                dE_tot = strain_inc[iSE];
            //                dE_el  = strain_inc[iSE];
            //                dE_pl  = strain_inc[iSE];
            //                dE_pwl = strain_inc[iSE];
            //                dE_exp = strain_inc[iSE];
            //                dE_lin = strain_inc[iSE];
            //                dE_gbs = strain_inc[iSE];
            //            }
            //            if (j_part==Nx-2 && i_part==Nz-2) {
            //                dE_tot = strain_inc[iSW];
            //                dE_el  = strain_inc[iSW];
            //                dE_pl  = strain_inc[iSW];
            //                dE_pwl = strain_inc[iSW];
            //                dE_exp = strain_inc[iSW];
            //                dE_lin = strain_inc[iSW];
            //                dE_gbs = strain_inc[iSW];
            //            }
            
            //            if (j_part==0 && i_part>0 && i_part<Nz-2) {
            //
            //                if (tag[iSE]!=30 && tag[iSE]!=31) {
            //                    dE_tot += (1.0-dzm/dz)   * strain_inc[iSE];
            //                    dE_el  += (1.0-dzm/dz)   * strain_inc_el[iSE];
            //                    dE_pl  += (1.0-dzm/dz)   * strain_inc_pl[iSE];
            //                    dE_pwl += (1.0-dzm/dz)   * strain_inc_pwl[iSE];
            //                    dE_exp += (1.0-dzm/dz)   * strain_inc_exp[iSE];
            //                    dE_lin += (1.0-dzm/dz)   * strain_inc_lin[iSE];
            //                    dE_gbs += (1.0-dzm/dz)   * strain_inc_gbs[iSE];
            //                    sumW +=  (1.0-dzm/dz);
            //                }
            //                if (tag[iNE]!=30 && tag[iNE]!=31) {
            //                    dE_tot += (0.0-dzm/dz)   * strain_inc[iNE];
            //                    dE_el  += (0.0-dzm/dz)   * strain_inc_el[iNE];
            //                    dE_pl  += (0.0-dzm/dz)   * strain_inc_pl[iNE];
            //                    dE_pwl += (0.0-dzm/dz)   * strain_inc_pwl[iNE];
            //                    dE_exp += (0.0-dzm/dz)   * strain_inc_exp[iNE];
            //                    dE_lin += (0.0-dzm/dz)   * strain_inc_lin[iNE];
            //                    dE_gbs += (0.0-dzm/dz)   * strain_inc_gbs[iNE];
            //                    sumW +=  (dzm/dz);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }
            //
            //            if (j_part==Nx-2 && i_part>0 && i_part<Nz-2) {
            //
            //                if (tag[iSW]!=30 && tag[iSW]!=31) {
            //                    dE_tot += (1.0-dzm/dz)   * strain_inc[iSW];
            //                    dE_el  += (1.0-dzm/dz)   * strain_inc_el[iSW];
            //                    dE_pl  += (1.0-dzm/dz)   * strain_inc_pl[iSW];
            //                    dE_pwl += (1.0-dzm/dz)   * strain_inc_pwl[iSW];
            //                    dE_exp += (1.0-dzm/dz)   * strain_inc_exp[iSW];
            //                    dE_lin += (1.0-dzm/dz)   * strain_inc_lin[iSW];
            //                    dE_gbs += (1.0-dzm/dz)   * strain_inc_gbs[iSW];
            //                    sumW +=  (1.0-dzm/dz);
            //                }
            //                if (tag[iNW]!=30 && tag[iNW]!=31) {
            //                    dE_tot += (0.0-dzm/dz)   * strain_inc[iNW];
            //                    dE_el  += (0.0-dzm/dz)   * strain_inc_el[iNW];
            //                    dE_pl  += (0.0-dzm/dz)   * strain_inc_pl[iNW];
            //                    dE_pwl += (0.0-dzm/dz)   * strain_inc_pwl[iNW];
            //                    dE_exp += (0.0-dzm/dz)   * strain_inc_exp[iNW];
            //                    dE_lin += (0.0-dzm/dz)   * strain_inc_lin[iNW];
            //                    dE_gbs += (0.0-dzm/dz)   * strain_inc_gbs[iNW];
            //                    sumW +=  (dzm/dz);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }
            
            //            if (i_part==0 && j_part>0 && j_part<Nx-2) {
            //                if (tag[iNW]!=30 && tag[iNW]!=31) {
            //                    dE_tot += (1.0-dxm/dx)   * strain_inc[iNW];
            //                    dE_el  += (1.0-dxm/dx)   * strain_inc_el[iNW];
            //                    dE_pl  += (1.0-dxm/dx)   * strain_inc_pl[iNW];
            //                    dE_pwl += (1.0-dxm/dx)   * strain_inc_pwl[iNW];
            //                    dE_exp += (1.0-dxm/dx)   * strain_inc_exp[iNW];
            //                    dE_lin += (1.0-dxm/dx)   * strain_inc_lin[iNW];
            //                    dE_gbs += (1.0-dxm/dx)   * strain_inc_gbs[iNW];
            //                    sumW += (1.0-dxm/dx);
            //                }
            //                if (tag[iNE]!=30 && tag[iNE]!=31) {
            //                    dE_tot += (dxm/dx)   * strain_inc[iNE];
            //                    dE_el  += (dxm/dx)   * strain_inc_el[iNE];
            //                    dE_pl  += (dxm/dx)   * strain_inc_pl[iNE];
            //                    dE_pwl += (dxm/dx)   * strain_inc_pwl[iNE];
            //                    dE_exp += (dxm/dx)   * strain_inc_exp[iNE];
            //                    dE_lin += (dxm/dx)   * strain_inc_lin[iNE];
            //                    dE_gbs += (dxm/dx)   * strain_inc_gbs[iNE];
            //                    sumW += (dxm/dx);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }
            //            if (i_part==Nz-2 && j_part>0 && j_part<Nx-2) {
            //                if (tag[iSW]!=30 && tag[iSW]!=31) {
            //                    dE_tot += (1.0-dxm/dx)   * strain_inc[iSW];
            //                    dE_el  += (1.0-dxm/dx)   * strain_inc_el[iSW];
            //                    dE_pl  += (1.0-dxm/dx)   * strain_inc_pl[iSW];
            //                    dE_pwl += (1.0-dxm/dx)   * strain_inc_pwl[iSW];
            //                    dE_exp += (1.0-dxm/dx)   * strain_inc_exp[iSW];
            //                    dE_lin += (1.0-dxm/dx)   * strain_inc_lin[iSW];
            //                    dE_gbs += (1.0-dxm/dx)   * strain_inc_gbs[iSW];
            //                    sumW += (1.0-dxm/dx);
            //                }
            //                if ( tag[iSE] != 30 && tag[iSE] != 31 ) {
            //                    dE_tot += (dxm/dx)   * strain_inc[iSE];
            //                    dE_el  += (dxm/dx)   * strain_inc_el[iSE];
            //                    dE_pl  += (dxm/dx)   * strain_inc_pl[iSE];
            //                    dE_pwl += (dxm/dx)   * strain_inc_pwl[iSE];
            //                    dE_exp += (dxm/dx)   * strain_inc_exp[iSE];
            //                    dE_lin += (dxm/dx)   * strain_inc_lin[iSE];
            //                    dE_gbs += (dxm/dx)   * strain_inc_gbs[iSE];
            //                    sumW += (dxm/dx);
            //                }
            //                if(sumW>1e-13) {
            //                    dE_tot /= sumW;
            //                    dE_el  /= sumW;
            //                    dE_pl  /= sumW;
            //                    dE_pwl /= sumW;
            //                    dE_exp /= sumW;
            //                    dE_lin /= sumW;
            //                    dE_gbs /= sumW;
            //                }
            //            }
            
            if(sumW>1e-13) {
                dE_tot /= sumW;
                dE_el  /= sumW;
                dE_pl  /= sumW;
                dE_pwl /= sumW;
                dE_exp /= sumW;
                dE_lin /= sumW;
                dE_gbs /= sumW;
            }
            
            
        }
        
        particles->strain[k]      += dE_tot;
        particles->strain_el[k]   += dE_el;
        particles->strain_pl[k]   += dE_pl;
        particles->strain_pwl[k]  += dE_pwl;
        particles->strain_exp[k]  += dE_exp;
        particles->strain_lin[k]  += dE_lin;
        particles->strain_gbs[k]  += dE_gbs;
    }
    
    //---------------------------------------------------------------//
    
    // Free
    DoodzFree(strain_inc);
    DoodzFree(strain_inc_el);
    DoodzFree(strain_inc_pl);
    DoodzFree(strain_inc_pwl);
    DoodzFree(strain_inc_exp);
    DoodzFree(strain_inc_lin);
    DoodzFree(strain_inc_gbs);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DeformationGradient ( grid mesh, scale scaling, params model, markers *particles ) {
    
    int k, l, cp, cu, cv, Nx, Nz;
    double dx, dz;
    double fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o;
    
    Nx = mesh.Nx;
    Nz = mesh.Nz;
    dx = mesh.dx;
    dz = mesh.dz;
    
    double *dudx, *dudz, *dvdx, *dvdz;
    DoodzFP *pdudx, *pdudz, *pdvdx, *pdvdz;
    
    dudx   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
    dvdz   = DoodzMalloc ((Nx-1)*(Nz-1)*sizeof(double));
    dudz   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
    dvdx   = DoodzMalloc ((Nx)*(Nz)*sizeof(double));
    pdudx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdudz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdvdx  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    pdvdz  = DoodzMalloc (particles->Nb_part*sizeof(DoodzFP));
    
    // Compute dudx and dvdz (cell centers)
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            cp = k  + l*(Nx-1);
            cu = k  + l*(Nx) + Nx;
            cv = k  + l*(Nx+1) + 1;
            dudx[cp] = 1.0/dx * ( mesh.u_in[cu+1]      - mesh.u_in[cu] );
            dvdz[cp] = 1.0/dz * ( mesh.v_in[cv+(Nx+1)] - mesh.v_in[cv] );
        }
    }
    
    // Compute dudx and dvdz (cell vertices)
    for (k=0; k<Nx; k++) {
        for (l=0; l<Nz; l++) {
            cp = k  + l*(Nx-1);
            cu = k  + l*(Nx);
            cv = k  + l*(Nx+1);
            dudz[cu] = 1.0/dz * ( mesh.u_in[cu+Nx] - mesh.u_in[cu] );
            dvdx[cu] = 1.0/dx * ( mesh.v_in[cv+1]  - mesh.v_in[cv] );
        }
    }
    
    // Interpolate from grid to particles
    Interp_Grid2P_centroids2( *(particles), pdudx, &mesh, dudx, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model  );
    Interp_Grid2P_centroids2( *(particles), pdvdz, &mesh, dvdz, mesh.xvz_coord,  mesh.zvx_coord,  mesh.Nx-1, mesh.Nz-1, mesh.BCp.type, &model );
    Interp_Grid2P( *(particles), pdvdx, &mesh, dvdx, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    Interp_Grid2P( *(particles), pdudz, &mesh, dudz, mesh.xg_coord,  mesh.zg_coord,  mesh.Nx,   mesh.Nz, mesh.BCg.type   );
    
#pragma omp parallel for shared ( particles ) private ( k, fxx, fxz, fzx, fzz, fxx_o, fxz_o, fzx_o, fzz_o ) firstprivate( model ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {
        
        fxx = 1.0 + pdudx[k]*model.dt;
        fxz = pdudz[k]*model.dt;
        fzz = 1.0 + pdvdz[k]*model.dt;
        fzx = pdvdx[k]*model.dt;
        fxx_o = particles->Fxx[k];
        fxz_o = particles->Fxz[k];
        fzx_o = particles->Fzx[k];
        fzz_o = particles->Fzz[k];
        particles->Fxx[k] = fxx*fxx_o + fxz*fzx_o;
        particles->Fxz[k] = fxx*fxz_o + fxz*fzz_o;
        particles->Fzx[k] = fzx*fxx_o + fzz*fzx_o;
        particles->Fzz[k] = fzx*fxz_o + fzz*fzz_o;
    }
    
    // Free
    DoodzFree( dudx );
    DoodzFree( dvdz );
    DoodzFree( dudz );
    DoodzFree( dvdx );
    DoodzFree( pdudx );
    DoodzFree( pdudz );
    DoodzFree( pdvdx );
    DoodzFree( pdvdz );
    
    printf("-----> Deformation gradient tensor Updated\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FiniteStrainAspectRatio ( grid *mesh, scale scaling, params model, markers *particles ) {
    
    int k, l, cp, cu, cv;
    double *FS_AR, CGxx, CGxz, CGzx, CGzz;
    double Tr, Det, U0xx, U0xz, U0zx, U0zz, s, t, e1, e2;
    int    cent=1, vert=0, prop=1, interp=0;
    
    FS_AR = DoodzCalloc (particles->Nb_part, sizeof(double));
    
    //#pragma omp parallel for shared ( particles ) private ( k, CGxx, CGxz, CGzx, CGzz, Tr, Det, U0xx, U0xz, U0zx, U0zz, s, t, e1, e2  ) schedule( static )
    for (k=0; k<particles->Nb_part; k++) {
        
        // Cauchy-Green
        CGxx = particles->Fxx[k]*particles->Fxx[k] + particles->Fzx[k]*particles->Fzx[k];
        CGxz = particles->Fxx[k]*particles->Fxz[k] + particles->Fzx[k]*particles->Fzz[k];
        CGzx = particles->Fxx[k]*particles->Fxz[k] + particles->Fzx[k]*particles->Fzz[k];
        CGzz = particles->Fxz[k]*particles->Fxz[k] + particles->Fzz[k]*particles->Fzz[k];
        
        // U0 = Sqrt(CG)
        Tr  = CGxx + CGzz;
        Det = CGxx*CGzz - CGxz*CGzx;
        s   = sqrt(Det);
        t   = sqrt(Tr + 2.0*s);
        U0xx = 1/t*(CGxx + s);
        U0xz = 1/t*CGxz;
        U0zx = 1/t*CGzx;
        U0zz = 1/t*(CGzz + s);
        
        // eigs(U0)
        Tr   = U0xx + U0zz;
        Det  = U0xx*U0zz - U0xz*U0zx;
        e1   = Tr/2.0 + sqrt( Tr*Tr/4.0 - Det );
        e2   = Tr/2.0 - sqrt( Tr*Tr/4.0 - Det );
        
        // aspect ratio
        FS_AR[k] = e1/e2;
    }
    
    P2Mastah( &model, *particles, FS_AR, mesh, mesh->FS_AR_n,   mesh->BCp.type,  1, 0, interp, cent, model.interp_stencil);
    P2Mastah( &model, *particles, FS_AR, mesh, mesh->FS_AR_s,   mesh->BCg.type,  1, 0, interp, vert, model.interp_stencil);
    
    DoodzFree(FS_AR);
    
    printf("-----> Finite strain aspect ratio updated\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AccumulatedStrain( grid* mesh, scale scaling, params model, markers* particles ) {
    
    double *strain_inc;
    int Nx, Nz, k, l, c1;
    DoodzFP *strain_inc_mark, *strain_inc_el, *strain_inc_pl, *strain_inc_pwl, *strain_inc_exp, *strain_inc_lin, *strain_inc_gbs;
    
    Nx = mesh->Nx;
    Nz = mesh->Nz;
    
    // Allocate
    //    exz_n           = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc      = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_el   = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_pl   = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_pwl  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_exp  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_lin  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_gbs  = DoodzMalloc(sizeof(double)*(Nx-1)*(Nz-1));
    strain_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    
    // Interpolate exz to cell centers
    for (k=0; k<Nx-1; k++) {
        for (l=0; l<Nz-1; l++) {
            c1 = k  + l*(Nx-1);
            
            strain_inc[c1]     = model.dt*(mesh->eII_pl[c1]+mesh->eII_pwl[c1]+mesh->eII_exp[c1]+mesh->eII_lin[c1]+mesh->eII_gbs[c1]+mesh->eII_el[c1]+mesh->eII_cst[c1]);
            if (strain_inc[c1]<0) {
                printf("negative strain increment\n");
                printf("%2.2e %2.2e %2.2e %2.2e %2.2e %2.2e\n" ,mesh->eII_pl[c1],mesh->eII_pwl[c1],mesh->eII_exp[c1],mesh->eII_lin[c1],mesh->eII_gbs[c1],mesh->eII_el[c1]);
                exit(0);
            }
            
            strain_inc_el[c1]  = model.dt*mesh->eII_el[c1];
            strain_inc_pl[c1]  = model.dt*mesh->eII_pl[c1];
            strain_inc_pwl[c1] = model.dt*mesh->eII_pwl[c1];
            strain_inc_exp[c1] = model.dt*mesh->eII_exp[c1];
            strain_inc_lin[c1] = model.dt*mesh->eII_lin[c1];
            strain_inc_gbs[c1] = model.dt*mesh->eII_gbs[c1];
        }
    }
    
    // Interp increments to particles
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_el, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_el, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_pl, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_pl, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_pwl, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_pwl, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_exp, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_exp, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_lin, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_lin, strain_inc_mark, particles->Nb_part );
    
    Interp_Grid2P_strain( *particles, strain_inc_mark, mesh, strain_inc_gbs, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type );
    ArrayPlusArray( particles->strain_gbs, strain_inc_mark, particles->Nb_part );
    
    // Free
    //    DoodzFree(exz_n);
    DoodzFree(strain_inc);
    DoodzFree(strain_inc_el);
    DoodzFree(strain_inc_pl);
    DoodzFree(strain_inc_pwl);
    DoodzFree(strain_inc_exp);
    DoodzFree(strain_inc_lin);
    DoodzFree(strain_inc_gbs);
    DoodzFree(strain_inc_mark);
    
    printf("-----> Accumulated strain updated\n");
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateMaxPT ( scale scaling, params model, markers *particles ) {
    
    // Update maximum pressure and temperature on markers
    int k;
    
#pragma omp parallel for shared ( particles ) private ( k ) firstprivate( model ) schedule( static )
    for(k=0; k<particles->Nb_part; k++) {
        
        if (particles->T[k] > particles->Tmax[k] ) particles->Tmax[k] = particles->T[k];
        if (particles->P[k] > particles->Pmax[k] ) particles->Pmax[k] = particles->P[k];
        
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleDensity( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *rho_inc_mark, *rho_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    rho_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    rho_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));

    //----------------------

    //    Interp_Grid2P_centroids2( *particles, particles->rho, mesh, mesh->rho_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
    //----------------------

    for (k=0;k<Ncx*Ncz;k++) {
        rho_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) rho_inc_grid[k] = mesh->rho_n[k] - mesh->rho0_n[k];
    }
    
    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, rho_inc_mark, mesh, rho_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
    // Increment temperature on particles
    ArrayPlusArray( particles->rho, rho_inc_mark, particles->Nb_part );
    
    DoodzFree(rho_inc_grid);
    DoodzFree(rho_inc_mark);

    //----------------------

    // for (int k=0; k<particles->Nb_part; k++ ) {
    //     particles->rho[k] = EvaluateDensity( particles->phase[k], particles->T[k], particles->P[k], particles->X[k],  particles->phi[k], &model, materials );
    // }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePhi( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *phi_inc_mark, *phi_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    phi_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    phi_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
    for (k=0;k<Ncx*Ncz;k++) {
        phi_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) phi_inc_grid[k] = mesh->phi_n[k] - mesh->phi0_n[k];
    }
    
    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, phi_inc_mark, mesh, phi_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
    // Increment temperature on particles
    ArrayPlusArray( particles->phi, phi_inc_mark, particles->Nb_part );
    
    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->phi[k] <= 0.0) particles->phi[k] = 0.0;
        if (particles->phi[k] >= 1.0) particles->phi[k] = 1.0;
    }
    
    DoodzFree(phi_inc_grid);
    DoodzFree(phi_inc_mark);
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleX( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
//     DoodzFP *T_inc_mark, *Tm0, dtm, *dTms, *dTgr, *dTmr, *rho_part;
//     double *dT, *Tg0, *dTgs, dx=model.dx, dz=model.dz, d=1.0;
//     int    cent=1, vert=0, prop=1, interp=0;
//     int Nx, Nz, Ncx, Ncz, k, c0, p;
//     Nx = mesh->Nx; Ncx = Nx-1;
//     Nz = mesh->Nz; Ncz = Nz-1;
    
//     // Allocations
//     Tm0        = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     Tg0        = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     dT        = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
//     T_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    
//     // Old temperature grid
// #pragma omp parallel for shared(mesh, Tg0) private(c0) firstprivate(Ncx,Ncz)
//     for ( c0=0; c0<Ncx*Ncz; c0++ ) {
//         if (mesh->BCt.type[c0] != 30) { 
//             Tg0[c0] = mesh->X0_n[c0];
//             dT[c0]  = mesh->X_n[c0] - mesh->X0_n[c0];
//         }
//     }
//     Interp_Grid2P_centroids2( *particles, Tm0, mesh, Tg0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
 
//     /* CASE WITHOUT SUBGRID DIFFUSION: INTERPOLATE INCREMENT DIRECTLY */
        
//         // Interp increments to particles
//         Interp_Grid2P_centroids2( *particles, T_inc_mark, mesh, dT, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

//         // Increment temperature on particles
//         ArrayPlusArray( particles->X, T_inc_mark, particles->Nb_part );
    
    
//     // Freedom
//     DoodzFree(dT);
//     DoodzFree(Tg0);
//     DoodzFree(Tm0);
//     DoodzFree(T_inc_mark);

    //--------------------------------------------------
    DoodzFP *X_inc_mark, *X_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    // // Total update (do not use!!!)
    // Interp_Grid2P_centroids2( *particles, particles->X, mesh, mesh->X_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
    
    // Incremental update
    X_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    X_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
#pragma omp parallel for shared ( X_inc_grid, mesh ) private( k ) firstprivate( Ncx, Ncz )
    for (k=0; k<Ncx*Ncz; k++) {
        X_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) X_inc_grid[k] = mesh->X_n[k] - mesh->X0_n[k];
    }
    
    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, X_inc_mark, mesh, X_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
    
    // Increment X on particles
    ArrayPlusArray( particles->X, X_inc_mark, particles->Nb_part );
    
    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->X[k] <= 0.0) particles->X[k] = 0.0;
        if (particles->X[k] >= 1.0) particles->X[k] = 1.0;
    }
    
    DoodzFree(X_inc_grid);
    DoodzFree(X_inc_mark);
}

//--------------------------------------------------
void UpdateParticleXpips( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {

    DoodzFP *X_inc_mark, *X_inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    //==================================================================
    // Diffuse Xreac array ---------------------------------------------
    //==================================================================
    
    //double diff = model.Reac_diff/(scaling->L*scaling->L/scaling->t);
    //double diff = 1.0e-15/(scaling.L*scaling.L/scaling.t);
    double DW, DE, DS, DN;
    double diff_time = model.dt;
    double dx = model.dx;
    double dz = model.dz;
    
    double X_wanted, valmin, valmax, deltaval, valt_test, Xtest;
    double d=1.0;
    int nb_it, nb_itmax;
    int i,c,l, c1, c2, c3;
    
    //    get min max kc
    double min_kc=1.0e30, max_kc=0.0;
    for ( k=0; k<materials->Nb_phases; k++) {
        if (materials->k_chem[k]<min_kc) min_kc = materials->k_chem[k];
        if (materials->k_chem[k]>max_kc) max_kc = materials->k_chem[k];
    }
    printf("--> min kc = %2.2e m2.s-1 - max kc = %2.2e m2.s-1 \n", min_kc*scaling.L*scaling.L/scaling.t, max_kc*scaling.L*scaling.L/scaling.t );
    
    
    double *knodes,*newkx,*newkz;
    
    knodes  = DoodzCalloc(Nx*Nz, sizeof(DoodzFP));
    newkx   = DoodzCalloc(Nx*(Nz+1), sizeof(DoodzFP));
    newkz   = DoodzCalloc((Nx+1)*Nz, sizeof(DoodzFP));
    
    int interp=0, vert=0;

    P2Mastah( &model, *particles, materials->k_chem, mesh, knodes, mesh->BCg.type,  0, 0, interp, vert, (&model)->interp_stencil);

    // Interp_P2N ( *particles, materials->k_chem,  mesh, knodes, mesh->xg_coord,  mesh->zg_coord, 0, 0, &model );
    MinMaxArray(knodes, scaling.L*scaling.L/scaling.t, Nx*Nz, "knodes" );
    
    // fill newkx
    for( l=0; l<Nz+1; l++) {
        for( k=0; k<Nx; k++) {
            
            c1 = k + l*Nx;     // current value
            c3 = k + (l-1)*Nx; // value below
            
            if (l==0)          newkx[c1] = knodes[c1];
            if (l==Nz)         newkx[c1] = knodes[c3];
            if (l>0 && l<Nz)   newkx[c1] = 0.5*(knodes[c1] + knodes[c3]);
            
        }
    }
    MinMaxArray(newkx, scaling.L*scaling.L/scaling.t, Nx*(Nz+1), "newkx" );
    
    // fill newkz
    for( l=0; l<Nz; l++) {
        for( k=0; k<Nx+1; k++) {
            
            c1 = k   + l*(Nx+1);     // current value
            c2 = k   + l*(Nx);       // value on nodes
            c3 = k   + l*(Nx) - 1 ;     // value on nodes left
            
            //            if (k==0)          newkz[c1] = knodes[c2];
            //            if (k==Nx)         newkz[c1] = knodes[c3];
            
            //          if (k==0 && model.isperiodic_x == 1)    newkz[c1] = knodes[c2 + Nx-1];
            if (k==0 && model.periodic_x == 1)    newkz[c1] = 0.5*(knodes[c2] + knodes[c2+(Nx-2)]);
            if (k==0 && model.periodic_x == 0)    newkz[c1] = knodes[c2];
            
            //          if (k==Nx && model.isperiodic_x == 1)   newkz[c1] = knodes[c3 - Nx-1];
            if (k==Nx && model.periodic_x == 1)   newkz[c1] = 0.5*(knodes[c3 - (Nx-2)] + knodes[c3]);
            if (k==Nx && model.periodic_x == 0)   newkz[c1] = knodes[c3];
            
            if (k>0 && k<Nz)   newkz[c1] = 0.5*(knodes[c2] + knodes[c3]);
            
        }
    }
    MinMaxArray(newkz, scaling.L*scaling.L/scaling.t, (Nx+1)*Nz, "newkz" );
    
    
    
    double dt   = 0.249*model.dx*model.dz/max_kc, time=0.0;
    double dt_special = 0.0;
    
    
    printf("\n--------\n");
    printf("--> model time  = %f \n", model.dt );
    printf("--> dt_explicit = %f \n", dt );
    
    if (dt>=diff_time) dt = diff_time;
    
    int nstep   = floor(diff_time/dt);
    
    double vN, vS, vE, vW, vC;
    
    double *temp;
    
    temp  = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
    //    printf("\n--------\n");
    printf("--> dt_chosen   = %f \n", dt );
    //    printf("--> K           = %f \n", diff );
    printf("--> nstep       = %d \n", nstep );
    //    printf("--> K           = %2.2e m2.s-1 \n", diff*scaling.L*scaling.L/scaling.t );
    
    
    if (diff_time - (nstep*dt) > 0.0){
        dt_special = diff_time - (nstep*dt);
    }
    
    for (i=0; i<=nstep; i++) {
        
        ArrayEqualArray(temp, mesh->X_n,Ncx * Ncz);
        
        if (i==nstep){
            time += dt_special;
            dt = dt_special;
        }
        else{
            time += dt;
        }
        
        for( l=0; l<Ncz; l++) {
            for( k=0; k<Ncx; k++) {
                
                c = k + l*Ncx;
                
                c1  = k   + (l+1)*Nx;
                c3  = k   + l*(Nx+1)+1;
                
                //                DW    = mesh->kc_x[c1]         ;
                //                DE    = mesh->kc_x[c1 + 1]     ;
                //                DS    = mesh->kc_z[c3]         ;
                //                DN    = mesh->kc_z[c3 + Nx+1]  ;
                
                DW    = newkx[c1]         ;
                DE    = newkx[c1 + 1]     ;
                DS    = newkz[c3]         ;
                DN    = newkz[c3 + Nx+1]  ;
                
                
                vC = temp[c];
                
                if (k==0 && model.periodic_x == 1)     vW = temp[c + Ncx-1];
                if (k==0 && model.periodic_x == 0)     vW = temp[c]        ;
                if (k>0)                                 vW = temp[c - 1]    ;
                
                if (k==Ncx-1 && model.periodic_x == 1) vE = temp[c-Ncx+1]  ;
                if (k==Ncx-1 && model.periodic_x == 0) vE = temp[c]        ;
                if (k<Ncx-1)                             vE = temp[c+1]      ;
                
                if (l==0) vS = temp[c];
                if (l>0)  vS = temp[c - Ncx];
                
                if (l==Ncz-1) vN = temp[c];
                if (l<Ncz-1)  vN = temp[c + Ncx];
                
                //mesh->X_n[c] += dt*diff*((vE-2.0*vC+vW)/(model.dx*model.dx)+(vN-2.0*vC+vS)/(model.dz*model.dz));
                //mesh->X_n[c] += dt*1.0e-15/(scaling.L*scaling.L/scaling.t)*((vE-2.0*vC+vW)/(model.dx*model.dx)+(vN-2.0*vC+vS)/(model.dz*model.dz));
                
                mesh->X_n[c] += dt*( DE*(vE-vC)/(dx*dx) - DW*(vC-vW)/(dx*dx) + DN*(vN-vC)/(dz*dz) - DS*(vC-vS)/(dz*dz));
                
            }
        }
    }
    
    printf("--> end of loops, time = %2.2e s \n", time*scaling.t );
    printf("--> model time  = %2.2e s \n", model.dt*scaling.t );
    printf("\n--------\n");
    MinMaxArrayTag( mesh->X_n,             1.0, Ncx*Ncz, "Xreac_n ", mesh->BCp.type );
    
    free(temp);
    free(knodes);
    free(newkx);
    free(newkz);
    //==================================================================
    //==================================================================
    
    // // Total update (do not use!!!)
    //    Interp_Grid2P( *particles, particles->X, mesh, mesh->X_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
    // Incremental update
    X_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    X_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
#pragma omp parallel for shared ( X_inc_grid, mesh ) private( k ) firstprivate( Ncx, Ncz )
    for (k=0; k<Ncx*Ncz; k++) {
        X_inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) X_inc_grid[k] = mesh->X_n[k] - mesh->X0_n[k];
    }
    
    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, X_inc_mark, mesh, X_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
    
    // Increment temperature on particles
    ArrayPlusArray( particles->X, X_inc_mark, particles->Nb_part );
    
    // Bounds: numerical capote
#pragma omp parallel for shared ( particles ) private( k )
    for (k=0; k<particles->Nb_part; k++) {
        if (particles->X[k] <= 0.0) particles->X[k] = 0.0;
        if (particles->X[k] >= 1.0) particles->X[k] = 1.0;
    }
    
    DoodzFree(X_inc_grid);
    DoodzFree(X_inc_mark);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleGrainSize( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
//     DoodzFP *d_inc_mark, *d_inc_grid;
//     int Nx, Nz, Ncx, Ncz, k;
    const int Nx = mesh->Nx; //Ncx = Nx-1;
    const int Nz = mesh->Nz; //Ncz = Nz-1;
    
//     d_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
//     d_inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    
//     for (k=0;k<Ncx*Ncz;k++) {
//         d_inc_grid[k] =0.0;
//         if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) d_inc_grid[k] = mesh->d_n[k] - mesh->d0_n[k];
//     }
    
//     // Interp increments to particles
//     Interp_Grid2P_centroids2( *particles, d_inc_mark, mesh, d_inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

//     MinMaxArrayPart( particles->d, scaling.L, particles->Nb_part, "d on markers", particles->phase ) ;

//     // Increment grain size on particles
//     ArrayPlusArray( particles->d, d_inc_mark, particles->Nb_part );

//     MinMaxArrayPart( particles->d, scaling.L, particles->Nb_part, "d on markers", particles->phase ) ;


// #pragma omp parallel for shared ( particles )
//     for (k=0;k<particles->Nb_part;k++) {
//         if (particles->d[k]<1.0e-16) particles->d[k] = 1.0e-16;
//     }
    
//     DoodzFree(d_inc_mark);
//     DoodzFree(d_inc_grid);

    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, particles->d, mesh, mesh->d_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleEnergy( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *T_inc_mark, *Tm0, dtm, *dTms, *dTgr, *dTmr, *rho_part;
    double *dT, *dTgs, dx=model.dx, dz=model.dz, d=1.0;
    int    cent=1, vert=0, prop=1, interp=0;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    // Allocations
    Tm0        = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    dT         = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
    T_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));

    // Old temperature grid
#pragma omp parallel for shared(mesh, dT) private(c0) firstprivate(Ncx,Ncz)
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        if (mesh->BCt.type[c0] != 30) {
            dT[c0]  = mesh->T[c0] - mesh->T0_n[c0];
        }
    }
    Interp_Grid2P_centroids2( *particles, Tm0, mesh, mesh->T0_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
    // SUBGRID
    if ( model.subgrid_diffusion >= 1 ) { /* CASE WITH SUBGRID DIFFUSION */
        
        printf("Subgrid diffusion for temperature update\n");
        dTgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        dTms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        dTmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        rho_part = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        
        // Get density
        Interp_Grid2P_centroids2( *particles, rho_part, mesh, mesh->rho_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
        
        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Tm0,dTms) private(k,p,dtm) firstprivate(materials,dx,dz,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p = particles->phase[k];
                dtm     = materials->Cp[p] * rho_part[k]/ (materials->k[p] * (1.0/dx/dx + 1.0/dz/dz));
                dTms[k] = -( particles->T[k] - Tm0[k] ) * (1.0-exp(-d*model.dt/dtm));
            }
        }
        
        // Subgrid temperature increments markers --> grid
        P2Mastah( &model, *particles, dTms,     mesh, dTgs,   mesh->BCp.type,  1, 0, interp, cent, 1);
        
        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dTgs, dTgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) dTgr[c0] = dT[c0] - dTgs[c0];
        }
        
        // Remaining temperature increments grid --> markers
        Interp_Grid2P_centroids2( *particles, dTmr, mesh, dTgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
        
        // Final temperature update on markers
#pragma omp parallel for shared(particles,dTms,dTmr,T_inc_mark) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) T_inc_mark[k]    = dTms[k] + dTmr[k];
            if (particles->phase[k] != -1) particles->T[k]  = particles->T[k] + T_inc_mark[k];
        }
        DoodzFree(dTms);
        DoodzFree(dTmr);
        DoodzFree(dTgs);
        DoodzFree(dTgr);
        DoodzFree(rho_part);
    }
    else {  /* CASE WITHOUT SUBGRID DIFFUSION: INTERPOLATE INCREMENT DIRECTLY */
        
        // Interp increments to particles
        Interp_Grid2P_centroids2( *particles, T_inc_mark, mesh, dT, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );

        // Increment temperature on particles
        ArrayPlusArray( particles->T, T_inc_mark, particles->Nb_part );
    }
    
    // Freedom
    DoodzFree(dT);
    DoodzFree(Tm0);
    DoodzFree(T_inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticlePressure( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *P_inc_mark;
    int Nx, Nz, Ncx, Ncz, k, c0, p;
    double d=1.0, dtm;
    int    cent=1, vert=0, prop=1, interp=0;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    // Compute increment
#pragma omp parallel for shared(mesh) private(c0) 
    for ( c0=0; c0<Ncx*Ncz; c0++ ) {
        mesh->dp[c0] = 0.0;
        if (mesh->BCp.type[c0] != 30 ) {
            mesh->dp[c0] = (mesh->p_in[c0]-mesh->p0_n[c0]);
        }
    }
    
    if ( model.subgrid_diffusion >= 2 ) {
        
        printf("Subgrid diffusion for pressure update\n");
        double *Pg0  = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *dPgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *dPgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
        double *Pm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *dPms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *dPmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        double *etam = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        
        Interp_Grid2P_centroids2( *particles, etam,  mesh, mesh->eta_phys_n, mesh->xvz_coord,  mesh->zvx_coord, Nx-1  , Nz-1 , mesh->BCp.type, &model);
        
        /* -------------- */
        // Old Pressure grid
#pragma omp parallel for shared(mesh, Pg0) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCt.type[c0] != 30) Pg0[c0] = mesh->p0_n[c0];
        }
        Interp_Grid2P_centroids2( *particles, Pm0, mesh, Pg0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        /* -------------- */
        
        //        Interp_Grid2P_centroids( *particles, Pm0, mesh, mesh->p0_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        
        
        
        // Old temperature grid --> markers
        //        Interp_Grid2P_centroids( *particles, Pm0, mesh, mesh->p0_n, mesh->xc_coord,  mesh->zc_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        
        // Compute subgrid temperature increments on markers
#pragma omp parallel for shared(particles,Pm0,dPms) private(k,p,dtm) firstprivate(materials,model,d)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                p         = particles->phase[k];
                dtm       = etam[k]/ (materials->bet[p]);
                dPms[k]   = -( particles->P[k] - Pm0[k]) * (1.0 - exp(-d*model.dt/dtm));
            }
        }
        
        // Subgrid temperature increments markers --> grid
        P2Mastah( &model, *particles, dPms,     mesh, dPgs,   mesh->BCp.type,  1, 0, interp, cent, 1);
        
        // Remaining temperature increments on the grid
#pragma omp parallel for shared(mesh, dPgs, dPgr) private(c0) firstprivate(Ncx,Ncz)
        for ( c0=0; c0<Ncx*Ncz; c0++ ) {
            if (mesh->BCp.type[c0] != 30 ) dPgr[c0] = mesh->dp[c0] - dPgs[c0];
        }
        
        // Remaining temperature increments grid --> markers
        Interp_Grid2P_centroids2( *particles, dPmr, mesh, dPgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        
        // Final temperature update on markers
#pragma omp parallel for shared(particles,dPms,dPmr) private(k)
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) particles->P[k] = particles->P[k] + dPms[k] + dPmr[k];
        }
        
        DoodzFree(Pg0);
        DoodzFree(Pm0);
        DoodzFree(dPms);
        DoodzFree(dPmr);
        DoodzFree(dPgs);
        DoodzFree(dPgr);
        DoodzFree(etam);
    }
    else {
        
        P_inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
        
        // Interp increments to particles
        Interp_Grid2P_centroids2( *particles, P_inc_mark, mesh, mesh->dp, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, &model  );
        
        // Increment pressure on particles
        ArrayPlusArray( particles->P, P_inc_mark, particles->Nb_part );
        
        DoodzFree(P_inc_mark);
        
    }
    
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleStress( grid* mesh, markers* particles, params* model, mat_prop* materials, scale *scaling ) {
    
    int k, l, c0, c1, c2, Nx, Nz, Ncx, Ncz, p, k1, c3;
    DoodzFP *mdsxxd, *mdszzd, *mdsxz, *dsxxd, *dszzd, *dsxz, d=1.0, dtaum;
    DoodzFP *dtxxgs, *dtzzgs, *dtxzgs, *dtxxgr, *dtzzgr, *dtxzgr, *txxm0, *tzzm0, *txzm0, *dtxxms, *dtzzms, *dtxzms, *dtxxmr, *dtzzmr, *dtxzmr, *etam;
    double *dudx_n, *dvdz_n, *dudz_s, *dvdx_s, *om_s, *om_n, *dudz_n, *dvdx_n, *dudx_s, *dvdz_s;
    double angle, tzz, txx, txz, dx, dz, dt;
    double *txz_n, *txx_s, *tzz_s, *dtxxg0, *dtzzg0, *dtxzg0;
    int    cent=1, vert=0, prop=1, interp=0;
    
    Nx = model->Nx;
    Nz = model->Nz;
    dx = model->dx;
    dz = model->dz;
    dt = model->dt;
    
    dudx_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dvdz_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dudz_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    dvdx_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    dudz_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dvdx_n = DoodzCalloc ((Nx-1)*(Nz-1), sizeof(double));
    dudx_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    dvdz_s = DoodzCalloc ((Nx-0)*(Nz-0), sizeof(double));
    
#pragma omp parallel for shared ( mesh, om_s, dudz_s, dvdx_s ) \
private ( k, l, k1, c1, c3 )                                   \
firstprivate( model )
    for ( k1=0; k1<Nx*Nz; k1++ ) {
        k  = mesh->kn[k1];
        l  = mesh->ln[k1];
        c1 = k + l*Nx;
        c3 = k + l*(Nx+1);
        if ( mesh->BCg.type[c1] != 30 ) {
            dudz_s[c1] = (mesh->u_in[c1+Nx] - mesh->u_in[c1])/dz;
            dvdx_s[c1] = (mesh->v_in[c3+1 ] - mesh->v_in[c3])/dx;
        }
    }
    
#pragma omp parallel for shared ( mesh, dudx_n, dvdz_n ) \
private ( k, l, k1, c0, c1, c2 )                         \
firstprivate( model )
    for ( k1=0; k1<(Nx-1)*(Nz-1); k1++ ) {
        k  = mesh->kp[k1];
        l  = mesh->lp[k1];
        c0 = k  + l*(Nx-1);
        c1 = k  + l*(Nx);
        c2 = k  + l*(Nx+1);
        if ( mesh->BCp.type[c0] != 30 && mesh->BCp.type[c0] != 31) {
            dudx_n[c0]  = (mesh->u_in[c1+1+Nx]     - mesh->u_in[c1+Nx] )/dx;
            dvdz_n[c0]  = (mesh->v_in[c2+1+(Nx+1)] - mesh->v_in[c2+1]  )/dz;
        }
    }
    
    InterpCentroidsToVerticesDouble( dudx_n, dudx_s, mesh, model );
    InterpCentroidsToVerticesDouble( dvdz_n, dvdz_s, mesh, model );
    InterpVerticesToCentroidsDouble( dudz_n, dudz_s, mesh, model );
    InterpVerticesToCentroidsDouble( dvdx_n, dvdx_s, mesh, model );

// Rotate director directly on particles
    if ( model->anisotropy == 1 && model->advection==1) {

#pragma omp parallel for shared( particles, mesh ) firstprivate( dt, model ) private( k )
        for ( k=0; k<particles->Nb_part; k++ ) {
            if (particles->phase[k] != -1) {
                double nx = particles->nx[k];
                double nz = particles->nz[k];
                double mdudx =  Centers2Particle( particles, dudx_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, k, model->periodic_x );
                double mdudz = Vertices2Particle( particles, dudz_s,     mesh->xg_coord,  mesh->zg_coord,  mesh->Nx-0, mesh->Nz-0, mesh->BCg.type, mesh->dx, mesh->dz, k );
                double mdvdz =  Centers2Particle( particles, dvdz_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, k, model->periodic_x );
                double mdvdx = Vertices2Particle( particles, dvdx_s,     mesh->xg_coord,  mesh->zg_coord,  mesh->Nx-0, mesh->Nz-0, mesh->BCg.type, mesh->dx, mesh->dz, k );
                particles->nx[k] += dt*(-(mdudx - mdvdz)*nx*nz - mdvdx*nz*nz + mdudz*nx*nx)*nz;
                particles->nz[k] += dt*( (mdudx - mdvdz)*nx*nz + mdvdx*nz*nz - mdudz*nx*nx)*nx;
                double norm = sqrt( pow(particles->nx[k],2) + pow(particles->nz[k],2));
                particles->nx[k] /= norm;
                particles->nz[k] /= norm;
            }
        }
    }

    // Marker stress update 
    if ( model->elastic==1 ) {
        
        Nx = mesh->Nx; Ncx = Nx-1;
        Nz = mesh->Nz; Ncz = Nz-1;
        
        if ( model->subgrid_diffusion == 2 ) {
            
            // Alloc
            dtxxgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtzzgs = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtxzgs = DoodzCalloc(Nx*Nz,   sizeof(DoodzFP));
            dtxxgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtzzgr = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));
            dtxzgr = DoodzCalloc(Nx*Nz,   sizeof(DoodzFP));
            txxm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            tzzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            txzm0  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxxms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtzzms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxzms = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxxmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtzzmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            dtxzmr = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            etam   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            
            // Old stresses grid --> markers
            Interp_Grid2P_centroids2( *particles, txxm0, mesh, mesh->sxxd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P_centroids2( *particles, tzzm0, mesh, mesh->szzd0, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P( *particles, txzm0, mesh, mesh->sxz0,       mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);
            Interp_Grid2P( *particles, etam,  mesh, mesh->eta_phys_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);
        

            MinMaxArray(mesh->eta_phys_s, scaling->eta, Nx*Nz, "eta grid  ");
            MinMaxArrayTag( mesh->eta_phys_s, scaling->eta, Nx*Nz,   "eta_s", mesh->BCg.type );
            MinMaxArray(etam, scaling->eta, particles->Nb_part, "eta phys part  ");

            // Interp_Grid2P_eta( *particles, etam,  mesh, mesh->eta_phys_s, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type);

            
            if ( model->subgrid_diffusion == 2 ) {
                
                printf("Subgrid diffusion for stress tensor component update\n");
                
                // Compute subgrid stress increments on markers
    #pragma omp parallel for shared(particles,txxm0,tzzm0,txzm0,dtxxms,dtzzms,dtxzms,etam) private(k,p,dtaum) firstprivate(materials,model,d)
                for ( k=0; k<particles->Nb_part; k++ ) {
                    if (particles->phase[k] != -1) {
                        p         = particles->phase[k];
                        dtaum     = etam[k] / materials->G[p];
                        dtxxms[k] = -( particles->sxxd[k] - txxm0[k]) * (1.0 - exp(-d*dt/dtaum));
                        dtzzms[k] = -( particles->szzd[k] - tzzm0[k]) * (1.0 - exp(-d*dt/dtaum));
                        dtxzms[k] = -( particles->sxz[k]  - txzm0[k]) * (1.0 - exp(-d*dt/dtaum));
                        if (isinf(dtxxms[k])) {
                            printf("Inf dtxxms[k]: sxxd=%2.2e txxm0=%2.2e exp(-d*dt/dtaum)=%2.2e d=%2.2e dt=%2.2e dtaum=%2.2e etam=%2.2e G=%2.2e\n", particles->sxxd[k], txxm0[k], exp(-d*dt/dtaum), d, dt, dtaum, etam[k]*scaling->eta, materials->G[p]*scaling->S);
                            printf("%2.2e %2.2e %2.2e %2.2e %2.2e", d, dt, dtaum, etam[k]*scaling->eta, materials->G[p]*scaling->S );
                            exit(1);
                        }
                        if (isnan(dtxxms[k])) {
                            printf("NaN dtxxms[k]: sxxd=%2.2e txxm0=%2.2e exp(-d*dt/dtaum)=%2.2e d=%2.2e dt=%2.2e dtaum=%2.2e etam=%2.2e G=%2.2e\n", particles->sxxd[k], txxm0[k], exp(-d*dt/dtaum), d, dt, dtaum, etam[k]*scaling->eta, materials->G[p]*scaling->S);
                            exit(1);
                        }
                    }
                }
                
                // Subgrid stress increments markers --> grid
                P2Mastah( model, *particles, dtxxms,     mesh, dtxxgs,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil);
                P2Mastah( model, *particles, dtzzms,     mesh, dtzzgs,   mesh->BCp.type,  1, 0, interp, cent, model->interp_stencil);
                P2Mastah( model, *particles, dtxzms,     mesh, dtxzgs,   mesh->BCg.type,  1, 0, interp, vert, model->interp_stencil);
                
                // Remaining stress increments on the grid
    #pragma omp parallel for shared(mesh,dtxxgs,dtxxgr,dtzzgs,dtzzgr) private(c0) firstprivate(Ncx,Ncz)
                for ( c0=0; c0<Ncx*Ncz; c0++ ) {
                    if (mesh->BCp.type[c0]!=30 && mesh->BCp.type[c0]!=31) dtxxgr[c0] = (mesh->sxxd[c0] - mesh->sxxd0[c0]) - dtxxgs[c0];
                    if (mesh->BCp.type[c0]!=30 && mesh->BCp.type[c0]!=31) dtzzgr[c0] = (mesh->szzd[c0] - mesh->szzd0[c0]) - dtzzgs[c0];
                }
    #pragma omp parallel for shared(mesh,dtxzgs,dtxzgr) private(c0) firstprivate(Nx,Nz)
                for ( c0=0; c0<Nx*Nz; c0++ ) {
                    if (mesh->BCg.type[c0]!=30) dtxzgr[c0] = (mesh->sxz[c0]-mesh->sxz0[c0]) - dtxzgs[c0];
                }
                
                // Remaining stress increments grid --> markers
                Interp_Grid2P_centroids2( *particles, dtxxmr, mesh, dtxxgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
                Interp_Grid2P_centroids2( *particles, dtzzmr, mesh, dtzzgr, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCp.type, model  );
                Interp_Grid2P(           *particles, dtxzmr, mesh, dtxzgr, mesh->xg_coord,  mesh->zg_coord, Nx  , Nz  , mesh->BCg.type         );
                
                // Final stresses update on markers
    #pragma omp parallel for shared(particles,dtxxms,dtzzms,dtxzms,dtxxmr,dtzzmr,dtxzmr) private(k) 
                for ( k=0; k<particles->Nb_part; k++ ) {
                    if (particles->phase[k] != -1) particles->sxxd[k]  = particles->sxxd[k] + dtxxms[k] + dtxxmr[k];
                    if (particles->phase[k] != -1) particles->szzd[k]  = particles->szzd[k] + dtzzms[k] + dtzzmr[k];
                    if (particles->phase[k] != -1) particles->sxz[k]   = particles->sxz[k]  + dtxzms[k] + dtxzmr[k];
                }
            }
            
            // Free
            DoodzFree( dtxxgs );
            DoodzFree( dtzzgs );
            DoodzFree( dtxzgs );
            DoodzFree( dtxxgr );
            DoodzFree( dtzzgr );
            DoodzFree( dtxzgr );
            DoodzFree( txxm0  );
            DoodzFree( tzzm0  );
            DoodzFree( txzm0  );
            DoodzFree( dtxxms );
            DoodzFree( dtzzms );
            DoodzFree( dtxzms );
            DoodzFree( dtxxmr );
            DoodzFree( dtzzmr );
            DoodzFree( dtxzmr );
            DoodzFree( etam   );  
        }

        if (model->subgrid_diffusion==0 || model->subgrid_diffusion==1){
            
            printf("No subgrid diffusion for stress tensor component update\n");
            
            // Alloc
            dsxxd  = DoodzCalloc((Nx-1)*(Nz-1), sizeof(double));
            dszzd  = DoodzCalloc((Nx-1)*(Nz-1), sizeof(double));
            dsxz   = DoodzCalloc((Nx)*(Nz), sizeof(double));
            mdsxxd = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            mdszzd = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            mdsxz  = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            
            // Cell: normal stress change
            for (k=0; k<Nx-1; k++) {
                for (l=0; l<Nz-1; l++) {
                    c0 = k  + l*(Nx-1);
                    if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dsxxd[c0] = mesh->sxxd[c0] - mesh->sxxd0[c0];
                    if (mesh->BCp.type[c0] !=30 && mesh->BCp.type[c0] !=31) dszzd[c0] = mesh->szzd[c0] - mesh->szzd0[c0];
                }
            }
            
            // Vertex: shear stress change
            for (k=0; k<Nx; k++) {
                for (l=0; l<Nz; l++) {
                    c1 = k  + l*(Nx);
                    if (mesh->BCg.type[c1] !=30 ) dsxz[c1] =  mesh->sxz[c1] - mesh->sxz0[c1];
                }
            }
            
            // Interpolate stress changes to markers
            Interp_Grid2P_centroids2( *particles, mdsxxd, mesh, dsxxd, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P_centroids2( *particles, mdszzd, mesh, dszzd, mesh->xvz_coord,  mesh->zvx_coord,  mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, model  );
            Interp_Grid2P( *particles, mdsxz,  mesh, dsxz,  mesh->xg_coord,  mesh->zg_coord,  mesh->Nx,   mesh->Nz,   mesh->BCg.type  );
            
            // Update marker stresses
            ArrayPlusArray(  particles->sxxd, mdsxxd, particles->Nb_part );
            ArrayPlusArray(  particles->szzd, mdszzd, particles->Nb_part );
            ArrayPlusArray(  particles->sxz,  mdsxz,  particles->Nb_part );

            // Free
            DoodzFree(dsxxd);
            DoodzFree(dszzd);
            DoodzFree(dsxz);
            DoodzFree(mdsxxd);
            DoodzFree(mdszzd);
            DoodzFree(mdsxz);
        }

        // Large strain: stress rotation/deformation
        if (model->elastic==1) { 

            double *wxzm   = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
            P2Mastah( model, *particles, wxzm,     mesh, mesh->wxz,   mesh->BCg.type,  1, 0, interp, vert, model->interp_stencil);

#pragma omp parallel for shared( particles, wxzm ) firstprivate( dt, model ) private( k, angle, txx, tzz, txz )
            for ( k=0; k<particles->Nb_part; k++ ) {
                if (particles->phase[k] != -1) {
                    double txx = particles->sxxd[k];
                    double tzz = particles->szzd[k];
                    double txz = particles->sxz[k];
                    if (model->stress_rotation==1) { // Jaumann rate
                        particles->sxxd[k] -=  dt*2.0*txz*wxzm[k];
                        particles->szzd[k] -= -dt*2.0*txz*wxzm[k];
                        particles->sxz[k]  -=  dt*(tzz-txx)*wxzm[k];
                    }
                    if (model->stress_rotation==2) { // Analytical rotation
                        double angle = dt*wxzm[k];
                        particles->sxxd[k] =  ( txx*cos(angle) + txz*sin(angle))*cos(angle) + ( txz*cos(angle) + tzz*sin(angle))*sin(angle);
                        particles->szzd[k] = -(-txx*sin(angle) + txz*cos(angle))*sin(angle) + (-txz*sin(angle) + tzz*cos(angle))*cos(angle);
                        particles->sxz[k]  = -( txx*cos(angle) + txz*sin(angle))*sin(angle) + ( txz*cos(angle) + tzz*sin(angle))*cos(angle);
                    }
                }
            }      
            DoodzFree(wxzm);
        }
    }

    DoodzFree(dudx_n);
    DoodzFree(dvdz_n);
    DoodzFree(dvdx_s);
    DoodzFree(dudz_s);
    DoodzFree(dudz_n);
    DoodzFree(dvdx_n);
    DoodzFree(dvdz_s);
    DoodzFree(dudx_s);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void UpdateParticleDivThermal( grid* mesh, scale scaling, params model, markers* particles, mat_prop* materials ) {
    
    DoodzFP *inc_mark, *inc_grid;
    int Nx, Nz, Ncx, Ncz, k;
    Nx = mesh->Nx; Ncx = Nx-1;
    Nz = mesh->Nz; Ncz = Nz-1;
    
    inc_mark = DoodzCalloc(particles->Nb_part, sizeof(DoodzFP));
    inc_grid = DoodzCalloc(Ncx*Ncz, sizeof(DoodzFP));    
    
    for (k=0;k<Ncx*Ncz;k++) {
        inc_grid[k] = 0.0;
        if (mesh->BCp.type[k] != 30 && mesh->BCp.type[k] != 31) inc_grid[k] = mesh->divth_n[k] - mesh->divth0_n[k];
    }
    
    // Interp increments to particles
    Interp_Grid2P_centroids2( *particles, inc_mark, mesh, inc_grid, mesh->xvz_coord,  mesh->zvx_coord, Nx-1, Nz-1, mesh->BCt.type, &model  );
    
    // Increment temperature on particles
    ArrayPlusArray( particles->divth, inc_mark, particles->Nb_part );
    
    DoodzFree(inc_grid);
    DoodzFree(inc_mark);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/