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

#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "mdoodz-private.h"
#include "time.h"

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num()  0
#define omp_get_num_threads() 1
#define omp_get_wtime() clock()/CLOCKS_PER_SEC
#endif

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Compute4Cell( int ic, int jc , int Nx, double** xc, grid *mesh, markers *particles, int k, int**npc, double** npvolc, double*** phc, int ith ) {
    double dxm, dzm, w;
    int ind;
    ind = ic + jc * Nx;
    dxm = fabs(        xc[ith][ic] - particles->x[k]);
    dzm = fabs( mesh->zc_coord[jc] - particles->z[k]);
    w   = (1.0-dxm/mesh->dx)*(1.0-dzm/mesh->dz);
    npc[ith][ind]++;
    npvolc[ith][ind]+=w;
    phc[ith][particles->phase[k]][ind]+=w;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void V2P( double *Vxm, double *Vzm, markers* particles, double*Vx,  double*Vz, double*xg, double*zg, double *zvx, double *xvz, int Nxg, int Nzg, int NzVx, int NxVz, char *tagVx, char *tagVz, double dx, double dz, int k, int TarasCorr ) {
    
    double val  = 0.0;
    double sumW = 0.0;
    double dxm, dzm, dst, vxm13, vxm24, vzm12, vzm34;
    int    i_part, j_part, iSW, iNW, iSE, iNE;
    int    Nx, Nz, iSEE, iNEE, iSWW, iNWW, iNNW, iNNE, iSSW, iSSE;
    double VSEE, VSWW, VNEE, VNWW, VSSW, VSSE, VNNW, VNNE, VSW, VSE, VNW, VNE;
    
    (*Vxm) = 0.0;
    (*Vzm) = 0.0;
    
    // Filter out particles that are inactive (out of the box)
    if (particles->phase[k] != -1) {
        
        Nx = Nxg; Nz = NzVx;
        
        dst    = (particles->x[k]-xg[0]);
        j_part = ceil((dst/dx)) - 1;
        if (j_part<0) {
            j_part = 0;
        }
        if (j_part>Nx-2) {
            j_part = Nx-2;
        }
        
        // Get the line:
        dst    = (particles->z[k]-zvx[0]);
        i_part = ceil((dst/dz)) - 1;
        if (i_part<0) {
            i_part = 0;
        }
        if (i_part>Nz-2) {
            i_part = Nz-2;
        }
        
        dxm = fabs(particles->x[k] - xg[j_part]);
        dzm = fabs(particles->z[k] - zvx[i_part]);
        
        iSW = j_part+i_part*Nx;
        iSE = j_part+i_part*Nx+1;
        iNW = j_part+(i_part+1)*Nx;
        iNE = j_part+(i_part+1)*Nx+1;
        
        VNW = Vx[iNW];
        VNE = Vx[iNE];
        VSW = Vx[iSW];
        VSE = Vx[iSE];
        
        if (tagVx[iSW] == 30) VSW = Vx[iNW];
        if (tagVx[iSE] == 30) VSE = Vx[iNE];
        if (tagVx[iNW] == 30) VNW = Vx[iSW];
        if (tagVx[iNE] == 30) VNE = Vx[iSE];

        // if (tagVx[iNE] == 30 && tagVx[iSE] == 30 && tagVx[iSW] == 30) exit(1130);
        // if (tagVx[iNW] == 30 && tagVx[iSW] == 30 && tagVx[iSE] == 30) exit(1131);
        
        if (tagVx[iNE] == 30 && tagVx[iSE] == 30) {
            VNE = VSW;
            VSE = VSW;
        }
        
        if (tagVx[iNW] == 30 && tagVx[iSW] == 30) {
            VNW = VSE;
            VSW = VSE;
        }
        
        vxm13 = VSW*(1.0-dxm/dx) + VSE*dxm/dx;
        vxm24 = VNW*(1.0-dxm/dx) + VNE*dxm/dx;
        
        // Compute correction
        if (TarasCorr==1) {
            if (dxm/dx>=0.5) {
                if(j_part<Nx-2) {
                    iSEE  = j_part+i_part*Nx+1 + 1;
                    iNEE  = j_part+(i_part+1)*Nx+1 + 1;
                    VSEE  = Vx[iSEE];
                    VNEE  = Vx[iNEE];
                    if (tagVx[iSEE] == 30) VSEE = Vx[iSW];
                    if (tagVx[iNEE] == 30) VNEE = Vx[iNW];
                    vxm13 = vxm13 + 0.5*pow((dxm/dx-0.5),2.0) * (VSW - 2.0*VSE + VSEE);
                    vxm24 = vxm24 + 0.5*pow((dxm/dx-0.5),2.0) * (VNW - 2.0*VNE + VNEE);
                }
            }
            else {
                if(j_part>0) {
                    iSWW  = j_part+i_part*Nx - 1;
                    iNWW  = j_part+(i_part+1)*Nx - 1;
                    VSWW  = Vx[iSWW];
                    VNWW  = Vx[iNWW];
                    if (tagVx[iSWW] == 30) VSWW = Vx[iSE];
                    if (tagVx[iNWW] == 30) VNWW = Vx[iNE];
                    vxm13 = vxm13 + 0.5*pow((dxm/dx-0.5),2.0) * (VSWW   -2.0*VSW   + VSE);
                    vxm24 = vxm24 + 0.5*pow((dxm/dx-0.5),2.0) * (VNWW   -2.0*VNW   + VNE);
                }
            }
        }

        (*Vxm) = (1.0-dzm/dz) * vxm13 + (dzm/dz) * vxm24;
        
        
//        (*Vxm) =
        
        //---------------------------------------------------//
        
        Nx = NxVz; Nz = Nzg;
        
        dst    = (particles->x[k]-xvz[0]);
        j_part = ceil((dst/dx)) - 1;
        if (j_part<0) {
            j_part = 0;
        }
        if (j_part>Nx-2) {
            j_part = Nx-2;
        }
        
        // Get the line:
        dst    = (particles->z[k]-zg[0]);
        i_part = ceil((dst/dz)) - 1;
        if (i_part<0) {
            i_part = 0;
        }
        if (i_part>Nz-2) {
            i_part = Nz-2;
        }
        
        dxm = fabs(particles->x[k] - xvz[j_part]);
        dzm = fabs(particles->z[k] - zg[i_part]);
        
        iSW = j_part+i_part*Nx;
        iSE = j_part+i_part*Nx+1;
        iNW = j_part+(i_part+1)*Nx;
        iNE = j_part+(i_part+1)*Nx+1;
        
        VSE = Vz[iSE];
        VSW = Vz[iSW];
        VNE = Vz[iNE];
        VNW = Vz[iNW];

        if (tagVz[iSW] == 30) VSW = Vz[iNE];
        if (tagVz[iSE] == 30) VSE = Vz[iNE];
        if (tagVz[iNW] == 30) VNW = Vz[iSW];
        if (tagVz[iNE] == 30) VNE = Vz[iSE];

        // if (tagVz[iNE] == 30 && tagVz[iSE] == 30 && tagVz[iSW] == 30) exit(1130);
        // if (tagVz[iNW] == 30 && tagVz[iSW] == 30 && tagVz[iSE] == 30) exit(1131);
        
        if (tagVz[iNE] == 30 && tagVz[iSE] == 30 ) {
            VNE = VSW;
            VSE = VSW;
        }
        
        if (tagVz[iNW] == 30 && tagVz[iSW] == 30) {
            VNW = VSE;
            VSW = VSE;
        }
        
        vzm12 = VSW*(1.0-dzm/dz) + VNW*dzm/dz;
        vzm34 = VSE*(1.0-dzm/dz) + VNE*dzm/dz;
        
        // Compute correction
        if (TarasCorr==1) {
            if(dzm/dz>=0.5) {
                if(i_part<Nz-2) {
                    iNNW  = j_part+(i_part+2)*Nx;
                    iNNE  = j_part+(i_part+2)*Nx + 1;
                    VNNW  = Vz[iNNW];
                    VNNE  = Vz[iNNE];
                    if (tagVz[iNNW] == 30) VNNW = Vz[iSW];
                    if (tagVz[iNNE] == 30) VNNE = Vz[iSE];
                    vzm12 = vzm12 + 0.5*pow((dzm/dz-0.5),2.0) * (VSW - 2.0*VNW + VNNW);
                    vzm34 = vzm34 + 0.5*pow((dzm/dz-0.5),2.0) * (VSE - 2.0*VNE + VNNE);
                }
            }
            else {
                if(i_part>0) {
                    iSSW  = j_part+(i_part-1)*Nx;
                    iSSE  = j_part+(i_part-1)*Nx + 1;
                    VSSW = Vz[iSSW];
                    VSSE = Vz[iSSE];
                    if (tagVz[iSSW] == 30) VSSW = Vz[iNW];
                    if (tagVz[iSSE] == 30) VSSE = Vz[iNE];
                    vzm12 = vzm12 + 0.5*pow((dzm/dz-0.5),2.0) * (VSSW - 2.0*VSW + VNW);
                    vzm34 = vzm34 + 0.5*pow((dzm/dz-0.5),2.0) * (VSSE - 2.0*VSE + VNE);
                }
            }
        }
        
        (*Vzm) = (1.0-dxm/dx)*vzm12 + (dxm/dx)*vzm34;
        
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


double Centers2Particle( markers* particles, double* NodeField, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag, double dx, double dz, int k, int periodix ) {
    
    double val  = 0.0;
    double sumW = 0.0;
    double dxm, dzm, dst;
    int    i_part, j_part, iSW, iNW, iSE, iNE, ic, jc;
    
    // Filter out particles that are inactive (out of the box)
    if (particles->phase[k] != -1) {

        // Check distance with regards to extended array
        dst    = (particles->x[k]-X_vect[0]);
        j_part = ceil((dst/dx)) - 1;
        
        // Get the line:
        dst    = (particles->z[k]-Z_vect[0]);
        i_part = ceil((dst/dz)) - 1;
        
        if (j_part>0 && j_part<Nx && i_part>0 && i_part<Nz) {
            ic  = i_part - 1;
            jc  = j_part - 1;
            iSW = jc + ic*Nx;
            iSE = jc + ic*Nx+1;
            iNW = jc + (ic+1)*Nx;
            iNE = jc + (ic+1)*Nx+1;
        }
        
        if (j_part==0 && i_part > 0 && i_part<Nz ) {
            if (periodix==1) {
                ic  = i_part-1;
                jc  = j_part;
                iSW = (Nx-1) + ic*Nx;
                iSE = jc     + ic*Nx;
                iNW = (Nx-1) + (ic+1)*Nx;
                iNE = jc     + (ic+1)*Nx;
            }
            else {
                // Copy solution along edge
                ic  = i_part-1;
                jc  = j_part;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx;
                iNW = jc + (ic+1)*Nx;
                iNE = jc + (ic+1)*Nx;
            }
        }

        if (j_part==Nx && i_part > 0 && i_part<Nz ) {
            if (periodix==1) {
                ic  = i_part-1;
                jc  = j_part-1;
                iSW = jc + ic*Nx;
                iSE = 0  + ic*Nx;
                iNW = jc + (ic+1)*Nx;
                iNE = 0  + (ic+1)*Nx;
            }
            else {
                // Copy solution along edge
                ic  = i_part-1;
                jc  = j_part-1;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx;
                iNW = jc + (ic+1)*Nx;
                iNE = jc + (ic+1)*Nx;
            }
        }
        
        if (i_part==0 && j_part > 0 && j_part<Nx) {
            // Copy solution along edge
            ic  = i_part;
            jc  = j_part-1;
            iSW = jc + ic*Nx;
            iSE = jc + ic*Nx + 1;
            iNW = jc + ic*Nx;
            iNE = jc + ic*Nx + 1;
        }

        if (i_part==Nz && j_part > 0 && j_part<Nx) {
            // Copy solution along edge
            ic  = i_part-1;
            jc  = j_part-1;
            iSW = jc + ic*Nx;
            iSE = jc + ic*Nx + 1;
            iNW = jc + ic*Nx;
            iNE = jc + ic*Nx + 1;
        }
        
        // SW
        if (i_part==0 && j_part==0) {

            if (periodix==1) {
                ic  = 0;
                jc  = 0;
                iSW = (Nx-1) + ic*Nx;
                iSE = jc     + ic*Nx;
                iNW = (Nx-1) + ic*Nx;
                iNE = jc     + ic*Nx;
            }
            else {
                ic  = i_part - 0;
                jc  = j_part - 0;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx;
                iNW = jc + ic*Nx;
                iNE = jc + ic*Nx;
            }
        }
        // SE
        if (i_part==0 && j_part==Nx) {

            if (periodix==1) {
                ic  = i_part - 0;
                jc  = j_part - 1;
                iSW = jc + ic*Nx;
                iSE = 0  + ic*Nx;
                iNW = jc + ic*Nx;
                iNE = 0  + ic*Nx;
            }
            else {
                ic  = i_part - 0;
                jc  = j_part - 1;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx;
                iNW = jc + ic*Nx;
                iNE = jc + ic*Nx;
            }
        }
        // NW
        if (i_part==Nz && j_part==0) {

            if (periodix==1) {
                ic  = i_part - 1;
                jc  = j_part - 0;
                iSW = (Nx-1) + ic*Nx;
                iSE = jc     + ic*Nx;
                iNW = (Nx-1) + ic*Nx;
                iNE = jc     + ic*Nx;
            }
            else {
                ic  = i_part - 1;
                jc  = j_part - 0;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx;
                iNW = jc + ic*Nx;
                iNE = jc + ic*Nx;
            }
        }
        // NE
        if (i_part==Nz && j_part==Nx) {

            if (periodix==1) {
                ic  = i_part - 1;
                jc  = j_part - 1;
                iSW = jc + ic*Nx;
                iSE = 0  + ic*Nx;
                iNW = jc + ic*Nx;
                iNE = 0  + ic*Nx;
            }
            else {
                ic  = i_part - 1;
                jc  = j_part - 1;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx;
                iNW = jc + ic*Nx;
                iNE = jc + ic*Nx;
            }
        }
    
        dxm = (particles->x[k] - X_vect[j_part]);
        dzm = (particles->z[k] - Z_vect[i_part]);

        if (tag[iSW]!=30 && tag[iSW]!=31) {
            val  += (1.0-dxm/dx) * (1.0-dzm/dz) * NodeField[iSW];
            sumW += (1.0-dxm/dx) * (1.0-dzm/dz);

        }
        if (tag[iSE]!=30 && tag[iSE]!=31) {
            val  += (dxm/dx) * (1.0-dzm/dz)  * NodeField[iSE];
            sumW += (dxm/dx) * (1.0-dzm/dz);
        }
        if (tag[iNW]!=30 && tag[iNW]!=31) {
            val  += (1.0-dxm/dx) * (dzm/dz)    * NodeField[iNW];
            sumW += (1.0-dxm/dx) * (dzm/dz);
        }
        if (tag[iNE]!=30 && tag[iNE]!=31) {
            val  += (dxm/dx) * (dzm/dz)    * NodeField[iNE];
            sumW += (dxm/dx) * (dzm/dz);
        }

        if(sumW>1e-13) val /= sumW;
        
    }
    
    return val;
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


double Vertices2Particle( markers* particles, double* NodeField, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag, double dx, double dz, int k ) {
    
    double val  = 0.0;
    double sumW = 0.0;
    double dxm, dzm, dst;
    int    i_part, j_part, iSW, iNW, iSE, iNE;
    
    // Filter out particles that are inactive (out of the box)
    if (particles->phase[k] != -1) {
        
        dst    = fabs(particles->x[k]-X_vect[0]);
        j_part = ceil((dst/dx)) - 1;
        if (j_part<0) {
            // printf("Should never be here I! (Vertices2Particle)\n");
            // printf("%2.10e %2.10e\n", particles->x[k],  X_vect[0]);
            // exit(1);
            j_part = 0;
        }
        if (j_part>Nx-2) {
            printf("Should never be here II! (Vertices2Particle)\n"); exit(1);
            j_part = Nx-2;
        }
        
        // Get the line:
        dst    = fabs(particles->z[k]-Z_vect[0]);
        i_part = ceil((dst/dz)) - 1;
        if (i_part<0) {
            // printf("particles->x[k] = %lf, zmin= %lf\n", particles->z[k], Z_vect[0]);
            // printf("Should never be here IV! (Vertices2Particle)\n"); exit(1);
            i_part = 0;
        }
        if (i_part>Nz-2) {
            // printf("Should never be here V! (Vertices2Particle)\n"); exit(1);
            i_part = Nz-2;
        }
        
        dxm = particles->x[k] - X_vect[j_part];
        dzm = particles->z[k] - Z_vect[i_part];
        
        iSW = j_part+i_part*Nx;
        iSE = j_part+i_part*Nx+1;
        iNW = j_part+(i_part+1)*Nx;
        iNE = j_part+(i_part+1)*Nx+1;
        
        if (tag[iSW]!=30 && tag[iSW]!=31) {
            val  += (1.0-dxm/dx) * (1.0-dzm/dz) * NodeField[iSW];
            sumW += (1.0-dxm/dx) * (1.0-dzm/dz);
        }
        if (tag[iSE]!=30 && tag[iSE]!=31) {
            val  += (dxm/dx) * (1.0-dzm/dz)  * NodeField[iSE];
            sumW += (dxm/dx) * (1.0-dzm/dz);
        }
        if (tag[iNW]!=30 && tag[iNW]!=31) {
            val  += (1.0-dxm/dx) * (dzm/dz)    * NodeField[iNW];
            sumW += (1.0-dxm/dx) * (dzm/dz);
        }
        if (tag[iNE]!=30 && tag[iNE]!=31) {
            val  += (dxm/dx) * (dzm/dz)    * NodeField[iNE];
            sumW += (dxm/dx) * (dzm/dz);
        }
        if(sumW>1e-13) val /= sumW;
        

//
//        val =  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iSW];
//        val += (dxm/dx)    * (1.0-dzm/dz) * NodeField[iSE];
//        val += (1.0-dxm/dx)* (dzm/dz)     * NodeField[iNW];
//        val += (dxm/dx)    * (dzm/dz)     * NodeField[iNE];
//
//        sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
//        sumW += (dxm/dx)* (1.0-dzm/dz);
//        sumW += (1.0-dxm/dx)* (dzm/dz);
//        sumW += (dxm/dx)* (dzm/dz);
//
//        if( sumW > 1e-13 ) val /= sumW;

        
    }
    else {
        printf("Should never be here VI ! (Vertices2Particle)\n"); exit(1);
    }
    return val;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void ParticleInflowCheck ( markers* particles, grid *mesh, MdoodzInput *input, surface topo, int flag, SetParticles_ff setParticles  ) {
    
    params model = input->model;
    int k, Nb_part=particles->Nb_part, npW=0, npE=0;
    double xW, xE, dx2=model.dx/2.0;
    int finite_strain = model.finite_strain;
    
    clock_t t_omp = (double)omp_get_wtime();

    // Count particles ready for recycling
    int np_reuse = 0;
    for (k=0; k<Nb_part; k++) {
        if ( particles->phase[k] == -1 ) np_reuse++;
    }
    printf( GREEN "Particles available for recycling: %04d\n" RESET, np_reuse);

    // Stores indices of particles to recycle
    int *ip_reuse = DoodzCalloc( np_reuse, sizeof(int));
    np_reuse = 0;
    for (k=0; k<Nb_part; k++) {
        if ( particles->phase[k] == -1 ) {
            ip_reuse[np_reuse] = k; 
            np_reuse++;
        }
    }
    
    // Check if particles are in the first west/east half-cell column
    if (flag == 0) {
        xW=mesh->xg_coord[0];
        xE=mesh->xg_coord[model.Nx-1];
        
#pragma omp parallel for shared (particles ) \
private( k ) \
firstprivate( Nb_part, xW, xE, dx2 ) \
reduction(+:npW,npE )
        for (k=0; k<Nb_part; k++) {
            
            particles->intag[k] = 0;
            
            if ( particles->phase[k] != -1 ) {
                
                // WEST COAST
                if (particles->x[k]>=xW && particles->x[k]<=xW+dx2) {
                    particles->intag[k] = 1;
                    npW += 1;
                }
                
                // EAST COAST
                if (particles->x[k]>=xE-dx2 && particles->x[k]<=xE) {
                    particles->intag[k] = 1;
                    npE += 1;
                }
            }
        }
        printf("Number of particles west boundary: %d\n", npW);
        printf("Number of particles east boundary: %d\n", npE);
    }
    
    if (flag == 1) {
        
        xW=mesh->xg_coord[0]          + dx2;
        xE=mesh->xg_coord[model.Nx-1] - dx2;

        int NewIndex;
        
        for (k=0; k<Nb_part; k++) {
            
            if ( particles->phase[k] != -1 && particles->intag[k] == 1 ) {
                
                // WEST COAST
                if (particles->x[k]>=xW && particles->x[k]<=xW+dx2) {
                    particles->intag[k] = 2;
                    npW += 1;
                    
                    // Add particles
                    if (np_reuse>0) {
                        // Recycle
                        NewIndex = ip_reuse[np_reuse-1];
                        np_reuse -= 1;
                    }
                    else if (Nb_part < particles->Nb_part_max) {
                        // Increase number of particles
                        NewIndex = Nb_part;
                        Nb_part++;
                    }
                    else {
                        printf("Too many particles in my mind %d, maximum %d\n", Nb_part, particles->Nb_part_max);
                    }

                    // Add particles
                    particles->x[NewIndex] = particles->x[k]-dx2;
                    particles->z[NewIndex] = particles->z[k];
                    // Assign new marker point properties
                    AssignMarkerPropertiesInflow ( particles, NewIndex, k, input, mesh, 1, setParticles );
                }
                
                // EAST COAST
                if (particles->x[k]>=xE-dx2 && particles->x[k]<=xE) {
                    particles->intag[k] = 2;
                    npE += 1;

                    // Add particles
                    if (np_reuse>0) {
                        // Recycle
                        NewIndex = ip_reuse[np_reuse-1];
                        np_reuse -= 1;
                    }
                    else if (Nb_part < particles->Nb_part_max) {
                        // Increase number of particles
                        NewIndex = Nb_part;
                        Nb_part++;
                    }
                    else {
                        printf("Too many particles in my mind %d, maximum %d\n", Nb_part, particles->Nb_part_max);
                    }
                    
                    // Add particles
                    particles->x[NewIndex] = particles->x[k]+dx2;
                    particles->z[NewIndex] = particles->z[k];
                    // Assign new marker point properties
                    AssignMarkerPropertiesInflow ( particles, NewIndex, k, input, mesh, 1, setParticles );
                }
            }
        }
        printf("Number of particles west boundary: %d\n", npW);
        printf("Number of particles east boundary: %d\n", npE);
        particles->Nb_part = Nb_part;
    }
    
    printf("** Time for particles inflow check = %lf sec\n",  (double)((double)omp_get_wtime() - t_omp) );
    
    DoodzFree(ip_reuse);

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AssignMarkerProperties (markers* particles, int new_ind, int min_index, params *model, grid* mesh, int direct_neighbour ) {
    
    particles->phase[new_ind]         = particles->phase[min_index];
    particles->dual[new_ind]          = particles->dual[min_index];
    if (particles->phase[min_index]==-1) {printf("New particle is has phase -1, error! (AssignMarkerProperties)\n" ); exit(99);}
    particles->Vx[new_ind]            = particles->Vx[min_index];
    particles->Vz[new_ind]            = particles->Vz[min_index];
    particles->strain[new_ind]        = particles->strain[min_index];
    particles->strain_el[new_ind]     = particles->strain_el[min_index];
    particles->strain_pl[new_ind]     = particles->strain_pl[min_index];
    particles->strain_pwl[new_ind]    = particles->strain_pwl[min_index];
    particles->strain_exp[new_ind]    = particles->strain_exp[min_index];
    particles->strain_lin[new_ind]    = particles->strain_lin[min_index];
    particles->strain_gbs[new_ind]    = particles->strain_gbs[min_index];
    particles->divth[new_ind]         = particles->divth[min_index]; // to be changed
    if ( direct_neighbour == 1 ) {
        particles->d[new_ind]             = particles->d[min_index];
        particles->T[new_ind]             = particles->T[min_index];
        particles->P[new_ind]             = particles->P[min_index];
        particles->phi[new_ind]           = particles->phi[min_index]; // to be changed
        particles->X[new_ind]             = particles->X[min_index];   // to be changed
        particles->noise[new_ind]         = particles->noise[min_index];   // to be changed
        particles->sxxd[new_ind]          = particles->sxxd[min_index];
        particles->szzd[new_ind]          = particles->szzd[min_index];
        particles->sxz[new_ind]           = particles->sxz[min_index];
    }
    else {
        particles->d[new_ind]             = Centers2Particle( particles, mesh->d_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->T[new_ind]             = Centers2Particle( particles, mesh->T,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->P[new_ind]             = Centers2Particle( particles, mesh->p_in,  mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->phi[new_ind]           = Centers2Particle( particles, mesh->phi_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->X[new_ind]             = Centers2Particle( particles, mesh->X_n,  mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->noise[new_ind]         = Centers2Particle( particles, mesh->noise_n,  mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->sxxd[new_ind]          = Centers2Particle( particles, mesh->sxxd,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->szzd[new_ind]          = Centers2Particle( particles, mesh->szzd,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->sxz[new_ind]           = Vertices2Particle( particles, mesh->sxz,     mesh->xg_coord,  mesh->zg_coord,  mesh->Nx-0, mesh->Nz-0, mesh->BCg.type, mesh->dx, mesh->dz, new_ind );
    }
    // particles->dsxxd[new_ind]         = particles->dsxxd[min_index];
    // particles->dszzd[new_ind]         = particles->dszzd[min_index];
    // particles->dsxz[new_ind]          = particles->dsxz[min_index];
    // particles->syy[new_ind]           = particles->syy[min_index];
    // particles->dsyy[new_ind]           = particles->dsyy[min_index];
     
    if (model->finite_strain == 1) {
        // do not set default to 0 beause then it can not accumulate, better to identify which markers are new and start to accumulate as we do for the general case (fxx=fyy=1, fxz=fzx=0).
        particles->Fxx[new_ind]           = particles->Fxx[min_index];
        particles->Fxz[new_ind]           = particles->Fxz[min_index];
        particles->Fzx[new_ind]           = particles->Fzx[min_index];
        particles->Fzz[new_ind]           = particles->Fzz[min_index];
    }
    if (model->track_T_P_x_z == 1) {
        particles->T0[new_ind]           = particles->T0[min_index];
        particles->P0[new_ind]           = particles->P0[min_index];
        particles->x0[new_ind]           = particles->x0[min_index];
        particles->z0[new_ind]           = particles->z0[min_index];
        particles->Tmax[new_ind]         = particles->Tmax[min_index];
        particles->Pmax[new_ind]         = particles->Pmax[min_index];
    }
    if (model->anisotropy == 1) {
        // !!!! difference with AssignMarkerPropertiesInflow:
        // inside the box AssignMarkerProperties does not need to know out inflow
        particles->nx[new_ind]           = particles->nx[min_index];
        particles->nz[new_ind]           = particles->nz[min_index];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AssignMarkerPropertiesInflow (markers* particles, int new_ind, int min_index, MdoodzInput *input, grid* mesh, int direct_neighbour, SetParticles_ff setParticles ) {
    
    params *model = &(input->model);  
    particles->phase[new_ind]         = particles->phase[min_index];
    particles->dual[new_ind]          = particles->dual[min_index];
    if (particles->phase[min_index]==-1) {printf("New particle is has phase -1, error! (AssignMarkerProperties)\n" ); exit(99);}
    particles->Vx[new_ind]            = particles->Vx[min_index];
    particles->Vz[new_ind]            = particles->Vz[min_index];
    particles->strain[new_ind]        = particles->strain[min_index];
    particles->strain_el[new_ind]     = particles->strain_el[min_index];
    particles->strain_pl[new_ind]     = particles->strain_pl[min_index];
    particles->strain_pwl[new_ind]    = particles->strain_pwl[min_index];
    particles->strain_exp[new_ind]    = particles->strain_exp[min_index];
    particles->strain_lin[new_ind]    = particles->strain_lin[min_index];
    particles->strain_gbs[new_ind]    = particles->strain_gbs[min_index];
    particles->divth[new_ind]         = particles->divth[min_index]; // to be changed
    if ( direct_neighbour == 1 ) {
        particles->d[new_ind]             = particles->d[min_index];
        particles->T[new_ind]             = particles->T[min_index];
        particles->P[new_ind]             = particles->P[min_index];
        particles->phi[new_ind]           = particles->phi[min_index]; // to be changed
        particles->X[new_ind]             = particles->X[min_index];   // to be changed
        particles->noise[new_ind]         = particles->noise[min_index];   // to be changed
        particles->sxxd[new_ind]          = particles->sxxd[min_index];
        particles->szzd[new_ind]          = particles->szzd[min_index];
        particles->sxz[new_ind]           = particles->sxz[min_index];
    }
    else {
        particles->d[new_ind]             = Centers2Particle( particles, mesh->d_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->T[new_ind]             = Centers2Particle( particles, mesh->T,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->P[new_ind]             = Centers2Particle( particles, mesh->p_in,  mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->phi[new_ind]           = Centers2Particle( particles, mesh->phi_n,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->X[new_ind]             = Centers2Particle( particles, mesh->X_n,  mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->noise[new_ind]         = Centers2Particle( particles, mesh->noise_n,  mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->sxxd[new_ind]          = Centers2Particle( particles, mesh->sxxd,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->szzd[new_ind]          = Centers2Particle( particles, mesh->szzd,     mesh->xvz_coord, mesh->zvx_coord, mesh->Nx-1, mesh->Nz-1, mesh->BCp.type, mesh->dx, mesh->dz, new_ind, model->periodic_x );
        particles->sxz[new_ind]           = Vertices2Particle( particles, mesh->sxz,     mesh->xg_coord,  mesh->zg_coord,  mesh->Nx-0, mesh->Nz-0, mesh->BCg.type, mesh->dx, mesh->dz, new_ind );
    }
    // particles->dsxxd[new_ind]         = particles->dsxxd[min_index];
    // particles->dszzd[new_ind]         = particles->dszzd[min_index];
    // particles->dsxz[new_ind]          = particles->dsxz[min_index];
    // particles->syy[new_ind]           = particles->syy[min_index];
    // particles->dsyy[new_ind]           = particles->dsyy[min_index];
     
    if (input->model.finite_strain == 1) {
        // do not set default to 0 beause then it can not accumulate, better to identify which markers are new and start to accumulate as we do for the general case (fxx=fyy=1, fxz=fzx=0).
        particles->Fxx[new_ind]           = particles->Fxx[min_index];
        particles->Fxz[new_ind]           = particles->Fxz[min_index];
        particles->Fzx[new_ind]           = particles->Fzx[min_index];
        particles->Fzz[new_ind]           = particles->Fzz[min_index];
    }
    if (input->model.track_T_P_x_z == 1) {
        particles->T0[new_ind]           = particles->T0[min_index];
        particles->P0[new_ind]           = particles->P0[min_index];
        particles->x0[new_ind]           = particles->x0[min_index];
        particles->z0[new_ind]           = particles->z0[min_index];
        particles->Tmax[new_ind]         = particles->Tmax[min_index];
        particles->Pmax[new_ind]         = particles->Pmax[min_index];
    }
    if (input->model.anisotropy == 1) {
        // add with William: the incoming markers does not need to inherit the nearest neighbour angle
           Coordinates coordinates = {
            .x = particles->x[new_ind],
            .z = particles->z[new_ind]};
        if (input->model.marker_aniso_angle && setParticles.SetAnisoAngle) {
            double predefined_angle = acos(particles->nx[min_index]);
            double angle = setParticles.SetAnisoAngle(input, coordinates, particles->phase[new_ind], predefined_angle);
            particles->nx[new_ind]           = cos(angle);
            particles->nz[new_ind]           = sin(angle);
        }
        else {
            particles->nx[new_ind]           = particles->nx[min_index];
            particles->nz[new_ind]           = particles->nz[min_index];
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/
//
//void PrintMarkerProperties (markers* particles, int new_ind, int min_index ) {
//    printf("xp = %2.6e\n", particles->x[new_ind]*30e4/1000         );
//    printf("zp = %2.6e\n", particles->z[new_ind]*30e4/1000           );
//    printf("xn = %2.6e\n", particles->x[min_index]*30e4/1000           );
//    printf("zn = %2.6e\n", particles->z[min_index]*30e4/1000           );
//    printf("%d\n", particles->phase[new_ind]        );
//    printf("%2.6e\n", particles->sxxd[new_ind]         );
//    printf("%2.6e\n", particles->sxz[new_ind]          );
//    printf("%2.6e\n", particles->Vx[new_ind]           );
//    printf("%2.6e\n", particles->Vz[new_ind]           );
//    printf("%2.6e\n", particles->strain[new_ind]       );
//    printf("%2.6e\n", particles->strain_el[new_ind]    );
//    printf("%2.6e\n", particles->strain_pl[new_ind]    );
//    printf("%2.6e\n", particles->strain_pwl[new_ind]   );
//    printf("%2.6e\n", particles->strain_exp[new_ind]   );
//    printf("%2.6e\n", particles->strain_lin[new_ind]   );
//    printf("%2.6e\n", particles->strain_gbs[new_ind]    );
//    printf("%2.6e\n", particles->d[new_ind]            );
//    printf("%2.6e\n", particles->T[new_ind]             );
//    printf("%2.6e\n", particles->P[new_ind]             );
//    printf("%2.6e\n", particles->phi[new_ind]           );
//    printf("%2.6e\n", particles->X[new_ind]             );
//    printf("%d\n", particles->generation[new_ind]    );
//}
//
///*--------------------------------------------------------------------------------------------------------------------*/
///*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
///*--------------------------------------------------------------------------------------------------------------------*/

void PartInit( markers *particles, params* model ) {

    // Initialise particle fields to zero.
    int k;

    #pragma omp parallel for shared ( particles ) private ( k ) schedule ( static )
    for(k=0; k<particles->Nb_part_max; k++) {

        particles->x[k]     = 0.0;
        particles->z[k]     = 0.0;
        particles->P[k]     = 0.0;

        // Viscoelasticity
        particles->sxxd[k]    = 0.0;
        particles->sxz[k]     = 0.0;

        // phase structure
        particles->phase[k]      = 0.0;
        particles->dual[k]       = 0.0;
        particles->strain[k]     = 0.0;
        particles->strain_el[k]  = 0.0;
        particles->strain_pl[k]  = 0.0;
        particles->strain_pwl[k] = 0.0;
        particles->strain_exp[k] = 0.0;
        particles->strain_lin[k] = 0.0;
        particles->strain_gbs[k] = 0.0;
        particles->T[k]          = 0.0;
        particles->phi[k]        = 0.0;
        particles->X[k]          = 0.0;

        // Reaction front progress
        particles->generation[k] = 0.0;
        particles->progress[k]   = 0.0;

        // Finite strain - deformation gradient tensor
        if ( model->finite_strain == 1 ) {
            particles->Fxx[k]   = 1.0;
            particles->Fxz[k]   = 0.0;
            particles->Fzx[k]   = 0.0;
            particles->Fzz[k]   = 1.0;
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FindClosestPhase( markers* particles, int ic, int jc, grid mesh, int *ind_list, int new_ind, int nb_neigh, params* model ) {

    // Find phase of a new particle according to its closest neighbour. This function works at the scale of a cell.
    // The array 'ind_list' is list of indices of particle already owing to that cell.

    double dst, min_dst;
    int ind, min_index=ind_list[0];


    min_dst = 20.0*sqrt( pow( mesh.dx, 2.0) + pow( mesh.dz, 2.0) );

    // Loop on the particles inside the current cell, find which one is the closest to the newly created one.
    for (ind=0; ind<nb_neigh; ind++) {

        dst = sqrt( pow (particles->x[new_ind] - particles->x[ind_list[ind]], 2.0) + pow(particles->z[new_ind] - particles->z[ind_list[ind]], 2.0) );

        if (dst<min_dst) {
            min_dst = dst;
            min_index = ind_list[ind];
        }
    }

    // Assign new marker point properties
    AssignMarkerProperties ( particles, new_ind, min_index, model, &mesh, model->direct_neighbour );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int FindClosestPhase2( markers* particles, double *newx, double* newz, int ic, int jc, grid mesh, int *ind_list, int new_ind, int nb_neigh, int sedim, int sed_phase,  params model ) {

    // Find phase of a new particle according to its closest neighbour. This function works at the scale of a cell.
    // The array 'ind_list' is list of indices of particle already owing to that cell.

    double dst, min_dst;
    int ind, min_index=ind_list[0];

    min_dst = (model.xmax-model.xmin);//*sqrt( pow( mesh.dx, 2.0) + pow( mesh.dz, 2.0) );

    if (sedim==0) {

        // Loop on the particles inside the current cell, find which one is the closest to the newly created one.
        for (ind=0; ind<nb_neigh; ind++) {

            dst = sqrt( pow (newx[new_ind] - particles->x[ind_list[ind]], 2.0) + pow(newz[new_ind] - particles->z[ind_list[ind]], 2.0) );

            if (dst<min_dst) {
                min_dst = dst;
                min_index = ind_list[ind];
            }
        }

    }
    else {
        //        printf("Looking for sediments\n");

        for (ind=0; ind<nb_neigh; ind++) {

            dst = sqrt( pow (newx[new_ind] - particles->x[ind_list[ind]], 2.0) + pow(newz[new_ind] - particles->z[ind_list[ind]], 2.0) );

            if (dst<min_dst && ( particles->phase[ind_list[ind]] ==  model.surf_ised1 || particles->phase[ind_list[ind]] ==  model.surf_ised2) ) {
                min_dst = dst;
                min_index = ind_list[ind];
            }
        }

    }

    return min_index;

    // Assign new marker point properties
    // AssignMarkerProperties ( particles, new_ind, min_index, &mesh, model->direct_neighbour );
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FindClosestPhaseVertex( markers* particles, int ic, int jc, grid mesh, int *ind_list, int new_ind, int nb_neigh, params* model ) {

    // Find phase of a new particle according to its closest neighbour.
    // This function works at the scale of a cell.
    // The array 'ind_list' is list of indices of particle already owing to that cell.

    double dst, min_dst;
    int ind, min_index=ind_list[0];

    //    min_dst = sqrt( pow( 4*mesh.dxgrid, 2) + pow( 4*mesh.dzgrid, 2) );
    min_dst = 20.0*sqrt( pow( mesh.dx, 2.0) + pow( mesh.dz, 2.0) );

    for (ind=0; ind<nb_neigh; ind++) {

        dst = sqrt( pow (particles->x[new_ind] - particles->x[ind_list[ind]], 2.0) + pow(particles->z[new_ind] - particles->z[ind_list[ind]], 2.0) );

        if (dst<min_dst) {
            min_dst = dst;
            min_index = ind_list[ind];
        }
    }

    if (particles->phase[min_index] <-1 || particles->phase[min_index]>50 ) {
        printf("Could not find matching phase for newly created particle, this is a bug (attribute phase id : %d)\n Error at node ic = %d jc = %d \nExiting...", particles->phase[min_index], ic, jc);
        exit(50);
    }

    // Assign new marker point properties
    AssignMarkerProperties ( particles, new_ind, min_index, model, &mesh, model->direct_neighbour );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartCell( markers *particles, grid mesh, int ic, int jc, int* ind_list, params model, int nb_neigh, int *ind_part_reuse, int *inc_reuse, int nb_part_reuse, surface topo  ) {

    // This functions add particles to cells that are tag as 'particle deficient'.
    // So far, it adds 4 new particles at different locations within the cell.
    // The new particles get the phase and stored stresses of the closest particle.
    int nxp=1, nzp=1;
    int Nb_add = nxp*nzp;
    double new_x, new_z, h=model.zmax;
    int new_ind;
    int k, l;

    // Reallocate particle arrays if no more space in arrays (TO BE IMPLEMENTED)
    if (particles->Nb_part + Nb_add > particles->Nb_part_max && Nb_add > nb_part_reuse) {
        printf("You have reached the maximum number of particles currently available (%d), please increase it...\n", particles->Nb_part_max );
        printf("Exiting...\n");
        exit(1);
    }

    for (l=0;l<nzp;l++) {
        for (k=0;k<nxp;k++) {

            //-------------------------------------------------------------------------------------------------------------------------//

            new_x = mesh.xg_coord[ic] + (1+k)*mesh.dx/(nxp+1) - 0.0*mesh.dx/10.0;
            new_z = mesh.zg_coord[jc] + (1+l)*mesh.dz/(nzp+1) - 0.0*mesh.dz/10.0;

            if ( new_x > model.xmin && new_z > model.zmin && new_z<h  ) {
                if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
                    new_ind = ind_part_reuse[*inc_reuse];
                    (*inc_reuse)++;
                }
                else {
                    new_ind = particles->Nb_part;
                    particles->Nb_part++;
                }

                // Add 4 particules
                particles->x[new_ind]      = new_x;
                particles->z[new_ind]      = new_z;
                FindClosestPhase( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh, &model  );
            }

        }
    }

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartCell2( int *pidx, int *Nb_part_thread, markers *particles, grid mesh, int ic, int jc, int* ind_list, params model, int nb_neigh, surface topo, double* xg , double* zg, int cell, int *nnewp, double *newx, double *newz, int *newi, int sed_phase, surface topo_ini ) {

    // This functions add particles to cells that are tag as 'particle deficient'.
    // So far, it adds 4 new particles at different locations within the cell.
    // The new particles get the phase and stored stresses of the closest particle.


    int nxp=1, nzp=1;
    double new_x, new_z, h=model.zmax, h_ini=model.zmax;
    int  k, l, sedim;

    //    // Reallocate particle arrays if no more space in arrays (TO BE IMPLEMENTED)
    //    if (Nb_part_thread + Nb_add > particles->Nb_part_max && Nb_add > nb_part_reuse) {
    //        printf("You have reached the maximum number of particles currently available (%d), please increase it...\n", particles->Nb_part_max );
    //        printf("Exiting...\n");
    //        exit(1);
    //    }

    for (l=0;l<nzp;l++) {
        for (k=0;k<nxp;k++) {

            //-------------------------------------------------------------------------------------------------------------------------//

            new_x = xg[0] + ic*mesh.dx/2.0 + mesh.dx/4.0;
            new_z = zg[0] + jc*mesh.dz/2.0 + mesh.dz/4.0;

            if ( model.free_surface==1 ) {
//                h0    =  topo.a0[cell]*new_x      + topo.b0[cell];
                h     =  topo.a[cell]*new_x      +  topo.b[cell];
                h_ini =  topo_ini.a[cell]*new_x  +  topo_ini.b[cell];
            }

            if ( new_x > model.xmin && new_z > model.zmin  && new_z<h ) { //
                (*Nb_part_thread)++;

                sedim = 0;
                if (model.surface_processes>0 && new_z > h_ini ) sedim = 1;

                // Add 4 particules
                newx[(*nnewp)] = new_x;
                newz[(*nnewp)] = new_z;
                newi[(*nnewp)] = FindClosestPhase2( particles, newx, newz, ic, jc, mesh, ind_list, (*nnewp), nb_neigh, sedim, sed_phase, model  );

                (*nnewp)++;
            }

        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartVert( markers *particles, grid mesh, int ic, int jc, int* ind_list, params model, int nb_neigh, int *ind_part_reuse, int *inc_reuse, int nb_part_reuse, surface topo ) {

    // This functions add particles to cells that are tag as 'particle deficients'.
    // So far, it adds 4 new particles at different locations within the cell.
    // The new particles get the phase and stored stress of the closest particle.

    int Nb_add = 4;
    double new_x, new_z, h=model.zmax;
    int new_ind;
    int cell;

    // Reallocate particle arrays if no more space in arrays (TO BE IMPLEMENTED)
    if (particles->Nb_part + Nb_add > particles->Nb_part_max && Nb_add > nb_part_reuse ) {
        printf("You have reached the maximum number of particles currently available (%d), please increase it...\n", particles->Nb_part_max );
        printf("Exiting...\n");
        exit(1);
    }

    // Add 4 particules
    if ( Nb_add == 4 ) {

        cell = ic;

        //-------------------------------------------------------------------------------------------------------------------------//
        // Get the particle index
        new_x = mesh.xg_coord[ic] - mesh.dx/3.0;
        new_z = mesh.zg_coord[jc] - mesh.dz/3.0;
        //        distance=(new_x - xmin);
        //        cell = ceil((distance/model.dx)+0.5) - 1;
        //        if (cell<0) cell = 0;
        //        if (cell>mesh.Nx[0]-2) cell = mesh.Nx[0]-2;
        if ( model.free_surface==1 ) h  = topo.a[cell]*new_x  + topo.b[cell];


        if ( new_x > model.xmin && new_z > model.zmin && new_z<h   ) {

            //            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
            //                new_ind = ind_part_reuse[*inc_reuse];
            //                (*inc_reuse)++;
            //            }
            //            else {
            new_ind = particles->Nb_part;
            particles->Nb_part++;
            //            }

            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh, &model );
        }

        //-------------------------------------------------------------------------------------------------------------------------//

        // Get the particle index
        new_x = mesh.xg_coord[ic] + mesh.dx/3.0;
        new_z = mesh.zg_coord[jc] - mesh.dz/3.0;
        if ( model.free_surface==1 ) h  = topo.a[cell]*new_x  + topo.b[cell];

        if ( new_x < model.xmax && new_z > model.zmin && new_z<h   ) {

            //            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
            //                new_ind = ind_part_reuse[*inc_reuse];
            //                (*inc_reuse)++;
            //            }
            //            else {
            new_ind = particles->Nb_part;
            particles->Nb_part++;
            //            }

            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh, &model );
        }

        //-------------------------------------------------------------------------------------------------------------------------//
        // Get the particle index
        new_x = mesh.xg_coord[ic] - mesh.dx/3.0;
        new_z = mesh.zg_coord[jc] + mesh.dz/3.0;
        if ( model.free_surface==1 ) h  = topo.a[cell]*new_x  + topo.b[cell];

        if ( new_z < model.zmax && new_x > model.xmin && new_z<h  ) {

            //            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
            //                new_ind = ind_part_reuse[(*inc_reuse)];
            //                (*inc_reuse)++;
            //            }
            //            else {
            new_ind = particles->Nb_part;
            particles->Nb_part++;
            //            }

            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh, &model );
        }

        //------------------------------------------------------------------------------------------------------------------------/
        // Get the particle index
        new_x = mesh.xg_coord[ic] + mesh.dx/3.0;
        new_z = mesh.zg_coord[jc] + mesh.dz/3.0;
        if ( model.free_surface==1 ) h  = topo.a[cell]*new_x  + topo.b[cell];


        if ( new_z < model.zmax && new_x < model.xmax && new_z<h ) {

            //            if ( (*inc_reuse) < nb_part_reuse && nb_part_reuse>0  ) {
            //                new_ind = ind_part_reuse[*inc_reuse];
            //                (*inc_reuse)++;
            //            }
            //            else {
            new_ind = particles->Nb_part;
            particles->Nb_part++;
            //            }

            particles->x[new_ind]      = new_x;
            particles->z[new_ind]      = new_z;
            FindClosestPhaseVertex( particles, ic, jc, mesh, ind_list, new_ind, nb_neigh, &model );
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void PutPartInBox( markers *particles, grid *mesh, params model, surface topo, scale scaling ) {

    // Set the original particle layout throughout the mesh.
    // The particles are set with regular spacing.
    int i, j, ki, kj, np;
    double dx_particles, dz_particles;
    int add_noise=model.initial_noise;
    double noise=0.1, random_x, random_z;
    
    // Compute the spacing between particles:
    dx_particles = mesh->dx/(double)(particles->Nx_part); // new
    dz_particles = mesh->dz/(double)(particles->Nz_part); // new
    printf("Initial particle spacing : dxm = %lf dzm = %lf m\n", dx_particles*scaling.L, dz_particles*scaling.L);
    
    // Loop over loop over loop over loop (on s'en branle, c'est du C)
    np = 0;
    double h, a, b, coord_x, coord_z;
    
    for (j=0; j<mesh->Nx-1; j++) {
        for (i=0; i<mesh->Nz-1; i++) {
            
            if (model.free_surface == 1) {
                a = topo.a[j];
                b = topo.b[j];
            }
            for (kj=0; kj<particles->Nx_part; kj++) {
                for (ki=0; ki<particles->Nz_part; ki++) {
                    
                    if (model.free_surface == 1) {
                        
//                        coord_x = mesh->xg_coord[j] + dx_particles*kj + dx_particles;
//                        coord_z = mesh->zg_coord[i] + dz_particles*ki + dz_particles;
                        
                        coord_x = mesh->xg_coord[j] + dx_particles*(kj+0.5);
                        coord_z = mesh->zg_coord[i] + dz_particles*(ki+0.5);
                        
                        h = b + a*coord_x;
                        
                        if ( coord_z < h ) {
                            particles->x[np]  = coord_x;
                            particles->z[np]  = coord_z;
                            
                            if ( add_noise == 1 ) {
                                
                                random_x = (double)rand()/(double)RAND_MAX;
                                random_z = (double)rand()/(double)RAND_MAX;
                                random_x = ((2*random_x)-1)*dx_particles*noise;
                                random_z = ((2*random_z)-1)*dz_particles*noise;
                                
                                particles->x[np]+= random_x;
                                particles->z[np]+= random_z;
                                
                                isoutPart( particles, &model, np );
                            }
            
                            np++;
                        }
                    }
                    
                    if (model.free_surface == 0) {
//                        coord_x = mesh->xg_coord[j] + dx_particles*kj + dx_particles;
//                        coord_z = mesh->zg_coord[i] + dz_particles*ki + dz_particles;
                        
                        coord_x = mesh->xg_coord[j] + dx_particles*(kj+0.5);
                        coord_z = mesh->zg_coord[i] + dz_particles*(ki+0.5);
                        particles->x[np]  = coord_x;
                        particles->z[np]  = coord_z;
                        
                        if ( add_noise == 1 ) {
                            
                            random_x = (double)rand()/(double)RAND_MAX;
                            random_z = (double)rand()/(double)RAND_MAX;
                            random_x = ((2*random_x)-1)*dx_particles*noise;
                            random_z = ((2*random_z)-1)*dz_particles*noise;
                            
                            particles->x[np]+= random_x;
                            particles->z[np]+= random_z;
                            
                            isoutPart( particles, &model, np );


                            
                        }
                        
                        np++;
                    }
                }
            }
        }
    }
    particles->Nb_part = np;
    printf("Initial number of particles = %d\n", particles->Nb_part);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double MarkerValue( DoodzFP* mat_prop, markers *particles, int k, int itp_type, int flag ) {
    
    double mark_val;
    
    if (flag==0) {
        mark_val = mat_prop[particles->phase[k]];
    }
    if (flag==1) {
        mark_val = mat_prop[k];
    }
    if (itp_type==1) {
        mark_val =  1.0/mark_val;
    }
    if (itp_type==2) {
        mark_val =  log(mark_val);
    }
    return mark_val;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Interp_Grid2P ( markers particles, DoodzFP* PartField, grid *mesh, double* NodeField, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag ) {

    int k;
    int i_part,j_part, iSW, iNW, iSE, iNE;
    double dx,dz,dxm,dzm;
    int  Nb_part;
    double distance, sumW;
    int periodix = 1;
//    int periodix = model->periodic_x;

    Nb_part=particles.Nb_part;
    double Xp, Lx = mesh->xg_coord[mesh->Nx-1] - mesh->xg_coord[0];
    int Nxg = mesh->Nx;
    
    dx=mesh->dx;
    dz=mesh->dz;
    
    int inSE, inSW, inNW, inNE;


#pragma omp parallel for shared ( particles, PartField, X_vect, Z_vect, NodeField, Nb_part, dx, dz, Nx, Nz, tag ) \
private ( k, dxm, dzm, j_part, i_part, distance, iSW, iNW, iSE, iNE, sumW, Xp ) firstprivate( Lx, Nxg )// schedule( static )
    for (k=0;k<Nb_part;k++) {

        PartField[k] = 0.0;

        // Filter out particles that are inactive (out of the box)
        if (particles.phase[k] != -1) {
            
            distance = (particles.x[k]-X_vect[0]);
            j_part   = ceil((distance/dx)) - 1;
            
            if (j_part<0) {
                j_part = 0;
            }

            if (j_part>Nx-2) {
                j_part = Nx-2;

                if (abs(Nx - (Nxg-1))>0) {
                    printf("j_part = %d --- Nx = %03d Nxc = %d --- px = %2.2e xg = %2.2e\n", j_part, Nx, Nxg-1, particles.x[k], X_vect[j_part]);
                    printf("i_part = %d --- Nz = %03d Nzc = %d --- pi = %2.2e ig = %2.2e\n", i_part, Nz, Nz-1, particles.z[k], Z_vect[i_part]);
                    exit(89);
                }
            }

            // Get the line:
            distance = (particles.z[k]-Z_vect[0]);
            i_part   = ceil((distance/dz)) - 1;
            if (i_part<0) {
                i_part = 0;
            }
            if (i_part>Nz-2) {
                i_part = Nz-2;
            }

            dxm = fabs(particles.x[k] - X_vect[j_part]);
            dzm = fabs(particles.z[k] - Z_vect[i_part]);

            iSW = j_part+i_part*Nx;
            iSE = j_part+i_part*Nx+1;
            iNW = j_part+(i_part+1)*Nx;
            iNE = j_part+(i_part+1)*Nx+1;

            PartField[k] = 0.0;
            sumW         = 0.0;
            
//            inSE = 0; inSW = 0; inNW = 0; inNE = 0;

            if (tag[iSW]!=30 && tag[iSW]!=31) {
                PartField[k] +=  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iSW];
                sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
            }
            else {
                PartField[k] +=  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iNW];
                sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
            }
            if (tag[iSE]!=30 && tag[iSE]!=31) {
                PartField[k] += (dxm/dx)  * (1.0-dzm/dz)  * NodeField[iSE];
                sumW += (dxm/dx)* (1.0-dzm/dz);
            }
            else {
                PartField[k] += (dxm/dx)  * (1.0-dzm/dz)  * NodeField[iNE];
                sumW += (dxm/dx)* (1.0-dzm/dz);
            }
            if (tag[iNW]!=30 && tag[iNW]!=31) {
                PartField[k] += (1.0-dxm/dx)* (dzm/dz)    * NodeField[iNW];
                sumW += (1.0-dxm/dx)* (dzm/dz);
            }
            else {
                PartField[k] += (1.0-dxm/dx)* (dzm/dz)    * NodeField[iSW];
                sumW += (1.0-dxm/dx)* (dzm/dz);
            }
            if (tag[iNE]!=30 && tag[iNE]!=31) {
                PartField[k] += (dxm/dx)  * (dzm/dz)    * NodeField[iNE];
                sumW += (dxm/dx)* (dzm/dz);
            }
            else {
                PartField[k] += (dxm/dx)  * (dzm/dz)    * NodeField[iSE];
                sumW += (dxm/dx)  * (dzm/dz) ;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Interp_Grid2P_centroids2 ( markers particles, DoodzFP* PartField, grid *mesh, double* NodeField, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag, params* model ) {

    int k;
    int i_part,j_part, iSW, iNW, iSE, iNE, ic, jc;
    double dx,dz,dxm,dzm;
    int  Nb_part;
    double distance, sumW;
    int periodix = model->periodic_x;

    Nb_part=particles.Nb_part;
    double Xp, Lx = mesh->xg_coord[mesh->Nx-1] - mesh->xg_coord[0];
    int Nxg = mesh->Nx;
    
    dx=mesh->dx;
    dz=mesh->dz;

#pragma omp parallel for shared ( particles, PartField, X_vect, Z_vect, NodeField, Nb_part, dx, dz, Nx, Nz, tag ) \
private ( k, dxm, dzm, j_part, i_part, distance, iSW, iNW, iSE, iNE, sumW, Xp, ic, jc ) firstprivate( Lx, Nxg, periodix )// schedule( static )
    for (k=0;k<Nb_part;k++) {

        PartField[k] = 0.0;

        // Filter out particles that are inactive (out of the box)
        if (particles.phase[k] != -1) {

            // Check distance with regards to extended array
            distance = (particles.x[k]-X_vect[0]);
            j_part   = ceil((distance/dx)) - 1;
            
            // Get the line:
            distance = (particles.z[k]-Z_vect[0]);
            i_part   = ceil((distance/dz)) - 1;
            
            if (j_part>0 && j_part<Nx && i_part>0 && i_part<Nz) {
                ic  = i_part - 1;
                jc  = j_part - 1;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx+1;
                iNW = jc + (ic+1)*Nx;
                iNE = jc + (ic+1)*Nx+1;
            }

            if (j_part==0 && i_part > 0 && i_part<Nz ) {
                if (periodix==1) {
                    ic  = i_part-1;
                    jc  = j_part;
                    iSW = (Nx-1) + ic*Nx;
                    iSE = jc     + ic*Nx;
                    iNW = (Nx-1) + (ic+1)*Nx;
                    iNE = jc     + (ic+1)*Nx;
                }
                else {
                    // Copy solution along edge
                    ic  = i_part-1;
                    jc  = j_part;
                    iSW = jc + ic*Nx;
                    iSE = jc + ic*Nx;
                    iNW = jc + (ic+1)*Nx;
                    iNE = jc + (ic+1)*Nx;
                }
            }

            if (j_part==Nx && i_part > 0 && i_part<Nz ) {
                if (periodix==1) {
                    ic  = i_part-1;
                    jc  = j_part-1;
                    iSW = jc + ic*Nx;
                    iSE = 0  + ic*Nx;
                    iNW = jc + (ic+1)*Nx;
                    iNE = 0  + (ic+1)*Nx;
                }
                else {
                    // Copy solution along edge
                    ic  = i_part-1;
                    jc  = j_part-1;
                    iSW = jc + ic*Nx;
                    iSE = jc + ic*Nx;
                    iNW = jc + (ic+1)*Nx;
                    iNE = jc + (ic+1)*Nx;
                }
            }


            if (i_part==0 && j_part > 0 && j_part<Nx) {
                // Copy solution along edge
                ic  = i_part;
                jc  = j_part-1;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx + 1;
                iNW = jc + ic*Nx;
                iNE = jc + ic*Nx + 1;
            }

            if (i_part==Nz && j_part > 0 && j_part<Nx) {
                // Copy solution along edge
                ic  = i_part-1;
                jc  = j_part-1;
                iSW = jc + ic*Nx;
                iSE = jc + ic*Nx + 1;
                iNW = jc + ic*Nx;
                iNE = jc + ic*Nx + 1;
            }

            // SW
            if (i_part==0 && j_part==0) {

                if (periodix==1) {
                    ic  = 0;
                    jc  = 0;
                    iSW = (Nx-1) + ic*Nx;
                    iSE = jc     + ic*Nx;
                    iNW = (Nx-1) + ic*Nx;
                    iNE = jc     + ic*Nx;
                }
                else {
                    ic  = i_part - 0;
                    jc  = j_part - 0;
                    iSW = jc + ic*Nx;
                    iSE = jc + ic*Nx;
                    iNW = jc + ic*Nx;
                    iNE = jc + ic*Nx;
                }
            }
            // SE
            if (i_part==0 && j_part==Nx) {

                if (periodix==1) {
                    ic  = i_part - 0;
                    jc  = j_part - 1;
                    iSW = jc + ic*Nx;
                    iSE = 0  + ic*Nx;
                    iNW = jc + ic*Nx;
                    iNE = 0  + ic*Nx;
                }
                else {
                    ic  = i_part - 0;
                    jc  = j_part - 1;
                    iSW = jc + ic*Nx;
                    iSE = jc + ic*Nx;
                    iNW = jc + ic*Nx;
                    iNE = jc + ic*Nx;
                }
            }
            // NW
            if (i_part==Nz && j_part==0) {

                if (periodix==1) {
                    ic  = i_part - 1;
                    jc  = j_part - 0;
                    iSW = (Nx-1) + ic*Nx;
                    iSE = jc     + ic*Nx;
                    iNW = (Nx-1) + ic*Nx;
                    iNE = jc     + ic*Nx;
                }
                else {
                    ic  = i_part - 1;
                    jc  = j_part - 0;
                    iSW = jc + ic*Nx;
                    iSE = jc + ic*Nx;
                    iNW = jc + ic*Nx;
                    iNE = jc + ic*Nx;
                }
            }
            // NE
            if (i_part==Nz && j_part==Nx) {

                if (periodix==1) {
                    ic  = i_part - 1;
                    jc  = j_part - 1;
                    iSW = jc + ic*Nx;
                    iSE = 0  + ic*Nx;
                    iNW = jc + ic*Nx;
                    iNE = 0  + ic*Nx;
                }
                else {
                    ic  = i_part - 1;
                    jc  = j_part - 1;
                    iSW = jc + ic*Nx;
                    iSE = jc + ic*Nx;
                    iNW = jc + ic*Nx;
                    iNE = jc + ic*Nx;
                }
            }
        
            dxm = (particles.x[k] - X_vect[j_part]);
            dzm = (particles.z[k] - Z_vect[i_part]);
            
            PartField[k] = 0.0;
            sumW         = 0.0;
            
            if (tag[iSW]!=30 && tag[iSW]!=31) {
                PartField[k] +=  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iSW];
                sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
            }
            else {
                PartField[k] +=  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iNW];
                sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
            }
            if (tag[iSE]!=30 && tag[iSE]!=31) {
                PartField[k] += (dxm/dx)  * (1.0-dzm/dz)  * NodeField[iSE];
                sumW += (dxm/dx)* (1.0-dzm/dz);
            }
            else {
                PartField[k] += (dxm/dx)  * (1.0-dzm/dz)  * NodeField[iNE];
                sumW += (dxm/dx)* (1.0-dzm/dz);
            }
            if (tag[iNW]!=30 && tag[iNW]!=31) {
                PartField[k] += (1.0-dxm/dx)* (dzm/dz)    * NodeField[iNW];
                sumW += (1.0-dxm/dx)* (dzm/dz);
            }
            else {
                PartField[k] += (1.0-dxm/dx)* (dzm/dz)    * NodeField[iSW];
                sumW += (1.0-dxm/dx)* (dzm/dz);
            }
            if (tag[iNE]!=30 && tag[iNE]!=31) {
                PartField[k] += (dxm/dx)  * (dzm/dz)    * NodeField[iNE];
                sumW += (dxm/dx)* (dzm/dz);
            }
            else {
                PartField[k] += (dxm/dx)  * (dzm/dz)    * NodeField[iSE];
                sumW += (dxm/dx)  * (dzm/dz) ;
            }

//            if (tag[iSW]!=30 && tag[iSW]!=31) {
//                PartField[k] +=  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iSW];
//                sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
//
//            }
//            if (tag[iSE]!=30 && tag[iSE]!=31) {
//                PartField[k] += (dxm/dx)  * (1.0-dzm/dz)  * NodeField[iSE];
//                sumW += (dxm/dx)* (1.0-dzm/dz);
//            }
//            if (tag[iNW]!=30 && tag[iNW]!=31) {
//                PartField[k] += (1.0-dxm/dx)* (dzm/dz)    * NodeField[iNW];
//                sumW += (1.0-dxm/dx)* (dzm/dz);
//            }
//            if (tag[iNE]!=30 && tag[iNE]!=31) {
//                PartField[k] += (dxm/dx)  * (dzm/dz)    * NodeField[iNE];
//                sumW += (dxm/dx)* (dzm/dz);
//            }

            if (sumW>1e-13) PartField[k] /= sumW;
            if (isinf(PartField[k])) {
                printf("%2.2e %2.2e %2.2e %2.2e \n", NodeField[iSW],NodeField[iSE], NodeField[iNW],NodeField[iNE]);
                printf("%2.2e\n", sumW);
                exit(1);
            }
        }
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void Interp_Grid2P_strain( markers particles, DoodzFP* PartField, grid *mesh, double* NodeField, double* X_vect, double* Z_vect, int Nx, int Nz, char *tag ) {

    int    k, i_part, j_part, iSW, iNW, iSE, iNE;
    double dx,dz,dxm,dzm;
    int    Nb_part;
    double distance, sumW;

    Nb_part=particles.Nb_part;

    dx=mesh->dx;
    dz=mesh->dz;

#pragma omp parallel for shared ( particles, PartField, X_vect, Z_vect, NodeField, Nb_part, dx, dz, Nx, Nz, tag ) \
private ( k, dxm, dzm, j_part, i_part, distance, iSW, iNW, iSE, iNE, sumW ) // schedule( static )
    for (k=0;k<Nb_part;k++) {

        PartField[k] = 0.0;

        // Filter out particles that are inactive (out of the box)
        if (particles.phase[k] != -1) {

            distance = fabs(particles.x[k]-X_vect[0]);
            j_part   = ceil((distance/dx)) - 1;
            if (j_part<0) {
                j_part = 0;
            }
            if (j_part>Nx-2) {
                j_part = Nx-2;
            }

            // Get the line:
            distance = fabs(particles.z[k]-Z_vect[0]);
            i_part   = ceil((distance/dz)) - 1;
            if (i_part<0) {
                i_part = 0;
            }
            if (i_part>Nz-2) {
                i_part = Nz-2;
            }

            dxm = (particles.x[k] - X_vect[j_part]);
            dzm = (particles.z[k] - Z_vect[i_part]);

            iSW = j_part+i_part*Nx;
            iSE = j_part+i_part*Nx+1;
            iNW = j_part+(i_part+1)*Nx;
            iNE = j_part+(i_part+1)*Nx+1;

            PartField[k] = 0.0;
            sumW         = 0.0;

            if (j_part>0 && j_part<Nx-2 && i_part>0 && i_part<Nz-2) {

                if (tag[iSW]!=30 && tag[iSW]!=31) {
                    PartField[k] +=  (1.0-dxm/dx)* (1.0-dzm/dz) * NodeField[iSW];
                    sumW += (1.0-dxm/dx)* (1.0-dzm/dz);
                }
                if (tag[iSE]!=30 && tag[iSE]!=31) {
                    PartField[k] += (dxm/dx)  * (1.0-dzm/dz)  * NodeField[iSE];
                    sumW += (dxm/dx)* (1.0-dzm/dz);
                }
                if (tag[iNW]!=30 && tag[iNW]!=31) {
                    PartField[k] += (1.0-dxm/dx)* (dzm/dz)    * NodeField[iNW];
                    sumW += (1.0-dxm/dx)* (dzm/dz);
                }
                if (tag[iNE]!=30 && tag[iNE]!=31) {
                    PartField[k] += (dxm/dx)  * (dzm/dz)    * NodeField[iNE];
                    sumW += (dxm/dx)* (dzm/dz);
                }
            }

            if (j_part==0    && i_part==0   ) PartField[k] =  NodeField[iNE];
            if (j_part==Nx-2 && i_part==0   ) PartField[k] =  NodeField[iNW];
            if (j_part==0    && i_part==Nz-2) PartField[k] =  NodeField[iSE];
            if (j_part==Nx-2 && i_part==Nz-2) PartField[k] =  NodeField[iSW];

            if (j_part==0 && i_part>0 && i_part<Nz-2) {

                if (tag[iSE]!=30 && tag[iSE]!=31) {
                    PartField[k] +=  (1.0-dzm/dz)  * NodeField[iSE];
                    sumW +=  (1.0-dzm/dz);
                }
                if (tag[iNE]!=30 && tag[iNE]!=31) {
                    PartField[k] +=  (dzm/dz)    * NodeField[iNE];
                    sumW +=  (dzm/dz);
                }
                if(sumW>1e-13) PartField[k] /= sumW;
            }

            if (j_part==Nx-2 && i_part>0 && i_part<Nz-2) {

                if (tag[iSW]!=30 && tag[iSW]!=31) {
                    PartField[k] +=  (1.0-dzm/dz) * NodeField[iSW];
                    sumW +=  (1.0-dzm/dz);
                }
                if (tag[iNW]!=30 && tag[iNW]!=31) {
                    PartField[k] +=  (dzm/dz)    * NodeField[iNW];
                    sumW +=  (dzm/dz);
                }
                if(sumW>1e-13) PartField[k] /= sumW;
            }

            if (i_part==0 && j_part>0 && j_part<Nx-2) {
                if (tag[iNW]!=30 && tag[iNW]!=31) {
                    PartField[k] += (1.0-dxm/dx)   * NodeField[iNW];
                    sumW += (1.0-dxm/dx);
                }
                if (tag[iNE]!=30 && tag[iNE]!=31) {
                    PartField[k] += (dxm/dx)      * NodeField[iNE];
                    sumW += (dxm/dx);
                }
                if(sumW>1e-13) PartField[k] /= sumW;
            }
            if (i_part==Nz-2 && j_part>0 && j_part<Nx-2) {
                if (tag[iSW]!=30 && tag[iSW]!=31) {
                    PartField[k] +=  (1.0-dxm/dx) * NodeField[iSW];
                    sumW += (1.0-dxm/dx);
                }
                if ( tag[iSE] != 30 && tag[iSE] != 31 ) {
                    PartField[k] += (dxm/dx)   * NodeField[iSE];
                    sumW += (dxm/dx);
                }
                if(sumW>1e-13) PartField[k] /= sumW;
            }
            if(sumW>1e-13) PartField[k] /= sumW;

        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Particles to reference nodes
void Interp_Phase2VizGrid ( markers particles, int* PartField, grid *mesh, char* NodeField, double* X_vect, double* Z_vect, int Nxg, int Nzg, params model, surface topo ) {

    int     k, i_part, j_part, Nb_part, Nx=Nxg-1, Nz=Nzg-1, ncx=mesh->Nx-1, ncz=mesh->Nz-1, ic, jc, flag;
    double  dx, dz, dxm, dzm, *Wm, distance;
    Nb_part=particles.Nb_part;
    dx=X_vect[1]-X_vect[0];
    dz=Z_vect[1]-Z_vect[0];

    double dx_std = model.dx;
    double dz_std = model.dz;

    double xmin = mesh->xg_coord[0] + dx_std/2.0;


    Wm   = DoodzCalloc (Nx*Nz,sizeof(double));    // allocate storage for the array

    //--------------------------------------------------------------
    // Compute Wm and BmWm
    //--------------------------------------------------------------

    // Set default values to -1
    for (k=0;k<Nz*Nx;k++) {
        NodeField[k] = -1;
    }

    // Get the column:
    distance         = (particles.x[k] - xmin);
    ic               = ceil((distance/dx_std)+0.5) - 1;
    if (ic<0)     ic = 0;
    if (ic>ncx-1) ic = ncx-1;

    // Compute topography
//    if ( model.free_surface == 1 ) h = topo.b[ic] + topo.a[ic]*particles.x[k];

    for (k=0;k<Nb_part;k++) {

        // Filter out particles that are inactive (out of the box)
        if ( particles.phase[k] != -1  ) { //&& particles.z[k] < h

            // Get the column:
            distance = ( particles.x[k] - (mesh->xg_coord[0]+0.5*dx) );
            j_part   = ceil( (distance/dx) + 0.5) - 1;

            if (j_part<0)    j_part = 0;
            if (j_part>Nx-1) j_part = Nx-1;

            // Get the line:
            distance = ( particles.z[k] - (mesh->zg_coord[0]+0.5*dz) );
            i_part   = ceil( (distance/dz) + 0.5) - 1;

            if (i_part<0)    i_part = 0;
            if (i_part>Nz-1) i_part = Nz-1;

            dxm = fabs(0.5*(X_vect[j_part]  + X_vect[j_part+1]) - particles.x[k]);
            dzm = fabs(0.5*(Z_vect[i_part+1]+ Z_vect[i_part])   - particles.z[k]);



            // Find index of the particle with regard to the computation grid
            // Then we can find the flag of the cell and figure out in we are below the free surface or not
            // x index
            distance         = particles.x[k] - mesh->xc_coord[0];
            ic               = ceil((distance/dx_std)+0.5) -1;
            if (ic<0)     ic = 0;
            if (ic>ncx-1) ic = ncx-1;

            // z index
            distance         = particles.z[k] - mesh->zc_coord[0];
            jc               = ceil((distance/dz_std)+0.5) -1;
            if (jc<0)     jc = 0;
            if (jc>ncz-1) jc = ncz-1;

            flag = mesh->BCp.type[ic+jc*(ncx)];

//            printf("ic %d jc %d \n", ic, jc);
//            if (flag==30 || flag==31) printf("OKAY\n");


            if ( (flag != 30 && flag != 31) && (1.0-(dxm/dx))*(1.0-(dzm/dz)) >  Wm[j_part+i_part*(Nx)] ) {
                Wm[j_part+i_part*(Nx)]   = (1.0-(dxm/dx))*(1.0-(dzm/dz)) ;

                NodeField[j_part+i_part*(Nx)] = (char)PartField[k];

//                // Update cell is flag indicate below free surface
//                if (flag != 30 && flag != 31) NodeField[j_part+i_part*(Nx)] = (char)PartField[k];
//                else NodeField[j_part+i_part*(Nx)] = -1;
//
//                if (flag == 30 || flag == 31) NodeField[j_part+i_part*(Nx)] = -1;
            }

        }
    }
    DoodzFree(Wm);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Particles to reference nodes
void P2Mastah ( params *model, markers particles, DoodzFP* mat_prop, grid *mesh, double* NodeField, char* NodeType, int flag, int avg, int prop, int centroid, int interp_stencil) {
    
    // flag == 0 --> interpolate from material properties structure
    // flag == 1 --> interpolate straight from the particle arrays
    // flag ==-1 --> specific for anisotropy (see Anisotropy_v2.ipynb)
    // flag ==-2 --> specific for anisotropy (see Anisotropy_v2.ipynb)
    // flag ==-3 --> specific for anisotropy: angle of director
    // avg  == 0 --> arithmetic distance-weighted average
    // avg  == 1 --> harmonic distance-weighted average
    // avg  == 2 --> geometric distance-weighted average
    
    int p, i, k, kp, ip, jp, Np, nthreads, thread_num, Nx, Nz;
    double dx, dz, dxm, dzm,  distance, mark_val, dx_itp, dz_itp;
    double *WM, *BMWM;
    double **Wm, **BmWm, ***Wm_ph;
    double *X_vect, *Z_vect;
    int peri;
    
    if (centroid==1) {
        Nx = mesh->Nx-1;
        Nz = mesh->Nz-1;
        X_vect = mesh->xc_coord;
        Z_vect = mesh->zc_coord;
    }
    if (centroid==0) {
        Nx = mesh->Nx;
        Nz = mesh->Nz;
        X_vect = mesh->xg_coord;
        Z_vect = mesh->zg_coord;
    }
    if (centroid==-1) {           // Vx nodes
        Nx = mesh->Nx;
        Nz = mesh->Nz+1;
        X_vect = mesh->xg_coord;
        Z_vect = mesh->zvx_coord;
    }
    if (centroid==-2) {           // Vy nodes
        Nx = mesh->Nx+1;
        Nz = mesh->Nz;
        X_vect = mesh->xvz_coord;
        Z_vect = mesh->zg_coord;
    }

    Np = particles.Nb_part;
    dx = mesh->dx;
    dz = mesh->dz;
    if (interp_stencil==1) dx_itp =     dx/2.0; // 1-cell
    if (interp_stencil==1) dz_itp =     dz/2.0; // 1-cell
    if (interp_stencil==9) dx_itp = 3.0*dx/2.0; // 9-cell
    if (interp_stencil==9) dz_itp = 3.0*dz/2.0; // 9-cell
    
    // Initialisation
    if ( prop == 1 ) {
#pragma omp parallel for shared ( mesh ) private( i, k, p ) firstprivate( Nx, Nz, nthreads, model, centroid ) schedule( static )
        for ( i=0; i<Nx*Nz; i++ ) {
            for (p=0; p<model->Nb_phases; p++) {
                if (centroid==0) mesh->phase_perc_s[p][i] = 0.0;
                if (centroid==1) mesh->phase_perc_n[p][i] = 0.0;
            }
        }
    }
    
#pragma omp parallel
    {
        nthreads = omp_get_num_threads();
    }
    
    //--------------------------------------------------------------
    // Initialize Wm and BmWm
    //--------------------------------------------------------------
    Wm    = DoodzCalloc ( nthreads, sizeof(double*));    // allocate storage for the array
    BmWm  = DoodzCalloc ( nthreads, sizeof(double*));    // allocate storage for the array
    Wm_ph = DoodzCalloc ( nthreads, sizeof(double**));
    
    for ( k=0; k<nthreads; k++ ) {
        Wm[k]    = DoodzCalloc ( Nx*Nz, sizeof(double));
        BmWm[k]  = DoodzCalloc ( Nx*Nz, sizeof(double));
        Wm_ph[k] = DoodzCalloc ( model->Nb_phases, sizeof(double*));
        for (p=0; p<model->Nb_phases; p++) Wm_ph[k][p] = DoodzCalloc ( Nx*Nz, sizeof(double));
    }
    
    WM   = DoodzCalloc ( Nx*Nz, sizeof(double));
    BMWM = DoodzCalloc ( Nx*Nz, sizeof(double));
    
    //--------------------------------------------------------------
    // Compute Wm and BmWm
    //--------------------------------------------------------------
    
    #pragma omp parallel for shared ( particles, BmWm, Wm, Wm_ph, X_vect, Z_vect )       \
    private ( k, kp, dxm, dzm, ip, jp, distance, mark_val, thread_num, p, peri   )       \
    firstprivate ( mat_prop, dx, dz, Np, Nx, Nz, mesh, flag, avg, dx_itp, dz_itp, prop)  schedule( static )
        
        for (k=0; k<Np; k++) {
            
            // Filter out particles that are inactive (out of the box)
            if (particles.phase[k] != -1) {
                
                thread_num = omp_get_thread_num();
                p          = particles.phase[k];
                
                // Get the column:
                distance = ( particles.x[k] - X_vect[0] );
                ip       = ceil( (distance/dx) + 0.5) - 1;
                if (ip<0   ) ip = 0;
                if (ip>Nx-1) ip = Nx-1;

                // Get the line:
                distance = ( particles.z[k] - Z_vect[0] );
                jp       = ceil( (distance/dz) + 0.5) - 1;
                if (jp<0   ) jp = 0;
                if (jp>Nz-1) jp = Nz-1;

                // Distance to node
                dxm = fabs( X_vect[ip] - particles.x[k]);
                dzm = fabs( Z_vect[jp] - particles.z[k]);
                
                if ( prop == 0 ) {
                    // Layer orientation
                    double nx, nz;
                    // Get material properties (from particules or mat_prop array)
                    switch (flag) {
                        case 0:
                            mark_val = mat_prop[particles.phase[k]];
                            break;
                            
                        case 1:
                            mark_val = mat_prop[k];
                            break;
                            
                        case -1: // Anisotropy: Anisotropy_v2.ipynb
                            nx = particles.nx[k];
                            nz = particles.nz[k];
                            mark_val =  2.0*pow(nx, 2.0)*pow(nz, 2.0);
                            break;
                            
                        case -2: // Anisotropy: Anisotropy_v2.ipynb
                            nx = particles.nx[k];
                            nz = particles.nz[k];
                            mark_val = nx*nz*(-pow(nx, 2.0) + pow(nz, 2.0));
                            break;

                        case -3: // Anisotropy: angle of director
                            mark_val = acos(particles.nx[k]);
                            break;
                            
                        default:
                            break;
                    }
                    //----------------------
                    if (avg==1) {
                        mark_val =  1.0/mark_val;
                    }
                    if (avg==2) {
                        mark_val =  log(mark_val);
                    }
                }
            
                // Center
                kp = ip + jp*Nx;
                Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                
                // --------------------------
                
                // Case for 9 Cell
                if ( interp_stencil==9 ) {
                
                    // N
                    if ( jp<Nz-1 ){
                        kp  = ip + (jp+1)*Nx;
                        dxm = fabs( X_vect[ip]      - particles.x[k]);
                        dzm = fabs( Z_vect[jp] + dz - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }

                    // S
                    if ( jp>0 ){
                        kp  = ip + (jp-1)*Nx;
                        dxm = fabs( X_vect[ip]      - particles.x[k]);
                        dzm = fabs( Z_vect[jp] - dz - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }

                    // E
                    peri = 0;
                    if (ip==Nx-1 && model->periodic_x==1) peri = 1;

                    if ( ip<Nx-1 || peri==1 ){
                        kp  = ip + jp*Nx + (1.0-peri)*1 - peri*(Nx-1);
                        dxm = fabs( X_vect[ip] + dx - particles.x[k]);
                        dzm = fabs( Z_vect[jp]      - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }

                    // W
                    peri = 0;
                    if (ip==0 && model->periodic_x==1) peri = 1;

                    if ( ip>0 || peri==1 ){
                        kp  = ip + jp*Nx - (1.0-peri)*1 + peri*(Nx-1);
                        dxm = fabs( X_vect[ip] - dx - particles.x[k]);
                        dzm = fabs( Z_vect[jp]      - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }
                    
                    // --------------------------
                    
                    // SW
                    peri = 0;
                    if (ip==0 && model->periodic_x==1) peri = 1;

                    if ( jp>0 && (ip>0 || peri==1) ){
                        kp = ip + (jp-1)*Nx - (1.0-peri)*1 + peri*(Nx-1);
                        dxm = fabs( X_vect[ip] - dx - particles.x[k]);
                        dzm = fabs( Z_vect[jp] - dz - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }

                    // NW
                    peri = 0;
                    if (ip==0 && model->periodic_x==1) peri = 1;

                    if ( jp<Nz-1 && (ip>0 || peri==1) ){
                        kp = ip + (jp+1)*Nx - (1.0-peri)*1 + peri*(Nx-1);
                        dxm = fabs( X_vect[ip] - dx - particles.x[k]);
                        dzm = fabs( Z_vect[jp] + dz - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }

                    // NE
                    peri = 0;
                    if (ip==Nx-1 && model->periodic_x==1) peri = 1;

                    if ( jp<Nz-1 && (ip<Nx-1 || peri==1) ){
                        kp = ip + (jp+1)*Nx + (1.0-peri)*1 - peri*(Nx-1);
                        dxm = fabs( X_vect[ip] + dx - particles.x[k]);
                        dzm = fabs( Z_vect[jp] + dz - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }

                    // SE
                    peri = 0;
                    if (ip==Nx-1 && model->periodic_x==1) peri = 1;

                    if ( jp>0 && (ip<Nx-1 || peri==1) ){
                        kp = ip + (jp-1)*Nx + (1.0-peri)*1 - peri*(Nx-1);
                        dxm = fabs( X_vect[ip] + dx - particles.x[k]);
                        dzm = fabs( Z_vect[jp] - dz - particles.z[k]);
                        Wm_ph[thread_num][p][kp]  +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        Wm[thread_num][kp]        +=          (1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                        BmWm[thread_num][kp]      += mark_val*(1.0-dxm/dx_itp)*(1.0-dzm/dz_itp);
                    }
           
                }
                
            }
        }

    // Final reduction
#pragma omp parallel for shared ( BmWm, Wm, BMWM, WM, Wm_ph, mesh ) private( i, k, p ) firstprivate( Nx, Nz, nthreads, model, centroid, prop) schedule( static )
    for ( i=0; i<Nx*Nz; i++ ) {
        for ( k=0; k<nthreads; k++ ) {
            WM[i]   += Wm[k][i];
            BMWM[i] += BmWm[k][i];
            if ( prop == 1 ) {
                for (p=0; p<model->Nb_phases; p++) {
                    if (centroid==0) mesh->phase_perc_s[p][i] += Wm_ph[k][p][i];
                    if (centroid==1) mesh->phase_perc_n[p][i] += Wm_ph[k][p][i];
                }
            }
        }
    }
    
//    for ( i=0; i<Nx*Nz; i++ ) {
//        printf("%lf\n", mesh->phase_perc_n[1][i]);
//    }

    //--------------------------------------------------------------
    // Get interpolated value on nodes
    //--------------------------------------------------------------
    
    
#pragma omp parallel for shared ( NodeField, BMWM, WM, Nx, Nz, mesh, NodeType ) private( i, p ) firstprivate ( avg ) schedule( static )
    for (i=0;i<Nx*Nz;i++) {
        
        if ( fabs(WM[i])<1e-30  || (NodeType[i]==30 || NodeType[i]==31) ) {
        }
        else {
            
            if ( prop == 0 ) {
                
                NodeField[i] = BMWM[i]/WM[i];
                if (avg==1) {
                    NodeField[i] =  1.0 / NodeField[i];
                }
                if (avg==2) {
                    NodeField[i] =  exp(NodeField[i]);
                }
            }
            if ( prop == 1 ) {
                for (p=0; p<model->Nb_phases; p++) {
                    if (centroid==0) mesh->phase_perc_s[p][i] /= WM[i];
                    if (centroid==1) mesh->phase_perc_n[p][i] /= WM[i];
                }
            }
        }
    }

    // Make periodic
    if (centroid==0 && prop==0 && model->periodic_x==1) {
        double av;
        for( int l=0; l<Nz; l++) {
            int c1 = l*Nx + Nx-1;
            av = 0.5*(NodeField[c1] + NodeField[l*Nx]);
            NodeField[c1] = av; NodeField[l*Nx] = av;
        }
    }
    
    // Clean up
    DoodzFree(WM);
    DoodzFree(BMWM);
    
    for ( k=0; k<nthreads; k++ ) {
        for (p=0; p<model->Nb_phases; p++) DoodzFree(Wm_ph[k][p]);
        DoodzFree(Wm[k]);
        DoodzFree(BmWm[k]);
        DoodzFree(Wm_ph[k]);
    }
    DoodzFree(Wm_ph);
    DoodzFree(Wm);
    DoodzFree(BmWm);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CountPartCell( markers* particles, grid *mesh, params model, surface topo, surface topo_ini, int reseed_markers, scale scaling ) {

    // First operation - compute phase proportions on centroids and vertices
    int cent=1, vert=0, prop=1, interp=0;
    P2Mastah ( &model, *particles, NULL, mesh, NULL, mesh->BCp.type,  0, 0, prop, cent, model.interp_stencil);
    P2Mastah ( &model, *particles, NULL, mesh, NULL, mesh->BCg.type,  0, 0, prop, vert, model.interp_stencil);

    // Let's consider a finer mesh of double resolution shifted by -dx/4 and -dz/4, this way both centroids and vertices are coverer
    const int nvx = 2*mesh->Nx-1, ncx = 2*mesh->Nx-2;  
    const int nvz = 2*mesh->Nz-1, ncz = 2*mesh->Nz-2;
    const double dx = model.dx/2.0;
    const double dz = model.dz/2.0;
    const double xmin = model.xmin, xmax = model.xmax;
    const double zmin = model.zmin, zmax = model.zmax;

    // (1) Let's count the number of active particle in each finer mesh cell
    int Nb_part = particles->Nb_part, ic, jc, kc, nb_part_reuse = 0, icoarse, jcoarse, Ncx=model.Nx-1, Ncz=model.Nz-1;
    double distance, x, z;

    // One array will store the number of particle in each of the cells
    int* nb_part_cell = DoodzCalloc(     ncx*ncz, sizeof(int));
    int* part_reuse   = DoodzCalloc(     Nb_part, sizeof(int));
    // char* Cell_flags  = DoodzCalloc(     ncx*ncz, sizeof(int));

    // Copy flags from coarse to fine mesh
    for (ic=0; ic<ncx; ic++) {
        for (jc=0; jc<ncz; jc++) {
            kc = ic + jc*ncx;

            // Get the column:
            x = xmin + dx/2.0 + ic*dx;
            distance = (x - (xmin + model.dx/2) );
            icoarse  = ceil( (distance/model.dx) + 0.5) - 1;
            if (icoarse<0    ) icoarse = 0;
            if (icoarse>Ncx-1) icoarse = Ncx-1;

            // Get the line:
            z = zmin + dz/2.0 + jc*dz;
            distance = (z - (zmin + model.dz/2));
            jcoarse  = ceil( (distance/model.dz) + 0.5) - 1;
            if (jcoarse<0    ) jcoarse = 0;
            if (jcoarse>Ncz-1) jcoarse = Ncz-1;

            // Cell_flags[kc] = mesh->BCt.type[icoarse + jcoarse*Ncx];
            mesh->BCt_fine.type[kc] = mesh->BCt.type[icoarse + jcoarse*Ncx];

        }
    }

    // Count particles in cells
    for (int k=0; k<Nb_part; k++) {
        if ( particles->phase[k] != -1 )  {
            // Get the column:
            distance = ( particles->x[k] - (xmin + dx/2) );
            ic       = ceil( (distance/dx) + 0.5) - 1;
            if (ic<0    ) ic = 0;
            if (ic>ncx-1) ic = ncx-1;

            // Get the line:
            distance = ( particles->z[k] - (zmin + dz/2) );
            jc       = ceil( (distance/dz) + 0.5) - 1;
            if (jc<0    ) jc = 0;
            if (jc>ncz-1) jc = ncz-1;

            // Center
            kc = ic + jc*ncx;

            // Add particle to finer mesh cell
            nb_part_cell[kc] += 1;
        }
        else {
            part_reuse[nb_part_reuse] = k;
            nb_part_reuse++;
        }
    }

    // // TBD
    //     int should_fill = 0;
    //     for (ic=0; ic<ncx; ic++) {
    //         for (jc=0; jc<ncz; jc++) {
    //             kc = ic + jc*ncx;
    //             if (nb_part_cell[kc] < 2 && Cell_flags[kc]<30) {
    //                 should_fill++;
    //                 printf("should_fill %d i=%d j=%d, nc = %d\n", should_fill, ic, jc, nb_part_cell[kc]);
    //             }
    //         }
    //     }

    if (reseed_markers==1) {

        // Allocate memory for an array that list particle indices in each cell
        int** part_cell = DoodzCalloc(ncx*ncz, sizeof(int*));
        for (kc=0; kc<ncx*ncz; kc++) {
            part_cell[kc] = DoodzCalloc(nb_part_cell[kc], sizeof(int));
        }

        // Store indices
        int kp;
        int* part_count = DoodzCalloc(ncx*ncz, sizeof(int));
        for (int k=0; k<Nb_part; k++) {
            if ( particles->phase[k] != -1 )  {

                // Get the column:
                distance = ( particles->x[k] - (xmin + dx/2) );
                ic       = ceil( (distance/dx) + 0.5) - 1;
                if (ic<0    ) ic = 0;
                if (ic>ncx-1) ic = ncx-1;

                // Get the line:
                distance = ( particles->z[k] - (zmin + dz/2) );
                jc       = ceil( (distance/dz) + 0.5) - 1;
                if (jc<0    ) jc = 0;
                if (jc>ncz-1) jc = ncz-1;

                // Center
                kc = ic + jc*ncx;

                part_cell[kc][part_count[kc]] = k;
                part_count[kc] += 1;
            }
        }
        DoodzFree(part_count);

        // Identify cells that need a new particle
        int nb_new_parts = 0, new_ind;
        int count_part_reuse = 0, count_created = 0;
        double new_x, new_z, h = model.zmax;

        for (ic=0; ic<ncx; ic++) {
            for (jc=0; jc<ncz; jc++) {
                kc = ic + jc*ncx;
                
                h = model.free_surface == 1 ? topo.height_finer_c[ic] : model.zmax+model.dz/2;
                double xc = xmin + dx/2 + ic*dx; // centroid of fine grid cell
                double zc = zmin + dz/2 + jc*dz; // centroid of fine grid cell

                if (nb_part_cell[kc] < 2 && mesh->BCt_fine.type[kc] != 30) {  // Cell_flags[kc]<30
                    // printf("Adding... i=%d j=%d, nb = %d\n", ic, jc, nb_part_cell[kc]);
                    int imin = ic>0     ? ic-1 : 0;
                    int imax = ic<ncx-1 ? ic+1 : ncx-1;
                    int jmin = jc>0     ? jc-1 : 0;
                    int jmax = jc<ncz-1 ? jc+1 : ncz-1;

                    if (ic>1    ) imin = ic-2; 
                    if (ic<ncx-2) imax = ic+2;
                    if (jc>1    ) jmin = jc-2; 
                    if (jc<ncz-2) jmax = jc+2;

                   // Number of particles in neigbour cells
                    int nb_neigh = 0;
                    for (int i=imin; i<=imax; i++) {
                        for (int j=jmin; j<=jmax; j++) {
                            int kk = i + j*ncx;
                            nb_neigh += nb_part_cell[kk];
                        }
                    }
                    if (nb_neigh==0) {
                        printf("Cell i=%d j=%d,flag=%d has 0 neighouring particles... Exiting\nimin=%d imax=%d jmin=%d jmax=%d\n", ic, jc, mesh->BCt_fine.type[kc], imin, imax, jmin, jmax);
                        for (int i=imin; i<=imax; i++) {
                            for (int j=jmin; j<=jmax; j++) {
                                printf("(%d,%d): %d", i, j, nb_part_cell[i+j*ncx]);
                            }
                        }
                        printf("xmin = %lf - xmax = %lf --- zmin = %lf - zmax = %lf\n", xmin-dx/2.0 + imin*dx, xmin-dx/2.0 + imax*dx, xmin-dz/2.0 + jmin*dz, xmin-dz/2.0 + jmax*dz );
                        exit(122);
                    }
                    // Collect particle indices in an array
                    int* neighbours = DoodzCalloc(nb_neigh, sizeof(int));
                    nb_neigh = 0;
                    for (int i=imin; i<=imax; i++) {
                        for (int j=jmin; j<=jmax; j++) {
                            kc = i + j*ncx;
                            for (int n=0; n<part_count[kc]; n++) {
                                neighbours[nb_neigh] = part_cell[kc][n];
                                nb_neigh++;
                            }
                        }
                    }
                    // Find which particle is the closest to the new position
                    new_x = xc; // centroid of fine grid cell
                    new_z = zc; // centroid of fine grid cell

                    if ( new_x >= model.xmin && new_x <= model.xmax && new_z >= model.zmin && (new_z <= model.zmax || new_z <= h) ) { 
                        // printf("Adding marker\n");
                        if ( count_part_reuse < nb_part_reuse && nb_part_reuse>0  ) {
                            new_ind = part_reuse[count_part_reuse];
                            count_part_reuse++;
                            // printf("reuse - nb_part_reuse = %d - count_part_reuse = %d\n", nb_part_reuse, count_part_reuse);
                        }
                        else {
                            // printf("create\n");
                            new_ind = particles->Nb_part;
                            if (particles->Nb_part+1>=particles->Nb_part_max) {
                                printf("Max number of particles reached, Exiting...");
                                exit(190);
                            }
                            particles->Nb_part++;
                            count_created++;
                        }

                        // Add 1 particule
                        particles->x[new_ind]          = new_x;
                        particles->z[new_ind]          = new_z;
                        particles->generation[new_ind] = 1;

                        double distance, closest_distance = model.xmax - model.xmin;
                        int ind, closest_neighbour = neighbours[0];

                        // Loop on the particles inside the current cell, find which one is the closest to the newly created one.
                        for (ind=0; ind<nb_neigh; ind++) {

                            distance = sqrt( pow (particles->x[new_ind] - particles->x[neighbours[ind]], 2.0) + pow(particles->z[new_ind] - particles->z[neighbours[ind]], 2.0) );

                            if (distance<closest_distance) {
                                closest_distance  = distance;
                                closest_neighbour = neighbours[ind];
                            }
                        }

                        // Closest neighbour was found, now assign properti of the new marker
                        AssignMarkerProperties ( particles, new_ind, closest_neighbour, &model, mesh, model.direct_neighbour );
                        // printf("new_x = %lf, new_z = %lf, ph = %d\n", new_x, new_z, particles->phase[new_ind]);
                    }
                    DoodzFree(neighbours);
                    nb_new_parts++;
                }
            }
        }
        printf("%d markers reused, %d markers created, out of %d new markers\n", count_part_reuse, count_created, nb_new_parts);

        // // TBD
        // for (int k=0; k<Nb_part; k++) {
        //     if ( particles->generation[k] == 1 )  {
        //         printf("   k = %d\n",  k);
        //     }
        // }

        // Freedom
        for (kc=0; kc<ncx*ncz; kc++) {
            DoodzFree(part_cell[kc]);
        }
        DoodzFree(part_cell);
    }
    // Infos
    MinMaxArrayI(particles->phase,          1, particles->Nb_part, "phase       ");
    MinMaxArrayTagChar(mesh->BCt_fine.type, 1, ncx*ncz           , "Cell_flags  ", mesh->BCt_fine.type );
    MinMaxArrayTagInt(nb_part_cell,         1, ncx*ncz           , "nb_part_cell", mesh->BCt_fine.type );
    DoodzFree(nb_part_cell);
    // DoodzFree(Cell_flags);
    DoodzFree(part_reuse);

    // nb_part_cell = DoodzCalloc( Ncx*Ncz, sizeof(int));

    // // Count particles in cells
    // for (int k=0; k<Nb_part; k++) {
    //     if ( particles->phase[k] != -1 )  {
    //         // Get the column:
    //         distance = ( particles->x[k] - (xmin + model.dx/2) );
    //         ic       = ceil( (distance/model.dx) + 0.5) - 1;
    //         if (ic<0   ) ic = 0;
    //         if (ic>=Ncx) ic = Ncx-1;

    //         // Get the line:
    //         distance = ( particles->z[k] - (zmin + model.dz/2) );
    //         jc       = ceil( (distance/model.dz) + 0.5) - 1;
    //         if (jc<0   ) jc = 0;
    //         if (jc>=Ncz) jc = Ncz-1;

    //         // Center
    //         kc = ic + jc*Ncx;

    //         // Add particle to finer mesh cell
    //         nb_part_cell[kc] += 1;

    //     }
    // }

    // for (ic=0; ic<Ncx; ic++) {
    //     for (jc=0; jc<Ncz; jc++) {
    //         kc = ic + jc*Ncx;
    //         if ( mesh->BCt.type[kc]!=30 && nb_part_cell[kc]==0) {
    //             printf("There is a hole @ x=%lf z=%lf!\n", mesh->xc_coord[ic]*scaling.L, mesh->zc_coord[jc]*scaling.L); //exit(1);
    //         }
    //     }
    // }

    // DoodzFree(nb_part_cell);

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CountPartCell_OLD( markers* particles, grid *mesh, params model, surface topo, surface topo_ini, int reseed_markers, scale scaling ) {

    // This function counts the number of particle that are currently in each cell of the domain.
    // The function detects cells that are lacking of particle and call the particle re-seeding routine.

    int k, ith, ic, jc, ip, kc, kd, l, nb, start;
    int nthreads=1, Nx=mesh->Nx, Nz=mesh->Nz;
    int Nb_part=particles->Nb_part;
    int Ncx=Nx-1, Ncz=Nz-1, nxl;
    double dx = model.dx, dz = model.dz, distance, weight, dxm, dzm;
    int finite_strain = model.finite_strain;
    int track_T_P_x_z   = model.track_T_P_x_z;

    int sed_phase=model.surf_ised1;
    int time_My = floor(model.time*scaling.t / (3600*365.25*24*1e6));
    if ( time_My % 2 > 0 ) sed_phase = model.surf_ised1;
    else                   sed_phase = model.surf_ised2;

    printf("USING NEW PART IN CELL\n");
    if (reseed_markers==1) printf("OLD NUMBER OF MARKERS = %02d\n", particles->Nb_part);
    
    // First operation - compute phase proportions on centroids and vertices
    int cent=1, vert=0, prop=1, interp=0;
    P2Mastah ( &model, *particles, NULL, mesh, NULL, mesh->BCp.type,  0, 0, prop, cent, model.interp_stencil);
    P2Mastah ( &model, *particles, NULL, mesh, NULL, mesh->BCg.type,  0, 0, prop, vert, model.interp_stencil);
    
    // Split the domain in N threads in x direction
// #pragma omp parallel
//     {
//         nthreads = omp_get_num_threads();
//     }

    nthreads = 1;

    int nx[nthreads], ncx[nthreads], nx_e[nthreads], ncx_e[nthreads], nz_e = Nz+2, ncz_e = Ncz+2;
    int num_new[nthreads], id_new[nthreads], npart[nthreads], npreuse[nthreads], nnewp[nthreads], npartr[nthreads], nb_new[nthreads];
    double **xg, **xc, **xg_e, **xc_e, *zg_e, *zc_e, ***phv, ***phc;
    double **npvolc, **npvolv;
    int **pidx, **npc, **mpc, **npv, **ipreuse, ***ipcell, **npcell, **ind_list, **pidxr;
    int neighs=0, oo=0;
    double **newx, **newz;
    int    **newi;
    int offc=0;
    int flag1;
    
    // Get the value of nx for each processor
    for (ith=0; ith<nthreads; ith++) {
        ncx[ith]  = (Ncx-Ncx%nthreads)/nthreads;
        if (nthreads>1 && ith==nthreads-1) {
            ncx[ith] += (Ncx%nthreads); // Add additional nodes
        }
        ncx_e[ith] = ncx[ith] + 2;
        nx[ith]    = ncx[ith] + 1;
        nx_e[ith]  = nx[ith]  + 2;
    }

    // Create extended zg array
    zg_e = DoodzCalloc( nz_e,  sizeof(double*));
    zc_e = DoodzCalloc( ncz_e, sizeof(double*));
    zg_e[0] = mesh->zg_coord[0] - dz;
    for (k=1; k<nz_e; k++) {
        zg_e[k] = zg_e[k-1] + dz;
    }

    // Create extended zc array
    for (k=0;   k<ncz_e; k++) zc_e[k]   = 0.5*( zg_e[k] + zg_e[k+1]  );

    // Allocate parallel grid
    pidx = DoodzCalloc( nthreads, sizeof(int*));
    pidxr= DoodzCalloc( nthreads, sizeof(int*));
    npc  = DoodzCalloc( nthreads, sizeof(int*));
    npvolv  = DoodzCalloc( nthreads, sizeof(double*));
    npvolc  = DoodzCalloc( nthreads, sizeof(double*));
    npv  = DoodzCalloc( nthreads, sizeof(int*));
    xg   = DoodzCalloc( nthreads, sizeof(double*));
    xc   = DoodzCalloc( nthreads, sizeof(double*));
    xg_e = DoodzCalloc( nthreads, sizeof(double*));
    xc_e = DoodzCalloc( nthreads, sizeof(double*));
    phv  = DoodzCalloc( nthreads, sizeof(double**));
    phc  = DoodzCalloc( nthreads, sizeof(double**));
    mpc  = DoodzCalloc( nthreads, sizeof(int*));

    ipreuse = DoodzCalloc( nthreads, sizeof(int*));
    ipcell  = DoodzCalloc( nthreads, sizeof(int**));
    npcell  = DoodzCalloc( nthreads, sizeof(int*));
    //    flag    = DoodzCalloc( nthreads, sizeof(char*));
    ind_list= DoodzCalloc( nthreads, sizeof(int*));
    newi    = DoodzCalloc( nthreads, sizeof(int*));
    newx    = DoodzCalloc( nthreads, sizeof(double*));
    newz    = DoodzCalloc( nthreads, sizeof(double*));

    for (ith=0; ith<nthreads; ith++) {
        xg[ith]     = DoodzCalloc(        nx[ith], sizeof(double));
        xc[ith]     = DoodzCalloc(       ncx[ith], sizeof(double));
        xg_e[ith]   = DoodzCalloc(      nx_e[ith], sizeof(double));
        xc_e[ith]   = DoodzCalloc(     ncx_e[ith], sizeof(double));

        npv[ith]    = DoodzCalloc(     nx[ith]*Nz, sizeof(int));
        npc[ith]    = DoodzCalloc( (ncx[ith])*Ncz, sizeof(int));
        npvolv[ith]    = DoodzCalloc(     nx[ith]*Nz, sizeof(double));
        npvolc[ith]    = DoodzCalloc( (ncx[ith])*Ncz, sizeof(double));
        phv[ith]    = DoodzCalloc( model.Nb_phases, sizeof(double*));
        phc[ith]    = DoodzCalloc( model.Nb_phases, sizeof(double*));

        mpc[ith]    = DoodzCalloc(4*(ncx_e[ith+0])*ncz_e, sizeof(int));
        npcell[ith] = DoodzCalloc(4*(ncx_e[ith]+0)*ncz_e, sizeof(int));
        ipcell[ith] = DoodzCalloc(4*(ncx_e[ith]+0)*ncz_e, sizeof(int*));
        //        flag[ith]   = DoodzCalloc(4*(ncx[ith])*Ncz, sizeof(char));
        for (k=0; k<model.Nb_phases; k++) {
            phv[ith][k]  = DoodzCalloc(     nx[ith]*Nz, sizeof(double));
            phc[ith][k]  = DoodzCalloc( (ncx[ith])*Ncz, sizeof(double));
        }
    }

    // start parallel section
//#pragma omp parallel private (ith, k, ip, distance, weight, dxm, dzm, ic, jc, kc, kd, l, neighs, oo, nb, offc, start, flag1, nxl) shared (ncx_e, ncz_e, nx_e, nz_e, xg, xg_e, xc, xc_e, dx, dz, nx, ncx, Nb_part, npart, particles, nthreads, pidx, phv, phc, npv, npc, npvolv, npvolc, npreuse, mesh, ind_list, nnewp, newi, newx, newz, pidxr, npartr, nb_new, id_new, num_new, topo, topo_ini) firstprivate(Nx, Nz, Ncx, Ncz, model, sed_phase )
             for (ith=0; ith<nthreads; ith++)
    {
        // ith = omp_get_thread_num();
        
        // Compute anchor for global cell index in x direction
        int igc, ii, ancrage_cell = 0;
        if (ith == 0) ancrage_cell = 0;
        if (ith  > 0) {
            for (ii=0;ii<ith;ii++) ancrage_cell += (nx[ii]-1);
        }

        nnewp[ith]   = 0;
        npart[ith]   = 0;
        npartr[ith]  = 0;
        npreuse[ith] = 0;
        num_new[ith] = 0;
        id_new[ith]  = 0;
        nb_new[ith]  = 0;

        // Create parallel grid coordinates: vertices
        offc=0;
        for (ic=1; ic<ith+1; ic++) offc += (nx[ic-1]-1);

        xg[ith][0] = mesh->xg_coord[offc];
        for (k=1; k<nx[ith]; k++) {
            xg[ith][k] = xg[ith][k-1] + dx;
        }

        xg_e[ith][0] = xg[ith][0] - dx;
        for (k=1; k<nx_e[ith]; k++) {
            xg_e[ith][k] = xg_e[ith][k-1] + dx;
        }

        // Create parallel grid coordinates: centroids
        for (k=0;   k<ncx[ith]; k++) xc[ith][k]   = 0.5*( xg[ith][k] + xg[ith][k+1]  );
        for (k=0; k<ncx_e[ith]; k++) xc_e[ith][k] = 0.5*( xg_e[ith][k] + xg_e[ith][k+1]  );

        // Count # of particles per processor
        if (nthreads>1) {
            for (k=0; k<Nb_part; k++) {
                if ( particles->phase[k] != -1 )  {
                if (particles->x[k] > xg_e[ith][0] && particles->x[k] < xg_e[ith][nx_e[ith]-1]) npart[ith]++;
                if (particles->x[k] > xg  [ith][0] && particles->x[k] < xg  [ith][nx  [ith]-1]) npartr[ith]++;
            }
            }
        }
        else {
            npart[ith]  = particles->Nb_part;
            npartr[ith] = particles->Nb_part;
        }

//        printf("ith = %d --> Npart  = %d tot = %d\n", ith,  npart[ith],  particles->Nb_part);
//        printf("ith = %d --> Npartr = %d\n", ith,  npartr[ith]);

        // Create list of particle indices per processor
        pidx[ith]    = DoodzCalloc( 2*npart[ith]  , sizeof(int)); // over-allocate
        pidxr[ith]   = DoodzCalloc( 2*npartr[ith]  , sizeof(int)); // over-allocate
        npart[ith]   = 0;
        npartr[ith]  = 0;
        npreuse[ith] = 0;

        if (nthreads>1) {
            // Multi core: particle distribution over the threads
            for (k=0; k<Nb_part; k++) {

                if ( particles->phase[k] != -1 )  {

                    if ( particles->x[k] > xg_e[ith][0] && particles->x[k] < xg_e[ith][nx_e[ith]-1]) {
                        pidx[ith][npart[ith]] = k;
                        npart[ith]++;
                    };

                    if ( particles->x[k] > xg  [ith][0] && particles->x[k] < xg  [ith][nx  [ith]-1]) {
                        pidxr[ith][npartr[ith]] = k;
                        npartr[ith]++;
                    }
                }
            }
        }
        else {
            // Single core: no particle distribution over the threads
            npart[ith] = particles->Nb_part;
            npartr[ith] = particles->Nb_part;
            for (ip=0; ip<npart[ith]; ip++)  {
                pidx[ith][ip]  = ip;
                pidxr[ith][ip] = ip;
            }
        }

//        printf("ith = %d --> Npart  = %d tot = %d\n", ith,  npart[ith],  particles->Nb_part);
//        printf("ith = %d --> Npartr = %d\n", ith,  npartr[ith]);

//        // Count number of particles per nodes and cell
//        for (ip=0; ip<npart[ith]; ip++) {
//
//            k = pidx[ith][ip];
//
//            if ( particles->phase[k] != -1 && particles->x[k] > xg[ith][0] - dx/2.0 && particles->x[k] < xg[ith][nx[ith]-1] + dx/2.0) {
//
//                // -------- Check NODES ---------- //
//                distance = (particles->x[k]-          xg[ith][0]);
//                ic       = ceil((distance/dx)+0.5) -1;
//                distance = (particles->z[k]-mesh->zg_coord[0]);
//                jc       = ceil((distance/dz)+0.5) -1;
//
//                //
////                printf("particles->x[k] = %2.2e ; (PT CENTRE X) = %2.2e; gauche  = %2.2e; xg droite = %2.2e\n", particles->x[k],xg[ith][ic],xg[ith][ic]-dx/2.0, xg[ith][ic]+dx/2.0 );
////                printf("particles->z[k] = %2.2e ; (PT CENTRE Z) = %2.2e; dessous = %2.2e; dessus    = %2.2e\n", particles->z[k],mesh->zg_coord[jc],mesh->zg_coord[jc]-dz/2.0, mesh->zg_coord[jc]+dz/2.0 );
//
//                weight = 1.0;
//
//                if (model.interp_stencil==1) {
//                    dxm = 1.0*fabs( xg[ith][ic] - particles->x[k]);
//                    dzm = 1.0*fabs( mesh->zg_coord[jc] - particles->z[k]);
//                    weight = (1.0-dxm/(dx/2.0))*(1.0-dzm/(dz/2.0));
//                }
//
////                if (dxm>dx/2.0) {
////                 printf("nodes Outside X !!!");
////                    exit(1);
////                }
////                if (dzm>dz/2.0) {
////                 printf("nodes Outside Z !!!");
////                    exit(1);
////                }
//
//                // If particles are around the nodes: count them
//                if ( particles -> phase[k]>=0 ) {
//                    npv[ith][ic + jc * (nx[ith])]++;
//                    npvolv[ith][ic + jc * (nx[ith])]+=weight;
//                    phv[ith][particles->phase[k]][ic + jc * (nx[ith])] +=weight;
////                    phv[ith][particles->phase[k]][ic + jc * (nx[ith])] +=1.0;
//
//                }
//
//                 //-------- Check CELLS ---------- //
//
//
//                // Filter out overlapping domain
//                if (particles->x[k] > xg[ith][0] && particles->x[k] < xg[ith][nx[ith]-1]) {
//                    distance = (particles->x[k]-          xc[ith][0]);
//                    ic       = ceil((distance/dx)+0.5) -1;
//                    if (ic<0)    ic = 0;
//                    if (ic>nx[ith]-2) ic = nx[ith]-2;
//                    distance = (particles->z[k]-mesh->zc_coord[0]);
//                    jc       = ceil((distance/dz)+0.5) -1;
//                    if (jc<0)    jc = 0;
//                    if (jc>Nz-2) jc = Nz-2;
//
//                    weight = 1.0;
//
//                    if (model.interp_stencil==1) {
//                        dxm    = 1.0*fabs( xc[ith][ic] - particles->x[k]);
//                        dzm    = 1.0*fabs( mesh->zc_coord[jc] - particles->z[k]);
//                        weight = (1.0-dxm/(dx/2.0))*(1.0-dzm/(dz/2.0));
//                    }
//                    // If particles are around the nodes: count them
//                    if ( particles -> phase[k]>=0 )  { //&& (model.interp_stencil==0 || model.interp_stencil==1
//                        npc[ith][ic + jc * (nx[ith]-1)]++;
//                        npvolc[ith][ic + jc * (nx[ith]-1)]+=weight;
//                        phc[ith][particles->phase[k]][ic + jc * (nx[ith]-1)]+=weight;
////                        phc[ith][particles->phase[k]][ic + jc * (nx[ith]-1)]+=1.0;
//                    }
//
////                    if ( particles -> phase[k]>=0 && model.interp_stencil==4) {
////
////                        // Global cell index in x direction
////                        igc = ancrage_cell + ic;
////
////                        // if particle lies in the SW quadrant of the cell -- reference to NE
////                        if ( particles->x[k] < xc[ith][ic] &&  particles->z[k] < mesh->zc_coord[jc] ) {
////                                                   Compute4Cell( ic  , jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NE
////                            if (igc>0            ) Compute4Cell( ic-1, jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NW
////                            if (igc>0   &&   jc>0) Compute4Cell( ic-1, jc-1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SW
////                            if (             jc>0) Compute4Cell( ic  , jc-1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SE
////                        }
////
////                        // if particle lies in the SE quadrant of the cell -- reference to NW
////                        if ( particles->x[k] > xc[ith][ic] &&  particles->z[k] < mesh->zc_coord[jc] ) {
////                                                   Compute4Cell( ic  , jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NW
////                            if (igc<Ncx          ) Compute4Cell( ic+1, jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NE
////                            if (igc<Ncx  &&  jc>0) Compute4Cell( ic+1, jc-1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SE
////                            if (             jc>0) Compute4Cell( ic  , jc-1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SW
////                        }
////                        // if particle lies in the NW quadrant of the cell -- reference to SE
////                        if ( particles->x[k] < xc[ith][ic] &&  particles->z[k] > mesh->zc_coord[jc] ) {
////                                                   Compute4Cell( ic  , jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SE
////                            if (igc>0            ) Compute4Cell( ic-1, jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SW
////                            if (           jc<Ncz) Compute4Cell( ic  , jc+1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NE
////                            if (igc>0  &&  jc<Ncz) Compute4Cell( ic-1, jc+1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NW
////
////                        }
////                        // if particle lies in the NE quadrant of the cell -- reference to SW
////                        if ( particles->x[k] > xc[ith][ic] &&  particles->z[k] > mesh->zc_coord[jc] ) {
////                                                   Compute4Cell( ic  , jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SW
////                            if (igc<Ncx          ) Compute4Cell( ic+1, jc  , nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> SE
////                            if (           jc<Ncz) Compute4Cell( ic  , jc+1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NW
////                            if (igc<Ncx && jc<Ncz) Compute4Cell( ic+1, jc+1, nx[ith]-1, xc, mesh, particles, k, npc,  npvolc, phc, ith );  // ---> NE
////                        }
////
////
////                    }
//
//                    if (dxm>dx/2.0) {
//                        printf("centers Outside X !!!");
//                        exit(1);
//                    }
//                    if (dzm>dz/2.0) {
//                        printf("centers Outside Z !!!");
//                        exit(1);
//                    }
//
//
//
//                }
//            }
//        }

        // Count number of particles available for re-use (no overlapping particles)
        for (ip=0; ip<npartr[ith]; ip++) {
            k = pidxr[ith][ip];
            if ( particles->phase[k] == -1 ) {
                // -------- Check What's OUT  ----------
                // If particles are taged as 'out' (-1): count them for re-use
                (npreuse[ith])++;
            }
        }

//        // Phase proportions
//        for (k=0; k<model.Nb_phases; k++) {
//
//            // Calculate phase proportions vertices
//            for (ip=0; ip<nx[ith]*Nz; ip++) {
//                //if (npv[ith][ip]>0) phv[ith][k][ip] /= npv[ith][ip];
//                if (npv[ith][ip]>0) phv[ith][k][ip] /= npvolv[ith][ip];
//            }
//
//            // Calculate phase proportions centroids
//            for (ip=0; ip<ncx[ith]*Ncz; ip++) {
//                //if (npc[ith][ip]>0) phc[ith][k][ip] /= npc[ith][ip];
//                if (npc[ith][ip]>0) phc[ith][k][ip] /= npvolc[ith][ip];
//            }
//
//            // Fill global arrays - centroids
//            for (ic=0; ic<ncx[ith]; ic++) {
//                for (jc=0; jc<Ncz; jc++) {
//                    ip =        ic + jc*ncx[ith];
//                    kc = offc + ic + jc*Ncx;
//                    mesh->nb_part_cell[kc]    = npc[ith][ip];
//                    if ( mesh->nb_part_cell[kc]>0 ) mesh->phase_perc_n[k][kc] = phc[ith][k][ip];
//                }
//            }
//
//            start=0;
//            if (ith>0) start = 1; // Do not include overlapping nodal values
//
//            // Fill global arrays - vertices
//            for (ic=start; ic<(nx[ith]); ic++) {
//                for (jc=0; jc<Nz; jc++) {
//                    ip =        ic + jc*nx[ith];
//                    kc = offc + ic + jc*Nx;
//                    mesh->nb_part_vert[kc]    = npv[ith][ip];
//                    if ( mesh->nb_part_vert[kc]>0 ) mesh->phase_perc_s[k][kc] = phv[ith][k][ip];
//                }
//            }
//
//        }
        
//        MinMaxArrayTag( mesh->phase_perc_n[0], 1.0, Ncx*Ncz, "phase % n", mesh->BCp.type );
        


        // ------------------------------------ RE-SEEDING ------------------------------------ //
        // Re-seeding on a double resolution grid
        if (reseed_markers==1) {

            newi[ith]   = DoodzCalloc( npart[ith]  , sizeof(int));
            newx[ith]   = DoodzCalloc( npart[ith]  , sizeof(double));
            newz[ith]   = DoodzCalloc( npart[ith]  , sizeof(double));

            // Count number of particle per re-seeding volume
            for (ip=0; ip<npart[ith]; ip++) {

                k = pidx[ith][ip];

                if ( particles->phase[k] != -1 ) {

                    // -------- Check CELLS ---------- //
                    distance                      = particles->x[k] - (xc_e[ith][0]-dx/4.0);
                    ic                            = ceil((distance/dx*2.0)+0.5) -1;
                    if (ic<0)                  ic = 0;
                    if (ic>2*(ncx_e[ith]+0)-1) ic = 2*(ncx_e[ith]+0)-1;
                    distance                      = particles->z[k] - (zc_e[0]-dz/4.0);
                    jc                            = ceil((distance/dz*2.0)+0.5) -1;
                    if (jc<0)                  jc = 0;
                    if (jc>2*ncz_e-1)          jc = 2*ncz_e-1;

                    // If particles are around the nodes: count them
                    if ( particles -> phase[k]>=0 ) {
                        kc = ic + jc * 2*(ncx_e[ith]+0);
                        mpc[ith][kc]++;

                        // Periodic
                        if (model.periodic_x == 1 ) {

                            // Count particles from East side to West side of first thread
                            if (ic==2*(ncx_e[ith]+0)-1 && ith==0) {
                                kc =  0 + jc * 2*(ncx_e[ith]+0);
                                mpc[ith][kc]++;
                            }
                            if (ic==2*(ncx_e[ith]+0)-2 && ith==0) {
                                kc =  1 + jc * 2*(ncx_e[ith]+0);
                                mpc[ith][kc]++;
                            }

                            // Count particles from West side to East side of last thread
                            if (ic==0 && ith==nthreads-1) {
                                kc =  2*(ncx_e[ith]+0)-1+ jc * 2*(ncx_e[ith]+0);
                                mpc[ith][kc]++;
                            }
                            if (ic==1 && ith==nthreads-1) {
                                kc =  2*(ncx_e[ith]+0)-2+ jc * 2*(ncx_e[ith]+0);
                                mpc[ith][kc]++;
                            }

                        }
                    }
                }
            }

            // Build list of particles to reuse
            ipreuse[ith] = DoodzCalloc( npreuse[ith], sizeof(int));
            ic           = 0;
            for (ip=0; ip<npartr[ith]; ip++) {
                k = pidxr[ith][ip];
                if ( particles->phase[k] == -1 ) {
                    ipreuse[ith][ic] = k;
//                    printf("ic = %d;k = %d\n", ic,k);
                    ic++;
                }
            }
            
            //printf("ic = %d and npreuse[ith] = %d \n", ic, npreuse[ith]); //=> Ic and npreuse looks ok (i.e. equal and pair)
            

            // Build list of particles index and particles number per cell
            for (k=0;k<4*(ncx_e[ith]+0)*ncz_e;k++) ipcell[ith][k] = DoodzCalloc( mpc[ith][k]  , sizeof(int));

            // Loop on particles, give their ID's to the cell that owns them
            for (ip=0; ip<npart[ith]; ip++) {

                k = pidx[ith][ip];

                if ( particles->phase[k] != -1 ) {
                    // -------- Check CELLS ---------- //
                    distance                    = particles->x[k] - (xc_e[ith][0]-dx/4.0);
                    ic                          = ceil((distance/dx*2.0)+0.5) -1;
                    if (ic<0)                ic = 0;
                    if (ic>2*(ncx_e[ith]+0)-1) ic = 2*(ncx_e[ith]+0)-1;
                    distance                    = particles->z[k] - (zc_e[0]-dz/4.0);
                    jc                          = ceil((distance/dz*2.0)+0.5) -1;
                    if (jc<0)                jc = 0;
                    if (jc>2*ncz_e-1)        jc = 2*ncz_e-1;
                    kc = ic + jc * 2*(ncx_e[ith]);

                    ipcell[ith][kc][npcell[ith][kc]] = k;
                    npcell[ith][kc] ++;


                    // Periodic
                    if (model.periodic_x == 1 ) {

                        // Count particles from East side to West side of first thread
                        if (ic==2*(ncx_e[ith]+0)-1 && ith==0) {
                            kc =  0 + jc * 2*(ncx_e[ith]+0);
                            ipcell[ith][kc][npcell[ith][kc]] = k;
                            npcell[ith][kc] ++;
                        }
                        if (ic==2*(ncx_e[ith]+0)-2 && ith==0) {
                            kc =  1 + jc * 2*(ncx_e[ith]+0);
                            ipcell[ith][kc][npcell[ith][kc]] = k;
                            npcell[ith][kc] ++;
                        }
                        // Count particles from West side to East side of last thread
                        if (ic==0 && ith==nthreads-1) {
                            kc =  2*(ncx_e[ith]+0)-1+ jc * 2*(ncx_e[ith]+0);
                            ipcell[ith][kc][npcell[ith][kc]] = k;
                            npcell[ith][kc] ++;

                        }
                        if (ic==1 && ith==nthreads-1) {
                            kc =  2*(ncx_e[ith]+0)-2+ jc * 2*(ncx_e[ith]+0);
                            ipcell[ith][kc][npcell[ith][kc]] = k;
                            npcell[ith][kc] ++;
                        }
                    }
                }
            }

            // Loop on cells and add particles
            for (k=2; k<2*ncx_e[ith]-2; k++) {
                for (l=2; l<2*ncz_e-2; l++) {

                    kd = k + l * 2*(ncx_e[ith]+0);

                    // ---  Get index of the corresponding coarse grid cell

                    // x index
                    distance         = (xg_e[ith][0]+dx/4.0+(k)*dx/2.0) - mesh->xc_coord[0];
                    ic               = ceil((distance/dx)+0.5) -1;
                    if (ic<0)     ic = 0;
                    if (ic>Ncx-1) ic = Ncx-1;

                    // z index
                    distance         = (zg_e[0]+dz/4.0+(l)*dz/2.0) - mesh->zc_coord[0];
                    jc               = ceil((distance/dz)+0.5) -1;
                    if (jc<0)     jc = 0;
                    if (jc>Ncz-1) jc = Ncz-1;


                    ip = ic + jc * Ncx;
                    flag1 = mesh->BCt.type[ip]; //mesh->BCp.type[0][ip];

                    // ---  Get index of the corresponding coarse grid cell

                    if (flag1 != 30) {

//                        printf("number = %d condition = %d\n",  mpc[ith][kd], (int)(ceil(particles->min_part_cell) / 4.0) );
                        if ( mpc[ith][kd] < (int)(ceil(particles->min_part_cell) / 4.0) ) {

                            neighs = 0;
                            oo     = 0;

                            // list the number of particles in the neighbouring cells
                            neighs += mpc[ith][kd];
                            nxl     = 2*(ncx_e[ith]+0);
                            // Direct neighbours
                            neighs += mpc[ith][kd-1];
                            neighs += mpc[ith][kd+1];
                            neighs += mpc[ith][kd-nxl];
                            neighs += mpc[ith][kd+nxl];
                            // Diagonal neighbours
                            neighs += mpc[ith][kd-1-nxl];
                            neighs += mpc[ith][kd+1-nxl];
                            neighs += mpc[ith][kd-1+nxl];
                            neighs += mpc[ith][kd+1+nxl];
                            // Extended neighbours
                            neighs += mpc[ith][kd-2];
                            neighs += mpc[ith][kd+2];
                            neighs += mpc[ith][kd-2*nxl];
                            neighs += mpc[ith][kd+2*nxl];
                            // Extended diagonal neighbours
                            neighs += mpc[ith][kd-2 - 2*nxl];
                            neighs += mpc[ith][kd+2 - 2*nxl];
                            neighs += mpc[ith][kd-2 + 2*nxl];
                            neighs += mpc[ith][kd+2 + 2*nxl];
                            // Missing guys
                            neighs += mpc[ith][kd-1 -2*nxl];
                            neighs += mpc[ith][kd+1 -2*nxl];
                            neighs += mpc[ith][kd-1 +2*nxl];
                            neighs += mpc[ith][kd+1 +2*nxl];
                            neighs += mpc[ith][kd-2 -1*nxl];
                            neighs += mpc[ith][kd-2 +1*nxl];
                            neighs += mpc[ith][kd+2 -1*nxl];
                            neighs += mpc[ith][kd+2 +1*nxl];

                            if (neighs == 0) {
                                printf("All the neighbouring CELLS of ix = %d iz = %d are empty, simulation will stop\n", k, l);
                                printf("fine   cell: xc = %2.2e zc = %2.2e flag = %d thread #%d nthreads=%d\n",(xg[ith][0]+dx/4.0+(k)*dx/2.0)*scaling.L, (mesh->zg_coord[0]+dz/4.0+(l)*dz/2.0)*scaling.L, flag1, ith, nthreads);
                                printf("coarse cell: xc = %2.2e zc = %2.2e flag = %d\n", mesh->xc_coord[ic]*scaling.L, mesh->zc_coord[jc]*scaling.L, mesh->BCt.type[ip]);
                                exit(333);
                            }

                            ind_list[ith] = DoodzMalloc( neighs * sizeof(int));

                            // add particles from current cell
                            for ( nb=0; nb<mpc[ith][kd]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd][nb];
                                oo++;
                            }

                            // Add particles from neigbours
                            for ( nb=0; nb<mpc[ith][kd-1]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-1][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+1]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+1][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd-nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-nxl][nb];
                                oo++;
                            }
                            for ( nb=0; nb<mpc[ith][kd+nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+nxl][nb];
                                oo++;
                            }

                            // Diagonal cells: add particles from neigbours
                            for ( nb=0; nb<mpc[ith][kd-1-nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-1-nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+1-nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+1-nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd-1+nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-1+nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+1+nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+1+nxl][nb];
                                oo++;
                            }

                            // Extended neighbours
                            for ( nb=0; nb<mpc[ith][kd-2]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-2][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+2]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+2][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd-2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+2*nxl][nb];
                                oo++;
                            }

                            // Extended diagonal neighbours
                            for ( nb=0; nb<mpc[ith][kd-2-2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-2-2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+2-2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+2-2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd-2+2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-2+2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+2+2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+2+2*nxl][nb];
                                oo++;
                            }

                            // Missing guys
                            for ( nb=0; nb<mpc[ith][kd-1 -2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-1 -2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+1 -2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+1 -2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd-1 +2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-1 +2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+1 +2*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+1 +2*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd-2 -1*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-2 -1*nxl][nb];
                                oo++;
                            }


                            for ( nb=0; nb<mpc[ith][kd-2 +1*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd-2 +1*nxl][nb];
                                oo++;
                            }


                            for ( nb=0; nb<mpc[ith][kd+2 -1*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+2 -1*nxl][nb];
                                oo++;
                            }

                            for ( nb=0; nb<mpc[ith][kd+2 +1*nxl]; nb++ ) {
                                ind_list[ith][oo] = ipcell[ith][kd+2 +1*nxl][nb];
                                oo++;
                            }

                            // Identify closer points
                            AddPartCell2( pidx[ith], &(npart[ith]), particles, *mesh, k-2, l-2, ind_list[ith], model, neighs, topo, xg[ith], mesh->zg_coord, ic, &(nnewp[ith]), newx[ith], newz[ith], newi[ith], sed_phase, topo_ini );
                            DoodzFree(ind_list[ith]);

                        }
                    }
                }
            }
        }
    }

    // Generate particle indices
    if (reseed_markers==1) {
        int dummy1 = 0;
        int dummy2 = 0;
        
        // Create indices in global particle arrays for newly created particles (if not reused indices!)
        for (ith=0; ith<nthreads; ith++) {
            num_new[ith] = 0;     // nb of particles using new indices
            id_new[ith]  = 0;     // new indices
            if (nnewp[ith] > npreuse[ith]) {
                num_new[ith] = nnewp[ith] - npreuse[ith]; // nb of particles using new indices
            }
        }
        for (ith=0; ith<nthreads; ith++) {
            if (nnewp[ith] > npreuse[ith]) {
                ip           = 0;
                for (ic=0; ic<ith; ic++) {  // removed the +1 => for (ic=1; ic<ith+1; ic++)
                    ip += num_new[ic];      // removed the -1 => ip += num_new[ic-1]
                }
                id_new[ith]  = Nb_part + ip;
            }
        }
        
        
//        printf("=============================================================\n");
//        ith = 0;
//        printf("num_new[ith] = %d;  nnewp[ith] = %d; dummy1 = %d; dummy2 = %d --- ith = %d --- nth = %d\n", num_new[ith], nnewp[ith], dummy1, dummy2, ith, nthreads);
//        printf("=============================================================\n");
        
        
        // Generate new particles in global array
    // #pragma omp parallel private ( ip, k, ith ) shared( particles, nb_new, newi, newx, newz, id_new, ipreuse, npreuse, nnewp, mesh ) firstprivate( model )
               for (ith=0; ith<nthreads; ith++)
        {

            // ith         = omp_get_thread_num();
            nb_new[ith] = 0;

            // Loop on all new particles
            for (ip=0; ip<nnewp[ith]; ip++) {
                
                if ( ip<npreuse[ith]) {
                    // Reuse indices
                    // printf("Recycle marker %d\n", ipreuse[ith][ip]);
                    k = ipreuse[ith][ip];
                }
                else {
                    // New indices
                    k = id_new[ith] + nb_new[ith];
                    nb_new[ith]++;
                }
                // Security check
                if (k>particles->Nb_part_max) {
                    printf ("Maximum number of particles exceeded!\n Exiting!\n");
                    exit(1);
                }
                else {
                    // Assign coordinates and properties
                    particles->x[k]          = newx[ith][ip];
                    particles->z[k]          = newz[ith][ip];
                    particles->generation[k] = 1;
                    AssignMarkerProperties( particles, k, newi[ith][ip], &model, mesh, model.direct_neighbour );
                }
            }
        }

    }

    for (ith=0; ith<nthreads; ith++) {
        particles->Nb_part += nb_new[ith];
    }
    
//    //________________________________________
//    for (k=0;k<particles->Nb_part;k++) {
//    isoutPart( particles, &model, k );
//    }
//    int dum1 = 0;
//
//    for (k=0;k<particles->Nb_part;k++) {
//        if ( particles->phase[k] == -1 ) dum1++;
//    }
//    printf("dum1 = %d\n", dum1);
//    //________________________________________
    
    if (reseed_markers==1) printf("NEW NUMBER OF MARKERS = %02d\n", particles->Nb_part);

    IsNanArray2DFP(mesh->phase_perc_s[0], Nx*Nz);
    IsInfArray2DFP(mesh->phase_perc_s[0], Nx*Nz);
    IsNanArray2DFP(mesh->phase_perc_n[0], Ncx*Ncz);
    IsInfArray2DFP(mesh->phase_perc_n[0], Ncx*Ncz);

    // Periodic - only for grid values
    double av;
    int c1;
    if (model.periodic_x==1 && (model.Nx-Nx==0) ) {
        for (ip=0;ip<model.Nb_phases;ip++) {
            for( l=0; l<Nz; l++) {
                c1 = l*Nx + Nx-1;
                av = 0.5*(mesh->phase_perc_s[ip][c1] + mesh->phase_perc_s[ip][l*Nx]);
                mesh->phase_perc_s[ip][c1] = av; mesh->phase_perc_s[ip][l*Nx] = av;
            }
        }
    }


//    // check
//    printf( "CHECK PARTICULE RESEEDING\n");
//    double a;
//    int sum = 0, sum_tot = 0, kn, ln;
//    for (kc=0; kc<Nx*Nz; kc++) {
////        kn = mesh->kn[kc];
////        ln = mesh->ln[kc];
//        a = 0.0;
//        sum_tot++;
//
//        for (k=0; k<model.Nb_phases; k++) {
//            a += mesh->phase_perc_s[k][kc];
//
//
//
//
//                if (k==model.Nb_phases-1 &&  mesh->BCg.type[kc] != 30 && a<1e-13) {
//                    kn = mesh->kn[kc];
//                    ln = mesh->ln[kc];
//                    printf("node = %d %2.2e %0d %0d %0d %0d %0f phase=%0d model.Nb_phases=%0d\n", kc, a, mesh->BCg.type[kc], kc, ln*Nx+kn, npv[0][kc], phv[0][k][kc], k, model.Nb_phases  );
//                    printf("line %0d column %0d xg = %2.2e zg = %2.2e\n", ln,kn, mesh->xg_coord[0][kn]*scaling.L, mesh->zg_coord[0][ln]*scaling.L );
//                    sum++;
//                    printf("\n");
//                }
//
////            if (kn == 135 && ln == 15) printf("a = %2.2e\n");
//
//        }
//
//    }
//    printf("detected %d empty vertices over %d vertices\n", sum, sum_tot);
//
//    ln = 172; kn = 163;
//    printf("topo.height[%0d] = %2.2e %d\n", kn, topo.height[kn]*scaling.L, topo.height[kn]>mesh->zg_coord[0][ln] );
//
    //or (k=0;k<4*(ncx_e[ith]+0)*ncz_e;k++) ipcell[ith][k] = DoodzCalloc( mpc[ith][k]  , sizeof(int));
    // Free parallel grid
    for (ith=0; ith<nthreads; ith++) {
        DoodzFree( xg[ith]);
        DoodzFree( xc[ith]);
        DoodzFree( xg_e[ith]);
        DoodzFree( xc_e[ith]);
        DoodzFree( pidx[ith]);
        DoodzFree( pidxr[ith]);
        DoodzFree( npv[ith]);
        DoodzFree( npc[ith]);
        DoodzFree( npvolv[ith]);
        DoodzFree( npvolc[ith]);
        DoodzFree( mpc[ith]);
        if (reseed_markers==1) DoodzFree( ipreuse[ith]);
        for (k=0; k<model.Nb_phases; k++) {
            DoodzFree( phv[ith][k] );
            DoodzFree( phc[ith][k] );
        }
        if (reseed_markers==1) {
            for (k=0;k<4*(ncx_e[ith]+0)*ncz_e;k++) { // This was a leak --- for (k=0; k<4*(ncx[ith]+1)*Ncz; k++) {
                DoodzFree( ipcell[ith][k] );
            }
        }
        DoodzFree( ipcell[ith] );
        DoodzFree( npcell[ith] );
        DoodzFree( phv[ith] );
        DoodzFree( phc[ith] );
        //        DoodzFree( flag[ith]);

        if (reseed_markers==1) DoodzFree( newi[ith]);
        if (reseed_markers==1) DoodzFree( newx[ith]);
        if (reseed_markers==1) DoodzFree( newz[ith]);
        //        DoodzFree( ind_list[ith]);
    }

    DoodzFree( xg );
    DoodzFree( xc );
    DoodzFree( xg_e );
    DoodzFree( xc_e );
    DoodzFree( pidx );
    DoodzFree( pidxr );
    DoodzFree( npv );
    DoodzFree( npc );
    DoodzFree( npvolv );
    DoodzFree( npvolc );
    DoodzFree( mpc );
    DoodzFree( phv );
    DoodzFree( phc );
    DoodzFree( ipreuse );
    DoodzFree( ipcell );
    DoodzFree( npcell );
    //    DoodzFree( flag );
    DoodzFree( ind_list );

    DoodzFree( newi );
    DoodzFree( newx );
    DoodzFree( newz );

    DoodzFree( zg_e );
    DoodzFree( zc_e );

}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void TransmutateMarkers(markers *particles, mat_prop *materials, double scaling_T) {
  int transmutated_particles = 0;
#pragma omp parallel for reduction(+:transmutated_particles) shared(particles)
  for (int k=0; k<particles->Nb_part-1; k++) {
    const int phase = particles->phase[k];
    const double temperature = particles->T[k];
    const double transmutation = materials->transmutation[phase];

    if ((transmutation == -1 && materials->transmutation_temperature[phase] >= temperature)
        || (transmutation == 1 && materials->transmutation_temperature[phase] <= temperature)) {
      particles->phase[k] = materials->transmutation_phase[phase];
      particles->dual[k] = materials->transmutation_phase[phase];
      transmutated_particles++;
    }
  }

  printf("Transmutated particles = %03d\n", transmutated_particles);
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/