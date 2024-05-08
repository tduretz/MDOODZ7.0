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
//---- M-Doodz header file
#include "mdoodz-private.h"

#ifdef _VG_
#define printf(...) printf("")
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AdvectFreeSurf( markers *topo_chain, params model, scale scaling ) {
    int k;
#pragma omp parallel for shared ( topo_chain ) firstprivate ( model ) private ( k )
    for (k=0; k<topo_chain->Nb_part; k++) {
        topo_chain->x[k] += model.dt*topo_chain->Vx[k];
        topo_chain->z[k] += model.dt*topo_chain->Vz[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SetTopoChainHorizontalCoords( surface *topo, markers *topo_chain, params model, grid Mmesh, scale scaling ) {

    int Nb_part_per_cell = 23;
    int k, counter;
    double dxm = model.dx/Nb_part_per_cell;
    topo_chain->Nb_part = (model.Nx-1)*Nb_part_per_cell;

#pragma omp parallel for shared ( topo_chain ) firstprivate ( model, dxm, scaling ) private ( k )
    for ( k=0; k<topo_chain->Nb_part; k++ ) {

        topo_chain->x[k]     = model.xmin + dxm/2 + k*dxm;
        topo_chain->z[k]     = 0.0/scaling.L;
        topo_chain->z0[k]    = 0.0/scaling.L;
        topo_chain->phase[k] = 0;
    }
    printf( "Topographic chain initialised with %d markers\n", topo_chain->Nb_part );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CorrectTopoIni( markers *particles, mat_prop materials, markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh) {
    
    int    k, Ncx=model.Nx-1;
    double distance, dx=model.dx;
    int in;
    double grid_topo;
    
    for (k=0;k<topo_chain->Nb_part;k++) {
        // Index of the coarse grid column
        distance        = (topo_chain->x[k]-model.xmin-dx/2.0);
        in              = ceil((distance/dx)+0.5) - 1;
        if (in<0)    in = 0;
        if (in>Ncx-1)in = Ncx-1;
        grid_topo       = (topo->b[in] + topo->a[in] * ( topo_chain->x[k] ));
        if ( topo_chain->z[k] > grid_topo ) {
            topo_chain->z[k]   = grid_topo;
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddPartSed( markers *particles, mat_prop materials, markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh) {
    
    int sed_phase = model.surf_ised1;
    int finite_strain = model.finite_strain;
    int track_T_P_x_z = model.track_T_P_x_z;
    int time_My = floor(model.time*scaling.t / (3600.0*365.25*24.0*1.0e6));
    if ( time_My % 2 > 0 ) sed_phase = model.surf_ised1;
    else                   sed_phase = model.surf_ised2;
    
    int    k, Ncx=model.Nx-1, res = 2, ip;
    double distance, dx=model.dx, xW, xE, xmin = model.xmin;
    int nb_part = 2, in, new_ind = particles->Nb_part;
    double xnew, znew, hm0, hm;
    
    for ( k=0; k<res*Ncx; k++ ) {
        
        // West/East cell boundaries
        xW = xmin + k*dx/res;
        xE = xmin + k*dx/res + dx/res;
        
        for ( ip=0; ip<nb_part; ip++ ) {
            
            if (ip == 0) xnew = xW + dx/res/nb_part/2.0;
            if (ip == 1) xnew = xE - dx/res/nb_part/2.0;
            
            // Index of the coarse grid column
            distance        = xnew - model.xmin - dx/2.0;
            in              = ceil((distance/dx)+0.5) - 1;
            if (in<0)    in = 0;
            if (in>Ncx-1)in = Ncx-1;
            
            
            // Topography computed from the coarse grid column
            hm0  = (topo->b0[in] + topo->a0[in] * ( xnew ));
            hm   = (topo->b[in]  + topo->a[in]  * ( xnew ));
            
            if ( hm > hm0 ) { //;# && (hm-hm0)>1e-6
                
                // Compute new altitude of the marker right in between the new and old position of the free surface
                znew                      = 0.5*(hm+hm0);
                particles->x[new_ind]     = xnew;
                particles->z[new_ind]     = znew;
                particles->phase[new_ind] = sed_phase;
                
                particles->sxxd[new_ind]          =  0.0;
                particles->sxz[new_ind]           =  0.0;
                particles->Vx[new_ind]            =  0.0;
                particles->Vz[new_ind]            =  0.0;
                particles->strain[new_ind]        =  0.0;
                particles->strain_el[new_ind]     =  0.0;
                particles->strain_pl[new_ind]     =  0.0;
                particles->strain_pwl[new_ind]    =  0.0;
                particles->strain_exp[new_ind]    =  0.0;
                particles->strain_lin[new_ind]    =  0.0;
                particles->strain_gbs[new_ind]    =  0.0;
                particles->d[new_ind]             =  materials.gs_ref[sed_phase];
                particles->T[new_ind]             =  zeroC/scaling.T;
                particles->P[new_ind]             =  0.0;
                
                particles->phi[new_ind]           =  0.0;
                particles->X[new_ind]             =  0.0;
                particles->sxxd[new_ind]          =  0.0;
                particles->sxz[new_ind]           =  0.0;
                
                if (finite_strain==1) {
                    particles->Fxx[new_ind]           = 1.0;
                    particles->Fxz[new_ind]           = 0.0;
                    particles->Fzx[new_ind]           = 1.0;
                    particles->Fzz[new_ind]           = 0.0;
                }
                
                if (track_T_P_x_z==1) {
                    particles->T0[new_ind]           = zeroC/scaling.T;
                    particles->P0[new_ind]           = 0.0;
                    particles->x0[new_ind]           = particles->x[new_ind];
                    particles->z0[new_ind]           = particles->z[new_ind];
                    particles->Tmax[new_ind]         = zeroC/scaling.T;
                    particles->Pmax[new_ind]         = 0.0;
                }
                
                new_ind++;
            }
        }
    }
    particles->Nb_part = new_ind;
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void RemeshMarkerChain( markers *topo_chain, surface *topo, params model, scale scaling, grid *mesh, int mode ) {
    
    int    k, Ncx=model.Nx-1,count=0, ip, ic, fact=4;
    double dxm=model.dx/(fact+1), xmin = model.xmin + model.dx/2.0;
    int    Nb_part0=topo_chain->Nb_part, *NumMarkCell;
    double distance, dx=model.dx, grid_topo, mismax = 0.01;
    int    in, minPartCell=4, NewInd, inc = 0;
    int    res = 2;
    
    printf("Remesh surface markers, step %d \n", mode);
  
    // Here the marker chain is not remeshed but new marker points are added in deficient locations
    NumMarkCell = DoodzCalloc( res*Ncx, sizeof(int) );
    
    if ( mode == 1 ) {
        
        int nout=0;
        // Find to which cell each marker contribute / find number of topo. markers per FINE cell column (DX/res)
        for (k=0;k<topo_chain->Nb_part;k++) {
            
            if (topo_chain->x[k]>model.xmax || topo_chain->x[k]<model.xmin  ) topo_chain->phase[k] = -1;
            else topo_chain->phase[k]=0;
            
            // Index of the fine grid column
            distance             = topo_chain->x[k] - (model.xmin + dx/2.0/res);
            in                   = ceil((distance/dx*res)+0.5) - 1;
            if (in<0        ) in = 0;
            if (in>res*Ncx-1) in = res*Ncx-1;
            if (topo_chain->phase[k]!=-1) NumMarkCell[in]++;
            // NEW
            if (topo_chain->phase[k]== -1) nout++;
        }
        
        // NEW
        int ii=0;
        int *reuse = DoodzCalloc( nout, sizeof(int) );
        for (k=0;k<topo_chain->Nb_part;k++) {
            if (topo_chain->phase[k]== -1) {
                reuse[ii] = k;
                ii++;
            }
        }
        int recycle = 0;
        
        ii = 0;
        if (nout>0) recycle = 1;
        printf("%d surface markers are out, so recycle is %d\n", nout, recycle);
        
        for ( k=0; k<res*Ncx; k++ ) {
            
            if ( NumMarkCell[k]<minPartCell ) {
                //                printf("Adding topo. markers in cell %d who owns %d markers --  xc = %2.4lf\n", k, NumMarkCell[k], (model.xmin + k*dx/res + dx/2/res)*scaling.L);
                // check for possibility of adding markers
                if( topo_chain->Nb_part+1<topo_chain->Nb_part_max && topo_chain->Nb_part+2<topo_chain->Nb_part_max) {
                    // Add one particle on the WEST side of the fine column
                    if (recycle==0) NewInd                = topo_chain->Nb_part;
                    if (recycle==0) topo_chain->Nb_part++;
                    if (recycle==1) NewInd                = reuse[ii];
                    if (recycle==1) ii++;
                    if (  ii>=nout) recycle=0;
                    // Local x coordinate:  1.0*dx/4.0/res
                    topo_chain->x[NewInd] = model.xmin + k*dx/res + dx/4.0/res;
                    // Index of the coarse grid column
                    distance        = (topo_chain->x[NewInd]-model.xmin-dx/2.0);
                    in              = ceil((distance/dx)+0.5) - 1;
                    if (in<0)    in = 0;
                    if (in>Ncx-1)in = Ncx-1;
                    // Topography computed from the coarse grid column
                    topo_chain->Vx[NewInd] = 0.5*(topo->vx[in] + topo->vx[in+1]);
                    topo_chain->Vz[NewInd] = 0.5*(topo->vz[in] + topo->vz[in+1]);
                    topo_chain->z[NewInd]  = topo->b[in]  + topo->a[in]  * topo_chain->x[NewInd];
                    topo_chain->z0[NewInd] = topo->b0[in] + topo->a0[in] * topo_chain->x[NewInd];

                    // Add one particle on the EAST side of the fine column
                    if (recycle==0) NewInd                = topo_chain->Nb_part;
                    if (recycle==0) topo_chain->Nb_part++;
                    if (recycle==1) NewInd                = reuse[ii];
                    if (recycle==1) ii++;
                    if (  ii>=nout) recycle=0;
                    // Local x coordinate:  3.0*dx/4.0/res
                    topo_chain->x[NewInd] = model.xmin + k*dx/res + 3.0*dx/4.0/res;
                    // Index of the coarse grid column
                    distance        = (topo_chain->x[NewInd]-model.xmin-dx/2.0);
                    in              = ceil((distance/dx)+0.5) - 1;
                    if (in<0)    in = 0;
                    if (in>Ncx-1)in = Ncx-1;
                    // Topography computed from the coarse grid column
                    topo_chain->Vx[NewInd] = 0.5*(topo->vx[in] + topo->vx[in+1]);
                    topo_chain->Vz[NewInd] = 0.5*(topo->vz[in+1] + topo->vz[in+2]);
                    topo_chain->z[NewInd]  = topo->b[in]  + topo->a[in]  * topo_chain->x[NewInd];
                    topo_chain->z0[NewInd] = topo->b0[in] + topo->a0[in] * topo_chain->x[NewInd];
                }
                else {
                    printf("The max. number of topographic particles (currently %d) needs to be increased (number of particles %d)\n", topo_chain->Nb_part_max, topo_chain->Nb_part);
                    exit(45);
                }
            }
        }
        DoodzFree(reuse);
    }
    
    if (mode==2) {
        // Check if marker topography is too different than grid topography - reset to topo
#pragma omp parallel for shared( topo_chain, topo ) private( k, distance, in, grid_topo ) firstprivate( model, dx, Ncx, mismax ) reduction(+ :inc)
        for (k=0;k<topo_chain->Nb_part;k++) {
            // Index of the coarse grid column
            distance        = (topo_chain->x[k]-model.xmin-dx/2.0);
            in              = ceil((distance/dx)+0.5) - 1;
            if (in<0)    in = 0;
            if (in>Ncx-1)in = Ncx-1;
            grid_topo       = (topo->b[in] + topo->a[in] * ( topo_chain->x[k] ));
            //if ( fabs((grid_topo-topo_chain->z[k]) / grid_topo) > mismax) {
                topo_chain->z[k]  = grid_topo;
                topo_chain->z0[k] = grid_topo;
                inc++;
            //}
        }
        printf("Had to correct %d marker topographies for a mismax of %lf\n", inc, mismax);
    }
    
    printf( "Surface remesher 1: old number of marker %d --> New number of markers %d \n", Nb_part0, topo_chain->Nb_part );
    DoodzFree( NumMarkCell );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6 TODO: find a more explicit name like `InterpolateTopographyMarker2Grid`
void ProjectTopography( surface *topo, markers *topo_chain, params model, grid mesh, scale scaling, double* X_vect, int itp_type ) { 

    int k, in, Nx=mesh.Nx;
    double dx=mesh.dx, distance, dxm, mark_val, *Xc_virtual, *Wm, *BmWm;
    Wm              = DoodzCalloc ( Nx-1, sizeof(double));
    BmWm            = DoodzCalloc ( Nx-1, sizeof(double));
    double *heightc = DoodzCalloc ( Nx-1, sizeof(double));
    double *Delta_heightc = DoodzCalloc ( Nx-1, sizeof(double));
    
    for (k=0; k<topo_chain->Nb_part; k++) {

        if ( topo_chain->phase[k] != -1 ) {
            distance  = ( topo_chain->x[k] - mesh.xc_coord[0] );
            in        = ceil( (distance/dx) + 0.5) - 1;
            dxm       = 2.0*fabs( mesh.xc_coord[in] - topo_chain->x[k]);
            if (model.topo_update==0) mark_val  = topo_chain->z[k];
            if (model.topo_update==1) mark_val  = topo_chain->z[k] - topo_chain->z0[k];
            Wm[in]   += 1.0 - dxm/dx;
            BmWm[in] += mark_val*(1.0-dxm/dx);
        }
    }

    for (k=0;k<Nx-1;k++) {
        if (Wm[k]<1e-30){
            Delta_heightc[k] = 0.0;
        }
        else {
            Delta_heightc[k] = BmWm[k]/Wm[k];
        }
    }

    for (k=1;k<Nx-1;k++) {
        if (model.topo_update==0) topo->height[k] = 0.5*(Delta_heightc[k] + Delta_heightc[k-1]);
        if (model.topo_update==1) topo->height[k] = topo->height0[k] + 0.5*(Delta_heightc[k] + Delta_heightc[k-1]);
    }
   topo->height[0]    = topo->height[1];
   topo->height[Nx-1] = topo->height[Nx-2];

    // Correct for sides is the box in case of inflow conditions
    for (k=0;k<Nx;k++) {
        if ( model.polar==0 && model.pure_shear_ALE <= 0 && k==Nx-1 ) topo->height[k]=topo->height[k-1];
        if ( model.polar==0 && model.pure_shear_ALE <= 0 && k==0    ) topo->height[k]=topo->height[k+1];
    }

    // Correct for sides when working in polar mode
    double Rad=6370e3/scaling.L, zW, zE;
    if ( model.polar==1 ) {
        zW = sqrt((Rad - mesh.xg_coord[   0] )*(Rad + mesh.xg_coord[   0]));
        zE = sqrt((Rad - mesh.xg_coord[Nx-1] )*(Rad + mesh.xg_coord[Nx-1]));
        topo->height[   0] = zW;
        topo->height[Nx-1] = zE;
    }

    // Free memory
    DoodzFree(Wm);
    DoodzFree(BmWm);
    DoodzFree(heightc);
    DoodzFree(Delta_heightc);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MarkerChainPolyFit( surface *topo, markers *topo_chain, params model, grid mesh ) {
    
    // Create a cell-based representation of the topography
    int ic;
    int ncx = model.Nx-1;
    double dx = model.dx;
    
    // Linear polynomial for each cell based on vertex values of topography (y = a.x + b)
    for ( ic=0; ic<ncx; ic++ ) {
        // Save old topography
        topo->a0[ic] =  (topo->height0[ic+1] - topo->height0[ic]) / dx;
        topo->b0[ic] =  topo->height0[ic] - (mesh.xg_coord[ic]) * topo->a0[ic];
        // Find slope
        topo->a[ic]  = (topo->height[ic+1] - topo->height[ic]) / dx;
        // Find origin value
        topo->b[ic]  =  topo->height[ic] - (mesh.xg_coord[ic]) * topo->a[ic];
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AllocateMarkerChain( surface *topo, markers* topo_chain, params model ) {
    
    topo_chain->Nb_part_max = 50*model.Nx;
    topo_chain->x           = DoodzMalloc( topo_chain->Nb_part_max*sizeof(DoodzFP) );
    topo_chain->z           = DoodzMalloc( topo_chain->Nb_part_max*sizeof(DoodzFP) );
    topo_chain->z0           = DoodzMalloc( topo_chain->Nb_part_max*sizeof(DoodzFP) );
    topo_chain->Vx          = DoodzCalloc( topo_chain->Nb_part_max, sizeof(DoodzFP) );
    topo_chain->Vz          = DoodzCalloc( topo_chain->Nb_part_max, sizeof(DoodzFP) );
    topo_chain->phase       = DoodzCalloc( topo_chain->Nb_part_max, sizeof(int) );
    // for (int i=1; i<topo_chain->Nb_part_max; i++) topo_chain->phase[i] = -1; 
    topo->height            = DoodzCalloc( (model.Nx),sizeof(DoodzFP) );
    topo->height0           = DoodzCalloc( (model.Nx),sizeof(DoodzFP) );
    topo->vx                = DoodzCalloc( (model.Nx),sizeof(DoodzFP) );
    topo->vz                = DoodzCalloc( (model.Nx+1),sizeof(DoodzFP) );
    topo->a                 = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->b                 = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->a0                = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->b0                = DoodzCalloc( (model.Nx-1),sizeof(DoodzFP) );
    topo->VertInd           = DoodzCalloc( model.Nx,sizeof(DoodzFP) );
    printf( "Marker chain for topography was allocated, %d\n", topo_chain->Nb_part_max );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void FreeMarkerChain( surface *topo, markers* topo_chain ) {
    
    DoodzFree( topo_chain->x );
    DoodzFree( topo_chain->z );
    DoodzFree( topo_chain->z0 );
    DoodzFree( topo_chain->Vx );
    DoodzFree( topo_chain->Vz );
    DoodzFree( topo_chain->phase );
    DoodzFree( topo->height );
    DoodzFree( topo->height0 );
    DoodzFree( topo->a  );
    DoodzFree( topo->b  );
    DoodzFree( topo->a0 );
    DoodzFree( topo->b0 );
    DoodzFree( topo->vx );
    DoodzFree( topo->vz );
    DoodzFree( topo->VertInd );
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double TopoFun0( double X, double H, double L, scale scaling ) {
    double Y = -10.0e3/scaling.L;
    Y = 0 + 2000/scaling.L*cos(X*M_PI/L*6.0);
    Y = -5000.0/scaling.L + 4500.0/scaling.L*cos(X*M_PI/L*1);
    return Y;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

double TopoFun( double X, int ic, surface topo, scale scaling ) {
    double Y;
    Y = topo.b[ic] + topo.a[ic] * X;
    return Y;
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CellFlagging( grid *mesh, params model, surface topo, scale scaling ) {
    
    int nx   = mesh->Nx;
    int nz   = mesh->Nz;
    int ncx  = nx-1;
    int ncz  = nz-1;
    int nxvz = nx+1;
    int nzvx = nz+1;
    int i, j, c1, c2;
    double h;
    int *PVtag0, *PVtag;
    PVtag0 = DoodzCalloc( (ncx+2)*(ncz+2), sizeof(int) );
    PVtag  = DoodzCalloc( (ncx+2)*(ncz+2), sizeof(int) );
    
    Initialise1DArrayChar(   mesh->BCp.type, ncx*ncz, -1  );
    Initialise1DArrayChar(   mesh->BCt.type,    ncx*ncz, -1  );
    Initialise1DArrayChar(   mesh->BCu.type, nx*nzvx, -1  );
    Initialise1DArrayChar(   mesh->BCv.type, nxvz*nz, -1  );
    Initialise1DArrayChar(   mesh->BCg.type,    nx *nz , -1  );
    
    Initialise1DArrayDouble( mesh->BCt.val,     ncx*ncz, 0.0 );
    Initialise1DArrayDouble( mesh->BCp.val,  ncx*ncz, 0.0 );
    Initialise1DArrayDouble( mesh->BCu.val,  nx*nzvx, 0.0 );
    Initialise1DArrayDouble( mesh->BCv.val,  nxvz*nz, 0.0 );
    Initialise1DArrayDouble( mesh->BCg.val,     nx *nz , 0.0  );
    
    //---------------------------- METHOD 1
    
    //------------- FLAG EXTENDED PRESSURE CELLS -------------//
    for( i=0; i<ncx+2; i++ ) {
        
        // Topographic function
        if (i==0) {
            h = TopoFun( mesh->xvz_coord[i], 0, topo, scaling );
        }
        
        if (i==ncx+1) {
            h = TopoFun( mesh->xvz_coord[i], ncx-1, topo, scaling );
        }
        
        if (i>0 && i<ncx+1) {
            h = TopoFun( mesh->xvz_coord[i], i-1, topo, scaling );
        }
        
        //-------------------------------------------------------------//
        
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
            PVtag0[c1] = -1;
            if (mesh->zvx_coord[j]>h) {
                // Deactivate pressure cells over the surface
                PVtag0[c1] = 30;
            }
            
        }
    }
    
    //---------------------------- METHOD 2
    for( i=0; i<ncx+2; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
            
            PVtag[c1] = PVtag0[c1];
            
            if ( (j>0 && j<ncz+1) && (PVtag0[c1] != PVtag0[c1-(ncx+2)]) ) {
                PVtag[c1] = 60;
            }
            
        }
    }
    
    ArrayEqualArrayI(PVtag0, PVtag, (ncx+2)*(ncz+2));
    
    for( i=1; i<ncx+1; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
            
            if ( PVtag0[c1]==-1 && (PVtag0[c1-1]==60 && PVtag0[c1+1]==60) ) {
                PVtag[c1] = 60;
            }
            
            // Correction at Lausanne Gare - this is important
            // It corrects for large slopes, Stokes matrix would become asymmetric otherwise
            if ( PVtag0[c1]==30 && PVtag0[c1+1]==-1 && (PVtag0[c1-1]==30 || PVtag0[c1-1]==60 ) ) {
                PVtag[c1] = 60;
            }
            
            if ( PVtag0[c1]==30 && PVtag0[c1-1]==-1 && (PVtag0[c1+1]==30 || PVtag0[c1+1]==60) ) {
                PVtag[c1] = 60;
            }
            
            // Additional fix on 28/02/19
            if ( PVtag0[c1]==30 && PVtag0[c1-1]==-1 && PVtag0[c1+1]==-1) {
                PVtag[c1] = 60;
            }
        }
    }
    
    // Sides
    for( j=0; j<ncz+2; j++ ) {
        
        c1 = 0 + j*(ncx+2);
        c1 = ncx+1 + j*(ncx+2);
        
        PVtag[c1] = PVtag[c1-1];
    }
    
    for( j=0; j<ncz; j++ ) {
        for( i=0; i<ncx; i++ ) {
            
            c1 = i + j*(ncx);
            c2 = i + 1 + (j+1)*(ncx+2);
            
            mesh->BCp.type[c1]     = PVtag[c2];
            mesh->BCt.type[c1]     = PVtag[c2];
            mesh->BCp_exp.type[c2] = PVtag[c2];

            
            if ( PVtag[c2] == 30 ) {
                mesh->BCt.type[c1]     = 30;
                mesh->T[c1]            = zeroC/scaling.T ;
                mesh->BCp_exp.type[c2] = 30;
            }
            
            // Above surface pressure nodes
            if ( PVtag[c2] == 60 ) {
                mesh->BCp.type[c1] = 31;
                mesh->BCp.val[c1]  = 0.0*7.0e6/scaling.S;
                
                // TEMPERATURE above THE SURFACE
                mesh->BCt.type[c1] = 30;
                mesh->T[c1]        = zeroC/scaling.T ;
                // activate first layer above surface for interpolation
                mesh->BCp_exp.type[c2] = -1;
                
            }
        }
    }
    
    //------------- FLAG Vx CELLS -------------//
    for ( i=0; i<nx; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            c2 = i + (j)*(nx);
            mesh->BCu.type[c2] = 30;
            mesh->BCu.val[c2]  = 0.0;
        }
    }
    
    for ( i=0; i<ncx+2; i++ ) {
        for( j=0; j<ncz+2; j++ ) {
            
            c1 = i + j*(ncx+2);
            c2 = i + (j)*(nx);
            
            if ( i<ncx+1 && PVtag[c1]==-1 ) {
                mesh->BCu.type[c2] = -1;
                mesh->BCu.val[c2]  = 0.0;
            }
            
            if ( i>0 && PVtag[c1]==-1) {
                mesh->BCu.type[c2-1] = -1;
                mesh->BCu.val[c2-1]  = 0.0;
            }
        }
    }
    
    for (j=0; j<nzvx; j++) {
        
        c1 = 0 + j*(nx);
        c2 = nx-1 + j*(nx);
        mesh->BCu.type[c1] = mesh->BCu.type[c1+1];
        mesh->BCu.type[c2] = mesh->BCu.type[c2-1];
    }
    
    //------------- FLAG Vz CELLS -------------//
    for (i=0; i<ncx+2; i++) {
        for (j=0; j<nz; j++) {
            
            c2 = i + j*(nxvz);
            mesh->BCv.type[c2] = 30;
            mesh->BCv.val[c2] = 0.0;
        }
    }
    
    for (i=0; i<ncx+2; i++) {
        for (j=0; j<ncz+2; j++) {
            
            c1 = i + j*(ncx+2);
            c2 = i + j*(nxvz);
            
            if ( j<ncz+1 && PVtag[c1]==-1 ) {
                mesh->BCv.type[c2] = -1;
                mesh->BCv.val[c2]  = 0.0;
            }
            
            if ( j>0  && PVtag[c1]==-1 ) {
                mesh->BCv.type[c2-nxvz] = -1;
                mesh->BCv.val[c2-nxvz]  = 0.0;
            }
        }
    }
    
    for (j=0; j<nz; j++) {
        
        c1 = 0 + j*(nxvz);
        c2 = nxvz-1 + j*(nxvz);
        
        mesh->BCv.type[c1] = mesh->BCv.type[c1+1];
        mesh->BCv.type[c2] = mesh->BCv.type[c2-1];
    }
    
    //------------- FLAG VERTICES -------------//
    
    for( j=0; j<nz; j++ ) {
        for ( i=0; i<nx; i++ ) {
            
            c1 = i + j*nx;
            mesh->BCg.type[c1] = -1;
            
            if (mesh->zg_coord[j] > topo.height[i]) {
                
                mesh->BCg.type[c1] = 30;
            }
        }
    }
    
    // set surface vertices to zero vel
    for ( i=0; i<nx; i++ ) {
        for( j=0; j<nz; j++ ) {
            
            c1 = i + j*nx;
            
            if ( j<nz-1 ) {
                
                if ( mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+nx] == 30 ) mesh->BCg.type[c1] = 30;
            }
        }
    }
    
    // set surface vertices to zero vel
    for ( i=0; i<nx; i++ ) {
        for( j=0; j<nz; j++ ) {
            
            c1 = i + j*nx;
            
            if ( i>0 ) {
                
                if ( mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1-1] == 30 ) mesh->BCg.type[c1] = 30;
            }
        }
    }
    
    // set surface vertices to zero vel
    for ( i=0; i<nx; i++ ) {
        for( j=0; j<nz; j++ ) {
            
            c1 = i + j*nx;
            
            if ( i<nx-1 ) {
                
                if ( mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+1] == 30 ) mesh->BCg.type[c1] = 30;
            }
        }
    }
    
    // This is suspected to create issues by forcing vertices to be activated whilie there are no surrounding particles (regions of intense cusping)
    // METHOD 2 Vertices
    for( j=0; j<nz; j++ ) {
        for ( i=0; i<nx; i++ ) {
            
            c1 = i + j*nx;
            c2 = i + j*ncx;
            
            if ( i>0 && i<nx-1 && j>0 && j<nz-1 && mesh->BCg.type[c1] == 30 ) {
                
                if ( mesh->BCp.type[c2-1] == -1  && mesh->BCp.type[c2-1-ncx] == -1 && mesh->BCp.type[c2] == -1 && mesh->BCp.type[c2-ncx] == -1 ) {
                    mesh->BCg.type[c1] = -1;
                }
                
            }
            
            if ( i==0 && i<nx-1 && j>0 && j<nz-1 && mesh->BCg.type[c1] == 30 ) {
                
                if (  mesh->BCp.type[c2] == -1 && mesh->BCp.type[c2-ncx] == -1 ) {
                    mesh->BCg.type[c1] = -1;
                }
                
            }
            
            if ( i==nx-1 && j>0 && j<nz-1 && mesh->BCg.type[c1] == 30 ) {
                
                if ( mesh->BCp.type[c2-1] == -1  && mesh->BCp.type[c2-1-ncx] == -1 ) {
                    mesh->BCg.type[c1] = -1;
                }
                
            }
            
            
        }
    }
    
    DoodzFree( PVtag  );
    DoodzFree( PVtag0 );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// MD6
void CleanUpSurfaceParticles( markers* particles, grid *mesh, surface topo, scale scaling ) {

    int k, ic, ncx=mesh->Nx-1;
    double h;
    double dx=mesh->dx, dz=mesh->dz;
    double xmin = mesh->xg_coord[0] + dx/2;
    double dst;

    int count = 0;

#pragma omp parallel for shared ( particles, topo ) private ( k, h, ic, dst ) firstprivate( ncx, xmin, dx )
    for ( k=0; k<particles->Nb_part; k++ ) {

        if ( particles->phase[k] != -1 ) {

            // Get the column:
            dst = fabs(particles->x[k] - xmin);
            ic  = ceil((dst/dx)+0.5) - 1;

            if (ic<0) ic = 0;
            if (ic>ncx-1) ic = ncx-1;

            h = topo.b[ic] + topo.a[ic]*particles->x[k];

            if ( particles->z[k]>h ) {
                particles->phase[k] = -1;
            }
        }
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/


void SurfaceDensityCorrection( grid *mesh, params model, surface topo, scale scaling ) {
    
    int nx   = mesh->Nx;
    int nz   = mesh->Nz;
    int ncx  = nx-1;
    int ncz  = nz-1;
    int i, j, c1;
    double h0, h, dz = fabs(mesh->zg_coord[1]-mesh->zg_coord[0]);
    
    // that's MD6 commented for testing
    // Density on cell centers
    for( j=0; j<ncz; j++ ) {
        for( i=0; i<ncx; i++ ) {
            c1 = i + j*(ncx);
            mesh->FreeSurfW_n[c1] = 1.0;
            
            if (mesh->BCp.type[c1] == -1 && mesh->BCp.type[c1+ncx] == 31 ) {
                h  = topo.b[i] + topo.a[i]*mesh->xc_coord[i];
                h0               = fabs(h - mesh->zc_coord[j]);
                // double rho_ini = mesh->rho_n[c1];
                mesh->rho_n[c1]  *= h0/dz;
                mesh->rho0_n[c1] *= h0/dz;
                mesh->FreeSurfW_n[c1] = h0/dz;
                // printf("rho_ini = %3.3lf --- rho = %3.3lf -- h = %3.3lf --- h0/dz = %3.3lf \n", rho_ini*scaling.rho, mesh->rho_n[c1]*scaling.rho, h*scaling.L, h0/dz);
            }
            if ( mesh->BCp.type[c1] == 30 || mesh->BCp.type[c1] == 31 ) {
                mesh->rho_n[c1]       = 1.0/scaling.rho;
                mesh->FreeSurfW_n[c1] = 0.0;
            }
        }
    }

    // Density on cell vertices
    for( j=0; j<nz; j++ ) {
        for( i=0; i<nx; i++ ) {

            c1 = i + j*(nx);
            mesh->FreeSurfW_s[c1] = 1.0;

            if (mesh->BCg.type[c1] == -1 && mesh->BCg.type[c1+nx] == 30) {
                if (i==0)    h  = topo.b[i] + topo.a[i]*mesh->xc_coord[i];
                if (i==nx-1) h  = topo.b[nx-2] + topo.a[nx-2]*mesh->xc_coord[nx-2];
                if (i>0 && i<nx-1) {
                    h  = 0.5*(topo.b[i] + topo.a[i]*mesh->xc_coord[i]);
                    h += 0.5*(topo.b[i-1] + topo.a[i-1]*mesh->xc_coord[i-1]);
                }
                h0               = fabs(h - mesh->zg_coord[j]);
                mesh->rho_s[c1] *= h0/dz;
                mesh->FreeSurfW_s[c1] = h0/dz;
            }
            if ( mesh->BCg.type[c1] == 30 ) {
                mesh->rho_s[c1]     = 1.0/scaling.rho;
                mesh->FreeSurfW_s[c1] = 0.0;
            }
        }
    }
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void DiffuseAlongTopography( grid *mesh, params model, scale scaling, double *array_ini, double *array, int size, double dummy,  double diff_time ) {
    
    // Explicit diffusion solver;
    int i, it;
    double dx   = fabs(mesh->xg_coord[1]-mesh->xg_coord[0]);
    double diff = model.surf_diff;
    double dt   = 0.4*dx*dx/diff, time=0.0, dtr;
    int nstep   = (int)(diff_time/dt + 1);
    double correct[size], s, e;
    double base_level = model.surf_baselev; // left side
    double sedi_rate  = model.surf_sedirate;
    double Wvalley    = model.surf_Winc;
    double Vinc       = -model.surf_Vinc, Vinc_num;
    
    printf("****** Surface processes ******");
    printf("Going to make %03d substeps for surface processes\n", nstep);
    printf("W valley   = %2.2e m\n", Wvalley*scaling.L);
    printf("Vincision  = %2.2e m.s-1\n", Vinc*scaling.V);
    printf("Kero       = %2.2e m2.s-1\n", diff*(pow(scaling.L,2.0)/scaling.t));
    printf("Sed. rate  = %2.2e m/y with base level: %2.2e m\n", model.surf_sedirate*scaling.V*3600.0*365.0*24.0, base_level*scaling.L);

    if ( model.surface_processes == 1 || model.surface_processes == 5 ) {

        if ( model.surface_processes == 5 ) {
        
            // Compute volume of cells in the valley region
            int ncell = 0;
            for (i=1; i<size-1; i++) {
                if (fabs(mesh->xg_coord[i]) < 0.5*Wvalley){
                    ncell = ncell + 1;
                }
            }

            // Recompute erosion rate to satisfy mass
            Vinc_num = Wvalley*dt*Vinc / (ncell*dx*dt);

            printf("We currently have %0d cell(s) within the valley region\n", ncell);
            printf("Real surface of eroded material should be: %2.2e\n", Wvalley*dt*Vinc);
            printf("Actual surface of eroded material is     : %2.2e\n", ncell*dx*dt*Vinc);
            printf("Corrected surface of eroded material is  : %2.2e\n", ncell*dx*dt*Vinc_num);
        }
        
        // Calculate timestep for diffusion sub-steps
        if (dt >diff_time) dtr = diff_time;
        if (dt<=diff_time) dtr = diff_time/nstep;
        
        // Sub-time loop
        for (it=0; it<nstep; it++) {
                        
            for (i=1; i<size-1; i++) {
                correct[i]  = 0.5*dtr/dx/dx*diff*(array[i-1]+array[i+1]-2.0*array[i]);
            }
            for (i=size-2; i==1; i--) {
                correct[i] += 0.5*dtr/dx/dx*diff*(array[i-1]+array[i+1]-2.0*array[i]);
            }
            for (i=1; i<size-1; i++) {
                // Activate source term (sedimentation) only below base level
                s = 0.0;
                e = 0.0;
                
                if ( model.surface_processes == 5 && fabs(mesh->xg_coord[i] ) < 0.5*Wvalley) {
                    e = dtr*Vinc_num;
                }
                
                if (array[i]<base_level) s = sedi_rate*dtr;
                array[i] = array_ini[i] + correct[i] + s + e;
                
            }
            time += dtr;
        }
        printf("Do %d topographic diffusion steps - whole time: %2.2e s - final time: %2.2e s - explicit dt: %2.2e s - diffusivity: %2.2e m^2/s\n", nstep, diff_time*scaling.t, time*scaling.t, dt*scaling.t, diff*pow(scaling.L,2)/scaling.t );
    }
    
    // Instantaneous basin filling
    if (model.surface_processes == 2) {
        
        for (i=0; i<size; i++) {
            if (array[i]<base_level)
                array[i]  = array_ini[i] + sedi_rate*model.dt;
        }
    }
    
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void SurfaceVelocity( grid *mesh, params model, surface *topo, markers* topo_chain, scale scaling ) {
    
    int nx   = mesh->Nx;
    int nz   = mesh->Nz;
    int nxvz = nx+1;
    int nzvx = nz+1;
    int i, j, c1, c2, k;
    double dx = fabs(mesh->xg_coord[1]-mesh->xg_coord[0]);
    double distance, dxm;
    int in;
    
    double sumvel=0.0, sumvelx=0.0;
    int    ncell, ncellx;
    
    // Build surface velocity vectors from the mesh
    
    // Vx on topography vertices
    for( j=0; j<nzvx; j++ ) {
        for( i=0; i<nx; i++ ) {
            c1 = i + j*(nx);
            if (mesh->BCu.type[c1] == 30 && mesh->BCu.type[c1-nx] != 30 ) {
                topo->vx[i] = mesh->u_in[c1-nx];
            }
        }
    }
   
    ncell  = 0;
    ncellx = 0;
    // Vz on topography vertices
    for( j=0; j<nz; j++ ) {
        for( i=0; i<nxvz; i++ ) {
            c2 = i + j*(nxvz);
            if (mesh->BCv.type[c2] == 30 && mesh->BCv.type[c2-nxvz] != 30 ) {
                topo->vz[i] = mesh->v_in[c2-2*nxvz];
                sumvel += topo->vz[i];
                ncell++;
            }
        }
    }
    
    sumvel = 0.0;
    ncell  = 0;
    ncellx = 0;

    //---------------------------
    // Interpolate velocities from topography nodes to topography markers
    for( k=0; k<topo_chain->Nb_part; k++ ) {
        
        topo_chain->Vx[k] = 0.0;
        distance          = (topo_chain->x[k]-mesh->xg_coord[0]);
        in                = ceil((distance/dx)) - 1;
        if (in<0)    in   = 0;
        if (in>nx-2) in   = nx-2;
        dxm               = topo_chain->x[k] - mesh->xg_coord[in];
        topo_chain->Vx[k] =  (1.0-dxm/dx)* topo->vx[in];
        topo_chain->Vx[k]+=  (dxm/dx)* topo->vx[in+1];
        
        sumvelx += topo_chain->Vx[k];
        ncellx ++;
        //---------------------------
        topo_chain->Vz[k] = 0.0;
        distance          =(topo_chain->x[k]-mesh->xvz_coord[0]);
        in                = ceil((distance/dx)) - 1;
        if (in<0)      in = 0;
        if(in>(nx+1)-2)in = (nx+1)-2;
        dxm               = topo_chain->x[k] - mesh->xvz_coord[in];
        topo_chain->Vz[k] =  (1.0-dxm/dx)* topo->vz[in];
        topo_chain->Vz[k]+=  (dxm/dx)* topo->vz[in+1];
        
        sumvel += topo_chain->Vz[k];
        ncell ++;
    }   
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void KeepZeroMeanTopo(params *model, surface *topo, markers *topo_chain ) {
    printf("Make zero mean topo\n");
    double mean_z = 0.0;
    for (int i=0; i<model->Nx; i++){
        mean_z += topo->height[i];
    }
    mean_z /= model->Nx;
    for (int i=0; i<model->Nx; i++) {
        topo->height[i] -= mean_z;
    }
    mean_z = 0.0;
    int nmark = 0;
    for (int i=0; i<topo_chain->Nb_part_max; i++){
        if (topo_chain->phase[i] > -1) {
        mean_z += topo_chain->z[i];
        nmark  ++;
        }
    }
    mean_z /= nmark;
    for (int i=0; i<topo_chain->Nb_part_max; i++){
        topo_chain->z[i] -= mean_z;
    }
}