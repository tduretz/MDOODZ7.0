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

#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "mdoodz-private.h"
#include "hdf5.h"
#include "zlib.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#ifdef _OMP_
#include "omp.h"
#else
#define omp_get_thread_num() 0
#endif

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void CreateDir(const char *dirName) {
  struct stat st = {0};
  if (stat(dirName, &st) == -1) {
#ifdef _WIN32
    mkdir(dirName);
#else
    mkdir(dirName, 0700);
#endif
  }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Scale parameters to physical scale
void ScaleBack(float* FieldF, double scale, int size) {
    int k;
    for (k=0; k<size; k++) {
        FieldF[k] = FieldF[k]*scale;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Scale parameters to physical scale
void ScaleBackD(double* FieldD, double scale, int size) {
    int k;
    for (k=0; k<size; k++) {
        FieldD[k] = FieldD[k]*scale;
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Function that cast double arrays to float (used for visualisation purposes only)
void DoubleToFloat(double* FieldD, float* FieldF, int size) {
    int k;
    for (k=0; k<size; k++) {
        FieldF[k] = (float)FieldD[k];
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// The following routines were originally writen by D. A. May, and used in I2VIS //
void CreateOutputHDF5( const char FileName[] )
{
	hid_t       file_id;   /* file identifier */

	/* Create a new file using default properties. */
	file_id = H5Fcreate( FileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	/* Terminate access to the file. */
	H5Fclose(file_id);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddGroupToHDF5( const char filename[], const char group[] )
{
	hid_t       file_id, group_id;  /* identifiers */
	char        *group_name;

	asprintf( &group_name, "/%s", group );

	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);

	/* Create group "ParticleGroup" in the root group using absolute name. */
	group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	/* Close group. */
	H5Gclose(group_id);

	/* Close the file. */
	H5Fclose(file_id);

	free( group_name );
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void AddFieldToGroup( const char filename[], const char group[], const char field[], char d_type, int np, void* data, int dim )
{
	hid_t   file_id, particle_group_id, coord_dataset_id, coord_dataspace_id;  /* identifiers */
	hsize_t length;
	char *dataset_name;
	char *group_name;
	hid_t    plist=0,SET_CREATION_PLIST;
	hsize_t  cdims[2] = {0,0};
	double percentage_chunk;
	double chunk_size;
	int deflation_level, quiet=1;
    int compress=_TRUE_;

	asprintf( &group_name, "/%s", group );
	asprintf( &dataset_name, "%s/%s", group, field );

	length = dim * np;

	/* Open exisiting file */
	file_id = H5Fopen( filename, H5F_ACC_RDWR, H5P_DEFAULT);

	/* Open group "ParticleGroup" */
	particle_group_id = H5Gopen(file_id, group_name,H5P_DEFAULT);

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	percentage_chunk = 5.0;
	chunk_size = percentage_chunk*((double)np)/((double)100.0 );
	if( chunk_size < 1 ) {
		cdims[0] = 1;
	}
	else {
		cdims[0] = (int)chunk_size;
	}

	plist  = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(plist, 1, cdims);
	deflation_level = 4;
    H5Pset_deflate( plist, deflation_level);

	/* Create the data space for the dataset. */
	coord_dataspace_id = H5Screate_simple( 1, &length, NULL);

	if(compress==_TRUE_) {
		SET_CREATION_PLIST = plist;
        if ( quiet == 0 ) {
            printf("*** Compression info *** \n");
            printf("  chunk_size = %f \n", chunk_size );
            printf("  deflation level = %d \n", deflation_level );
        }
	}
	else {
		SET_CREATION_PLIST = H5P_DEFAULT;
	}

	if( d_type == 'd' ) {
		/* Create a dataset within "group". */
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_DOUBLE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);

		/* Write the particle dataset. */
		H5Dwrite(coord_dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );
	}

	else if( d_type == 'c' ) {
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_CHAR, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
        H5Dwrite(coord_dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );

	}

	else if( d_type == 'i' ) {
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_STD_I32BE, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
        H5Dwrite(coord_dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );

	}
	else if( d_type == 'f' ) {
		coord_dataset_id = H5Dcreate( file_id, dataset_name, H5T_NATIVE_FLOAT, coord_dataspace_id, H5P_DEFAULT,SET_CREATION_PLIST,H5P_DEFAULT);
        H5Dwrite(coord_dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data );

	}

	else {
		printf("ERROR: Only know how to write doubles (d), ints (i), or chars (c) \n");
		exit(1);
	}

	/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

	/* Close the data space for this particle dataset. */
	//status = H5Sclose(attr_dataspace_id);

	/* Close the particle coord dataset. */
	H5Dclose(coord_dataset_id);

	/* free(dataset_nameP); */
	/* } */

	/* Close datatspace. */
	H5Sclose(coord_dataspace_id);

	/* Close list. */
	H5Pclose(plist);

	/* Close group. */
	H5Gclose(particle_group_id);

	/* Close the file. */
	H5Fclose(file_id);

	free(group_name);
	free(dataset_name);
}
// The above routines were originally writen by D. A. May, and used in I2VIS //

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void WriteOutputHDF5( grid *mesh, markers *particles, surface *topo, markers* topo_chain, params model, Nparams Nmodel, char *txtout, mat_prop materials, scale scaling ) {

    char *FileName;
    double A[8];
    int k;
    double *strain, *strain_el, *strain_pl, *strain_pwl, *strain_exp, *strain_lin, *strain_gbs, *X;
    float *Crho_s, *Crho_n, *Ceta_s, *Ceta_n, *CVx, *CVz, *CP, *Csxxd, *Cszzd, *Csxz, *Cexxd, *Cezzd, *Cexz, *Cstrain, *Cstrain_el, *Cstrain_pl, *Cstrain_pwl, *Cstrain_exp, *Cstrain_lin, *Cstrain_gbs, *CT, *Cd;
    float *Cxg_coord, *Czg_coord, *Cxc_coord, *Czc_coord, *Czvx_coord, *Cxvz_coord;
    float *CeII_el, *CeII_pl, *CeII_pwl, *CeII_exp, *CeII_lin, *CeII_gbs, *CX;
    double *Fxx, *Fxz, *Fzx, *Fzz, *nx, *nz;
    float *CFxx, *CFxz, *CFzx, *CFzz, *Cnx, *Cnz;
    double *T0, *P0, *x0, *z0, *Tmax, *Pmax;
    float *CT0, *CP0, *Cx0, *Cz0, *CTmax, *CPmax;
    float *CXreac;
    float *COverS, *Cdivu, *Cdivu_el, *Cdivu_pl, *Cdivu_th, *Cdivu_r;
        
    int    cent=1, vert=0, prop=1, interp=0;
    int    res_fact = 1;
    int    nxviz, nzviz, nxviz_hr, nzviz_hr;
    char  *compo, *compo_hr;
    float *Cxviz, *Czviz, *Cxviz_hr, *Czviz_hr, *Cxtopo, *Cztopo, *Cheight, *Ctopovx, *Ctopovz, *Ctopovx_mark, *Ctopovz_mark;
    double *P_total;
    float  *Ccohesion, *Cfriction, *Cani_fac;

    P_total  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));

    // Build total pressure
    for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
        P_total[k] = mesh->p_in[k];//(mesh->p_in[k] -  mesh->p0_n[k]);
    }

    // int N= 5;//(mesh->Nx-0)*(mesh->Nz-0);
    // for (int c1=0; c1<N; c1++) {
    //      printf("eta = %2.2e\n", mesh->sxz[c1]/2.0/mesh->exz[c1]);
    // }

    // ---------------------------------------------------------
    // Genrate phase map with normal resolution
    res_fact = 1;
    nxviz = res_fact*(mesh->Nx-1) + 1;
    nzviz = res_fact*(mesh->Nz-1) + 1;
    double xviz[nxviz], zviz[nzviz];
    // Define visual grid
    xviz[0] = mesh->xg_coord[0];
    zviz[0] = mesh->zg_coord[0];
    for (k=1; k<nxviz; k++) xviz[k] = xviz[k-1] + mesh->dx/res_fact;
    for (k=1; k<nzviz; k++) zviz[k] = zviz[k-1] + mesh->dz/res_fact;
    compo  = DoodzMalloc( sizeof(char)*(nxviz-1)*(nzviz-1));
    // Closest point interpolation: marker phase --> visual nodes
    Interp_Phase2VizGrid ( *particles, particles->dual, mesh, compo, xviz, zviz, nxviz, nzviz, model, *topo );

    // ---------------------------------------------------------
    // Genrate phase map with double resolution
    res_fact = 2;
    nxviz_hr = res_fact*(mesh->Nx-1) + 1;
    nzviz_hr = res_fact*(mesh->Nz-1) + 1;
    double xviz_hr[nxviz_hr], zviz_hr[nzviz_hr];
    // Define visual grid
    xviz_hr[0] = mesh->xg_coord[0];
    zviz_hr[0] = mesh->zg_coord[0];
    for (k=1; k<nxviz_hr; k++) xviz_hr[k] = xviz_hr[k-1] + mesh->dx/res_fact;
    for (k=1; k<nzviz_hr; k++) zviz_hr[k] = zviz_hr[k-1] + mesh->dz/res_fact;
    
    compo_hr  = DoodzMalloc( sizeof(char)*(nxviz_hr-1)*(nzviz_hr-1));
    // Closest point interpolation: marker phase --> visual nodes
    Interp_Phase2VizGrid ( *particles, particles->dual, mesh, compo_hr, xviz_hr, zviz_hr, nxviz_hr, nzviz_hr, model, *topo );
    // ---------------------------------------------------------
    // Smooth rheological contrasts
    if (model.diffuse_X) P2Mastah( &model, *particles, particles->X,     mesh, mesh->X_n,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    // ---------------------------------------------------------
    // Cast grid arrays
    Crho_s  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->rho_s, Crho_s, model.Nx*model.Nz);
    ScaleBack( Crho_s, scaling.rho, model.Nx*model.Nz );

    Crho_n  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->rho_n, Crho_n, (model.Nx-1)*(model.Nz-1));
    ScaleBack( Crho_n, scaling.rho, (model.Nx-1)*(model.Nz-1));

    Ceta_s  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->eta_phys_s, Ceta_s, model.Nx*model.Nz);
    ScaleBack( Ceta_s, scaling.eta, model.Nx*model.Nz );

    Ceta_n  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eta_phys_n, Ceta_n, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Ceta_n, scaling.eta, (model.Nx-1)*(model.Nz-1) );

    CVx  = DoodzMalloc( sizeof(float)*model.Nx*(model.Nz+1));
    DoubleToFloat( mesh->u_in, CVx, model.Nx*(model.Nz+1) );
    ScaleBack( CVx, scaling.V, model.Nx*(model.Nz+1) );

    CVz  = DoodzMalloc( sizeof(float)*(model.Nx+1)*model.Nz);
    DoubleToFloat( mesh->v_in, CVz, (model.Nx+1)*model.Nz );
    ScaleBack( CVz, scaling.V, (model.Nx+1)*model.Nz );

    CP  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( P_total, CP, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CP, scaling.S, (model.Nx-1)*(model.Nz-1) );

    Csxxd  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->sxxd, Csxxd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Csxxd, scaling.S, (model.Nx-1)*(model.Nz-1) );
    
    Cszzd  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->szzd, Cszzd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cszzd, scaling.S, (model.Nx-1)*(model.Nz-1) );

    Csxz  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->sxz, Csxz, model.Nx*model.Nz );
    ScaleBack( Csxz, scaling.S, model.Nx*model.Nz );

    Cexxd  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->exxd, Cexxd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cexxd, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cezzd  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->ezzd, Cezzd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cezzd, scaling.E, (model.Nx-1)*(model.Nz-1) );

    Cexz  = DoodzMalloc( sizeof(float)*model.Nx*model.Nz);
    DoubleToFloat( mesh->exz, Cexz, model.Nx*model.Nz );
    ScaleBack( Cexz, scaling.E, model.Nx*model.Nz );

    CeII_el  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_el, CeII_el, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_el, scaling.E, (model.Nx-1)*(model.Nz-1) );

    CeII_pl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_pl, CeII_pl, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_pl, scaling.E, (model.Nx-1)*(model.Nz-1) );

    CeII_pwl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_pwl, CeII_pwl, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_pwl, scaling.E, (model.Nx-1)*(model.Nz-1) );

    CeII_exp  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_exp, CeII_exp, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_exp, scaling.E, (model.Nx-1)*(model.Nz-1) );

    CeII_lin  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_lin, CeII_lin, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_lin, scaling.E, (model.Nx-1)*(model.Nz-1) );

    CeII_gbs  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->eII_gbs, CeII_gbs, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CeII_gbs, scaling.E, (model.Nx-1)*(model.Nz-1) );

    Cd        = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->d_n, Cd, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cd, scaling.L, (model.Nx-1)*(model.Nz-1) );

    CX        = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->X_n, CX, (model.Nx-1)*(model.Nz-1) );

    COverS    = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->OverS_n, COverS, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( COverS, scaling.S, (model.Nx-1)*(model.Nz-1) );
    
    Cdivu        = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->div_u, Cdivu, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cdivu, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cdivu_el     = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->div_u_el, Cdivu_el, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cdivu_el, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cdivu_pl     = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->div_u_pl, Cdivu_pl, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cdivu_pl, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cdivu_th     = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->divth0_n, Cdivu_th, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cdivu_th, scaling.E, (model.Nx-1)*(model.Nz-1) );
    
    Cdivu_r      = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->div_u_r , Cdivu_r , (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Cdivu_r , scaling.E, (model.Nx-1)*(model.Nz-1) );

    //---------------------------------------------------
    CT  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->T, CT, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( CT, scaling.T, (model.Nx-1)*(model.Nz-1) );

    //---------------------------------------------------
    Cxg_coord = DoodzMalloc( sizeof(float)*model.Nx);
    DoubleToFloat( mesh->xg_coord, Cxg_coord, model.Nx );
    ScaleBack( Cxg_coord, scaling.L, model.Nx );

    Czg_coord = DoodzMalloc( sizeof(float)*model.Nz);
    DoubleToFloat( mesh->zg_coord, Czg_coord, model.Nz );
    ScaleBack( Czg_coord, scaling.L, model.Nz );

    Cxc_coord = DoodzMalloc( sizeof(float)*(model.Nx-1));
    DoubleToFloat( mesh->xc_coord, Cxc_coord, model.Nx-1 );
    ScaleBack( Cxc_coord, scaling.L, model.Nx-1 );

    Czc_coord = DoodzMalloc( sizeof(float)*(model.Nz-1));
    DoubleToFloat( mesh->zc_coord, Czc_coord, model.Nz-1 );
    ScaleBack( Czc_coord, scaling.L, model.Nz-1 );

    Czvx_coord = DoodzMalloc( sizeof(float)*(model.Nz+1));
    DoubleToFloat( mesh->zvx_coord, Czvx_coord, model.Nz+1 );
    ScaleBack( Czvx_coord, scaling.L, model.Nz+1 );

    Cxvz_coord = DoodzMalloc( sizeof(float)*(model.Nx+1));
    DoubleToFloat( mesh->xvz_coord, Cxvz_coord, model.Nx+1 );
    ScaleBack( Cxvz_coord, scaling.L, model.Nx+1 );

    Cxviz = DoodzMalloc( sizeof(float)*nxviz);
    DoubleToFloat( xviz, Cxviz, nxviz );
    ScaleBack( Cxviz, scaling.L, nxviz );

    Czviz = DoodzMalloc( sizeof(float)*nzviz);
    DoubleToFloat( zviz, Czviz, nzviz );
    ScaleBack( Czviz, scaling.L, nzviz );

    Cxviz_hr = DoodzMalloc( sizeof(float)*nxviz_hr);
    DoubleToFloat( xviz_hr, Cxviz_hr, nxviz_hr );
    ScaleBack( Cxviz_hr, scaling.L, nxviz_hr );

    Czviz_hr = DoodzMalloc( sizeof(float)*nzviz_hr);
    DoubleToFloat( zviz_hr, Czviz_hr, nzviz_hr );
    ScaleBack( Czviz_hr, scaling.L, nzviz_hr );

    // Total strain
    strain  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles, particles->strain,     mesh, strain,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain, Cstrain, (model.Nx-1)*(model.Nz-1) );

    // Elastic strain
    strain_el  = DoodzCalloc((model.Nx-1)*(model.Nz-1),sizeof(double));
    P2Mastah( &model, *particles, particles->strain_el,     mesh, strain_el,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain_el  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_el, Cstrain_el, (model.Nx-1)*(model.Nz-1) );

    // Plastic strain
    strain_pl  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles, particles->strain_pl,     mesh, strain_pl,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain_pl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_pl, Cstrain_pl, (model.Nx-1)*(model.Nz-1) );

    // Power-law strain
    strain_pwl  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles, particles->strain_pwl,     mesh, strain_pwl,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain_pwl  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_pwl, Cstrain_pwl, (model.Nx-1)*(model.Nz-1) );

    // Exponential flow strain
    strain_exp  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles, particles->strain_exp,     mesh, strain_exp,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain_exp  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_exp, Cstrain_exp, (model.Nx-1)*(model.Nz-1) );

    // Linear flow strain
    strain_lin  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles, particles->strain_lin,     mesh, strain_lin,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain_lin  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_lin, Cstrain_lin, (model.Nx-1)*(model.Nz-1) );

    // GBS flow strain
    strain_gbs  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
    P2Mastah( &model, *particles, particles->strain_gbs,     mesh, strain_gbs,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
    Cstrain_gbs  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( strain_gbs, Cstrain_gbs, (model.Nx-1)*(model.Nz-1) );

    if ( model.fstrain == 1 ) {
        // Fxx
        Fxx  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->Fxx,     mesh, Fxx,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        CFxx  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fxx, CFxx, (model.Nx-1)*(model.Nz-1) );

        // Fxz
        Fxz  = DoodzCalloc((model.Nx-1)*(model.Nz-1),sizeof(double));
        P2Mastah( &model, *particles, particles->Fxz,     mesh, Fxz,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        CFxz  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fxz, CFxz, (model.Nx-1)*(model.Nz-1) );

        // Fzx
        Fzx  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->Fzx,     mesh, Fzx,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        CFzx  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fzx, CFzx, (model.Nx-1)*(model.Nz-1) );

        // Fzz
        Fzz  = DoodzCalloc((model.Nx-1)*(model.Nz-1),sizeof(double));
        P2Mastah( &model, *particles, particles->Fzz,     mesh, Fzz,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        CFzz = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Fzz, CFzz, (model.Nx-1)*(model.Nz-1) );
    }

    if ( model.aniso == 1 ) {

        // nx
        nx  = DoodzCalloc((model.Nx-1)*(model.Nz-1),sizeof(double));
        P2Mastah( &model, *particles, particles->nx,     mesh, nx,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        Cnx = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( nx, Cnx, (model.Nx-1)*(model.Nz-1) );

        // nz
        nz  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->nz,     mesh, nz,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        Cnz = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( nz, Cnz, (model.Nx-1)*(model.Nz-1) );

        // aniso factor on centroids
        Cani_fac = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( mesh->aniso_factor_n, Cani_fac, (model.Nx-1)*(model.Nz-1) );

    }

    if ( model.rec_T_P_x_z == 1 ) {

        // T0
        T0  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->T0,     mesh, T0,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) T0[k]   = zeroC/scaling.T;
        }
        CT0  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( T0, CT0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CT0, scaling.T, (model.Nx-1)*(model.Nz-1) );

        // P0
        P0  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->P0,     mesh, P0,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) P0[k]   = 0.0;
        }
        CP0  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( P0, CP0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CP0, scaling.S, (model.Nx-1)*(model.Nz-1) );

        // Tmax
        Tmax  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->Tmax,     mesh, Tmax,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) Tmax[k] = zeroC/scaling.T;
        }
        CTmax  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Tmax, CTmax, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CTmax, scaling.T, (model.Nx-1)*(model.Nz-1) );

        // Pmax
        Pmax  = DoodzCalloc((model.Nx-1)*(model.Nz-1),sizeof(double));
        P2Mastah( &model, *particles, particles->Pmax,     mesh, Pmax,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        for (k=0; k<(mesh->Nx-1)*(mesh->Nz-1); k++) {
            if ( mesh->BCt.type[k] == 30 ) Pmax[k] = 0.0;
        }
        CPmax  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( Pmax, CPmax, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( CPmax, scaling.S, (model.Nx-1)*(model.Nz-1) );


        // x0
        x0  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->x0,     mesh, x0,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        Cx0  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( x0, Cx0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( Cx0, scaling.L, (model.Nx-1)*(model.Nz-1) );

        // z0
        z0  = DoodzCalloc((model.Nx-1)*(model.Nz-1), sizeof(double));
        P2Mastah( &model, *particles, particles->z0,     mesh, z0,   mesh->BCp.type,  1, 0, interp, cent, model.itp_stencil);
        Cz0 = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
        DoubleToFloat( z0, Cz0, (model.Nx-1)*(model.Nz-1) );
        ScaleBack( Cz0, scaling.L, (model.Nx-1)*(model.Nz-1) );

    }

    Cfriction  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->fric_n, Cfriction, (model.Nx-1)*(model.Nz-1) );
    Ccohesion  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
    DoubleToFloat( mesh->C_n, Ccohesion, (model.Nx-1)*(model.Nz-1) );
    ScaleBack( Ccohesion, scaling.S, (model.Nx-1)*(model.Nz-1) );
    
//    // Get X from particles
//    X  = DoodzCalloc((model.Nx-1)*(model.Nz-1),sizeof(double));
//    Interp_P2C ( *particles,  particles->X, mesh, X, mesh->xg_coord, mesh->zg_coord, 1, 0 );
//    CX  = DoodzMalloc( sizeof(float)*(model.Nx-1)*(model.Nz-1));
//    DoubleToFloat( X, CX, (model.Nx-1)*(model.Nz-1) );


    //---------------------------------------------------

    // Topography
    if ( model.free_surf == 1 ) {

        Cxtopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->x, Cxtopo, topo_chain->Nb_part );
        ScaleBack( Cxtopo, scaling.L, topo_chain->Nb_part );

        Cztopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->z, Cztopo, topo_chain->Nb_part );
        ScaleBack( Cztopo, scaling.L, topo_chain->Nb_part );

        Cheight = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->height, Cheight, model.Nx );
        ScaleBack( Cheight, scaling.L, model.Nx );

        Ctopovx = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->vx, Ctopovx, model.Nx );
        ScaleBack( Ctopovx, scaling.V, model.Nx );

        Ctopovz = DoodzMalloc( sizeof(float)*(model.Nx+1));
        DoubleToFloat( topo->vz, Ctopovz, model.Nx+1 );
        ScaleBack( Ctopovz, scaling.V, model.Nx+1 );
        
        Ctopovx_mark = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->Vx, Ctopovx_mark, topo_chain->Nb_part );
        ScaleBack( Ctopovx_mark, scaling.V, topo_chain->Nb_part );
        
        Ctopovz_mark = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->Vz, Ctopovz_mark, topo_chain->Nb_part );
        ScaleBack( Ctopovz_mark, scaling.V, topo_chain->Nb_part );

    }

    // Generate file name
    asprintf( &FileName, "%s%05d%s",txtout, model.step, ".gzip.h5");
    if (model.writerSubfolder && strcmp(model.writerSubfolder, "")) {
      CreateDir(model.writerSubfolder);
      asprintf( &FileName, "%s/%s", model.writerSubfolder, FileName);
    }
    CreateOutputHDF5( FileName );

    // Add groups
    AddGroupToHDF5( FileName, "Model" );
    AddGroupToHDF5( FileName, "Vertices" );
    AddGroupToHDF5( FileName, "Centers" );
    AddGroupToHDF5( FileName, "VxNodes" );
    AddGroupToHDF5( FileName, "VzNodes" );
    AddGroupToHDF5( FileName, "Particles" );
    AddGroupToHDF5( FileName, "VizGrid" );
    AddGroupToHDF5( FileName, "Topo" );
    AddGroupToHDF5( FileName, "Flags" );
    AddGroupToHDF5( FileName, "Iterations");
    AddGroupToHDF5( FileName, "TimeSeries" );

    AddFieldToGroup( FileName, "Model", "Description"   , 'c', 500 * sizeof(char),  model.description, 1 );

    // Model Parameters
    A[0] = (double)(model.time) * scaling.t;
    A[1] = (double)(model.xmax - model.xmin) * scaling.L;
    A[2] = (double)(model.zmax - model.zmin) * scaling.L;
    A[3] = (double)model.Nx;
    A[4] = (double)model.Nz;
    A[5] = (double)model.dx * scaling.L;
    A[6] = (double)model.dz * scaling.L;
    A[7] = (double)model.dt * scaling.t;

    // Parameter array
    AddFieldToGroup( FileName, "Model", "Params"   , 'd', 8,  A, 1 );

    // Add time series
    AddFieldToGroup(  FileName, "TimeSeries", "Time_time"      , 'd', model.Nt+1,  mesh->Time_time      , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "Short_time"     , 'd', model.Nt+1,  mesh->Short_time     , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "Work_time"      , 'd', model.Nt+1,  mesh->Work_time      , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "Uthermal_time"  , 'd', model.Nt+1,  mesh->Uthermal_time  , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "Uelastic_time"  , 'd', model.Nt+1,  mesh->Uelastic_time  , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "T_mean_time"    , 'd', model.Nt+1,  mesh->T_mean_time    , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "P_mean_time"    , 'd', model.Nt+1,  mesh->P_mean_time    , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "sxxd_mean_time" , 'd', model.Nt+1,  mesh->sxxd_mean_time , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "szzd_mean_time" , 'd', model.Nt+1,  mesh->szzd_mean_time , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "sxz_mean_time"  , 'd', model.Nt+1,  mesh->sxz_mean_time  , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "Tii_mean_time"  , 'd', model.Nt+1,  mesh->Tii_mean_time  , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "exxd_mean_time" , 'd', model.Nt+1,  mesh->exxd_mean_time , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "ezzd_mean_time" , 'd', model.Nt+1,  mesh->ezzd_mean_time , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "exz_mean_time"  , 'd', model.Nt+1,  mesh->exz_mean_time  , 1 );
    AddFieldToGroup(  FileName, "TimeSeries", "Eii_mean_time"  , 'd', model.Nt+1,  mesh->Eii_mean_time  , 1 );


    // Grid coordinate arrays
    AddFieldToGroup( FileName, "Model", "xg_coord" , 'f', model.Nx,    Cxg_coord,  1 );
    AddFieldToGroup( FileName, "Model", "zg_coord" , 'f', model.Nz,    Czg_coord,  1 );
    AddFieldToGroup( FileName, "Model", "xc_coord" , 'f', model.Nx-1,  Cxc_coord,  1 );
    AddFieldToGroup( FileName, "Model", "zc_coord" , 'f', model.Nz-1,  Czc_coord,  1 );
    AddFieldToGroup( FileName, "Model", "xvz_coord", 'f', model.Nx+1,  Cxvz_coord, 1 );
    AddFieldToGroup( FileName, "Model", "zvx_coord", 'f', model.Nz+1,  Czvx_coord, 1 );

    // Visualisation grid
    AddFieldToGroup( FileName, "VizGrid", "xviz"    , 'f', nxviz, Cxviz,     1 );
    AddFieldToGroup( FileName, "VizGrid", "zviz"    , 'f', nzviz, Czviz,     1 );
    AddFieldToGroup( FileName, "VizGrid", "xviz_hr" , 'f', nxviz_hr, Cxviz_hr,  1 );
    AddFieldToGroup( FileName, "VizGrid", "zviz_hr" , 'f', nzviz_hr, Czviz_hr,  1 );
    AddFieldToGroup( FileName, "VizGrid", "compo"   , 'c', (nxviz-1)*(nzviz-1), compo,    1 );
    AddFieldToGroup( FileName, "VizGrid", "compo_hr", 'c', (nxviz_hr-1)*(nzviz_hr-1), compo_hr,    1 );

    // Add casted grid fields
    AddFieldToGroup( FileName, "Vertices", "rho_s", 'f', model.Nx*model.Nz,         Crho_s, 1 );
    AddFieldToGroup( FileName, "Centers" , "rho_n", 'f', (model.Nx-1)*(model.Nz-1), Crho_n, 1 );

    AddFieldToGroup( FileName, "Vertices", "eta_s", 'f', model.Nx*model.Nz,         Ceta_s, 1 );
    AddFieldToGroup( FileName, "Centers" , "eta_n", 'f', (model.Nx-1)*(model.Nz-1), Ceta_n, 1 );
    AddFieldToGroup( FileName, "VxNodes" , "Vx"   , 'f', model.Nx*(model.Nz+1),     CVx, 1 );
    AddFieldToGroup( FileName, "VzNodes" , "Vz"   , 'f', (model.Nx+1)*model.Nz,     CVz, 1 );
    AddFieldToGroup( FileName, "Centers" , "P"    , 'f', (model.Nx-1)*(model.Nz-1), CP, 1 );
    AddFieldToGroup( FileName, "Centers" , "sxxd" , 'f', (model.Nx-1)*(model.Nz-1), Csxxd, 1 );
    AddFieldToGroup( FileName, "Centers" , "szzd" , 'f', (model.Nx-1)*(model.Nz-1), Cszzd, 1 );

    AddFieldToGroup( FileName, "Centers" , "strain",    'f', (model.Nx-1)*(model.Nz-1), Cstrain, 1 );
    AddFieldToGroup( FileName, "Centers" , "strain_el", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_el, 1 );
    AddFieldToGroup( FileName, "Centers" , "strain_pl", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_pl, 1 );
    AddFieldToGroup( FileName, "Centers" , "strain_pwl", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_pwl, 1 );
    AddFieldToGroup( FileName, "Centers" , "strain_exp", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_exp, 1 );
    AddFieldToGroup( FileName, "Centers" , "strain_lin", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_lin, 1 );
    AddFieldToGroup( FileName, "Centers" , "strain_gbs", 'f', (model.Nx-1)*(model.Nz-1), Cstrain_gbs, 1 );
    AddFieldToGroup( FileName, "Centers" , "T",     'f', (model.Nx-1)*(model.Nz-1), CT, 1 );
    AddFieldToGroup( FileName, "Vertices", "sxz"  , 'f', model.Nx*model.Nz,         Csxz, 1 );
    AddFieldToGroup( FileName, "Centers" , "exxd" , 'f', (model.Nx-1)*(model.Nz-1), Cexxd, 1 );
    AddFieldToGroup( FileName, "Centers" , "ezzd" , 'f', (model.Nx-1)*(model.Nz-1), Cezzd, 1 );
    AddFieldToGroup( FileName, "Centers" , "divu" , 'f', (model.Nx-1)*(model.Nz-1), Cdivu, 1 );
    AddFieldToGroup( FileName, "Centers" , "divu_el" , 'f', (model.Nx-1)*(model.Nz-1), Cdivu_el, 1 );
    AddFieldToGroup( FileName, "Centers" , "divu_pl" , 'f', (model.Nx-1)*(model.Nz-1), Cdivu_pl, 1 );
    AddFieldToGroup( FileName, "Centers" , "divu_th" , 'f', (model.Nx-1)*(model.Nz-1), Cdivu_th, 1 );
    AddFieldToGroup( FileName, "Centers" , "divu_r"  , 'f', (model.Nx-1)*(model.Nz-1), Cdivu_r , 1 );
    AddFieldToGroup( FileName, "Vertices", "exz"  , 'f', model.Nx*model.Nz,         Cexz, 1 );
    AddFieldToGroup( FileName, "Centers" , "eII_el" , 'f', (model.Nx-1)*(model.Nz-1), CeII_el, 1 );
    AddFieldToGroup( FileName, "Centers" , "eII_pl" , 'f', (model.Nx-1)*(model.Nz-1), CeII_pl, 1 );
    AddFieldToGroup( FileName, "Centers" , "eII_pwl" , 'f', (model.Nx-1)*(model.Nz-1), CeII_pwl, 1 );
    AddFieldToGroup( FileName, "Centers" , "eII_exp" , 'f', (model.Nx-1)*(model.Nz-1), CeII_exp, 1 );
    AddFieldToGroup( FileName, "Centers" , "eII_lin" , 'f', (model.Nx-1)*(model.Nz-1), CeII_lin, 1 );
    AddFieldToGroup( FileName, "Centers" , "eII_gbs" , 'f', (model.Nx-1)*(model.Nz-1), CeII_gbs, 1 );
    AddFieldToGroup( FileName, "Centers" , "d" , 'f', (model.Nx-1)*(model.Nz-1), Cd, 1 );
    AddFieldToGroup( FileName, "Centers" , "X" , 'f', (model.Nx-1)*(model.Nz-1), CX, 1 );
    AddFieldToGroup( FileName, "Centers" , "OverS",'f', (model.Nx-1)*(model.Nz-1), COverS, 1 );

    if ( model.free_surf == 1 ) {
        AddFieldToGroup( FileName, "Topo", "z_grid" , 'f', (model.Nx), Cheight, 1 );
        AddFieldToGroup( FileName, "Topo", "Vx_grid" , 'f', (model.Nx), Ctopovx, 1 );
        AddFieldToGroup( FileName, "Topo", "Vz_grid" , 'f', (model.Nx+1), Ctopovz, 1 );
        AddFieldToGroup( FileName, "Topo", "x_mark" , 'f', topo_chain->Nb_part, Cxtopo, 1 );
        AddFieldToGroup( FileName, "Topo", "z_mark" , 'f', topo_chain->Nb_part, Cztopo, 1 );
        AddFieldToGroup( FileName, "Topo", "Vx_mark" , 'f', topo_chain->Nb_part, Ctopovx_mark, 1 );
        AddFieldToGroup( FileName, "Topo", "Vz_mark" , 'f', topo_chain->Nb_part, Ctopovz_mark, 1 );
        AddFieldToGroup( FileName, "Topo", "phase_mark" , 'i', topo_chain->Nb_part, topo_chain->phase, 1 );
    }

    AddFieldToGroup( FileName, "Centers" , "friction", 'f', (model.Nx-1)*(model.Nz-1), Cfriction, 1 );
    AddFieldToGroup( FileName, "Centers" , "cohesion", 'f', (model.Nx-1)*(model.Nz-1), Ccohesion, 1 );

    if (model.fstrain == 1) {
        AddFieldToGroup( FileName, "Centers" , "Fxx", 'f', (model.Nx-1)*(model.Nz-1), CFxx, 1 );
        AddFieldToGroup( FileName, "Centers" , "Fxz", 'f', (model.Nx-1)*(model.Nz-1), CFxz, 1 );
        AddFieldToGroup( FileName, "Centers" , "Fzx", 'f', (model.Nx-1)*(model.Nz-1), CFzx, 1 );
        AddFieldToGroup( FileName, "Centers" , "Fzz", 'f', (model.Nx-1)*(model.Nz-1), CFzz, 1 );
    }

    if (model.aniso == 1) {
        AddFieldToGroup( FileName, "Centers" , "nx", 'f', (model.Nx-1)*(model.Nz-1), Cnx, 1 );
        AddFieldToGroup( FileName, "Centers" , "nz", 'f', (model.Nx-1)*(model.Nz-1), Cnz, 1 );
        AddFieldToGroup( FileName, "Centers" , "ani_fac", 'f', (model.Nx-1)*(model.Nz-1), Cani_fac, 1 );
    }

    if (model.rec_T_P_x_z == 1) {
        AddFieldToGroup( FileName, "Centers" , "T0", 'f', (model.Nx-1)*(model.Nz-1), CT0, 1 );
        AddFieldToGroup( FileName, "Centers" , "P0", 'f', (model.Nx-1)*(model.Nz-1), CP0, 1 );
        AddFieldToGroup( FileName, "Centers" , "Tmax", 'f', (model.Nx-1)*(model.Nz-1), CTmax, 1 );
        AddFieldToGroup( FileName, "Centers" , "Pmax", 'f', (model.Nx-1)*(model.Nz-1), CPmax, 1 );
        AddFieldToGroup( FileName, "Centers" , "x0", 'f', (model.Nx-1)*(model.Nz-1), Cx0, 1 );
        AddFieldToGroup( FileName, "Centers" , "z0", 'f', (model.Nx-1)*(model.Nz-1), Cz0, 1 );
    }

    // Add cell flags for debugging 
    AddFieldToGroup( FileName, "Centers" , "BCp" , 'c', (model.Nx-1)*(model.Nz-1), mesh->BCp.type, 1 );
    AddFieldToGroup( FileName, "Flags", "tag_s" , 'c', (model.Nx-0)*(model.Nz-0), mesh->BCg.type,  1 );
    AddFieldToGroup( FileName, "Flags", "tag_n" , 'c', (model.Nx-1)*(model.Nz-1), mesh->BCp.type,  1 );
    AddFieldToGroup( FileName, "Flags", "tag_u" , 'c', (model.Nx-0)*(model.Nz+1), mesh->BCu.type,  1 );
    AddFieldToGroup( FileName, "Flags", "tag_v" , 'c', (model.Nx+1)*(model.Nz-0), mesh->BCv.type,  1 );

    // Add non-iteration log
    AddFieldToGroup( FileName, "Iterations", "rx_abs"          , 'd', Nmodel.nit_max+1, Nmodel.rx_abs         ,  1 );
    AddFieldToGroup( FileName, "Iterations", "rz_abs"          , 'd', Nmodel.nit_max+1, Nmodel.rz_abs         ,  1 );
    AddFieldToGroup( FileName, "Iterations", "rp_abs"          , 'd', Nmodel.nit_max+1, Nmodel.rp_abs         ,  1 );
    AddFieldToGroup( FileName, "Iterations", "LogIsNewtonStep" , 'i', Nmodel.nit_max+1, Nmodel.LogIsNewtonStep,  1 );
    AddFieldToGroup( FileName, "Iterations", "NumberSteps"     , 'i',                1, &(Nmodel.nit)         ,  1 );

    // Freedom
    free( FileName );
    //------------//
    DoodzFree( Crho_s );
    DoodzFree( Crho_n );
    DoodzFree( Ceta_s );
    DoodzFree( Ceta_n );
    DoodzFree( CVx );
    DoodzFree( CVz );
    DoodzFree( CP );
    DoodzFree( Csxxd );
    DoodzFree( Cszzd );
    DoodzFree( Csxz );
    DoodzFree( Cexxd );
    DoodzFree( Cezzd );
    DoodzFree( Cexz );
    DoodzFree( CeII_el  );
    DoodzFree( CeII_pl  );
    DoodzFree( CeII_pwl );
    DoodzFree( CeII_exp );
    DoodzFree( CeII_lin );
    DoodzFree( CeII_gbs );

    DoodzFree( Cxg_coord );
    DoodzFree( Czg_coord );
    DoodzFree( Cxc_coord );
    DoodzFree( Czc_coord );
    DoodzFree( Czvx_coord );
    DoodzFree( Cxvz_coord );
	DoodzFree( Cxviz );
	DoodzFree( Czviz );
    DoodzFree( Cxviz_hr );
    DoodzFree( Czviz_hr );

    DoodzFree( strain     );
    DoodzFree( strain_el  );
    DoodzFree( strain_pl  );
    DoodzFree( strain_pwl );
    DoodzFree( strain_exp );
    DoodzFree( strain_lin );
    DoodzFree( strain_gbs );
    DoodzFree( Cstrain     );
    DoodzFree( Cstrain_el  );
    DoodzFree( Cstrain_pl  );
    DoodzFree( Cstrain_pwl );
    DoodzFree( Cstrain_exp );
    DoodzFree( Cstrain_lin );
    DoodzFree( Cstrain_gbs );
    DoodzFree( CX );
    DoodzFree( CT );
    DoodzFree( Cd );
    DoodzFree( COverS );
    DoodzFree( Cdivu  );
    DoodzFree( Cdivu_el );
    DoodzFree( Cdivu_pl );
    DoodzFree( Cdivu_th );
    DoodzFree( Cdivu_r );

    if ( model.free_surf == 1 ) {
        DoodzFree( Cxtopo );
        DoodzFree( Cztopo );
        DoodzFree( Cheight );
        DoodzFree( Ctopovx );
        DoodzFree( Ctopovz );
        DoodzFree( Ctopovx_mark );
        DoodzFree( Ctopovz_mark );
    }

    DoodzFree( Cfriction );
    DoodzFree( Ccohesion );

    DoodzFree( compo    );
    DoodzFree( compo_hr );

    if ( model.fstrain == 1 ) {
        DoodzFree( Fxx  );
        DoodzFree( Fxz  );
        DoodzFree( Fzx  );
        DoodzFree( Fzz  );
        DoodzFree( CFxx );
        DoodzFree( CFxz );
        DoodzFree( CFzx );
        DoodzFree( CFzz );
    }

     if ( model.aniso == 1 ) {
         DoodzFree( nx  );
         DoodzFree( nz  );
         DoodzFree( Cnx );
         DoodzFree( Cnz );
         DoodzFree( Cani_fac );
     }

    if ( model.rec_T_P_x_z == 1 ) {
        DoodzFree( T0 );
        DoodzFree( P0 );
        DoodzFree( Tmax );
        DoodzFree( Pmax );
        DoodzFree( x0 );
        DoodzFree( z0 );
        DoodzFree( CT0 );
        DoodzFree( CP0 );
        DoodzFree( CTmax );
        DoodzFree( CPmax );
        DoodzFree( Cx0 );
        DoodzFree( Cz0 );
    }

    DoodzFree(P_total);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void WriteOutputHDF5Particles( grid *mesh, markers *particles, surface *topo, markers* topo_chain, surface *topo_ini, markers* topo_chain_ini, params model, char *txtout, mat_prop materials, scale scaling ) {

    char *FileName;
    double A[8];
    char  *part_ph, *part_gen;
    float *part_x, *part_z, *part_Vx, *part_Vz, *part_T, *part_P, *part_sxxd, *part_sxz;
    int d_part=1, ind=0, Nb_part_viz=particles->Nb_part, *part_index;
    float *Cxtopo, *Cztopo, *Cvxtopo, *Cvztopo, *Cheight, *Cvxsurf, *Cvzsurf;
    float *Cxtopo_ini, *Cztopo_ini, *Cvxtopo_ini, *Cvztopo_ini, *Cheight_ini, *Cvxsurf_ini, *Cvzsurf_ini;
    int k;

//    // Only save a give number of particles
//    if(particles->Nb_part<Nb_part_viz) {
        Nb_part_viz = particles->Nb_part;
//    }
    
//    // Tracer: find number of particles after filtering
//    Nb_part_viz = 0;
//    for (k=0; k<particles->Nb_part; k++) {
//        if (particles->phase[k]==1) { // filter the mantle - VALID FOR MY MODEL ONLY
//            Nb_part_viz++;
//        }
//    }
//
//    printf("Saving %d particles of crust\n",  Nb_part_viz);
    
//    d_part      = particles->Nb_part / (Nb_part_viz);
//    d_part      = particles->Nb_part / (Nb_part_viz);
    part_x      = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_z      = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_Vx     = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_Vz     = DoodzMalloc( sizeof(float)*Nb_part_viz );
    part_ph     = DoodzMalloc( sizeof(char) *Nb_part_viz );
    part_gen    = DoodzMalloc( sizeof(char) *Nb_part_viz );
    part_T      = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_P      = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_sxxd   = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_sxz    = DoodzMalloc( sizeof(float) *Nb_part_viz );
    part_index  = DoodzMalloc( sizeof(int) *Nb_part_viz );  // Tracer: allocate
    
//    // Tracer: store selected particle
//    for (k=0; k<particles->Nb_part; k++) {
//        if (particles->phase[k]==1) {
//        part_x[ind]     = (float)particles->x[k];
//        part_z[ind]     = (float)particles->z[k];
//        part_Vx[ind]    = (float)particles->Vx[k];
//        part_Vz[ind]    = (float)particles->Vz[k];
//        part_P[ind]     = (float)particles->P[k];
//        part_T[ind]     = (float)particles->T[k];
//        part_sxxd[ind]  = (float)particles->sxxd[k];
//        part_sxz[ind]   = (float)particles->sxz[k];
//        part_ph[ind]    = (char)particles->phase[k];
//        part_gen[ind]   = (char)particles->generation[k];
//        part_index[ind] = k;
//        ind ++;
//        }
//    }

    for (k=0; k<Nb_part_viz; k++) {
        part_x[k]     = (float)particles->x[ind];
        part_z[k]     = (float)particles->z[ind];
        part_Vx[k]    = (float)particles->Vx[ind];
        part_Vz[k]    = (float)particles->Vz[ind];
        part_P[k]     = (float)particles->P[ind];
        part_T[k]     = (float)particles->T[ind];
        part_sxxd[k]  = (float)particles->sxxd[ind];
        part_sxz[k]   = (float)particles->sxz[ind];
        part_ph[k]    = (char)particles->phase[ind];
        part_gen[k]   = (char)particles->generation[ind];
        part_index[k] = k;
        ind += d_part;
    }

    ScaleBack( part_x,     scaling.L, Nb_part_viz );
    ScaleBack( part_z,     scaling.L, Nb_part_viz );
    ScaleBack( part_Vx,    scaling.V, Nb_part_viz );
    ScaleBack( part_Vz,    scaling.V, Nb_part_viz );
    ScaleBack( part_T,     scaling.T, Nb_part_viz );
    ScaleBack( part_P,     scaling.S, Nb_part_viz );
    ScaleBack( part_sxxd,  scaling.S, Nb_part_viz );
    ScaleBack( part_sxz,   scaling.S, Nb_part_viz );

    // Topography
    if ( model.free_surf == 1 ) {

        // Real topo
        Cxtopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->x, Cxtopo, topo_chain->Nb_part );
        ScaleBack( Cxtopo, scaling.L, topo_chain->Nb_part );

        Cztopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->z, Cztopo, topo_chain->Nb_part );
        ScaleBack( Cztopo, scaling.L, topo_chain->Nb_part );

        Cheight = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->height, Cheight, model.Nx );
        ScaleBack( Cheight, scaling.L, model.Nx );

        Cvxtopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->Vx, Cvxtopo, topo_chain->Nb_part );
        ScaleBack( Cvxtopo, scaling.V, topo_chain->Nb_part );

        Cvztopo = DoodzMalloc( sizeof(float)*topo_chain->Nb_part);
        DoubleToFloat( topo_chain->Vz, Cvztopo, topo_chain->Nb_part );
        ScaleBack( Cvztopo, scaling.V, topo_chain->Nb_part );

        Cvxsurf = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo->vx, Cvxsurf, model.Nx );
        ScaleBack( Cvxsurf, scaling.V, model.Nx );

        Cvzsurf = DoodzMalloc( sizeof(float)*(model.Nx+1));
        DoubleToFloat( topo->vz, Cvzsurf, model.Nx+1);
        ScaleBack( Cvzsurf, scaling.V, model.Nx+1 );

        // Advected initial topo
        Cxtopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->x, Cxtopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cxtopo_ini, scaling.L, topo_chain_ini->Nb_part );

        Cztopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->z, Cztopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cztopo_ini, scaling.L, topo_chain_ini->Nb_part );

        Cheight_ini = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo_ini->height, Cheight_ini, model.Nx );
        ScaleBack( Cheight_ini, scaling.L, model.Nx );

        Cvxtopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->Vx, Cvxtopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cvxtopo_ini, scaling.V, topo_chain_ini->Nb_part );

        Cvztopo_ini = DoodzMalloc( sizeof(float)*topo_chain_ini->Nb_part);
        DoubleToFloat( topo_chain_ini->Vz, Cvztopo_ini, topo_chain_ini->Nb_part );
        ScaleBack( Cvztopo_ini, scaling.V, topo_chain_ini->Nb_part );

        Cvxsurf_ini = DoodzMalloc( sizeof(float)*(model.Nx));
        DoubleToFloat( topo_ini->vx, Cvxsurf_ini, model.Nx );
        ScaleBack( Cvxsurf_ini, scaling.V, model.Nx );

        Cvzsurf_ini = DoodzMalloc( sizeof(float)*(model.Nx+1));
        DoubleToFloat( topo_ini->vz, Cvzsurf_ini, model.Nx+1);
        ScaleBack( Cvzsurf_ini, scaling.V, model.Nx+1 );

    }

    // Generate file name
    asprintf( &FileName, "%s%05d%s",txtout, model.step, ".gzip.h5");
    CreateOutputHDF5( FileName );

    // Add groups
    AddGroupToHDF5( FileName, "Model" );
    AddGroupToHDF5( FileName, "Vertices" );
    AddGroupToHDF5( FileName, "Centers" );
    AddGroupToHDF5( FileName, "VxNodes" );
    AddGroupToHDF5( FileName, "VzNodes" );
    AddGroupToHDF5( FileName, "Particles" );
    AddGroupToHDF5( FileName, "VizGrid" );
    AddGroupToHDF5( FileName, "Topo" );
    AddGroupToHDF5( FileName, "Topo_ini" );

    // Model Parameters
    A[0] = (double)(model.time) * scaling.t;
    A[1] = (double)(model.xmax - model.xmin) * scaling.L;
    A[2] = (double)(model.zmax - model.zmin) * scaling.L;
    A[3] = (double)model.Nx;
    A[4] = (double)model.Nz;
    A[5] = (double)model.dx * scaling.L;
    A[6] = (double)model.dz * scaling.L;
    A[7] = (double)model.dt * scaling.t;

    // Parameter array
    AddFieldToGroup( FileName, "Model", "Params"   , 'd', 8,  A, 1 );

    // Add casted partial particle fields
    AddFieldToGroup( FileName, "Particles", "x"    , 'f', Nb_part_viz, part_x,     1 );
    AddFieldToGroup( FileName, "Particles", "z"    , 'f', Nb_part_viz, part_z,     1 );
    AddFieldToGroup( FileName, "Particles", "phase", 'c', Nb_part_viz, part_ph,    1 );
    
    // Tracer: write to file
    AddFieldToGroup( FileName, "Particles", "index"    , 'i', Nb_part_viz, part_index,     1 );

    AddFieldToGroup( FileName, "Particles", "generation", 'c', Nb_part_viz, part_gen,    1 );


    AddFieldToGroup( FileName, "Particles", "T", 'f', Nb_part_viz, part_T,     1 );
    AddFieldToGroup( FileName, "Particles", "P", 'f', Nb_part_viz, part_P,    1 );


    AddFieldToGroup( FileName, "Particles", "Vx", 'f', Nb_part_viz, part_Vx,    1 );
    AddFieldToGroup( FileName, "Particles", "Vz", 'f', Nb_part_viz, part_Vz,    1 );

    AddFieldToGroup( FileName, "Particles", "sxxd", 'f', Nb_part_viz, part_sxxd,    1 );
    AddFieldToGroup( FileName, "Particles", "sxz",  'f', Nb_part_viz, part_sxz,    1 );

    if ( model.free_surf == 1 ) {
        AddFieldToGroup( FileName, "Topo", "height" , 'f', (model.Nx), Cheight, 1 );
        AddFieldToGroup( FileName, "Topo", "vxsurf" , 'f', (model.Nx), Cvxsurf, 1 );
        AddFieldToGroup( FileName, "Topo", "vzsurf" , 'f', (model.Nx+1), Cvzsurf, 1 );
        AddFieldToGroup( FileName, "Topo", "x" , 'f', topo_chain->Nb_part, Cxtopo, 1 );
        AddFieldToGroup( FileName, "Topo", "z" , 'f', topo_chain->Nb_part, Cztopo, 1 );
        AddFieldToGroup( FileName, "Topo", "vx" , 'f', topo_chain->Nb_part, Cvxtopo, 1 );
        AddFieldToGroup( FileName, "Topo", "vz" , 'f', topo_chain->Nb_part, Cvztopo, 1 );
        AddFieldToGroup( FileName, "Topo", "phase" , 'i', topo_chain->Nb_part, topo_chain->phase, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "height" , 'f', (model.Nx), Cheight_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "vxsurf" , 'f', (model.Nx), Cvxsurf_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "vzsurf" , 'f', (model.Nx+1), Cvzsurf_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "x" , 'f', topo_chain_ini->Nb_part, Cxtopo_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "z" , 'f', topo_chain_ini->Nb_part, Cztopo_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "vx" , 'f', topo_chain_ini->Nb_part, Cvxtopo_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "vz" , 'f', topo_chain_ini->Nb_part, Cvztopo_ini, 1 );
        AddFieldToGroup( FileName, "Topo_ini", "phase" , 'i', topo_chain_ini->Nb_part, topo_chain_ini->phase, 1 );
    }

    // Freedom
    free( FileName );
    //------------//
    DoodzFree( part_x );
    DoodzFree( part_z );
    DoodzFree( part_Vx );
    DoodzFree( part_Vz );
    DoodzFree( part_ph );
    DoodzFree( part_gen );
    DoodzFree( part_T );
    DoodzFree( part_P );
    DoodzFree( part_sxxd );
    DoodzFree( part_sxz  );
    DoodzFree( part_index  ); // Tracer: clean
    //------------//

    if ( model.free_surf == 1 ) {
        DoodzFree( Cxtopo );
        DoodzFree( Cztopo );
        DoodzFree( Cvxtopo );
        DoodzFree( Cvxsurf );
        DoodzFree( Cvzsurf );
        DoodzFree( Cvztopo );
        DoodzFree( Cheight );

        DoodzFree( Cxtopo_ini );
        DoodzFree( Cztopo_ini );
        DoodzFree( Cvxtopo_ini );
        DoodzFree( Cvxsurf_ini );
        DoodzFree( Cvzsurf_ini );
        DoodzFree( Cvztopo_ini );
        DoodzFree( Cheight_ini );
    }
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void WriteResiduals( grid Mmesh, params model, Nparams Nmodel, scale scaling ) {

    char *FileName;
    double A[8];

    asprintf( &FileName, "Residuals%05d%s", Nmodel.nit, ".gzip.h5");
    CreateOutputHDF5( FileName );


    // Model Parameters
    A[0] = (double)(model.time) * scaling.t;
    A[1] = (double)(model.xmax - model.xmin) * scaling.L;
    A[2] = (double)(model.zmax - model.zmin) * scaling.L;
    A[3] = (double)model.Nx;
    A[4] = (double)model.Nz;
    A[5] = (double)model.dx * scaling.L;
    A[6] = (double)model.dz * scaling.L;
    A[7] = (double)model.dt * scaling.t;

    // Add groups
    AddGroupToHDF5( FileName, "Model" );
    AddGroupToHDF5( FileName, "Vertices" );
    AddGroupToHDF5( FileName, "Centers" );
    AddGroupToHDF5( FileName, "VxNodes" );
    AddGroupToHDF5( FileName, "VzNodes" );
    AddGroupToHDF5( FileName, "Particles" );
    AddGroupToHDF5( FileName, "VizGrid" );
    AddGroupToHDF5( FileName, "Topo" );

    // Scaling
    ArrayTimesScalar( Mmesh.eta_phys_n, scaling.eta, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.rho_n,  scaling.rho, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.eta_phys_s, scaling.eta, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rho_s,  scaling.rho, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rp,      scaling.E,   (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.ru,      scaling.F,   (model.Nx)*(model.Nz+1) );
    ArrayTimesScalar( Mmesh.rv,      scaling.F,   (model.Nx+1)*(model.Nz) );

    // Parameter array
    AddFieldToGroup( FileName, "Model", "Params"   , 'd', 8,  A, 1 );

    AddFieldToGroup( FileName, "VxNodes" , "ru"  , 'd', model.Nx*(model.Nz+1),     Mmesh.ru, 1 );
    AddFieldToGroup( FileName, "VzNodes" , "rv"  , 'd', (model.Nx+1)*model.Nz,     Mmesh.rv, 1 );
    AddFieldToGroup( FileName, "Centers" , "rp"  , 'd', (model.Nx-1)*(model.Nz-1), Mmesh.rp, 1 );
    AddFieldToGroup( FileName, "Centers" , "rho"  , 'd', (model.Nx-1)*(model.Nz-1), Mmesh.rho_n, 1 );
    AddFieldToGroup( FileName, "Vertices" , "rho"  , 'd', (model.Nx)*(model.Nz), Mmesh.rho_s, 1 );
    AddFieldToGroup( FileName, "Centers" , "eta"  , 'd', (model.Nx-1)*(model.Nz-1), Mmesh.eta_phys_n, 1 );
    AddFieldToGroup( FileName, "Vertices" , "eta"  , 'd', (model.Nx)*(model.Nz), Mmesh.eta_phys_s, 1 );

    // Scaling
    ArrayTimesScalar( Mmesh.eta_phys_n, 1.0/scaling.eta, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.rho_n,  1.0/scaling.rho, (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.eta_phys_s, 1.0/scaling.eta, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rho_s,  1.0/scaling.rho, (model.Nx)*(model.Nz) );
    ArrayTimesScalar( Mmesh.rp,      1.0/scaling.E,   (model.Nx-1)*(model.Nz-1) );
    ArrayTimesScalar( Mmesh.ru,      1.0/scaling.F,   (model.Nx)*(model.Nz+1) );
    ArrayTimesScalar( Mmesh.rv,      1.0/scaling.F,   (model.Nx+1)*(model.Nz) );

    free( FileName );

}
