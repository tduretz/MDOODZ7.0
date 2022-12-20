#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include "stdio.h"
#include "stdbool.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "time.h"

// #ifdef _OMP_
#include "omp.h"
// #else
// #define omp_get_thread_num()  0
// #define omp_get_num_threads() 1
// #define omp_get_wtime() clock()/CLOCKS_PER_SEC
// #endif

// Materials
typedef struct {
   double *Exx, *Ezz, *Exz;
   double *Txx, *Tzz, *Txz;
} grid;

// Materials
typedef struct {
   double eta_v;
   double ani_fac_v;
   double n;
} mat_params;

// Tensor 2D
typedef struct {
   double xx, zz, yy, xz;
} Tensor2D;

// Vector 2D
typedef struct {
   double x, z;
} Vector2D;

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Input: material params, strain rate components
// Outut: stress tensor components 
void StressEffectiveViscosity( double *Txx, double *Tzz, double *Txz, double exx, double ezz, double exz, mat_params* materials ) {
    *Txx = 0.; *Tzz = 0., *Txz = 0.;
    const double eta_v = materials->eta_v;
    const double n     = materials->n;
    // Compute isotropic stress
    *Txx = 2*eta_v* pow(exx, n);
    *Tzz = 2*eta_v* pow(ezz, n);
    *Txz = 2*eta_v* pow(exz, n);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Input: material params, strain rate components
// Outut: stress tensor components 
void StressEffectiveViscosity_TensorInterface( Tensor2D* Tau, Tensor2D* Eps, mat_params* materials ) {
    Tau->xx = 0.; Tau->zz = 0.; Tau->xz = 0.;
    // const double eta_v = materials->eta_v;
    // const double n     = materials->n;
    // // Compute isotropic stress
    // Tau->xx = 2*eta_v* pow(Eps->xx, n);
    // Tau->zz = 2*eta_v* pow(Eps->zz, n);
    // Tau->xz = 2*eta_v* pow(Eps->xz, n);
    //  const double eta_v = materials->eta_v;
    // const double n     = materials->n;
    // Compute isotropic stress
    Tau->xx = 2*materials->eta_v* pow(Eps->xx, materials->n);
    Tau->zz = 2*materials->eta_v* pow(Eps->zz, materials->n);
    Tau->xz = 2*materials->eta_v* pow(Eps->xz, materials->n);
}

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int main( int argc, char *argv[] ) {

    // Kinematics
    const double Exx = -1.0;
    const double Ezz = -Exx;
    const double Exz = 0.0;

    // Materials
    mat_params materials;
    materials.eta_v     = 0.5;
    materials.ani_fac_v = 1.0;
    materials.n         = 10.;

    // Deviatoric stress
    double Txx, Tzz, Txz;
    
    if (argc<2) {
        printf("Error: please define number of points\n"); 
        exit(1);
    }
    const int Nx = atoi(argv[1]);
    const int Nz = Nx;
    grid mesh;
    mesh.Exx = calloc( Nx*Nz, sizeof(double) );
    mesh.Ezz = calloc( Nx*Nz, sizeof(double) );
    mesh.Exz = calloc( Nx*Nz, sizeof(double) );
    mesh.Txx = calloc( Nx*Nz, sizeof(double) );
    mesh.Tzz = calloc( Nx*Nz, sizeof(double) );
    mesh.Txz = calloc( Nx*Nz, sizeof(double) );
    printf("Nx = %03d, Nz = %03d\n", Nx, Nz);

    // Initialise strain rates (input for rheology computation)
    for (int k=0; k<Nx*Nz; k++) {
        mesh.Exx[k] = Exx;
        mesh.Ezz[k] = Ezz;
        mesh.Exz[k] = Exz;
    }

    Tensor2D Tau, Eps;

    int nt = 10, n_warmup = 11;
    clock_t t_omp;
    double  t_diff, t_mean;

    // // ------------------------------------- TEST 1 ------------------------------------- //
    // printf("No tensor interface\n");
    // t_mean=0.;

    // for (int it=0; it<nt; it++) {
    //     t_omp = (double)omp_get_wtime();

    //     // Compute rheology at all points
    //     for (int k=0; k<Nx*Nz; k++) {
    //         StressEffectiveViscosity( &Txx, &Tzz, &Txz, mesh.Exx[k], mesh.Ezz[k], mesh.Exz[k], &materials);
    //         mesh.Txx[k] = Txx;
    //         mesh.Tzz[k] = Tzz;
    //         mesh.Txz[k] = Txz;
    //     }
    //     t_diff = (double)((double)omp_get_wtime() - t_omp);
    //     printf("Elapsed time: %f s\n", t_diff);
    //     t_mean += t_diff;
    // }
    // t_mean /=  nt;
    // printf("Mean Elapsed time: %f s\n", t_mean);

    // // ------------------------------------- TEST 2 ------------------------------------- //
    // printf("With tensor interface\n");        
    // t_mean=0.;

    // for (int it=0; it<nt; it++) {
    //     t_omp = (double)omp_get_wtime();

    //     // Compute rheology at all points
    //     for (int k=0; k<Nx*Nz; k++) {
    //         Eps.xx = mesh.Exx[k];
    //         Eps.zz = mesh.Ezz[k];
    //         Eps.xz = mesh.Exz[k];
    //         StressEffectiveViscosity_TensorInterface( &Tau, &Eps, &materials);
    //         mesh.Txx[k] = Tau.xx;
    //         mesh.Tzz[k] = Tau.zz;
    //         mesh.Txz[k] = Tau.xz;
    //     }
    //     t_diff = (double)((double)omp_get_wtime() - t_omp);
    //     printf("Elapsed time: %f s\n", t_diff);
    //     t_mean += t_diff;
    // }
    // t_mean /=  nt;
    // printf("Mean Elapsed time: %f s\n", t_mean);

      // ------------------------------------- TEST 1 ------------------------------------- //
#pragma omp parallel 
    { 
        if (omp_get_thread_num()==0) ("No tensor interface on %02d threads\n", omp_get_num_threads() ); 
    }
    t_mean=0.;

    for (int it=0; it<nt+n_warmup; it++) {
        t_omp = (double)omp_get_wtime();

        // Compute rheology at all points
#pragma omp parallel for private(Txx, Tzz, Txz) shared(mesh) firstprivate(materials) 
        for (int k=0; k<Nx*Nz; k++) {
            StressEffectiveViscosity( &Txx, &Tzz, &Txz, mesh.Exx[k], mesh.Ezz[k], mesh.Exz[k], &materials);
            mesh.Txx[k] = Txx;
            mesh.Tzz[k] = Tzz;
            mesh.Txz[k] = Txz;
        }
        if (it>n_warmup) {
            t_diff = (double)((double)omp_get_wtime() - t_omp);
            printf("Elapsed time: %f s\n", t_diff);
            t_mean += t_diff;
        }
    }
    t_mean /=  nt;
    printf("Mean Elapsed time: %f s\n", t_mean);

    //---------------------------------------//

    #pragma omp parallel 
    { 
        if (omp_get_thread_num()==0) ("No tensor interface on %02d threads\n", omp_get_num_threads() ); 
    }
    t_mean=0.;

    for (int it=0; it<nt+n_warmup; it++) {
        t_omp = (double)omp_get_wtime();

        // Compute rheology at all points
#pragma omp parallel for private(Txx, Tzz, Txz) shared(mesh) firstprivate(materials) 
        for (int k=0; k<Nx*Nz; k++) {
            StressEffectiveViscosity( &Txx, &Tzz, &Txz, mesh.Exx[k], mesh.Ezz[k], mesh.Exz[k], &materials);
            mesh.Txx[k] = Txx;
            mesh.Tzz[k] = Tzz;
            mesh.Txz[k] = Txz;
        }
        if (it>n_warmup) {
            t_diff = (double)((double)omp_get_wtime() - t_omp);
            printf("Elapsed time: %f s\n", t_diff);
            t_mean += t_diff;
        }
    }
    t_mean /=  nt;
    printf("Mean Elapsed time: %f s\n", t_mean);

    //---------------------------------------//

    #pragma omp parallel 
    { 
        if (omp_get_thread_num()==0) ("No tensor interface on %02d threads\n", omp_get_num_threads() ); 
    }
    t_mean=0.;

    for (int it=0; it<nt+n_warmup; it++) {
        t_omp = (double)omp_get_wtime();

        // Compute rheology at all points
#pragma omp parallel for private(Txx, Tzz, Txz) shared(mesh) firstprivate(materials) 
        for (int k=0; k<Nx*Nz; k++) {
            StressEffectiveViscosity( &Txx, &Tzz, &Txz, mesh.Exx[k], mesh.Ezz[k], mesh.Exz[k], &materials);
            mesh.Txx[k] = Txx;
            mesh.Tzz[k] = Tzz;
            mesh.Txz[k] = Txz;
        }
        if (it>n_warmup) {
            t_diff = (double)((double)omp_get_wtime() - t_omp);
            printf("Elapsed time: %f s\n", t_diff);
            t_mean += t_diff;
        }
    }
    t_mean /=  nt;
    printf("Mean Elapsed time: %f s\n", t_mean);

    //---------------------------------------//

    #pragma omp parallel 
    { 
        if (omp_get_thread_num()==0) ("No tensor interface on %02d threads\n", omp_get_num_threads() ); 
    }
    t_mean=0.;

    for (int it=0; it<nt+n_warmup; it++) {
        t_omp = (double)omp_get_wtime();

        // Compute rheology at all points
#pragma omp parallel for private(Txx, Tzz, Txz) shared(mesh) firstprivate(materials) 
        for (int k=0; k<Nx*Nz; k++) {
            StressEffectiveViscosity( &Txx, &Tzz, &Txz, mesh.Exx[k], mesh.Ezz[k], mesh.Exz[k], &materials);
            mesh.Txx[k] = Txx;
            mesh.Tzz[k] = Tzz;
            mesh.Txz[k] = Txz;
        }
        if (it>n_warmup) {
            t_diff = (double)((double)omp_get_wtime() - t_omp);
            printf("Elapsed time: %f s\n", t_diff);
            t_mean += t_diff;
        }
    }
    t_mean /=  nt;
    printf("Mean Elapsed time: %f s\n", t_mean);


    // ------------------------------------- TEST 2 ------------------------------------- //
// #pragma omp parallel 
//     { 
//         if (omp_get_thread_num()==0) printf("With tensor interface on %02d threads\n", omp_get_num_threads() ); 
//     }    
//     t_mean=0.;

//     for (int it=0; it<nt+n_warmup; it++) {
//         t_omp = (double)omp_get_wtime();
//         // Compute rheology at all points
// // #pragma omp parallel for shared( mesh ) private( Tau, Eps ) firstprivate( materials )
// #pragma omp parallel for private(Txx, Tzz, Txz) shared(mesh) firstprivate(materials) 
//         for (int k=0; k<Nx*Nz; k++) {

//             StressEffectiveViscosity( &Txx, &Tzz, &Txz, mesh.Exx[k], mesh.Ezz[k], mesh.Exz[k], &materials);
//             mesh.Txx[k] = Txx;
//             mesh.Tzz[k] = Tzz;
//             mesh.Txz[k] = Txz;

//             // Eps.xx = mesh.Exx[k];
//             // Eps.zz = mesh.Ezz[k];
//             // Eps.xz = mesh.Exz[k];
//             // StressEffectiveViscosity_TensorInterface( &Tau, &Eps, &materials);
//             // mesh.Txx[k] = Tau.xx;
//             // mesh.Tzz[k] = Tau.zz;
//             // mesh.Txz[k] = Tau.xz;
//         }
//         if (it>n_warmup) {
//             t_diff = (double)((double)omp_get_wtime() - t_omp);
//             printf("Elapsed time: %f s\n", t_diff);
//             t_mean += t_diff;
//         }
//     }
//     t_mean /=  nt;
//     printf("Mean Elapsed time: %f s\n", t_mean);


    //----------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------//
    free( mesh.Exx );
    free( mesh.Ezz );
    free( mesh.Exz );
    free( mesh.Txx );
    free( mesh.Tzz );
    free( mesh.Txz );
}