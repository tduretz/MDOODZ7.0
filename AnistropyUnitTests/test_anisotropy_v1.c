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
   double *Txx, *Tzz, *Txz, *Tii;
} grid;

// Materials
typedef struct {
   double eta_v;
   double ani_fac_v;
   double n;
} mat_params;

// Tensor 2D
typedef struct {
   double xx, zz, yy, xz, ii;
} Tensor2D;

// Vector 2D
typedef struct {
   double x, z;
} Vector2D;

/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

void MinMaxArrayVal( double* array, int size, double* min, double* max ) { // ALREADY IN MISCFUNCTIONS.c

    // Search min/max values of an array and prints it to standard out
    *min=array[0];
    *max=array[0];
    int k;

#pragma omp parallel
    {
        double pmin=array[0], pmax=array[0];
#pragma omp for
        for (k=0; k<size; k++) {
            if (array[k]>pmax) pmax = array[k];
            if (array[k]<pmin) pmin = array[k];
        }
#pragma omp flush (max)
        if (pmax>*max) {
#pragma omp critical
            {
                if (pmax>*max) *max = pmax;
            }
        }

#pragma omp flush (min)
        if (pmin<*min) {
#pragma omp critical
            {
                if (pmin<*min) *min = pmin;
            }
        }
    }

}

void Invii( Tensor2D *T ) {
    T->ii = sqrt( 0.5*( pow(T->xx,2) + pow(T->zz,2)) + pow(T->xz,2) );
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

// Input: material params, strain rate components
// Outut: stress tensor components 
void StressEffectiveViscosity_TensorInterface( Tensor2D* Tau, Tensor2D* Eps, mat_params* materials, double nx2, double nxnz ) {
    // Initialise
    Tau->xx = 0.; Tau->zz = 0.; Tau->xz = 0.;
    const double eta_v     = materials->eta_v;
    const double n         = materials->n;
    const double ani_fac_v = materials->ani_fac_v;
    Tensor2D Eps_rot, Tau_rot;
    // Transform strain rate: Q*E*Qt
    const double nz2 = 1.0-nx2;
    Eps_rot.xx =   nx2*Eps->xx +  nz2*Eps->zz +   2.*nxnz*Eps->xz;
    Eps_rot.zz =   nz2*Eps->xx +  nx2*Eps->zz -   2.*nxnz*Eps->xz;
    Eps_rot.xz = -nxnz*Eps->xx + nxnz*Eps->zz + (nx2-nz2)*Eps->xz;
    // Compute stress in principal plane
    Tau_rot.xx = 2.*eta_v* Eps_rot.xx;
    Tau_rot.zz = 2.*eta_v* Eps_rot.zz;
    Tau_rot.xz = 2.*eta_v* Eps_rot.xz / ani_fac_v;
    // Backtransform stress: Qt*E*Q
    Tau->xx =   nx2*Tau_rot.xx +  nz2*Tau_rot.zz -   2.*nxnz*Tau_rot.xz;
    Tau->zz =   nz2*Tau_rot.xx +  nx2*Tau_rot.zz +   2.*nxnz*Tau_rot.xz;
    Tau->xz =  nxnz*Tau_rot.xx - nxnz*Tau_rot.zz + (nx2-nz2)*Tau_rot.xz;
    Invii( Tau );
}


/*--------------------------------------------------------------------------------------------------------------------*/
/*------------------------------------------------------ M-Doodz -----------------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------------------*/

int main( int argc, char *argv[] ) {

    // Kinematics
    const double Exx = -0.0;
    const double Ezz = -Exx;
    const double Exz = -1.0;

    // Materials
    mat_params materials;
    materials.eta_v     = 0.5;
    materials.ani_fac_v = 4.0;
    materials.n         = 10.;

    FILE        *GNUplotPipe;

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
    mesh.Tii = calloc( Nx*Nz, sizeof(double) );
    printf("Nx = %03d, Nz = %03d\n", Nx, Nz);

    // Initialise strain rates (input for rheology computation)
    for (int k=0; k<Nx*Nz; k++) {
        mesh.Exx[k] = Exx;
        mesh.Ezz[k] = Ezz;
        mesh.Exz[k] = Exz;
    }

    Tensor2D Tau, Eps;

    int nt = 10, n_warmup = 11;
    int nang = 11; 
    clock_t t_omp;
    double  t_diff, t_mean;
    double dir_ang = 0.0, max_ang = M_PI/2., dang = max_ang/(nang-1), a_dir, b_dir, nx2, nxnz; 
    double Txx_vec[nang][2], Tzz_vec[nang][2], Txz_vec[nang][2], Tii_vec[nang][2];
    double min_val, max_val;

    // ------------------------------------- Variable angle ------------------------------------- //

    for (int iang=0; iang<nang; iang++) {

        dir_ang = iang*dang;
        nx2     = pow( cos(dir_ang), 2);
        nxnz    = cos(dir_ang)*sin(dir_ang);
        printf("Angle: %2.3f\n", dir_ang*180./M_PI);

    #pragma omp parallel 
        { 
            if (omp_get_thread_num()==0) printf("With tensor interface on %02d threads\n", omp_get_num_threads() ); 
        }    

        // Compute rheology at all points
    #pragma omp parallel for shared( mesh ) private( Tau, Eps ) firstprivate( materials )
        for (int k=0; k<Nx*Nz; k++) {

            Eps.xx = mesh.Exx[k];
            Eps.zz = mesh.Ezz[k];
            Eps.xz = mesh.Exz[k];
            Invii(&Eps);
            StressEffectiveViscosity_TensorInterface( &Tau, &Eps, &materials, nx2, nxnz);
            mesh.Txx[k] = Tau.xx;
            mesh.Tzz[k] = Tau.zz;
            mesh.Txz[k] = Tau.xz;
            mesh.Tii[k] = Tau.ii;
        }
        MinMaxArrayVal( mesh.Txx, Nx*Nz, &min_val, &max_val ); 
        Txx_vec[iang][1] = min_val;
        Txx_vec[iang][2] = max_val;
        printf("min Txx = %2.2f --- max Txx = %2.2f\n", Txx_vec[iang][1], Txx_vec[iang][2]);
        MinMaxArrayVal( mesh.Tzz, Nx*Nz, &min_val, &max_val ); 
        Tzz_vec[iang][1] = min_val;
        Tzz_vec[iang][2] = max_val;
        printf("min Tzz = %2.2f --- max Tzz = %2.2f\n", Tzz_vec[iang][1], Tzz_vec[iang][2]);
        MinMaxArrayVal( mesh.Txz, Nx*Nz, &min_val, &max_val ); 
        Txz_vec[iang][1] = min_val;
        Txz_vec[iang][2] = max_val;
        printf("min Txz = %2.2f --- max Txz = %2.2f\n", Txz_vec[iang][1], Txz_vec[iang][2]);
        MinMaxArrayVal( mesh.Tii, Nx*Nz, &min_val, &max_val ); 
        Tii_vec[iang][1] = min_val;
        Tii_vec[iang][2] = max_val;
        printf("min Tii = %2.2f --- max Tii = %2.2f\n", Tii_vec[iang][1], Tii_vec[iang][2]);
    }

    //-------------------

    FILE* file=fopen("data4plot.txt", "wt");
            fprintf(file, "%s %s %s %s %s\n", "angle", "Txx", "Tzz", "Txz", "Tii"); //Write the data to a temporary file
        for (int iang=0; iang<nang; iang++) {
        fprintf(file, "%lf %lf %lf %lf %lf\n", (double)iang*dang*180./M_PI, Txx_vec[iang][1], Tzz_vec[iang][1], Txz_vec[iang][1], Tii_vec[iang][1]); //Write the data to a temporary file
    }
    fclose(file);

    system("gnuplot -p plot.gp");

    //----------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------//
    //----------------------------------------------------------------------------------//
    free( mesh.Exx );
    free( mesh.Ezz );
    free( mesh.Exz );
    free( mesh.Txx );
    free( mesh.Tzz );
    free( mesh.Txz );
    free( mesh.Tii );
}