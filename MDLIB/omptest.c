#include <omp.h>
 
#include <stdio.h>
#include <stdlib.h>
 
int main(int argc, char* argv[])
{
    printf("Num threads = %03d\n", omp_get_num_threads()); 
    // Beginning of parallel region
    #pragma omp parallel
    {

    printf("Num threads = %03d\n", omp_get_num_threads()); 
 
        printf("Hello World... from thread = %d\n",
               omp_get_thread_num());
    }
    // Ending of parallel region
}
