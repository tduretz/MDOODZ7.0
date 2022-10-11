#include <stdio.h>
#include <omp.h>

int main(int argc, char** argv){
    printf("Hello from process: %d\n", omp_get_num_threads());
    printf("Hello from process: %d\n", omp_get_thread_num());

    return 0;
}
