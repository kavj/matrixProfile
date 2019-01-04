#include<iostream>
#include<algorithm>
#include<ctime>
#include<omp.h>
#include "alloc.h"
#include "primitive_print_funcs.h"
#include "pearson.h"

int main(int argc, char* argv[]){
    if(argc < 4){
        printf("check input arguments\n");
        exit(0);
    }
    long long len = atoll(argv[2]);
    FILE* f = fopen(argv[1],"r");
    long long sublen = atoll(argv[3]);
    
    if(f == NULL){
        perror("fopen");
        exit(1);
    }
    int mlen = len - sublen + 1;
    double* ts = static_cast<double*>(alloc_aligned_buffer(len * sizeof(double)));

    if(ts == nullptr){
       std::cerr << "failed to allocate memory" << std::endl;
       exit(1);
    }
    for(int i = 0; i < len; i++){
       int throwaway = fscanf(f, "%lf\n", &ts[i]);
    }
    fclose(f);
    #if defined(_OPENMP)
    double t1 = omp_get_wtime(); 
    #else
    clock_t t1 = clock();
    #endif
    double* mp = static_cast<double*>(alloc_aligned_buffer(mlen * sizeof(double)));
    long long* mpi = static_cast<long long*>(alloc_aligned_buffer(mlen * sizeof(long long)));
    std::fill(mp, mp + mlen, -1.0);
    std::fill(mpi, mpi + mlen, -1); 
    int e = nautocorr_reduc(ts, mp, mpi, sublen, sublen);
    pearson2zned(mp, mlen, sublen);    
    #if defined(_OPENMP)
    double t2 = omp_get_wtime();
    std::cout << "time: " << t2 - t1 << std::endl;
    #else
    clock_t t2 = clock();
    std::cout << "time: " << static_cast<double>(t2 - t1)/CLOCKS_PER_SEC << std::endl;
    #endif
    if(e != errs::none){
       printf("miscellaneous error (this is a debugging file anyway)\n");
    }
    writeDoubles("cppoutput/mp", mp, mlen);
    writeLongs("cppoutput/mpi", mpi, mlen);
    return 0;
}
