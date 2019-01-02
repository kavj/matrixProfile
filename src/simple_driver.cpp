#include<cstdio>
#include<ctime>
#include<omp.h>
#include "utils/descriptors.h"
#include "utils/primitive_print_funcs.h"
#include "solvers/pearson.h"
//#define prefalign 64

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
    dbuf ts(len, 1);   
 
    if(!ts.valid()){
       printf("failed to allocate memory\n");
    }
    for(int i = 0; i < len; i++){
       fscanf(f, "%lf\n", ts(i));
    }
    fclose(f);
    #if defined(_OPENMP)
    double t1 = omp_get_wtime(); 
    #else
    clock_t t1 = clock();
    #endif
    dbuf mp(mlen, 1, -1.0);
    ibuf mpi(mlen, 1, -1); 
    int e = nautocorr_reduc(ts, mp, mpi, sublen, sublen);
    pearson2zned(mp(0), mlen, sublen);    
    #if defined(_OPENMP)
    double t2 = omp_get_wtime();
    printf("time: %lf\n", t2 - t1);
    #else
    clock_t t2 = clock();
    printf("time: %lf\n", static_cast<double>((t2 - t1))/CLOCKS_PER_SEC);
    #endif
    if(e != errs::none){
       printf("miscellaneous error (this is a debugging file anyway)\n");
    }
    writeDoubles("cppoutput/mp", mp(0), mlen);
    writeLongs("cppoutput/mpi", mpi(0), mlen);
    return 0;
}
