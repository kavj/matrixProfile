#include<cstdio>
#include<ctime>
#include "utils/descriptors.h"
#include "utils/primitive_print_funcs.h"
#include "solvers/pearson.h"


int main(int argc, char* argv[]){
    if(argc < 4){
        printf("check input arguments\n");
        exit(0);
    }
    int n = atoi(argv[2]);
    FILE* f = fopen(argv[1],"r");
    int m = atoi(argv[3]);
    
    if(f == NULL){
        perror("fopen");
        exit(1);
    }
    int mlen = n-m+1;
    primbuf<double> ts(n);   
 
    if(!ts.valid()){
       printf("failed to allocate memory for time series\n");
    }
    double* t = ts(0);
    if(t == nullptr){
       printf("null pointer returned\n");
       exit(1);
    }  
    for(int i = 0; i < n; i++){
       fscanf(f, "%lf\n", t+i);
    }
    fclose(f);
    clock_t t1 = clock();
    primbuf<double> mp(ts.len-m+1);
    primbuf<long long>mpi(ts.len-m+1);

    pearson_pauto_reduc(ts, mp, mpi, m, m);
    writeDoubles("mp", mp.dat, n - m + 1);
    writeLongs("mpi", mpi.dat, n - m + 1);
    clock_t t2 = clock();
    printf("time: %lf\n", static_cast<double>((t2 - t1))/CLOCKS_PER_SEC);
    return 0;
}
