#include<cstdio>
#include<ctime>
#include<cstdlib>
#include "utils/descriptors.h"
#include "utils/primitive_print_funcs.h"
#include "solvers/pearson.h"
#define prefalign 64

int main(int argc, char* argv[]){
    if(argc < 4){
        printf("check input arguments\n");
        exit(0);
    }
    int len = atoi(argv[2]);
    FILE* f = fopen(argv[1],"r");
    int sublen = atoi(argv[3]);
    
    if(f == NULL){
        perror("fopen");
        exit(1);
    }
    int mlen = len - sublen + 1;
    primbuf<double> ts(len);   
 
    if(!ts.valid()){
       printf("failed to allocate memory for time series\n");
    }
    double* t = ts(0);
    if(t == nullptr){
       printf("null pointer returned\n");
       exit(1);
    }  
    for(int i = 0; i < len; i++){
       fscanf(f, "%lf\n", t + i);
    }
    fclose(f);
    clock_t t1 = clock();
    primbuf<double> mp(mlen);
    primbuf<long long>mpi(mlen);

    pearson_pauto_reduc(ts, mp, mpi, sublen, sublen);
    writeDoubles("mp", mp.dat, mlen);
    writeLongs("mpi", mpi.dat, mlen);
    clock_t t2 = clock();
    printf("time: %lf\n", static_cast<double>((t2 - t1))/CLOCKS_PER_SEC);
    return 0;
}
