#include "alloc.h"
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include "descriptors.h"
#include "exact_Pearson.h"
#include "../utils/primitive_print_funcs.h"



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
    stridedbuf<double> ts(n);   
 
    if(!ts.isvalid()){
       printf("failed to allocate memory for time series\n");
    }
    double* t = ts(0);
    if(t == nullptr){
       printf("null pointer returned\n");
       exit(1);
    }  
    
    for(int i = 0; i < n; i++){
       fscanf(f,"%lf\n",t+i);
    }
    fclose(f);
    clock_t t1 = clock();
    stridedbuf<double> mp(ts.len-m+1);
    stridedbuf<int>mpi(ts.len-m+1);
    maxpearson_partialauto<double,int>(ts,mp,mpi,m,m);
    printf("check2\n");
    writeDoubles("testoutputs/mp",mp.dat,n-m+1);
    writeInts("testoutputs/mpi",mpi.dat,n-m+1);
    clock_t t2 = clock();
    printf("time: %lf\n",static_cast<double>((t2-t1))/CLOCKS_PER_SEC);
    return 0;
}

