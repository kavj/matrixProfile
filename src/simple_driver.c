#ifndef __STDC_WANT_LIB_EXT1__ 
#define __STDC_WANT_LIB_EXT1__ 1
#endif
#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include "utils/primitive_print_funcs.h"
#include "utils/alloc.h"
#include "solvers/pearson.h"


int main(int argc, char* argv[]){
    if(argc < 4){
        printf("check input arguments\n");
        exit(0);
    }
    long long len = atoll(argv[2]);
    int sublen = atoll(argv[3]);
    FILE* f;
    fopen(&f, argv[1], "r");
    if(f == NULL){
        perror("fopen_s");
        exit(1);
    }
    double* ts = alloc_buff(paddedlen(len));   
    if(ts == NULL){
       printf("failed to allocate memory\n");
    }
    for(int i = 0; i < len; i++){
       fscanf(f, "%lf\n", &ts[i]);
    }
    fclose(f);
    clock_t t1 = clock();
    mprofile p = pearson_pauto_reduc(ts, len, sublen, sublen);
    pearson_to_normalized_euclidean(p.matrixProfile, p.len, sublen);
    clock_t t2 = clock();
    writeDoubles("coutput/mp", p.matrixProfile, p.len);
    writeLongs("coutput/mpi", p.matrixProfileIndex, p.len);
    dealloc_buff(ts);
    printf("time: %lf\n", (double)((t2 - t1))/CLOCKS_PER_SEC);
    return 0;
}
