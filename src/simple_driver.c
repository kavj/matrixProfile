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
    int len = atoi(argv[2]);
    FILE* f;
    fopen_s(&f, argv[1], "r");
    int sublen = atoi(argv[3]);
    if(f == NULL){
        perror("fopen_s");
        exit(1);
    }
    int mlen = len - sublen + 1;
    double* ts = alloc_buff(paddedlen(len));   
    if(ts == NULL){
       printf_s("failed to allocate memory\n");
    }
    for(int i = 0; i < len; i++){
       fscanf_s(f, "%lf\n", &ts[i]);
    }
    fclose(f);
    clock_t t1 = clock();
    mprofile p = pearson_pauto_reduc(ts, len, sublen, sublen);
    clock_t t2 = clock();
    //writeDoubles("coutput/mp", p.mp, p.len);
    //writeLongs("coutput/mpi", p.mpi, p.len);
    dealloc_buff(ts);
    printf_s("time: %lf\n", (double)((t2 - t1))/CLOCKS_PER_SEC);
    return 0;
}
