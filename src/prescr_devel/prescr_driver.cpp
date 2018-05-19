#include "alloc.h"
#include<cstdio>
#include<cstdlib>
#include<ctime>
#include "descriptors.h"
#include "auto_Pearson.h"
#include "exact_Pearson.h"



void writeDoubles(const char* name, double* t, int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%.15lf\n",t[i]);
    }
    fclose(f);
}

void writeInts(const char* name, long* t, int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%d\n",t[i]);
    }
    fclose(f);
}


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
    //const int dstride = 65536-m+1;
    stridedbuf<double> ts(1,n,n);
    
    // these should be count, stride, length so update later
    //std::fill(
    if(!ts.isvalid()){
       printf("failed to allocate memory for time series\n");
    //   exit(1);
    }
    double* t = ts(0);
    if(t == nullptr){
       printf("null pointer returned\n");
       exit(1);
    }  
//    printf("length: %d address:%lu\n",n,t);
//    printf("%lf %lf \n",ts.dat[0],ts.dat[10]);
    
    for(int i = 0; i < n/10; i++){
       fscanf(f,"%lf\n",&(t[i]));
    }
    printf("check1\n");

    fclose(f);
    clock_t t1 = clock();
    maxpearson_partialauto<double,long long>(ts,m,mlen,m);
    clock_t t2 = clock();
    printf("time: %lu\n",t2-t1);
    
    return 0;
}
