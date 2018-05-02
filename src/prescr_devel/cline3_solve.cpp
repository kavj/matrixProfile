#define _POSIX_C_SOURCE 200809L
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include<omp.h>
#include<immintrin.h>
#include "../kernel/1NN_reference.h"
#include "prescr_Pearson.h"
#include "../arch/avx256.h"
#include "../utils/xprec_math.h"
#include "../utils/reg.h"
#include<unistd.h>
typedef double dt;
typedef double dtype;
typedef __m256d vtype;

using namespace avx256_t;
extern void  symm_pearson_reduct_kern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp);

extern void reference_pearson_reduc(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp);



//extern void bleh(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp);

//extern void bleh2(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int offsetr, int offsetc, int offsetmp);
//extern void bleh2(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp);
//template<typename dtype,typename vtype>
void batchcov(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   const int unroll = 8;
   const int simlen = sizeof(vtype)/sizeof(dtype);
   const int stride = simlen*unroll;
   const int block_count = count/stride;
   const int last = count - stride;
   //printf("last : %d\n",last);
   //exit(1);
   for(int i = 0; i < last; i+=stride){
      block<vtype> c;
      block<vtype> m;
      for(int j = 0; j < unroll; j++){
         m(j) = aload(mu,i+j*simlen);
      }
      for(int j = i; j < i+sublen; j++){
         for(int k = 0; k < unroll; k++){
            c(k+unroll) = uload(ts,j+simlen*k) - m(k+simlen*k);
         }
         vtype q = brdcst(query,j);
         for(int k = 0; k < unroll; k++){
            c(k) = mul_add(q,c(k+unroll),c(k));
         }
      }
      for(int j = 0; j < unroll; j++){
         astore(c(j),cov,j*simlen+i);
      }
   }
   //for(int j = 0; 

}



//extern void initcov(const dt* __restrict__ a, dt* __restrict__ cx, const dt* __restrict__ mu, int mindiag, int len, int sublen);

//This file is a somewhat ugly leftover hack that I will replace later

/* Use wc -l <filename> to get the number of lines in a single column csv. Pass it as argument 2 here.*/

//extern void  accumTest4_7_10(double* cx,double* dx, double* df, double*s, double* a, double* output, long*;outputi,int n, int m);

typedef double dt;
typedef long dti;
extern void solvempref(dt *a, dt *cx, dt *mu, dt *df, dt *dx, dt *s, dt *mp, dti *mpi, int len, int diagmin, int sublen);



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
    double* T = NULL;
    posix_memalign((void**)&T,64,2*n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf\n",&T[i]);
    }
    printf("check1\n");
    fclose(f);
    double* cx = NULL;
    double* mu = NULL;
    double* buffer = NULL;
    double* sigmaInv = NULL;
    double* dX = NULL;
    double* dF = NULL;
    int*   bufferI2 = NULL;
    posix_memalign((void**)&bufferI2,64,n*sizeof(long));
 
    long*   bufferI = NULL;
    posix_memalign((void**)&bufferI,64,n*sizeof(long));
    double* q = NULL;
    posix_memalign((void**)&q,64,2*m*sizeof(double));
    posix_memalign((void**)&cx,64,n*sizeof(double));
    posix_memalign((void**)&mu,64,n*sizeof(double));
    posix_memalign((void**)&buffer,64,(2*n)*sizeof(double));
    posix_memalign((void**)&sigmaInv,64,n*sizeof(double));
    posix_memalign((void**)&dX,64,n*sizeof(double));
    posix_memalign((void**)&dF,64,n*sizeof(double));

    xmean_windowed(T,mu,n,m); 
    xsInv(T, mu, sigmaInv, n, m);   
    init_dfdx(T,mu, dF, dX, m, n);
  
    printf("check2\n");
    for(int i = 0; i < n-m+1; i++){
       buffer[i] = -1.0;
    }       
    clock_t t1 = clock();
    batchcov(T,cx,T, mu, m,n-m+1,m);

    printf("starting\n");
    clock_t t2 = clock();

    for(int i = 0; i < 65281*64; i++){ //  i < 2064544; i++){
       //reference_pearson_reduc(cx,dF,dX,sigmaInv,buffer,bufferI2,m,m,m);
       ///symm_pearson_reduct_kern(cx,dF,dX,sigmaInv,buffer,bufferI2,m,m,m);
    }   
    clock_t t3 = clock();
    
    //printf("%lf %lf %lf %lf\n",T[n-1],dX[n-m-1],dF[n-m-1],sigmaInv[n-m-1]);
    clock_t t4 = clock();
    printf("done\n");
    clock_t t5 = clock();
    clock_t t6 = clock();
    printf("test:  %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    printf("test:  %lf\n",(double)(t3-t2)/CLOCKS_PER_SEC);
    //printf("test blocking in 4s:  %lf\n",(double)(t4-t3)/CLOCKS_PER_SEC);
    //printf("test blocking in 8s:  %lf\n",(double)(t6-t5)/CLOCKS_PER_SEC);
    free(T);
    free(cx);
    free(buffer);
    free(bufferI);
    free(mu);
    free(sigmaInv);
    free(dX);
    free(dF);
    free(q);
    return 0;
}
