#define _POSIX_C_SOURCE 200809L
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include<omp.h>
#include "../utils/xprec_math.hpp"
#include "../utils/repack.hpp"


//This file is a somewhat ugly leftover hack that I will replace later

/* Use wc -l <filename> to get the number of lines in a single column csv. Pass it as argument 2 here.*/

extern void  accumTest4_7_10(double* cx,double* dx, double* df, double*s, double* a, double* output, long* outputi,int n, int m);



/*static inline void pack4s(double *src, double *dst, int i){
   __m256d aa = aload(src,i);
   __m256d ab = preshuffle(aa,aload(src,i+4));
   __m256d ac = shift2(aa,ab);
           ab = shift1(aa,ab);
   __m256d ad = uload(src,i+3);
   astore(aa,dst,4*i);
   astore(ab,dst,4*(i+1));
   astore(ac,dst,4*(i+2));
   astore(ad,dst,4*(i+3));
}
*/


void initDXDF(double* x,double* mu, double* dF, double* dX,int n, int m){    
    for(int i = 0; i < n-m; i++){
        dX[i] = (1/2)*x[i+m]-x[i];
        dF[i] = (x[i+m]-mu[i+1]) + (x[i]-mu[i]);
    }
}


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

void writeInts(const char* name, int* t, int n){
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
    posix_memalign((void**)&T,64,n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf\n",&T[i]);
    }    
    printf("check1\n");
    fclose(f);
    double* a = NULL;
    double* mu = NULL;
    double* buffer = NULL;
    double* sigmaInv = NULL;
    double* dX = NULL;
    double* dF = NULL;
    long*   bufferI = NULL;
    posix_memalign((void**)&a,64,n*sizeof(double));
    posix_memalign((void**)&mu,64,n*sizeof(double));
    posix_memalign((void**)&buffer,64,n*sizeof(double));
    posix_memalign((void**)&sigmaInv,64,n*sizeof(double));
    posix_memalign((void**)&dX,64,n*sizeof(double));
    posix_memalign((void**)&dF,64,n*sizeof(double));
    posix_memalign((void**)&bufferI,64,n*sizeof(long));

    xmean_windowed(T,mu,n,m); 
    xsInv(T, mu, sigmaInv, n, m);   
    initDXDF((double*)T, (double*)mu, (double*)dF, (double*)dX, n, m);
#ifdef __AVX2__ 
    double *ak = NULL;
    double *muk= NULL;
    double *bufferk = NULL;
    double *sigmaInvk = NULL;
    double *dFk = NULL;
    double *dXk = NULL;   
 
    posix_memalign((void**)&ak,64,4*n*sizeof(double));
    posix_memalign((void**)&muk,64,4*n*sizeof(double));
    posix_memalign((void**)&bufferk,64,4*n*sizeof(double));
    posix_memalign((void**)&sigmaInvk,64,4*n*sizeof(double));
    posix_memalign((void**)&dXk,64,4*n*sizeof(double));
    posix_memalign((void**)&dFk,64,4*n*sizeof(double));
 
    unfold(a,ak,n);
    unfold(mu,muk,n); 
    unfold(buffer,bufferk,n);
    unfold(sigmaInv,sigmaInvk,n);
    unfold(dX,dXk,n);
    unfold(dF,dFk,n);
    printf("check2\n");
    clock_t t1 = clock();
    accumTest4_7_10(T,dXk,dFk,sigmaInvk,ak,bufferk,bufferI,n,256);
    clock_t t2 = clock();
#else

#endif

    printf("%lf %lf %lf %lf\n",T[n-1],dX[n-m-1],dF[n-m-1],sigmaInv[n-m-1]);
    clock_t t3 = clock();
    clock_t t4 = clock();
    printf("done\n");
    clock_t t5 = clock();
    clock_t t6 = clock();
    printf("test:  %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    printf("test blocking in 4s:  %lf\n",(double)(t4-t3)/CLOCKS_PER_SEC);
    printf("test blocking in 8s:  %lf\n",(double)(t6-t5)/CLOCKS_PER_SEC);
    writeDoubles("meancheck2.csv",muk,4*(n-m+1));
    free(T);
    free(buffer);
    free(bufferI);
    free(mu);
    free(sigmaInv);
    free(dX);
    free(dF);
#ifdef __AVX2__
    free(ak);
    free(muk);
    free(bufferk);
    free(sigmaInvk);
#endif
    return 0;

}
