#define _POSIX_C_SOURCE 200112L 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>



void  accumTest3(double* Cxy, double* dX, double* dF, double*s, double* output,int n,int m){
    long q = 0;
    double a[128];
    //printf("%d  %d \n",n,m);
    for(int i = 0; i < (n-m); i +=4){
        for(int j = i; j < n-m;j+=4){
            __m256d c = _mm256_load_pd(Cxy+j);
            for(int k = 0; k < 4; k++){
                //q += 4;
                int p = 4*k;
                __m256d op1 = _mm256_broadcast_sd(dF+i+p);
                __m256d op2 = _mm256_loadu_pd(dX+i+k);
                __m256d op3 = _mm256_broadcast_sd(dF+j+k+m);
                __m256d op4 = _mm256_loadu_pd(dX+j+k+m);
                c = _mm256_fmadd_pd(op1,op2,c);
                c = _mm256_fmadd_pd(op3,op4,c);
                _mm256_store_pd(a+p,c);
                _mm256_store_pd(a+p,_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k))); 
                _mm256_store_pd(a+p,_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k+m)));
                __m256d f = _mm256_load_pd(output+p);
                __m256d g = _mm256_load_pd(a+p);
                __m256d h = _mm256_loadu_pd(a+k+m);
                __m256d d1 = _mm256_cmp_pd(g,f,13);
                __m256d d2 = _mm256_cmp_pd(g,h,13);         
                g = _mm256_blendv_pd(g,f,d1);
                h = _mm256_blendv_pd(h,f,d2);
                _mm256_store_pd(output+p,g);
                _mm256_store_pd(output+p+m,h);
            }
        
        }
    }
   // prin`tf("%d\n",q/n);
}


#define buffLen 32
#define buffSteps 8
#define buffWid 4
#define blockLen 128
void accumTest5(double* Cxy, double* dX, double* dF, double*s, double* output,int n,int m){
    double a[buffWid*buffWid];
    for(int i = 0; i < n-m; i+= buffLen){
        for(int j = 0; j < blockLen; j+= 3){
            __m256d c1 = _mm256_load_pd(Cxy+j);
            __m256d c2 = _mm256_load_pd(Cxy+j+4);
            __m256d c3 = _mm256_load_pd(Cxy);
            for(int k = 0; k < buffSteps; k++){
                __m256d op1 = _mm256_broadcast_sd(dF+i+k);
                __m256d op2 = _mm256_loadu_pd(dX+i+j+k);
                __m256d op3 = _mm256_loadu_pd(dF+i+j+k);
                __m256d op4 = _mm256_loadu_pd(dX+i+k);
                //c = _mm256_fmadd_pd(op1,op2,c);
            }
        }
    }
}


void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output,int n,int m){
    double a[32];
    //printf("%d  %d \n",n,m);
    for(int i = 0; i < (n-m); i +=32){
        for(int j = i; j < n-m;j+=4){
            __m256d c = _mm256_load_pd(Cxy+j);
            for(int k = 0; k < 32; k+=4){
                __m256d op1 = _mm256_broadcast_sd(dF+i+k);
                __m256d op2 = _mm256_loadu_pd(dX+i+k);
                __m256d op3 = _mm256_broadcast_sd(dF+j+k+m);
                __m256d op4 = _mm256_loadu_pd(dX+j+k+m);
                c = _mm256_fmadd_pd(op1,op2,c);
                c = _mm256_fmadd_pd(op3,op4,c);
                _mm256_store_pd(a+k,c);
                _mm256_store_pd(a+k,_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k))); 
                _mm256_store_pd(a+k,_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k+m)));
                __m256d f = _mm256_load_pd(output+k);
                __m256d g = _mm256_load_pd(a+k);
                __m256d h = _mm256_loadu_pd(a+k+m);
                __m256d d1 = _mm256_cmp_pd(g,f,13);
                __m256d d2 = _mm256_cmp_pd(g,h,13);         
                g = _mm256_blendv_pd(g,f,d1);
                h = _mm256_blendv_pd(h,f,d2);
                _mm256_store_pd(output+k,g);
                _mm256_store_pd(output+k+m,h);
            }
        }
    }
}


 
