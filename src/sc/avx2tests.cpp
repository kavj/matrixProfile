#include "avx2arith.h"
#include<stdio.h>
#include<stdlib.h>

// It would be nice to cast everything to __m256d and use arithmetic operators where possible. 
// Unfortunately certain compilers (particularly visual C++) do not implement arithmetic operators for vector types
// Given that we're after fairly strict logic and that we avoid heavy reliance on cross lane permutations and masked loads and stores
// it should be sufficient to just include the appropriate implementation and data type definitions at build time, using a generic simd data type
// We're using a separate scalar kernel to avoid cases where the compiler might inappropriately use branch instructions over conditional moves 
// (which have been notably faster here).

// loads and stores

void printDArray(double *a){
   printf("%lf %lf %lf %lf\n",a[0],a[1],a[2],a[3]);
}

void printIntArray(int *a){
   printf("%d %d %d %d\n",a[0],a[1],a[2],a[3]);
}

void printm256I(__m256i a){
   long b[4];
   _mm256_store_si256((__m256i*)b,a);
   printf("%d %d %d %d\n",b[0],b[1],b[2],b[3]);
}

void printm256D(__m256d a){
   double b[4];
   _mm256_store_pd(b,a);
   printf("%lf %lf %lf %lf\n",b[0],b[1],b[2],b[3]);
}

void makeshiftmultadd(double* a, double* b, double* c){
    for(int i = 0; i < 3; i++){
       c[i] = a[i]*b[i] + a[i];
    }
}

void makeshiftmultsub(double* a, double* b, double* c){
    for(int i = 0; i < 3; i++){
       c[i] = a[i] - a[i]*b[i];
    } 
}

void makeshiftmult(double* a, double* b, double* c){
   for(int i = 0; i < 3; i++){
      c[i] = a[i]*b[i];
   }
}

void makeshiftadd(double* a, double* b, double* c){
   for(int i = 0; i < 3; i++){
      c[i] = a[i]+b[i];
   }
}

void makeshiftcompare(double* a, double* b, double* c){
   for(int i = 0; i < 3; i++){
      c[i] = b[i] > a[i];
   }
}

void expect(void){
   printf("expected output\n");
}

void actual(void){
   printf("actual output\n");
}

int main(void){
   double a[4] = {1.0,5.2,3.0,7.0};
   double b[4] = {1.2,2.0,10.0,11.0};
   double f[4];
   long h[4] = {10,22,43,54};
   long z[4] = {3, 4, 6, 7};
   __m256d c = loada(a,0);
   __m256d d = loadu(b,0);
   
   expect();
   printDArray(a);
   printDArray(b);
   actual();
   printm256D(c);
   printm256D(d);
   expect();
   printDArray(a);
   printDArray(b);
   storea(c,a);
   storeu(d,b);
   actual();
   printDArray(a);
   printDArray(b);
   
   __m256d g = bcast(a[1]);
   __m256i j = bcast(h[2]);
   printf("testing broadcast");
   printf("expecting %lf\n",a[1]);
   printf("expecting %lu\n",h[2]);
   print256D(g);
   print256I(j);

   __m256d k = mult_add(c,d,c);
   __m256d m = mult_sub(c,d,c);
   double buf[4];
   makeshiftmultadd(a,b,buf);
   expect();
   printDArray(buf);
   makeshiftmultsub(a,b,buf);
   printDArray(buf);
   actual();
   print256D(k);
   print256D(m);
   expect();
   makeshiftmult(a,b,buf);
   actual();
   k = mult(c,d);
   printm256D(k);
   makeshiftadd(a,b,buf);
   k = add(c,d);
   expect();
   printDArray(buf);
   actual();
   print256D(k);
   expect();
   makeshiftcompare(a,b,buf);
   printDArray(buf);
   k = cmpgtr(c,d);
   actual();
   print256D(k);
   expect();
   printf("comparing blends (swizzles)\n");
   printf("array A: ");
   printDArray(a);
   printf("array B: ");
   printDArray(b);
   printf("select a 0 only\n");
   k = select1(c,d);
   print256D(k);
   printf("select a 0 and 1\n");
   k = select12(c,d);
   print256D(k);
   printf("select a 0 through 2\n");
   k = select123(c,d);
   print256D(k);
   //printf("array H: ");
   //printArrayI(h);
   //printf("array Z: ");
   //printArrayI(z);
   
   
}
