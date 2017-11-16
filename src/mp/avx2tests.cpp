#include "avx2arith.hpp"
#include<cstdio>
#include<stdlib.h>

// It would be nice to cast everything to __m256d and use arithmetic operators where possible. 
// Unfortunately certain compilers (particularly visual C++) do not implement arithmetic operators for vector types
// Given that we're after fairly strict logic and that we avoid heavy reliance on cross lane permutations and masked loads and stores
// it should be sufficient to just include the appropriate implementation and data type definitions at build time, using a generic simd data type
// We're using a separate scalar kernel to avoid cases where the compiler might inappropriately use branch instructions over conditional moves 
// (which have been notably faster here).

// loads and stores

static void printDArray(double *a){
   printf("%lf %lf %lf %lf\n",a[0],a[1],a[2],a[3]);
}

static void printIntArray(long *a){
   printf("%d %d %d %d\n",a[0],a[1],a[2],a[3]);
}

static void printm256I(__m256i a){
   long b[4];
   _mm256_store_si256((__m256i*)b,a);
   printf("%d %d %d %d\n",b[0],b[1],b[2],b[3]);
}

static void printm256D(__m256d a){
   double b[4];
   _mm256_store_pd(b,a);
   printf("%lf %lf %lf %lf\n",b[0],b[1],b[2],b[3]);
}

static void printboolarray(double * a, double* b, double* mask){
   for(int i = 0; i < 4; i++){
      if(mask[i])
         printf("%lf ",a[i]);
      else
         printf("%lf ",b[i]);
   }
   printf("\n");
}



static void makeshiftmultadd(double* a, double* b, double* c){
    for(int i = 0; i < 4; i++){
       c[i] = a[i]*b[i] + a[i];
    }
}

static void makeshiftmultsub(double* a, double* b, double* c){
    for(int i = 0; i < 4; i++){
       c[i] = a[i] - a[i]*b[i];
    } 
}

static void makeshiftmult(double* a, double* b, double* c){
   for(int i = 0; i < 4; i++){
      c[i] = a[i]*b[i];
   }
}

static void makeshiftadd(double* a, double* b, double* c){
   for(int i = 0; i < 4; i++){
      c[i] = a[i]+b[i];
   }
}

static void makeshiftcompare(double* a, double* b, double* c){
   for(int i = 0; i < 4; i++){
      c[i] = a[i] > b[i];
   }
}

static void expect(void){
   printf("expected output\n");
}

static void actual(void){
   printf("actual output\n");
}




using namespace vmth;

int main(void){
   double a[4] = {1.0,5.2,3.0,7.0};
   double b[4] = {1.2,2.0,10.0,11.0};
   double f[4];
   long h[4] = {10,22,43,54};
   long z[4] = {3, 4, 6, 7};
   __m256d c = loada(a,0);
   __m256d d = loadu(b,0);
   printf("\n\ntesting aligned and unaligned load instructions\n");
   expect();
   printDArray(a);
   actual();
   printm256D(c);
   expect();
   printDArray(b);
   actual(); 
   printm256D(d);

   printf("\n\ntesting aligned and unaligned store instructions\n");  // this should probably test an acutal unaligned store as well, not just the instruction
   expect();
   printDArray(a);
   actual();
   storea(c,&a[0],0);
   printDArray(a);
   expect();
   printDArray(b);
   actual();
   storeu(d,&b[0],0);
   printDArray(b);

   printf("\n\ntesting broadcasts\n");
   __m256d g = bcast(a,1);
   __m256i j = bcast(h[2]);

   printf("\nexpecting %lf\n",a[1]);
   actual();
   printm256D(g);

   printf("\nexpecting %lu\n",h[2]);
   actual();
   printm256I(j);

   printf("\n\ntesting multiply add and multiply subtract\n");
   double buf[4];
   expect();
   makeshiftmultadd(a,b,buf); // naive multiply add
   printDArray(buf);
   actual();
   __m256d k = mult_add(c,d,c);
   printm256D(k);

   expect();
   makeshiftmultsub(a,b,buf); // and naive multiply subtract
   printDArray(buf);
   actual();
   __m256d m = mult_sub(c,d,c);
   printm256D(m);

   printf("\n\ntesting standalone packed multiply\n");
   expect();
   makeshiftmult(a,b,buf);
   printDArray(buf);
   actual();
   k = mult(c,d);
   printm256D(k);

   printf("\n\n testing standalone packed addition\n");
   expect();
   makeshiftadd(a,b,buf);
   printDArray(buf);
   actual();
   k = add(c,d);
   printm256D(k);

   printf("\n\ntesting comparison instructions\n");
   expect();
   makeshiftcompare(a,b,buf);
   printf("this one is slightly different,in the simd version, an integer output of -1 indicates that all bits are set over that operand\n"); 
   printDArray(buf);
   actual();
   k = cmpgtr(c,d);
   __m256i p = _mm256_castpd_si256(k);
   printm256I(p); // this is a mask rather than a value. Using a floating point interpretation, it won't produce a reasonable readable output.

   printf("\n\ntesting blends\n");

   printf("testing the results of the prior comparison\n");
   expect();
   printboolarray(a,b,buf);
   actual();
   __m256d q = blend(c,d,k);
   printm256D(q);
   printf("testing the use of double precision floating point comparators on 64 bit integer (long) vectors\n"); 
   expect();
   long ibuff[4];
   for(int i = 0; i < 4; i++){
       if(buf[i]){
          ibuff[i] = h[i];
       }
       else{
          ibuff[i] = z[i];
       }
   }
   printIntArray(ibuff);
   actual(); 
   p = loada(h,0);
   j = loada(z,0);
   __m256i r = blend(p,j,k);
   printm256I(r); 

   printf("\n\narray A: ");
   printDArray(a);
   printf("\n\narray B: ");
   printDArray(b);
   printf("\n\nselect a 0 only\n");
   expect();
   double blbuf1[4] = {a[0],b[1],b[2],b[3]};
   printDArray(blbuf1);
   actual();
   k = select1(c,d);
   printm256D(k);
   printf("\n\nselect a 0 and 1\n");
   expect();
   double blbuf2[4] = {a[0],a[1],b[2],b[3]};
   printDArray(blbuf2);
   actual();
   k = select12(c,d);
   printm256D(k);
   printf("\n\nselect a 0 through 2\n");
   expect();
   double blbuf3[4] = {a[0],a[1],a[2],b[3]};
   printDArray(blbuf3);
   actual();
   k = select123(c,d);
   printm256D(k);
   double f1[4] = {1.0,2.0,3.0,4.0};
   double f2[4] = {5.0,6.0,7.0,8.0};
   c = loada(f1,0);
   d = loada(f2,0);
   printf("testing register shift merge\n");
   printf("original 1\n");
   printm256D(c);
   printf("original 2\n");
   printm256D(d);
   expect();
   printf("%d %d %d %d\n",c[1],c[2],c[3],d[0]);
   actual();
   d = shiftmerge1(c,d);
   printm256D(d);
   h = shiftmerge1(h,z);
   printm256I(h);
   //printf("array H: ");
   //printArrayI(h);
   //printf("array Z: ");
   //printArrayI(z);
      
}

