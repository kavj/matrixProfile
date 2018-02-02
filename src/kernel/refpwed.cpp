#include<cstdio>
#include<cmath>
#include <algorithm>
//#include "refpwed.h" 
#define tilesz 128 
#define largetile 512 
using namespace std;

static inline void mpkern(double* __restrict__ cx, const double*  __restrict__ dx, const double* __restrict__ xf, const double* __restrict__ s, double* __restrict__ buf, int diag, int offset) __attribute__ ((always_inline));

/*
static inline void mpkern(double* __restrict__ cx, const double* __restrict__ dx, const double* __restrict__ xf, double* __restrict__ s, double* __restrict__ mp, long* __restrict__ mpi, int diag, int offset){
   cx = (double*)__builtin_assume_aligned(cx,64);
   dx = (double*)__builtin_assume_aligned(dx,64);
   xf = (double*)__builtin_assume_aligned(xf,64);
   s  = (double*)__builtin_assume_aligned(s,64);
   mp =(double*)__builtin_assume_aligned(mp,64);
   
   const int m = diag+offset;
   for(int i = diag; i < offset+tilesz; i+=8){
      double r1 = cx[i];
      double r2 = cx[i+1];
      double r3 = cx[i+2];
      double r4 = cx[i+3];
      double r5 = cx[i+4];
      double r6 = cx[i+5];
      double r7 = cx[i+6];
      double r8 = cx[i+7];       
      for(int j = offset; j < diag+tilesz; j++){
         int k = i+j;
         double dxi = dx[i];
         r1 = fma(dxi,xf[k],r1);
         r2 = fma(dxi,xf[k+1],r2);
         r3 = fma(dxi,xf[k+2],r3);
         r4 = fma(dxi,xf[k+3],r4);
         r5 = fma(dxi,xf[k+4],r5);
         r6 = fma(dxi,xf[k+5],r6);
         r7 = fma(dxi,xf[k+6],r7);
         r8 = fma(dxi,xf[k+7],r8);

         dxi = xf[i];
         r1 = fma(dxi,dx[k],r1);
         r2 = fma(dxi,dx[k+1],r2);
         r3 = fma(dxi,dx[k+2],r3);
         r4 = fma(dxi,dx[k+3],r4);
         r5 = fma(dxi,dx[k+4],r5);
         r6 = fma(dxi,dx[k+5],r6);
         r7 = fma(dxi,dx[k+6],r7);
         r8 = fma(dxi,dx[k+7],r8);

         double t1 = r1*s[k];
         double t2 = r2*s[k+1];
         double t3 = r3*s[k+2];
         double t4 = r4*s[k+3];
         double t5 = r1*s[k+4];
         double t6 = r2*s[k+5];
         double t7 = r3*s[k+6];
         double t8 = r4*s[k+7];

         t1 = r1*s[k];
         t2 = r2*s[k+1];
         t3 = r3*s[k+2];
         t4 = r4*s[k+3];
         t5 = r1*s[k+4];
         t6 = r2*s[k+5];
         t7 = r3*s[k+6];
         t8 = r4*s[k+7];
        
         mp[k] = max(t1,mp[k]);
         mp[k+1] = max(t2,mp[k+1]);
         mp[k+2] = max(t3,mp[k+2]);
         mp[k+3] = max(t4,mp[k+3]);
         mp[k+4] = max(t5,mp[k+4]);
         mp[k+5] = max(t6,mp[k+5]);
         mp[k+6] = max(t7,mp[k+6]);
         mp[k+7] = max(t8,mp[k+7]);

         t1 = max(t1,t2);
         t3 = max(t3,t4);
         t5 = max(t5,t6);
         t7 = max(t7,t8);
         t1 = max(t1,t3);
         t5 = max(t5,t7);
         t1 = max(t1,mp[i]);
         mp[i] = max(t1,t5);
      }
   }
}*/




static inline void mpkern(double* __restrict__ cx, const double* __restrict__ dx, const double* __restrict__ xf, double* __restrict__ s, double* __restrict__ mp, long* __restrict__ mpi, int diag, int offset){
   cx = (double*)__builtin_assume_aligned(cx,64);
   dx = (double*)__builtin_assume_aligned(dx,64);
   xf = (double*)__builtin_assume_aligned(xf,64);
   s  = (double*)__builtin_assume_aligned(s,64);
   mp =(double*)__builtin_assume_aligned(mp,64);
   
   const int m = diag+offset;
   for(int i = diag; i < offset+tilesz; i++){
      double r1 = cx[i];
      double r2 = cx[i+1];
      double r3 = cx[i+2];
      double r4 = cx[i+3];
      double r5 = cx[i+4];
      double r6 = cx[i+5];
      double r7 = cx[i+6];
      double r8 = cx[i+7];       
      for(int j = offset; j < diag+tilesz; j++){
         int k = i+j;
         double dxi = dx[i];
         r1 = fma(dxi,xf[k],r1);
         r2 = fma(dxi,xf[k+1],r2);
         r3 = fma(dxi,xf[k+2],r3);
         r4 = fma(dxi,xf[k+3],r4);
         r5 = fma(dxi,xf[k+4],r5);
         r6 = fma(dxi,xf[k+5],r6);
         r7 = fma(dxi,xf[k+6],r7);
         r8 = fma(dxi,xf[k+7],r8);

         dxi = xf[i];
         r1 = fma(dxi,dx[k],r1);
         r2 = fma(dxi,dx[k+1],r2);
         r3 = fma(dxi,dx[k+2],r3);
         r4 = fma(dxi,dx[k+3],r4);
         r5 = fma(dxi,dx[k+4],r5);
         r6 = fma(dxi,dx[k+5],r6);
         r7 = fma(dxi,dx[k+6],r7);
         r8 = fma(dxi,dx[k+7],r8);

         double t1 = r1*s[k];
         double t2 = r2*s[k+1];
         double t3 = r3*s[k+2];
         double t4 = r4*s[k+3];
         double t5 = r1*s[k+4];
         double t6 = r2*s[k+5];
         double t7 = r3*s[k+6];
         double t8 = r4*s[k+7];

         t1 = r1*s[k];
         t2 = r2*s[k+1];
         t3 = r3*s[k+2];
         t4 = r4*s[k+3];
         t5 = r1*s[k+4];
         t6 = r2*s[k+5];
         t7 = r3*s[k+6];
         t8 = r4*s[k+7];
        
         if(mp[k] < t1){
            mp[k] = t1;
            mpi[i] = k;
         }
         if(mp[k+1] < t2){
            mp[k+1] = t2;
            mpi[k+1] = j;
         }
         if(mp[k+2] < t3){
            mp[k+2] = t3;
            mpi[k+2] = j;
         }
         if(mp[k+3] < t4){
            mp[k+3] = t4;
            mpi[k+3] = j;
         }
         if(mp[k+4] < t5){
            mp[k+4] = t5;
            mp[k+4] = j;
         }
         if(mp[k+5] < t6){
            mp[k+5] = t6;
            mpi[k+5] = j;
         }
         if(mp[k+6] < t7){
            mp[k+6] = t7;
            mpi[k+6] = j;
         }
         if(mp[k+7] < t8){
            mp[k+7] = t8; 
            mpi[k+7] = j;
         }
         k -= m;

         int col1 = offset;
         int col2 = offset+2;
         int col3 = offset+4;
         int col4 = offset+6;

         if(t1 < t2){
            t1 = t2;
            col1++;
         }
         if(t3 < t4){
            t3 = t4;
            col2++;
         }
         if(t5 < t6){
            t5 = t6;
            col3++;
         }
         if(t7 < t8){
            t7 = t8;
            col4++;
         }
         if(t1 < t3){
            t1 = t3;
            col1 = col2;
         }
         if(t5 < t7){
            t5 = t7;
            col3 = col4;
         }
         if(t1 < t5){
            t1 = t5;
            col1 = col3;
         }
         if(t1 < mp[offset]){
            mp[offset] = t1;
            mpi[offset] = col1;
         }
      
      }
   }
}


static inline void xcorr_kern(double* __restrict__ cx, double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double* __restrict__ mp, long* __restrict__ mpi, int diag, int offset){
   for(int subdiag = diag; subdiag < diag + largetile; subdiag += tilesz){
      for(int suboff = offset; suboff < offset+largetile; suboff+= tilesz){
         mpkern(cx, dx, xf, s, mp, mpi, diag, offset);
      }
   }
}

//template<typename dtype, typename itype>
void  accumTest4_7_10(double* cx,double* dx, double* xf, double*s, double* buf, double* output, long* outputi,int n, int m){
   unsigned long t= 0;
   #pragma omp parallel for  
   for(int diag = 0; diag < n-16*m; diag += largetile){
      for(int offset = 0; offset < n-diag-16*m; offset+=largetile){
         //mpkern(cx, dx, xf, s, output, outputi, diag, offset);
         xcorr_kern(cx,dx,xf,s,output,outputi,diag,offset);
         t+= largetile*largetile;
      }
   }   
   printf("iterations: %lu \n",t);
}


