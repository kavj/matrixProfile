#include<stdio.h>
#include<cmath>
#include "../utils/reg.h"
// this could implement pairwise later

// It could also be optimized to do sets of 8



/*
void pears_reduct_avx(double* __restrict__ cov,  const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, int cindoffset){
   cov = (double*) __builtin_assume_aligned(cov,32);
   df = (const double*) __builtin_assume_aligned(df,32);
   dx = (const double*) __builtin_assume_aligned(dx,32);
   invn = (const double*) __builtin_assume_aligned(invn,32);
   xcorr = (double*) __builtin_assume_aligned(xcorr,32);
   cindex = (int*) __builtin_assume_aligned(cindex,32); 
   
   const int unroll = 8;
   const int count = 32;
   const int aligned = count;
   block<__m256d> cov_r;

   for(int i = 0; i < aligned; i++){
      for(int j = 0; j < unroll; j++){
         cov_r(j) = aload(cov,j*4); 
      } 
      double q0 = brdcst(dx,i);
      double q1 = brdcst(df,i);
      for(int j = 0; j < unroll; j++){
         cov_r(j) = mul_add(q0,uload(df,i+4*j+cindoffset));
         cov_r(j) = mul_add(q1,uload(dx,i+4*j+cindoffset));
      }
      q0 = brdcst(uload(invn,i));
      for(int j = 0; j < unroll; j++){
         astore(cov_r(j),cov,i+4*j);
         cov_r(j) *= q0;
      }
      for(int j = 0; j < unroll; j++){
         corr(j) = cov_r(j)*uload(invn,i+j);
      }
      block<__m256d> cmpmsk;
      for(int j = 0; j < unroll; j++){
         cmpmsk(j) = x
         if(uload(xcorr,i+4*j) < corr(j)){
             = corr(j);
            cindex[i+j] = i+cindoffset;
         }
      }
      block<int> cind;
      for(int j = 0; j < 4; j++){
         if(corr(2*j) < corr(2*j+1)){
            corr(2*j) = corr(2*j+1);
            cind(j) = 2*j+1;
         }
         else{
            cind(j) = 2*j;
         }
      }
      for(int j = 0; j < 2; j++){
         if(corr(4*j) < corr(4*j+2)){
            corr(4*j) = corr(4*j+2);
            cind(4*j) = cind(4*j+2);
         }
      }
      if(corr(0) < corr(4)){
         if(xcorr[i] < corr(4)){
            xcorr[i] = corr(4);
            cindex[i] = cind(4) + cindoffset;
         }
      }
      else if(xcorr[i+cindoffset] < corr(0)){
         xcorr[i] = corr(0);
         cindex[i] = cind(0) + i;
      }
   }
}
*/



void pears_reduct(double* __restrict__ cov,  const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, int cindoffset){
   cov = (double*) __builtin_assume_aligned(cov,32);
   df = (const double*) __builtin_assume_aligned(df,32);
   dx = (const double*) __builtin_assume_aligned(dx,32);
   invn = (const double*) __builtin_assume_aligned(invn,32);
   xcorr = (double*) __builtin_assume_aligned(xcorr,32);
   cindex = (int*) __builtin_assume_aligned(cindex,32); 
   
   const int unroll = 8;
   const int count = 32;
   const int aligned = count;
   block<double> cov_r;
   for(int i = 0; i < unroll; i++){
      cov_r(i) = cov[i];
   }
   for(int i = 0; i < aligned; i+= unroll){
      double q0 = dx[i];
      double q1 = df[i];
      for(int j = 0; j < unroll; j++){
         cov_r(j) += q0*df[i+j+cindoffset];
         cov_r(j) += q1*dx[i+j+cindoffset];
      }
      block<double> corr;
      q0 = invn[i];
      for(int j = 0; j < unroll; j++){
         corr(j) = cov_r(j)*q0;
      }
      for(int j = 0; j < unroll; j++){
         corr(j) = cov_r(j)*invn[i+j];
         cov[i+j] = corr(j);
      }
      for(int j = 0; j < unroll; j++){
         if(xcorr[i+j] < corr(j)){
            xcorr[i+j] = corr(j);
            cindex[i+j] = i+cindoffset;
         }
      }
      block<int> cind;
      for(int j = 0; j < 4; j++){
         if(corr(2*j) < corr(2*j+1)){
            corr(2*j) = corr(2*j+1);
            cind(j) = 2*j+1;
         }
         else{
            cind(j) = 2*j;
         }
      }
      for(int j = 0; j < 2; j++){
         if(corr(4*j) < corr(4*j+2)){
            corr(4*j) = corr(4*j+2);
            cind(4*j) = cind(4*j+2);
         }
      }
      if(corr(0) < corr(4)){
         if(xcorr[i] < corr(4)){
            xcorr[i] = corr(4);
            cindex[i] = cind(4) + cindoffset;
         }
      }
      else if(xcorr[i+cindoffset] < corr(0)){
         xcorr[i] = corr(0);
         cindex[i] = cind(0) + i;
      }
   }
}



   //  maxpearson_reference(double*,                  double const*,                 double const*,                 double const*,                   double*,                  int*,                  int,     int,         int)'

void maxpearson_reference(double*  __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ mp, int* __restrict__ mpi, int len, int diagmin, int sublen){
   for(int i = diagmin; i < len; i++){
      double c = cov[i];
      // if I don't update at initialization, I would need to do so here
      for(int j = i; j < len; j++){
         int k = j - i;
         c = c + df[j]*dx[k] + df[k]*dx[j];
         double cor = c*invn[j]*invn[k];
         if(mp[j] < cor){
            mp[j] = cor;
            mpi[j] = k;
         }
         if(mp[k] < cor){
            mp[k] = cor;
            mpi[k] = j;
         }
      }
   }
}

