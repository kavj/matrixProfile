#include<algorithm>
#include"../utils/max_reduce.h"
#include "prescr_Pearson.h"
#include "auto_Pearson.h"
#include "../utils/reg.h"


#define tsz 64
#define simlen 4
#define unroll 8


struct rpair rescaled_max_reduct(double* __restrict__ cov,  const double* __restrict__ invn, double* __restrict__ xcorr, long long* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset){
   cov = (double*)__builtin_assume_aligned(cov,32);
   invn = (const double*)__builtin_assume_aligned(invn,32);
   xcorr = (double*)__builtin_assume_aligned(xcorr,32);
   cindex = (long long*)__builtin_assume_aligned(cindex,32);
   for(int i = 0; i < tsz; i++){
      block<__m256d> cov_r;
      __m256d q = brdcst(qinvn);
      for(int k = 0; k < unroll; k++){
         cov_r(k) = q*aload(cov,simlen*k);
      }
      for(int k = 0; k < unroll; k++){
         cov_r(k) *= uload(invn,i+simlen*k);
      }
      block<__m256i> mask;
      for(int k = 0; k < unroll; k++){
         mask(k) = cov_r(k) > uload(xcorr,i+simlen*k);
      }
      for(int k = 0; k < unroll; k++){
         __m256i r = brdcst(i);
         if(testnz(mask(k))){
            maskstore(cov_r(k),mask(k),xcorr+i+simlen*k);
            maskstore(r,mask(k),cindex+i+simlen*k);
         }
      }
      struct rpair r = max_reduce_8x1(cov_r(0),cov_r(1),cov_r(2),cov_r(3),cov_r(4),cov_r(5),cov_r(6),cov_r(7));
      // should just convert q to a vector
      __m256i msk = r.val > brdcst(xcorr,i);
      if(testnz(msk)){
         
         maskstore(r.val,msk,xcorr+tsz);
         maskstore(r.index,msk,cindex+tsz);
      }
   }
}


struct query_stat rescaled_max_reduct_ref(double* __restrict__ cov,  const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset){
   cov = (double*)__builtin_assume_aligned(cov,32);
   invn = (const double*) __builtin_assume_aligned(invn,32);
   xcorr = (double*) __builtin_assume_aligned(xcorr,32);
   cindex = (int*) __builtin_assume_aligned(cindex,32);
   int ind = -1;
   double a[tsz];
   for(int i = 0; i < tsz; i++){
      a[i] = cov[i]*qinvn*invn[i];
   }
   for(int i = 0; i < tsz; i++){
      if(xcorr[i] < a[i]){
         xcorr[i] = a[i];
         cindex[i] = qbaseind;
      }
   }
   int j = 0;
   for(int i = 0; i < tsz; i+=4){
      if(a[j] < a[i]){
         j = i;
      }
   }
   struct query_stat s;
   s.qcov = cov[j];
   s.qcorr = a[j];
   s.qind = j+cindoffset;
   return s;
}



/* qcov follows the partial cross covariance formula without rescaling applied. This form assumes that consecutive terms in qcov correspond to strided queries, where stride indicates the distance between consecutive entries*/

void maxpearson_ext_auto(const double* __restrict__ qcov, const double* __restrict__ invn, const double* __restrict__ df, const double* __restrict__ dx, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len){
   for(int i = 0; i < count; i++){
      int ci = qind[i];  // need to know the index 
      int qi = i*stride;
      int lim = std::min(extraplen,std::min(ci,qi));
      double c = qcov[i];
      for(int j = 1; j < lim; j++){
         c -= dx[ci-j+1]*df[qi-j+1];
         c -= dx[qi-j+1]*df[ci-j+1];
         double corr = c*invn[qi-j]*invn[ci-j];
         if(mp[ci-j] < corr){
            mp[ci-j] = corr;
            mpi[ci-j] = qi-j;
         }
         if(mp[qi] < corr){
            mp[qi-j] = corr;
            mpi[qi-j] = ci-j;
         }
      }
      lim = std::min(extraplen,len-std::max(ci,qi));
      c = qcov[i];
      for(int j = 1; j < lim; j++){
         c += dx[ci+j]*df[qi+j];
         c += dx[qi+j]*df[ci+j]; 
         double corr = c*invn[qi+j]*invn[ci+j];
         if(mp[ci+j] < corr){
            mp[ci+j] = corr;
            mpi[ci+j] = qi+j;
         }
         if(mp[qi+j] < corr){
            mp[qi+j] = corr;
            mpi[qi+j] = ci+j;
         }
      }
   }
}


