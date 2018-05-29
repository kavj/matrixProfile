#include<algorithm>
#include "../utils/max_reduce.h"
#include "../utils/reg.h"
#include "../utils/cov.h"
#include "descriptors.h"
//#define prefalign 32 
#define klen 64 
#define step 32
#define unroll 8
#define simlen 4

// rename and reorganize later. This should be in a different namespace or compilation unit
static inline void  pauto_pearson_AVX_kern (double* __restrict__ cov, double* __restrict__ mp, long long* __restrict__ mpi, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, int offsetr, int offsetc){
   cov = (double*)__builtin_assume_aligned(cov,prefalign);
   df  = (const double*)__builtin_assume_aligned(df,prefalign);
   dx  = (const double*)__builtin_assume_aligned(dx,prefalign);
   invn   = (const double*)__builtin_assume_aligned(invn,prefalign);
   mp  = (double*)__builtin_assume_aligned(mp,prefalign);
   mpi = (long long*)__builtin_assume_aligned(mpi,prefalign);
   for(int i = 0; i < klen; i++){
      for(int j = 0; j < klen; j+=step){
         block<__m256d> cov_r;
         for(int k = 0; k < unroll; k++){
            cov_r(k) = aload(cov,simlen*(j+k));
         }
         __m256d q = brdcst(dx,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(df,i+j+simlen*k),cov_r(k));
         }
         q = brdcst(df,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(dx,i+j+simlen*k),cov_r(k));
            astore(cov_r(k),cov,simlen*(j+k));
         }
         q = brdcst(invn,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= q;
         }
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= uload(invn,i+j+simlen*k);
         }
         block<__m256i> mask;
         for(int k = 0; k < unroll; k++){
            mask(k) = cov_r(k) > uload(mp,i+j+simlen*k);
         }
         __m256i s = brdcst(i);
         for(int k = 0; k < unroll; k++){
            if(testnz(mask(k))){
               maskstore(cov_r(k),mask(k),mp+i+j+simlen*k);
               maskstore(s,mask(k),mpi+i+simlen*k);
            }
         }
         struct rpair r = max_reduce_8x1(cov_r(0),cov_r(1),cov_r(2),cov_r(3),cov_r(4),cov_r(5),cov_r(6),cov_r(7));
         //Todo: technically this is incorrect, as we're writing back 4 operands rather than 1. It would be possible to do a final comparison against a temporary buffer and make one sweep at the end
         // alternatively use the if statement to determine if we need to do a horizontal max (since that is quite expensive)
         __m256i msk = r.val > brdcst(mp,i+j);
         if(testnz(msk)){
            maskstore(r.val,msk,mp+i+j);
            maskstore(r.index+brdcst(i+j+offsetc),msk,mpi+i+j);
         }
      }
   }
}



static inline  __attribute__((always_inline)) void pauto_pearson_refkern_initial (
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int offsetr,
   const int offsetc)
{
    cov =  (double*)__builtin_assume_aligned(cov,prefalign);
    mp =   (double*)__builtin_assume_aligned(mp,prefalign);
    mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
    df =   (const double*)__builtin_assume_aligned(df,prefalign);
    dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
    invn = (const double*)__builtin_assume_aligned(invn,prefalign);

    // Note this goes by row then diagonal/column, so the inner block 
    // determines the column. This messes with our normal variable naming convention 
    // it's necessary so that values of cov may be updated independently

   double corr[klen];
   for(int j = 0; j < klen; j++){
      corr[j] = cov[j]*invn[j+offsetc];
   }
   for(int j = 0; j < klen; j++){
      corr[j] *= invn[0];
   }
   for(int i = 0; i < klen; i++){
      for(int j = 1; j < klen; j++){
         cov[j] += dg[i]*df[i+j+offsetc];
      }
      for(int j = 1; j < klen; j++){
         cov[j] += df[i]*dg[i+j+offsetc];
      } 
      for(int j = 1; j < klen; j++){
         corr[j] = cov[j]*invn[i+j+offsetc];
      }
      for(int j = 1; j < klen; j++){
         corr[j] *= invn[i];
      }
      for(int j = 1; j < klen; j++){
         if(mp[i] < corr[j]){
            mp[i] = corr[j];
            mpi[i] = i+j+offsetr+offsetc;
         }
      }
      for(int j = 1; j < klen; j++){
         if(mp[i+j+offsetc] < corr[j]){
            mp[i+j+offsetc] = corr[j];
            mpi[i+j+offsetc] = i+offsetr;
         }
      }
   }
}




static inline  __attribute__((always_inline)) void pauto_pearson_refkern (
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int offsetr,
   const int offsetc)
{
    cov =  (double*)__builtin_assume_aligned(cov,prefalign);
    mp =   (double*)__builtin_assume_aligned(mp,prefalign);
    mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
    df =   (const double*)__builtin_assume_aligned(df,prefalign);
    dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
    invn = (const double*)__builtin_assume_aligned(invn,prefalign);

    // Note this goes by row then diagonal/column, so the inner block 
    // determines the column. This messes with our normal variable naming convention 
    // it's necessary so that values of cov may be updated independently

    for(int i = 0; i < klen; i++){
      for(int j = 0; j < klen; j++){
         cov[j] += dg[i]*df[i+j+offsetc];
      }
      for(int j = 0; j < klen; j++){
         cov[j] += df[i]*dg[i+j+offsetc];
      } 
      double corr[klen];
      for(int j = 0; j < klen; j++){
         corr[j] = cov[j]*invn[i+j+offsetc];
      }
      for(int j = 0; j < klen; j++){
         corr[j] *= invn[i];
      }
      for(int j = 0; j < klen; j++){
         if(mp[i] < corr[j]){
            mp[i] = corr[j];
            mpi[i] = i+j+offsetr+offsetc;
         }
      }
      for(int j = 0; j < klen; j++){
         if(mp[i+j+offsetc] < corr[j]){
            mp[i+j+offsetc] = corr[j];
            mpi[i+j+offsetc] = i+offsetr;
         }
      }
   }
}


// annoying discrepancy between int and long long. May need to template these
// The problem is that int is generally ideal unless long long happens to match the size of a particular mask
static inline void pauto_pearson_edgekern(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   int*          __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int offsetr, 
   int offsetc, 
   int bound)
{
    cov =  (double*)__builtin_assume_aligned(cov,prefalign);
    mp =   (double*)__builtin_assume_aligned(mp,prefalign);
    mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
    df =   (const double*)__builtin_assume_aligned(df,prefalign);
    dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
    invn = (const double*)__builtin_assume_aligned(invn,prefalign);
   
 
   int dlim = std::min(klen,bound);
   for(int i = 0; i < dlim; i++){
      double c = cov[i];
      for(int j = i; j < bound-i; j++){
         c += df[j]*dg[j+offsetc];
         c += df[j+offsetc]*dg[j];
         double corr = c*invn[i]*invn[i+j];
         if(corr > mp[i]){
            mp[j] = corr;
            mpi[j] = j+offsetr+offsetc;
         }
         if(corr > mp[i+j+offsetc]){
            mp[j+offsetc] = corr;
            mpi[j+offsetc] = j+offsetr;
         }
      }
      cov[i] = c;
   }
}

//Todo: make preambles compiler agnostic



void pauto_pearson(
  double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ ts,
   const double* __restrict__ mu,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int mlen,
   const int sublen,
   const int minlag)
{
   mu  = (double*) __builtin_assume_aligned(mu,prefalign);
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mp =   (double*)__builtin_assume_aligned(mp,prefalign);
   mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   const int tlen = 65536; // add in dynamic formula later
   const int kcount = tlen/klen;
   int tilesperdim = (mlen-minlag)/tlen;
   int fringe = mlen - minlag - tilesperdim*tlen;
   int tailkcount = kcount;
   if(fringe != 0){
      tailkcount = fringe/klen;
      fringe -= tailkcount*klen;
      ++tilesperdim;
   }
   stridedbuf<double> q(tlen,tilesperdim);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts+i*tlen, mu+i*tlen, q(i), sublen);
   }
   for(int d = 0; d < tilesperdim; d += tlen){
      #pragma omp parallel for
      for(int r = 0; r < tilesperdim - d; r += tlen){
         bool inner = d + r + 2 < tilesperdim;
         bool aligned_edge = d + r + 2 == tilesperdim;
         int sd_aligned = (inner || aligned_edge) ? kcount : tailkcount;
         if(inner || (aligned_edge && (tailkcount != 0))){
            // regular initialization 
         }
         else{
            // use tapered initialization
         }
         for(int sd = 0; sd < sd_aligned; sd++){
            int sr_aligned = (inner || (aligned_edge && (sd < tailkcount))) ? kcount : tailkcount - sd;
            for(int sr = 1; sr < srmx; sr++){
               pauto_pearson_inner();
            } 
            if((fringe != 0) && !(inner || (aligned_edge && (sd < tailcount)))){
               pauto_pearson_edge();
            }
         }
      }
   }
}



void pauto_pearson_reftest(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   int minlag,
   int mlen)
{

   for(int i = minlag; i < mlen; i++){
      for(int j = 0; j < mlen-i; j++){
         cov[i] += df[j]*dg[i+j];
         cov[i] += df[i+j]*dg[j];
         double corr = cov[i]*invn[j]*invn[i+j];
         if(mp[j] < corr){
            mp[j] = corr;
            mpi[j] = i+j;
         }
         if(mp[i+j] < corr){
            mp[i+j] = corr;
            mpi[i+j] = i;
         }
      }
   }
}

