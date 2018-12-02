#include<stdbool.h>
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "../utils/alloc.h"
#include <omp.h>
#include "pearson.h"
const long long klen  = 256; 

// Note: The use of long long here allows better compiler generated simd vectorization on architectures where sizeof(double) == sizeof(long long).

static void dfdg_init(const double* restrict ts, const double* restrict mu, double* restrict df, double* restrict dg, long long len, long long sublen){
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - ts[i])/2.0;
      dg[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
   }
}

static inline void pauto_pearson_kern(
   double*       restrict  cov,
   double*       restrict  mp,
   long long*    restrict mpi,
   const double* restrict df,
   const double* restrict dg,
   const double* restrict invn,
   const long long ofr,
   const long long ofc,
   const long long iters)
{
   df   = (const double*)__builtin_assume_aligned(df, prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (long long*) __builtin_assume_aligned(mpi, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);
   for(long long r = 0; r < iters; r++){
      double cr[klen];
      if(r != 0){
         for(long long d = 0; d < klen; d++){ 
            cov[d] += df[r] * dg[r + d + ofc];
            cov[d] += df[r + d + ofc] * dg[r];
         }
      }
      for(long long d = 0; d < klen; d++){
         cr[d] = cov[d] * invn[r] * invn[r + d + ofc];
         if(cov[d] * invn[r] * invn[r + d + ofc] > mp[r + d + ofc]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         }
      }
      long long ci[klen/2];
      for(long long i = 0; i < 128; i++){
         ci[i] = cr[i] > cr[i + 128] ? i : i + 128;
         cr[i] = cr[i] > cr[i + 128] ? cr[i] : cr[i + 128];
      }
      for(long long i = 0; i < 64; i++){ // gcc can generally simplify conditional writes using if statements to masked writes, but in this case
                                         // it generates all scalar code. For now I rewrote it using ternary operators. I'll file a bug report later.
         ci[i] = cr[i] > cr[i + 64] ? ci[i] : ci[i + 64];
         cr[i] = cr[i] > cr[i + 64] ? cr[i] : cr[i + 64];
      }
      for(long long i = 0; i < 32; i++){
         ci[i] = cr[i] > cr[i + 32] ? ci[i] : ci[i + 32];
         cr[i] = cr[i] > cr[i + 32] ? cr[i] : cr[i + 32];
      } 
      for(long long i = 0; i < 16; i++){
         ci[i] = cr[i] > cr[i + 16] ? ci[i] : ci[i + 16];
         cr[i] = cr[i] > cr[i + 16] ? cr[i] : cr[i + 16];
      }
      for(long long i = 0; i < 8; i++){
         ci[i] = cr[i] > cr[i + 8] ? ci[i] : ci[i + 8];
         cr[i] = cr[i] > cr[i + 8] ? cr[i] : cr[i + 8];
      }
      for(long long i = 0; i < 4; i++){
         ci[i] = cr[i] > cr[i + 4] ? ci[i] : ci[i + 4];
         cr[i] = cr[i] > cr[i + 4] ? cr[i] : cr[i + 4];
      }
      for(long long i = 0; i < 2; i++){
         ci[i] = cr[i] > cr[i + 2] ? ci[i] : ci[i + 2];
         cr[i] = cr[i] > cr[i + 2] ? cr[i] : cr[i+2];
      }
      ci[0] = cr[0] > cr[1] ? ci[0] : ci[1];
      cr[0] = cr[0] > cr[1] ? cr[0] : cr[1];
      mpi[r] = mp[r] > cr[0] ? mpi[r] : ci[0] + r + ofr + ofc;
      mp[r] =  mp[r] > cr[0] ? mp[r] : cr[0];
   }
}

static inline void pauto_pearson_edge(
   double*       restrict cov, 
   double*       restrict mp,  
   long long*    restrict mpi, 
   const double* restrict df,  
   const double* restrict dg, 
   const double* restrict invn, 
   const long long ofr, 
   const long long ofd, 
   const long long dlim,
   const long long clim,
   const bool init)
{
   for(long long r = ofr; r < clim - ofd; r++){
      long long constraineddlim = (clim - r) < dlim ? (clim - r) : dlim; 
      for(long long d = ofd; d < constraineddlim; d++){
         if(!init || r != ofr){
            cov[d - ofd] += df[r] * dg[r + d];
            cov[d - ofd] += df[r + d] * dg[r];
         }
         double cr = cov[d - ofd] * invn[r] * invn[r + d];
         if(mp[r] < cr){
            mp[r] = cr;
            mpi[r] = r + d;
         }
         if(mp[r + d] < cr){
            mp[r + d] = cr;
            mpi[r + d] = r;
         }
      }
   }
}

mprofile pearson_pauto_reduc(const double* ts, long long len, long long minlag, long long sublen){
   if(ts == NULL){ 
      //printf("invalid arguments\n");
      exit(1);
   }
   else if(sublen + minlag >= len){
      //printf("time series is too short relative to chosen subsequence length and minimum acyclic lag factor\n");
      exit(1);
   }
   const long long mlen = len - sublen + 1;
   const long long tlen = (2 << 14) > (4 * sublen - (4 * sublen) % klen) ? (2 << 14) : (4 * sublen - (4 * sublen) % klen);
   const long long tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   const long long qstride = paddedlen(sublen);
   double* mu = alloc_buff(mlen);
   double* invn = alloc_buff(mlen);
   double* df = alloc_buff(mlen);
   double* dg = alloc_buff(mlen);
   double* cov = alloc_buff(mlen);
   double* mp = alloc_buff(mlen);
   long long* mpi = alloc_buff(mlen);
   double* q = alloc_buff(tilesperdim * qstride);
   // update with a more sensible return
   if((mu == NULL) || (invn == NULL) || (df == NULL) || (dg == NULL) || (cov == NULL) || (mp == NULL) || (mpi == NULL) || (q == NULL)){
      //printf("unable to allocate memory\n");
      exit(1);
   }

   xmean_windowed(ts, mu, len, sublen);
   xInvn(ts, mu, invn, len, sublen);
   dfdg_init(ts, mu, df, dg, len, sublen);
  
   #pragma omp parallel for
   for(long long i = 0; i < tilesperdim; i++){
      center_query(&ts[i * tlen], &mu[i * tlen], &q[i * qstride], sublen); 
   }
   for(long long diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for 
      for(long long ofst = 0; ofst < tilesperdim - diag; ofst++){
         const long long di = diag * tlen + minlag;
         const long long ofi = ofst * tlen;
         const long long dlim = di + tlen < mlen - ofi ? (di + tlen) : (mlen - ofi); 
         batchcov(&ts[di + ofi], &mu[di + ofi], &q[ofst * qstride], &cov[ofi], dlim - di, sublen);
         for(long long d = di; d < dlim; d += klen){
            if(d + klen <= dlim){
               const long long r = tlen < (mlen - d - ofi - klen) ? tlen : (mlen - d - ofi - klen);
	       const long long ral = r > 0 ? r : 0;
               pauto_pearson_kern(&cov[ofi + d - di], &mp[ofi], &mpi[ofi], &df[ofi], &dg[ofi], &invn[ofi], ofi, d, ral);
               if(ral < tlen){
                  pauto_pearson_edge(&cov[ofi + d - di], mp, mpi, df, dg, invn, ofi + ral, d, d + klen, mlen, false);
               }
            }
            else{
               pauto_pearson_edge(&cov[ofi + d - di], mp, mpi, df, dg, invn, ofi, d, dlim, mlen, true);
            }
         }
      }
   }
   dealloc_buff(mu);
   dealloc_buff(invn);
   dealloc_buff(df);
   dealloc_buff(dg);
   dealloc_buff(cov);
   dealloc_buff(q);
   mprofile p = {mp, mpi};
   return p;
}

