#include <iostream>
//#include "../utils/xprec.h"
//#include "../utils/cov.h"
#include "../utils/alloc.h"
#include <omp.h>
#include <array>
#include "pearson.h"
#include "../utils/moments.h"
constexpr long long klen  = 256; 

static void dfdg_init(const double* __restrict__ ts, const double* __restrict__ mu, double* __restrict__ df, double* __restrict__ dg, long long len, long long sublen){
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - ts[i])/2.0;
      dg[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
   }
}

static inline void pauto_pearson_kern(
   double*       __restrict__  cov,
   double*       __restrict__  mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
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
      std::array<double, klen> cr;
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
      std::array<long long, klen/2> ci;
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
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   const long long ofr, 
   const long long ofd, 
   const long long dlim,
   const long long clim,
   const bool init)
{
   for(long long r = ofr; r < clim - ofd; r++){
      for(long long d = ofd; d < std::min(clim - r, dlim); d++){
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

int partialauto(bufd& ts, bufd& mp, bufi& mpi, long long minlag, long long sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   const int qstride = paddedlen(sublen, prefalign);
   const long long mlen = ts.len - sublen + 1;
   const long long tlen = std::max(static_cast<long long>(2 << 14), 4 * sublen - (4 * sublen) % klen);        
   const long long tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   bufd mu(mlen); bufd invn(mlen); bufd df(mlen);  
   bufd dg(mlen); bufd cov(mlen);  bufd q(tilesperdim * qstride);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   
   sw_mean(ts(), mu(), ts.len, sublen);
   sw_inv_meancentered_norm(ts(), mu(), invn(), ts.len, sublen);

   dfdg_init(ts(), mu(), df(), dg(), ts.len, sublen);
   #pragma omp parallel for
   for(long long i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i * qstride), sublen); 
   }
   for(long long diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for 
      for(long long ofst = 0; ofst < tilesperdim - diag; ofst++){
         const long long di = diag * tlen + minlag;
         const long long ofi = ofst * tlen;
         const long long dlim = std::min(di + tlen, mlen - ofi);
         autocov(ts(di + ofi), mu(di + ofi), q(ofst * qstride), cov(ofi), dlim - di, sublen);
         for(long long d = di; d < dlim; d += klen){
            if(d + klen <= dlim){
               const long long ral = std::min(tlen, mlen - d - ofi - klen + 1); 
               pauto_pearson_kern(cov(ofi + d - di), mp(ofi), mpi(ofi), df(ofi), dg(ofi), invn(ofi), ofi, d, ral);
               if(ral < tlen){
                  pauto_pearson_edge(cov(ofi + d - di), mp(), mpi(), df(), dg(), invn(), ofi + ral, d, d + klen, mlen, false);
               }
            }
            else{
               pauto_pearson_edge(cov(ofi + d - di), mp(), mpi(), df(), dg(), invn(), ofi, d, dlim, mlen, true);
            }
         }
      }
   }
   return errs::none;
}

