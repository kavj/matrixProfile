#include <iostream>
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "../utils/alloc.h"
#include <omp.h>
#include <array>
#include "pearson.h"

// at the moment the kernel functions rely on this being 256, since that allows auto-vectorization across most up to date compilers without writing or generating lots of intrinsics based code.

constexpr long long klen  = 256; 


void pearson2zned(double* __restrict mp, long long len, long long sublen){
   mp = (double*) __builtin_assume_aligned(mp, prefalign);
   double scale = 2 * sublen;
   #pragma omp parallel for
   for(long long i = 0; i < len; i++){
      mp[i] = sqrt(scale * (1.0 - mp[i]));
   }
}

static void dfdg_init(const double* __restrict ts, const double* __restrict mu, double* __restrict df, double* __restrict dg, long long len, long long sublen){
   df   = (double*)__builtin_assume_aligned(df, prefalign);
   dg   = (double*)__builtin_assume_aligned(dg, prefalign);
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - ts[i])/2.0;
      dg[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
   }
}

static inline void nautocorr_reduc_kern(
   double*       __restrict  cov,
   double*       __restrict  mp,
   long long*    __restrict mpi,
   const double* __restrict df,
   const double* __restrict dg,
   const double* __restrict invn,
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

static inline void nautocorr_reduc_edge(
   double*       __restrict cov, 
   double*       __restrict mp,  
   long long*    __restrict mpi, 
   const double* __restrict df,  
   const double* __restrict dg, 
   const double* __restrict invn, 
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

int nautocorr_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, long long minlag, long long sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   const long long mlen = ts.len - sublen + 1;
   const long long tlen = std::max(8 * klen, 4 * sublen - (4 * sublen) % klen);        
   const long long tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   
   sw_mean(ts(0), mu(0), ts.len, sublen);
   sw_inv_meancentered_norm(ts(0), mu(0), invn(0), ts.len, sublen);
   //xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(long long i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i), sublen); 
   }
   for(long long diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for 
      for(long long ofst = 0; ofst < tilesperdim - diag; ofst++){
    	 const long long di = diag * tlen + minlag;
         const long long ofi = ofst * tlen;
         const long long dlim = std::min(di + tlen, mlen - ofi);
         batchcov(ts(di + ofi), mu(di + ofi), q(ofst), cov(ofi), dlim - di, sublen);
         for(long long d = di; d < dlim; d += klen){
            if(d + klen <= dlim){
               const long long ral = std::max(static_cast<long long>(0), std::min(tlen, mlen - d - ofi - klen)); // stupid compiler
               nautocorr_reduc_kern(cov(ofi + d - di), mp(ofi), mpi(ofi), df(ofi), dg(ofi), invn(ofi), ofi, d, ral);
               if(ral < tlen){
                  nautocorr_reduc_edge(cov(ofi + d - di), mp(0), mpi(0), df(0), dg(0), invn(0), ofi + ral, d, d + klen, mlen, false);
               }
            }
            else{
               nautocorr_reduc_edge(cov(ofi + d - di), mp(0), mpi(0), df(0), dg(0), invn(0), ofi, d, dlim, mlen, true);
            }
         }
      }
   }
   return errs::none;
}


/*
int ncrosscorr_rowwise_reduc(dsbuf &a, dsbuf &b, dsbuf &mp, dsbuf &mpi, long long sublen){
   const long long mlen = std::min(a.len, b.len) - sublen + 1;
   const long long tlen = std::max(8 * klen, 4 * sublen - (4 * sublen) % klen);
   const long long tilesperdim = mlen / tlen + (mlen % tlen ? 1 : 0);
   dsbuf mua(mlen); dsbuf invna(mlen); dsbuf dfa(mlen); dsbuf dga(mlen); 
   dsbuf mub(mlen); dsbuf invnb(mlen); dsbuf dfb(mlen); dsbuf dgb(mlen);
   // check valid or more likely factor out common steps?
   //
   dfdg_init(a(0), mua(0), dfa(0), dga(0), a.len, sublen);
   dfdg_init(b(0), mub(0), dfb(0), dgb(0), b.len, sublen);

   #pragma omp parallel for 
   for(long long i = 0; i < tilesperdim; i++){
      
   }
}



int ncrosscorr_colwise_reduc(dsbuf &a, dsbuf &b, dsbuf &mp, dsbuf &mpi, long long sublen){
   const long long mlen = a.len - sublen + 1;
   const long long tlen = std::max(8 * klen, 4 * sublen - (4 * sublen) % klen);
   const long long tilesperdim = mlen / tlen + (mlen % tlen ? 1 : 0);
   dsbuf mua(mlen); dsbuf invna(mlen); dsbuf dfa(mlen); dsbuf dga(mlen); 
   dsbuf mub(mlen); dsbuf invnb(mlen); dsbuf dfb(mlen); dsbuf dgb(mlen);
   // check valid or more likely factor out common steps?
   //
   #pragma omp parallel for 
   for(long long i = 0; i < tilesperdim; i++){

   }
}

*/
