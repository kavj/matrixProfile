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


// The seemingly repeated loops in the following section are to enable partial auto - vectorization over a very expensive region. 
// If they're structured with nesting instead, most compilers are unable to vectorize them. Unfortunately the STL can't really do vectorized reductions, particularly not ones like this.

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
      for(long long i = 0; i < 64; i++){ 
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

// the low level parts with pointer nonsense will always be ugly, but the naming could be improved,
// and we can hide more tiling calculations with the available striding methods. The buffer objects themselves could probably be improved somewhat.
// Will think on this

int nautocorr_reduc(dbuf& ts, dbuf& mp, ibuf& mpi, long long minlag, long long sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   long long mlen = ts.len - sublen + 1;
   long long tlen = std::max(8 * klen, 4 * sublen - (4 * sublen) % klen);        
   ts.set_stride(tlen);
   mp.set_stride(tlen);
   mpi.set_stride(tlen);
   long long tiles = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   long long qstride = paddedlen(sublen);
   dbuf mu(mlen, tlen); dbuf invn(mlen, tlen); dbuf df(mlen, tlen); 
   dbuf dg(mlen, tlen); dbuf cov(mlen, tlen);  dbuf q(tiles * qstride, qstride);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   
   sw_mean(ts(0), mu(0), ts.len, sublen);
   sw_inv_meancentered_norm(ts(0), mu(0), invn(0), ts.len, sublen);

   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);

   #pragma omp parallel for
   for(long long i = 0; i < tiles; i++){
      center_query(ts(i), mu(i), q(i), sublen); 
   }
   for(long long diag = 0; diag < tiles; diag++){
      #pragma omp parallel for 
      for(long long ofst = 0; ofst < tiles - diag; ofst++){
         long long dlim = std::min(tlen, mlen - (diag + ofst) * tlen - minlag);
         // shouldn't batchcov's ts param be offset by minlag as well? 
	 batchcov(ts(diag + ofst, minlag), mu(diag + ofst, minlag), q(ofst), cov(ofst), dlim, sublen);
	 for(long long d = 0; d < dlim; d += klen){
            if(d + klen <= dlim){
               long long ral = std::min(tlen, mlen  - (diag + ofst) * tlen - minlag - d - klen);
	       ral = (ral > 0) ? ral : 0;   
               nautocorr_reduc_kern(cov(ofst, d), mp(ofst), mpi(ofst), df(ofst), dg(ofst), invn(ofst), ofst * tlen, d, ral);
               if(ral < tlen){
                  nautocorr_reduc_edge(cov(ofst, d), mp(0), mpi(0), df(0), dg(0), invn(0), ofst * tlen + ral, diag * tlen + minlag + d, diag * tlen + minlag + d + klen, mlen, false);
               }
            }
            else{
               nautocorr_reduc_edge(cov(ofst, d), mp(0), mpi(0), df(0), dg(0), invn(0), ofst * tlen, diag * tlen + minlag + d, dlim, mlen, true);
            }
         }
      }
   }
   return errs::none;
}


int ncrosscorr_rowwise_reduc(dbuf &a, dbuf &b, dbuf &mp, dbuf &mpi, long long sublen){
   long long mlena = a.len - sublen + 1;
   long long mlenb = b.len - sublen + 1;
   long long tlen = std::max(8 * klen, 4 * sublen - (4 * sublen) % klen);
   long long tilesa = mlena / tlen + (mlena % tlen ? 1 : 0);
   long long tilesb = mlenb / tlen + (mlenb % tlen ? 1 : 0);

   dbuf cov(std::max(mlena, mlenb), tlen);
   dbuf mua(mlena, tlen); dbuf invna(mlena, tlen); dbuf dfa(mlena, tlen); dbuf dga(mlena , tlen); 
   dbuf mub(mlenb, tlen); dbuf invnb(mlenb, tlen); dbuf dfb(mlenb, tlen); dbuf dgb(mlenb, tlen);
   dbuf q(paddedlen(sublen) * std::max(tilesa, tilesb), sublen); 

   sw_mean(a(0), mua(0), a.len, sublen);
   sw_mean(b(0), mub(0), b.len, sublen);

   sw_inv_meancentered_norm(a(0), mua(0), invna(0), a.len, sublen);
   sw_inv_meancentered_norm(b(0), mub(0), invnb(0), b.len, sublen);

   dfdg_init(a(0), mua(0), dfa(0), dga(0), a.len, sublen);
   dfdg_init(b(0), mub(0), dfb(0), dgb(0), b.len, sublen);
   
   #pragma omp parallel for
   for(long long i = 0; i < tilesb; i++){
      center_query(b(i * tlen), mub(i), q(i), sublen);
   }
   for(long long ia = 0; ia < tilesa; ia++){
      //const long long mx = std:min(tilesperdima - ia + 1, tilesperdimb);
      #pragma omp parallel for 
      for(long long ib = 0; ib < tilesb; ib++){
         long long alim = std::min(mlena - ia * tlen, tlen);
	 batchcov(a(ia), mua(ia), q(ib), cov(ia), alim, sublen);
         long long aligned_alim = alim - alim % klen;
	 for(long long d = ia * tlen; d < aligned_alim; d += klen){
            
	 }
         // call function for aligned_alim through alim
      }
   }
}


