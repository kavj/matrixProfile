#include <array>
#include<immintrin.h>
#include "pearson.h"
#include "alloc.h"
#include "moments.h"
constexpr int klen  = 64; 

static void dfdg_init(const double* __restrict__ ts, const double* __restrict__ mu, double* __restrict__ df, double* __restrict__ dg, int len, int sublen){
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - ts[i])/2.0;
      dg[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
   }
}


// These sections tend to be difficult for an optimizing compiler, thus the code duplication. If I use a cleaner nested loop structure,
// gcc and clang are unable to eliminate some redundnant loop comparisons. The use of inline functions also 
// extends the scope of the temporary arrays, resulting in a typical loss of 8% or so on the entire runtime (not just this section).
// Since the aligned kernel sections have a sufficiently regular structure, I may consider the use of jited code, particularly 
// as this is close to being compliant with C99 or C11 anyway.


 void pavx_kern(
   double*       __restrict cv,
   double*       __restrict mp,
   const double* __restrict df,
   const double* __restrict dg,
   const double* __restrict invn,
   const int ofr,
   const int colOf,
   const int iters,
   const int tlen)
{

   cv  =  (double*)      __builtin_assume_aligned(cv, prefalign);
   mp   = (double*)      __builtin_assume_aligned(mp, prefalign);
   df   = (const double*)__builtin_assume_aligned(df, prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);

   __m256d c0 = _mm256_load_pd(cv);
   __m256d c1 = _mm256_load_pd(cv + 4);
   __m256d c2 = _mm256_load_pd(cv + 8);
   __m256d c3 = _mm256_load_pd(cv + 12);
   __m256d c4 = _mm256_load_pd(cv + 16);
   __m256d c5 = _mm256_load_pd(cv + 20);
   __m256d c6 = _mm256_load_pd(cv + 24);
   __m256d c7 = _mm256_load_pd(cv + 28);
   #pragma omp simd aligned(cv, df, dg, invn : 64)
   for(int i = 0; i < iters; i++){
      __m256d df_ = _mm256_broadcast_sd(df + i);
      if(i > 0){
         c0 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf), c0);
         c1 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 4), c1);
         c2 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 8), c2);
         c3 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 12), c3);
         c4 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 16), c4);
         c5 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 20), c5);
         c6 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 24), c6);
         c7 = _mm256_fmadd_pd(df_, _mm256_loadu_pd(dg + i + colOf + 28), c7);

         __m256d dg_ = _mm256_broadcast_sd(dg + i);

         c0 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf), c0);
         c1 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 4), c1);
         c2 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 8), c2);
         c3 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 12),c3);
         c4 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 16),c4);
         c5 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 20),c5);
         c6 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 24),c6);
         c7 = _mm256_fmadd_pd(dg_, _mm256_loadu_pd(df + i + colOf + 28),c7);
      }
   
      __m256d invn_ = _mm256_broadcast_sd(invn + i);
      __m256d f0 = _mm256_mul_pd(c0, _mm256_loadu_pd(invn + i + colOf));
      __m256d f1 = _mm256_mul_pd(c1, _mm256_loadu_pd(invn + i + colOf + 4));
      __m256d f2 = _mm256_mul_pd(c2, _mm256_loadu_pd(invn + i + colOf + 8));
      __m256d f3 = _mm256_mul_pd(c3, _mm256_loadu_pd(invn + i + colOf + 12));
      __m256d f4 = _mm256_mul_pd(c4, _mm256_loadu_pd(invn + i + colOf + 16));
      __m256d f5 = _mm256_mul_pd(c5, _mm256_loadu_pd(invn + i + colOf + 20));
      __m256d f6 = _mm256_mul_pd(c6, _mm256_loadu_pd(invn + i + colOf + 24));
      __m256d f7 = _mm256_mul_pd(c7, _mm256_loadu_pd(invn + i + colOf + 28));

      __m256d g0 = _mm256_mul_pd(f0, invn_);
      __m256d g1 = _mm256_mul_pd(f1, invn_);
      __m256d g2 = _mm256_mul_pd(f2, invn_);
      __m256d g3 = _mm256_mul_pd(f3, invn_);
      __m256d g4 = _mm256_mul_pd(f4, invn_);
      __m256d g5 = _mm256_mul_pd(f5, invn_);
      __m256d g6 = _mm256_mul_pd(f6, invn_);
      __m256d g7 = _mm256_mul_pd(f7, invn_);

      __m256d h0 = _mm256_max_pd(g0, _mm256_loadu_pd(mp + i + colOf));
      __m256d h1 = _mm256_max_pd(g1, _mm256_loadu_pd(mp + i + colOf + 4));
      __m256d h2 = _mm256_max_pd(g2, _mm256_loadu_pd(mp + i + colOf + 8));
      __m256d h3 = _mm256_max_pd(g3, _mm256_loadu_pd(mp + i + colOf + 12));
      __m256d h4 = _mm256_max_pd(g4, _mm256_loadu_pd(mp + i + colOf + 16));
      __m256d h5 = _mm256_max_pd(g5, _mm256_loadu_pd(mp + i + colOf + 20));
      __m256d h6 = _mm256_max_pd(g6, _mm256_loadu_pd(mp + i + colOf + 24));
      __m256d h7 = _mm256_max_pd(g7, _mm256_loadu_pd(mp + i + colOf + 28));

      _mm256_storeu_pd(mp + i + colOf, h0);
      _mm256_storeu_pd(mp + i + colOf + 4, h1);
      _mm256_storeu_pd(mp + i + colOf + 8, h2);
      _mm256_storeu_pd(mp + i + colOf + 12, h3);
      _mm256_storeu_pd(mp + i + colOf + 16, h4);
      _mm256_storeu_pd(mp + i + colOf + 20, h5);
      _mm256_storeu_pd(mp + i + colOf + 24, h6);
      _mm256_storeu_pd(mp + i + colOf + 28, h7);

      __m256d j0 = _mm256_max_pd(g0, g1);
      __m256d j1 = _mm256_max_pd(g2, g3);
      __m256d j2 = _mm256_max_pd(g4, g5);
      __m256d j3 = _mm256_max_pd(g6, g7);

      __m256d L0 = _mm256_max_pd(j0, j1);
      __m256d L1 = _mm256_max_pd(j2, j3);
   
      __m256d L2 = _mm256_max_pd(L0, L1);

      double q[4];
      _mm256_store_pd(&q[0], L2);
      double q0 = q[0] > q[1] ? q[0] : q[1];
      double q1 = q[2] > q[3] ? q[2] : q[3];
      double q2 = q0 > q1 ? q0 : q1;
      mp[i] = q2 > mp[i] ? q2 : mp[i];
   }
   if(tlen > iters){
      _mm256_store_pd(cv, c0);
      _mm256_store_pd(cv + 4, c1);
      _mm256_store_pd(cv + 8, c2);
      _mm256_store_pd(cv + 12, c3);
      _mm256_store_pd(cv + 16, c4);
      _mm256_store_pd(cv + 20, c5);
      _mm256_store_pd(cv + 24, c6);
      _mm256_store_pd(cv + 28, c7);
   }
}


static inline void pauto_pearson_kern(
   double*       __restrict  cov,
   double*       __restrict  mp,
   int*          __restrict mpi,
   const double* __restrict df,
   const double* __restrict dg,
   const double* __restrict invn,
   const int ofr,
   const int ofc,
   const int iters)
{
   df   = (const double*)__builtin_assume_aligned(df, prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign);
   mp   = (double*)      __builtin_assume_aligned(mp, prefalign);
   mpi  = (int*)         __builtin_assume_aligned(mpi, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);
   cov  = (double*)      __builtin_assume_aligned(cov, prefalign);
   for(int r = 0; r < iters; r++){
      std::array<double, klen> cr;
      if(r != 0){
         for(int d = 0; d < klen; d++){ 
            cov[d] += df[r] * dg[r + d + ofc];
            cov[d] += df[r + d + ofc] * dg[r];
         }
      }
      for(int d = 0; d < klen; d++){
         cr[d] = cov[d] * invn[r] * invn[r + d + ofc];
         if(cov[d] * invn[r] * invn[r + d + ofc] > mp[r + d + ofc]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         }
      }
      std::array<int, klen/2> ci;
      auto maxReduction = [&cr, &ci](int offset){
         for(int i = 0; i < offset; i++){
            ci[i] = cr[i] >= cr[i + offset] ? cr[i] : cr[i + offset];
            cr[i] = cr[i] >= cr[i + offset] ? cr[i] : cr[i + offset];
         }
      };
      for(int i = klen/4; i > 0; i /= 2){
         maxReduction(i);
      }
      mpi[r] = mp[r] >= cr[0] ? mpi[r] : ci[0] + r + ofr + ofc;
      mp[r] =  mp[r] >= cr[0] ? mp[r] : cr[0];
   }
}

static inline void pauto_pearson_edge(
   double*       __restrict cov, 
   double*       __restrict mp,  
   int*          __restrict mpi, 
   const double* __restrict df,  
   const double* __restrict dg, 
   const double* __restrict invn, 
   const int ofr, 
   const int ofd, 
   const int dlim,
   const int clim,
   const bool init)
{
   for(int r = ofr; r < clim - ofd; r++){
      for(int d = ofd; d < std::min(clim - r, dlim); d++){
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

int pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   const int mlen = ts.len - sublen + 1;
   const int tlen = std::max(static_cast<int>(2 << 14), 4 * sublen - (4 * sublen) % klen);        
   const int tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   
   windowedMean(ts(0), mu(0), ts.len, sublen);
   windowedInverseCenteredNorm(ts(0), mu(0), invn(0), ts.len, sublen);
   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      meanCenter(ts(i * tlen), q(i), *mu(i * tlen), sublen); 
   }
   for(int diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for 
      for(int ofst = 0; ofst < tilesperdim - diag; ofst++){
         const int di = diag * tlen + minlag;
         const int ofi = ofst * tlen;
         const int dlim = std::min(di + tlen, mlen - ofi);
         crossCov(ts(di + ofi), mu(di + ofi), q(ofst), cov(ofi), dlim - di, sublen);
         for(int d = di; d < dlim; d += klen){
            const int ral = std::min(tlen, mlen - d - ofi - klen + 1);
            if(ral > 0){
               pavx_kern(cov(ofi + d - di), mp(ofi), df(ofi), dg(ofi), invn(ofi), ofi, d, ral,tlen);
             //  pauto_pearson_kern(cov(ofi + d - di), mp(ofi), mpi(ofi), df(ofi), dg(ofi), invn(ofi), ofi, d, ral);
               if(ral < tlen){
                  pauto_pearson_edge(cov(ofi + d - di), mp(0), mpi(0), df(0), dg(0), invn(0), ofi + ral, d, d + klen, mlen, false);
               }
            }
            else{
               pauto_pearson_edge(cov(ofi + d - di), mp(0), mpi(0), df(0), dg(0), invn(0), ofi, d, dlim, mlen, true);
            }
         }
      }
   }
   return errs::none;
}

