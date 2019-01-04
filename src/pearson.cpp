#include <iostream>
#include <omp.h>
#include <array>
#include <cmath>
#include "pearson.h"
#include "alloc.h"
#include "moments.h"


// at the moment the kernel functions rely on this being 256, since that allows auto-vectorization across most up to date compilers without writing or generating lots of intrinsics based code.

constexpr long long klen  = 256; 
constexpr int maxunroll = 256;

void pearson2zned(double* __restrict mp, long long len, long long sublen){
   mp = (double*) __builtin_assume_aligned(mp, prefalign);
   double scale = 2 * sublen;
   #pragma omp parallel for simd
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


void nautocorr_reduc_kern(
   double*       __restrict cov, 
   double*       __restrict mp,  
          int*   __restrict mpi, 
   const double* __restrict df,  
   const double* __restrict dg, 
   const double* __restrict invn, 
   const long long ofr, 
   const long long ofd, 
   const long long dlim,
   const long long clim)
   {
   cov = (double*) __builtin_assume_aligned(cov,prefalign);
   mp = (double*) __builtin_assume_aligned(mp,prefalign);
   mpi = (int*) __builtin_assume_aligned(mpi,prefalign);
   df = (double*) __builtin_assume_aligned(df,prefalign);
   dg = (double*) __builtin_assume_aligned(dg,prefalign);
   invn = (double*) __builtin_assume_aligned(invn,prefalign);

   for(long long r = ofr; r < clim - ofd; r++){
      std::array<double, maxunroll> cr;
      #pragma omp simd safelen(maxunroll) aligned(cov : prefalign)
      for(long long d = 0; d < maxunroll; d++){
         cov[d] += df[r] * dg[r + d + ofd];
         cov[d] += df[r + d + ofd] * dg[r];
         cr[d] = cov[r + d] * invn[r] * invn[r + d + ofd];
      }
      for(long long d = 0; d < maxunroll; d++){
         if(mp[r + d + ofd] < cr[d]){
            mp[r + d + ofd] = cr[d];
	    mpi[r + d + ofd] = r;
         }
      }
      std::array<int, maxunroll> cri;
      #pragma omp simd safelen(maxunroll/2) 
      for(long long d = 0; d < maxunroll/2; d++){
         cri[d] = (cr[d] > cr[d + maxunroll/2]) ? d + ofd : d + maxunroll/2 + ofd;
         cr[d] =  (cr[d] > cr[d + maxunroll/2]) ? cr[d] : cr[d + maxunroll/2];
      }
      for(long long reduce = maxunroll/4; reduce > 0; reduce /= 2){
	 #pragma omp simd  
         for(long long d = 0; d < reduce; d++){
            cri[d] = (cr[d] > cr[d + reduce]) ? cri[d] : cri[d + reduce];
            cr[d] =  (cr[d] > cr[d + reduce]) ? cr[d] : cr[d + reduce];
         }
      }
      if(mp[r] < cr[0]){
         mp[r] = cr[0];
	 mpi[r] = cri[0];
      }
   }
}


int nautocorr_reduc(double* ts, double* mp, int* mpi, int len, int minlag, int sublen){
   int basestride = 256;
   int tilecount = 1;
   int qstride = paddedlen(sublen);
   double* cv = static_cast<double*>(alloc_aligned_buffer((len - sublen + 1) * sizeof(double)));
   double* mu = static_cast<double*>(alloc_aligned_buffer((len - sublen + 1) * sizeof(double)));
   double* invn = static_cast<double*>(alloc_aligned_buffer((len - sublen + 1) * sizeof(double)));
   double* df = static_cast<double*>(alloc_aligned_buffer((len - sublen + 1) * sizeof(double)));
   double* dg = static_cast<double*>(alloc_aligned_buffer((len - sublen + 1) * sizeof(double)));
   double* query = static_cast<double*>(alloc_aligned_buffer(qstride * sizeof(double)));
   if((ts == nullptr) || (mp == nullptr) || (mpi == nullptr) || (mu == nullptr) ){
      return errs::mem_error; 
   }

   sw_mean(ts, mu, len, sublen);
   sw_inv_meancentered_norm(ts, mu, invn, len, sublen);
   dfdg_init(ts, mu, df, dg, len, sublen);
   int mlen = len - sublen + 1;
   #pragma omp parallel for
   for(int i = 0; i * basestride < mlen; i++){
      //center_query(ts, mu, query, i * basestride, i * qstride, sublen); // update func to allow offsets to be passed in
   }

   for(int diag = 0; diag < mlen; diag += basestride){
      #pragma omp parallel for
      for(int offset = 0; offset < mlen - diag; offset += basestride){
         //batchcov(ts, mu, query, cv, count, winlen); // should be updated for an offset range too
	 //nautocorr_reduc_kern(cv, df, dg, invn, mp, mpi, ...);
      }
   }

   dealloc_aligned_buffer(cv);
   dealloc_aligned_buffer(mu);
   dealloc_aligned_buffer(invn);
   dealloc_aligned_buffer(df);
   dealloc_aligned_buffer(dg);
   dealloc_aligned_buffer(query);
   return errs::none;
}


/*
int ncrosscorr_rowwise_reduc(double* a, double* b, double* mp, double* mpi, int lena, int lenb, int sublen){
   long long mlena = len - sublen + 1;
   long long mlenb = len - sublen + 1;
   long long tlen = std::max(8 * klen, 4 * sublen - (4 * sublen) % klen);
   long long tilesa = mlena / tlen + (mlena % tlen ? 1 : 0);
   long long tilesb = mlenb / tlen + (mlenb % tlen ? 1 : 0);

   dbuf cov(std::max(mlena, mlenb), tlen);
   dbuf mua(mlena, tlen); dbuf invna(mlena, tlen); dbuf dfa(mlena, tlen); dbuf dga(mlena , tlen); 
   dbuf mub(mlenb, tlen); dbuf invnb(mlenb, tlen); dbuf dfb(mlenb, tlen); dbuf dgb(mlenb, tlen);
   dbuf q(paddedlen(sublen) * std::max(tilesa, tilesb), sublen); 

   sw_mean(a, mua, lena, sublen);
   sw_mean(b, mub, lenb, sublen);

   sw_inv_meancentered_norm(a, mua, invna, lena, sublen);
   sw_inv_meancentered_norm(b, mub, invnb, lenb, sublen);

   dfdg_init(a, mua, dfa, dga, lena, sublen);
   dfdg_init(b, mub, dfb, dgb, lenb, sublen);
   
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

*/
