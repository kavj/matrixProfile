#include <cstdio>
#include <algorithm>
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "../utils/alloc.h"
#include "../utils/primitive_print_funcs.h"
#include "tiled_pearson.h"
#include "pearson.h"
#define klen 128 

static void dfdg_init(const double* __restrict__ ts, const double* __restrict__ mu, double* __restrict__ df, double* __restrict__ dg, long long len, long long sublen){
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
      dg[i + 1] = (ts[i + sublen] - ts[i])/2.0;
   }
}
 
static inline void pauto_pearson_kern(
   double*       cov,
   double*       mp,
   long long*    mpi,
   const double* df,
   const double* dg,
   const double* invn,
   const long long ofr,
   const long long ofc)
{
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);
   df   = (const double*)__builtin_assume_aligned(df, prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (long long*) __builtin_assume_aligned(mpi, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);
   for(long long r = 0; r < klen; r++){
      double cr[klen];
      for(long long d = 0; d < klen; d++){ 
         cov[d] += df[r] * dg[r + d + ofc];
         cov[d] += df[r + d + ofc] * dg[r];
      }
      for(long long d = 0; d < klen; d++){
         cr[d] = cov[d] * invn[r] * invn[r + d + ofc];
      }
      int ci[64];
      for(long long i = 0; i < 64; i++){ // if we write this in the obvious way using conditional statements, gcc 7.3 is unable to generate
                                   // a max reduction via masked writes or blend operations. I'll file a bug report at some point.
         ci[i] = cr[i] > cr[i + 64] ? i : i + 64;
         cr[i] = cr[i] > cr[i + 64] ? cr[i] : cr[i + 64];
      }
      for(int i = 0; i < 32; i++){
         ci[i] = cr[i] > cr[i + 32] ? ci[i] : ci[i + 32];
         cr[i] = cr[i] > cr[i + 32] ? cr[i] : cr[i + 32];
      } 
      for(int i = 0; i < 16; i++){
         ci[i] = cr[i] > cr[i + 16] ? ci[i] : ci[i + 16];
         cr[i] = cr[i] > cr[i + 16] ? cr[i] : cr[i + 16];
      }
      for(int i = 0; i < 8; i++){
         ci[i] = cr[i] > cr[i + 8] ? ci[i] : ci[i + 8];
         cr[i] = cr[i] > cr[i + 8] ? cr[i] : cr[i + 8];
      }
      for(int i = 0; i < 4; i++){
         ci[i] = cr[i] > cr[i + 4] ? ci[i] : ci[i + 4];
         cr[i] = cr[i] > cr[i + 4] ? cr[i] : cr[i + 4];
      }
      for(int i = 0; i < 2; i++){
         ci[i] = cr[i] > cr[i + 2] ? ci[i] : ci[i + 2];
         cr[i] = cr[i] > cr[i + 2] ? cr[i] : cr[i+2];
      }
      ci[0] = cr[0] > cr[1] ? ci[0] : ci[1];
      cr[0] = cr[0] > cr[1] ? cr[0] : cr[1];
      mpi[r] = mp[r] > cr[0] ? mpi[r] : ci[0] + ofr + ofc;
      mp[r] =  mp[r] > cr[0] ? mp[r] : cr[0];
   }
}


static inline void pauto_pearson_edge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*          __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   long long ofr, 
   long long ofc, 
   long long dlim,
   long long clim)
{
   for(int d = 0; d < dlim; d++){
      for(int r = 0; r < clim - d; r++){
         cov[d] += df[r] * dg[r + d + ofc];
         cov[d] += df[r + d + ofc]*dg[r];
         if(cov[d] * invn[r] * invn[r + d + ofc] > mp[r]){
            mp[r] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r] = r + d + ofc + ofr;
         }
         if(cov[d] * invn[r] * invn[r + d + ofc] > mp[r + d + ofc]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         }
      }
   }
}


static inline void pauto_pearson_corner(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   long long ofr, 
   long long ofc, 
   long long dlim,
   long long clim)
{
   for(int d = 0; d < dlim; d++){
      if(cov[0] * invn[0] * invn[d + ofc] > mp[0]){
         mp[0] = cov[0] * invn[0] * invn[d + ofc];
         mpi[0] = d + ofc + ofr;
      }
      if(cov[d] * invn[0] * invn[d + ofc] > mp[d + ofc]){
         mp[d + ofc] = cov[0] * invn[0] * invn[d + ofc];
         mpi[d + ofc] = ofr;
      }
   }
   pauto_pearson_edge(cov, mp + 1, mpi + 1, df + 1, dg + 1, invn + 1, ofr + 1, ofc, dlim, clim - 1);
}


static inline void pauto_pearson_init(
   double* __restrict__ cov, 
   const double* __restrict__ df,
   const double* __restrict__ dg,
   double*       __restrict__ mp,  
   long long*          __restrict__ mpi, 
   const double* __restrict__ invn, 
   long long ofr, 
   long long ofc)
{
   cov  = (double*)__builtin_assume_aligned(cov,prefalign);
   df   = (const double*) __builtin_assume_aligned(df, prefalign);
   dg   = (const double*) __builtin_assume_aligned(dg, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp,prefalign);
   mpi  = (long long*) __builtin_assume_aligned(mpi,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   double cbuf[klen];
   for(long long d = 0; d < klen; d++){
      cbuf[d] = cov[d] * invn[0] * invn[d + ofc];
   }
   for(long long d = 0; d < klen; d++){
      if(mp[0] < cbuf[d]){
         mp[0] = cov[d] * invn[0] * invn[d + ofc];
         mpi[0] = d + ofc + ofr;
      }
   }
   for(long long d = 0; d < klen; d++){
      if(mp[d + ofc] < cbuf[d]){
         mp[d+ofc] = cbuf[d];
         mpi[d+ofc] = ofr;
      }
   }
   for(long long r = 1; r < klen; r++){
      for(int d = 0; d < klen; d++){
         cov[d] += df[r] * dg[r + d + ofc];
         cov[d] += df[r + d + ofc] * dg[r];
      }
      for(long long d = 0; d < klen; d++){
         cbuf[d] = cov[d] * invn[r] * invn[r + d + ofc];
      }
      for(long long d = 0; d < klen; d++){
         if(mp[r] < cbuf[d]){
            mp[r] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r] = r + d + ofc + ofr;
         }
      }
      for(long long d = 0; d < klen; d++){
         if(mp[r + d + ofc] < cbuf[d]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         } 
      }
   }
}


// Todos: This appears to have a resurfaced range bug, but the loop ranges are quite a bit clearer due to alignment checks being restricted to the innermost tile.
//        At this point we never process more than kernel len in either direction in a single loop block,
//        and we determine alignment based only on that condition. This makes it easier to either step through in a 
//        debugger or print the indices per starting position to verify start and end conditions. 
//
//
//        Once that's resolved we should swap in the intrinsics based kernel anytime we're not in the initialization phase.
//        Then we should eventually make it possible to run this on Windows by supporting their aligned allocator and Visual C++ syntax for restrict and force inline
//        under the names noalias and stronginline or something similar. 

int pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // need to write the supporting namespace for this
   }
   const int mlen = ts.len - sublen + 1;
   const int tlen = std::max(2 << 14, 4 * sublen - (4 * sublen) % klen);        
   const int tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i), sublen); 
   }
   const int kal = mlen - (mlen - minlag) % klen;  // kal is the first column index such that column_max - kal < kernel_length
   for(int diag = minlag; diag < mlen; diag += tlen){
      #pragma omp parallel for
      for(int ofst = 0; ofst < mlen - diag; ofst += tlen){
         const int dlim = std::min(tlen, mlen - diag - ofst);
         batchcov(ts(diag + ofst), mu(diag + ofst), q(ofst/tlen), cov(ofst), dlim, sublen);
         const int dalim = std::max(0, std::min(tlen, kal - diag - ofst - klen));
         for(int d = diag; d < diag + dalim; d += klen){
            pauto_pearson_init(cov(d - diag + ofst), df(ofst), dg(ofst), mp(ofst), mpi(ofst), invn(ofst), ofst, d);
            const int ralim = std::min(tlen, kal - d - ofst); 
            //printf("dalim: %d ralim %d\n",dalim, ralim);
            for(int r = ofst + klen; r < ofst + ralim; r += klen){
               pauto_pearson_kern(cov(d - diag + ofst), mp(ofst), mpi(ofst), df(ofst), dg(ofst), invn(ofst), ofst, d);
            }
            if(ralim < tlen){
               //printf("edge: offset: %d, row lim: %d\n", ofst + ralim, mlen - d - ofst - ralim);
               pauto_pearson_edge(cov(d - diag + ofst + ralim), 
                                  mp(ofst + ralim), 
                                  mpi(ofst + ralim), 
                                  df(ofst + ralim), 
                                  dg(ofst + ralim), 
                                  invn(ofst + ralim), 
                                  ofst + ralim, 
                                  d, 
                                  klen, 
                                  mlen - d - ofst - ralim);
            }   
         }
         if(dalim < tlen){
            //printf("corner: diag lim: %d column lim: %d\n", dlim - dalim, mlen - diag - dalim - ofst); 
            pauto_pearson_corner(cov(dalim + ofst), 
                                 mp(ofst), 
                                 mpi(ofst), 
                                 df(ofst), 
                                 dg(ofst), 
                                 invn(ofst), 
                                 ofst, 
                                 diag + dalim, 
                                 dlim - dalim, 
                                 mlen - diag - ofst - dalim); 
         }
      }
   }
   return errs::none;
}

