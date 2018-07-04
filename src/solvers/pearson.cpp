#include <cstdio>
#include <algorithm>
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "../utils/alloc.h"
#include "../utils/primitive_print_funcs.h"
#include "tiled_pearson.h"
#include "pearson.h"
#define klen 64

static void init_dfdg(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ df, dtype* __restrict__ dg, int len, int sublen){
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
      dg[i + 1] = (ts[i + sublen] - ts[i])/2.0;
   }
}


template<int ofs>
static inline void pauto_pearson_kern(
   double*       cov,
   double*       mp,
   int*          mpi,
   const double* df,
   const double* dg,
   const double* invn,
   const int ofr,
   const int ofc)
{
   cov  = (double*)__builtin_assume_aligned(cov, prefalign, ofs);
   df   = (const double*)__builtin_assume_aligned(df, prefalign, ofs);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign, ofs);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign, ofs);
   mpi  = (int*) __builtin_assume_aligned(mpi, prefalign, ofs);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign, ofs);
   
   for(int r = ofs; r < klen; r++){
      for(int d = 0; d < klen; d++){ 
         cov[d] += df[r] * dg[r + d + ofc];
         cov[d] += df[r + d + ofc] * dg[r];
         if(mp[r] < cov[d] * invn[r] * invn[r + d + ofc]){
            mp[r] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r] = r + d + ofr + ofc;
	 }
	 if(mp[r + d + ofc] < cov[d] * invn[r] * invn[r + d + ofc]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         }
      }
   }
};


auto pauto_pearson_init_edge = [&](
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   int*          __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int ofr, 
   int ofc, 
   int dlim,
   int clim)
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
            mp[r+d+ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r+d+ofc] = r + ofr;
         }
      }
   }
};


static inline void pauto_pearson_edge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   int*          __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int ofr, 
   int ofc, 
   int dlim,
   int clim)
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
            mp[r+d+ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r+d+ofc] = r + ofr;
         }
      }
   }
}


static inline void pauto_pearson_init(
   const double* __restrict__ cov, 
   double*       __restrict__ mp,  
   int*          __restrict__ mpi, 
   const double* __restrict__ invn, 
   int ofr, 
   int ofc, 
   int dlim)
{
   cov  = (const double*)__builtin_assume_aligned(cov,prefalign);
   mp   = (double*) __builtin_assume_aligned(mp,prefalign);
   mpi  = (int*) __builtin_assume_aligned(mpi,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   int fringe = dlim % klen;
   for(int d = 0; d < dlim - fringe; d += klen){
      for(int sd = d; sd < d + klen; sd++){
         if(cov[d] * invn[0] * invn[d + ofc] > mp[0]){
            mp[0] = cov[d] * invn[0] * invn[d + ofc];
            mpi[0] = d + ofc + ofr;
         }
         if(cov[d] * invn[0] * invn[d + ofc] > mp[d + ofc]){
            mp[d+ofc] = cov[d] * invn[0] * invn[d + ofc];
            mpi[d+ofc] = ofr;
         }
      }
   }
   for(int d = dlim - fringe; d < dlim; d++){
      if(cov[d] * invn[0] * invn[d + ofc] > mp[0]){
         mp[0] = cov[d] * invn[0] * invn[d + ofc];
         mpi[0] = d + ofc + ofr;
      }
      if(cov[d] * invn[0] * invn[d + ofc] > mp[d + ofc]){
         mp[d+ofc] = cov[d] * invn[0] * invn[d + ofc];
         mpi[d+ofc] = ofr;
      }
   }
}


int pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_inputs;
   }
   const int mlen = ts.len - sublen + 1;
   const int tlen = std::max(16384, 4 * sublen - (4 * sublen) % klen);  
   const int tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   init_dfdg(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i), sublen); 
   }
   const int kal = mlen - minlag - (mlen - minlag) % klen;
   for(int diag = minlag; diag < mlen; diag += tlen){
      #pragma omp parallel for
      for(int ofst = 0; ofst < mlen - diag; ofst += tlen){
         const int dlim = std::min(diag + tlen, mlen - ofst);
         batchcov(ts(diag + ofst), mu(diag + ofst), q(ofst/tlen), cov(ofst), dlim - diag, sublen);
         for(int d = diag; d < dlim; d += klen){
            const int rlim = std::min(ofst + tlen, mlen - diag);
            const int ral = std::min(rlim, kal - d);
            pauto_pearson_init(cov(ofst + d - diag), mp(ofst), mpi(ofst), invn(ofst), ofst, d, dlim - diag);
            if(ofst + klen < ral){
               // first section
               for(int r = ofst + klen; ofst < ral; r += klen){
                  pauto_pearson_kern<0>(cov(ofst + d - diag), mp(r), mpi(r), df(r), dg(r), invn(r), r, d);
               }
            }
            if(ral < rlim){
               pauto_pearson_edge(cov(ofst + d - diag), mp(ral), mpi(ral), df(ral), dg(ral), invn(ral), ral, d, dlim, mlen);
            }
         }
      }
   }
   return errs::none;
}


