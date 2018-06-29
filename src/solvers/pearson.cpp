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


auto pauto_pearson_kern = [&](
   double*       cov,
   double*       mp,
   int*          mpi,
   const double* df,
   const double* dg,
   const double* invn,
   const int ofr,
   const int ofc)
{
   for(int r = 0; r < klen; r++){
      for(int d = 0; d < klen; d++){ 
         if(r > 0){
            cov[d] += df[r] * dg[r + d + ofc];
            cov[d] += df[r + d + ofc] * dg[r];
         }
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


auto pauto_pearson_edge = [&](
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


auto pauto_pearson_init = [&](
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
};


void pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      printf("bad inputs\n");
      exit(1);
   }
   const int mlen = ts.len - sublen + 1;
   const int tlen = std::max(16384, 4 * sublen - (4 * sublen) % klen);  
   const int tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      printf("error allocating memory\n");
      exit(1);  // Todo: Exceptions should be thrown on memory allocation issues, or we could add a return value
   }
   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   init_dfdg(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i), sublen); 
   }
   const int fringe = (mlen - minlag) % klen;
   for(int diag = minlag; diag < mlen; diag += tlen){
      #pragma omp parallel for
      for(int ofst = 0; ofst < mlen - diag; ofst += tlen){
         const int dmx = std::min(diag + tlen, mlen - ofst - fringe);
         batchcov(ts(diag + ofst), mu(diag + ofst), q(ofst/tlen), cov(ofst), dmx - diag, sublen);
         pauto_pearson_init(cov(ofst), mp(ofst), mpi(ofst), invn(ofst), ofst, diag, dmx - diag);
         for(int d = diag; d < dmx; d += klen){
            const int rmx = std::min(ofst + tlen, mlen - fringe - d);
            //Todo: need a truncated kernel here          
            for(int r = ofst + klen; r < rmx; r += klen){
               pauto_pearson_kern(cov(ofst + d - diag), mp(r), mpi(r), df(r), dg(r), invn(r), r, d);
            }
         }
         if(diag + ofst + 2 * tlen > mlen){
            // tile fringe 
         }
      }
   }
}


