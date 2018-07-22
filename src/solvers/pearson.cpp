#include <cstdio>
#include <algorithm>
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "../utils/alloc.h"
#include "../utils/primitive_print_funcs.h"
#include "tiled_pearson.h"
#include "pearson.h"
#define klen 64

static void dfdg_init(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ df, dtype* __restrict__ dg, int len, int sublen){
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
   int*          mpi,
   const double* df,
   const double* dg,
   const double* invn,
   const int ofr,
   const int ofc)
{
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);
   df   = (const double*)__builtin_assume_aligned(df, prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (int*) __builtin_assume_aligned(mpi, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);
   
   for(int r = 0; r < klen; r++){
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
}


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
   double c = *cov;
   for(int d = 0; d < dlim; d++){
      for(int r = 0; r < clim - d; r++){
         c += df[r] * dg[r + d + ofc];
         c += df[r + d + ofc]*dg[r];
         if(c * invn[r] * invn[r + d + ofc] > mp[r]){
            mp[r] = c * invn[r] * invn[r + d + ofc];
            mpi[r] = r + d + ofc + ofr;
         }
         if(c * invn[r] * invn[r + d + ofc] > mp[r + d + ofc]){
            mp[r + d + ofc] = c * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         }
      }
   }
   *cov = c;
}


static inline void pauto_pearson_corner(
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
      if(cov[0] * invn[0] * invn[d + ofc] > mp[r]){
         mp[0] = cov[0] * invn[0] * invn[d + ofc];
         mpi[0] = d + ofc + ofr;
      }
      if(cov[d] * invn[0] * invn[d + ofc] > mp[d + ofc]){
         mpr[d + ofc] = cov[0] * invn[0] * invn[d + ofc];
         mpi[d + ofc] = ofr;
      }
   }
   pauto_pearson_edge(cov, mp + 1, mpi + 1, df + 1, dg + 1, invn + 1, ofr + 1, ofc, dlim, clim - 1);
}


// changing this to be aligned init only
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

// Note to self: The std::min expressions can be used to organize various boundary conditions in a way that they can be quickly verified by hand
// Later on it would probably be more readable if they're factored out into temporary variables. 


// Note: The prior function was extremely difficult to test. Since the fixed sized kernels run klen * klen updates, we can only worry about alignment here.
// This means if diagonal + offset + 2 * klen <= mlen, then we can process a full kernel. The extra statements just handle two other cases, the one where we have an initialized block
// where some updates may be unaligned and the case where we can't even initialize a full section. The former processes a block of size klen in one direction and iterates over the other to either a tile boundary
// or an edge column. The other basically folds in the initialization step and processes however many diagonals remain then runs the same basic process to the last column.

// This was really really annoying to work out, but it's much easier to verify boundary conditions now due to the lack of obfuscation.

int pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return -1;
      //return errs::bad_input;  // need to write the supporting namespace for this
   }
   const int mlen = ts.len - sublen + 1;
   const int tlen = std::max(2 << 15, 4 * sublen - (4 * sublen) % klen);        
   const int tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
       return -1;
      //return errs::mem_error;
   }
   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i), sublen); 
   }
   const int kal = mlen - (mlen - minlag) % klen;  // Column ranges are from minlag to mlen. Since minlag may not be a multiple of our preferred simd alignment, we align by row range 0 to mlen - diag
                                                   // so our boundary cases for aligned sections are based on (mlen - minlag) modulo tile length rather than mlen modulo tile length
   for(int diag = minlag; diag < mlen; diag += tlen){
      #pragma omp parallel for
      for(int ofst = 0; ofst < mlen - diag; ofst += tlen){
         const int dlim = std::min(tlen, mlen - diag - ofst);
         batchcov(ts(diag + ofst), mu(diag + ofst), q(ofst/tlen), cov(ofst), dlim, sublen);
         const int dalim = std::min(tlen, kal - diag - ofst - klen);
         for(int d = diag; d < diag + dalim; d += klen){
            pauto_pearson_init(cov(d - diag + ofst), mp(ofst), mpi(ofst), invn(ofst), ofst, d, klen);
            const int ralim = std::min(tlen, kal - d - ofst); 
            for(int r = ofst + klen; r < ofst + ralim; r += klen){
               
            }
            if(d + ofst + klen + tlen > mlen){
               // This shouldn't be inline, but these are the appropriate semantics
               // should go from column == kal to min(mlen, row + tlen)
               // that's why kal - d is a starting point for row val
               for(int sd = d; sd < d + klen; sd++){
                  const int re = std::min(tlen, mlen - sd - ofst);
                  for(int r = kal - d; r < re; r++){
                     
                  }
               }
            }
         }
         for(int d = dalim; d < dlim; d++){ 
            // corners, basically single diagonal 

         }
      }
   }
   return 0;
   //return errs::none;
}

