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

// It is probably a better idea to externally allocate these buffers, since this is generating a ridiculous amount of stack traffic 
// and the current method prevents refactoring of pairwise max pooling. Unfortunately it drives the fused kernel and reference further apart, and it would cause a divergence in function signatures. Until the implementation stabilizes, we will just live with the lack of refactoring.

static void pauto_pearson_simple(
   double*       __restrict__  cov,
   double*       __restrict__  mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const long long ofr,
   const long long ofc,
   const long long iters,
   const bool init)
{
   for(long long r = 0; r < iters; r++){
      //double cr[128];
      if(!init || r != 0){
         for(long long d = 0; d < 128; d++){ 
            cov[d] += df[r] * dg[r + d + ofc];
            cov[d] += df[r + d + ofc] * dg[r];
         }
      }
      for(int d = 0; d < 128; d++){
         double cr = cov[d] * invn[r] * invn[r + d + ofc];
         if((cr > 1) || (cr < -1)){
            printf("simple makes me sad\n");
         }
         if(mp[r] < cr){
            mp[r] = cr;
            mpi[r] = r + d + ofr + ofc;
         }
         if(mp[r + d + ofc] < cr){
            mp[r + d + ofc] = cr;
            mpi[r + d + ofc] = r + ofr;
         }
      }
   }
}

void pauto_pearson_kern(
   double*       __restrict__  cov,
   double*       __restrict__  mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const long long ofr,
   const long long ofc,
   const long long iters,
   const bool init)
{
   df   = (const double*)__builtin_assume_aligned(df, prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (long long*) __builtin_assume_aligned(mpi, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);
   for(long long r = 0; r < iters; r++){
      double cr[128];
      if(!init || r != 0){
         for(long long d = 0; d < 128; d++){ 
            cov[d] += df[r] * dg[r + d + ofc];
            cov[d] += df[r + d + ofc] * dg[r];
         }
      }
      for(long long d = 0; d < 128; d++){
         cr[d] = cov[d] * invn[r] * invn[r + d + ofc];
      }
      long long ci[64];
      for(long long i = 0; i < 64; i++){ // if we write this in the obvious way using conditional statements, gcc 7.3 is unable to generate
                                   // a max reduction via masked writes or blend operations. I'll file a bug report at some point.
         ci[i] = cr[i] > cr[i + 64] ? i : i + 64;
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
      if(cr[0] > 1.0){
         printf("this makes me sad\n");
      }
      mpi[r] = mp[r] > cr[0] ? mpi[r] : ci[0] + ofr + ofc;
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
   const long long ofc, 
   const long long dlim,
   const long long clim,
   const bool init)
{
   for(long long d = 0; d < dlim; d++){
      for(long long r = 0; r < clim - d; r++){
         if(!init || r != 0){
            cov[d] += df[r] * dg[r + d + ofc];
            cov[d] += df[r + d + ofc] * dg[r];
         }
         if(cov[d] * invn[r] * invn[r + d + ofc] > mp[r]){
            mp[r] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r] = r + d + ofc + ofr;
         }
         if(cov[d] * invn[r] * invn[r + d + ofc] > mp[r + d + ofc]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = r + ofr;
         }
         if(cov[d] * invn[r] * invn[r + d + ofc] > 1.0){
            printf("This edge makes me sad, row: %d col: %d cov: %lf \n", r + ofr, r + d + ofr + ofc, cov[d]);
         }
      }
   }
}



void pearson_pauto_reference_solve(double*       __restrict__ cov, 
                                   double*       __restrict__ mp, 
                                   long long*    __restrict__ mpi, 
                                   const double* __restrict__ df, 
                                   const double* __restrict__ dg, 
                                   const double* __restrict__ invn,
                                   const long long minlag, 
                                   const long long mlen, 
                                   const long long sublen){
   if((cov == nullptr) || (mp == nullptr) || (mpi == nullptr) || 
      (df == nullptr)  || (dg == nullptr) || (invn == nullptr)){
      printf("problem in solver input\n");
      exit(1);
   }

   for(int diag = minlag; diag < mlen; diag++){
      double c = cov[diag - minlag];
      if(mp[0] < c * invn[0] * invn[diag]){
         mp[0] = c * invn[0] * invn[diag];
         mpi[0] = diag;
      }
      if(mp[diag] < c * invn[0] * invn[diag]){
         mp[diag] = c * invn[0] * invn[diag];
         mpi[diag] = 0;
      }
      for(int ofst = 0; ofst < mlen - diag; ofst++){
         c += df[ofst] * dg[diag + ofst];
         c += df[diag + ofst] * dg[ofst];
         if(mp[ofst] < (c * invn[ofst] * invn[diag + ofst])){
            mp[ofst] = c * invn[ofst] * invn[diag + ofst];
            mpi[ofst] = diag + ofst;
         }
         if(mp[diag + ofst] < c * invn[ofst] * invn[diag + ofst]){
            mp[diag + ofst] = c * invn[ofst] * invn[diag + ofst];
            mpi[diag + ofst] = ofst;
         }
      }
   }
}

int pearson_pauto_reduc_ref(dsbuf& ts, dsbuf& mp, lsbuf& mpi, long long minlag, long long sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // need to write the supporting namespace for this
   }
   const long long mlen = ts.len - sublen + 1;
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(1, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   center_query(ts(0), mu(0), q(0), sublen); 
   batchcov(ts(minlag), mu(minlag), q(0), cov(0), mlen - minlag, sublen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/q",q(0),sublen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/cov",cov(0),mlen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/mu",mu(0),mlen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/df",df(0),mlen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/dg",dg(0),mlen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/invn",invn(0),mlen);
   
   pearson_pauto_reference_solve(cov(0), mp(0), mpi(0), df(0), dg(0), invn(0), minlag, mlen, sublen);  
   writeLongs("/home/kkamg001/matlabscripts/cppoutput/mpi",mpi(0),mlen);
   writeDoubles("/home/kkamg001/matlabscripts/cppoutput/mp",mp(0),mlen);
 
   return errs::none;
}
  

int pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, long long minlag, long long sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   const long long mlen = ts.len - sublen + 1;
   const long long tlen = std::max(static_cast<long long>(2 << 14), 4 * sublen - (4 * sublen) % klen);        
   const long long tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   dfdg_init(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   #pragma omp parallel for
   for(long long i = 0; i < mlen; i+= tlen){
      center_query(ts(i), mu(i), q(i/tlen), sublen); 
   }
   for(long long diag = minlag; diag < mlen; diag += tlen){
      #pragma omp parallel for
      for(long long ofst = 0; ofst < mlen - diag; ofst += tlen){
         const long long dlim = std::min(tlen, mlen - diag - ofst);
         batchcov(ts(diag + ofst), mu(diag + ofst), q(ofst/tlen), cov(ofst), dlim, sublen);
         for(long long d = diag; d < dlim; d += klen){
            const long long rlim = std::min(tlen, mlen - d - ofst);
            const long long rfringe = rlim % klen;
            pauto_pearson_simple(cov(ofst + d - diag), 
                               mp(ofst), 
                               mpi(ofst), 
                               df(ofst), 
                               dg(ofst), 
                               invn(ofst), 
                               ofst, 
                               d, 
                               rlim - rfringe,
                               true);
            if(rfringe > 0){
               pauto_pearson_edge(cov(ofst + d - diag), 
                                  mp(ofst), 
                                  mpi(ofst), 
                                  df(ofst), 
                                  dg(ofst), 
                                  invn(ofst), 
                                  ofst, 
                                  diag + d, 
                                  klen, 
                                  rfringe,
                                  rlim - rfringe - ofst == 0);
            }
         }
      }
   }
   return errs::none;
}

