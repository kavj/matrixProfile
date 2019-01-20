#include <iostream>
#include <cmath>
#include "../utils/alloc.h"
#include <omp.h>
#include <array>
#include "pearson.h"
#include "../utils/moments.h"

using namespace pearson;

constexpr int klen  = 256; 
struct argmax{
   double val;
   int ind;
};


void pearson::tonormalizedeuclidean(double* mp, int len, int sublen){
   #pragma omp parallel for
   for(int i = 0; i < len - sublen + 1; i++){
      mp[i] = sqrt(2 * sublen * (1 - mp[i]));
   }
}

static void dfdg_init(const double* __restrict ts, const double* __restrict mu, double* __restrict df, double* __restrict dg, int len, int sublen){
   df[0] = 0;
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - ts[i])/2.0;
      dg[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
   }
}

static inline struct argmax pw_reduc(std::array<double, klen>& cr){
   std::array<int, klen/2> ci;
   for(int i = 0; i < 128; i++){
      ci[i] = cr[i] > cr[i + 128] ? i : i + 128;
      cr[i] = cr[i] > cr[i + 128] ? cr[i] : cr[i + 128];
   }
   for(int i = 0; i < 64; i++){ 
      ci[i] = cr[i] > cr[i + 64] ? ci[i] : ci[i + 64];
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
   struct argmax r = {cr[0], ci[0]};
   return r;
}


static inline void partialcross_kern(
   double*       __restrict cov,
   double*       __restrict mp,
   int*          __restrict mpi,
   const double* __restrict dfa,
   const double* __restrict dga,
   const double* __restrict invna,
   const double* __restrict dfb,
   const double* __restrict dgb,
   const double* __restrict invnb,
   const int amx, const int bmx)
{
   dfa   = (const double*)__builtin_assume_aligned(dfa, prefalign);
   dga   = (const double*)__builtin_assume_aligned(dga, prefalign);
   invna = (const double*)__builtin_assume_aligned(invna, prefalign);
   dfb   = (const double*)__builtin_assume_aligned(dfb, prefalign);
   dgb   = (const double*)__builtin_assume_aligned(dgb, prefalign);
   invnb = (const double*)__builtin_assume_aligned(invnb, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (int*) __builtin_assume_aligned(mpi, prefalign);
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);

   for(int ia = 0; ia < amx; ia++){
      double cr = cov[ia] * invna[ia] * invnb[0];
      if(cr > mp[ia]){
         mp[ia] = cr;
         mpi[ia] = ia;
      }
   }
   for(int ia = 0; ia < amx; ia++){
      int mx = std::min(amx - ia, bmx);
      for(int ib = 1; ib < mx; ib++){
         cov[ia] += dfa[ib + ia] * dgb[ib];
         cov[ia] += dfb[ib] * dga[ib + ia];
         double cr = cov[ia] * invna[ib + ia] * invnb[ib];
         if(cr > mp[ib + ia]){
            mp[ib + ia] = cr;
            mpi[ib + ia] = ib;
         }
      }
   }
}


static inline void partialcross_transkern(
   double*       __restrict cov,
   double*       __restrict mp,
   int*          __restrict mpi,
   const double* __restrict dfa,
   const double* __restrict dga,
   const double* __restrict invna,
   const double* __restrict dfb,
   const double* __restrict dgb,
   const double* __restrict invnb,
   const int amx, const int bmx)
{
   dfa   = (const double*)__builtin_assume_aligned(dfa, prefalign);
   dga   = (const double*)__builtin_assume_aligned(dga, prefalign);
   invna = (const double*)__builtin_assume_aligned(invna, prefalign);
   dfb   = (const double*)__builtin_assume_aligned(dfb, prefalign);
   dgb   = (const double*)__builtin_assume_aligned(dgb, prefalign);
   invnb = (const double*)__builtin_assume_aligned(invnb, prefalign);
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (int*) __builtin_assume_aligned(mpi, prefalign);
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);

   for(int ib = 0; ib < bmx; ib++){
      double cr = cov[ib] * invna[0] * invnb[ib];
      if(cr > mp[0]){
         mp[0] = cr;
         mpi[0] = ib;
      }
   }
   for(int ib = 0; ib < bmx; ib++){
      int mx = std::min(bmx - ib + 1, amx);
      for(int ia = 0; ia < mx; ia++){
         cov[ib] += dfa[ia] * dgb[ia + ib];
         cov[ib] += dfb[ia + ib] * dga[ia];
         double cr = cov[ib] * invna[ia] * invnb[ia + ib];
         if(cr > mp[ia]){
            mp[ia] = cr;
            mpi[ia] = ia + ib;
         }
      }
   }
}


int pearson::partialcross(bufd& a, bufd& b, bufd& mp, bufi& mpi, int sublen){
   if(!(a.valid() && b.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   const int qstride = paddedlen(sublen, prefalign);
   const int amx = a.len - sublen + 1;
   const int bmx = b.len - sublen + 1;
   const int tlen = std::max(static_cast<int>(2 << 14), 4 * sublen - (4 * sublen) % klen);        
   const int tilesperdima = amx/tlen +  amx % tlen ? 1 : 0;
   const int tilesperdimb = bmx/tlen + bmx % tlen ? 1 : 0;
   bufd mua(amx); bufd invna(amx); bufd dfa(amx); bufd dga(amx); bufd mub(bmx); bufd invnb(bmx);
   bufd dfb(bmx); bufd dgb(bmx);   bufd cov(bmx); bufd q(std::max(tilesperdima, tilesperdimb) * qstride);
   if(!(mua.valid() && invna.valid() && dfa.valid() && dga.valid() && mub.valid() && invnb.valid()) && invnb.valid() && dfb.valid() && dgb.valid() && cov.valid()){
      return errs::mem_error;
   }
   sw_mean(a(), mua(), a.len, sublen);
   sw_mean(b(), mub(), b.len, sublen);
   sw_inv_meancentered_norm(a(), mua(), invna(), a.len, sublen);
   sw_inv_meancentered_norm(b(), mub(), invnb(), b.len, sublen);
   dfdg_init(a(), mua(), dfa(), dga(), a.len, sublen);
   dfdg_init(b(), mub(), dfb(), dgb(), b.len, sublen);

   #pragma omp parallel for
  /* for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i * qstride), sublen); 
   }*/
   
   center_query(b(), mub(), q(), sublen); 
   crosscov(a(), mua(), q(), cov(), amx, sublen); 
   partialcross_kern(cov(), mp(), mpi(), dfa(), dga(), invna(), dfb(), dgb(), invnb(), amx, bmx); 
   center_query(a(), mua(), q(), sublen);
   crosscov(b(), mub(), q(), cov(), bmx, sublen);
   partialcross_kern(cov(), mp(), mpi(), dfa(), dga(), invna(), dfb(), dgb(), invnb(), amx, bmx); 

   return pearson::errs::none;
}



static inline void partialauto_kern(
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
   mp   = (double*) __builtin_assume_aligned(mp, prefalign);
   mpi  = (int*) __builtin_assume_aligned(mpi, prefalign);
   invn = (const double*)__builtin_assume_aligned(invn, prefalign);
   cov  = (double*)__builtin_assume_aligned(cov, prefalign);
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
      struct argmax p = pw_reduc(cr);
      if(p.val > mp[r]){
         mp[r] = p.val;
         mpi[r] = p.ind + r + ofr + ofc;
      }
   }
}

static inline void partialauto_edge(
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



int pearson::partialauto(bufd& ts, bufd& mp, bufi& mpi, int minlag, int sublen){
   if(!(ts.valid() && mp.valid() && mpi.valid())){
      return errs::bad_input;  // Todo: Build a real set of error checking functions 
   }
   const int qstride = paddedlen(sublen, prefalign);
   const int mlen = ts.len - sublen + 1;
   const int tlen = std::max(static_cast<int>(2 << 14), 4 * sublen - (4 * sublen) % klen);        
   const int tilesperdim = (mlen - minlag)/tlen + ((mlen - minlag) % tlen ? 1 : 0);
   bufd mu(mlen); bufd invn(mlen); bufd df(mlen);  
   bufd dg(mlen); bufd cov(mlen);  bufd q(tilesperdim * qstride);
   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      return errs::mem_error;
   }
   
   sw_mean(ts(), mu(), ts.len, sublen);
   sw_inv_meancentered_norm(ts(), mu(), invn(), ts.len, sublen);

   dfdg_init(ts(), mu(), df(), dg(), ts.len, sublen);
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i * tlen), mu(i * tlen), q(i * qstride), sublen); 
   }
   for(int diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for 
      for(int ofst = 0; ofst < tilesperdim - diag; ofst++){
         const int di = diag * tlen + minlag;
         const int ofi = ofst * tlen;
         const int dlim = std::min(di + tlen, mlen - ofi);
         crosscov(ts(di + ofi), mu(di + ofi), q(ofst * qstride), cov(ofi), dlim - di, sublen);
         for(int d = di; d < dlim; d += klen){
            const int aligned_rlim = std::min(tlen, mlen - d - ofi - klen + 1); 
            partialauto_kern(cov(ofi + d - di), mp(ofi), mpi(ofi), df(ofi), dg(ofi), invn(ofi), ofi, d, aligned_rlim);
            if(aligned_rlim > 0 && aligned_rlim < tlen){
               partialauto_edge(cov(ofi + d - di), mp(), mpi(), df(), dg(), invn(), ofi + aligned_rlim, d, d + klen, mlen, false);
            }
            else if(aligned_rlim == 0){
               partialauto_edge(cov(ofi + d - di), mp(), mpi(), df(), dg(), invn(), ofi, d, dlim, mlen, true);
            }
         }
      }
   }
   return errs::none;
}


