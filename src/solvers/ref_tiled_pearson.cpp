#include<algorithm>
#include "../utils/reg.h"
#include "../utils/xprec.h"
#include "../utils/reduce.h"
#define prefalign 64  // Todo: put this somewhere in build settings, since its used by both allocators and anything that may apply auto vectorization

auto pauto_pearson_init_naive = [&](
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
   const double* __restrict__ invn, 
   int ofr, 
   int ofc, 
   int dlim)
{
   for(int d = 0; d < dlim; d++){
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

auto pauto_pearson_update_naive = [&](
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
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

/*
static void pauto_pearson_edge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int tlen,
   int ofr, 
   int ofc, 
   int clim)
{
   cov =  (double*)      __builtin_assume_aligned(cov,prefalign);
   mp =   (double*)      __builtin_assume_aligned(mp,prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   const int tail = clim%klen;
   const int dalgn = std::max(0,std::min(clim - tail - 2*klen, tlen));
   for(int d = 0; d < dalgn; d+= klen){
      const int ralgn = std::min(dalgn + klen - d, tlen);
      // optimized init
      for(int r = klen; r < ralgn; r += klen){
         pauto_pearson_AVX_kern(cov+d,mp+r,mpi+r,df+r,dg+r,invn+r,ofr+r,ofc+d);
      } 
   }
   // We run the fringe sections after doing everything where simd is possible
   // This is particularly import for AVX where the compiler may not call zero_upper for every scalar section
   if(tail != 0){
      for(int d = 0; d < dalgn; d+= klen){
         const int ral = std::max(0,std::min(dalgn + klen - d, tlen));
         if(ral < tlen){
            pauto_pearson_update_naive(cov+d,mp+ral,mpi+ral,df+ral,dg+ral,invn+ral,ofr,ofc+dalgn,klen,clim-d);
         }
      }
   }
   if(clim < tlen){
      pauto_pearson_init_naive(cov+dalgn,mp,mpi,invn,ofr,ofc+dalgn,clim-dalgn);
      pauto_pearson_update_naive(cov+dalgn,mp+1,mpi+1,df+1,dg+1,invn+1,ofr+1,ofc+dalgn,clim-dalgn-1,clim-dalgn-1); 
   }
}*/



void pauto_pearson_edge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int tlen,
   int ofr, 
   int ofc, 
   int bound)
{
   cov =  (double*)      __builtin_assume_aligned(cov, prefalign);
   mp =   (double*)      __builtin_assume_aligned(mp,  prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi, prefalign);
   df =   (const double*)__builtin_assume_aligned(df,  prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,  prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);
 
   for(int d = 0; d < std::min(tlen,bound); d++){
      if(cov[d] * invn[0] * invn[d + ofc] > mp[0]){
         mp[0] = cov[d] * invn[0] * invn[d + ofc];  
         mpi[0] = d + ofr + ofc;
      }
      if(cov[d] * invn[0] * invn[d + ofc] > mp[d + ofc]){
         mp[d + ofc] = cov[d] * invn[0] * invn[d + ofc]; 
         mpi[d + ofc] = ofr;
      }
      for(int r = 1; r < std::min(tlen, bound - d); r++){
         cov[d] += df[r] * dg[r + d + ofc];
         cov[d] += df[r + d + ofc] * dg[r];
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

/*
void pauto_pearson_inner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int ofr,
   const int ofc)
{
   cov  = (double*)      __builtin_assume_aligned(cov,prefalign);
   mp   = (double*)      __builtin_assume_aligned(mp,prefalign);
   mpi  = (long long*)   __builtin_assume_aligned(mpi,prefalign);
   df   = (const double*)__builtin_assume_aligned(df,prefalign);
   dg   = (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int d = 0; d < tlen; d += klen){
      // call init pass here
      for(int r = klen; r < tlen; r += klen){
         pauto_pearson_AVX_kern(cov+d,mp+r,mpi+r,df+r,dg+r,invn+r,ofr+r,ofc+d);
      }
   }
}*/


void pauto_pearson_inner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int ofr,
   const int ofc)
{
   cov =  (double*)      __builtin_assume_aligned(cov, prefalign);
   mp =   (double*)      __builtin_assume_aligned(mp,  prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi, prefalign);
   df =   (const double*)__builtin_assume_aligned(df,  prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,  prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int d = 0; d < tlen; d++){
      if(mp[0] < cov[d] * invn[0] * invn[d + ofc]){
         mp[0] = cov[d] * invn[0] * invn[d + ofc];
         mpi[0] = static_cast<long long>(d + ofr + ofc);
      }
      if(mp[d + ofc] < cov[d] * invn[0] * invn[d + ofc]){
         mp[d + ofc] = cov[d] * invn[0] * invn[d + ofc];
         mpi[d + ofc] = ofr;
      }
   }
   for(int r = 1; r < tlen; r++){
      for(int d = 0; d < tlen; d++){ 
         cov[d] += df[r] * dg[r + d + ofc];
         cov[d] += df[r + d + ofc] * dg[r];
         if(mp[r] < cov[d] * invn[r] * invn[r + d + ofc]){
            mp[r] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r] = static_cast<long long>(r + d + ofr + ofc);
	 }
	 if(mp[r + d + ofc] < cov[d] * invn[r] * invn[r + d + ofc]){
            mp[r + d + ofc] = cov[d] * invn[r] * invn[r + d + ofc];
            mpi[r + d + ofc] = static_cast<long long>(r + ofr);
         }
      }
   }
}

