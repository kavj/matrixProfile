#include<algorithm>
//#include "../utils/max_reduce.h"
#include "../utils/reg.h"
#include "../utils/cov.h"
#include "../utils/xprec_math.h"
#include "../utils/max_reduce.h"
#include "descriptors.h"
#include "../utils/primitive_print_funcs.h"
#ifdef prefalign
#undef prefalign
#endif
#define prefalign 64 
#define unroll 8
//#define prefalign 32 
#define klen 64 
#define simlen 4  // move to header later
#define unroll 8
#define step 32

// rename and reorganize later. This should be in a different namespace or compilation unit



static inline void  pauto_pearson_AVX_kern (double*       __restrict__ cov, 
                                            double*       __restrict__ mp, 
                                            long long*    __restrict__ mpi, 
                                            const double* __restrict__ df, 
                                            const double* __restrict__ dg, 
                                            const double* __restrict__ invn, 
                                            int ofr, 
                                            int ofc){

   cov = (double*)__builtin_assume_aligned(cov,prefalign);
   df  = (const double*)__builtin_assume_aligned(df,prefalign);
   dg  = (const double*)__builtin_assume_aligned(dg,prefalign);
   invn   = (const double*)__builtin_assume_aligned(invn,prefalign);
   mp  = (double*)__builtin_assume_aligned(mp,prefalign);
   mpi = (long long*)__builtin_assume_aligned(mpi,prefalign);

   for(int r = 0; r < klen; r++){
      for(int d = 0; d < klen; d+=step){
         block<__m256d> cov_r;
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = aload(cov,simlen*(d+sd));
         }
         __m256d q = brdcst(dg,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = mul_add(q,uload(df,r+d+simlen*sd),cov_r(sd));
         }
         q = brdcst(df,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = mul_add(q,uload(dg,r+d+simlen*sd),cov_r(sd));
            astore(cov_r(sd),cov,simlen*(d+sd));
         }
         q = brdcst(invn,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) *= q;
         }
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) *= uload(invn,r+d+simlen*sd);
         }
         block<__m256i> mask;
         for(int sd = 0; sd < unroll; sd++){
            mask(sd) = cov_r(sd) > uload(mp,r+d+simlen*sd);
         }
         __m256i s = brdcst(r);
         for(int sd = 0; sd < unroll; sd++){
            if(testnz(mask(sd))){
               maskstore(cov_r(sd),mask(sd),mp+r+d+simlen*sd);
               maskstore(s,mask(sd),mpi+r+simlen*sd);
            }
         }
         struct rpair v = max_reduce_8x1(cov_r(0),cov_r(1),cov_r(2),cov_r(3),cov_r(4),cov_r(5),cov_r(6),cov_r(7));
         //Todo: technically this is incorrect, as we're writing back 4 operands rather than 1. It would be possible to do a final comparison against a temporary buffer and make one sweep at the end
         // alternatively use the if statement to determine if we need to do a horizontal max (since that is quite expensive)
         __m256i msk = v.val > brdcst(mp,r+d);
         if(testnz(msk)){
            maskstore(v.val,msk,mp+r+d);
            maskstore(v.index+brdcst(r+d+ofc),msk,mpi+r+d);
         }
      }
   }
}


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
      if(cov[d]*invn[0]*invn[d+ofc] > mp[0]){
         mp[0] = cov[d]*invn[0]*invn[d+ofc];
         mpi[0] = d+ofc+ofr;
      }
      if(cov[d]*invn[0]*invn[d+ofc] > mp[d+ofc]){
         mp[d+ofc] = cov[d]*invn[0]*invn[d+ofc];
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
      for(int r = 0; r < clim-d; r++){
         cov[d] += df[r]*dg[r+d+ofc];
         cov[d] += df[r+d+ofc]*dg[r];
         if(cov[d]*invn[r]*invn[r+d+ofc] > mp[r]){
            mp[r] = cov[d]*invn[r]*invn[r+d+ofc];
            mpi[r] = r+d+ofc+ofr;
         }
         if(cov[d]*invn[r]*invn[r+d+ofc] > mp[r+d+ofc]){
            mp[r+d+ofc] = cov[d]*invn[r]*invn[r+d+ofc];
            mpi[r+d+ofc] = r+ofr;
         }
      }
   }
};

void pauto_pearson_xedge(
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
}



void pauto_pearson_naive_edge(
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
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mp =   (double*)__builtin_assume_aligned(mp,prefalign);
   mpi =  (long long*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);
 
   for(int d = 0; d < std::min(tlen,bound); d++){
      if(cov[d]*invn[0]*invn[d+ofc] > mp[0]){
         mp[0] = cov[d]*invn[0]*invn[d+ofc];  
         mpi[0] = d+ofr+ofc;
      }
      if(cov[d]*invn[0]*invn[d+ofc] > mp[d+ofc]){
         mp[d+ofc] = cov[d]*invn[0]*invn[d+ofc]; 
         mpi[d+ofc] = ofr;
      }
      for(int r = 1; r < std::min(tlen,bound-d); r++){
         cov[d] += df[r]*dg[r+d+ofc];
         cov[d] += df[r+d+ofc]*dg[r];
         if(cov[d]*invn[r]*invn[r+d+ofc] > mp[r]){
            mp[r] = cov[d]*invn[r]*invn[r+d+ofc];
            mpi[r] = r+d+ofc+ofr;
         }
         if(cov[d]*invn[r]*invn[r+d+ofc] > mp[r+d+ofc]){
            mp[r+d+ofc] = cov[d]*invn[r]*invn[r+d+ofc];
            mpi[r+d+ofc] = r+ofr;
         }
      }
   }
}


void pauto_pearson_xinner(
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
}


void pauto_pearson_naive_inner(
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
   cov =  (double*)      __builtin_assume_aligned(cov,prefalign);
   mp =   (double*)      __builtin_assume_aligned(mp,prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int d = 0; d < tlen; d++){
      if(mp[0] < cov[d]*invn[0]*invn[d+ofc]){
         mp[0] = cov[d]*invn[0]*invn[d+ofc];
         mpi[0] = d+ofr+ofc;
      }
      if(mp[d+ofc] < cov[d]*invn[0]*invn[d+ofc]){
         mp[d+ofc] = cov[d]*invn[0]*invn[d+ofc];
         mpi[d+ofc] = ofr;
      }
   }
   for(int r = 1; r < tlen; r++){
      for(int d = 0; d < tlen; d++){ 
         cov[d] += df[r]*dg[r+d+ofc];
         cov[d] += df[r+d+ofc]*dg[r];
         if(mp[r] < cov[d]*invn[r]*invn[r+d+ofc]){
            mp[r] = cov[d]*invn[r]*invn[r+d+ofc];
            mpi[r] = static_cast<long long>(r+d+ofr+ofc);
	 }
	 if(mp[r+d+ofc] < cov[d]*invn[r]*invn[r+d+ofc]){
            mp[r+d+ofc] = cov[d]*invn[r]*invn[r+d+ofc];
            mpi[r+d+ofc] = static_cast<long long>(r+ofr);
         }
      }
   }
}



