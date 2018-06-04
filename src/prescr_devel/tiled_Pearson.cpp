#include<algorithm>
//#include "../utils/max_reduce.h"
#include "../utils/reg.h"
#include "../utils/cov.h"
#include "../utils/xprec_math.h"
#include "descriptors.h"
#include "../utils/primitive_print_funcs.h"
#ifdef prefalign
#undef prefalign
#endif
#define prefalign 64 

//#define prefalign 32 
#define klen 64 

// rename and reorganize later. This should be in a different namespace or compilation unit
/*
static inline void  pauto_pearson_AVX_kern (double* __restrict__ cov, double* __restrict__ mp, long long* __restrict__ mpi, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, int offsetr, int offsetc){
   cov = (double*)__builtin_assume_aligned(cov,prefalign);
   df  = (const double*)__builtin_assume_aligned(df,prefalign);
   dx  = (const double*)__builtin_assume_aligned(dx,prefalign);
   invn   = (const double*)__builtin_assume_aligned(invn,prefalign);
   mp  = (double*)__builtin_assume_aligned(mp,prefalign);
   mpi = (long long*)__builtin_assume_aligned(mpi,prefalign);
   for(int i = 0; i < klen; i++){
      for(int j = 0; j < klen; j+=step){
         block<__m256d> cov_r;
         if( (j != 0)){
         for(int k = 0; k < unroll; k++){
            cov_r(k) = aload(cov,simlen*(j+k));
         }
         __m256d q = brdcst(dx,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(df,i+j+simlen*k),cov_r(k));
         }
         q = brdcst(df,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(dx,i+j+simlen*k),cov_r(k));
            astore(cov_r(k),cov,simlen*(j+k));
         }
         q = brdcst(invn,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= q;
         }
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= uload(invn,i+j+simlen*k);
         }
         block<__m256i> mask;
         for(int k = 0; k < unroll; k++){
            mask(k) = cov_r(k) > uload(mp,i+j+simlen*k);
         }
         __m256i s = brdcst(i);
         for(int k = 0; k < unroll; k++){
            if(testnz(mask(k))){
               maskstore(cov_r(k),mask(k),mp+i+j+simlen*k);
               maskstore(s,mask(k),mpi+i+simlen*k);
            }
         }
         struct rpair r = max_reduce_8x1(cov_r(0),cov_r(1),cov_r(2),cov_r(3),cov_r(4),cov_r(5),cov_r(6),cov_r(7));
         //Todo: technically this is incorrect, as we're writing back 4 operands rather than 1. It would be possible to do a final comparison against a temporary buffer and make one sweep at the end
         // alternatively use the if statement to determine if we need to do a horizontal max (since that is quite expensive)
         __m256i msk = r.val > brdcst(mp,i+j);
         if(testnz(msk)){
            maskstore(r.val,msk,mp+i+j);
            maskstore(r.index+brdcst(i+j+offsetc),msk,mpi+i+j);
         }
      }
   }
}
}*/


/*
void pauto_pearson_xinner(
   double*      __restrict__ cov,
   double*      __restrict__ mpr,
      int*      __restrict__ mpri,
   double*      __restrict__ mpc,
      int*      __restrict__ mpci,
const double*   __restrict__ df,
const double*   __restrict__ dg,
const double*   __restrict__ dx,
const double*   __restrict__ dy,
const double*   __restrict__ invnf,
const double*   __restrict__ invnx,
int tlen,
int offsetr,
int offsetc,
int bound)
{
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mpr =   (double*)__builtin_assume_aligned(mpr,prefalign);
   mpri =  (int*)__builtin_assume_aligned(mpri,prefalign);
   mpc =  (double*) __builtin_assume_aligned(mpc,prefalign);
   mpci = (int*) __builtin_assume_aligned(mpci,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   dx =   (const double*)__builtin_assume_aligned(dx,prefalign);
   dy =   (const double*)__builtin_assume_aligned(dy,prefalign);
   invnf = (const double*)__builtin_assume_aligned(invnf,prefalign);
   invnx = (const double*)__builtin_assume_aligned(invnx,prefalign);
    
   for(int d = 0; d < tlen; d++){
      if(cov[d]*invnf[0]*invnx[d+offsetc] > mpr[0]){
         mpr[0] = cov[d]*invnf[0]*invnx[d+offsetc];
         mpri[0] = d+offsetr+offsetc;
      }
      if(cov[d]*invnf[0]*invnx[d+offsetc] > mpc[d+offsetc]){
         mpc[d] = cov[d]*invnf[0]*invnx[d+offsetc];
         mpci[d] = offsetr;
      }
   }
   for(int r = 0; r < tlen; r++){
      
   }

}
*/



auto accum = [&] (auto mp, auto mpi,auto cov,auto df,auto dg,auto invn, auto offsetr, auto offsetc, auto itersr, auto itersc){
   for(auto r = 0; r < itersr; r++){
      for(auto i = 0; i < itersc; i++){
         cov[i] += df[i]*dg[i+offsetc];
         cov[i] += df[i+offsetc]*dg[i];
      }
      for(auto i = 0; i < itersc; i++){
         if(cov[i]*invn[i]*invn[i+offsetc] > mp[i]){
            mp[i] = cov[i]*invn[i]*invn[i+offsetc];
            mpi[i] = i+offsetr+offsetc;
         }
      }
      for(auto i = 0; i < itersc; i++){
         if(cov[i]*invn[i]*invn[i+offsetc] > mp[i]){
            mp[i+offsetc] = cov[i]*invn[i]*invn[i+offsetc];
            mpi[i+offsetc] = i+offsetr;
         }
      }
   }
};

auto reduce = [&] (auto mp, auto mpi,auto cov, auto invn, auto offsetr, auto offsetc, auto iters){
   for(auto i = 0; i < iters; i++){
      if(cov[i]*invn[i]*invn[i+offsetc] > mp[i]){
         mp[i] = cov[i]*invn[i]*invn[i+offsetc];
         mpi[i] = i+offsetr+offsetc;
      }
      if(cov[i]*invn[i]*invn[i+offsetc] > mp[i]){
         mp[i+offsetc] = cov[i]*invn[i]*invn[i+offsetc];
         mpi[i+offsetc] = i+offsetr;
      }
   }
};

void pauto_pearson_colproj(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   int*          __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int tlen,
   int offsetr, 
   int offsetc, 
   int bound)
{
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mp =   (double*)__builtin_assume_aligned(mp,prefalign);
   mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int c = 0; c < c + tlen; c++){
      for(int d = 0; d < c; d++){
         cov[d] += df[c+offsetc]*dg[c-d];
         cov[d] += df[c-d]*dg[c+offsetc];
      }
      for(int d = 0; d < c; d++){
         if(cov[d]*invn[c-d]*invn[c+offsetc] > mp[c-d]){
            mp[c-d] = cov[d]*invn[c-d]*invn[c+offsetc];
            mpi[c-d] = c+offsetc+offsetr;
         }
      }
      for(int d = 0; d < c; d++){
         if(cov[d]*invn[c-d]*invn[c+offsetc] > mp[c+offsetc]){
            mp[c+offsetc] = cov[d]*invn[c-d]*invn[c+offsetc];
            mpi[c+offsetc] = c-d+offsetr;
         }
      } 
   } 
}



void pauto_pearson_edge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   int*          __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int tlen,
   int offsetr, 
   int offsetc, 
   int bound)
{
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mp =   (double*)__builtin_assume_aligned(mp,prefalign);
   mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);
 
   for(int d = 0; d < std::min(tlen,bound); d++){
      if(cov[d]*invn[0]*invn[d+offsetc] > mp[0]){
         mp[0] = cov[d]*invn[0]*invn[d+offsetc];  
         mpi[0] = d+offsetr+offsetc;
      }
      if(cov[d]*invn[0]*invn[d+offsetc] > mp[d+offsetc]){
         mp[d+offsetc] = cov[d]*invn[0]*invn[d+offsetc]; 
         mpi[d+offsetc] = offsetr;
      }
      for(int r = 1; r < std::min(tlen,bound-d); r++){
         cov[d] += df[r]*dg[r+d+offsetc];
         cov[d] += df[r+d+offsetc]*dg[r];
         if(cov[d]*invn[r]*invn[r+d+offsetc] > mp[r]){
            mp[r] = cov[d]*invn[r]*invn[r+d+offsetc];
            mpi[r] = r+d+offsetc+offsetr;
         }
         if(cov[d]*invn[r]*invn[r+d+offsetc] > mp[r+d+offsetc]){
            mp[r+d+offsetc] = cov[d]*invn[r]*invn[r+d+offsetc];
            mpi[r+d+offsetc] = r+offsetr;
         }
      }
   }
}

void pauto_pearson_basic_inner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc)
{
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mp =   (double*)__builtin_assume_aligned(mp,prefalign);
   mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int d = 0; d < tlen; d++){
      if(mp[0] < cov[d]*invn[0]*invn[d+offsetc]){
         mp[0] = cov[d]*invn[0]*invn[d+offsetc];
         mpi[0] = d+offsetr+offsetc;
      }
      if(mp[d] < cov[d]*invn[0]*invn[d]){
         mp[d] = cov[d]*invn[0]*invn[d];
         mpi[d] = offsetr;
      }
   }
   for(int r = 1; r < tlen; r++){
      for(int d = 0; d < tlen; d++){ 
         cov[d] += df[r]*dg[r+d+offsetc];
         cov[d] += df[r+d+offsetc]*dg[r];
         if(mp[r] < cov[d]*invn[r]*invn[r+d+offsetc]){
            mp[r] = cov[d]*invn[r]*invn[r+d+offsetc];
            mpi[r] = r+d+offsetr+offsetc;
	 }
	 if(mp[r+d+offsetc] < cov[d]*invn[r]*invn[r+d+offsetc]){
            mp[r+d+offsetc] = cov[d]*invn[r]*invn[r+d+offsetc];
            mpi[r+d+offsetc] = r+offsetr;
         }
      }
   }
}

void pauto_pearson_inner_alt(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc)
{
   cov =  (double*)__builtin_assume_aligned(cov,prefalign);
   mp =   (double*)__builtin_assume_aligned(mp,prefalign);
   mpi =  (int*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int d = 0; d < tlen; d += klen){
      reduce(mp,mpi,cov+d,invn,offsetr,d+offsetc,klen);
      accum(mp,mpi,cov+d,df,dg,invn,0,d+offsetc,klen-1,klen);
      for(int r = klen; r < tlen; r+=klen){
         accum(cov+d,mp+r,mpi+r,df+r,dg+r,invn+r,r,d+offsetc,klen,klen);
      }
   }
}


