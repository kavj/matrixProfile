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




static inline void  pauto_pearson_AVX_kern (double* __restrict__ cov, double* __restrict__ mp, long long* __restrict__ mpi, const double* __restrict__ df, const double* __restrict__ dg, const double* __restrict__ invn, int offsetr, int offsetc){
   cov = (double*)__builtin_assume_aligned(cov,prefalign);
   df  = (const double*)__builtin_assume_aligned(df,prefalign);
   dg  = (const double*)__builtin_assume_aligned(dg,prefalign);
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
            __m256d q = brdcst(dg,i);
            for(int k = 0; k < unroll; k++){
               cov_r(k) = mul_add(q,uload(df,i+j+simlen*k),cov_r(k));
            }
            q = brdcst(df,i);
            for(int k = 0; k < unroll; k++){
               cov_r(k) = mul_add(q,uload(dg,i+j+simlen*k),cov_r(k));
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
}


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



void pauto_pearson_xedge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
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
   mpi =  (long long*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);
   
   const int dmx = std::min(bound,tlen);
   // this should be initialize, 
 
   for(int d = 0; d < dmx; d += klen){
      const int rmx = std::min(bound-d,tlen);
      const int ralgn = rmx - rmx%klen - 2*klen;
      for(int r = klen; r < ralgn; r += klen){
         pauto_pearson_AVX_kern(cov, mp, mpi, df, dg, invn, offsetr+r, offsetc+d);
      }     
      for(int r = ralgn; r < rmx; r++){  
         for(int subd = d; subd < std::min(bound-r,d+klen); subd++){
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
}



void pauto_pearson_edge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
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
   mpi =  (long long*)__builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);
   printf("offsetr:%d offsetc:%d\n",offsetr,offsetc);
 
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


void pauto_pearson_xinner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc)
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
         pauto_pearson_AVX_kern(cov+d,mp+r,mpi+r,df+r,dg+r,invn+r,offsetr+r,offsetc+d);
      }
   }
}


void pauto_pearson_basic_inner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc)
{
   cov =  (double*)      __builtin_assume_aligned(cov,prefalign);
   mp =   (double*)      __builtin_assume_aligned(mp,prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi,prefalign);
   df =   (const double*)__builtin_assume_aligned(df,prefalign);
   dg =   (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int d = 0; d < tlen; d++){
      if(mp[0] < cov[d]*invn[0]*invn[d+offsetc]){
         mp[0] = cov[d]*invn[0]*invn[d+offsetc];
         mpi[0] = d+offsetr+offsetc;
      }
      if(mp[d] < cov[d]*invn[0]*invn[d+offsetc]){
         mp[d] = cov[d]*invn[0]*invn[d+offsetc];
         mpi[d] = offsetr;
      }
      if(cov[d]*invn[0]*invn[d+offsetc] > 1.0){
         printf("r: %d d: %d\n",offsetr,d+offsetr+offsetc);
         exit(0);
      }
      else if(cov[d]*invn[0]*invn[d+offsetc] < -1.0){
         printf("r: %d d: %d\n",offsetr,d+offsetr+offsetc);
         exit(0);
      }
   }
   for(int r = 1; r < tlen; r++){
      for(int d = 0; d < tlen; d++){ 
         cov[d] += df[r]*dg[r+d+offsetc];
         cov[d] += df[r+d+offsetc]*dg[r];
         if(cov[d]*invn[r]*invn[r+d+offsetc] > 1.0){
            printf("r: %d d: %d\n",r+offsetr,d+r+offsetr+offsetc);
            exit(0);
         }
         else if(cov[d]*invn[r]*invn[r+d+offsetc] < -1.0){
            printf("r: %d d: %d\n",r+offsetr,d+r+offsetr+offsetc);
            exit(0);
         }
         if(mp[r] < cov[d]*invn[r]*invn[r+d+offsetc]){
            mp[r] = cov[d]*invn[r]*invn[r+d+offsetc];
            mpi[r] = static_cast<long long>(r+d+offsetr+offsetc);
	 }
	 if(mp[r+d+offsetc] < cov[d]*invn[r]*invn[r+d+offsetc]){
            mp[r+d+offsetc] = cov[d]*invn[r]*invn[r+d+offsetc];
            mpi[r+d+offsetc] = static_cast<long long>(r+offsetr);
         }
      }
   }
}



