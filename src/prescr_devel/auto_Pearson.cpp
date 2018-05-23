#include "../utils/max_reduce.h"
#include "../arch/avx256.h"
#include "../utils/reg.h"
using namespace avx256_t;

#define tsz 65536
#define unroll 8
#define simlen 4
#define step 32

// Todo: separate different kernel implementations by namespace. This makes them selectable at compile time without messy ifdef nonsense






void   pauto_pearson_kern (double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ mp, long long* __restrict__ mpi, int offsetr, int offsetc){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df  = (const double*)__builtin_assume_aligned(df,32);
   dx  = (const double*)__builtin_assume_aligned(dx,32);
   invn   = (const double*)__builtin_assume_aligned(invn,32);
   mp  = (double*)__builtin_assume_aligned(mp,32);
   mpi = (long long*)__builtin_assume_aligned(mpi,32);
   for(int i = 0; i < tsz; i++){
      for(int j = 0; j < tsz; j+=step){
         block<__m256d> cov_r;
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


void pauto_pearson_edgekern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ mp, long long* __restrict__ mpi, int offsetr, int offsetc, int bound){
   for(int i = 0; i < bound; i++){
      double c = cov[i];
      for(int j = i; j < bound; j++){
         c += df[j]*dx[j+offsetc];
         c += df[j+offsetc]*dx[j];
         double corr = c*invn[i]*invn[i+j];
         if(corr > mp[i]){
            mp[j] = corr;
            mpi[j] = j+offsetr+offsetc;
         }
         if(corr > mp[i+j+offsetc]){
            mp[j+offsetc] = corr;
            mpi[j+offsetc] = j+offsetr;
         }
      }
   }
}


// If we don't have an optimized intrinsic version for the target
// we can usually get reasonable code if we guarantee alignment conditions, tile shape, and lack of aliasing

void pauto_pearson_refkern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);

   for(int i = 0; i < tsz; i++){
      for(int j = 0; j < tsz; j++){
         cov[j] += dx[i]*df[i+j+offsetc];
      }
      double a[tsz];
      for(int j = 0; j < tsz; j++){
         cov[j] += df[i]*dx[i+j+offsetc];
         a[j] = cov[j]*s[i]*s[i+j+offsetc];
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i] < a[j]){
            mp[i] = a[j];
            mpi[i] = i+j+offsetc;
         }
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i+j+offsetc] < a[j]){
            mp[i+j+offsetc] = a[j];
            mpi[i+j+offsetc] = i+j+offsetr;
         }
      }
   }
}


