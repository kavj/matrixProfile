#include "../arch/avx256.h"
#include "../utils/reg.h"
using namespace avx256_t;

struct rpair { __m256d val; __m256i index;};

 
static inline struct rpair __attribute__((always_inline)) reduce(const __m256d& r0, const __m256d& r1, const __m256d& r2, const __m256d& r3,  const __m256d& r4, const __m256d& r5, const __m256d& r6, const __m256d& r7){

   __m256i mask1, mask2, mask3, mask4, mask5, mask6;
   __m256d  op1, op2, op3, op4;

   mask1 = r0 > r1;
   mask2 = r2 > r3;
   mask3 = r4 > r5;
   mask4 = r6 > r7;

   static const __m256i four = brdcst(4);
   static const __m256i eight = brdcst(8);
  
   op1 = vmax(r0,r1);
   op2 = vmax(r2,r3);
   op3 = vmax(r4,r5);
   op4 = vmax(r6,r7);
  
   __m256i g = set(0,1,2,3); 
   __m256i h = g + four;
   mask1 = blend(g,h,mask1);
   g += eight;
   h += eight;
   mask2 = blend(g,h,mask2);
   g += eight;
   h += eight;
   mask3 = blend(g,h,mask3);
   g += eight;
   h += eight;
   mask4 = blend(g,h,mask4);
   
   mask5 = op1 > op2;
   mask6 = op3 > op4;
   
   op1 = vmax(op1,op2);
   op3 = vmax(op3,op4);

   mask1 = blend(mask1,mask2,mask5);
   mask3 = blend(mask3,mask4,mask6);
   
   struct rpair x;
   x.val = vmax(op1,op3);
   x.index = blend(mask1,mask3,op1 > op3);
   return x;
}


#define tsz 64
#define unroll 8
#define simlen 4


void  symm_pearson_reduct_kern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);

   for(int i = 0; i < tsz; i++){
      for(int j = 0; j < tsz; j+=32){
         block<__m256d> cov_r;
         for(int k = 0; k < unroll; k++){
            cov_r(k) = aload(cov,simlen*(j+k));
         }
         __m256d q = brdcst(dx,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(df,i+j+4*k),cov_r(k));
         }
         q = brdcst(df,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(dx,i+j+simlen*k),cov_r(k));
            astore(cov_r(k),cov,simlen*(j+k));
         }
         q = brdcst(s,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= q;
         }
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= uload(s,i+j+simlen*k);
         }
         block<__m256i> mask;
         for(int k = 0; k < unroll; k++){
            mask(k) = cov_r(k) > uload(mp,i+j+simlen*k);
         }
         for(int k = 0; k < unroll; k++){
            __m256i r = brdcst(i);
            if(testnz(mask(k))){
               maskstore(cov_r(k),mask(k),mp+i+j+simlen*k);
               maskstore(r,mask(k),mpi+i+simlen*k);
            }
         }
         block<__m256d> c = cov_r;
         struct rpair r = reduce(c(0),c(1),c(2),c(3),c(4),c(5),c(6),c(7));
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



// If we don't have an optimized intrinsic version for the target
// we can usually get reasonable code if we guarantee alignment conditions, tile shape, and lack of aliasing

#define tsz 64 
void reference_pearson_reduc(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
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
      for(int j = 0; j < tsz; j++){
         cov[j] += df[i]*dx[i+j+offsetc];
      }
      double a[64];
      for(int j = 0; j < tsz; j++){
         a[j] = cov[j]*s[i]*s[i+j+offsetc];
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i] < a[j]){
            mp[i] = a[j];
            mpi[i] = i+offsetc+j;
         }
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i+j+offsetc] < a[j]){
            mp[i+j+offsetc] = a[j];
            mpi[i+j+offsetc] = i+j+offsetc;
         }
      }
   }
}


