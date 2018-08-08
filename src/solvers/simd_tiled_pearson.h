#include "../utils/reg.h"
#include "../utils/reduce.h"
#include "../utils/vector_print_funcs.h"
#define simlen 4
#define unroll 8
using namespace avx256_t;

// can't recall what initially provided the struct definition, but I'll probably remove it soon anyway


void  pauto_pearson_AVX_kern (double*       __restrict__ cov, 
                              double*       __restrict__ mp, 
                              long long*    __restrict__ mpi, 
                              const double* __restrict__ df, 
                              const double* __restrict__ dg, 
                              const double* __restrict__ invn, 
                              int ofr, 
                              int ofc,
                              int iters){

   cov =  (double*)      __builtin_assume_aligned(cov,prefalign);
   mp  =  (double*)      __builtin_assume_aligned(mp,prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi,prefalign);
   df  =  (const double*)__builtin_assume_aligned(df,prefalign);
   dg  =  (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);
   
  
   for(int r = 0; r < iters; r++){
      double t1[simlen] = {-1.0, -1.0, -1.0, -1.0}; // This is an alternative to the use of horizontal operations, since we only need a scalar result.
      long long t2[simlen] = {0, 0, 0, 0};         
      struct rpair v;   
      for(int d = 0; d < 128; d += simlen * unroll){
         block<__m256d> cov_r;
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = aload(cov, d + simlen * sd);
         }
         __m256d q;
         if(r != 0){
            q = brdcst(dg,r);
            for(int sd = 0; sd < unroll; sd++){
               cov_r(sd) = mul_add(q, uload(df, r + d + ofc + simlen * sd), cov_r(sd));
            }
            q = brdcst(df, r);
            for(int sd = 0; sd < unroll; sd++){
               cov_r(sd) = mul_add(q, uload(dg, r + d + ofc + simlen * sd), cov_r(sd));
               prvec("cov_r", cov_r(sd));
               astore(cov_r(sd), cov, d + simlen * sd);
            }
         }
         q = brdcst(invn,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) *= q;
         }
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) *= uload(invn, r + d + ofc + simlen * sd);
         }
         block<__m256i> mask;
         for(int sd = 0; sd < unroll; sd++){
            mask(sd) = cov_r(sd) < uload(mp, r + d + ofc + simlen * sd);
         }
         __m256i s = brdcst(r + ofr);
         for(int sd = 0; sd < unroll; sd++){
            if(testnz(mask(sd))){
               maskstore(cov_r(sd), mask(sd), mp + r + d + simlen * sd + ofc);
               maskstore(s, mask(sd), mpi + r + d + ofc + simlen * sd);
            }
         }
         struct rpair v = max_reduce_8x1(cov_r(0), cov_r(1), cov_r(2), cov_r(3), cov_r(4), cov_r(5), cov_r(6), cov_r(7));
         __m256i msk = v.val > aload(&t1[0],0);
         if(testnz(msk)){ 
            maskstore(v.val, msk, &t1[0]);
            maskstore(v.index + brdcst(ofr + ofc), msk, &t2[0]);
         }
      }
      __m256d t3 = aload(&t1[0],0) > brdcst(mp,r);
      if(testnz(t3)){  
         double r0, r1;
         double i0, i1;
         i0 = t1[0] > t1[2] ? t1[0] : t1[2];
         i1 = t1[1] > t1[3] ? t1[1] : t1[3];
         r0 = t1[0] > t1[2] ? t1[0] : t1[2];
         r1 = t1[1] > t1[3] ? t1[1] : t1[3];
         mpi[r] = r0 > r1 ? i0 : i1;
         mp[r] = r0 > r1 ? r0 : r1;
      } 
   }
}



