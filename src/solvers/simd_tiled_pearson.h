#include <immintrin.h>
constexpr long long simlen  = 4;
constexpr long long unroll = 8;

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
         fblock cov_r;
         if(r != 0){
            __m256d q = brdcst(dg, r);
            fblock cv;
            for(int sd = 0; sd < unroll; sd++){
               cv(sd) = mul_add(q, uload(df, r + d + ofc + simlen * sd), aload(cov, d + simlen * sd));
            }
            q = brdcst(df, r);
            for(int sd = 0; sd < unroll; sd++){
               cov_r(sd) = mul_add(q, uload(dg, r + d + ofc + simlen * sd), cv(sd));
            }
            for(int sd = 0; sd < unroll; sd++){
               astore(cov_r(sd), cov, d + simlen * sd);
            }
         }
         else{
            for(int sd = 0; sd < unroll; sd++){
               cov_r(sd) = aload(cov, d + simlen * sd);
            }
         }
         __m256d q = brdcst(invn, r);
         fblock corr;
         for(int sd = 0; sd < unroll; sd++){
            corr(sd) = cov_r(sd) * q * uload(invn, r + d + ofc + simlen * sd);
         }
         iblock mask;
         for(int sd = 0; sd < unroll; sd++){
            mask(sd) = _mm256_castpd_si256(corr(sd) > uload(mp, r + d + ofc + simlen * sd));
         }
         __m256i s = brdcst(r + ofr);
         for(int sd = 0; sd < unroll; sd++){
            if(testz(mask(sd)) == 0){
               maskstore(corr(sd), mask(sd), mp + r + d + simlen * sd + ofc);
               maskstore(s, mask(sd), mpi + r + d + ofc + simlen * sd);
            }
         }
         struct rpair v = max_reduce_8x1(corr(0), corr(1), corr(2), corr(3), corr(4), corr(5), corr(6), corr(7));
         __m256i msk = _mm256_castpd_si256(v.val > aload(&t1[0],0));
         if(testz(msk) == 0){ 
            maskstore(v.val, msk, &t1[0]);
            maskstore(v.index + brdcst(r + d + ofr + ofc), msk, &t2[0]);
         }
      }
      __m256d t3 = aload(&t1[0],0) > brdcst(mp,r);
      if(testz(t3) == 0){  
         double r0, r1;
         long long i0, i1;
         i0 = t1[0] >= t1[1] ? t2[0] : t2[1];
         i1 = t1[2] >= t1[3] ? t2[2] : t2[3];
         r0 = t1[0] >= t1[1] ? t1[0] : t1[1];
         r1 = t1[2] >= t1[3] ? t1[2] : t1[3];
        /// if(r0 > mp[r] || r1 > mp[r]){
         mpi[r] = r0 >= r1 ? i0 : i1;
         mp[r] = r0 >= r1 ? r0 : r1;
        // }
      } 
   }
}






