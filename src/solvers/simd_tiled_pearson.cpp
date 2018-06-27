
static inline void  pauto_pearson_AVX_kern (double*       __restrict__ cov, 
                                            double*       __restrict__ mp, 
                                            long long*    __restrict__ mpi, 
                                            const double* __restrict__ df, 
                                            const double* __restrict__ dg, 
                                            const double* __restrict__ invn, 
                                            int ofr, 
                                            int ofc){

   cov =  (double*)      __builtin_assume_aligned(cov,prefalign);
   mp  =  (double*)      __builtin_assume_aligned(mp,prefalign);
   mpi =  (long long*)   __builtin_assume_aligned(mpi,prefalign);
   df  =  (const double*)__builtin_assume_aligned(df,prefalign);
   dg  =  (const double*)__builtin_assume_aligned(dg,prefalign);
   invn = (const double*)__builtin_assume_aligned(invn,prefalign);

   for(int r = 0; r < klen; r++){
      for(int d = 0; d < klen; d+=step){
         block<__m256d> cov_r;
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = aload(cov, simlen *(d + sd));
         }
         __m256d q = brdcst(dg,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = mul_add(q,uload(df,r + d + simlen * sd),cov_r(sd));
         }
         q = brdcst(df,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) = mul_add(q,uload(dg, r + d + simlen * sd),cov_r(sd));
            astore(cov_r(sd), cov, simlen * (d + sd));
         }
         q = brdcst(invn,r);
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) *= q;
         }
         for(int sd = 0; sd < unroll; sd++){
            cov_r(sd) *= uload(invn, r + d + simlen * sd);
         }
         block<__m256i> mask;
         for(int sd = 0; sd < unroll; sd++){
            mask(sd) = cov_r(sd) > uload(mp, r + d + simlen * sd);
         }
         __m256i s = brdcst(r);
         for(int sd = 0; sd < unroll; sd++){
            if(testnz(mask(sd))){
               maskstore(cov_r(sd), mask(sd), mp + r + d + simlen * sd);
               maskstore(s,mask(sd), mpi + r + simlen * sd);
            }
         }
         struct rpair v = max_reduce_8x1(cov_r(0), cov_r(1), cov_r(2), cov_r(3), cov_r(4), cov_r(5), cov_r(6), cov_r(7));
         //Todo: technically this is incorrect, as we're writing back 4 operands rather than 1. It would be possible to do a final comparison against a temporary buffer and make one sweep at the end
         // alternatively use the if statement to determine if we need to do a horizontal max (since that is quite expensive)
         __m256i msk = v.val > brdcst(mp,r + d);
         if(testnz(msk)){
            maskstore(v.val, msk, mp + r + d);
            maskstore(v.index + brdcst(r + d + ofc), msk, mpi + r + d);
         }
      }
   }
}




