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
      double t0[simlen] = {-1.0, -1.0, -1.0, -1.0}; 
      long long t1[simlen] = {0, 0, 0, 0};         
      for(int d = 0; d < 128; d += 8 * simlen){
         __m256d cv0, cv1, cv2, cv3, cv4, cv5, cv6, cv7;
         if(r != 0){
            __m256d q0 = _mm256_broadcast_sd(df + r);
            
            __m256d r0 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 0 * simlen + ofc), _mm256_load_pd(cov + d + 0 * simlen));
            __m256d r1 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 1 * simlen + ofc), _mm256_load_pd(cov + d + 1 * simlen));
            __m256d r2 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 2 * simlen + ofc), _mm256_load_pd(cov + d + 2 * simlen));
            __m256d r3 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 3 * simlen + ofc), _mm256_load_pd(cov + d + 3 * simlen));
            __m256d r4 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 4 * simlen + ofc), _mm256_load_pd(cov + d + 4 * simlen));
            __m256d r5 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 5 * simlen + ofc), _mm256_load_pd(cov + d + 5 * simlen));
            __m256d r6 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 6 * simlen + ofc), _mm256_load_pd(cov + d + 6 * simlen));
            __m256d r7 = _mm256_fmadd_pd(q0, _mm256_loadu_pd(dg + r + d + 7 * simlen + ofc), _mm256_load_pd(cov + d + 7 * simlen));
            
            __m256d q1 = _mm256_broadcast_sd(dg + r); 
            
            cv0 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 0 * simlen + ofc), r0);
            cv1 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 1 * simlen + ofc), r1);
            cv2 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 2 * simlen + ofc), r2);
            cv3 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 3 * simlen + ofc), r3);
            cv4 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 4 * simlen + ofc), r4);
            cv5 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 5 * simlen + ofc), r5);
            cv6 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 6 * simlen + ofc), r6);
            cv7 = _mm256_fmadd_pd(q1, _mm256_loadu_pd(df + r + d + 7 * simlen + ofc), r7);
            
            _mm256_store_pd(cov + d + 0 * simlen, cv0);
            _mm256_store_pd(cov + d + 1 * simlen, cv1);
            _mm256_store_pd(cov + d + 2 * simlen, cv2);
            _mm256_store_pd(cov + d + 3 * simlen, cv3);
            _mm256_store_pd(cov + d + 4 * simlen, cv4);
            _mm256_store_pd(cov + d + 5 * simlen, cv5);
            _mm256_store_pd(cov + d + 6 * simlen, cv6);
            _mm256_store_pd(cov + d + 7 * simlen, cv7);          
         }  
         else{
            cv0 = _mm256_load_pd(cov + d + 0 * simlen);
            cv1 = _mm256_load_pd(cov + d + 1 * simlen);
            cv2 = _mm256_load_pd(cov + d + 2 * simlen);
            cv3 = _mm256_load_pd(cov + d + 3 * simlen);
            cv4 = _mm256_load_pd(cov + d + 4 * simlen);
            cv5 = _mm256_load_pd(cov + d + 5 * simlen);
            cv6 = _mm256_load_pd(cov + d + 6 * simlen);
            cv7 = _mm256_load_pd(cov + d + 7 * simlen);
         }
         __m256d q0 = _mm256_broadcast_sd(invn + r);

         __m256d cr0 = _mm256_mul_pd(_mm256_mul_pd(cv0, q0), _mm256_loadu_pd(invn + r + d + 0 * simlen + ofc));
         __m256d cr1 = _mm256_mul_pd(_mm256_mul_pd(cv1, q0), _mm256_loadu_pd(invn + r + d + 1 * simlen + ofc));
         __m256d cr2 = _mm256_mul_pd(_mm256_mul_pd(cv2, q0), _mm256_loadu_pd(invn + r + d + 2 * simlen + ofc));
         __m256d cr3 = _mm256_mul_pd(_mm256_mul_pd(cv3, q0), _mm256_loadu_pd(invn + r + d + 3 * simlen + ofc));
         __m256d cr4 = _mm256_mul_pd(_mm256_mul_pd(cv4, q0), _mm256_loadu_pd(invn + r + d + 4 * simlen + ofc));
         __m256d cr5 = _mm256_mul_pd(_mm256_mul_pd(cv5, q0), _mm256_loadu_pd(invn + r + d + 5 * simlen + ofc));
         __m256d cr6 = _mm256_mul_pd(_mm256_mul_pd(cv6, q0), _mm256_loadu_pd(invn + r + d + 6 * simlen + ofc));
         __m256d cr7 = _mm256_mul_pd(_mm256_mul_pd(cv7, q0), _mm256_loadu_pd(invn + r + d + 7 * simlen + ofc));
        
         __m256i msk0 = _mm256_castpd_si256(cr0 > _mm256_loadu_pd(mp + r + d + 0 * simlen + ofc));
         __m256i msk1 = _mm256_castpd_si256(cr1 > _mm256_loadu_pd(mp + r + d + 1 * simlen + ofc));
         __m256i msk2 = _mm256_castpd_si256(cr2 > _mm256_loadu_pd(mp + r + d + 2 * simlen + ofc));
         __m256i msk3 = _mm256_castpd_si256(cr3 > _mm256_loadu_pd(mp + r + d + 3 * simlen + ofc));
         __m256i msk4 = _mm256_castpd_si256(cr4 > _mm256_loadu_pd(mp + r + d + 4 * simlen + ofc));
         __m256i msk5 = _mm256_castpd_si256(cr5 > _mm256_loadu_pd(mp + r + d + 5 * simlen + ofc));
         __m256i msk6 = _mm256_castpd_si256(cr6 > _mm256_loadu_pd(mp + r + d + 6 * simlen + ofc));
         __m256i msk7 = _mm256_castpd_si256(cr7 > _mm256_loadu_pd(mp + r + d + 7 * simlen + ofc));

         __m256i ind = _mm256_set1_epi64x(r);
 
         if(_mm256_testnzc_si256(msk0, msk0)){
            _mm256_maskstore_pd(mp + r + d + 0 * simlen + ofc, msk0, cr0);
            _mm256_maskstore_epi64(mpi + r + d + 0 * simlen + ofc, msk0, ind);
         }
         if(_mm256_testnzc_si256(msk1, msk1)){
            _mm256_maskstore_pd(mp + r + d + 1 * simlen + ofc, msk1, cr1);
            _mm256_maskstore_epi64(mpi + r + d + 1 * simlen + ofc, msk1, ind);
         }
         if(_mm256_testnzc_si256(msk2, msk2)){
            _mm256_maskstore_pd(mp + r + d + 2 * simlen + ofc, msk2, cr2);
            _mm256_maskstore_epi64(mpi + r + d + 2 * simlen + ofc, msk2, ind);
         }
         if(_mm256_testnzc_si256(msk3, msk3)){
            _mm256_maskstore_pd(mp + r + d + 3 * simlen + ofc, msk3, cr3);
            _mm256_maskstore_epi64(mpi + r + d + 3 * simlen + ofc, msk3, ind);
         }
         if(_mm256_testnzc_si256(msk4, msk4)){
            _mm256_maskstore_pd(mp + r + d + 4 * simlen + ofc, msk4, cr4);
            _mm256_maskstore_epi64(mpi + r + d + 4 * simlen + ofc, msk4, ind);
         }
         if(_mm256_testnzc_si256(msk5, msk5)){
            _mm256_maskstore_pd(mp + r + d + 5 * simlen + ofc, msk5, cr5);
            _mm256_maskstore_epi64(mpi + r + d + 5 * simlen + ofc, msk5, ind);
         }
         if(_mm256_testnzc_si256(msk6, msk6)){
            _mm256_maskstore_pd(mp + r + d + 6 * simlen + ofc, msk6, cr6);
            _mm256_maskstore_epi64(mpi + r + d + 6 * simlen + ofc, msk6, ind);
         }
         if(_mm256_testnzc_si256(msk7, msk7)){
            _mm256_maskstore_pd(mp + r + d + 7 * simlen + ofc, msk7, cr7);
            _mm256_maskstore_epi64(mpi + r + d + 7 * simlen + ofc, msk7, ind);
         }
         
         __m256i rd00msk = _mm256_castpd_si256(cr1 > cr0);
         __m256i rd01msk = _mm256_castpd_si256(cr3 > cr2);
         __m256i rd02msk = _mm256_castpd_si256(cr5 > cr4);
         __m256i rd03msk = _mm256_castpd_si256(cr7 > cr6);

         __m256d rd00 = _mm256_max_pd(cr1, cr0);
         __m256d rd01 = _mm256_max_pd(cr3, cr3);
         __m256d rd02 = _mm256_max_pd(cr5, cr4);
         __m256d rd03 = _mm256_max_pd(cr7, cr6);

         __m256i rd00ind = _mm256_blendv_epi8(_mm256_set1_epi64x(r + d + 0 * simlen + ofc), _mm256_set1_epi64x(r + d + 1 * simlen + ofc), rd00msk);
         __m256i rd01ind = _mm256_blendv_epi8(_mm256_set1_epi64x(r + d + 2 * simlen + ofc), _mm256_set1_epi64x(r + d + 3 * simlen + ofc), rd01msk);
         __m256i rd02ind = _mm256_blendv_epi8(_mm256_set1_epi64x(r + d + 4 * simlen + ofc), _mm256_set1_epi64x(r + d + 5 * simlen + ofc), rd02msk);
         __m256i rd03ind = _mm256_blendv_epi8(_mm256_set1_epi64x(r + d + 6 * simlen + ofc), _mm256_set1_epi64x(r + d + 7 * simlen + ofc), rd03msk);       

         __m256i rd10msk = _mm256_castpd_si256(rd00 > rd01);
         __m256i rd11msk = _mm256_castpd_si256(rd02 > rd03);

         __m256d rd10 = _mm256_max_pd(rd00, rd01);
         __m256d rd11 = _mm256_max_pd(rd02, rd03);

         __m256i rd10ind = _mm256_blendv_epi8(rd01ind, rd00ind, rd10msk);
         __m256i rd11ind = _mm256_blendv_epi8(rd02ind, rd03ind, rd11msk);

         __m256i rd20msk = _mm256_castpd_si256(rd10 > rd11);
         __m256d rd20 = _mm256_max_pd(rd10, rd11);
         __m256i rd20ind = _mm256_blendv_epi8(rd11ind, rd10ind, rd20msk);
         
         __m256i tmp = _mm256_castpd_si256(rd20 > _mm256_load_pd(mp + r));
         if(_mm256_testnzc_si256(tmp, tmp)){
            __m256i rdcmp = _mm256_castpd_si256(rd20 > _mm256_loadu_pd(t0));
            _mm256_maskstore_pd(t0, rdcmp, rd20);
            _mm256_maskstore_epi64(t1, rdcmp, rd20ind); 
         }
      }
      // __m256i t3 = _mm256_castpd_si256(_mm256_load_pd(&t0[0]));
      __m256d t3 = _mm256_loadu_pd(&t0[0]);
      if(_mm256_testnzc_pd(t3,t3)){  
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






