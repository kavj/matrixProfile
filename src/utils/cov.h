#include "../arch/avx256.h"
#define prefalign 64
// Todo: Split reference and simd types


template<typename dtype, typename itype>
void center_query(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ q, itype sublen){
   q = (double*)__builtin_assume_aligned(q,prefalign);
   for(itype i = 0; i < sublen; i++){
      q[i] = ts[i] - mu[0];
   }
}

template<typename dtype, typename itype>
void batchcov(const dtype* __restrict__ ts, const dtype* __restrict__ mu, const dtype* __restrict__ query, dtype* __restrict__ cov, itype count, itype sublen){
   for(itype i = 0; i < count; i++){
      cov[i] = 0; 
      for(itype j = 0; j < sublen; j++){
         cov[i] += (ts[i + j] - mu[i]) * query[j];
      }
   }
}

