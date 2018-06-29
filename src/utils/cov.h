#include "../arch/avx256.h"
#define prefalign 64
// Todo: Split reference and simd types


template<typename dtype>
void center_query(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ q, int sublen){
   q = (double*)__builtin_assume_aligned(q,prefalign);
   for(int i = 0; i < sublen; i++){
      q[i] = ts[i] - mu[0];
   }
}

template<typename dtype>
void batchcov(const dtype* __restrict__ ts, const dtype* __restrict__ mu, const dtype* __restrict__ query, dtype* __restrict__ cov, int count, int sublen){
   static int k = 0;
   for(int i = 0; i < count; i++){
      cov[i] = 0; 
      for(int j = 0; j < sublen; j++){
         cov[i] += (ts[i + j] - mu[i]) * query[j];
      }
   }
}

