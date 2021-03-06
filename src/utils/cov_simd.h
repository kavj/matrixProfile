#include "../arch/avx256.h"
#include "reg.h"
#define prefalign 64
// Todo: Split reference and simd types

template<typename dtype>
void center_query(const dtype* __restrict__ ts, const dtype* __restrict__  mu,  dtype* __restrict__ q, int sublen){
   q = (dtype*) __builtin_assume_aligned(q,prefalign);
   const int aligned = sublen - sublen%prefalign; 
   const int unroll = 8;
   const int simlen = 4;
   const int loopwid = prefalign;
   for(int j = 0; j < aligned; j+=loopwid){
      __m256d m = brdcst(mu,0);
      block<__m256d> op;
      for(int k = 0; k < unroll; k++){
         op(k) = uload(ts, j + k * simlen) - m;
      }
      for(int k = 0; k < unroll; k++){
         astore(op(k), q, j + k * simlen);
      }
   }
   dtype m = mu[0];
   for(int j = aligned; j < sublen; j++){
      q[j] = ts[j] - m;
   }
}

template<typename dtype>
void batchcov(const dtype* __restrict__ ts, const dtype* __restrict__ mu, const dtype* __restrict__ query, dtype* __restrict__ cov, int count, int sublen){
   cov =   (dtype*)       __builtin_assume_aligned(cov,prefalign);
   query = (const dtype*) __builtin_assume_aligned(query,prefalign);
   mu =    (const dtype*) __builtin_assume_aligned(mu,prefalign);
   const int simlen = 4;
   const int loopwid = prefalign;
   const int unroll = 8;
   #pragma omp parallel for
   for(int i = 0; i < count; i += loopwid){
      block<dtype> c;
      block<dtype> m;
      for(int j = 0; j < unroll; j++){
         m(j) = uload(mu, i + j * simlen);
      }
      for(int j = 0; j < sublen; j++){
         for(int k = 0; k < unroll; k++){
            c(k+unroll) = uload(ts, i + j + k * simlen) - m(j); 
         }
         dtype q = brdcst(query,j);
         for(int k = 0; k < unroll; k++){
            c(k) = mul_add(q, c(k + unroll), c(k));
         }
      }
      for(int j = 0; j < unroll; j++){
         astore(c(j), cov, i + j * simlen);
      }
   }   // fringe should be done here. it would be possible to use a while loop countdown, like std::min(stride,count-i), and i need to add either a backshift or scalar loop that won't break
   if(count % loopwid){
      int fringe = count - count % loopwid;
      batchcov_ref(ts + fringe, cov + fringe, query, mu + fringe, count - fringe, sublen);
   }
}


