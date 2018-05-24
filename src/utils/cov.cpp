#include<cstdio>
#include <algorithm>
#include <cmath>
#include "../arch/avx256.h"
#include "reg.h"
#include "checkArray.h"

using namespace avx256_t;
typedef double dtype;
typedef __m256d vtype;

// todo: move reference to separate file, make scalar namespace

//static inline void batchcov_kern(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int unroll, int sublen) __attribute__ ((always_inline));

void center_query_ref(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ q, int sublen){
   for(int j = 0; j < sublen; j++){
      q[j] = ts[j] - mu[0];
   }
}

void center_query(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ q,  int sublen){
   q = (dtype*) __builtin_assume_aligned(q,32);
   const int aligned = sublen - sublen%32; 
   const int unroll = 8;
   const int simlen = 4;
   const int loopwid = 32;
   for(int j = 0; j < aligned; j+=loopwid){
      __m256d m = brdcst(mu,0);
      block<__m256d> op;
      for(int k = 0; k < unroll; k++){
         op(k) = uload(ts,j+k*simlen) - m;
      }
      for(int k = 0; k < unroll; k++){
         astore(op(k),q,j+k*simlen);
      }
   }
   dtype m = mu[0];
   for(int j = aligned; j < sublen; j++){
      q[j] = ts[j] - m;
   }
}



void batchcov_ref(const dtype* __restrict__ ts, dtype* cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int count, int sublen){
   for(int i = 0; i < count; i++){
      cov[i] = 0; 
      for(int j = 0; j < sublen; j++){
         cov[i] += (ts[i+j] - mu[i])*query[j];
      }
   }
}

//Todo: cleanup edge handling
//template<typename dtype,typename vtype>

void batchcov(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int count, int sublen){
   cov = (dtype*) __builtin_assume_aligned(cov,32);
   query = (const dtype*) __builtin_assume_aligned(query,32);
   mu = (const dtype*) __builtin_assume_aligned(mu,32);
   const int simlen = 4;
   const int loopwid = 32;
   const int unroll = 8;
   #pragma omp parallel for
   for(int i = 0; i < count; i+=loopwid){
      block<vtype> c;
      block<vtype> m;
      for(int j = 0; j < unroll; j++){
         m(j) = uload(mu,i+j*simlen);
      }
      for(int j = 0; j < sublen; j++){
         for(int k = 0; k < unroll; k++){
            c(k+unroll) = uload(ts,i+j+k*simlen) - m(j); 
         }
         vtype q = brdcst(query,j);
         for(int k = 0; k < unroll; k++){
            c(k) = mul_add(q,c(k+unroll),c(k));
         }
      }
      for(int j = 0; j < unroll; j++){
         astore(c(j),cov,i+j*simlen);
      }
   }   // fringe should be done here. it would be possible to use a while loop countdown, like std::min(stride,count-i), and i need to add either a backshift or scalar loop that won't break
   if(count%loopwid){
      int fringe = count - count%loopwid;
      batchcov_ref(ts+fringe,cov+fringe,query,mu+fringe,count-fringe,sublen);
   }
}


