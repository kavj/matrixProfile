#include <cmath>
#include "../arch/avx256.h"
#include "../utils/reg.h"

using namespace avx256_t;
typedef double dtype;
typedef __m256d vtype;



static inline void batchcov_kern_scalar(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int sublen) __attribute__ ((always_inline));
static inline void cov_kern_simd(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int unroll, int sublen) __attribute__ ((always_inline));


//template<typename dtype>
void batchcov_reference(const dtype* __restrict__ ts, dtype* cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   for(int i = 0; i < count; i++){
      dtype c = 0;
      dtype m = mu[i];
      for(int j = 0; j < sublen; j++){
         c += query[j]*(ts[i+j]-m);
      }
      cov[i] = c;
   }
}

//template<typename dtype, typename vtype>
static inline void cov_kern_simd(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int unroll, int sublen){
   block<vtype> c;
   block<vtype> m;
   block<vtype> aux;
   for(int i = 0; i < unroll; i++){
      m(i) = aload(mu,offset+simlen*i);
   }
   for(int i = offset; i < offset+sublen; i++){
      for(int j = 0; j < unroll; j++){
         aux(j) = uload(ts,i+simlen*j) - m(j);
      }
      vtype q = brdcst(query,i-offset);
      for(int j = 0; j < unroll; j++){
         c(j) = mul_add(q,aux(j),c(j));
      }
   }
   for(int i = 0; i < unroll; i++){
      astore(c(i),cov,offset+simlen*i);
   }
}

//template<typename dtype>
static inline void batchcov_kern_scalar(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int sublen){
   const int block_lim = offset + (count/unroll)*unroll;   
   block<dtype> c;
   block<dtype> m;
   block<dtype> aux;
   for(int i = 0; i < unroll; i++){
      m(i) = mu[offset+i];
   }
   for(int i = offset; i < offset+sublen; i++){
      for(int j = 0; j < unroll; j++){
         aux(j) = ts[i+j]-m(j);
      }
      dtype q = query[i-offset];
      for(int j = 0; j < unroll; j++){
         c(j) = fma(q,aux(j),c(j));
      }
   } 
   for(int i = 0; i < unroll; i++){
      cov[i+offset] = c(i);
   }
}

// not sure I want to keep this, probably only want the AVX type, otherwise just use fft assuming stability doesn't suck
//template<typename dtype>
void batchcov_scalar(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   const int block_lim = offset + (count/unroll)*unroll;   
   const int block_count = count/unroll;
   if(block_count > 0){
      for(int i = offset; i < block_lim + offset; i+= unroll){
         batchcov_kern_scalar(ts,cov,query,mu,i,sublen);
      }
      batchcov_kern_scalar(ts,cov,query,mu,offset+count-unroll,sublen);
   }
   else{
      batchcov_reference(ts,cov,query,mu,offset,count,sublen);
   }
}


//template<typename dtype,typename vtype>
void batchcov_simd(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   const int simlen = sizeof(vtype)/sizeof(dtype);
   const int stride = simlen*unroll;
   const int simd_block_count = count/stride;
   const int last = count - simd_block_count*stride;   
   if(simd_block_count > 0){
      for(int i = offset; i < block_lim + offset; i+=stride){  // <-- I got rid of block_lim, so this needs to be fixed
         cov_kern_simd(ts,cov,query,mu,i,unroll,sublen);
      }
      if(block_lim+offset < last){
         cov_kern_simd(ts,cov,query,mu,block_lim,last-(block_lim+offset),sublen); 
      }
   }
   else{
      batch_cov_reference(ts,cov,query,mu,offset,count,sublen);
   }
}


