#include "../arch/avx256.h"
#include "../utils/reg.h"

using namespace avx256_t;
typedef double dtype;
typedef __m256d vtype;


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


//template<typename dtype>
// not sure I want to keep this, probably only want the AVX type, otherwise just use fft assuming stability doesn't suck

/*void batchcov_scalar(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   const int unroll = 8;
   const int block_lim = offset + (count/unroll)*unroll;   
   for(int i = offset; i < block_lim; i+=unroll){
      block<dtype> c;
      block<dtype> m;
      block<dtype> aux;
      for(int j = 0; j < unroll; j++){
         m(j) = mu[i+j];
      }
      for(int j = 0; j < sublen; j++){
         dtype q = query[j];
         for(int k = 0; k < unroll; k++){
            aux(k) = ts[i+j]-m(j);
         }
         for(int k = 0; k < unroll; k++){
            c(k) += q*aux(k);
         }
      } 
      for(int j = 0; j < unroll; j++){
         cov[i+j] = c(j);
      }
   }
   batchcov_reference(ts,cov,query,mu,offset+block_lim,count-block_lim,sublen);
}*/


static inline void cov_kern_simd(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int unroll, int sublen){
   block<vtype> c;
   block<vtype> m;
   block<vtype> aux;
   for(int j = 0; j < unroll; j++){
      m(j) = aload(mu,simlen*j);
   }
   for(int j = 0; j < sublen; j++){
      for(int k = 0; k < unroll; k++){
         aux(k) = uload(ts,i+j+simlen*k) - m(k);
      }
      vtype q = brdcst(query,i);
      for(int k = 0; k < unroll; k++){
         c(k) = mul_add(q,aux(k),c(k));
      }
   }
   for(int j = 0; j < unroll; j++){
      astore(c(j),cov,i+simlen*j);
   }
}


//template<typename dtype, typename vtype>

// It's easiest to just use one unroll width here. 
// Just allow it to repeat some work in the name of latency hiding. If there are fewer subsequences than the unroll width, then simply call reference instead

void batchcov_simd(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   const int unroll = 8;
   const int simlen = sizeof(vtype)/sizeof(dtype);
   const int stride = simlen*unroll;
   const int simd_block_count = count/stride;
   const int last = count - simd_block_count*stride;   
   // check whether it's too small here
   if(simd_block_count > 0){
      for(int i = offset; i < block_lim + offset; i+=stride){
         cov_kern_simd(ts, cov, query, mu, i, unroll,sublen);
      }
      if(block_lim+offset < last){
         cov_kern_simd(ts, cov, query, mu, block_lim,6,sublen); 
      }
   }
   else{
      batch_cov_reference(ts,cov,query,mu,offset,count,sublen);
   }
}

