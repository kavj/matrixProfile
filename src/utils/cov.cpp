#include <cmath>
#include "../arch/scalar.h"
#include "../utils/reg.h"

using namespace scalar;//avx256_t;
typedef double dtype;
typedef double vtype;//__m256d vtype;

// todo: move reference to separate file, make scalar namespace

static inline void batchcov_kern(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int unroll, int sublen) __attribute__ ((always_inline));


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

template<typename dtype, typename vtype,int unroll>
static inline void batchcov_kern(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int sublen){
   const int simlen = sizeof(vtype)/sizeof(dtype);
   block<vtype> c;
   for(int i = offset; i < offset+sublen; i++){
      vtype q = brdcst(mu,i-offset);
      for(int j = 0; j < unroll; j++){
         c(j+unroll) = uload(ts,i+simlen*j) - q;
      }
      q = brdcst(query,i);
      for(int j = 0; j < unroll; j++){
         c(j) = mul_add(q,c(j+unroll),c(j));
      }
   }
   vtype q = brdcst(1.0);
   for(int i = 0; i < unroll; i++){
      astore(q/c(i),cov,offset+simlen*i);
   }
}

//template<typename dtype,typename vtype>
void batchcov(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen){
   const int unroll = 8;
   const int simlen = sizeof(vtype)/sizeof(dtype);
   const int stride = simlen*unroll;
   const int simd_block_count = count/stride;
   const int last = count - simd_block_count*stride;   
   if(simd_block_count > 0){
      for(int i = offset; i < offset+count; i+=stride){  
         batchcov_kern<dtype,vtype,unroll>(ts,cov,query,mu,i,sublen);
      }
      if(offset+count > last){
        //last-(count+offset)
         batchcov_kern<dtype,vtype,1>(ts,cov,query,mu,last,sublen); 
      }
   }
   else{
      batchcov_reference(ts,cov,query,mu,offset,count,sublen);
   }
}
