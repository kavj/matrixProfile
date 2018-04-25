#include <cmath>
#include "../arch/avx256.h"
#include "../utils/reg.h"


using namespace avx256_t;
typedef double dtype;
typedef __m256d vtype;

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


//Todo: cleanup edge handling
template<typename dtype,typename vtype>
void batchcov(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int count, int sublen){
   const int simlen = sizeof(vtype)/sizeof(dtype);
   const int stride = simlen*unroll;
   const int block_count = count/stride;
   const int lastaligned = count - count%stride;
   for(int i = 0; i < count; i+=stride){
      block<vtype> c;
      block<vtype> m;
      for(int j = 0; j < unroll; j++){
         m(j) = mu(i+j*simlen);
      }
      for(int j = i; j < i+sublen; j++){
         for(int k = 0; k < unroll; k++){
            c(k+unroll) = uload(ts,j) - m(k);
         }
         for(int k = 0; k < unroll; k++){
            c(k) = mul_add(uload(query,j),c(k+unroll),c(k));
         }
      }
      for(int j = 0; j < unroll; j++){
         ustore(c(k),cov,i+j*simlen);
      }
   }
   // fringe should be done here. It would be possible to use a while loop countdown, like std::min(stride,count-i), and I need to add either a backshift or scalar loop that won't break
}


