#include "../arch/avx256.h"
#include "../utils/reg.h"

using namespace avx256_t;
typedef double dtype;
typedef __m256d vtype;





//template<typename dtype>
void batchcov_ref(const dtype* __restrict__ ts, 
                        dtype* __restrict__ cov, 
                  const dtype* __restrict__ query, 
                  const dtype* __restrict__ mu, 
                  int count, 
                  int sublen);

void batchcov(const dtype* __restrict__ ts, 
                    dtype* __restrict__ cov, 
              const dtype* __restrict__ query, 
              const dtype* __restrict__ mu, 
                int count, 
                int sublen);

void center_query(const dtype* __restrict__ ts, 
                  const dtype* __restrict__ mu, 
                        dtype* __restrict__ q,  
                    int sublen);


void center_query_ref(const dtype* __restrict__ ts, 
                      const dtype* __restrict__ mu, 
                            dtype* __restrict__ q, 
                        int sublen);

