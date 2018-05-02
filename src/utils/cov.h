#include "../arch/avx256.h"
#include "../utils/reg.h"

using namespace avx256_t;
typedef double dtype;
typedef __m256d vtype;





//template<typename dtype>
void batchcov_reference(const dtype* __restrict__ ts, dtype* cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen);

void batchcov_simd(const dtype* __restrict__ ts, dtype* __restrict__ cov, const dtype* __restrict__ query, const dtype* __restrict__ mu, int offset, int count, int sublen);

