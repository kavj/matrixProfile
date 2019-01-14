#include "../arch/avx256.h"
#define prefalign 64 
// Todo: Split reference and simd types
#define wid 32
#define unroll 8

void center_query(const double* __restrict__ ts, const double* __restrict__ mu, double* __restrict__ q, int sublen);

void batchcov(const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ query, double* __restrict__ cov, int count, int sublen);
