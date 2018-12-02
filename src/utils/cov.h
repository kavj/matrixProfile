#define prefalign 64 
// Todo: Split reference and simd types
#define wid 32
#define unroll 8

void center_query(const double* restrict ts, const double* restrict mu, double* restrict q, int sublen);

void batchcov(const double* restrict ts, const double* restrict mu, const double* restrict query, double* restrict cov, int count, int sublen);
