#ifndef __MOMENTS__
#define __MOMENTS__

void center_query(const double* __restrict ts, const double* __restrict mu, double* __restrict q, int sublen);

void batchcov(const double* __restrict ts, const double* __restrict mu, const double* __restrict query, double* __restrict cov, int count, int sublen);

void sw_mean(double* __restrict a, double* __restrict mu, int len, int winlen);

void sw_inv_meancentered_norm(double* __restrict a, double* __restrict mu, double* __restrict invn, int len, int winlen);

#endif
