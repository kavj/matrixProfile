#ifndef __MOMENTS__
#define __MOMENTS__

void centerQuery(const double* __restrict ts, double* __restrict qBuffer, double mu, int sublen);

void crossCov(const double* __restrict ts, const double* __restrict mu, const double* __restrict query, double* __restrict cov, int count, int sublen);

void windowedMean(double* __restrict a, double* __restrict mu, int len, int winlen);

void windowedInverseCenteredNorm(double* __restrict a, double* __restrict mu, double* __restrict invn, int len, int winlen);

#endif
