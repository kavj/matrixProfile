
#ifndef PAUTO
#define PAUTO

void  pauto_pearson_kern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ mp, long long* __restrict__ mpi, int offsetr, int offsetc);


void pauto_pearson_edgekern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, double* __restrict__ mp, long long* __restrict__ mpi, int offsetr, int offsetc, int bound);


void pauto_pearson_refkern(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp);

#endif
