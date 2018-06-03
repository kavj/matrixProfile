


void pauto_pearson_inner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc);

void pauto_pearson_edge(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc,
   const int bound);



void pauto_pearson_reftest(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   int minlag,
   int mlen);


