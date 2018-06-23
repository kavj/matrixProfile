

static inline void  pauto_pearson_AVX_kern(
   double* __restrict__ cov, 
   double* __restrict__ mp, 
   long long* __restrict__ mpi, 
   const double* __restrict__ df, 
   const double* __restrict__ dx, 
   const double* __restrict__ invn, 
   int offsetr,
   int offsetc);


void pauto_pearson_xinner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc);


void pauto_pearson_basic_inner_alt(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc);



void pauto_pearson_basic_inner(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc);

void pauto_pearson_xedge(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   long long*    __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int tlen,
   int offsetr, 
   int offsetc, 
   int bound);


void pauto_pearson_edge(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   long long*    __restrict__ mpi,
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


