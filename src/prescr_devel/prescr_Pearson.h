
#ifndef FFS
#define FFS
#include "../utils/max_reduce.h"
/*struct rpair{
   __m256d val;
   __m256i ind;
}*/

struct query_stat{
  // inline query_stat(double qcov, int qind) : qcov(qcov), qind(qind){}
   double qcov;
   double qcorr;
   int qind;
};

struct rpair rescaled_max_reduct(double* __restrict__ cov,  const double* __restrict__ invn, double* __restrict__ xcorr, long long* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset);
struct query_stat rescaled_max_reduct_ref(double* __restrict__ cov,  const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset);

#endif
