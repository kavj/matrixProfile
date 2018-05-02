
#ifndef FFS
#define FFS

#include<algorithm>
#include<cstdint>
#include<cmath>
#include "../utils/reg.h"
//#include "descriptors.h"
/*
struct query_stat{
   inline query_stat();
   inline query_stat(double qcov, int qind) : qcov(qcov), qind(qind){}
   double qcov;
   int qind;
};
*/

struct rpair rescaled_max_reduct(double* __restrict__ cov,  const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset);
//struct query_stat mx_reduct(const double* __restrict__ cov, const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset);
//struct query_stat mx_reduct(const double* __restrict__ cov, const double* __restrict__ invn, double* __restrict__ xcorr, int* __restrict__ cindex, double qinvn, double qcorr, int qbaseind, int cindoffset);

#endif
