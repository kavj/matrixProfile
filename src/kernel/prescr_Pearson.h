

// mmay not be needed

void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int count, int sublen, int qstart, int qstride, int qbufstride);


// This is probably the best I can do without avx
// I'm not completely happy with this, because it relies on updating a non-const reference

//struct query_stat rescaled_max_reduct(double* __restrict__ cov, double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, int qbaseind, int cindoffset, int count);
struct query_stat rescaled_max_reduct(double* __restrict__ cov,  double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, double qinvn, int qbaseind, int count, int cindoffset);

// This can probably be merged with the normal version and an unroll size of 1

// need to add structs in here
// This might be reworked later to remove temporary buffers or to wrap the messy formulas in small inline functions, not sure really but it's too hard to read with inline formulas
void maxpearson_ext_auto(const double* __restrict__ qcov, const double* __restrict__ invn, const double* __restrict__ df, const double* __restrict__ dx, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len);




