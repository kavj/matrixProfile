#include "../utils/descriptors.h"

typedef strided_buffer<double> dbuf; 
typedef strided_buffer<long long> ibuf;

int nautocorr_reduc(dbuf& ts, dbuf& mp, ibuf& mpi, long long minlag, long long sublen);

void pearson2zned(double* mp, long long len, long long sublen);


namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};
