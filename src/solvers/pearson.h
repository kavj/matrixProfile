#include "../utils/descriptors.h"

typedef primbuf<double> dsbuf;
typedef primbuf<long long> lsbuf;  
typedef multibuf<double> mdsbuf; 
typedef multibuf<long long> mdsibuf;
typedef double dtype;

int nautocorr_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, long long minlag, long long sublen);

void pearson2zned(double* mp, long long len, long long sublen);


namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};
