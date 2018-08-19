#include "../utils/descriptors.h"

typedef primbuf<double> dsbuf;
typedef primbuf<long long> lsbuf;  
typedef multibuf<double> mdsbuf; 
typedef multibuf<long long> mdsibuf;
typedef double dtype;

int pearson_pauto_reduc(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         long long minlag, 
                         long long sublen);

int pearson_pauto_reduc_ref(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         long long minlag, 
                         long long sublen);


int pearson_pauto_tileref(dsbuf& ts, 
                          dsbuf& mp, 
                          lsbuf& mpi, 
                          long long minlag, 
                          long long sublen);

namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};
