#include "../utils/descriptors.h"

typedef primbuf<double> dsbuf;
typedef primbuf<int> lsbuf;  
typedef multibuf<double> mdsbuf; 
typedef multibuf<int> mdsibuf;
typedef double dtype;

int pearson_pauto_reduc(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         int minlag, 
                         int sublen);

int pearson_pauto_reduc_ref(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         int minlag, 
                         int sublen);


int pearson_pauto_tileref(dsbuf& ts, 
                          dsbuf& mp, 
                          lsbuf& mpi, 
                          int minlag, 
                          int sublen);

namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};
