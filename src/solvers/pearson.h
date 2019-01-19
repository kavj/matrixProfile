#include "../utils/descriptors.h"

int partialauto(bufd& ts, bufd& mp, bufi& mpi, long long minlag, long long sublen);

// cross
//int partialauto(bufd& ts, bufd& mp, bufi& mpi, long long minlag, long long sublen);


namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};
