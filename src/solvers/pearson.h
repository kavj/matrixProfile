#include "../utils/descriptors.h"

namespace pearson{


void tonormalizedeuclidean(double* mp, int len, int sublen);

int partialauto(bufd& ts, bufd& mp, bufi& mpi, int minlag, int sublen);

int partialcross(bufd& a, bufd& b, bufd& mp, bufi& mpi, int sublen);

namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};

}
