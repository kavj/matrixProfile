#include "../utils/descriptors.h"

typedef primbuf<double> dsbuf;
typedef primbuf<int> lsbuf;      // Todo: We probably need to decide how to deal with the int vs long long issue. In the manually tuned version, it's faster to match datatype widths, but 
typedef multibuf<double> mdsbuf; //       compiler auto vectorization uses swizzles either way and tends to do better if it touches less memory
typedef double dtype;

void pearson_pauto_reduc(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         int minlag, 
                         int sublen);

