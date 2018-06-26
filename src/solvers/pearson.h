#include "descriptors.h"
#define prefalign 64

typedef stridedbuf<double> dsbuf;
typedef stridedbuf<long long> lsbuf;

void pearson_pauto_reduc(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         int minlag, 
                         int sublen);

