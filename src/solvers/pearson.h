#include "../utils/descriptors.h"
#define prefalign 64

typedef primbuf<double> dsbuf;
typedef primbuf<long long> lsbuf;
typedef multibuf<double> mdsbuf;

void pearson_pauto_reduc(dsbuf& ts, 
                         dsbuf& mp, 
                         lsbuf& mpi, 
                         int minlag, 
                         int sublen);

