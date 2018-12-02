
#include<math.h>

void pearson_to_normalized_euclidean(double* restrict a, long long len, long long sublen){
   double n = (double) (2 * sublen);
   for(long long i = 0; i < len; i++){
      a[i] = sqrt(n * (1.0 - a[i]));
   }
}
