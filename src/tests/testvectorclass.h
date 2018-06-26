#include<immintrin.h>


class vtf{
   protected:
      __m256d vd;   

   public:
   vtf(double* a, int i){
      vd = _mm256_load_pd(a+i);
   }

   vtf operator __m256d () {
      return vd;
   }

   vtf & operator = (__m256d &a){
      vd = a;
   }
  
   vtf & operator = (vtf &a){
      vd = a.vd;
   }

   static inline vtf & operator + (vtf &a, vtf &b){
      return _mm256_add_pd(a,b);
   }
};
