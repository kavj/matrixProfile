#include "../mp/avx2arith_2.hpp"





inline void incsum(double &a, double &b, double &err){
   double x = a + b;
   double z = x - a;
   double y = ((a - (x - z) + (b - z));
   
}

