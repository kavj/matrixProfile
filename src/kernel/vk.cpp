#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 16 
#define tileSz 128 
#define alignedInner (innerWid/4-1)

using namespace vmth;
#define ksz 64
typedef  __m256d vtf;
typedef  __m256i vti;



/* self join with unaligned loads removed */
void symmkern(double* Cxy, double* dF, double* dX, double* S, double* mp, long* mpI, int diag, int offs){
   double a[innerWid];
   for(int j = diag; j < diag + innerWid; j += simdWid){
      vtf csum = loada(Cxy(diag));
      for(int iter = 0; iter < innerWid; iter+=simdWid){
         vtf opF = bcast(dF,offs);
         vtf opG = loada(dF,j+offs);
         vtf opG2 = preshuffle(opG,loada(dF,j+offs+simdWid));
         vtf opX = bcast(dX,offs);
         vtf opY = loada(dX,j+offs);
         csum = mult_add(opG,opY,mult_add(opF,opX,csum));
         storea(a,iter); 
         vtf opG3 = shift1(opG,opG2);
          
   }

}
