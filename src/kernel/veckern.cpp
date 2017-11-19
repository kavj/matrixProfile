#include<cstdio>
#include "../mp/avx2arith.hpp"


using namespace vmth;
#define ksz 64
typedef  __m256d vtf;
typedef  __m256i vti;




#define innerWid 40
#define tileSz   512

void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    //double a[tileSz];
    unsigned long t3 = 0;
    for(int i = 0; i < n; i +=tileSz){
        for(int j = i; j < n;j+=4){
            for(int k = 0; k < tileSz; k+=innerWid){
                vtf a1 = bcast(dF,i+k); 
                vtf a2 = bcast(dF,j+k);
                for(int l = 0; l < innerWid; l+=4){
                   storea(mult_add(a2,loada(dX,j+l),mult_add(a1,loada(dX,i+l),loada(Cxy,i+l))),Cxy,i+l);
                }
                t3+=40;
          }
       }
    }
    printf("%lu %d \n",t3,n);
}






