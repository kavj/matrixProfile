#include<cstdio>
#include "../mp/avx2arith.hpp"
#define innerWid 8 
#define tileSz  512 

using namespace vmth;
#define ksz 64
typedef  __m256d vtf;
typedef  __m256i vti;
/*
void accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    
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
}*/

void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    double a[innerWid];
    double b[innerWid];
    double e[innerWid];
    double f[innerWid];
    unsigned long t3 = 0;
    for(int i = 0; i < n; i +=tileSz){
        for(int j = i; j < n;j+=4){
            vtf a1 = bcast(dF,i); 
            vtf a2 = bcast(dF,j);
            vtf s2 = bcast(s,i);
            __m256i i1 = loada(outputI,i);
            for(int k = 0; k < tileSz; k+=innerWid){
                for(int l = 0; l < innerWid; l+=4){
                    storea(mult_add(a1,loada(dX,i+k+l),loada(Cxy,i+k+l)),a,l);
             #ifdef split 
                }
                for(int l = 0; l < innerWid; l+=4){
                    storea(mult_add(a2,loada(dX,j+k+l),loada(a,l)),b,l);
                }
                for(int l = 0; l < innerWid; l+=4){
                    storea(mult(s2,loada(s,j+k+l)),e,l);
                }
                for(int l = 0; l < innerWid; l+=4){
                    storea(mult(loada(e,l),loada(b,l)),f,l);
                }
                for(int l = 0; l < innerWid; l+=4){
                    storea(max(loada(f,k+l),loada(output,i+k+l)),output,i+k+l);
             #else
                    storea(mult_add(a2,loada(dX,j+k+l),loada(a,l)),b,l);
                    storea(mult(s2,loada(s,j+k+l)),e,l);
                    storea(mult(loada(e,l),loada(b,l)),f,l);
                    storea(max(loada(f,k+l),loada(output,i+k+l)),output,i+k+l);
             #endif
                    t3 += 4;
                }
           }
            
        }
    }
    printf("%lu %d \n",t3,n);
}


