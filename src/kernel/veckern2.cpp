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



void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   double a[tileSz*simdWid];
   double b[tileSz*simdWid];
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += tileSz){
      for(int offset = 0; offset < diag; offset += tileSz){
         for(int suboff = offset; suboff < offset + tileSz; suboff++){
            vtf s1 = bcast(s,suboff);
            vtf dFi = bcast(dF,suboff);
            vtf dXi = bcast(dX,suboff);
            vtf mpi = bcast(output,suboff);
            vtf mpiI = bcast(output,suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = loada(Cxy,subdiag);
               c1 = mult_add(loadu(dF,suboff+subdiag),dXi,mult_add(dFi,loadu(dX,suboff+subdiag),c1));
               vtf c2 = mult(c1,s1);
               c2 = mult(c2,loadu(s,suboff+subdiag));
               mpi = max(c2,mpi);
               storea(c1,Cxy,subdiag-diag);
               t3 += 4;
            }
            storeu(mpi,output,suboff);
         }
      }   
   }
   printf("%lu %d \n",t3,n);
}




void  accumTest4_2(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   double a[tileSz*simdWid];
   double b[tileSz*simdWid];
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
            vtf c1 = loada(Cxy,subdiag);
            for(int suboff = offset; suboff < offset + innerWid; suboff++){
               c1 = mult_add(bcast(dF,suboff),loadu(dX,subdiag),c1);
               c1 = mult_add(bcast(dX,suboff),loadu(dF,subdiag),c1);
               //storea(c1,a,simdWid*(suboff-offset));
               vtf s1 = mult(bcast(s,subdiag),loadu(s,suboff));
               vtf s2 = mult(c1,s1);
               storea(s2,b,simdWid*(suboff-offset));
               t3 += 4;
            }
            storea(c1,Cxy,subdiag);
            c1 = loada(output,subdiag);
            
            for(int i = 0; i < innerWid; i+=simdWid){
               c1 = max(c1,loada(b,i));
            }
            storea(c1,output,subdiag);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}


void  accumTest4_3(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
            vtf c1 = loada(Cxy,subdiag);
            vtf c2 = loada(output,subdiag);
            for(int suboff = 0; suboff < innerWid; suboff++){
               c1 = mult_add(bcast(dF,suboff+offset),loadu(dX,subdiag),c1);
               c1 = mult_add(bcast(dX,suboff+offset),loadu(dF,subdiag),c1);
               vtf s1 = mult(bcast(s,subdiag),loadu(s,suboff+offset));
               vtf s2 = mult(c1,s1);
               c2 = max(s2,c2);
               storea(c2,output,4*suboff+subdiag);
               t3 += 4;
            }
            storea(c1,Cxy,subdiag);

         }
      }
   }
   printf("%lu %d \n",t3,n);
}



void  accumTest4_4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff++){
            vtf c2 = bcast(output,suboff);
            vtf s1 = bcast(s,suboff);
            vtf dx1 = bcast(dX,suboff);
            vtf df1 = bcast(dF,suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = mult_add(df1,loadu(dX,suboff),loada(Cxy,subdiag));
               c1 = mult_add(dx1,loadu(dF,suboff),c1);
               vtf s2 = mult(s1,loadu(s,suboff+offset));
               s2 = mult(c1,s2);
               c2 = max(s2,c2);
               t3 += 4;
               storea(s2,a,innerWid*suboff+(subdiag-diag));
               storea(c1,Cxy,subdiag);
            }
            storea(c2,output,4*suboff);
         }
         for(int i = 0; i < innerWid; i+= simdWid){
            vtf c1 = loada(output,diag);
            c1 = max(c1,loada(a,i));
            c1 = max(c1,loadu(a,i+innerWid+1));
            c1 = max(c1,loadu(a,i+2*(innerWid+1)));
            c1 = max(c1,loadu(a,i+3*(innerWid+1)));
            storea(c1,output,diag);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}




void  accumTest4_5(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff++){
            vtf c2 = bcast(output,suboff);
            vtf s1 = bcast(s,suboff);
            vtf dx1 = bcast(dX,suboff);
            vtf df1 = bcast(dF,suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = mult_add(df1,loadu(dX,suboff),loada(Cxy,subdiag));
               c1 = mult_add(dx1,loadu(dF,suboff),c1);
               vtf s2 = mult(s1,loadu(s,suboff+offset));
               s2 = mult(c1,s2);
               c2 = max(s2,c2);
               t3 += 4;
               storea(c1,Cxy,subdiag);
            }
            storea(c2,output,4*suboff);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}



