#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;



void  accumTest4_7_9(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n,int m){
   unsigned long t3 = 0;
   for(int diag = 0; diag < n-4096; diag += 32){
      vtf z = setzero();
  for(int offset = 0; offset < n-diag-32*128; offset+=4){
          vtf c1 = setzero();
          vtf c2 = setzero();
          vtf c3 = setzero();
          vtf c4 = setzero();
          vtf c5 = setzero();
          vtf c6 = setzero();
          vtf c7 = setzero();
          vtf c8 = setzero();
         for(int suboff = 0; suboff < 128; suboff+=4){
            vtf t1 = _mm256_set1_pd(0x7fffffff);
            vtf b1 = bcast(dx,suboff + offset);
            vtf d1 = max(z,flabs(loada(df,suboff+offset)-loada(dx,diag+suboff))*loada(s,suboff+offset) - b1);
            vtf d2 = max(z,flabs(loada(df,suboff+offset+4)-loada(dx,diag+suboff+4))*loada(s,suboff+offset+4) - b1);
            vtf d3 = max(z,flabs(loada(df,suboff+offset+8)-loada(dx,diag+suboff+8))*loada(s,suboff+offset+8) - b1);
            vtf d4 = max(z,flabs(loada(df,suboff+offset+12)-loada(dx,diag+suboff+12))*loada(s,suboff+offset+12) - b1);
            vtf d5 = max(z,flabs(loada(df,suboff+offset+16)-loada(dx,diag+suboff+16))*loada(s,suboff+offset+16) - b1);
            vtf d6 = max(z,flabs(loada(df,suboff+offset+20)-loada(dx,diag+suboff+20))*loada(s,suboff+offset+20) - b1);
            vtf d7 = max(z,flabs(loada(df,suboff+offset+24)-loada(dx,diag+suboff+24))*loada(s,suboff+offset+24) - b1);
            vtf d8 = max(z,flabs(loada(df,suboff+offset+28)-loada(dx,diag+suboff+28))*loada(s,suboff+offset+28) - b1);

            c1 = mult_add(d1,d1,c1);
            c2 = mult_add(d2,d2,c2);
            c3 = mult_add(d3,d3,c3);
            c4 = mult_add(d4,d4,c4);
            c5 = mult_add(d5,d5,c5);
            c6 = mult_add(d6,d6,c6);
            c7 = mult_add(d7,d7,c7);
	    c8 = mult_add(d8,d8,c8);
         }
         t3+= 32;
      
         storea(c1,output,offset);
         storea(c2,output,offset+4);
         storea(c3,output,offset+8);
         storea(c4,output,offset+12);
         storea(c5,output,offset+16);
         storea(c6,output,offset+20);
         storea(c7,output,offset+24);
         storea(c8,output,offset+28);
     }
   }
   printf("%lu %d \n",t3,n);
}


