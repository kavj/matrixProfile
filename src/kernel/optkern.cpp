#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;



static inline vtf ker_help(double* dx,double* df,double* s, int diag, int offset, int m){
    vtf z = setzero();
    vtf c1 = setzero();
    for(int suboff = 0; suboff < m; suboff += 4){
         vtf b1 = bcast(dx,suboff + offset);
         vtf d1 = max(z,flabs(loada(df,suboff+offset)-loada(dx,diag+suboff))*loada(s,suboff+offset) - b1);
         c1 = mult_add(d1,d1,c1);
    }
    return c1;
}


static inline void ker(double* dx, double* df, double* s, double* output, int diag, int offset, int m){
   for(int k = 0; k < 32; k+=4){
      vtf c1 = ker_help(dx,df,s,diag,offset,m);
     /* for(int suboff = 0; suboff < m; suboff += 4){
         vtf b1 = bcast(dx,suboff + offset);
         vtf b2 = bcast(dx,suboff + diag);
         vtf d1 = max(z,flabs(loada(df,suboff+offset)-b1)*loada(s,suboff+offset) - b2);
         c1 = mult_add(d1,d1,c1);
      }*/
      storea(c1,output,offset+k);
   }
}


/*
static inline void ker(double* dx, double* df, double* s, double* output, int diag, int offset, int m){
   vtf z = setzero();
   vtf t1 = _mm256_set1_pd(0x7fffffff);
   for(int k = 0; k < 32; k+=4){
      vtf c1 = setzero();
      for(int suboff = 0; suboff < m; suboff += 4){
         vtf b1 = bcast(dx,suboff + offset);
         vtf d1 = max(z,flabs(loada(df,suboff+offset)-loada(dx,diag+suboff))*loada(s,suboff+offset) - b1);
         c1 = mult_add(d1,d1,c1);
      }
      storea(c1,output,offset+k);
   }
}
*/




/*
static inline void actest_10_help(double* cxy, double* dx, double* df, double* s, double* output, int diag, int offset, int m){
   for(int i = 0; i < 32; i+= 4){
      for(int j = 0; j < m; j++){
         
      }
   }
}
*/

void  accumTest4_7_9(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n,int m){
   unsigned long t3 = 0;
   for(int diag = 0; diag < n-4096; diag += 32){
      for(int offset = 0; offset < n-diag-32*128; offset+=4){
         ker(dx,df,s,output,diag,offset,128);
         t3+= 32;      
     }
   }
   printf("%lu %d \n",t3,n);
}



