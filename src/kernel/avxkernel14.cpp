#include<cstdio>
#include "../mp/avx2arith_2.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;


/*
static inline void fma_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8, double *d, int i){
   c1 = fma(x1,bcast(d,i),c1);
   c2 = fma(x2,bcast(d,i+4),c2);
   c3 = fma(x3,bcast(d,i+8),c3);
   c4 = fma(x4,bcast(d,i+12),c4);
   c5 = fma(x5,bcast(d,i+16),c5);
   c6 = fma(x6,bcast(d,i+20),c6);
   c7 = fma(x7,bcast(d,i+24),c7);
   c8 = fma(x8,bcast(d,i+28),c8);
}*/


static inline void fma_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8, double *d, int i){
   c1 = fma(x1,aload(d,i),c1);
   c2 = fma(x2,aload(d,i+4),c2);
   c3 = fma(x3,aload(d,i+8),c3);
   c4 = fma(x4,aload(d,i+12),c4);
   c5 = fma(x5,aload(d,i+16),c5);
   c6 = fma(x6,aload(d,i+20),c6);
   c7 = fma(x7,aload(d,i+24),c7);
   c8 = fma(x8,aload(d,i+28),c8);
}



static inline void prod_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   c1 *= x1;
   c2 *= x2;
   c3 *= x3;
   c4 *= x4;
   c5 *= x5;
   c6 *= x6;
   c7 *= x7;
   c8 *= x8;
}

static inline void prod_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   c1 *= aload(d,i);
   c2 *= aload(d,i+4);
   c3 *= aload(d,i+8);
   c4 *= aload(d,i+12);
   c5 *= aload(d,i+16);
   c6 *= aload(d,i+20);
   c7 *= aload(d,i+24);
   c8 *= aload(d,i+28);
}

static inline void load_8x1(vtf &c1, vtf &c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   c1 = aload(d,i);
   c2 = aload(d,i+4);
   c3 = aload(d,i+8);
   c4 = aload(d,i+12);
   c5 = aload(d,i+16);
   c6 = aload(d,i+20);
   c7 = aload(d,i+24);
   c8 = aload(d,i+28);
}

static inline void shuffle_fwd_8(vtf &c1, vtf &c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   c8 = c7;
   c7 = c6;
   c6 = c5;
   c5 = c4;
   c4 = c3;
   c3 = c2;
   c2 = c1;
   c1 = aload(d,i);
}

static inline void store_8(vtf &c1, vtf &c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   astore(c1,d,i);
   astore(c2,d,i+4);
   astore(c3,d,i+8);
   astore(c4,d,i+12);
   astore(c5,d,i+16);
   astore(c6,d,i+20);
   astore(c7,d,i+24);
   astore(c8,d,i+28);
}


static inline void fma_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &x1, vtf &x2, vtf &x3, vtf &x4, double *d, int i){
   c1 = fma(x1,aload(d,i),c1);
   c2 = fma(x2,aload(d,i+4),c2);
   c3 = fma(x3,aload(d,i+8),c3);
   c4 = fma(x4,aload(d,i+12),c4);
}

static inline void prod_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   c1 *= x1;
   c2 *= x2;
   c3 *= x3;
   c4 *= x4;
}

static inline void prod_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, double *d, int i){
   c1 *= aload(d,i);
   c2 *= aload(d,i+4);
   c3 *= aload(d,i+8);
   c4 *= aload(d,i+12);
}

static inline void load_4x4(vtf &c1, vtf &c2, vtf &c3, vtf &c4, double *d, int i){
   c1 = aload(d,i);
   c2 = aload(d,i+4);
   c3 = aload(d,i+8);
   c4 = aload(d,i+12);
}

static inline void shuffle_fwd_4(vtf &c1, vtf &c2, vtf &c3, vtf &c4, double*d, int i){
   c4 = c3;
   c3 = c2;
   c2 = c1;
   c1 = aload(d,i);
}

static inline void store_4(vtf &c1, vtf &c2, vtf &c3, vtf &c4, double *d, int i){
   astore(c1,d,i);
   astore(c2,d,i+4);
   astore(c3,d,i+8);
   astore(c4,d,i+12);
}

void  accumTest4_7_9(double* cx,double* dx,double* df, double*s, double* a, double* output, long* outputi,int n, int m){
   unsigned long t = 0;
   for(int diag = 0; diag < n-32*m; diag += 32){
      vtf c1, c2, c3, c4;
      load_4x4(c1,c2,c3,c4,cx,diag);
      //vtf a[128];
     for(int offset = 0; offset < n-32*m; offset += 16){
        vtf x1, x2, x3, x4, f1, f2, f3, f4, s1, s2, s3, s4;
        load_4x4(x1,x2,x3,x4,dx,diag+offset);
        load_4x4(f1,f2,f3,f4,df,diag+offset);
        load_4x4(s1,s2,s3,s4,s,diag+offset);
        for(int k = 0; k < 64; k+= 16){
            int p = offset+diag+k;
            fma_4x4(c1,c2,c3,c4,x1,x2,x3,x4,df,offset+k);
            fma_4x4(c1,c2,c3,c4,f1,f2,f3,f4,dx,offset+k);
            prod_4x4(c1,c2,c3,c4,s1,s2,s3,s4);
            prod_4x4(c1,c2,c3,c4,s,offset+k);
            store_4(c1,c2,c3,c4,a,offset+k);
            shuffle_fwd_4(c1,c2,c3,c4,cx,diag+k);
            t += 16;
         }
         for(int k = 0; k < 128; k += 32){
            
         }
      }
   }     
   printf("iterations: %lu \n",t);
}

 
void  accumTest4_7_10(double* cx,double* dx, double* df, double*s, double* a, double* output, long* outputi,int n, int m){
   unsigned long t= 0;
   for(int diag = 0; diag < n-32*m; diag += 32){
      vtf c1, c2, c3, c4, c5, c6, c7, c8;
      load_8x1(c1,c2,c3,c4,c5,c6,c7,c8,cx,diag);
      //vtf a[128];
     for(int offset = 0; offset < n-32*m; offset += 64){
        vtf x1, x2, x3, x4, x5, x6, x7, x8,
            f1, f2, f3, f4, f5, f6, f7, f8,
            s1, s2, s3, s4, s5, s6, s7, s8;
        load_8x1(x1,x2,x3,x4,x5,x6,x7,x8,dx,diag+offset);
        load_8x1(f1,f2,f3,f4,f5,f6,f7,f8,df,diag+offset);
        load_8x1(s1,s2,s3,s4,s5,s6,s7,s8,s,diag+offset);
        for(int k = 0; k < 64; k+= 32){
            int p = offset+diag+k;
            fma_8x1(c1,c2,c3,c4,c5,c6,c7,c8,
                    x1,x2,x3,x4,x5,x6,x7,x8,df,offset+k);
            fma_8x1(c1,c2,c3,c4,c5,c6,c7,c8,
                    f1,f2,f3,f4,f5,f6,f7,c8,dx,offset+k);
            prod_8x1(c1,c2,c3,c4,c5,c6,c7,c8,
                     s1,s2,s3,s4,s5,s6,s7,s8);
            prod_8x1(c1,c2,c3,c4,c5,c6,c7,c8,s,offset+k);
            store_8(c1,c2,c3,c4,c5,c6,c7,c8,a,offset+k);
            shuffle_fwd_8(c1,c2,c3,c4,c5,c6,c7,c8,cx,diag+k);
            t += 32;
         }
         for(int k = 0; k < 128; k += 32){
            
         }
      }
   }     
   printf("iterations: %lu\n",t);
}

 
