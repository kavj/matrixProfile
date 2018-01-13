#include<cstdio>
#include "../arch/avx2arith_2.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
#define vwidth 4
typedef  __m256d vtf;
typedef  __m256i vti;



static inline void fma_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &x1, vtf &x2, vtf &x3, vtf &x4, double *d, int i){
   c1 = fma(x1,aload(d,i),c1);
   c2 = fma(x2,aload(d,i+4),c2);
   c3 = fma(x3,aload(d,i+8),c3);
   c4 = fma(x4,aload(d,i+12),c4);
}


template<typename dtype, typename vt>
static inline void fma_8x1(vt &c1, vt& c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                           vt &x1, vt &x2, vt &x3, vt &x4, vt &x5, vt &x6, vt &x7, vt &x8, dtype *d, int i){

   fma_4x4(c1,c2,c3,c4,x1,x2,x3,x4,d,i);
   fma_4x4(c5,c6,c7,c8,x5,x6,x7,x8,d,i+16);

}


static inline void seti_4x4(vti &c1, vti& c2, vti &c3, vti &c4, 
                            vtf &x1, vtf &x2, vtf &x3, vtf &x4, int t){
   
   vti u = bcast(4);
   vti v = set(t,t+1,t+2,t+3);
   c1 = blend(u,c1,x1);
   u += v;
   c2 = blend(u,c2,x2);
   u += v;
   c3 = blend(u,c3,x3);
   u += v;
   c4 = blend(u,c4,x4);
}



static inline void seti_8x8(vti &c1, vti& c2, vti &c3, vti &c4, vti &c5, vti &c6, vti &c7, vti &c8, 
                            vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8, int t){
   seti_4x4(c1,c2,c3,c4,x1,x2,x3,x4,t);
   seti_4x4(c5,c6,c7,c8,x5,x6,x7,x8,t+16);
}

static inline void cmp_4x4(vtf &c1, vtf &c2, vtf &c3, vtf &c4,
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   c1 = cmpgtr(c1,x1);
   c2 = cmpgtr(c2,x2);
   c3 = cmpgtr(c3,x3);
   c4 = cmpgtr(c4,x4);
}


static inline void cmp_8x8(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   cmp_4x4(c1,c2,c3,c4,x1,x2,x3,x4);
   cmp_4x4(c5,c6,c7,c8,x5,x6,x7,x8);
}


static inline void max_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   c1 = max(c1,x1);
   c2 = max(c2,x2);
   c3 = max(c3,x3);
   c4 = max(c4,x4);
}


static inline void max_8x8(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   max_4x4(c1,c2,c3,c4,x1,x2,x3,x4);
   max_4x4(c5,c6,c7,c8,x5,x6,x7,x8);
}


static inline void prod_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                            vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   x1 = c1*x1;
   x2 = c2*x2;
   x3 = c3*x3;
   x4 = c4*x4;
   x5 = c5*x5;
   x6 = c6*x6;
   x7 = c7*x7;
   x8 = c8*x8;
}


static inline void prod_8x1(vtf &d1, vtf &d2, vtf &d3, vtf &d4, vtf &d5, vtf &d6, vtf &d7, vtf &d8, 
                            vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                            vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   d1 = c1*x1;
   d2 = c2*x2;
   d3 = c3*x3;
   d4 = c4*x4;
   d5 = c5*x5;
   d6 = c6*x6;
   d7 = c7*x7;
   d8 = c8*x8;
}

static inline void prod_8x1(vtf &d1, vtf &d2, vtf &d3, vtf &d4, vtf &d5, vtf &d6, vtf &d7, vtf &d8, 
                            vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   d1 = c1 * aload(d,i);
   d2 = c2 * aload(d,i+4);
   d3 = c3 * aload(d,i+8);
   d4 = c4 * aload(d,i+12);
   d5 = c5 * aload(d,i+16);
   d6 = c6 * aload(d,i+20);
   d7 = c7 * aload(d,i+24);
   d8 = c8 * aload(d,i+28);
}

static inline void prod_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   x1 *= c1;
   x2 *= c2;
   x3 *= c3;
   x4 *= c4;
}

static inline void prod_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, double *d, int i){
  c1 = c1 * aload(d,i);
  c2 = c2 * aload(d,i+4);
  c3 = c3 * aload(d,i+8);
  c4 = c4 * aload(d,i+12);
}


static inline void prod_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   prod_4x4(c1,c2,c3,c4,d,i);
   prod_4x4(c5,c6,c8,c8,d,i+16);
}

template<typename dtype, typename vt>
static inline void load_4x4(vt &c1, vt &c2, vt &c3, vt &c4, dtype *d, int i){
   c1 = aload(d,i);
   c2 = aload(d,i+4);
   c3 = aload(d,i+8);
   c4 = aload(d,i+12);
}

template<typename dtype, typename vt>
static inline void load_8x1(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, dtype *d, int i){
   load_4x4(c1,c2,c3,c4,d,i);
   load_4x4(c1,c2,c3,c4,d,i+16);
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


template<typename dtype, typename vt>
static inline void store_4(vt &c1, vt &c2, vt &c3, vt &c4, dtype *d, int i){
   astore(c1,d,i);
   astore(c2,d,i+4);
   astore(c3,d,i+8);
   astore(c4,d,i+12);
}

template<typename dtype, typename vt>
static inline void store_8(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, dtype *d, int i){
   store_4(c1,c2,c3,c4,d,i);
   store_4(c5,c6,c7,c8,d,i+16);
}


static inline void shuffle_fwd_4(vtf &a, vtf &b, vtf &c, vtf &d, double*e, int i){
   d = c;
   c = b;
   b = a;
   a = aload(d,i);
}


static inline void fma_2(vtf &c1, vtf& c2, vtf &x1, vtf &x2, double *d, int i){
   c1 = fma(x1,aload(d,i),c1);
   c2 = fma(x2,aload(d,i+4),c2);
}

static inline void prod_2(vtf &c1, vtf& c2, vtf &x1, vtf &x2){
   c1 *= x1;
   c2 *= x2;
}

static inline void prod_2(vtf &c1, vtf& c2, double *d, int i){
   c1 *= aload(d,i);
   c2 *= aload(d,i+4);
}

static inline void load_2(vtf &c1, vtf &c2, double *d, int i){
   c1 = aload(d,i);
   c2 = aload(d,i+4);
}

static inline void shuffle_fwd_2(vtf &c1, vtf &c2, double*d, int i){
   c2 = c1;
   c1 = aload(d,i);
}

static inline void store_2(vtf &c1, vtf &c2, double *d, int i){
   astore(c1,d,i);
   astore(c2,d,i+4);
}

// with high register count variations, we should try using memory operands for s on both passes.

static inline void xcorr_kern(double *cx, double *dx, double *df, double *s, double *buf, int diag, int offset, int n, int m){
    vtf ca,cb,cc,cd,ce,cf,cg,ch;
    load_8x1(ca,cb,cc,cd,ce,cf,cg,ch,cx,diag);
    vtf xa, xb, xc, xd, xe, xf, xg, xh,
        fa, fb, fc, fd, fe, ff, fg, fh;
    load_8x1(xa,xb,xc,xd,xe,xf,xg,xh,dx,diag+offset);
    load_8x1(fa,fb,fc,fd,fe,ff,fg,fh,df,diag+offset);
    for(int k = 0; k < 128; k+= 4){
        int p = offset+diag+k;
        int q = offset+k;
        fma_8x1(ca,cb,cc,cd,ce,cf,cg,ch,
                xa,xb,xc,xd,xe,xf,xg,xh,df,q);
        fma_8x1(ca,cb,cc,cd,ce,cf,cg,ch,
                fa,fb,fc,fd,fe,ff,fg,fh,dx,q);
        vtf d1,d2,d3,d4,d5,d6,d7,d8;
        load_8x1(d1,d2,d3,d4,d5,d6,d7,d8,s,q);
        prod_8x1(ca,cb,cc,cd,ce,cf,cg,ch,d1,d2,d3,d4,d5,d6,d7,d8);
        prod_8x1(d1,d2,d3,d4,d5,d6,d7,d8,s,p);
        store_8(d1,d2,d3,d4,d5,d6,d7,d8,buf,offset+8*k);   
        astore(ch,cx,diag+k);
        shuffle_fwd_8(ca,cb,cc,cd,ce,cf,cg,ch,cx,diag+k);
    }
    store_8(ca,cb,cc,cd,ce,cf,cg,ch,cx,diag+128);
}

static void reduce_kern(double *a, double *output, long *outputi, int diag, int offset,int n, int m){
    vtf m1, m2, m3, m4,
        ca, cb, cc, cd;
    vti i1, i2, i3, i4; 
    load_4x4(m1,m2,m3,m4,output,offset);
    load_4x4(i1,i2,i3,i4,outputi,offset);
    for(int k = 0; k < 512; k += 16){
        load_4x4(ca,cb,cc,cd,a,offset+k);
        max_4x4(m1,m2,m3,m4,ca,cb,cc,cd);
        cmp_4x4(ca,cb,cc,cd,m1,m2,m3,m4);
        seti_4x4(i1,i2,i3,i4,ca,cb,cc,cd,offset+8*k);
    }
    store_4(m1,m2,m3,m4,output,offset);
    store_4(i1,i2,i3,i4,outputi,offset);
}

void  accumTest4_7_10(double* cx,double* dx, double* df, double*s, double* a, double* output, long* outputi,int n, int m){
   unsigned long t= 0;
   int offset = 0;
   for(int diag = 0; diag < n-32*m; diag += 32){
      for(int offset = 0; offset < n-32*m; offset += 64){
         xcorr_kern(cx,dx,df,s,a,diag,offset,n,m);
         reduce_kern(a,output,outputi,diag,offset,n,m);
      }
   }     
   printf("iterations: %lu\n",t);
}

