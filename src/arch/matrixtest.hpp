#include<cstdio>
#include "../arch/avx2arith_2.hpp"
using namespace vmth;
//typedef  __m256d vt;
typedef  __m256i vti;

// not portabue, which is why it uses raw intrinsics, this should probably be templated by instruction set 
/*

static inline void load_4(__m256i &ca,__m256i &cb,__m256i &cc,__m256i &cd, long *d,int index){
   d += index;
   ca = _mm256_load_si256((__m256i*)d);
   cb = _mm256_permute2x128_si256(ca,_mm256_load_si256((__m256i*)(d+4)),33);    
   cc = _mm256_alignr_epi8(cb,ca,8);
   cb = _mm256_alignr_epi8(cb,ca,16);
   cd = _mm256_loadu_si256((__m256i*)(d+3));
}

static inline void load_4(__m256d &ca,__m256d &cb,__m256d &cc,__m256d &cd, double *d,int index){
   d += index;
   ca = _mm256_load_pd(d);
   cb = _mm256_permute2f128_pd(ca,_mm256_load_pd(d+4),33);
   cc = _mm256_castsi_mm256_alignr_epi8(_mm256_castpd_si256(cb),_mm256_castpd_si256(ca),8);
   cb = _mm256_alignr_epi8(_mm256_castpd_si256(cb),_mm256_castpd_si256(ca),16);
   cd = _mm256_loadu_pd(d+3);
}
*/


/*
static inline void load_8_contig(vt &ca, vt &cb, vt &cc, vt &cd, vt &ce, vt &cf, vt &cg, vt &ch, double *cx, int index){
  load_4_contig(ca,cb,cc,cd,cx,index);
  load_4_contig(ce,cf,cg,ch,cx,index+4); 
} */


// fma3 semantics overwrite one register operand, so I only wrote in place versions for these.
template<typename vt,int stride>
static inline void fma_4x4(vt &c1, vt& c2, vt &c3, vt &c4, vt &h double* dg, int i, int j){
   c1 = fma(bcast(dg,i),h,c1);
   c2 = fma(bcast(dg,i),h,c2);
   c3 = fma(bcast(dg,i),h,c3);
   c4 = fma(bcast(dg,i),h,c4);
}

template<typename dtype, typename vt>
static inline void fma_8x1(vt &c1, vt& c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                           dtype *dg, dtype *dh, int i, int j){

   fma_4x4(c1,c2,c3,c4,dg,dh,i,j);
   fma_4x4(c5,c6,c7,c8,dg,dh,i+4*stride,j+4*stride);

}

template<typename dtype, typename vt>
static inline void fma_8x1(vt &c1, vt& c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                           vt &x1, vt &x2, vt &x3, vt &x4, vt &x5, vt &x6, vt &x7, vt &x8, dtype *d, int i){

   fma_4x4(c1,c2,c3,c4,x1,x2,x3,x4,d,i);
   fma_4x4(c5,c6,c7,c8,x5,x6,x7,x8,d,i+4*stride);

}

/*
template<typename vt, typename vti>
static inline void seti_4x4(vti &ca, vti& cb, vti &cc, vti &cd, 
                            vt &xa, vt &xb, vt &xc, vt &xd, int t){
   
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
*/
template<typename vt,typename vti>
static inline void seti_8x8(vti &ca, vti& cb, vti &cc, vti &cd, vti &ce, vti &cf, vti &cg, vti &ch, 
                            vt &xa, vt &xb, vt &xc, vt &xd, vt &xe, vt &xf, vt &xg, vt &xh, int t){
   seti_4x4(ca,cb,cc,cd,xa,xb,xc,xd,t);
   seti_4x4(ce,cf,cg,ch,xe,xf,xg,xh,t+4*stride);
}

template<typename vt>
static inline void gt_4ip(vt &c1, vt &c2, vt &c3, vt &c4,
                           vt &x1, vt &x2, vt &x3, vt &x4){
   c1 = c1 > x1;
   c2 = c2 > x2;
   c3 = c3 > x3;
   c4 = c4 > x4;
}

template<typename vt>
static inline void lt_4ip(vt &c1, vt &c2, vt &c3, vt &c4,
                           vt &x1, vt &x2, vt &x3, vt &x4){
   c1 = c1 < x1;
   c2 = c2 < x2;
   c3 = c3 < x3;
   c4 = c4 < x4;
}

template<typename vt>
static inline void gt_8ip(vt &c1, vt& c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                          vt &x1, vt &x2, vt &x3, vt &x4, vt &x5, vt &x6, vt &x7, vt &x8){
   gt_4ip(c1,c2,c3,c4,x1,x2,x3,x4);
   gt_4ip(c5,c6,c7,c8,x5,x6,x7,x8);
}

template<typename vt>
static inline void lt_8ip(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                          vt &x1, vt &x2, vt &x3, vt &x4, vt &x5, vt &x6, vt &x7, vt &x8){
   lt_4ip(c1,c2,c3,c4,x1,x2,x3,x4);
   lt_4ip(c5,c6,c7,c8,x5,x6,x7,x8);
}


template<typename vt>
static inline void max_4ip(vt &ca, vt &cb, vt &cc, vt &cd, 
                           vt &xa, vt &xb, vt &xc, vt &xd){
   ca = max(ca,xa);
   cb = max(cb,xb);
   cc = max(cc,xc);
   cd = max(cd,xd);
}

template<typename vt>
static inline void max_8ip(vt &ca, vt& cb, vt &cc, vt &cd, vt &ce, vt &cf, vt &cg, vt &ch, 
                           vt &xa, vt &xb, vt &xc, vt &xd, vt &xe, vt &xf, vt &xg, vt &xh){
   max_4x4(ca,cb,cc,cd,xa,xb,xc,xd);
   max_4x4(ce,cf,cg,ch,xe,xf,xg,xh);
}

template<typename vt>
static inline void prod_4ip(vt &ca, vt& cb, vt &cc, vt &cd, vt &xa, vt &xb, vt &xc, vt &xd){
   ca *= xc;
   cb *= xb;
   cc *= xc;
   cd *= xd;
}

template<typename vt>
static inline void prod_4(vt &a1, vt &a2, vt &a3, vt &a4, vt &c1, vt &c2, vt &c3, vt &c4, vt &x1, vt &x2, vt &x3, vt &x4){
   a1 = x1 * c1;
   a2 = x2 * c2;
   a3 = x3 * c3;
   a4 = x4 * c4;
}

template <typename vt, int stride>
static inline void prod_4ip(vt &c1, vt &c2, vt &c3, vt &c4, double *d, int i){
  c1 *= aload(d,i);
  c2 *= aload(d,i+stride);
  c3 *= aload(d,i+2*stride);
  c4 *= aload(d,i+3*stride);
}

template <typename vt, int stride>
static inline void uprod_4(vt &a1, vt &a2, vt &a3, vt &a4, vt &c1, vt &c2, vt &c3, vt &c4, double *d, int i){

  a1 = c1 * uload(d,i);
  a2 = c2 * uload(d,i+stride);
  a3 = c3 * uload(d,i+2*stride);
  a4 = c4 * uload(d,i+3*stride);
}


template <typename vt,int stride>
static inline void uprod_4ip(vt &c1, vt &c2, vt &c3, vt &c4, double *d, int i){
  c1 *= uload(d,i);
  c2 *= uload(d,i+stride);
  c3 *= uload(d,i+2*stride);
  c4 *= uload(d,i+3*stride);
}

template <typename vt>
static inline void prod_8x1(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, double *d, int i){
   prod_4ip(c1,c2,c3,c4,d,i);
   prod_4ip(c5,c6,c8,c8,d,i+4*stride);
}


// this should just use a stride template parameter with explicit semantics for special cases, such as a stride of 1 with vector types
template<typename dtype, typename vt, int stride>
static inline void load_4strided(vt &c1, vt &c2, vt &c3, vt &c4, dtype *d, int i){
   c1 = aload(d,i);
   c2 = aload(d,i+stride);
   c3 = aload(d,i+2*stride);
   c4 = aload(d,i+3*stride);
}

template<typename dtype, typename vt, int stride>
static inline void load_8strided(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, dtype *d, int i){
   load_4strided<dtype,vt,stride>(c1,c2,c3,c4,d,i);
   load_4strided<dtype,vt,stride>(c1,c2,c3,c4,d,i+4*stride);
}

template<typename dtype, typename vt, int stride>
static inline void store_4(vt &c1, vt &c2, vt &c3, vt &c4, dtype *d, int i){
   astore(c1,d,i);
   astore(c2,d,i+stride);
   astore(c3,d,i+2*stride);
   astore(c4,d,i+3*stride);
}

template<typename dtype, typename vt, int stride>
static inline void store_8(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, dtype *d, int i){
   store_4(c1,c2,c3,c4,d,i);
   store_4(c5,c6,c7,c8,d,i+4*stride);
}


