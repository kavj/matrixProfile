#include<cstdio>
#include "../arch/avx2arith_2.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
#define vwidth 4
typedef  __m256d vtf;
typedef  __m256i vti;

// not portabue, which is why it uses raw intrinsics, this should probably be templated by instruction set 
static inline load_4_contig(vtf &ca, vtf &cb, vtf &cc, vtf &cd, double *d,int index){
   d += index;
   ca = _mm256_load_pd(d);
   cb = _mm256_permute2x128_si256(ca,__mm256_load_pd(d+4),33);
   cc = _mm256_alignr_epi8(_mm256_castpd_si256(cb),_mm256_cast_pd(ca),8);
   cb = _mm256_alignr_epi8(_mm256_castpd_si256(cb),_mm256_cast_pd(ca),16);
   cd = _mm256_loadu_pd(d+3);

 _mm256_permute2x128_si256(a,b,33); 
_mm256_load_pd(a+offset);
 _mm256_mul_pd(a,b);
 _mm256_fmadd_pd(a,b,c);

}
static inline __m256i shiftmerge1(__m256i a, __m256i b){
   return _mm256_alignr_epi8(_mm256_permute2x128_si256(a,b,33),a,8);
}

static inline __m256d shiftmerge1(__m256d a, __m256d b){
   return _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(_mm256_permute2f128_pd(a,b,33)),_mm256_castpd_si256(a),8));
}

static inline __m256i shiftmerge2(__m256i a, __m256i b){
   return _mm256_alignr_epi8(_mm256_permute2x128_si256(a,b,43),a,16);
}

static inline __m256d preshuffle(__m256d a, __m256d b){
   return _mm256_permute2f128_pd(a,b,33);
}

static inline __m256i preshuffle(__m256i a, __m256i b){
   return _mm256_permute2x128_si256(a,b,33);
}

static inline __m256d shift1(__m256d lowerIndex, __m256d permuted){
   return _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(permuted),_mm256_castpd_si256(lowerIndex),8));
}

static inline __m256d shift2(__m256d lowerIndex, __m256d permuted){
   return _mm256_castsi256_pd(_mm256_alignr_epi8(_mm256_castpd_si256(permuted),_mm256_castpd_si256(lowerIndex),8));
}

static inline __m256i shift1(__m256i lowerIndex,__m256i permute){
   return _mm256_alignr_epi8(permute,lowerIndex,8);
}

static inline __m256i shift2(__m256i lowerIndex, __m256i permute){
   return _mm256_alignr_epi8(permute,lowerIndex,8);


// fma3 semantics overwrite one register operand, so I only wrote in place versions for these.
static inline void fmadd_4x4(vtf &c1, vtf& c2, vtf &c3, vtf &c4, double* dg, double *dh, int i, int j){
   c1 = fmadd(bcast(dg,i),aload(dh,j),c1);
   c2 = fmadd(bcast(dg,i),aload(dh,j+4),c2);
   c3 = fmadd(bcast(dg,i),aload(dh,j+8),c3);
   c4 = fmadd(bcast(dg,i),aload(dh,j+12),c4);
}


template<typename dtype, typename vt>
static inline void fmadd_8x1(vt &c1, vt& c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                           dtype *dg, dtype *dh, int i, int j){

   fma_4x4(c1,c2,c3,c4,dg,dh,i,j);
   fma_4x4(c5,c6,c7,c8,dg,dh,i+16,j+16);

}



template<typename dtype, typename vt>
static inline void fma_8x1(vt &c1, vt& c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, 
                           vt &x1, vt &x2, vt &x3, vt &x4, vt &x5, vt &x6, vt &x7, vt &x8, dtype *d, int i){

   fma_4x4(c1,c2,c3,c4,x1,x2,x3,x4,d,i);
   fma_4x4(c5,c6,c7,c8,x5,x6,x7,x8,d,i+16);

}

/*
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
}*/

static inline void gt_4ip(vtf &c1, vtf &c2, vtf &c3, vtf &c4,
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   c1 = c1 > x1;
   c2 = c2 > x2;
   c3 = c3 > x3;
   c4 = c4 > x4;
}

static inline void lt_4ip(vtf &c1, vtf &c2, vtf &c3, vtf &c4,
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   c1 = c1 < x1;
   c2 = c2 < x2;
   c3 = c3 < x3;
   c4 = c4 < x4;
}

static inline void gt_8ip(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   gt_4ip(c1,c2,c3,c4,x1,x2,x3,x4);
   gt_4ip(c5,c6,c7,c8,x5,x6,x7,x8);
}

static inline void lt_8ip(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   lt_4ip(c1,c2,c3,c4,x1,x2,x3,x4);
   lt_4ip(c5,c6,c7,c8,x5,x6,x7,x8);
}


static inline void max_4ip(vtf &c1, vtf& c2, vtf &c3, vtf &c4, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   c1 = max(c1,x1);
   c2 = max(c2,x2);
   c3 = max(c3,x3);
   c4 = max(c4,x4);
}


static inline void max_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                           vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){
   max_4x4(c1,c2,c3,c4,x1,x2,x3,x4);
   max_4x4(c5,c6,c7,c8,x5,x6,x7,x8);
}


static inline void prod_4x4_ip(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &x1, vtf &x2, vtf &x3, vtf &x4){
   x1 *= c1;
   x2 *= c2;
   x3 *= c3;
   x4 *= c4;
}

static inline void prod_4x4_ip(vtf &c1, vtf& c2, vtf &c3, vtf &c4, double *d, int i){
  c1 *= aload(d,i);
  c2 *= aload(d,i+4);
  c3 *= aload(d,i+8);
  c4 *=  aload(d,i+12);
}


static inline void prod_8ip(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, 
                            vtf &x1, vtf &x2, vtf &x3, vtf &x4, vtf &x5, vtf &x6, vtf &x7, vtf &x8){

   prod_4ip(c1,c2,c3,c4,x1,x2,x3,x4);
   prod_4ip(c5,c6,c7,c8,x5,x6,x7,x8);
}


static inline void prod_8x1(vtf &c1, vtf& c2, vtf &c3, vtf &c4, vtf &c5, vtf &c6, vtf &c7, vtf &c8, double *d, int i){
   prod_4x4(c1,c2,c3,c4,d,i);
   prod_4x4(c5,c6,c8,c8,d,i+16);
}


// this should just use a stride template parameter with explicit semantics for special cases, such as a stride of 1 with vector types
template<typename dtype, typename vt>
static inline void load_4strided(vt &c1, vt &c2, vt &c3, vt &c4, dtype *d, int i){
   c1 = aload(d,i);
   c2 = aload(d,i+4);
   c3 = aload(d,i+8);
   c4 = aload(d,i+12);
}

template<typename dtype, typename vt>
static inline void load_8strided(vt &c1, vt &c2, vt &c3, vt &c4, vt &c5, vt &c6, vt &c7, vt &c8, dtype *d, int i){
   load_4strided(c1,c2,c3,c4,d,i);
   load_4strided(c1,c2,c3,c4,d,i+16);
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
   a = aload(e,i);
}

