#ifdef __AVX__
#ifndef STRIDED
#define STRIDED

#include<immintrin.h>
namespace avx256_t{
#ifdef __GNUC__
#define wrapper __attribute__((always_inline,artificial))
#else
#define always_inline 
#endif


static inline int wrapper testnz(__m256i a){
  _mm256_testnzc_si256(a,a);
}

static inline int wrapper testnz(__m256d a){
   _mm256_testnzc_pd(a,a);
}

static inline void wrapper maskstore(__m256d a, __m256i msk, double* b){
   _mm256_maskstore_pd(b,msk,a);
}

static inline void wrapper maskstore(__m256i a, __m256i msk, int* b){
   _mm256_maskstore_epi32(b,msk,a);
}

// intel intrinsics guide specifies __int64, but using int64_t generates an error using gcc
static inline void wrapper maskstore(__m256i a, __m256i msk, long long* b){
   _mm256_maskstore_epi64(b,msk,a);
}

static inline __m256i wrapper blend(__m256i compar, __m256i base, __m256i mask){
   return _mm256_blendv_epi8(base,compar,mask);
}
/*
static inline __m256d wrapper blend(__m256d compar, __m256d base, __m256i mask){
   return _mm256_blendv_pd(base,compar,mask);
}

*/

static inline __m256d wrapper blend(__m256d compar, __m256d base, __m256d mask){
   return _mm256_blendv_pd(base,compar,mask);
}
/*
// In the case of indexing, we often set integer values using the result of a floating point comparison. 
static inline __m256i wrapper blend(__m256i compar, __m256i base, __m256d mask){
   return _mm256_blendv_epi8(base,compar,_mm256_castpd_si256(mask));
}
*/
static inline __m256d wrapper select1L(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x0E);
}

static inline __m256d wrapper aload(const double *a, int offset){
   return _mm256_load_pd(a+offset);
}

static inline __m256d wrapper uload(const double *a, int offset){
   return _mm256_loadu_pd(a+offset);
}

static inline void wrapper astore(const __m256d &b,double *a, int offset){
   _mm256_store_pd(a+offset,b);
}

static inline void wrapper ustore(const __m256d& b, double* a,int offset){
   _mm256_storeu_pd(a+offset,b);
}

// mul_add needs an ifdef in cases where fma intrinsics won't be loaded

static inline __m256d wrapper mul_add(const __m256d& a,const __m256d &b,const __m256d &c){
   return  _mm256_fmadd_pd(a,b,c);
}

static inline __m256d wrapper mul_sub(const __m256d &a, const __m256d &b,const  __m256d &c){
   return _mm256_fmsub_pd(a,b,c);
}

static inline __m256d wrapper mul_nadd(const __m256d &a, const __m256d &b, const __m256d &c){
   return _mm256_fnmadd_pd(a,b,c);
}

static inline __m256d wrapper vmax(const __m256d &a, const __m256d &b){
   return _mm256_max_pd(a,b);
}

static inline __m256d wrapper vmin(const __m256d &a, const __m256d &b){
   return _mm256_min_pd(a,b);
}

// long long is needed for masked stores on gcc
static inline __m256i wrapper aload(const long long* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i wrapper aload(const int* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i wrapper uload(const long* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline __m256i wrapper uload(const int* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline void wrapper astore(const __m256i &oper, long* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void wrapper astore(const __m256i &oper, int* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void wrapper ustore(const __m256i &oper, long* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline void wrapper ustore(const __m256i &oper, int* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline __m256d wrapper brdcst(const double* a,int offset){
   return _mm256_broadcast_sd(a+offset);
}

static inline __m256i wrapper brdcst(const long a){
   return _mm256_set1_epi64x(a);
}

static inline __m256i wrapper brdcst(const int a){
   return _mm256_set1_epi32(a);
}

static inline __m256d wrapper brdcst(const double a){
   return _mm256_set1_pd(a);
}

static inline __m256i wrapper set(const long i,const long j,const long k,const long l){
   return _mm256_set_epi64x(l,k,j,i);
}

static inline __m256d wrapper setzero(void){
   return _mm256_setzero_pd();
}

}


/*static inline __m256d wrapper mul(const __m256d& a,const  __m256d& b){
   return _mm256_mul_pd(a,b);
}

static inline __m256d wrapper add(const __m256d &a,const __m256d& b){
   return _mm256_add_pd(a,b);
}

static inline __m256i add(__m256i a, __m256i b){
   return _mm256_add_epi64(a,b);
}

static inline __m256d wrapper sub(const __m256d &a, const __m256d &b){
   return _mm256_sub_pd(a,b);
}

static inline __m256d wrapper div(const __m256d &a,const __m256d &b){
   return _mm256_div_pd(a,b);
}

// compares
static inline __m256d wrapper cmpgtr(const __m256d &compar,const __m256d &base){
   return _mm256_cmp_pd(compar,base,_CMP_GT_OS);
}

static inline __m256d wrapper cmplt(const __m256d &compar,const __m256d &base){
   return _mm256_cmp_pd(compar,base,_CMP_LT_OS);
}*/


#endif
#endif
