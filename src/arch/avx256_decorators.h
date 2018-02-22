#include<immintrin.h>
#ifdef __AVX__
#ifndef STRIDED
#define STRIDED

namespace avx256_t{

#define stride 4
// double precision
//
static inline __m256d aload(double *a, int offset){
   return _mm256_load_pd(a+offset);
}

static inline __m256d uload(double *a, int offset){
   return _mm256_loadu_pd(a+offset);
}

static inline void astore(__m256d b, double *a, int offset){
   _mm256_store_pd(a+offset,b);
}

static inline void ustore(__m256d b, double* a,int offset){
   _mm256_storeu_pd(a+offset,b);
}

static inline __m256d mul(__m256d a, __m256d b){
   return _mm256_mul_pd(a,b);
}

static inline __m256d add(__m256d a, __m256d b){
   return _mm256_add_pd(a,b);
}

static inline __m256i add(__m256i a, __m256i b){
   return _mm256_add_epi64(a,b);
}

static inline __m256d sub(__m256d a, __m256d b){
   return _mm256_sub_pd(a,b);
}

static inline __m256d div(__m256d a, __m256d b){
   return _mm256_div_pd(a,b);
}

// compares
static inline __m256d cmpgtr(__m256d compar, __m256d base){
   return _mm256_cmp_pd(compar,base,_CMP_GT_OS);
}

static inline __m256d cmplt(__m256d compar, __m256d base){
   return _mm256_cmp_pd(compar,base,_CMP_LT_OS);
}

static inline __m256d mul_add(__m256d a, __m256d b, __m256d c){
   return  _mm256_fmadd_pd(a,b,c);
}

static inline __m256d mul_sub(__m256d a, __m256d b, __m256d c){
   return _mm256_fmsub_pd(a,b,c);
}

static inline __m256d mul_nadd(__m256d a, __m256d b, __m256d c){
   return _mm256_fnmadd_pd(a,b,c);
}

static inline __m256d vmax(__m256d a, __m256d b){
   return _mm256_max_pd(a,b);
}

static inline __m256d vmin(__m256d a, __m256d b){
   return _mm256_min_pd(a,b);
}

static inline __m256i aload(long* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i aload(int* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i uload(long* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline __m256i uload(int* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline void astore(__m256i oper, long* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void astore(__m256i oper, int* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void ustore(__m256i oper, long* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline void ustore(__m256i oper, int* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline __m256d brdcst(double* a,int offset){
   return _mm256_broadcast_sd(a+offset);
}

static inline __m256i brdcst(long a){
   return _mm256_set1_epi64x(a);
}

/*static inline __m256i brdcst(int a){
   return _mm256_set1_epi32(a);
}*/

static inline __m256d brdcst(double a){
   return _mm256_set1_pd(a);
}

static inline __m256i set(long i, long j, long k, long l){
   return _mm256_set_epi64x(l,k,j,i);
}

static inline __m256d setzero(void){
   return _mm256_setzero_pd();
}
/*
static inline __m128d extractlo(__m256d a){
   return _mm256_extract128_pd(a,0);
}*/
/*
static inline __m256d testnz(__m256d a){

}
*/



// SINGLE PRECISION


static inline __m256 aload(float* a, int offset){
   return _mm256_load_ps(a+offset);
}

static inline __m256 uload(float* a, int offset){
   return _mm256_loadu_ps(a+offset);
}

static inline void astore(__m256 b, float* a, int offset){
   _mm256_store_ps(a+offset,b);
}

static inline void ustore(__m256 , float* a,int offset){
   _mm256_storeu_ps(a+offset,b);
}

static inline __m256 mul(__m256 a, __m256 b){
   return _mm256_mul_ps(a,b);
}

static inline __m256d add(__m256 a, __m256d b){
   return _mm256_add_ps(a,b);
}

static inline __m256i add(__m256i a, __m256i b){
   return _mm256_add_epi32(a,b);
}

static inline __m256d sub(__m256 a, __m256 b){
   return _mm256_sub_ps(a,b);
}

static inline __m256d div(__m256 a, __m256 b){
   return _mm256_div_ps(a,b);
}

// compares
static inline __m256d cmpgtr(__m256 compar, __m256 base){
   return _mm256_cmp_ps(compar,base,_CMP_GT_OS);
}

static inline __m256d cmplt(__m256 compar, __m256 base){
   return _mm256_cmp_ps(compar,base,_CMP_LT_OS);
}

static inline __m256d mul_add(__m256 a, __m256 b, __m256 c){
   return  _mm256_fmadd_ps(a,b,c);
}

static inline __m256d mul_sub(__m256 a, __m256 b, __m256 c){
   return _mm256_fmsub_ps(a,b,c);
}

static inline __m256d mul_nadd(__m256 a, __m256 b, __m256 c){
   return _mm256_fnmadd_ps(a,b,c);
}

static inline __m256d vmax(__m256 a, __m256 b){
   return _mm256_max_ps(a,b);
}

static inline __m256d vmin(__m256d a, __m256 b){
   return _mm256_min_ps(a,b);
}

static inline __m256i aload(long* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i aload(int* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i uload(long* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline __m256i uload(int* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline void astore(__m256i oper, long* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void astore(__m256i oper, int* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void ustore(__m256i oper, long* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline void ustore(__m256i oper, int* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline __m256d brdcst(double* a,int offset){
   return _mm256_broadcast_sd(a+offset);
}

/*static inline __m256i brdcst(long a){
   return _mm256_set1_epi64x(a);
}*/

static inline __m256i brdcst(int a){
   return _mm256_set1_epi32(a);
}

static inline __m256d brdcst(float a){
   return _mm256_set1_pd(a);
}

static inline __m256i set(int i,int j,int k,int l){ // should we use stdint here?
   return _mm256_set_epi64x(l,k,j,i);
}


}


#endif
#endif
