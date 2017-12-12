#include <immintrin.h>

#ifdef __FMA__            // should work on clang and gcc, update for visual C++ later
#define FMA_SUPPORTED 1
#endif

//a simd class might be prettier, but it leaves too much room for subtle bugs
//This is as close to a minimal interface as I can get. 

namespace vmth{
// loads and stores
static inline __m256d loada(double* a, int offset){
   return _mm256_load_pd(a+offset);
}

static inline __m256i loada(long* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256i loada(int* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256d loadu(double* a, int offset){
   return _mm256_loadu_pd(a+offset);
}

static inline __m256i loadu(long* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline __m256i loadu(int* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

static inline void storea(__m256d oper, double* a, int offset){
   return _mm256_store_pd(a+offset,oper);
}

static inline void storea(__m256i oper, long* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void storea(__m256i oper, int* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void storeu(__m256d oper, double* a,int offset){
   _mm256_storeu_pd(a+offset,oper);
}

static inline void storeu(__m256i oper, long* a,int offset){
   _mm256_storeu_si256((__m256i*)a+offset,oper);
}

static inline void storeu(__m256i oper, int* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

static inline __m256d bcast(double* a,int offset){
   return _mm256_broadcast_sd(a+offset);
}


static inline __m256i bcast(long a){
   return _mm256_set1_epi64x(a);
}

static inline __m256i bcast(int a){
   return _mm256_set1_epi32(a);
}

static inline __m256i set(long i, long j, long k, long l){
   return _mm256_set_epi64x(l,k,j,i);
}

static inline __m256d setzero(void){
   return _mm256_setzero_pd();
}

// temporary
/*static inline __m256d set(double i, double j, double k, double l){
   return (__m256d)  _mm256_set_epi64x(l,k,j,i);
}*/


// arithmetic
//
// prefer fused multiply add whenever possible. I might need to adjust flags for other compilers
#ifdef FMA_SUPPORTED
static inline __m256d mult_add(__m256d a, __m256d b, __m256d c){
   return  _mm256_fmadd_pd(a,b,c);
}

static inline __m256d mult_sub(__m256d a, __m256d b, __m256d c){
   return _mm256_fnmadd_pd(a,b,c);
}

#else
static inline __m256d mult_add(__m256d a, __m256d b, __m256d c){
   return  _mm256_add_pd(_mm256_mul_pd(a,b),c);
}

static inline __m256d mult_sub(__m256d a, __m256d b, __m256d c){
   return _mm256_sub_pd(c,_mm256_mul_pd(a,b));
}
#endif

static inline __m256d mult(__m256d a, __m256d b){
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

// compares
static inline __m256d cmpgtr(__m256d compar, __m256d base){
   return _mm256_cmp_pd(compar,base,_CMP_GT_OS);
}

// swizzles
static inline __m256d blend(__m256d compar, __m256d base, __m256d mask){
   return _mm256_blendv_pd(base,compar,mask);
}

// In the case of indexing, we often set integer values using the result of a floating point comparison. 
static inline __m256i blend(__m256i compar, __m256i base, __m256d mask){
   return _mm256_blendv_epi8(base,compar,_mm256_castpd_si256(mask));
}

static inline __m256d select1(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x0E);
}

static inline __m256i select1(__m256i a, __m256i b){
   return _mm256_blend_epi32(a,b,0x0E);
}

static inline __m256d select12(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x0C);
}

static inline __m256i select12(__m256i a, __m256i b){
   return _mm256_blend_epi32(a,b,0x0C);
}

static inline __m256d select123(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x08);
}

static inline __m256i select123(__m256i a, __m256i b){
   return _mm256_blend_epi32(a,b,0x08);
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
}

/*
static inline __m256d shiftOnly2(__m256d lowerIndex, __m256d permuted){
   return _mm256_alignr_epi8(permuted,lowerIndex,16);
}*/


/*
static inline __m256i packIntMask(__m256d a, __m256d b){
   
}*/

static inline __m256d max(__m256d a, __m256d b){
   return _mm256_max_pd(a,b);
}

// should go in a scalar header, here for now
static inline double max(double a, double b){
   return a > b ? a : b;
}

static inline __m256d min(__m256d a, __m256d b){
   return _mm256_min_pd(a,b);
}

}
