#include <immintrin.h>
#define msk1 0x01
#define msk12 0x03
#define msk123 0x07
#ifdef __FMA__            // should work on clang and gcc, update for visual C++ later
#define FMA_SUPPORTED 1
#endif

//a simd class might be prettier, but it leaves too much room for subtle bugs
//This is as close to a minimal interface as I can get. 


// loads and stores
static inline __m256d loada(double* a, int offset){
   return _mm256_load_pd(a+offset);
}

static inline __m256i loada(int* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

static inline __m256d loadu(double* a, int offset){
   return _mm256_loadu_pd(a+offset);
}

static inline __m256i loadu(int* a, int offset){
   return _mm256_loadu_si256((__m256i*)(a+offset));
}

static inline void storea(__m256d oper, double* a, int offset){
   return _mm256_store_pd(a+offset,oper);
}

static inline void storea(__m256i oper, int* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

static inline void storeu(__m256d oper, double* a,int offset){
   _mm256_storeu_pd(a+offset,oper);
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

static inline __m256i set(int i, int j, int k, int l){
   return _mm256_set_epi64x(l,k,j,i);
}

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

// compares
static inline __m256d cmpgtr(__m256d compar, __m256d base){
   return _mm256_cmp_pd(compar,base,_CMP_GT_OS);
}

// swizzles
static inline __m256d blend(__m256d compar, __m256d base, __m256d mask){
   return _mm256_blendv_pd(compar,base,mask);
}

// In the case of indexing, we often set integer values using the result of a floating point comparison. 
static inline __m256i blend(__m256i compar, __m256i base, __m256d mask){
   return _mm256_blendv_epi8(compar,base,_mm256_castpd_si256(mask));
}

static inline __m256d select1(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x01);
}

static inline __m256i select1(__m256i a, __m256i b){
   return _mm256_blend_epi32(a,b,0x01);
}

static inline __m256d select12(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x03);
}

static inline __m256i select12(__m256i a, __m256i b){
   return _mm256_blend_epi32(a,b,0x03);
}

static inline __m256d select123(__m256d a, __m256d b){
   return _mm256_blend_pd(a,b,0x07);
}

static inline __m256i select123(__m256i a, __m256i b){
   return _mm256_blend_epi32(a,b,0x07);
}


