#include<immintrin.h>
#ifdef __AVX__
#ifndef STRIDED
#define STRIDED


template<typename __m256d, int stride>
static inline void strided_store_x8(__m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4, __m256d &r5, __m256d &r6, __m256d &r7, __m256d &r8, double * dat, int offset) __attribute__ ((always_inline));

template<typename __m256d, int stride>
static inline void strided_load_x8(__m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4, __m256d &r5, __m256d &r6, __m256d &r7, __m256d &r8, double* dat, int offset) __attribute__ ((always_inline)); 

template<typename __m256d, int stride>
static inline void strided_load_x8(__m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4, __m256d &r5, __m256d &r6, __m256d &r7, __m256d &r8, double* dat, int offset){
   
   r1 = _mm256_load_pd(dat+offset);
   r2 = _mm256_load_pd(dat+offset+stride);
   r3 = _mm256_load_pd(dat+offset+2*stride);
   r4 = _mm256_load_pd(dat+offset+3*stride);
   r5 = _mm256_load_pd(dat+offset+4*stride);
   r6 = _mm256_load_pd(dat+offset+5*stride);
   r7 = _mm256_load_pd(dat+offset+6*stride);
   r8 = _mm256_load_pd(dat+offset+7*stride);
}

template<typename __m256d, int stride>
static inline void strided_store_x8(__m256d &r1, __m256d &r2, __m256d &r3, __m256d &r4, __m256d &r5, __m256d &r6, __m256d &r7, __m256d &r8, double * dat, int offset){
   _mm256_store_pd(dat+offset,r1);
   _mm256_store_pd(dat+offset+stride,r2);
   _mm256_store_pd(dat+offset+2*stride,r3);
   _mm256_store_pd(dat+offset+3*stride,r4);
   _mm256_store_pd(dat+offset+4*stride,r5);
   _mm256_store_pd(dat+offset+5*stride,r6);
   _mm256_store_pd(dat+offset+6*stride,r7);
   _mm256_store_pd(dat+offset+7*stride,r8);
}

template<typename __m256d>
static inline __m256d aload(double *a, int offset){
   return _mm256_load_pd(a+offset);
}

template<typename __m256d>
static inline __m256d uload(double *a, int offset){
   return _mm256_loadu_pd(a+offset);
}

template<typename __m256d> 
static inline void astore(__m256d b, double *a, int offset){
   _mm256_store_pd(a+offset,b);
}

template<typename __m256d>
static inline void ustore(__m256d b, double* a,int offset){
   _mm256_storeu_pd(a+offset,b);
}

template<typename __m256d>
static inline __m256d vfma(__m256d a, __m256d b, __m256d c){
   return  _mm256_fmadd_pd(a,b,c);
}

template<typename __m256d>
static inline __m256d vfms(__m256d a, __m256d b, __m256d c){
   return _mm256_fmsub_pd(a,b,c);
}

template<typename __m256d>
static inline __m256d vfmna(__m256d a, __m256d b, __m256d c){
   return _mm256_fnmadd_pd(a,b,c);
}

template<typename __m256i, long>
static inline __m256i aload(long* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

template<typename __m256i, int>
static inline __m256i aload(int* a, int offset){
   return  _mm256_load_si256((__m256i*)(a+offset));
}

template<typename __m256i, long>
static inline __m256i uload(long* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

template<typename __m256i, int>
static inline __m256i uload(int* a, int offset){
   return _mm256_loadu_si256((__m256i const*)(a+offset));
}

template<typename __m256i,long>
static inline void astore(__m256i oper, long* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

template<typename __m256i,int>
static inline void astore(__m256i oper, int* a, int offset){
   _mm256_store_si256((__m256i*)(a+offset),oper);
}

template<typename __m256i,long>
static inline void ustore(__m256i oper, long* a,int offset){
   _mm256_storeu_si256((__m256i*)a+offset,oper);
}

template<typename __m256i,int>
static inline void ustore(__m256i oper, int* a,int offset){
   _mm256_storeu_si256((__m256i*)(a+offset),oper);
}

#endif
#endif
//static inline __m256d bcast(double* a,int offset){