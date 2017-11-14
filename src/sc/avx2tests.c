#include "avx2arith.h"

// It would be nice to cast everything to __m256d and use arithmetic operators where possible. 
// Unfortunately certain compilers (particularly visual C++) do not implement arithmetic operators for vector types
// Given that we're after fairly strict logic and that we avoid heavy reliance on cross lane permutations and masked loads and stores
// it should be sufficient to just include the appropriate implementation and data type definitions at build time, using a generic simd data type
// We're using a separate scalar kernel to avoid cases where the compiler might inappropriately use branch instructions over conditional moves 
// (which have been notably faster here).

// loads and stores

void printDArray(double *a){
   printf("%lf %lf %lf %lf\n",a[0],a[1],a[2],a[3]);
}

void printIntArray(int *a){
   printf("%d %d %d %d\n",a[0],a[1],a[2],a[3]);
}

void printm256I(__m256i a){
   int b[4];
   __m256_store_pd(b,a);
   printf("%d %d %d %d\n",b[0],b[1],b[2],b[3]);
}

void printm256D(__m256d a){
   double b[4];
   __m256_store_pd(b,a);
   printf("%lf %lf %lf %lf\n",b[0],b[1],b[2],b[3]);
}

void makeshiftmultadd(double* a, double* b, double* c){
    for(int i = 0; i < 3; i++){
       c[i] = a[i]*b[i] + a[i];
    }
}

void makeshiftmultsub(double* a, double* b, double* c){
    for(int i = 0; i < 3; i++){
       c[i] = a[i] - a[i]*b[i];
    } 
}

void expect(void){
   printf("expected output\n");
}

void actual(void){
   printf("actual output\n");
}

int main(void){
   double a[4] = {1.0,5.2,3.0,7.0};
   double b[4] = {1.2,2.0,10.0,11.0};
   double f[4];
   int h[4] = {1,2,3,4};
   __m256d c = alignloadD(a);
   __m256d d = unalignloadD(b);
   
   expect();
   printDArray(a);
   printDArray(b);
   actual();
   printm256D(c);
   printm256D(d);
   expect();
   printDArray(a);
   printDArray(b);
   alignstoreD(c,a);
   unalignstoreD(d,b);
   actual();
   printDArray(a);
   printDArray(b);
   
   __m256d g = bcastD(a[1]);
   __m256i j = bcastI(h[2]);
   printf("testing broadcast");
   printf("expecting %lf\n",a[1]);
   printf("expecting %lu\n",h[2]);
   print256D(g);
   print256I(j);

   __m256d k = multadd(c,d,c);
   __m256d m = multsub(c,d,c);
   double buf[4];
   makeshiftmultadd(a,b,buf);
   expect();
   printDArray(buf);
   makeshiftmultsub(a,b,buf);
   printDArray(buf);
   actual();
   print256D(k);
   print256D(m);
   
   
}



static inline __m256i setInts(int i, int j, int k, int l){
    return _mm256_set_epi64x(l,k,j,i);
}

// arithmetic
//
// prefer fused multiply add whenever possible. I might need to adjust flags for other compilers

static inline __m256d mult(__m256d a, __m256d b){
    return _mm256_mul_pd(a,b);
}

static inline __m256i add(__m256i a, __m256i b){
    return _mm256_add_epi64(a,b);
}

// compares
static inline __m256d cmpgt(__m256d compar, __m256d base){
    return _mm256_cmp_pd(compar,base,_CMP_GT_OS);
}

// swizzles
static inline __m256d blend(__m256d compar, __m256d base, __m256d mask){
    return _mm256_blendv_pd(compar,base,mask);
}

static inline __m256d blendD1and234(__m256d a, __m256d b){
    return _mm256_blend_pd(a,b,msk1);
}

static inline __m256d blendD12and34(__m256d a, __m256d b){
    return _mm256_blend_pd(a,b,msk12);
}

static inline __m256d blendD123and4(__m256d a, __m256d b){
    return _mm256_blend_pd(a,b,msk123);
}

static inline __m256d blendI1and234(__m256i a, __m256i b){
    return _mm256_blend_pd(a,b,msk1);
}

static inline __m256d blendI12and34(__m256i a, __m256i b){
    return _mm256_blend_pd(a,b,msk12);
}

static inline __m256d blendI123and4(__m256i a, __m256i b){
    return _mm256_blend_pd(a,b,msk123);
}

static inline __m256 blendInt(__m256i compar, __m256i base, __m256i mask){
    return _mm256_blendv_epi8(compar,base,mask);
}

// casts 
static inline __m256i d2int(__m256d a){
    return _mm256_castpd_si256(a);
}

