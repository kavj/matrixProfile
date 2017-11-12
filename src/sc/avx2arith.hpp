

// loads and stores

static inline __m256d alignloadD(double* a, int offset){
    return __m256_load_pd(a+offset);
}

static inline __m256d unalignloadD(double* a, int offset){
    return __m256_loadu_pd(a+offset);
}

static inline void alignstoreD(__m256d oper, double* a){
    _mm256_store_pd(a,oper);
}

static inline void unalignstoreD(__m256d oper, double* a){
    _mm256_storeu_pd(a,oper);
}

static inline void alignstoreI(__m256i oper, int* a){

}

static inline __m256d bcast(double* a, int offset){
    return __m256_broadcast_sd(a+offset);
}

// arithmetic
static inline __m256d fma(__m256d a, __m256d b, __m256d c){
   return  _mm256_fmadd_pd(a,b,c);
}

static inline __m256d fms(__m256d a, __m256d b, __m256d c){
    return _mm256_fnmadd_pd(a,b,c);
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

}

static inline __m256d blendD12and34(__m256d a, __m256d b){

}

static inline __m256d blendD123and4(__m256d a, __m256d b){

}

static inline __m256d blendI1and234(__m256i a, __m256i b){

}

static inline __m256d blendI12and34(__m256i a, __m256i b){

}

static inline __m256d blendI123and4(__m256i a, __m256i b){

}

static inline __m256 blendInt(__m256i compar, __m256i base, __m256i mask){
    return _mm256_blendv_epi8(compar,base,mask);
}

static inline __m256d blend1and234(__m256d a, __m256d b){

}

static inline __m256d blend12and34(__m256d a, __m256d b){

}

static inline __m256d blend123and4(__m256d a, __m256d b){

}

// casts 
static inline __m256i d2int(__m256d a){
    return _mm256_castpd_si256(a);
}


static inline __m256d aandb(__m256d a, __m256d mask){
   // _mm256_
}
