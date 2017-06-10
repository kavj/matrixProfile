#define _POSIX_C_SOURCE  200112L
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdint.h>

// use mersenne twister file
void shuffle(double* i){

}

// scalar code for windowed dot product, m terms are considered during the first loop, m-1 terms are updated during the second
// Since SIMD versions require some kind of loop peeling, they should be written separately. Otherwise this could be refactored later,
// but it actually expresses the logic quite well. Do not merge it with other functions.
static inline void windot(const double* T1, const double* T2, double* d, const int m){
    T1 = __builtin_assume_aligned(T1,64);
    d  = __builtin_assume_aligned(d,64);
    double accum = 0;
    for(int i = m-1; i >= 0; i--){
        accum += T1[i]*T2[i];
        d[i] = accum;
    }
    __builtin_prefetch(&T1[m]);
    __builtin_prefetch(&T2[m]);
    accum = 0;
    for(int i = m; i < 2*m; i++){
        accum += T1[i]*T2[i];
        d[i-m+1] += accum;
    }
}

static int isunaligned(void* v){
    return (~((uintptr_t) v) & 0x40);
}

static void winsum(const double* T, double* d, const int m){
    double accum = 0;
    for(int i = m-1; i >= 0; i--){
        accum += T[i];
        d[i] = accum;
    }
    accum = 0;
    for(int i = m; i < 2*m; i++){
        accum += T[i];
        d[i]  += accum;
    }
}


/* This should be updated to use masked writebacks  */
static inline void minwb(const double* d, double* mp, int* mpi, const int offset, const int m){
    for(int i = 0; i < m; i++){
        if(d[i] < mp[i]){
            mp[i] = d[i];
            mpi[i] = offset;
        }
    }
}

static inline double onedot(const double* T1, const double* T2, const int m){
    double d = 0;
    for(int i = 0; i < m; i++){
        d += T1[i]*T2[i];
    }
    return d;
}

void winmean(const double* T, double* d, const int m, const int n){
    for(int i = 0; i < n; i+=m){
        windowSum(T,&d[i],m);
        //test if prefetching is needed later
        for(int j = i; j < i+m; j++){
            d[j] /= m;
        }
    }
}

void winstd(const double* T, double* d, const int m){
    
}

void  matrixProfile(const double* T, double* mp, const double n, const double m){
    double* d;
    int err = posix_memalign((void**)&d,64,2*m);
    d = __builtin_assume_aligned(d,64);
    if(err){
        perror("posix_memalign");
        exit(1);
    }
    
    for(int offset = m; offset < m; offset++){
        for(int i = 0; i < n-m-offset; i++){
            windot(&T[i],&T[i+offset],d,m);
            minwb(d,mp,m);
        }
    }
    free(d);
}
