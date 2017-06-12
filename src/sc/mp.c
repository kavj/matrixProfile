#define _POSIX_C_SOURCE  200112L
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<stdint.h>

// use mersenne twister file
void shuffle(double* i){

}

/* Computes dot products sum(a[i:i+winLen-1],b[i:i+winLen-1]) for 0 <= i < winLen 
 * Accumulating to a different variable *//

/* This accumulates a series of dot products over a fixed window length. */
void windot(const double* a, const double* b, double* d, const int winLen){
    double accum = 0;
    for(int i = winLen-1; i >= 0; i--){
        accum += a[i]*b[i];
        d[i]   = accum;
    }
    accum = 0;
    for(int i = 1; i < winLen; i++){
        accum += a[i+winLen]*b[i+winLen];
        d[i]  += accum;
    }
}

static inline double onedot(const double* a, const double* b, const int winLen){
    double accum = 0;
    for(int i = 0; i < m; i++){
        accum += a[i]*b[i];
    }
    return accum;
}

/*
static int aligned_256(void* v){
    return (~((uintptr_t) v) & 0x40);
}

static int aligned_512(void* v){
    return 0;
}

static int aligned_128(void* v){
    return 0;
}*/


/* This can be used to compute a windowed mean or other weighting */
static void winscaledsum(const double* a, double* d, const int m){
    double accum = 0;
    for(int i = winLen-1; i >= 0; i--){
        accum += T[i];
        d[i] = accum;
    }
    accum = 0;
    for(int i = 1; i < m; i++){
        accum += a[i+winLen];
        d[i] = (d[i]+accum)/winLen;
    }
}


/* This should be updated to use masked writebacks  */
static inline void minwb(const double* d, double* mp, int* mpi, const int base, const int m){
    
    for(int i = 0; i < m; i++){
        if(d[i] < mp[i]){
            mp[i] = d[i];
            mpi[i] = offset;
        }
    }
}

/* should use welford's method */
void winstd(const double* T, double* d, const int m){
    
}

int tsJoin(const double* a, double* mp, int n, int ssLen){

}

void  matrixProfile(const double* T, double* mp,int n, int  m){
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
