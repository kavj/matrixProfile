


/* This section should be cleaned up. It is supposed to check whether AVX is appropriate. It may be desirable to have 
 * runtime checks for this.
 */


/*static inline double initSS(const double* T, double M, int m, int offset){
    double s = 0;
    for(int i = offset; i < offset + m; i++){
        s += (T[i] - M) * (T[i] - M);
    }
    return s;
}*/



/*
 * An anytime version would need some way to determine which lag indices to run. I'll use this later if necessary.
 * */




/*int iterCnt(int* x, int xOffs, int n, int m, double perc){
    int w = n - m + 1;
    long numReq = (int)ceil(((double)(w*(w+1)/2))*perc);
    if(numReq + xOffs > w){
        numReq = w - xOffs;
    }
    int k = 0;
    int i = 0;
    while(i < w && k < numReq){
        k += w - x[i];
    }
    return k;
}*/


static inline double shiftMean(double mu, double xi, double xe, double m){
    return mu + (xi-xe)/m;
}

/* This section just organizes arithmetic operations using appropriate SIMD extensions if available
 * this way the source code for the basic algorithms don't directly commit to a particular SIMD width, and their 
 * source is still readable
 */

#ifdef __AVX2__

static inline __m128d decSX(__m128d xe, __m128d s){
    return _mm_fnmadd_sd(xe,xe,s);
}

static inline __128d incSX(__m128d xi, __m128d s){
    return _mm_fmadd_sd(xi,xi,s);
}
 
static inline __m256d decSXY(__m256d x, __m256d y, __m256d c){
    return _mm_fnmadd_pd(x,y,c);
}

static inline __m256d incSXY(__m256d x, __m256d y, __m256d c){
    return _mm_fmadd_pd(x,y,c);
}

#else

static inline double decSX(double xe, double s){
    return s - xe*xe;
}

static inline double incSX(double xi, double s){
    return s + xi*xi;
}

static inline double shiftSXY(double x, double y, double xprev, double yprev, double mux, double muyprev){
    return (x - mux)*(y - muyprev) - (xprev - mux)*(yprev - muyprev);
}

#endif



