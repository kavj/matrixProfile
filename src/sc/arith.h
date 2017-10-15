


/* This section should be cleaned up. It is supposed to check whether AVX is appropriate. It may be desirable to have 
 * runtime checks for this.
 */

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

#ifdef __AVX2__
//#ifdef 
static inline __m256d shiftMean(__m256d mu, __m256d xi, __m256d xe, __m256i m){

}

static inline __m256d shiftSX(__m256d s, __m256d mu){

}

static inline __m256d decSXY(__m256d c, __m256d x, __m256d y){
    return _mm_fnmadd_pd(x,y,c);
}

static inline __m256d incSXY(__m256d c, __m256d x, __m256d y){
    return _mm_fmadd_pd(x,y,c);
}

#else

static inline double shiftMean(double mu, double a, double b, double m){
    return mu + (a-b)/m;
}

static inline double shiftSSS(double s, double mu, double muprev, double x, double xprev){
    return s + ((x - mu) * (x - muprev) - (xprev - mu) * (xprev - muprev));
}

static inline double shiftSXY(double x, double y, double xprev, double yprev, double mux, double muyprev){
    return (x - mux)*(y - muyprev) - (xprev - mux)*(yprev - muyprev);
}

#endif


/*static inline double initSS(const double* T, double M, int m, int offset){
    double s = 0;
    for(int i = offset; i < offset + m; i++){
        s += (T[i] - M) * (T[i] - M);
    }
    return s;
}*/


