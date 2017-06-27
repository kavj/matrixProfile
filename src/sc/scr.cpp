//#include<ctgmath>
#include<algorithm>
#include<cmath>
#define err 0x40
#define chunkSz 4096  // We have to keep several arrays in L2

struct mpStats{
    const double* T;
    double* mu;
    double* sig;
    int n;
    int m;
};

typedef struct mpStats tsdesc;


/* allocate memory and setup structures */
/* still in debate whether I should overload this for single precision. I need to test stability first */
tsdesc* sc_init(const double* T, int n, int m){
    tsdesc* t = (tsdesc*)malloc(sizeof(tsdesc));
    if(t != NULL){
        t->T = T;
        t->n = n;
        t->m = m;
        t->mu = (double*)malloc(n*sizeof(double));
        t->sig = (double*)malloc(n*sizeof(double));
        if(t->mu == NULL || t->sig == NULL){
            free(t->mu);
            free(t->sig);
            free(t);
            t = NULL;
        }
    }
    return t;
}

/*static inline int alignedSeqLen(int n, int m){
    return (n - m + 1) - (n - m + 1) % m;
}*/


/* This may use a low quality random number generator, so it should be replaced with something more robust*/
/* n is the upper bound of the range of values to be indexed*/
int* buildIndex(int n){
    for(int i = 0; i < n; i++){
            index[i] = i;
        }
    shuffle(index,index+n);
    return index;
}


int iterCnt(int* x, int xOffs, int n, int m, double perc){
    long xLen = n-m+1;
    // this should be replaced by a static helper function
    long numReq = (int)ceil(((double)(xLen*(xLen+1)/2))*perc);
    // Add correction for upper triangular count
    if(numReq + xOffs > xLen){
        numReq = xLen - xOffs;
    }
    // Since segments vary in length, we have to check against both the completion and array length
    int k = 0;
    int i = 0;
    while(i < xLen && k < numReq){
        k += xLen - x[i];
    }
    return k;
}

/* This is based on Higham's approach to the rolling mean and variance, which I think is based on Welford's method for online variance
 * It's easy enough to simultaneously accumulate the means on a pass through the data, so I went ahead and did that. winmean only computes means. 
 * It's orphaned for now
 */
int winmeansig(const double* T, double* mu, double* sig, int n, int m){
    double p = T[0];
    double M = T[0];
    double Q = 0;
    for(int i = 1; i < m; i++){
        p += T[i];
        double t = T[i] - M;
        Q += (i-1)*t*t/i;
        M += t/i;
    }
    mu[0] = p/m;
    for(int i = m; i < alignedLen; i += m){
        if(i > n - m){
            m = n - i;
        }
        double q = p;
        p = 0; 
        for(int j = i; i < j+m; j++){
            q -= T[j-m];
            p += T[j];
            double t = T[i] - M;
            mu[j] = (p + q)/m;
            Q += (i-1)*t*t/i;
            M += t/i;
            sig[i] = sqrt(M/i);
        }
    }
    return 0;
}


/* This should always take aligned commands */
static void sjoincomp(const double* T, const double* mu, const double* sig,  double* mp, int* mpI, int lag, int n, int m){
    double x = 0;
    for(int j = 0; j < m; j++){
        x += T[j]*T[j+lag];
    }
    double z = (x-mu[0]*mu[lag])/(sig[0]*sig[lag]);
    if(z < mp[0]){
        mp[0] = z;
        mpI[0] = lag;
    }
    if(z < mp[lag]){
        mp[lag] = z;
        mpI[lag] = 0;
    }
    for(int i = m; i < n; i+=m){
        if(i > n - m){
            m = n - i;
        }
        double y = x;
        x = 0;
        for(int j = i; j < i+m; j++){
            int k = j + lag;
            y -= T[j-m]*T[k-m];
            x += T[j]  *T[k];
            z = (x + y - mu[j]*mu[k])/(sig[j]*sig[k]);
            if(z < mp[j-m]){
                mp[j-m] = z;
                mpI[j-m] = j-m+lag;
            }
            if(z < mp[j-m+lag]){
                mp[j-m+lag] = z;
                mpI[j-m+lag] = j-m;
            }
        }
    }
}
