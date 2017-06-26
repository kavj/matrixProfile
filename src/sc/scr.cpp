//#include<ctgmath>
#include<algorithm>
#define err 0x40

struct mpStats{
    const double* T;
    double* mu;
    double* sig;
    int n;
    int m;
};

typedef tsdesc mpStats;


/* allocate memory and setup structures */
/* still in debate whether I should overload this for single precision. I need to test stability first */
tsdesc* sc_init(const double* T, int n, int m){
    tsdesc* t = malloc(sizeof(tsdesc));
    if(t == NULL){
        goto done;
    }
    t->T = T;
    t->n = n;
    t->m = m;
    t->mu = malloc(n*sizeof(double));
    t->sig = malloc(n*sizeof(double));
    if(t->mu == NULL || t->sig == NULL){
        free(t->mu);
        free(t->sig);
        free(t);
        t = NULL;
    }
    done:
    return t;
}

static inline int alignedSeqLen(int n, int m){
    return (n - m + 1) - (n - m + 1) % m;
}



/* This is based on Higham's approach to the rolling mean and variance, which I think is based on Welford's method for online variance
 * It's easy enough to simultaneously accumulate the means on a pass through the data, so I went ahead and did that. winmean only computes means. 
 * It's orphaned for now
 */
int winmeansig(const double* T, double* mu, double* sig, int n, int m){
    double p = 0;
    double M = T[0];
    double Q = 0;
    int alignedLen = alignedSeqLen(n,m);
    if(alignedLen < m){
        return err;
    }
    for(int i = 0; i < m; i++){
        p += T[i];
        
    }
    mu[0] = p/m;
    for(int i = m; i < alignedLen; i += m){
        double q = p;
        p = 0; 
        for(int j = i; i < j+m; j++){
            q -= T[j-m];
            p += T[j];
            mu[j] = (p + q)/m;
            Q += (i-1)*(T[i]-M)/i;
            M += (T[i] - M)/i;
            sig[i] = sqrt(M/i);
            
        }
    }
    for(int i = alignedLen; i < n - m + 1; i++){

    }
    return 0;
}

/*
 * manually inline this
static void builtinShuffle(int* x,int n){
    shuffle(x,x+n);
}*/


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

/*
 * manually inline this
static void sortIndex(int* x, int offs, int k){
    sort(x+offs,x+offs+k);
}*/


/*
void abwindot(const tsdesc* a, const tsdesc* b, int n, int m, int lag){
    double x = 0;
    int alignedLen = alignedSeqLen(n-lag,m);
    for(int j = 0; j < m; j++){
        x += a->T[j];
    } 
    for(int i = m; i < alignedLen; i += m){
        double y = x;
        x = 0; 
        for(int j = i; j < i + m; j++){
            y -= a->T[j-m] * b->T[j-m];
            x += a->T[j]   * b->T[j];
            double z = (x + y - a->mu[j]*b->mu[j])/(a->sig[j]*b->sig[j]);
   
        }
    }
}*/

static inline void updateMP(double* mp, int* mpI, double z, int j, int k){
    if(z < mp[j]){
        mp[j] = z;
        mpI[j] = k;
    }
}

/* Since this algorithm doesn't lend itself to vectorization due to unknown alignment, conditional writebacks with different sized operands, etc, I opted to eliminate as many
 *  memory references as possible 
 *  We can eliminate the branching via a switch from int to long*/

static void sjoinComputePartial(const double* T, const double* mu, const double* sig, double* mp, double* mpi, const int cycles, const int m, const int lag){
    double x = 0;
    int alignedLen = alignedSeqLen(n-lag,m);
    for(int j = 0; j < m; j++){
        x += T[j]*T[j+lag];
    }
    updateMP(mp,mpI,(x - mu[0]*mu[lag])/(sig[0]*sig[lag]));

    for(int i = m; i < alignedLen; i+=m){
        double y = x;
        x = 0;
        for(int j = i; j < i+m; j++){
            int k = j + lag;
            y -= T[j-m]*T[k-m];
            x += T[j]  *T[k];
            double z = (x + y - mu[j]*mu[k])/(sig[j]*sig[k]);
            updateMP(mp,mpI,z,j-m,k);
            updateMP(mp,mpI,z,k-m,j);
        }
    }
}

// Replace arguments with a simpler struct that can be unpacked here. R
void sjoin(const double* T, double* mu, double* sig,  double* mp, int* mpI, int n, int m, int lag){
    double x = 0;
    int alignedLen = alignedSeqLen(n-lag,m);
    sjoinComputePartial(T,mu,sig,mp,mpi,alignedLen/m,m,lag);
    //if(alignedLen < n-m+1){
        sjoinComputePartial(&T[n-m],mu,sig,mp,mpi,2,lag);
    }
}
