#include<ctgmath>
#include<algorithm>

#define err 0x40

struct mpStats{
    double* T;
    double* mu;
    double* sig;
    int n;
    int m; 
};



/* allocate memory and setup structures */
void sc_init(int n, int m){
    

}

static inline int alignedSeqLen(int n, int m){
    return (n - m + 1) - (n - m + 1) % m;
}


/* Since we need to compute the means anyway, the two pass method is probably a better option. */
void winmean(const double* T, double* mu, int n, int m){
    int aligned = alignedSeqLen(n,m);
    double x = 0;
    for(int i = 0; i < m; i++){
        x += T[i];
    }
    mu[0] = x/m;
    for(int i = m; i < n-m+1; i+=m){
        double y = x;
        x = 0;
        for(int j = i; j < i+m; j++){
            y -= T[j-m];
            x += T[j];
            mu[j] = (x+y)/m;
        }
    }
}



/* This is based on Higham's approach to the rolling mean and variance, which I think is based on Welford's method for online variance
 * It's easy enough to simultaneously accumulate the means on a pass through the data, so I went ahead and did that. winmean only computes means. 
 * It's orphaned for now
 */
int rollmeanst(const double* T, double* mu, double* sig, int n, int m){
    double p = 0;
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
        double M = T[0];
        double Q = 0;
        for(int j = i; i < j+m; j++){
            q -= T[j-m];
            p += T[j];
            
            Q += (i-1)*(T[i]-M)/i;
            //
            M += (T[i] - M)/i;
            sig[i] = sqrt(M/i);
            // this will be something like sig[i]
        }
    }
    return 0;
}

/*
 * manually inline this
static void builtinShuffle(int* x,int n){
    shuffle(x,x+n);
}*/


int numItersReq(int* x, int xOffs, int n, int m, double perc){
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


void windotAB(ts* a, ts* b, int m){
    double x = 0;
    int len = a->len > b->len ? a->len : b->len;
    aligned = alignedSeqLen(len,m);

    for(int j = 0; j < m; j++){
        
    }
    for(int i = 0; i < aligned; i+=m){
        
    }
    for(int j = 0; j < len; j++){
    
    }
}


/* Since this algorithm doesn't lend itself to vectorization due to unknown alignment, conditional writebacks with different sized operands, etc, I opted to eliminate as many
 *  memory references as possible */
void windotself(const double* T, const double* mu, const double* sig, double* mp, int* mpI, const int n, const int m, int lag){
    double x = 0;
    int aligned = alignedSeqLen(n-lag,m);
    for(int j = 0; j < m; j++){

    }
    for(int i = m; i < aligned; i+=m){
        double y = x;
        x = 0;
        for(int j = i; j < i+m; j++){
            x += T[j]*T[j+lag];
            y -= T[j-m]*T[j-m+lag];
            double z = (x + y - mu[j]*mu[j+lag])/(std[j]*sig[j+lag]);
            // If we change the matrix profile index from int to long, we could eliminate the use of branching here. Typically the algorithm would have more writes early on,
            // so it may help for an anytime algorithm
            if(z <  mp[j]){
                mp[j] = z;
                mpI[j] = j+lag;
            }
            if(z < mp[i+lag]){
                mp[j+lag] = z;
                mpI[j+lag] = j;
            }
        }
    }
    for(int j = aligned; j < (n - lag ); j++){

    }
}

