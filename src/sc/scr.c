//#include<cstdlib>
//#include<cmath>
#include<stdlib.h>
#include<math.h>
//#include<random>
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



/*
 * We can add a correction to get the count considering diagonals only over the upper triangular portion.
 * We need a global bound. 
 * 
 */
int iterCnt(int* x, int xOffs, int n, int m, double perc){
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
}


static double decM(double M, double x, int k){
    return (k*M-x)/(k-1);
}

static double incM(double M, double x, int k){
    return M + (x-M)/k;
}

static double decQ(double Q, double M, double x, int k){
    return Q - (x-M)*(x-M)/k;
}

static double incQ(double Q, double M, double x, int k){
    return Q + (x-M)*(x-M)/k;
}


/* This is based on Higham's approach to the rolling mean and variance, which I think is based on Welford's method for online variance
 * It's easy enough to simultaneously accumulate the means on a pass through the data, so I went ahead and did that. 
 */
int winmeansig(const double* T, double* mu, double* sigma, int n, int m){
    double M = T[0];
    double Q = 0;
    for(int i = 1; i < m; i++){
        Q = incQ(Q,M,T[i],i+1);
        M = incM(M,T[i],i+1);
    }
    mu[0] = M;
    sigma[0] = sqrt(Q/m);
    for(int i = 0; i < n-m+1; i+=m){
        int w = (i < n-2*m+1)? i+m : n;
        double Q_prev = Q;
        double M_prev = M;
        for(int j = i; j < w; j++){

        }

    }
    return 0;
}

/*  We don't use an explicit loop kernel. Since relative alignment is completely unknown due to T being passed in externally and compared against an arbitrary lag and subsequence length,
 *  it's too difficult to guarantee a lack of protection faults. Instead I felt it might be better to limit the amount of memory accessed. In general the compiler will turn these temporary variables
 *  into register variables and this will run over a subset of T on any given call. The general assumption is that these calls will be blocked against L2, which will accommodate most subsequence lengths.
 *  This doesn't attempt to do any kind of vectorization. I even removed it from the min operation, because this eliminates a data dependency. It won't wait on mpI if the branch is not taken. */


