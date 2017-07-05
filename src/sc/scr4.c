#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<unistd.h>
#define err 0x40
#define chunkSz 0x10000  // We have to keep several arrays in L2

#ifndef NDEBUG
#define NDEBUG
#endif


/* allocate memory and setup structures */
/* still in debate whether I should overload this for single precision. I need to test stability first */
tsdesc* sc_init(const double* T, int n, int m){
    tsdesc* t = (tsdesc*)malloc(sizeof(tsdesc));
    if(t != NULL){
        t->T = T;
        t->n = n;
        t->m = m;
        t->mu = (double*)malloc(n*sizeof(double));
        t->sigma = (double*)malloc(n*sizeof(double));
        if(t->mu == NULL || t->sigma == NULL){
            free(t->mu);
            free(t->sigma);
            free(t);
            t = NULL;
        }
    }
    return t;
}

matrixProfileObj* mp_init(int n, int m){
    matrixProfileObj* matp = malloc(sizeof(matrixProfileObj));
    if(matp != NULL){
        matp->mp = malloc((n-m+1)*sizeof(double));
        matp->mpI = malloc((n-m+1)*sizeof(double));
        if(matp->mp == NULL || matp->mpI == NULL){
            free(matp->mp);
            free(matp->mpI);
            free(matp);
            matp = NULL;
        }
    }
    return matp;
}

void sc_destroy(tsdesc* t){
    free(t->mu);
    free(t->sigma);
    free(t);
}

void mp_destroy(matrixProfileObj* matp){
    free(matp->mp);
    free(matp->mpI);
    free(matp);
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


/* This uses a variation on Welford's method for sample variance to compute the mean and standard deviation. 
 * It resets the summation once the term under consideration does not share any terms with the last exact summation. */
void winmeansig(const double* T, double* mu, double* sigma, int n, int m){
    double M = T[0];
    double Q = 0;
    for(int i = 1; i < m; i++){
        Q = incQ(Q,M,T[i],i+1);
        M = incM(M,T[i],i+1);
    }
    mu[0] = M;
    sigma[0] = sqrt(Q/m);
    for(int i = 0; i < n-m+1; i+=m){
        int w = (i < n-2*m+1) ? i+m : n-m;
        double Q_prev = Q;
        double M_prev = M;
        Q = 0;
        M = T[i+1];
        for(int j = i; j < w; j++){
            if(j != i){
                Q = incQ(Q,M,T[j+m],j-i+1);
                M = incM(M,T[j+m],j-i+1);
            }
            M_prev = decM(M_prev,T[j],m);
            Q_prev = decQ(Q_prev,M_prev,T[j],m);
            Q_prev = incQ(Q_prev,M_prev,T[j+m],m);
            M_prev = incM(M_prev,T[j+m],m);
            mu[j+1] = M_prev;
            sigma[j+1] = sqrt(Q_prev/m);
        }

    }
}


static inline void scSetupBlock(const double* T, const double* mu, const double* sigma, double* mp, int* mpI, int n, int m, int lag, int offset){
}

void scBlockSolver(tsdesc* t, matrixProfileObj* matp){
    int m = t->m;
    int n = t->n;
    //printf("base n: %d\n",n);
    assert(t != NULL);
    assert(matp != NULL);
    for(int i = 0; i < n; i+= chunkSz){    
       // printf("i: %d\n",i);
        for(int lag = m; lag < n-i-m; lag++){ // technically could be lag++, we iterate over the same portion with different lag
            scSetupBlock(t->T,t->mu,t->sigma,matp->mp,matp->mpI,n,m,lag,i);
        } 
    }
} 

