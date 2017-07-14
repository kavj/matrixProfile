#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include "scr.h"
#define err 0x40
#define chunkSz 0x1000  // We have to keep several arrays in L2


/* allocate memory and setup structures */
/* still in debate whether I should overload this for single precision. I need to test stability first */


/* constructor for a 1D time series descriptor. */
tsdesc* sc_init(double* T, int n, int m){
    tsdesc* t = (tsdesc*)malloc(sizeof(tsdesc));
    if(t != NULL){
        t->T = T;
        t->n = n;
        t->m = m;
        t->mu = (double*)malloc(n*sizeof(double));
        t->sigmaInv = (double*)malloc(n*sizeof(double));
        if(t->mu == NULL || t->sigmaInv == NULL){
            free(t->mu);
            free(t->sigmaInv);
            free(t);
            t = NULL;
        }
    }
    return t;
}


/* somewhat generic 1D Matrix Profile Constructor*/
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
        for(int i = 0; i < n-m+1; i++){
            matp->mp[i] = (double) -1.0;
            matp->mpI[i] = 0;
        }
    }
    return matp;
}

void sc_destroy(tsdesc* t){
    free(t->mu);
    free(t->sigmaInv);
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


static double shiftMean(double mu, double a, double b, double m){
    return M + dMu/m;
}

static double incrementMean(double mu, double a, double m){
    return mu + (a - mu)/m;
}

static double shiftSSS(double s, double mu, double mup, double a, double ap){
    return s + ((a - mu) * (a - mup) - (ap - mu) * (ap - mup));
}

static double twoPassSSkern(double s, double a, double mu){
    return s + (a - mu)*(a - mu);
}

static double initSS(const double* T, double M, int m, int offset){
    double s = 0;
    for(int i = offset+1; i < offset + m; i++){
        twoPassSSkern(s,T[i],M);
    }
    return s;
}

static double initMean(const double* T, int m, int offset){
    double M = T[offset];
    for(int i = offset+1, i < offset+m; i++){
        M = incrementMean(M,T[i],i-offset+1);
    }
}


// helper function to interpolate additional points
static void muSigInterp(const double* T, double* mu, double* SS, int m, int offset){
    double M = mu[offset];
    double s = SS[offset];  
    for(int i = offset+1; i < offset + m; i++){
        double Mprev = M;
        M = shiftMean(M,T[i+m-1],T[i-1]);
        s = shiftSSS(s,M,Mprev,T[i+m-1],T[i-1]);
        mu[i]= M;
        SS[i]= s;
    }
}

static void covarInterp(){

}

/* This uses a variation on Welford's method for sample variance to compute the mean and standard deviation. 
 * It resets the summation once the term under consideration does not share any terms with the last exact summation. */
/* This should check for exceptions*/
/* I should make this clear somewhere that this isn't a full standard deviation function. It could be refactored into one,
 * but it should preserve the m factor difference given that this is used to cancel another similar factor */
void winmeansig(const double* T, double* mu, double* sigmaInv, int n, int m){
    int alignedBound = (n-m+1) - (n-m+1) % m;
    sigmaInv(0) = s;
    for(int i = 0; i < alignedBound; i += m){
        double M = initMean(T,m,i);
        double s = initSS(T,M,m,i);
        mu[i] = M;
        sigmaInv[i] = s;
        muSigInterp(T,mu,sigmaInv,m,i);
        for(int j = i; j < i+m; j++){
            sigmaInv[i] = 1/sqrt(sigmaInv[i]);
        }
    }
    /* compute unaligned portion */
    double M = initMean(T,m,alignedBound);
    double s = initSS(T,M,m,alignedBound);
    mu[alignedBound] = M;
    sigmaInv[alignedBound] = s;
    muSigInterp(T,mu,sigmaInv,n-m-alignedBound+1,alignedBound);
    for(int j = alignedBound; j < n-m+1; j++){
        sigmaInv[j] = 1/sqrt(sigmaInv[j]);
    }
}


void  sccomp(const double* T, const double* sigmaInv, double* mp, int* mpI, double Mx, double My, int k, int m, int base, int lag){
    for(int i = base; i < k; i += m){
        double corr = twoPass
    } 

}


void corrToDist(double* mp, int n, int m){
    for(int i = 0; i < n-m+1; i++){
        mp[i] = sqrt(2*m*(1-mp[i]));
    }
}

static inline void scSetupBlock(double* T, double* sigmaInv, double* mu, double* mp, int* mpI, int n, int m, int offset, int lag){

}

void scBlockSolver(tsdesc* t, matrixProfileObj* matp){
    int m = t->m;
    int n = t->n;
    printf("base n: %d\n",n);
    for(int i = 0; i < n; i+= chunkSz){       
        int lag = m;
        while(lag < n-i-m){
            scSetupBlock(t->T,t->sigmaInv,t->mu,matp->mp,matp->mpI,n,m,i,lag);
            lag++;
        }
    }
}


