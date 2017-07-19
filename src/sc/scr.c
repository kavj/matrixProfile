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
    return mu + (a-b)/m;
}

static double incrementMean(double mu, double a, double m){
    return mu + (a - mu)/m;
}

static double shiftSSS(double s, double mu, double muprev, double x, double xprev){
    return s + ((x - mu) * (x - muprev) - (xprev - mu) * (xprev - muprev));
}

static double shiftSXY(double x, double y, double xprev, double yprev, double mux, double muyprev){
    return (x - mux)*(y - muyprev) - (xprev - mux)*(yprev - muyprev);
}

static double initSS(const double* T, double M, int m, int offset){
    double s = 0;
    for(int i = offset; i < offset + m; i++){
        s += (T[i] - M) * (T[i] - M);
    }
    return s;
}

static double centeredSS(const double* X, const double* Y, double Mx, double My, int winLen){
    double sxy = 0;
    for(int i = 0; i < winLen; i++){
        sxy += (X[i] - Mx)*(Y[i] - My);
    }
    return sxy;
}


static double initMean(const double* T, int m, int offset){
    double M = T[offset];
    for(int i = offset+1; i < offset+m; i++){
        M = incrementMean(M,T[i],i-offset+1);
    }
    return M;
}

/* This uses a variation on Welford's method for sample variance to compute the mean and standard deviation. 
 * It resets the summation once the term under consideration does not share any terms with the last exact summation. */
/* This should check for exceptions*/
/* I should make this clear somewhere that this isn't a full standard deviation function. It could be refactored into one,
 * but it should preserve the m factor difference given that this is used to cancel another similar factor */
void winmeansig(const double* T, double* mu, double* sigma, int n, int m,int normConstant){
    int alignedBound = (n-m+1) - (n-m+1) % m;
    //printf("n: %d m: %d aligned: %d",n,m,alignedBound);
    for(int i = 0; i < alignedBound; i += m){
        double M    = initMean(T,m,i);
        double s    = initSS(T,M,m,i);
        //printf("%lf\n",M);
        mu[i]       = M;
        sigma[i]    = sqrt(s/normConstant);   

        for(int j = i+1; j < i+m; j++){
            double Mprev = M;
            M        = shiftMean(M,T[j+m-1],T[j-1],m);
            s        = shiftSSS(s,M,Mprev,T[j+m-1],T[j-1]);
            mu[j]    = M;
            sigma[j] = sqrt(s/normConstant);
        }

    }

    /* compute unaligned portion */
    double M               = initMean(T,m,alignedBound);
    double s               = initSS(T,M,m,alignedBound);
    mu[alignedBound]       = M;
    sigma[alignedBound]    = s;
    for(int i = alignedBound; i < n-m+1; i++){
        sigma[i] = sqrt(s/normConstant);
    }
}


static void corrToDist(double* mp, int n, int m){
    for(int i = 0; i < n-m+1; i++){
        mp[i] = sqrt(2*m*(1-mp[i]));
    }
}

static void computeBlock(const double* T, const double* mu, const double* sigmaInv, double* mp, int* mpI, int n, int m, int base, int offset){
    double Mx = mu[base];
    double My = mu[offset];
    double covxy = centeredSS(&T[base],&T[offset],Mx,My,m);
    double corrxy = covxy*sigmaInv[base]*sigmaInv[offset];
    int upperBound = (chunkSz - m + 1) < (n - m - offset + 1) ? (chunkSz - m + 1) : (n - m - offset + 1);
    if(mp[base] < corrxy){
        mp[base]  = corrxy;
        mpI[base] = offset;
    }
    if(mp[offset] < corrxy){
        mp[offset]  = corrxy;
        mpI[offset] = base;
    }
    for(int i = 1; i < n - offset - m + 1; i++){        
        double Myprev = My;
        Mx = mu[base+i];
        My = mu[offset+i];
        covxy = shiftSXY(T[base+i+m-1],T[offset+i+m-1],T[base+i-1],T[offset+i-1],Mx,Myprev);
        corrxy = covxy*sigmaInv[base+i]*sigmaInv[offset+i];
        if(mp[base+i] < corrxy){
            mp[base+i]  = corrxy;
            mpI[base+i] = offset + i;
        }
        if(mp[offset+i] < corrxy){
            mp[offset+i]  = corrxy;
            mpI[offset+i] = base + i;
        }   
    }
}


static void blockupdateED(tsdesc* t, matrixProfileObj* matp, int blockSt, int blockLen){
    int minlag = t->m;
    int count  = t->n - t->m + 1;
    double* T = t->T;
    double* sigmaInv = t->sigmaInv;
    for(int offs = blockSt + minlag; offs < count-t->m; offs++){ // t->m here is debug code
       computeBlock(T,t->mu,sigmaInv,matp->mp,matp->mpI,t->n,t->m,blockSt,chunkSz);
    }
}


void sjbs(tsdesc* t, matrixProfileObj* matp){

    int minlag  = t->m;
    int count   = t->n - minlag - t->m + 1;
    int step    = chunkSz - t->m + 1;
    int aligned = t->n - t->n%chunkSz;
 
    for(int base = 0; base < aligned; base += step){
        blockupdateED(t,matp,base,chunkSz);
    }
    if(aligned < count){
        blockupdateED(t,matp,aligned,t->n-aligned+1);
    }
    corrToDist(matp->mp,t->n,t->m);
}

/* This is the anytime version */
void mpIndexBlkSolver(tsdesc* t, matrixProfileObj* matp, int* permI, int strt, int perc){
// get the number of indices required.
// sort indices
// call block solve for each index
// finish writing this at some point 
}
