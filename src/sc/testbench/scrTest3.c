#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include "scr.h"
#define err 0x40
#define chunkSz 0x1000  // We have to keep several arrays in L2


/* allocate memory and setup structures */
/* still in debate whether I should overload this for single precision. I need to test stability first */

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











void corrToDist(double* mp, int n, int m){
    for(int i = 0; i < n-m+1; i++){
        mp[i] = sqrt(2*m*(1-mp[i]));
    }
}

static inline void scSetupBlock(double* T, double* sigmaInv, double* mu, double* mp, int* mpI, int n, int m, int offset, int lag){
   int k = chunkSz < n - offset ? chunkSz : n - offset;
   sccomp(T,sigmaInv,mp,mpI,k,m,offset,lag);
   //sccomp(&T[offset],&sigmaInv[offset],&mp[offset],&mpI[offset],mu[offset],mu[offset+lag],k,m,offset,lag);
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
       // printf("lag: %d\n",lag);
    }
}


