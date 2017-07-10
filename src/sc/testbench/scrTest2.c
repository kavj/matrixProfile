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

static double decM(double M, double x, int k){
    return (k*M-x)/(k-1);
}

static double incM(double M, double x, int k){
    return M + (x-M)/k;
}

static double decS(double S, double M, double x, int k){
    return S - (k-1)*(x-M)*(x-M)/k;
}


/* M here indicates the mean of k-1 values, not k values. You can just decrement M first*/
static double incS(double S, double M, double x, int k){
    return S + (k-1)*(x-M)*(x-M)/k;
}


/* This uses a variation on Welford's method for sample variance to compute the mean and standard deviation. 
 * It resets the summation once the term under consideration does not share any terms with the last exact summation. */
/* This should check for exceptions*/
void winmeansig(const double* T, double* mu, double* sigmaInv, int n, int m){
    double M = T[0];
    double Q = 0;
    for(int i = 1; i < m; i++){
        Q = incS(Q,M,T[i],i+1);
        M = incM(M,T[i],i+1);
    }
    mu[0] = M;
    sigmaInv[0] = sqrt(Q/m);
    for(int i = 0; i < n-m; i+=m){
        int w = (i < n-2*m+1)? i+m : n-m;
        double Q_prev = Q;
        double M_prev = M;
        Q = 0;
        M = T[i+m];
        for(int j = i; j < w; j++){
            if(j != i){
                Q = incS(Q,M,T[j+m],j-i+1);
                M = incM(M,T[j+m],j-i+1);
            }
            M_prev = decM(M_prev,T[j],m);
            Q_prev = decS(Q_prev,M_prev,T[j],m);
            Q_prev = incS(Q_prev,M_prev,T[j+m],m);
            M_prev = incM(M_prev,T[j+m],m);
            mu[j+1] = M_prev;
            sigmaInv[j+1] = sqrt(Q_prev/m);
        }

    }
}



/* These "effectively" add and remove terms from mean and sume of products calculations. */
/* C means data is already centered.*/

void sccomp(double* T, double* mu, double* sigmaInv, double* mp, int* mpI, int n, int m, int offset, int lag){
   double Mx = mu[offset]; 
   double My = mu[offset+lag];
   double Sx = 0;
   double Sy = 0;
   double C  = 0;
   for(int i = offset; i < offset+m; i++){
       double x = T[i] - Mx;
       double y = T[i+lag] - My;
       //Sx += x*x;
       //Sy += y*y;
       C  += x*y;
   }
   double c = C/sqrt(Sy*Sx);
   if(c < mp[offset]){
       mp[offset] = c;
       mpI[offset] = offset+lag;
   }
   if(c < mp[offset+lag]){
       mp[offset+lag] = c;
       mpI[offset+lag] = offset;
   }
   for(int i = offset+1; i < n-m; i++){
       double Mxp = Mx;  // Mx previous
       double Myp = My;
       //double S = sigmaInv[i]*sigmaInv[i+lag];
       Mx = mu[i];
       My = mu[i+lag];
       double xp = T[i-1];
       double x  = T[i-1+m];   
       double yp = T[i-1+offset];
       double y  = T[i-1+offset+m];
       //Sx += (x - Mx)*(x - Mxp) - (xp - Mx)*(xp - Mxp);
       //Sy += (y - My)*(y - Myp) - (yp - My)*(yp - Myp);
       double S = sqrt(Sx*Sy);
       
       C  += (x - Mxp)*(y - My) - (xp - Mxp)*(yp - My);
       S = C/S;
       if(S < mp[i]){
          mp[i] = S;
          mpI[i] = i+lag;
       } 
       if(S < mp[i+lag]){
          mp[i+lag] = S;
          mpI[i+lag] = i;
       }
   }

}


static double incCS(double S, double xp, double y, double m){
   return S + xp*y/m;
}



void corrToDist(double* mp, int n, int m){
    for(int i = 0; i < n-m+1; i++){
        mp[i] = sqrt(2*m*(1-mp[i]));
    }
}

static inline void scSetupBlock(double* T, double* sigmaInv, double* mu, double* mp, int* mpI, int n, int m, int offset, int lag){
   int ch = chunkSz >= 2*m ? chunkSz : 2*m;

   int k = ch < n - offset ? ch : n - offset;
   sccomp(T,mu,sigmaInv,mp,mpI,k,m,offset,lag);
}

void scBlockSolver(tsdesc* t, matrixProfileObj* matp){
    int m = t->m;
    int n = t->n;
    printf("base n: %d\n",n);
    int ch = chunkSz >= 2*m ? chunkSz : 2*m;
    for(int i = 0; i < n-m; i+= ch){       
        int lag = m;
        while(lag < n-i-m-ch){
            scSetupBlock(t->T,t->sigmaInv,t->mu,matp->mp,matp->mpI,n,m,i,lag);
            lag++;
        }
       // printf("lag: %d\n",lag);
    }
}


