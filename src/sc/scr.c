//#include<cstdlib>
//#include<cmath>
#include<stdlib.h>
#include<math.h>
//#include<random>
#define err 0x40
#define chunkSz 4096  // We have to keep several arrays in L2



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


/* This uses a variation on Welford's method for sample variance to compute the mean and standard deviation. 
 * It resets the summation once the term under consideration does not share any terms with the last exact summation. */
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
            M_prev = incM(M_prev,T[j+subLen],m);
            mu[j+1] = M_prev;
            sigma[j+1] = sqrt(Q_prev/m);
        }

    }
    return 0;
}


