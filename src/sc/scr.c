#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define err 0x40
#define chunkSz 4096  // We have to keep several arrays in L2


static long B;

struct desc{
    double* T;
    double* mu;
    double* sigmaInv;
    int n;
    int m;
    int chsz;
};

struct mpObj{
    double* mp;
    int* mpI;
};


typedef struct desc tsdesc;
typedef struct mpObj matrixProfileObj;





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
    return S - (x-M)*(x-M)/k;
}


/* M here indicates the mean of k-1 values, not k values. You can just decrement M first*/
static double incS(double S, double M, double x, int k){
    return S + (x-M)*(x-M)/k;
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
    sigmaInv[0] = sqrt(m/Q);
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
            sigmaInv[j+1] = sqrt(m/Q_prev);
        }

    }
}


/* These "effectively" add and remove terms from mean and sume of products calculations. */
/* C means data is already centered.*/
static double decCM(double M, double x, double invkk, int k){
    return (k*M-x)*invkk;
}

static double incCM(double M, double x, double invk){
    return M + x*invk;
}

static double decCSxy(double Sxy, double x, double y, double kkinvk){
    return Sxy - x*y*kkinvk;
}


/* M here indicates the mean of k-1 values, not k values. You can just decrement M first*/
static double incCSxy(double Sxy, double x, double y, double kkinvk){
   return Sxy + x*y*kkinvk;
}



/* The base/offset part is annoying, but its used to correctly set the index*/
/* At a high level, this is a somewhat optimized sum of squares formula, used to estimate the normalized difference between subsequences as a function of their correlation*/
/* I should probably precompute these constants and just pass a struct, depending on whether the compiler behaves*/
void  sccomp(const double* T, const double* sigmaInv, double* mp, int* mpI, double Mx, double My, int n, int m, int base, int lag){
    double Sxy = 0;
    double Cxy = sigmaInv[0]*sigmaInv[lag];
    double invm = (double)1/m;
    double invmm = (double)(1/(m-1));
    double mminv = (double)(m-1)/m;

    for(int i = 0; i < m; i++){ 
        Sxy += (T[i]-Mx)*(T[i+lag]-My);
    }
    Cxy *=  Sxy;
    if(mp[0] < Cxy){  /* Branching tends to be faster than unnecessary writes here */
        mp[0] = Cxy;
        mpI[0] = base+lag;
    }
    if(mp[base+lag] < Cxy){
        mp[lag] = Cxy;
        mpI[lag] = base;
    }
    for(int i = 0; i < n-m-lag; i++){
        Cxy = sigmaInv[i+1]*sigmaInv[i+lag+1];

        Mx = decCM(Mx,T[i],invmm,m);
        My = decCM(My,T[i+lag],invmm,m);

        double x = T[i+m] - Mx;
        double y = T[i+m+lag] - My;

        Sxy = incCSxy(Sxy,x,y,invmm);
        Cxy *= Sxy;

        Mx = incCM(Mx,x,invm);
        My = incCM(My,y,invm);

        if(mp[i+1] < Cxy){
            mp[i+1] = Cxy;
            mpI[i+1] = i+base+lag;
        }
        if(mp[i+lag+1] < Cxy){
            mp[i+lag+1] = Cxy;
            mpI[i+lag+1] = i+base;
        }
    }
}


void corrToDist(double* mp, int n, int m){
    for(int i = 0; i < n-m+1; i++){
        mp[i] = sqrt(2*m*(1-mp[i]));
    }
}

static inline void scSetupBlock(double* T, double* mu,  double* sigmaInv, double* mp, int* mpI, int n, int m, int lag, int offset){
   int k = chunkSz < n - offset ? chunkSz : n - offset;
   sccomp(&T[offset],&sigmaInv[offset],&mp[offset],&mpI[offset],mu[offset],mu[offset+lag],k,m,offset,lag);
}

void scBlockSolver(tsdesc* t, matrixProfileObj* matp){
    int m = t->m;
    int n = t->n;
    printf("base n: %d\n",n);
    for(int i = 0; i < n; i+= chunkSz){       
        int lag = 0;
        while(lag < n-i-m){
            scSetupBlock(t->T,t->mu,t->sigmaInv,matp->mp,matp->mpI,n,m,lag,i);
            lag++;
        }
        printf("lag: %d\n",lag);
    }
}


