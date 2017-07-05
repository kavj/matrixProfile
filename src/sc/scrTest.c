#include<stdio.h>
#include<time.h>
#include<stdlib.h>
#include<assert.h>
#include<math.h>
#include<unistd.h>
#define err 0x40
#define chunkSz 0x1000  // We have to keep several arrays in L2

#ifndef NDEBUG
#define NDEBUG
#endif

static long B;

struct mpStats{
    double* T;
    double* mu;
    double* sigma;
    int n;
    int m;
};

struct mpObj{
    double* mp;
    int* mpI;
};


typedef struct mpStats tsdesc;
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

static inline double distInc(double a, double b, double c){
    return a + b*c;
}

static inline double distDec(double a, double b, double c){
    return a - b*c;
}

static inline double znorm(double p, double mu_s, double mu_l, double sigma_s, double sigma_l){
    return (p - mu_s*mu_l)/(sigma_s*sigma_l);
}

void mpSelf(const double* T, const double* mu, const double* sigma, double* mp, int* mpI, const int n, const int m, const int lag){
    double x = 0;
    for(int i = 0; i < m; i++){
        x += T[i]*T[i+lag];
    }
    double z = znorm(x,mu[0],mu[lag],sigma[0],sigma[lag]);

    for(int i = 1; i < n-m-lag; i+=m){
        double y = x;
        x = 0;
        int w = (i+2*m+lag < n) ? i+m : n-lag-m;
        for(int j = i; j < w; j++){
            B++;
            assert(j+lag+m-1 < n);
            x = distInc(x,T[j+m-1],T[lag+j+m-1]);
            y = distDec(y,T[j-1],T[lag+j-1]);
            z = znorm(x,mu[j],mu[lag+j],sigma[j],sigma[lag+j]);
            if(z > mp[j]){
                mp[j] = z;
                mpI[j] = lag+j;
            }
            if(z > mp[lag+j]){
                mp[lag+j] = z;
                mpI[lag+j] = j;
            }
        }
       // printf("n: %d m: %d lag: %d  w: %d\n",n,m,lag,w);
    }  
}

static void testMPAlg(double* T, double* mp, int* mpI, double* mu, double* sig, int n, int m, int lag){
   double mux = mu[0];
   double muy = mu[lag];
   double cov = 0;
   double invm = 1/m;
   for(int i = 0; i < m; i++){
      /* double dx = T[i] - mux;
       double dy = T[i+lag] - muy;
       mux += mux*invm;
       muy += muy*invm;*/
       cov += (T[i]-mux)*(T[i+lag]-muy);
   }
   for(int i = 0; i < n - m - lag; i++){
       
       double dx = T[i] - mux;
       double dy = T[i+lag] - muy;
       mux -= dx*invm; //mux -= dx/m;
       muy -= dx*(T[lag+i]-muy)*invm; // /m;
       cov -= dx*(T[lag+i]-muy);
       dx = T[i+m] - mux;
       mux += dx*invm;  // /m;
       muy += dx*(T[lag+i+m]-muy)*invm; // /m;
       cov += dx*(T[lag+i+m]-muy);
       double r = cov*(sig[i+1]*sig[lag+i+1]);
       if(r > mp[i]){
           mp[i] = r;
           mpI[i] = i;
       }
       if(r > mp[i+lag]){
           mp[i+lag] = r;
           mpI[i+lag] = r;
       }
   }
}


static inline void scSetupBlock(double* T, double* mu,  double* sigma, double* mp, int* mpI, int n, int m, int lag, int offset){
   int k = chunkSz < n - offset ? chunkSz : n - offset;
   testMPAlg(&T[offset],&mp[offset],&mpI[offset],&mu[offset],&sigma[offset],k,m,lag);
   //mpSelf(&T[offset],&mu[offset],&sigma[offset],&mp[offset],&mpI[offset],k,m,lag);
}

void scBlockSolver(tsdesc* t, matrixProfileObj* matp){
    int m = t->m;
    int n = t->n;
    printf("base n: %d\n",n);
    assert(t != NULL);
    assert(matp != NULL);
    for(int i = 0; i < n; i+= chunkSz){    
    //   printf("i: %d\n",i);
       // printf("%d start\n",i);
        int lag = 0;
        while(lag < n-i-m){
            //for(int lag = m; lag < n-i-m; lag++){ // technically could be lag++, we iterate over the same portion with different lag
            scSetupBlock(t->T,t->mu,t->sigma,matp->mp,matp->mpI,n,m,lag,i);
            lag++;
        } 
        printf("lag: %d\n",lag);
    }
} 




int main(void){
    B = 0;
    const int n = 131072*4;
    const int m = 400;
    srand(395);
    double* x = (double*)malloc(n*sizeof(double));
    double* b = (double*)malloc(n*sizeof(double));
    matrixProfileObj* matp = mp_init(n,m);
    if((x == NULL) || (matp == NULL)){
        printf("malloc failed me again\n");
        free(x);
        free(matp);
        exit(1);
    }    
    for(int i = 0; i < n-m+1; i++){
        x[i] = (double)rand()/RAND_MAX;
    }
    printf("end random\n");
    tsdesc* t = sc_init(x,n,m);
    winmeansig(x,t->mu,t->sigma,n,m);
    clock_t t1 = clock();
    scBlockSolver(t,matp);
        // testMPAlg(x,matp->mp,matp->mpI,t->mu,t->sigma,n,m,lag);
    clock_t t2 = clock();
    printf("time: %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    mp_destroy(matp);
    free(x);
    free(b);
    printf("%lu\n",B);
}
