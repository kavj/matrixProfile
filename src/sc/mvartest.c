#define _POSIX_C_SOURCE 200112L 
#include<stdio.h> 
#include<time.h>
#include<stdlib.h>
#include<math.h>
#include<assert.h>

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


static void winmeansig(const double* T, double* mu, double* sigma, int n, int m){
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

static void mpSelf(const double* T, const double* mu, const double* sigma, double* mp, int* mpI, const int n, const int m, const int lag){
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
    }
}


static void mvar(double* t,double* mp, long* mpI,int n, int m, int lag){
    double Mx = t[0];
    double My = t[lag];
    double Sx = 0;
    double Sy = 0;
    double Sxy = 0;

    for(int i = 1; i < m; i++){
        int k = i+lag;
        double v = t[i] - Mx;
        Mx += v/m;
        double w = t[i]-Mx;
        Sx += v*w;
        v = t[k] - My;
        Sxy += v*w;
        My += v/m;
        Sy += v*(t[k]-My);
        w = m*Sxy/(Sx*Sy);
        if(mp[i] > w){
            mp[i] = w;
            mpI[i] = k;      
        }
        if(mp[k] > w){
            mp[k] = w;
            mpI[k] = i;
        } 
    }   
}


static double controlSeq(double* T, int n, int m){
    double a = 0;
    for(int i = 0; i < n; i++){
        a += T[i]*T[i+m];
    }
    return a;
}

int main(void){
    const int n = 16*16*131072;
    const int m = n-2048;
    srand(395);
    double* x =  malloc(n*sizeof(double));
    double* mp = malloc(n*sizeof(double));
    long* mpI = malloc(n*sizeof(long)); 
    double* mu = malloc(n*sizeof(double));
    double* sigma = malloc(n*sizeof(double));
    
    for(int i = 0; i < n-m+1; i++){
        x[i] = (double)rand()/RAND_MAX;
    }
    //winmeansig(x,mu,sigma,n,m); 
    printf("end random\n");
    printf("%lf\n",mp[m-1]);
    double a = 0;
    clock_t t1 = clock();
    
    for(int i = 0; i < 100; i++){
        mvar(x,mp,mpI,n,m,1024);
       //mpSelf(x,mu,sigma,mp,mpI,n,m,m);  
        //a += controlSeq(x,n,m);
    }
    clock_t t2 = clock();
    
    printf("%lf\n",a);
    printf("time: %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    printf("%lf\n",mp[m-1]);
    free(x);
    free(mp);
    free(mpI);
    free(mu);
    free(sigma);
}




