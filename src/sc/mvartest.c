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
            sigma[j+1] = sqrt(m/Q_prev);  // we save the reciprocal instead of the standard deviation to avoid division later. It makes a big difference.
        }

    }
}


static void mvar(double* t, double* sigmaInv, double* mp, double Mx, double My, int* mpI,int n, int m, int lag){
    double Cxy = 0;
    double normFact = 1/m;
    double buffer[64];
    for(int i = 0; i < m; i++){
        Cxy += (t[i]-Mx)*(t[i+lag]-My);
    }
    double Cr = Cxy*sigmaInv[0]*sigmaInv[lag]*normFact;
    if(mp[0] < Cr){
        mp[0] = Cr;
        mpI[0] = 0;
    }
    if(mp[lag] < Cr){
        mp[lag] = Cr;
        mpI[0] = lag;
    }
    for(int i = 1; i < m; i++){
        int w = i+16 < m ? i+16: m-i;
        for(int j = i; j < w; j++){
           // buffer[j] = 
        }
        for(int j = 0; j < w; j++){
            
        }
        for(int j = 0; j < w; j++){
            /*if(mp[i] < Cr){
                mp[i] = Cr;
                //mpI[i] = k;
            }
            if(mp[k] < Cr){
                mp[k] = Cr;
                mpI[k] = i;
            }*/
        }
    }
}


static void mvar2(double* t,double* mp, int* mpI,int n, int m, int lag){
    double Mx = t[0];
    double My = t[lag];
    double Sx = 0;
    double Sy = 0;
    double Sxy = 0;
    double invm = 1/m;
    for(int i = 1; i < n-m-lag; i++){
        int k = i+lag;
        double v = t[i] - Mx;
        Mx += v*invm;
        double w = t[i]-Mx;
        Sx += v*w;
        v = t[k] - My;
        Sxy += v*w;
        My += v*invm;
        Sy += v*(t[k]-My);
        double z = 1/(Sx*Sy);
        w = m*Sxy;
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
    const int n = 131072;
    const int m = 3000;
    srand(395);
    double* x =  malloc(n*sizeof(double));
    double* mp = malloc(n*sizeof(double));
    int* mpI = malloc(n*sizeof(int)); 
    double* mu = malloc(n*sizeof(double));
    double* sigmaInv = malloc(n*sizeof(double));
    
    for(int i = 0; i < n-m+1; i++){
        x[i] = (double)rand()/RAND_MAX;
    }
    winmeansig(x,mu,sigmaInv,n,m); 
    printf("end random\n");
    printf("%lf\n",mp[m-1]);
    double a = 0;
    clock_t t1 = clock();
    
    for(int i = 0; i < n; i++){
         mvar2(x,mp,mpI,n,m,i);
      // mvar(x,mp,sigmaInv,mpI,n,m,1024);
    }

    clock_t t2 = clock();
    
    printf("%lf\n",a);
    printf("time: %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    printf("%lf\n",mp[m-1]);
    free(x);
    free(mp);
    free(mpI);
    free(mu);
    free(sigmaInv);
}
