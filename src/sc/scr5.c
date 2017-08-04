#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<time.h>
#include "scr.h"
#define err 0x40
//#define chunkSz 0x1F40  // We have to keep several arrays in L2


/* allocate memory and setup structures */
/* still in debate whether I should overload this for single precision. I need to test stability first */




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


static inline double shiftMean(double mu, double a, double b, double m){
    return mu + (a-b)/m;
}

static inline double incrementMean(double mu, double a, double m){
    return mu + (a - mu)/m;
}

static inline double shiftSSS(double s, double mu, double muprev, double x, double xprev){
    return s + ((x - mu) * (x - muprev) - (xprev - mu) * (xprev - muprev));
}

static inline double shiftSXY(double x, double y, double xprev, double yprev, double mux, double muyprev){
    return (x - mux)*(y - muyprev) - (xprev - mux)*(yprev - muyprev);
}

static inline double initSS(const double* T, double M, int m, int offset){
    double s = 0;
    for(int i = offset; i < offset + m; i++){
        s += (T[i] - M) * (T[i] - M);
    }
    return s;
}

static inline double centeredSS(const double* X, const double* Y, double Mx, double My, int winLen){
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
    for(int i = 0; i < alignedBound; i += m){
        double M    = initMean(T,m,i);
        double s    = initSS(T,M,m,i);
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


static inline void computeBlock(const double* T, const double* mu, const double* sigmaInv, double* mp, int* mpI, int baseSt, int chunkLen, int subLen, int offset){
    double Mx = mu[base];
    double My = mu[offset];
    double corrxy = sigmaInv[base]*sigmaInv[offset];
    double covxy = centeredSS(&T[base],&T[offset],Mx,My,m);
    corrxy *= covxy;
    if(mp[base] < corrxy){
        mp[base]  = corrxy;
        mpI[base] = offset;
    }
    if(mp[offset] < corrxy){
        mp[offset]  = corrxy;
        mpI[offset] = base;
    }
    for(int i = 1; i < chunkLen; i++){
        double Myprev = My;
        Mx = mu[base+i];
        My = mu[offset+i];
        double corrxy = sigmaInv[base+i]*sigmaInv[offset+i];
        covxy = shiftSXY(T[base+i+m-1],T[offset+i+m-1],T[base+i-1],T[offset+i-1],Mx,Myprev);
        corrxy *= covxy;
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





/* This part could use some cleanup. Unpacking structs is a bit messy here */
static inline void blockupdateED(double* T, double* mp, double* mu, double* sigmaInv, int* mpI, int base, int offs, int chunkLen, int subLen){
    int step = chunkLen - m + 1; 
    int stagger = 4; // sentinel
    while(base + chunkLen - 1 < n){
        for(int offset = base + m; offset < n-m+1; offset += stagger){
            computeBlock(T,mu,sigmaInv,mp,mpI,base,chunkLen,subLen,offset);
        }
        base += step;
    }
}


/*void sjbs(double* T, double* Mu, double* sigmaInv, int n, int m){
    int chunkLen = 2048; // sentinel for now
    int count = n - 2*m + 1;
    int step  = 3*m + 1;
    int stagg = 4;       // unroll factor, used to hide latency, if each is computed independently, then it's 1
    int aligned = n - n%(4*>m);
    for(base = 0; base < n - ; base += step){
        for(int offset = base + m; offset < n - chunkLen + 1; offset += stagg){
            
        }
    } */
    /*for(int base = 0; base < aligned; base += step){
       blockupdateED(t,matp,base);
    }*/
    /*if(aligned < count){
        blockupdateED(t,matp,aligned,t->n-aligned+1);
    }*/
/*    corrToDist(mp,n,m);
}*/

/* This is the anytime version */
void mpIndexBlkSolver(tsdesc* t, matrixProfileObj* matp, int* permI, int strt, int perc){
// get the number of indices required.
// sort indices
// call block solve for each index
// finish writing this at some point 
}
