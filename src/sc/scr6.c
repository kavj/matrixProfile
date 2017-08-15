#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<unistd.h>
#include<time.h>
#include "scr.h"
#define err 0x40
#define kernWid 4
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

static inline double initMean(const double* T, int m, int offset){
    double M = T[offset];
    for(int i = offset+1; i < offset+m; i++){
        M += (T[i] - M)/(i-offset+1);
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

static inline double sstpupdate(double x, double y, double Mx, double My){
    return (x - Mx)*(y - My);
}

static inline void updateMP(double* mp, int* mpI, double corr, int base, int offset){
    if(mp[base] < corr){
        mp[base] = corr;
        mpI[base] = offset;
    }
    if(mp[offset] < corr){
        mp[offset] = corr;
        mpI[offset] = base; 
    }
}

void computeBlock(const double* T, const double* mu, const double* sigmaInv, double* mp, int* mpI, int base, int chunkLen, int subLen, int offset){
    //#define kernWid 0x04
    double Mx;
    double My0 = mu[offset];
    double My1 = mu[offset+1];
    double My2 = mu[offset+2];
    double My3 = mu[offset+3];  
    double c0  = 0;
    double c1  = 0;
    double c2  = 0;
    double c3  = 0;   
 
    for(int i = 0; i < subLen; i++){
        c0 += (T[base+i]-Mx)*(T[offset+i]-My0);
        c1 += (T[base+i]-Mx)*(T[offset+i]-My1);
        c2 += (T[base+i]-Mx)*(T[offset+i]-My2);
        c3 += (T[base+i]-Mx)*(T[offset+i]-My3);
    }
    // rename variables to be less misleading with the current design
    c0 *= sigmaInv[base]*sigmaInv[offset];
    c1 *= sigmaInv[base]*sigmaInv[offset+1];
    c2 *= sigmaInv[base]*sigmaInv[offset+2];
    c3 *= sigmaInv[base]*sigmaInv[offset+3];
    // it's worth noting that this part could be improved using horizontal comparison for the "base" operand
    updateMP(mp,mpI,c0,base,offset);
    updateMP(mp,mpI,c1,base,offset+1);
    updateMP(mp,mpI,c2,base,offset+2);
    updateMP(mp,mpI,c3,base,offset+3);
    
    for(int i = 1; i < chunkLen; i++){
        c0 += shiftSXY(T[base+i+subLen],T[offset+i+subLen],T[base+i-1],T[offset+i-1],mu[base+i],mu[offset+i-1]);
        c1 += shiftSXY(T[base+i+subLen],T[offset+i+subLen],T[base+i-1],T[offset+i],mu[base+i],mu[offset+i]);
        c2 += shiftSXY(T[base+i+subLen],T[offset+i+subLen],T[base+i-1],T[offset+i+1],mu[base+i],mu[offset+i+1]);
        c3 += shiftSXY(T[base+i+subLen],T[offset+i+subLen],T[base+i-1],T[offset+i+2],mu[base+i],mu[offset+i+2]);
        updateMP(mp,mpI,c0,base+i,offset+i);
        updateMP(mp,mpI,c1,base+i,offset+i+1);
        updateMP(mp,mpI,c2,base+i,offset+i+2);
        updateMP(mp,mpI,c3,base+i,offset+i+3);
    }
}

/* This part could use some cleanup. Unpacking structs is a bit messy here */
void blockupdateED(double* T, double* mp, double* mu, double* sigmaInv, int* mpI, int base, int offs, int seqLen, int chunkLen, int subLen){
    int step = chunkLen - subLen + 1; 
    //int stagger = 4; // sentinel
    while(base + chunkLen - 1 < seqLen){
        for(int offset = base + subLen; offset < seqLen-subLen+1; offset += kernWid){
            computeBlock(T,mu,sigmaInv,mp,mpI,base,chunkLen,subLen,offset);
        }
        base += step;
    }
}

/* This is the anytime version */
void mpIndexBlkSolver(tsdesc* t, matrixProfileObj* matp, int* permI, int strt, int perc){
}
