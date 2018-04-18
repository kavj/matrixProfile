#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"

// better idea would be to make this neater then allow for MKL backend

// I will clean up the __restrict__ and other nonsense later. They're mainly used by the compiler during auto-vectorization. 
// There are a lot of different ways to set up this section in order to accommodate extremely long or short subsequence lengths and other things. 

// Queries are interleaved by some unrolling factor. This way we can stride continuously through memory when computing initial covariance with unrolling
// We would naturally use the upper triangular optimization so that it's easier to set exclusion zones. This means that for row 1, we would start with 1+excl given the symmetry of the whole thing
// otherwise it gets really nasty

// doing this in a more naive manner confused gcc in that it was seemingly worried whether things might overlap
void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, int qstart, int count, int sublen, int step){
   #pragma omp simd  // <-- check whether syntax is valid stacked, possibly inline this without restrict quantifieers
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qind = i*sublen;
      int tsind = i*step+qstart;
      double qm = mu[qind];
      for(int j = 0; j < sublen; j++){
         qbuf[qind] = ts[tsind]-qm;
         ++qind;
         ++tsind;
      }
   }
}

// prescrimp uses cross correlation, so it's only possible to fuse normalization and comparison
void max_normalize_reduction(double* __restrict__ cov, double* __restrict__ invn, int* __restrict__ blah, double qcorr, double qinvn, int qind, int offset, int count){
   #define unroll 8
   int aligned = count - count%unroll + offset;
   for(int i = offset; i < aligned; i+= unroll){
      block<double> corr;
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*qinvn;
      }
      for(int j = 0; j < unroll; j++){
         corr(j) *= invn[i+j];
      }
      block<int> cind;
      for(int j = 0; j < 4; j++){
         if(corr(2*j) < corr(2*j+1)){
            corr(2*j) = corr(2*j+1);
            cind(j) = 2*j+1;
         }
         else{
            cind(j) = 2*j;
         }
      }
      for(int j = 0; j < 2; j++){
         if(corr(4*j) < corr(4*j+2)){
            corr(4*j) = corr(4*j+2);
            cind(4*j) = cind(4*j+2);
         }
      }
      if(corr(0) < corr(4)){
         if(qcorr < corr(4)){
            qcorr = corr(4);
            qind = cind(4);
         }
      }
      else if(qcorr < corr(0)){
         qcorr = corr(0);
         qind = cind(0);
      }
   }
   *cov = qcorr;
   *blah = qind; 
   for(int i = aligned; i < offset+count; i++){
      double corr = cov[i]*qinvn*invn[i];
      if(qcorr < corr){
         qcorr = corr;
         qind =  i;
      } 
   }
   *blah = qind;  // <-- temporary so these aren't optimized out
   *cov = qcorr;
}


/*
void pxc_extrap_unifstride(double* __restrict__ qmcov, double* __restrict__ qcorr, double* __restrict__ qbuf, int64_t* __restrict__ qind, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ df, const double* __restrict__ dg, const double* __restrict__ invn, int qstart, int qcount, int step, int mxlen){
   for(int i = qstart; i < qcount; i++){
      //int jmx = std::min(mxlen,i+step);
      int jinit = qind[i];
      //double c = 
      for(int j = jinit; j < jinit+step; j++){
         
      }
   } 
}
*/



