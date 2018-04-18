#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"

// better idea would be to make this neater then allow for MKL backend

// I will clean up the __restrict__ and other nonsense later. They're mainly used by the compiler during auto-vectorization. 
// There are a lot of different ways to set up this section in order to accommodate extremely long or short subsequence lengths and other things. 

// Queries are interleaved by some unrolling factor. This way we can stride continuously through memory when computing initial covariance with unrolling
// We would naturally use the upper triangular optimization so that it's easier to set exclusion zones. This means that for row 1, we would start with 1+excl given the symmetry of the whole thing
// otherwise it gets really nasty

/*
void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, int qstart, int unroll, int sublen, int step){
   // hoist means 
   #pragma omp parallel for
   for(int i = 0; i < unroll; i++){
      double m = mu[i*step+qstart];
      for(int j = 0; j < sublen; j++){
         qbuf[i*step+j] = ts[i*step+j+qstart] - m;
      }
   }
}
*/

// doing this in a more naive manner confused gcc in that it was seemingly worried whether things might overlap
void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, int qstart, int count, int sublen, int step){
   //block<double> qm;
   // hoist means 
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

/*
// it would be best to do this using an upper triangular pattern. It avoids a lot of nonsense
static inline void partial_xcorr_kern(double* __restrict__ qmcov, double* __restrict__ corr, double* __restrict__ qbuf, int64_t* __restrict__ ind, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int start, int unroll, int sublen, int step, int excl){
   block<double> cov;
   int span = 0; // dummy var while I figure out what I want to do post changes
   for(int i = start; i < span; i++){ // <-- needs an arbitrary start point
      double m = mu[i];
      for(int j = 0; j < sublen; j++){  
         double c = ts[i] - m;
         for(int k = 0; k < unroll; k++){
            cov(k) += qbuf[j*unroll+k]*c;
         } 
      }
      m = invn[i];
      for(int j = 0; j < unroll; j++){
         cov(j+unroll) = cov(j)*m;
      }
      for(int j = 0; j < unroll; j++){
         cov(j+unroll) *= invn[start+j*step];
      }
      double qopt = -1.0;
      int qopt_i = -1;
      for(int k = 0; k < unroll; k++){
         int qi = start+k*step;
         if(cov(k+unroll) > corr[qi]){
            corr[qi] = cov(k+unroll);
            qmcov[qi]= cov(k);
         }
         if(qopt < cov(k+unroll)){ // <-- hello stall cycles
            qopt = cov(k+unroll);
            qopt_i = qi;
         }
      }
      if(qopt > corr[i]){
         corr[i] = qopt;
         ind[i] = qopt_i;
      }
   }
}*/

/*
// This takes a query range and preallocated buffers
// qcorr could technically be replaced with temporary variables. It's just rescaled from qcov
// In general we only care about the maxima here, which helps cut down on excess bus traffic
void prescr_partial_xcorr(double* __restrict__ qmcov, double* __restrict__ qcorr, double* __restrict__ qbuf, int64_t* __restrict__ qind, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int qstart, int qcount, int len, int sublen, int step, int excl){
   int rem = qcount;
   int maxqs = 8; // <-- Should be tuned using subsequence length and either detected cache size or dummy estimate

   while(rem > 0){
      int unroll = std::min(8,rem); // <-- this needs 2 steps of unrolling if we can reasonably stack enough queries to justify parallelizing this part
      int start = qstart+(qcount-rem)*step;
      for(int i = 0; i < unroll; i++){
         batch_normalize(qbuf,ts,mu,start+i*unroll,unroll,sublen,step);  
      }
      #pragma omp parallel for
      for(int i = 0; i < unroll; i++){  // <-- dummy setup
         partial_xcorr_kern(qmcov,qcorr,qbuf,qind,ts,mu,invn,start,unroll,sublen,step,excl); 
      }
      rem -= unroll;
   }
}*/


// This is a slightly exotic reduction using hoisted blocks of 8. I suspect it may be hard to beat even with vectorization due to the reduced throughput of blend operations and much higher register pressure
// This kind of thing often confuses gcc's scheduler

void max_reduction_8x1(double* __restrict__ cov, double* __restrict__ invn, int* __restrict__ blah, double qcorr, double qinvn, int qind, int offset, int count){
   //const int def_count = 8; // <-- temp val
   #define unroll 8
   for(int i = offset; i < offset+count; i+= unroll){
      block<double> corr;
      /*for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*qinvn;
      }*/
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
   *blah = qind; // <-- temporaries
}

//void max_reduction(double* __restrict__
// if count > 8 call 8x1 version  possibly inlined without restrict keywords then include logic for the reference here?
// should use peeling here

void max_reduction(double* __restrict__ cov, double* __restrict__ invn, int* __restrict__ blah, double qcorr, double qinvn, int qind, int offset, int count){
   block<double> corr;
   //int unaligned = 
   for(int i = offset; i < offset+count; i++){
      
   }

}

void pxc_extrap_unifstride(double* __restrict__ qmcov, double* __restrict__ qcorr, double* __restrict__ qbuf, int64_t* __restrict__ qind, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ df, const double* __restrict__ dg, const double* __restrict__ invn, int qstart, int qcount, int step, int mxlen){
   for(int i = qstart; i < qcount; i++){
      //int jmx = std::min(mxlen,i+step);
      int jinit = qind[i];
      //double c = 
      for(int j = jinit; j < jinit+step; j++){
         
      }
   } 
}


