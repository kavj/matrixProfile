#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"

// I will clean up the __Restrict__ and other nonsense later. They're mainly used by the compiler during auto-vectorization. 
// There are a lot of different ways to set up this section in order to accommodate extremely long or short subsequence lengths and other things. 

// Queries are interleaved by some unrolling factor. This way we can stride continuously through memory when computing initial covariance with unrolling
// We would naturally use the upper triangular optimization so that it's easier to set exclusion zones. This means that for row 1, we would start with 1+excl given the symmetry of the whole thing
// otherwise it gets really nasty

void batch_normalize(double* __Restrict__ qbuf,  const double* __Restrict__ ts, const double* __Restrict__ mu, int qstart, int qcount, int sublen, int step){
   int back = len-sublen+1;
   block<double> qm;
   // hoist means 
   for(int i = 0; i < qcount; i++){
      qm(i) = mu[i*step+qstart];
   } 
   for(int i = 0; i < sublen; i++){
      for(int j = 0; j < qcount; j++){
         q[i*qcount+j] = ts[i+qstart+j*step]-qm(j);
      }
   }
}



// it would be best to do this using an upper triangular pattern. It avoids a lot of nonsense

static inline void partial_xcorr_kern(double* __Restrict__ qmcov, double* __Restrict__ corr, double* __Restrict__ qbuf, int64_t* __Restrict__ ind, const double* __Restrict__ ts, const double* __Restrict__ mu, const double* __Restrict__ invn, int qstart, int unroll, int sublen, int qstart, int step, int excl){
   block<double> cov;
   for(int i = 0; i < len; i++){ // <-- needs an arbitrary start point
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
         cov(j+unroll) *= invn[qstart+j*step];
      }
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
}


// This takes a query range and preallocated buffers
// qcorr could technically be replaced with temporary variables. It's just rescaled from qcov
// In general we only care about the maxima here, which helps cut down on excess bus traffic
void prescr_partial_xcorr(double* __Restrict__ qmcov, double* __Restrict__ qcorr, double* __Restrict__ qbuf, int64_t* __Restrict__ qind, const double* __Restrict__ ts, const double* __Restrict__ mu, const double* __Restrict__ invn, int qstart, int qcount, int len, int sublen, int step, int excl){
   int rem = qcount;
   int maxqs = 8; // <-- Should be tuned using subsequence length and either detected cache size or dummy estimate

   while(rem > 0){
      int unroll = std::min(8,rem); // <-- this needs 2 steps of unrolling if we can reasonably stack enough queries to justify parallelizing this part
      int start = qstart+(qcount-rem)*step;
      #pragma omp parallel for
      for(int i = 0; i < unroll; i++){
         batch_normalize(ts,qbuf,mu,start+i*unroll,unroll,sublen,step);  
      }
      #pragma omp parallel for
      for(int i = 0; i < unroll; i++){  // <-- dummy setup
         partial_xcorr_kern(qmcov,qcorr,qbuf,qind,ts,mu,invn,start,step,excl); 
      }
      rem -= unroll;
   }
}

// this may need prefetching due to weird memory access
void max_reduc(double* __Restrict__ qcov, double* __Restrict__ qcorr, int* __Restrict__ qind, double* __Restrict__ corr, int* __Restrict__ ind, int qstart, int count, int stride,int passes){
   for(int i = qstart; i < qstart+count; i++){
      double qcr = qcorr[i];
      double qi  = qind[i];
      for(int j = 0; j < passes; j++){
         if(qcorr[i+j*stride] > qcr){
            qcr = qcorr[i+j*stride];
            qi = qind[i+j*stride];
         }
      }
   }
}


void pxc_extrap_unifstride(double* __Restrict__ qmcov, double* __Restrict__ qcorr, double* __Restrict__ qbuf, int64_t* __Restrict__ qind, const double* __Restrict__ ts, const double* __Restrict__ mu, const double* __Restrict__ df, const double* __Restrict__ dg, const double* __Restrict__ invn, int qstart, int qcount, int step, int mxlen){
   for(int i = qstart; i < qcount; i++){
      //int jmx = std::min(mxlen,i+step);
      int jinit = qind[i];
      double c = 
      for(int j = jinit; j < jinit+step; j++){
         
      }
   } 
}


