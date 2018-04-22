#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"
#include "descriptors.h"




void fast_invcn(double* __restrict__ invn, const double* __restrict__ ts, const double* __restrict__ mu, int len, int sublen){
   double a = 0;
   for(int i = 0; i < sublen; i++){
      double term = ts[i] - mu[0];
      a += term*term;
   }
   invn[0] = 1.0/a;
   for(int i = sublen; i < len-sublen+1; i++){
      double b = ts[i-sublen];
      double c = ts[i];
      a += ((b - mu[i-sublen]) + (c - mu[i])) * (b - c); 
      invn[i-sublen] = 1.0/a;
   }
}

// mmay not be needed
void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int qstart, int count, int sublen, int step){
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qind = i*sublen;
      int tsind = i*step+qstart;
      double qm = mu[qind];
      double qinvn = invn[qind];
      for(int j = 0; j < sublen; j++){
         qbuf[qind] = (ts[tsind]-qm)*qinvn;
         ++qind;
         ++tsind;
      }
   }
}


// This is probably the best I can do without avx
void max_reduction(double* __restrict__ cov,  double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, double qinvn, double qcov, double qcorr, int qind, int qbaseind, int offset, int count){
   const int unroll = 8;
   int aligned = count - count%unroll + offset;
   for(int i = offset; i < aligned; i+= unroll){
      block<double> corr;
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*qinvn;;
      }
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*invn[i+j];
      }
      for(int j = 0; j < unroll; j++){
         if(xcorr[i+j] < corr(j)){
            xcorr[i+j] = corr(j);
            cindex[i+j] = qbaseind;
         }
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
            qind = i+cind(4);
            qcov = cov[i+cind(4)];
         }
      }
      else if(qcorr < corr(0)){
         qcorr = corr(0);
         qind = cov[i+cind(0)];
      }
   }
   for(int i = aligned; i < offset+count; i++){
      double corr = cov[i]*qinvn*invn[i];
      if(qcorr < corr){
         qcorr = corr;
         qind =  i;  //<-- this way it's only forward iteration. Doing both consecutively may partially invalidate  hardware prefetching
         qcov = cov[i];
      } 
      if(xcorr[i] < corr){
         xcorr[i] = corr;
         cindex[i] = qbaseind;
      }
   }
   return;

}


// need to add structs in here
// This might be reworked later to remove temporary buffers or to wrap the messy formulas in small inline functions, not sure really but it's too hard to read with inline formulas
void maxpearson_extrap_partialauto(const double* __restrict__ qcov, const double* __restrict__ invn, const double* __restrict__ df, const double* __restrict__ dx, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len){
   for(int i = 0; i < count; i++){
      int qi = qind[i];
      int ci = i*stride;
      int lim = std::min(extraplen,std::min(ci,qi));
      double c = qcov[i];
      for(int j = 1; j < lim; j++){
         c -= dx[ci]*df[qi];
         c -= dx[qi]*df[ci];
         double corr = c*invn[qi]*invn[ci];
         if(mp[ci-j] < corr){
            mp[ci-j] = corr;
            mpi[ci-j] = qi-j;
         }
         if(mp[qi] < corr){
            mp[qi-j] = corr;
            mpi[qi-j] = ci-j;
         }
      }
      lim = std::min(extraplen,len-std::max(ci,qi));
      c = qcov[i];
      for(int j = 1; j < lim; j++){
         c += dx[ci]*df[qi];
         c += dx[qi]*df[ci]; 
         double corr = c*invn[qi]*invn[ci];
         if(mp[ci+j] < corr){
            mp[ci+j] = corr;
            mpi[ci+j] = qi+j;
         }
         if(mp[qi+j] < corr){
            mp[qi+j] = corr;
            mpi[qi+j] = ci+j;
         }
      }
   }
}


// if it's not parallel, we allocate fewer buffers assume parallel for interactivity
void prescr_exec_partialauto(const struct corr_desc* __restrict__ crd, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn){

  // initialize corr descriptors
   for(int i = 0; i < crd->qbufcount; i++){  // replace with qtotal
      #pragma omp parallel for
      for(int j = 0; j < crd->qbufcount; j++){
         double qm = mu[j*crd->qstride];
         for(int k = 0; k < crd->querylen; k++){
           // crd->qbuf[j*qmemstride+k] = ts[j*qstride+k] - qm;
         }
      }
      #pragma omp parallel for
      for(int j = 0; j < crd->ccount; ){
         for(int k = 0; k < crd->qbufcount; k++){
           //
           // corr(ts+j*kstride,
           // reduce
            
         }
      }
   }
   // refine stage
}


