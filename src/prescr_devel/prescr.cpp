#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"


// This could probably use some cleanup
// centered unit normalization
// this is probably an acceptable way to do it, 
// I was debating whether I should loop over singles in a parallel region
// Perhaps should consult elsewhere?


// somehow I'm not familiar with convention here on default initialization.


struct mkl_descs{
   // mkl corr descriptor array here

};

struct corr_databuffers{
   double* qbuf;  //
   double* qcov;  //
   double* qcorr; // need to know this
   double* qmu;   // may be helpful when 
   double* qinvn; //
   int* qbase;    //
   int* qmatch;   //
   double* cbuf;  // 
   double* xcorr; // cross correlation 
   int* xcorrind; // cross correlation index
};


struct corr_desc{
   int qbufcount; // number of query buffers
   int qtotalcount;//number of batch queries 
   int qlen;       //query / subsequence length
   int qstride;    //stride with respect to time series (distribution of queries)
   int qbufstride; //query stride with respect to query buffer
   int cbufcount;  //number of active correlation buffers
   int ctotalcount;
   int clen;       //total number of cross correlations
   int cstride;    //stride with respect to time series
   int cbufstride; //stride with respect to conv buffer
  
   // example, functions should be like this
   // nullptr could signify an error, but I'm not sure whether this is the best idea
   inline int qbuffer(int i){
      return i*qbufstride;
   }
   inline int cbuffer(int i){
      return i*cbufstride;
   }
};


struct prescr_desc{
   //prescr_desc() : ts(nullptr),qcov(nullptr),qcorr(nullptr),qind(nullptr),invn(nullptr),mp(nullptr),mpi(nullptr),len(0),sublen(0),qcount(0) {};
   double* ts;
   double* invn;

   void alloc(){

   }
   void free_buffers(){

   }
   void free_all(){

   }
  

};


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
         corr(j) = cov[i+j]*invn[i+j];
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
      double corr = cov[i]*invn[i];
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
   return ; // qcorr, qcov, qind, 
}




// The time complexity of this is linear or close to linear, so we can get it working and ignore
// I changed it to a zigzag pattern. 
//
void maxpearson_extrap_partialauto(const double* __restrict__ qcov, const double* __restrict__ invn, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len){
   for(int i = 0; i < count; i++){
      int qi = qind[i];
      int ci = i*stride;
      int lim = std::min(ci,qi);
      lim = std::min(extraplen,lim);
      double c = qcov[i];
      for(int j = 0; j < lim; j++){
      //   c -= dx[ci-j]*df[qi-j];
      //   c -= dx[qi-j]*df[ci-j];
         double scale = invn[qi]*invn[ci];
         double corr = c*scale;
         if(mp[ci] < corr){
            mp[ci] = corr;
            mpi[ci] = qi;
         }
         if(mp[qi] < corr){
            mp[qi] = corr;
            mpi[qi] = ci;
         }
         --ci;
         --qi;
      }
      lim = std::max(ci,qi);
      lim = std::min(extraplen,len-lim);
      for(int j = 0; j < lim; j++){
      // use the original formulas instead, this section doesn't justify more buffers
 
      //   c += dx[ci]*df[qi];
     //    c += dx[qi]*df[ci];
         double scale = invn[qi]*invn[ci];
         double corr = c*scale;
         if(mp[ci] < corr){
            mp[ci] = corr;
            mpi[ci] = qi;
         }
         if(mp[qi] < corr){
            mp[qi] = corr;
            mpi[qi] = ci;
         }
         ++ci; 
         ++qi;
      }
   }
}



// if it's not parallel, we allocate fewer buffers assume parallel for interactivity
void prescr_exec_partialauto(const struct corr_desc* __restrict__ crd, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn){
//   int qmemstride = 64; // sentinel value, should be on an appropriate boundary for MKL
//   int qbufcount = 10;
//   int kstride = 65536;
//   int kmemstride = 65536; // since correlation op outputs n + m - 1, we may need to round things or maybe just use separate buffers? not sure here really what is optimal
   // initialize corr descriptors
   for(int i = 0; i < crd->qbufcount; i++){  // replace with qtotal
      #pragma omp parallel for
      for(int j = 0; j < crd->qbufcount; j++){
         double qm = mu[j*qstride];
         for(int k = 0; k < sublen; k++){
            crd->qbuf[j*qmemstride+k] = ts[j*qstride+k] - qm;
         }
      }
      #pragma omp parallel for
      for(int j = 0; j < ccount; ){
         for(int k = 0; k < qbufcount; k++){
           //
           // corr(ts+j*kstride,
           // reduce
            
         }
      }
   }
   // refine stage
}







