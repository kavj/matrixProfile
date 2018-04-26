#include<cstdio>
#include<unistd.h>
#include "mkl.h"



static inline int paddedlen(int buflen, int unitsz, int alignmt){
   return  (buflen*unitsz) + (buflen*unitsz)%alignment ? alignnment - (bufferlen*unitsz)%alignmt : 0;
}

// standardized in case we want to use different allocators

void* init_buffer(int bufsz, int alignmt){
   void* buf;
   posix_memalign(buf,buflen*bufct,alignmt);
   return buf;
}

// I'll eventually build in a way to propagate errors from low level libraries to high level bindings. For now I'm mostly ignoring it.
//
VSLCorrTaskPtr* init_taskptrs(const double* ts, int len, int blkct,int blkstride,int sublen){
   VSLCorrTaskPtr* covtsks = init_buffer(1,blockct*sizeof(VSLCorrTaskPtr)*blkct,sizeof(VSLCorrTaskPtr));
   for(int i = 0; i < blockct; i++){
      int status = vsldCorrNewTaskX1D(covtsks+i,VSL_CORR_MODE_FFT,blocklen,sublen,blocklen,ts+i*stride,1);
      if(status != VSL_STATUS_OK){
         // Build up real error checking once we have a less ad hoc architecture in place
         perror("could not initialize tasks");
         return -1; 
      }
   }
   return covtsks;
}

int destroy_taskptrs(VSLCorrTaskPtr* covtsks, int blkct){
   for(int i = 0; i < blockct; i++){
      //int status = 
 
   }
   return 0;
}

template<typename dtype>
struct buf_set{
   buf_set(int blocklen, int blockcount, int alignmt){
      count = blockcount;
      len = blocklen;
      stride = paddedlen(blocklen,sizeof(dtype),alignmt);
      dat = init_buffer(stride*count,alignmt);
   }
   ~buf_set(){
      free(dat);
   }
   
  // inline dtype* operator()(int i) __attribute__((always_inline));
   inline dtype* operator()(int i){ return i < count ? dat + i*stride : nullptr;} __attribute__((always_inline))
   dtype* dat;
   int len;
   int count;
   int stride;
};


//template<typename dtype>
struct pcorrbuf{ 
   pcorrbuf(int querylen, int qbufct, int qstatslen, int qstatsct, int databuflen, int databufct,  int alignmt){
      covbufs = 
   }
   struct buf_set<double> q;
   struct buf_set<double> qcov;
   struct buf_set<double> qcorr;
   struct buf_set<double> covbufs;
   struct buf_set<int> qmatch; 
   VSLCorrTaskPtr* covtsks; // descriptor for MKL
};


struct p_autocorr{
   double* cov;
   double* xcorr;
   int* xind;  
   double* mu;
   double* invn;
   double* df;
   double* dx;
   int len;       //length of reduced cross correlation vector
};


