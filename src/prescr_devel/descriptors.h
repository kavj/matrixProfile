//#include<cstdio>
//#include<unistd.h>
#include "alloc.h"
#include "mkl.h"



static inline int paddedlen(int buflen, int unitsz, int alignmt){
   return  (buflen*unitsz) + (buflen*unitsz)%alignment ? alignnment - (bufferlen*unitsz)%alignmt : 0;
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
   
   inline dtype* operator()(int i){ return i < count ? dat + i*stride : nullptr;} __attribute__((always_inline))
   
   dtype* dat;
   int blen;   
   int bcount;
   int bstride; 
};


//template<typename dtype>
struct corr_auxbuf{ 
   corr_auxbuf(int querylen, int qbufct, int qstatslen, int qstatsct, int databuflen, int databufct,  int alignmt){
      
   }
   struct buf_set<double> q;
   struct buf_set<double> qcov;
   struct buf_set<double> qcorr;
   struct buf_set<double> covbufs;
   struct buf_set<int> qmatch; 
   VSLCorrTaskPtr* covtsks; // descriptor for MKL
   int querycount;   // total queries
   int querystride;  // indicates distance between queries
   int unalsegment;  // <-- need a better name, but this would be the unaligned segment
                     // this can be something like p_autocorr.len - blen*(bcount - 1) 
   // for covbufs, need way to identify the fringe component. Perhaps  
   
};


struct p_autocorr{
   double* ts;
   double* cov;
   double* xcorr;
   int* xind;  
   double* mu;
   double* invn;
   double* df;
   double* dx;
   int len;  
   int xlen; 
};

struct query_stat{
   double qcov;
   double qcorr;
   int qind;
   double qinvn;
   int qbaseind;
};

