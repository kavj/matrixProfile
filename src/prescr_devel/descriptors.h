//#include<cstdio>
//#include<unistd.h>
#include "alloc.h"
#include "mkl.h"

// This should probably be temporary. Since I end up having to distribute stuff frequently, it would be a good idea to generalize the striding somewhat for 1 and 2D strides


static inline int paddedlen(int buflen, int unitsz, int alignmt){ return  (buflen*unitsz) + ((buflen*unitsz)%alignment ? alignnment - (bufferlen*unitsz)%alignmt : 0); } __attribute__((always_inline))

template<typename dtype>
struct twod_buf{
   twod_buf(int blocklen, int blockcount, int alignmt){
      count = blockcount;
      len = blocklen;
      stride = paddedlen(blocklen,sizeof(dtype),alignmt);
      //dat = init_buffer(stride*count,alignmt);
   }

   ~buf_set(){
      free(dat);
   }
   
   inline dtype* operator()(int i){ return i < count ? dat + i*stride : nullptr;} __attribute__((always_inline))
   // inline int blocklen(int i){return i < bcount-1 ? blen : blen*bcount
   // could add overall length here to implement this, or we could just precompute last section len and make it a struct member  
   dtype* dat;
   int blen;   
   int bcount;
   int bstride; 
};



//template<typename dtype>
struct corr_auxbuf{ 
   corr_auxbuf(int qlen, int qbufct, int qstatslen, int qstatsbufct, int databuflen, int databufct,  int alignmt){
      
   }
   struct twod_buf<double> q;
   struct twod_buf<double> qcov;
   struct twod_buf<double> qcorr;
   struct twod_buf<double> covbufs;
   struct twod_buf<int> qmatch; 
   VSLCorrTaskPtr* covtsks; // descriptor for MKL
   int querycount;   // total queries
   int qbasestride;  // indicates distance between queries with respect to a time series, can be set to -1 if these are not uniformly sampled 
   int tailcount;  // <-- need a better name, but this would be the unaligned segment
                       // this is REM(length(time_series)/length(buffer))
   // for covbufs, need way to identify the fringe component. Perhaps  
   inline int isinitialized(){return (q.dat != nullptr) && (qcov.dat != nullptr) && (qcorr.dat != nullptr) && (covbufs.dat != nullptr) && (qmatch.dat != nullptr);} __attribute__((always_inline))
   
/*
 *  Ultimately need to be able to do something like reduce remaining from length(data) to 0 somewhere. 
 *  
*/
 
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
   int bstride;
   int bcount;
   int blen;
   int taillen;
   // Todo: the resource handle stuff should be much more opaque or it might cause problems later. This is fine for now given that strides are uniform. It is important to note that sections might overlap here on read only memory
   // due to the use of overlap and save style filtering

   inline int operator()(int i){ return i < bcount ? i*bstride : nullptr; } __attribute__((always_inline))
   inline int isinitialized(){return (ts != nullptr) && (cov != nullptr) && (xcorr != nullptr) && (xind != nullptr) && (mu != nullptr)  && (inv != nullptr) && (df != nullptr) && (dx != nullptr); } __attribute__((always_inline))

};

struct query_stat{
   double qcov;
   double qcorr;
   int qind;
   double qinvn;
   int qbaseind;
};

