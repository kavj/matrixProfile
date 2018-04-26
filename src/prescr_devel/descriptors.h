#include<cstdio>
#include<unistd.h>
#include "mkl.h"




// API may still need work. qbuf_desc describes a set of strided query buffers where each query is currently a 1D subsequence but that may change. We assume the use of padding to maintain SIMD alignment conditions when necessary.
// The first uses duplicated buffers for batch query normalization. This way multiple threads can share read access.
// The second structure provides reusable per thread vector valued variables via striding. In general I don't anticipate conflicts here, as its limited to qbuf_desc*blockct*max_threadcount .. If single threaded, it would be best to just update in place
// The third structure contains all persistent arrays. These are maintained using subdivided access. Note how qcorr isn't needed there. It's because the intended access pattern is to access normalized queries once, process over all buffers, then discard. 
// It would be possible to add an extrapolation pass prior to the next parallel section, thus removing qcov and qmatch completely. It depends on the performance impact. Matches might be scattered in memory, while the queries are somewhat localized due to striding.

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

int dest_taskptrs(VSLCorrTaskPtr* covtsks, int blkct){
   for(int i = 0; i < blockct; i++){
     //check this part
     // int status = vsldCorr
     //
   }
   return 0;
}

template<typename dtype>
struct buf_set{
   buf_set(int blocklen, int blockcount, int alignmt){
      count = blockcount;
      len = paddedlen(blocklen,sizeof(dtype),alignmt);
      dat = init_buffer(len*count,alignmt);
   }
   ~buf_set(){
      free(dat);
   }
   
   inline dtype* operator()(int i) __attribute((always_inline))__;

   inline dtype* operator()(int i){ return i < count ? dat + i*stride : nullptr;}

   dtype* dat;
   int len;
   int count;
   int stride;
};


// this should have a constructor and destructor. It can be passed as a constant to multithreaded code sections
struct qbuf{
   qbuf(int qlen, int qct, alignmt){
      
      blklen  = paddedlen(qlen,sizeof(double),alignmt);
      q = init_buffer(blklen*qct,alignmt);
      blockct = qct;
   }
   ~qbuf(){
      free(q);
   }
   double* gq(int i){ return i < qct ? q + i*blkstrd : nullptr;}
   double* q;  
   int blklen;
   int blkstrd;
   int blockct;
};


/*
 * this stores information for multiple partial correlation operations. 
 * this doesn't indicate whether qcorr is persistent. I'll deal with that later.
 * I assume yes for now unless I can discover that nothing depends on it
 *
 *
 * This needs so much cleanup git --facedesk
 * The objects should be able to initialize memory structures. I might change this to use lazy instantiation later, but for now it should at least be a bit cleaner
 */
//template<typename dtype>
struct pcorrbuf{ 
   pcorrbuf(int dct, int dbufl, int qct, int ql, int alignmt){
      qlen = ql;
      dct = dbufl;
      qbufstrd = paddedlen(ql,sizeof(double),alignmt);
      matchstrd = paddedlen(ql,sizeof(int),alignmt);
      bufstrd = paddedlen(dbufl,sizeof(double),alignmt);
      qcov = init_buffer(qbufstrd*qlen,alignmt);
      qcorr = init_buffer(qbufstrd*qlen,alignmt);
      qmatch = init_buffer(matchstrd*qlen,alignmt);
      covtsks = nullptr; 
   }

   ~pcorrbuff(){
      free(qbufstrd);
      free(matchstrd);
      free(buffstrd);
      free(qcov);
      free(qcorr);
      free(qmatch);
      /// delete taskptrs
   }

   // let's see if the compiler can inline this

   inline double gqcov(int i){ return    i < qbufct ? qcov    + i*qbufstrd  : nullptr; }
   inline double gqcorr(int i){ return   i < qbufct ? qcorr   + i*qbufstrd  : nullptr; }
   inline double gcovtsks(int i){ return i < qbufct ? covtsks + i           : nullptr; }
   inline double gqmatch(int i){ return  i < qbufct ? qcov    + i*matchstrd : nullptr; }

   double* qcov;             // nearest neighbor to query 
   double* qcorr;            // nearest neighbor candidates for each query
   int* qmatch;              // nearest neighbor index for each query
   VSLCorrTaskPtr* covtsks; // descriptor for MKL
   int matchstrd;
   int qbufstrd;             // stride between consecutive query buffer instances
   int qbufct;               // number of buffer instances per query
   int qlen;                 // query length
   int blklen;               //untruncated length of each correlation block section
   int blkct;                //number of queries that may be buffered at any given time
   int dstrd;                //stride between data blocks. We can't simply use buffer length, because we need to truncate edges
   int dbufstrd;              //stride between consecutive packed buffers   
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


