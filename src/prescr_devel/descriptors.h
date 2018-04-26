#include<cstdio>
#include<unistd.h>
#include "mkl.h"




// API may still need work. qbuf_desc describes a set of strided query buffers where each query is currently a 1D subsequence but that may change. We assume the use of padding to maintain SIMD alignment conditions when necessary.
// The first uses duplicated buffers for batch query normalization. This way multiple threads can share read access.
// The second structure provides reusable per thread vector valued variables via striding. In general I don't anticipate conflicts here, as its limited to qbuf_desc*blockct*max_threadcount .. If single threaded, it would be best to just update in place
// The third structure contains all persistent arrays. These are maintained using subdivided access. Note how qcorr isn't needed there. It's because the intended access pattern is to access normalized queries once, process over all buffers, then discard. 
// It would be possible to add an extrapolation pass prior to the next parallel section, thus removing qcov and qmatch completely. It depends on the performance impact. Matches might be scattered in memory, while the queries are somewhat localized due to striding.



void* init_buffer(int bufct, int buflen, int alignmt){
   int padding = buflen%alignment ? alignnment -  (bufferlen%alignmt) : 0;
   buflen += alignmt;
   void* buf;
   posix_memalign(buf,buflen*bufct,alignmt);
   // use posix_memalign here or something else?
   // VC++ uses unaligned instructions anyway
   // allocate all using aligned allocator if available
   // undecided whether we should use an ad hoc one otherwise
   return buf;
}


// I'll eventually build in a way to propagate errors from low level libraries to high level bindings. For now I'm mostly ignoring it.
//
int init_taskptrs(const double* ts, int len, int blkct,int blkstride,int sublen){
   VSLCorrTaskPtr* covdescs = init_buffer(1,blockct*sizeof(VSLCorrTaskPtr)*blkct,sizeof(VSLCorrTaskPtr));
   for(int i = 0; i < blockct; i++){
      int status = vsldCorrNewTaskX1D(covdescs+i,VSL_CORR_MODE_FFT,blocklen,sublen,blocklen,ts+i*stride,1);
      if(status != VSL_STATUS_OK){
         // Build up real error checking once we have a less ad hoc architecture in place
         perror("could not initialize tasks");
         return -1; 
      }
   }
   return 0;
}


// this should have a constructor and destructor. It can be passed as a constant to multithreaded code sections
struct qbuf_desc{
   double* q;  
   int blklen;
   int blkstrd;
   int blockct;
};


/*
 * this stores information for multiple partial correlation operations. 
 * this doesn't indicate whether qcorr is persistent. I'll deal with that later.
 * I assume yes for now unless I can discover that nothing depends on it
 */
struct partial_corr_buf_desc{ 
   double* qcov;              
   double* qcorr;            // nearest neighbor candidates for each query
   int* qmatch;              // nearest neighbor index for each query
   int qblkstrd;             // stride between consecutive queries
   int qblkct;               // number of queries which may be queued
   int qlen;                 // query length
   VSLCorrTaskPtr* covdescs; // descriptor for MKL
   int blklen;               //untruncated length of each correlation block section
   int blkct;                //number of queries that may be buffered at any given time
   int dstrd;                //stride between data blocks. We can't simply use buffer length, because we need to truncate edges
   int bufstrd;              //stride between consecutive packed buffers   
};


struct p_ac_desc{
   double* cov;
   double* xcorr;
   int* xind;
   int len;       //length of reduced cross correlation vector
};


