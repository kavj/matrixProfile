#include<cstdio>
#include<unistd.h>
#include "mkl.h"




// API may still need work. qbuf_desc describes a set of strided query buffers where each query is currently a 1D subsequence but that may change. We assume the use of padding to maintain SIMD alignment conditions when necessary.
// The first uses duplicated buffers for batch query normalization. This way multiple threads can share read access.
// The second structure provides reusable per thread vector valued variables via striding. In general I don't anticipate conflicts here, as its limited to qbuf_desc*blockct*max_threadcount .. If single threaded, it would be best to just update in place
// The third structure contains all persistent arrays. These are maintained using subdivided access. Note how qcorr isn't needed there. It's because the intended access pattern is to access normalized queries once, process over all buffers, then discard. 
// It would be possible to add an extrapolation pass prior to the next parallel section, thus removing qcov and qmatch completely. It depends on the performance impact. Matches might be scattered in memory, while the queries are somewhat localized due to striding.



int init(int buffercount, int bufferlen, int alignment){
   int padding = bufferlen - (bufferlen - alignment%sizeof(double));
   // I should sanitize input somewhere so that this can't overflow
   if(bufferlen%alignment){
      bufferlen += alignment;
      bufferlen -= bufferlen%alignment;
   }
   // allocate all using aligned allocator if available
   // undecided whether we should use an ad hoc one otherwise
}

// this should have a constructor and destructor. It can be passed as a constant to multithreaded code sections
struct qbuf_desc{
   double* q;  
   int blockstride;
   int blockct;
   int blockstride; 
   int bufct;
};


struct p_autocorrbuf_desc{  
  double* qcov;  // covariance of optimal query match, multiple copies are used in the case of multiple queries
   double* qcorr; // need to know this
   int* qmatch;   // nearest neighbor index for each query
   VSLCorrTaskPtr* covdescs;
   int blocklen;    //length of each correlation block section
   int blockct;
   int blockstride; //stride with respect to correlation buffer, a stride in memory accounts for fringe portions and elements used to pad memory alignment
   // should clean up task pointers as necessary ?

   int init_taskptrs(const double* ts, int len, int blockct,int stride,int sublen){
      covdescs = (VSLCorrTaskPtr*) mkl_malloc(blockct*sizeof(VSLCorrTaskPtr),64);
      for(int i = 0; i < blockct; i++){
         int status = vsldCorrNewTaskX1D(covdescs+i,VSL_CORR_MODE_FFT,blocklen,sublen,blocklen,ts+i*stride,1);
         if(status != VSL_STATUS_OK){
            perror("could not initialize tasks");
            return -1; 
         }
      }
      return 0;
   }  
};


struct p_autocorr_desc{
   double* qcov;
   int* qmatch;
   double* cov;
   double* xcorr;
   int* xind;
   int len;       //length of reduced cross correlation vector
   int sublen;
};
