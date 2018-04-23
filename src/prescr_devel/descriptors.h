#include<cstdio>
#include<unistd.h>
#include "mkl.h"



int init(int buffercount, int bufferlen, int alignment){
   int padding = bufferlen - (bufferlen - alignment%sizeof(double));
   // I should sanitize input somewhere so that this can't overflow
   if(bufferlen%alignment){
      bufferlen += alignment;
      bufferlen -= bufferlen%alignment;
   }
   // allocate all using aligned allocator
}

// this should have a constructor and destructor. It can be passed as a constant to multithreaded code sections
struct qbuf{
   double* q;  
   int len;
   int blockct;
   int blockstride;
};

struct qstats{
   double* qcov;  // covariance of optimal query match, multiple copies are used in the case of multiple queries
   double* qcorr; // need to know this
   int* qmatch;   // nearest neighbor index for each query
   int blockstride;
};

struct p_autocorr_desc{  
   double* xcorr;
   int* xind;
   VSLCorrTaskPtr* covdescs;
   double* cov;
   int blockct;
   int len;       //length of reduced cross correlation vector
   int blocklen;    //length of each correlation block section
   int blockstride; //stride with respect to correlation buffer, a stride in memory accounts for fringe portions and elements used to pad memory alignment
   int taillen;   // <-- sections can be long so this is unpadded
   int sublen;
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


