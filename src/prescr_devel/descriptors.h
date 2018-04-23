#include<cstdio>
#include<unistd.h>
#include "mkl.h"





struct qbuf{
   // Todo: explicit constructor and destructor. These can be assigned per block correlation section, so each section tracks a set of best so far reductions for the query set.
   // We require as many descriptors as max threads
   double* q;  
   int qbufcount;
   int qbufstride;
   double* qcov;  // covariance of optimal query match, multiple copies are used in the case of multiple queries
   double* qcorr; // need to know this
   int* qbase;    // stores the base index of each query with respect to the time series from which it was sampled
   int* qmatch;   // nearest neighbor index for each query
   int numcopies;    // since a lot of this is effectively read only (could enforce that) it is just strided, We don't allocate or deallocate in a shared section, and we avoid any cache line overlaps
};

struct xcorr_desc{
   double* xcorr;
   int* xind;
   VSLCorrTaskPtr* covdescs;
   int seccount;
   int len;       //length of reduced cross correlation vector
   int seclen;    //length of each correlation block section
   int bufstride; //stride with respect to correlation buffer, a stride in memory accounts for fringe portions and elements used to pad memory alignment
   int taillen;
   int subseclen;
   // should clean up task pointers as necessary ?
};


int init_taskptrs(struct xcorr_desc xcd, const double* ts){
   int count = xcd.seccount;;
   VSLCorrTaskPtr* covdescs = (VSLCorrTaskPtr*) mkl_malloc(count*sizeof(VSLCorrTaskPtr),64);
   for(int i = 0; i < count; i++){
      int status = vsldCorrNewTaskX1D(covdescs+i,VSL_CORR_MODE_FFT,xcd.seclen,xcd.subseclen,xcd.seclen,ts+i*xcd.bufstride,1);
      if(status != VSL_STATUS_OK){
         perror("could not initialize tasks");
         return -1; 
      }
   }
   xcd.covdescs = covdescs;
   return 0;
}


