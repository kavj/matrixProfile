#include<cstdio>
#include<unistd.h>
#include "mkl.h"


// Todo: figure out a reasonable way to turn factoor out separate back ends for mkl vs quasi Cook-Toom filtration


struct qdesc{
   int qind;
   double qcorr;
   double qcov;
};

struct qbuf{
   // this should probably have an explicit constructor and destructor

   qbuf() : q(nullptr), qbufcount(0), qbufstride(0){}

   double* q;  
   int qbufcount;
   int qbufstride;
  
   // should have an explicit destructor

   /*const double* getquery(int i){
      if(i <= qbufcount){
         return q+i*qbufstride;
      }
      else{
         return nullptr;
      }
   }*/

};

struct buffer{
   

};


struct qstatbuf{
   // Todo: explicit constructor and destructor. These can be assigned per block correlation section, so each section tracks a set of best so far reductions for the query set.
   // We require as many descriptors as max threads
   
   double* qcov;  // covariance of optimal query match, multiple copies are used in the case of multiple queries
   double* qcorr; // need to know this
   int* qbase;    // stores the base index of each query with respect to the time series from which it was sampled
   int* qmatch;   // nearest neighbor index for each query

};

struct xcorr_desc{
   double* xcorr;
   int* xcorrind;
   int ccount;     //number of elements in correlation vector, in our cases (we truncate aliased portions) length(time_series) - length(query) + 1
   int clen;       //length of reduced cross correlation vector
   int cseclen;    //length of each correlation block section
   int cbufstride; //stride with respect to correlation buffer, a stride in memory accounts for fringe portions and elements used to pad memory alignment
   int sublen;
   VSLCorrTaskPtr* covdescs;


   // should clean up task pointers as necessary ?

};




int init_taskptrs(VSLCorrTaskPtr* covdescs, double* ts, int clen, int seclen, int sublen){
   int count = clen/cseclen;
   if(count*cseclen < clen){
      ++cseclen;
   }
   covdescs = (VSLCorrTaskPtr*) mkl_malloc(count*sizeof(VSLCorrTaskPtr),64);
   for(int i = 0; i < ccount; i++){
      int status = vsldCorrNewTaskX1D(covdescs+i,VSL_CORR_MODE_FFT,cseclen,sublen,cseclen,ts+i*cbufstride,1);
      if(status != VSL_STATUS_OK){
         perror("could not initialize tasks");
         return -1; 
      }
   }
   return 0;
}




