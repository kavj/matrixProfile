#include<algorithm>
#include "../utils/xprec_math.h"
#include "../utils/cov.h"
//#include "prescr_Pearson.h"
#include "auto_Pearson.h"
#define tsz 65536 
#define ksz 248
#define prefalign 64

   // The rest needs to be a decision between full and partial tile. The full tile should instantiate a customizable size at runtime
   // The partial tile has a shared constraint. 
   // This is to be versioned as follows.
   // Symmetric auto, symmetric cross (requires equal length), rc_auto, rc_cross. This can be passed in as compile time options.
   // This function is the front end and should only deal with preprocessing and dispatch. 
   // It should hand off to a solver for the low level work. At the lowest level we have kernels which should be fully unrolled. 
 

template<typename dtype, typename itype>
int maxpearson_partialauto(stridedbuf<dtype>& ts, int minlag, int mlen, int sublen);

template<typename dtype, typename itype>
int maxpearson_partialauto(stridedbuf<dtype>& ts, int minlag, int mlen, int sublen){
   // allocate query buffers, basically covariance and normalized query
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }
   int tilesperdim = (mlen-minlag)/tsz; // find max tiles in 1 direction
   int taillen = mlen - minlag - tsz*tilesperdim;
   if(taillen > 0){
      tilesperdim++;
   }

   stridedbuf<dtype>mu(mlen); stridedbuf<dtype>invn(mlen); stridedbuf<dtype>df(mlen);  stridedbuf<dtype>dx(mlen);
   stridedbuf<dtype>mp(mlen); stridedbuf<dtype>cov(mlen);  stridedbuf<itype>mpi(mlen); stridedbuf<dtype>q(tilesperdim,sublen);

   if(!(q.isvalid() && cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dx.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");
      return 0;
   } 
  
   ts.setstride(tsz), mu.setstride(tsz), df.setstride(tsz), dx.setstride(tsz), invn.setstride(tsz),
   mp.setstride(tsz), cov.setstride(tsz), mpi.setstride(tsz);

   xmean_windowed(ts(0),mu(0),ts.len,sublen);

   xsInv(ts(0),mu(0),invn(0),ts.len,sublen);   
   init_dfdx(ts(0), mu(0), df(0), dx(0),sublen,ts.len);
   std::fill(mp(0),mp(0)+mp.len,-1.0);
   std::fill(mpi(0),mpi(0)+mpi.len,-1); 

   const int step = tsz/ksz;;
   const int tstride = 65536;
   const int tcount = ts.len/tstride;  
   
 
   #pragma omp parallel for
   for(int i = 0; i < tcount; i++){
      center_query(ts(i), mu(i), q(i), sublen);
   }
   for(int i = 0; i < tilesperdim; i++){
      #pragma omp parallel for
      for(int j = 0; j < tilesperdim-i; j++){
         int offset = (i+j)*tsz;
         batchcov(ts(j),cov(j),q(j),mu(i+j),tsz,sublen);
         
      }
   }
}
 
