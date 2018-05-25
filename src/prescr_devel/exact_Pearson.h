#include <cstdio>
#include<algorithm>
#include "descriptors.h"
#include "../utils/xprec_math.h"
#include "../utils/cov.h"
#include "tiled_Pearson.h"
#include "../utils/primitive_print_funcs.h"
#include "../utils/checkArray.h"
#define prefalign 64

   // The rest needs to be a decision between full and partial tile. The full tile should instantiate a customizable size at runtime
   // The partial tile has a shared constraint. 
   // This is to be versioned as follows.
   // Symmetric auto, symmetric cross (requires equal length), rc_auto, rc_cross. This can be passed in as compile time options.
   // This function is the front end and should only deal with preprocessing and dispatch. 
   // It should hand off to a solver for the low level work. At the lowest level we have kernels which should be fully unrolled. 
 


template<typename dtype, typename itype>
void maxpearson_partialauto(stridedbuf<dtype>& ts, stridedbuf<dtype>& mp, stridedbuf<itype>& mpi, int minlag, int sublen){
   // allocate query buffers, basically covariance and normalized query
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }

   int mlen = ts.len - sublen + 1;
   const int tlen = std::max(paddedlen(4*sublen,prefalign),131072); // this should incorporate the overall size
   int tilesperdim = (mlen-minlag)/tlen; // find max tiles in 1 direction
   int taillen = mlen - minlag - tlen*tilesperdim;

   if(taillen > 0){
      tilesperdim++;
   }

   stridedbuf<dtype>mu(mlen); stridedbuf<dtype>invn(mlen); stridedbuf<dtype>df(mlen);  stridedbuf<dtype>dg(mlen);
   stridedbuf<dtype>q(tilesperdim,sublen); stridedbuf<dtype>cov(mlen);

   if(!(q.isvalid() && cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dg.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");
      return;
   } 
  
   ts.setstride(tlen), mu.setstride(tlen), df.setstride(tlen), dg.setstride(tlen), invn.setstride(tlen),
   mp.setstride(tlen), cov.setstride(tlen), mpi.setstride(tlen);

   xmean_windowed(ts(0),mu(0),ts.len,sublen);
   xsInv(ts(0),mu(0),invn(0),ts.len,sublen);   
   
//   fast_invcn(invn(0),ts(0),mu(0),len,sublen);
   init_dfdx(ts(0), mu(0), df(0), dg(0),sublen,ts.len);

   std::fill(mp(0),mp(0)+mlen,-1.0);
   std::fill(mpi(0),mpi(0)+mlen,-1); 

   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i), mu(i), q(i), sublen);
   }
  
   for(int i = 0; i < tilesperdim; i++){
      // this should contain a check for if ! last iteration
      #pragma omp parallel for
      for(int j = 0; j < tilesperdim-i; j++){
         int mx = mlen-minlag-(i+j)*tlen;
         int maxperdim = mx > tlen ? tlen : mx;
         int upperbound = mx > 2*tlen ? 2*tlen : mx;
         batchcov_ref(ts(i+j)+minlag,cov(j),q(j),mu(j),maxperdim,sublen);
         pauto_pearson(cov(j),mp(j),mpi(j),df(j),dg(j),invn(j),maxperdim,j*tlen,i*tlen+minlag,upperbound);
      }
   }
}
 
