#include<algorithm>
#include "../utils/xprec_math.h"
#include "../utils/cov.h"
//#include "prescr_Pearson.h"
#include "auto_Pearson.h"
#define tsz 65536 
#define ksz 248
#define prefalign 64

template<typename dtype, typename itype>
int maxpearson_partialauto(stridedbuf<dtype>& ts, int minlag, int mlen, int sublen);

template<typename dtype, typename itype>
int maxpearson_partialauto(stridedbuf<dtype>& ts, int minlag, int mlen, int sublen){
   // allocate query buffers, basically covariance and normalized query
   int tilesperdim = mlen/ksz + (mlen%ksz ? 1 : 0); // find max tiles in 1 direction
   if(!ts.isvalid()){
      printf("invalid time series\n");
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
   
   // currently segfaults after this point
 
   #pragma omp parallel for
   for(int i = 0; i < tcount; i++){
      center_query(ts(i), mu(i), q(i), sublen);
   }
   for(int i = 0; i < tilesperdim; i++){
      #pragma omp parallel for
      for(int j = 0; j < tilesperdim-i; j++){
         int offset = (i+j)*tsz;
         if((i+j+2) < tilesperdim){
            batchcov(ts(j),cov(j),q(j),mu(i+j),tsz,sublen);
            for(int k = 0; k < step; k++){
               for(int l = 0; l < step; l++){
                  pauto_pearson_kern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz); 
               }
            }
         }
         else{
            batchcov(ts(j),cov(j),q(j),mu(i+j),tsz,sublen);
            int dlim = 0; 
            for(int k = 0; k < tsz; k++){
               for(int l = 0; l < tsz; l++){
                  pauto_pearson_edgekern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz,mlen-(offset+l)*ksz);
               }
            }
         }
      }
   }
}
 
