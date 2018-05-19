#include<algorithm>
#include "../utils/xprec_math.h"
#include "../utils/cov.h"
//#include "prescr_Pearson.h"
#include "auto_Pearson.h"
#define tsz 16384
#define ksz 248
#define prefalign 64

// should be using namespace whatever here to determine whether an optimized kernel may be used
// // dtype and itype should be cleaned up too. The rules are reference kernel itype is always base integer
// // optimized kernel and double means itype needs to be int64_t to avoid extra swizzles, and single precision 
// // itype should always be standard integer type. We don't really need to check further, because the compiler is the only thing
// // that can vectorize reference kernels. For hand tuned ones, we know the size of the data types anyway.
//
// // We can index here in terms of smaller tiles, where we use the fact that x number of smaller ones fit one large

//struct query_stat{double val; long long ind};

// This should be tested. If initialization overhead isn't an issue and we're using double precision,


template<typename dtype, typename itype>
int maxpearson_partialauto(stridedbuf<dtype>& ts, int minlag, int mlen, int sublen);

template<typename dtype, typename itype>
int maxpearson_partialauto(stridedbuf<dtype>& ts, int minlag, int mlen, int sublen){
   // allocate query buffers, basically covariance and normalized query
   const int step = tsz/ksz;
   int tilesperdim = mlen/ksz + (mlen%ksz ? 1 : 0); // find max tiles in 1 direction
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }
   int n = ts.len;
   stridedbuf<dtype> mu(tilesperdim,ksz,mlen);
   stridedbuf<dtype> invn(tilesperdim,ksz,mlen);
   stridedbuf<dtype> df(tilesperdim,ksz,mlen);
   stridedbuf<dtype> dx(tilesperdim,ksz,mlen);
   stridedbuf<dtype> mp(tilesperdim,ksz,mlen);
   stridedbuf<dtype>cov(ksz,n,n);
   stridedbuf<dtype>q(tilesperdim,paddedlen(sublen,prefalign));
   stridedbuf<itype> mpi(tilesperdim,ksz,mlen);
   if(!(q.isvalid() && cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dx.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");   
   } 
   xmean_windowed(ts(0),mu(0),ts.len,sublen);
   xsInv(ts(0),mu(0),invn(0),ts.len,sublen);   
   init_dfdx(ts(0), mu(0), df(0), dx(0),sublen,ts.len);
   std::fill(mp(0),mp(0)+mp.len,-1.0);
   std::fill(mpi(0),mpi(0)+mpi.len,-1); 
    
   //Todo: fix centering of query        
   #pragma omp parallel for
   for(int i = 0; i < tilesperdim; i++){
      center_query(ts(i), mu(i), q(i), sublen);
   }
   return 0;
   for(int i = 0; i < tilesperdim; i++){
      #pragma omp parallel for
      for(int j = 0; j < tilesperdim-i; j++){
         if((i+j+2) < tilesperdim){
             batchcov(ts(j),cov(j),q(j),mu(i+j),tsz,sublen);
            int cindoffset = i*tsz+minlag;
            int qbaseind = j*tsz; 
            double qinvn = invn(j)[0];
            double qcorr = -1.0; 
            itype qmatch = -1;
            int offset = (i+j)*tsz;
            for(int k = 0; k < tsz; k++){
               for(int l = 0; l < tsz; l++){
                  pauto_pearson_kern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz); 
               }
            }
         }
         else{
            continue;
            int offset = (i+j)*tsz;
            batchcov(ts(j),cov(j),q(j),mu(i+j),tsz,sublen);
            int cindoffset = i*tsz+minlag;
            int qbaseind = j*tsz; 
            double qinvn = invn(j)[0];
            double qcorr = -1.0; 
            itype qmatch = -1;
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
 
