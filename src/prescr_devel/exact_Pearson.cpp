#include<cstdio>
#include<algorithm>
#include "descriptors.h"
#include "exact_Pearson.h"


//template<>
//int maxpearson_partialauto<double,long>(const stridedbuf<double>& ts, int minlag, int mlen, int sublen);


template<typename dtype, typename itype>
int maxpearson_partialauto(const stridedbuf<dtype>& ts, int minlag, int mlen, int sublen){
   // allocate query buffers, basically covariance and normalized query
   int tilesperdim = mlen/tsz + (mlen%tsz ? 1 : 0); // find max tiles in 1 direction
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }
   int n = ts.len;
   stridedbuf<dtype> mu();
   stridedbuf<dtype> invn();
   stridedbuf<dtype> df();
   stridedbuf<dtype> dx();
   stridedbuf<dtype> mp();
   stridedbuf<dtype>cov(tsz,n,n);
   stridedbuf<dtype>q(tilesperdim,paddedlen(sublen,prefalign));
   stridedbuf<itype> mpi();
   if(!(q.isvalid() && cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dx.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");   
   } 
   xmean_windowed(ts(0),mu(0),ts.len,sublen);
   xsinv(ts(0),mu(0),invn(0),ts.len,sublen);   
   init_dfdx(ts(0), mu(0), df(0), dx(0),sublen,ts.len);
   std::fill(mp(0),mp(0)+mp.len,-1.0);
   std::fill(mpi(0),mpi(0)+mpi.len,-1);          
   const int step = tsz/ksz;
   #pragma omp parallel for
   for(int i = 0; i < mlen; i+= tsz){
      center_query(ts(i), mu(i), q(i), sublen);
   }
   for(int i = 0; i < tilesperdim; i++){
      #pragma omp parallel for
      for(int j = 0; j < tilesperdim-i; j++){
         batchcov(ts(j),cov(j),q(j),mu(i+j),tsz,sublen);
         int cindoffset = i*tsz+minlag;
         int qbaseind = j*tsz; 
         double qinvn = invn(j)[0];
         int prefixlen = (i+j)*tsz;
         for(int k = 0; k < step; k++){
            auto r = rescaled_max_reduct(cov(j+k),invn(j+k),mp(j+k),mpi(j+k),qinvn,mp.dat[0],qbaseind,cindoffset);
            if(r.val > mp.dat[0]){
               mp.dat[0] = r.val; 
               mpi.dat[0] = r.ind;
            }
            if(prefixlen+(k+1)*ksz <= mlen){
               pauto_pearson_kern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz); 
            }
            else{
               pauto_pearson_edgekern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz,mlen-(offset+l)*ksz);
            }
         }
         for(int k = 1; k < step; k++){  // do we need an edge check here? 
            for(int l = 0; l < step; l++){
               int offset = j*step+k;   // double check this
               if(prefixlen+(k+1)*ksz <= mlen){
                  pauto_pearson_kern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz); 
               }
               else{
                  pauto_pearson_edgekern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz,mlen-(offset+l)*ksz);
               }
            }
         }
      }
   }
}


//template<>
//int maxpearson_partialauto<double,long>(const stridedbuf<double>& ts, int minlag, int mlen, int sublen);
//template<float,int>
//int maxpearson_partialauto(const stridedbuf<float>& ts, int minlag, int mlen, int sublen);
//template<double,int>
//int maxpearson_partialauto(const stridedbuf<double>& ts, int minlag, int mlen, int sublen);
