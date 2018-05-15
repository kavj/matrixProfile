#include<cstdio>
               pauto_pearson_edgekern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz,mlen-(offset+l)*ksz);
               pauto_pearson_edgekern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz,mlen-(offset+l)*ksz);
#include<cmath>
#include<algorithm>
#include "descriptors.h"
#include "exact_Pearson.h"




template<typename dtype, typename itype>
int maxpearson_partialauto(const stridedbuf<dtype>& ts, int minlag, int mlen, int sublen){
   // allocate query buffers, basically covariance and normalized query
   int tilesperdim = mlen/tsz + (mlen%tsz ? 1 : 0); // find max tiles in 1 direction
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }
   //int n = ts.len;
   stridedbuf<dtype> mu(mlen,tsz,tilesperdim);
   stridedbuf<dtype> invn(mlen,tsz,tilesperdim);
   stridedbuf<dtype> df(mlen,tsz,tilesperdim);
   stridedbuf<dtype> dx(mlen,tsz,tilesperdim);
   stridedbuf<dtype> mp(mlen,tsz,tilesperdim);
   stridedbuf<dtype>cov(mlen,tsz,tilesperdim);
   stridedbuf<dtype>q(tilesperdim,paddedlen(sublen,prefalign));
   stridedbuf<itype> mpi(mlen,tsz,tilesperdim);
   if(!(q.isvalid() && cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dx.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");   
   } 
   xmean_windowed(ts(0),mu(0),ts.len,sublen);
   xsinv(ts(0),mu(0),invn(0),ts.len,sublen);   
   init_dfdx(ts(0), mu(0), df(0), dx(0),sublen,ts.len);
   std::fill(mp(0),mp(0)+mp.len,-1.0);
   std::fill(mpi(0),mpi(0)+mpi.len,-1);          
   //const int step = tsz/ksz;
   
   #pragma omp parallel for
   for(int i = 0; i < mlen; i+= tsz){
      center_query(ts(i), mu(i), q(i), sublen);
   }
   for(int i = 0; i < tilesperdim-minlag; i++){
      #pragma omp parallel for
      for(int j = 0; j < mlen-i; j++){
         int len = std::min(tsz,mlen-i-minlag-2*tsz);
         batchcov(ts(j),cov(j),q(j),mu(i+j),len,sublen);
         int cindoffset = i*tsz+minlag;
         int qbaseind = j*tsz; 
         double qinvn = invn(j)[0];
         int mxcolstep = std::min(step,(mlen-minlag)/ksz);
         for(int k = 0; k < step; k++){
            for(int l = 0; l < step; l++){
               if((i+j)*tsz + (k+l+1)*ksz < mlen){
                  pauto_pearson_kern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz); 
               }
               else{
                  pauto_pearson_edgekern(cov(offset),df(offset),dx(offset),invn(offset),mp(offset),mpi(offset),j*tsz,(offset+l)*ksz,mlen-(offset+l)*ksz);
               }
            }
            // when did I write this? This should be strided 
            auto r = rescaled_max_reduct(cov(j+k),invn(j+k),mp(j+k),mpi(j+k),qinvn,mp.dat[0],qbaseind,cindoffset);
            if(r.val > mp.dat[0]){
               mp.dat[0] = r.val; 
               mpi.dat[0] = r.ind;
            }
         }
      }
   }
}


