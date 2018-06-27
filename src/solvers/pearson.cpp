#include <cstdio>
#include<algorithm>
#include "../utils/descriptors.h"
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "tiled_pearson.h"
#include "pearson.h"



static void init_dfdx(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ df, dtype* __restrict__ dx, int w, int n){
   df[0] = 0; 
   dx[0] = 0;
   for(int i = 0; i < n-w; i++){
      df[i+1] = (ts[i+w] - mu[i+1]) + (ts[i] - mu[i]);
      dx[i+1] = (ts[i+w] - ts[i])/2.0;
   }
}


void pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){

   int mlen = ts.len - sublen + 1;
   const int tlen = 16384;
   int tail = (mlen - minlag)%tlen;
   int tilesperdim = (mlen - minlag - tail)/tlen + (tail ? 1 : 0);
 
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim,sublen);

   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   

   init_dfdx(ts(0), mu(0), df(0), dg(0), sublen, ts.len);
   std::fill(mp(0), mp(mlen), -1.0);
   std::fill(mpi(0), mpi(mlen), -1); 

   for(int i = 0; i < tilesperdim; i++){
      center_query_ref(ts(i * tlen), mu(i * tlen), q(i), sublen);
   }

   int aligned = std::max(0, tilesperdim - 1 - (tail > 0 ? 1 : 0));

   for(int diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for
      for(int ofst = 0; ofst < tilesperdim-diag; ofst++){
         int rofst = ofst * tlen;
         int cofst = (diag + ofst) * tlen;
         int initofst = cofst + minlag;
         if(ofst < tilesperdim - diag){
            batchcov_ref(ts(initofst), mu(initofst), q(ofst), cov(cofst), tlen, sublen);
            pauto_pearson_inner(cov(rofst), mp(rofst), mpi(rofst), df(rofst), dg(rofst), invn(rofst), tlen, rofst, cofst);
         }
         else{
            int mxofst = mlen - (diag + ofst)*tlen;
            int width = std::min(tlen, mxofst);
            batchcov_ref(ts(initofst), mu(initofst), q(ofst), cov(cofst), width, sublen);
            pauto_pearson_edge(cov(rofst), mp(rofst), mpi(rofst), df(rofst), dg(rofst), invn(rofst), tlen, rofst, cofst, mxofst);
         }
      }
   }
}
