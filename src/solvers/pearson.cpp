#include <cstdio>
#include <algorithm>
#include "../utils/xprec.h"
#include "../utils/cov.h"
#include "../utils/alloc.h"
#include "tiled_pearson.h"
#include "pearson.h"


static void init_dfdg(const dtype* __restrict__ ts, const dtype* __restrict__ mu, dtype* __restrict__ df, dtype* __restrict__ dg, int len, int sublen){
   df[0] = 0; 
   dg[0] = 0;
   for(int i = 0; i < len - sublen; i++){
      df[i + 1] = (ts[i + sublen] - mu[i + 1]) + (ts[i] - mu[i]);
      dg[i + 1] = (ts[i + sublen] - ts[i])/2.0;
   }
}


void pearson_pauto_reduc(dsbuf& ts, dsbuf& mp, lsbuf& mpi, int minlag, int sublen){

   if(!(ts.valid() && mp.valid() && mpi.valid())){
      printf("bad inputs\n");
      exit(1);
   }

   int mlen = ts.len - sublen + 1;
   const int tlen = 16384;
   int tail = (mlen - minlag) % tlen;
   int tilesperdim = (mlen - minlag - tail)/tlen + (tail ? 1 : 0);
 
   dsbuf mu(mlen); dsbuf invn(mlen); dsbuf df(mlen);  
   dsbuf dg(mlen); dsbuf cov(mlen);  mdsbuf q(tilesperdim, sublen);

   xmean_windowed(ts(0), mu(0), ts.len, sublen);
   xsInv(ts(0), mu(0), invn(0), ts.len, sublen);   
   init_dfdg(ts(0), mu(0), df(0), dg(0), ts.len, sublen);
   std::fill(mp(0), mp(0) + mlen, -1.0);
   std::fill(mpi(0), mpi(0) + mlen, -1); 

   if(!(mu.valid() && df.valid() && dg.valid() && invn.valid())){
      printf("error allocating memory\n");
      exit(1);
   }
   for(int i = 0; i < tilesperdim; i++){
      center_query_ref(ts(i * tlen), mu(i * tlen), q(i), sublen);
   }
   for(int diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for
      for(int ofst = 0; ofst < tilesperdim - diag; ofst++){
         int rofst = ofst * tlen;
         int cofst = diag * tlen;
         int initofst = (diag + ofst) * tlen + minlag;
         if(((diag + ofst + 2) * tlen + minlag) < mlen){
            batchcov_ref(ts(initofst), mu(initofst), q(ofst), cov(rofst), tlen, sublen);
            pauto_pearson_inner(cov(rofst), mp(rofst), mpi(rofst), df(rofst), dg(rofst), invn(rofst), tlen, rofst, cofst);
         }
         else{
            int mxofst = mlen - ((diag + ofst) * tlen + minlag);
            int width = std::min(tlen, mxofst);
            batchcov_ref(ts(initofst), mu(initofst), q(ofst), cov(rofst), width, sublen);
            pauto_pearson_edge(cov(rofst), mp(rofst), mpi(rofst), df(rofst), dg(rofst), invn(rofst), tlen, rofst, cofst, mxofst);
         }
      }
   }
}
