#include <cstdio>
#include<algorithm>
#include "descriptors.h"
#include "../utils/xprec_math.h"
#include "../utils/cov.h"
#include "tiled_Pearson.h"
#include "../utils/primitive_print_funcs.h"
#define prefalign 64


// redesigning this part
// It will handle basic preprocessing only in cases where some things are not externally supplied, particularly since we 
// can get the core of the algorithms down to sufficiently simple single functions

typedef double dtype;
typedef long long itype;

void maxpearson_partialauto(stridedbuf<dtype>& ts, stridedbuf<dtype>& mp, stridedbuf<itype>& mpi, int minlag, int sublen){
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }

   int mlen = ts.len - sublen + 1;
   const int tlen = 16384;
   int tail = (mlen - minlag)%tlen;
   int tilesperdim = (mlen - minlag - tail)/tlen + (tail ? 1 : 0);
 
   stridedbuf<dtype>mu(mlen); stridedbuf<dtype>invn(mlen); stridedbuf<dtype>df(mlen);  
   stridedbuf<dtype>dg(mlen); stridedbuf<dtype>cov(mlen);  
   multibuf<dtype> q(tilesperdim,sublen);

   if(!(cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dg.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");
      return;
   } 

   ts.setstride(tlen); cov.setstride(tlen); mu.setstride(tlen); df.setstride(tlen); 
   dg.setstride(tlen); invn.setstride(tlen); mp.setstride(tlen); mpi.setstride(tlen);

   xmean_windowed(ts(0),mu(0),ts.len,sublen);
   xsInv(ts(0),mu(0),invn(0),ts.len,sublen);   
   init_dfdx(ts(0), mu(0), df(0), dg(0),sublen,ts.len);
   std::fill(mp(0),mp(0)+mlen,-1.0);
   std::fill(mpi(0),mpi(0)+mlen,-1); 

   for(int i = 0; i < tilesperdim; i++){
      center_query_ref(ts(i),mu(i),q(i),sublen);
   }

   for(int diag = 0; diag < tilesperdim; diag++){
      #pragma omp parallel for
      for(int offst = diag; offst < tilesperdim-diag; offst++){
         printf("diag: %d, offset: %d diag+offset+2*tlen == %d\n",diag*tlen,offst*tlen,(diag+offst+2)*tlen);
         if(minlag+(diag+offst+2)*tlen <= mlen){
            printf("on inner\n");
            batchcov_ref(ts(diag+offst)+minlag,cov(offst),q(offst),mu(offst),tlen,sublen);
            pauto_pearson_basic_inner(cov(offst),mp(offst),mpi(offst),df(offst),dg(offst),invn(offst),tlen,offst*tlen,diag*tlen+minlag);
         }
         else{
            printf("on edge\n");
           // batchcov_ref(ts(offst)+minlag,cov(offst),q(offst),mu(offst),std::min(tlen,mlen-minlag-(diag+offst)*tlen),sublen);
           // pauto_pearson_edge(cov(offst),mp(offst),mpi(offst),df(offst),dg(offst),invn(offst),tlen,offst*tlen,diag*tlen+minlag,mlen-minlag-(diag+offst)*tlen);
         }
      }
   }
}

