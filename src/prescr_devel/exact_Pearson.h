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
typedef int itype;

void maxpearson_partialauto(stridedbuf<dtype>& ts, stridedbuf<dtype>& mp, stridedbuf<itype>& mpi, int minlag, int sublen){
   if(!ts.isvalid()){
      printf("invalid time series\n");
   }

   int mlen = ts.len - sublen + 1;

   stridedbuf<dtype>mu(mlen); stridedbuf<dtype>invn(mlen); stridedbuf<dtype>df(mlen);  
   stridedbuf<dtype>dg(mlen); stridedbuf<dtype>cov(mlen);

   if(!(cov.isvalid() && mu.isvalid() && invn.isvalid() && df.isvalid()  && dg.isvalid() && mp.isvalid() && mpi.isvalid())){
      printf("could not assign objects\n");
      return;
   } 
  

   xmean_windowed(ts(0),mu(0),ts.len,sublen);
   xsInv(ts(0),mu(0),invn(0),ts.len,sublen);   
   init_dfdx(ts(0), mu(0), df(0), dg(0),sublen,ts.len);
   std::fill(mp(0),mp(0)+mlen,-1.0);
   std::fill(mpi(0),mpi(0)+mlen,-1); 
   pauto_pearson(cov.dat,mp.dat,mpi.dat,ts.dat,mu.dat,df.dat,dg.dat,invn.dat,mlen,sublen,minlag); 

}


