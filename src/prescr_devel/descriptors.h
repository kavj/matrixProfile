//#include<cstdio>
//#include<unistd.h>
#include "alloc.h"
#include "mkl.h"

// This should probably be temporary. Since I end up having to distribute stuff frequently, it would be a good idea to generalize the striding somewhat for 1 and 2D strides

// rearrange attributes
static inline int   __attribute__((always_inline)) paddedlen(int buflen, int unitsz, int alignmt){ return  (buflen*unitsz) + ((buflen*unitsz)%alignmt ? alignmt - (buflen*unitsz)%alignmt : 0);}

template<typename dtype> struct buf_strided{
   buf_strided();
   buf_strided(int blklen, int blkcount, int alignmt) : bcount(blkcount), blen(blklen), bstride(paddedlen(blklen,sizeof(dtype),alignmt)), dat(nullptr){}
 
   ~buf_strided(){
      free(dat);
   }
   inline dtype*  __attribute__((always_inline)) operator()(int i){ return i < bcount ? dat + i*bstride : nullptr;}
   dtype* dat;
   int blen;   
   int bcount;
   int bstride; 
};



// indexed and continuous subdivision
template<typename dtype>
struct buf_indexed{
   buf_indexed() : dat(nullptr), len(0), blen(0), bstride(0){}
   buf_indexed(dtype* dat, int len, int blocklen) : dat(dat), len(len), blen(blocklen){}
   dtype* dat;
   int len;
   int blen;
   int bstride;
};


//template<typename dtype>
template<typename dtype>
struct corr_auxbuf{ 
   corr_auxbuf(int qlen, int qbufct, int qstatslen, int qstatsbufct, int databuflen, int databufct,  int alignmt) {}//: //ts(nullptr), cov(nullptr), xcorr(nullptr), xind(nullptr), mu(nullptr), invn(nullptr), df(nullptr), dx(nullptr), q(nullptr), 
                                                                                                                    //qcov(nullptr), qcorr(nullptr), covbufs(nullptr), qmatch(nullptr), covtsks(nullptr){}
   dtype* ts;
   dtype* cov;
   dtype* xcorr;
   dtype* xind;
   dtype* mu;
   dtype* invn;
   dtype* df;
   dtype* dx;;

   struct buf_strided<double> q;
   struct buf_strided<double> qcov;
   struct buf_strided<double> qcorr;
   struct buf_strided<double> covbufs;
   struct buf_strided<int> qmatch; 
   VSLCorrTaskPtr* covtsks; // descriptor for MKL
   int querycount;   // total queries required
   int qbasestride;  // indicates distance between queries with respect to a time series. This is used in indexing queries
   int tailcount;   // this is just however many queries  
                   // this is REM(length(time_series)/length(buffer))
                   // for covbufs, need way to identify the fringe component. Perhaps  

   inline dtype*  __attribute__((always_inline)) gts   (int i) { return i < bcount ? ts    + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) gcov  (int i) { return i < bcount ? cov   + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) gxcorr(int i) { return i < bcount ? xcorr + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) gxind (int i) { return i < bcount ? xind  + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) gmu   (int i) { return i < bcount ? mu    + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) ginvn (int i) { return i < bcount ? invn  + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) gdf   (int i) { return i < bcount ? df    + i*bstride : nullptr;}
   inline dtype*  __attribute__((always_inline)) gdx   (int i) { return i < bcount ? dx    + i*bstride : nullptr;}
   
   inline int __attribute__((always_inline)) isinitialized(){return (q.dat != nullptr) && (qcov.dat != nullptr) && (qcorr.dat != nullptr) && (covbufs.dat != nullptr) && (qmatch.dat != nullptr);} 
   inline int __attribute__((always_inline)) required_passes(){return (q.bcount/querycount) + ((q.bcount%querycount) ? 1 : 0);} 

  //  Ultimately need to be able to do something like reduce remaining from length(data) to 0 somewhere. 
  //  inline double  __attribute__((always_inline)) operator()(int i){ return i < bcount ? i*bstride : -1; }
  //  inline int   __attribute__((always_inline)) isinitialized(){return (ts != nullptr) && (cov != nullptr) && (xcorr != nullptr) && (xind != nullptr) && (mu != nullptr)  && (invn != nullptr) && (df != nullptr) && (dx != nullptr); } 
 
};


struct query_stat{
   double qcov;
   int qind;
};



