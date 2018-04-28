//#include<cstdio>
//#include<unistd.h>
#include "alloc.h"
#include "mkl.h"

// This should probably be temporary. Since I end up having to distribute stuff frequently, it would be a good idea to generalize the striding somewhat for 1 and 2D strides

// rearrange attributes
static inline int   __attribute__((always_inline)) paddedlen(int buflen, int unitsz, int alignmt){ return  (buflen*unitsz) + ((buflen*unitsz)%alignmt ? alignmt - (buflen*unitsz)%alignmt : 0);}

template<typename dtype> struct buf_strided{
   inline __attribute__((always_inline)) buf_strided(){}
   inline __attribute__((always_inline)) buf_strided(int blklen, int blkcount, int alignmt) : bcount(blkcount), blen(blklen), bstride(paddedlen(blklen,sizeof(dtype),alignmt)), dat(nullptr){}
 
   inline __attribute__((always_inline)) ~buf_strided(){
      if(ownsmem){
         free(dat);
      }
   }
   inline dtype*  __attribute__((always_inline)) operator()(int i){ return i < bcount ? dat + i*bstride : nullptr;}
   dtype* dat;
   int blen;   
   int bcount;
   int bstride; 
   private:
   bool ownsmem;
};



// indexed and continuous subdivision
template<typename dtype>
struct buf_overlapped{
   inline __attribute__((always_inline)) buf_overlapped() : dat(nullptr), len(0), blen(0), bstride(0){}
   buf_overlapped(dtype* dat, int len, int blocklen) : dat(dat), len(len), blen(blocklen){}
   ~buf_overlapped(){
      if(ownsmem){
         dtype* dat;
      }
   }
   dtype* dat;
   int len;
   int blen;
   int bstride;
   private:
   bool ownsmem;
};




// I'm a little conflicted on this part. The custom accessors are probably sufficient for now, but this might not be scalable compared to a more sophisticated buffer interface
// There's also the issue of memory lifetime. In cases where memory is externally allocated and deallocated, we need to avoid deallocation by the destructor
//

template<typename dtype>
struct acorr_desc{ 
  inline __attribute__((always_inline)) acorr_desc(int qlen, int qbufct, int qstatslen, int qstatsbufct, int databuflen, int databufct,  int alignmt) {}//: //ts(nullptr), cov(nullptr), xcorr(nullptr), xind(nullptr), mu(nullptr), invn(nullptr), df(nullptr), dx(nullptr), q(nullptr), 
                                                                                                                    //qcov(nullptr), qcorr(nullptr), covbufs(nullptr), qmatch(nullptr), covtsks(nullptr){}
   dtype* ts;
   dtype* cov;
   dtype* xcorr;
   dtype* xind;
   dtype* mu;
   dtype* invn;
   dtype* df;
   dtype* dx;

   struct buf_strided<double> q;
   struct buf_strided<double> qcov;
   struct buf_strided<double> qcorr;
   struct buf_strided<double> covbufs;
   struct buf_strided<int> qmatch; 
   VSLCorrTaskPtr* covtsks; // descriptor for MKL
   int querycount;   // total queries required
   int qbasestride;  // indicates distance between queries with respect to a time series. This is used in indexing queries
   int taillen;
   int tailqcount;   // this is just however many queries  
   int xcorrlen;
   int sublen;
   int bcount;
   int blen;
   int bstride;    // normal block stride

   // This allows for slightly more arbitrary indexing with a bounds check on the first point.
   // I am still somewhat undecided. The code itself might be better if I simply used structs to impl

   inline dtype*  __attribute__((always_inline)) sts   (int i, int stride) { return (i*stride) < xcorrlen+sublen-1 ? ts   + i*stride : nullptr;} 
   inline dtype*  __attribute__((always_inline)) scov  (int i, int stride) { return (i*stride) < xcorrlen ? cov  + i*stride : nullptr;}
   inline dtype*  __attribute__((always_inline)) sxcorr(int i, int stride) { return (i*stride) < xcorrlen ? xcorr+ i*stride : nullptr;}
   inline dtype*  __attribute__((always_inline)) sxind (int i, int stride) { return (i*stride) < xcorrlen ? xind + i*stride : nullptr;}
   inline dtype*  __attribute__((always_inline)) smu   (int i, int stride) { return (i*stride) < xcorrlen ? mu   + i*stride : nullptr;}
   inline dtype*  __attribute__((always_inline)) sinvn (int i, int stride) { return (i*stride) < xcorrlen ? invn + i*stride : nullptr;}
   inline dtype*  __attribute__((always_inline)) sdf   (int i, int stride) { return (i*stride) < xcorrlen ? df   + i*stride : nullptr;}
   inline dtype*  __attribute__((always_inline)) sdx   (int i, int stride) { return (i*stride) < xcorrlen ? dx   + i*stride : nullptr;}


   
   inline int __attribute__((always_inline)) isinitialized(){return (q.dat != nullptr) && (qcov.dat != nullptr) && (qcorr.dat != nullptr) && (covbufs.dat != nullptr) && (qmatch.dat != nullptr);} 
   inline int __attribute__((always_inline)) required_passes(){return (q.bcount/querycount) + ((q.bcount%querycount) ? 1 : 0);} 

};


struct query_stat{
   query_stat(double qcov, int qind) : qcov(qcov), qind(qind){}
   double qcov;
   int qind;
};



