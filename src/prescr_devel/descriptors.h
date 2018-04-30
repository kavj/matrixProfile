//#include<cstdio>
//#include<unistd.h>
#include "alloc.h"
#include "mkl.h"


// rearrange attributes
static inline int   __attribute__((always_inline)) paddedlen(int buflen, int unitsz, int alignmt){ return  (buflen*unitsz) + ((buflen*unitsz)%alignmt ? alignmt - (buflen*unitsz)%alignmt : 0);}

template<typename dtype> struct buf_strided{
   
   buf_strided(){} 
   inline __attribute__((always_inline)) buf_strided(int blklen, int blkcount, int alignmt) : bcount(blkcount), blen(blklen), bstride(paddedlen(blklen,sizeof(dtype),alignmt)), ownsmem(true){
       dat = reinterpret_cast<dtype>(init_buffer(blklen,blkcount,alignmt));
   }
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
struct buf_indexed{
   inline __attribute__((always_inline)) buf_indexed() : dat(nullptr), len(0), blen(0), ownsmem(false){}
   buf_indexed(dtype* dat, int len, int blocklen) : dat(dat), len(len), blen(blocklen), ownsmem(false) {}
   ~buf_indexed(){
      if(ownsmem){
         free(dat);
      }
   }
   inline dtype* __attribute__((always_inline)) operator()(int i)          { return (i*blen < len) ? dat + i*blen : nullptr;}
   inline dtype* __attribute__((always_inline)) bystride(int i,int stride) { return (i*stride < len) ? dat + i*stride : nullptr;}
   inline int    __attribute__((always_inline)) offset(int i)              { return (i*blen < len) ? i*blen : -1;}
   dtype* dat;
   int len;
   int blen;
   private:
   bool ownsmem;
};


template<typename dtype>
struct acorr_desc{ 
  inline __attribute__((always_inline)) acorr_desc(int qlen, int qbufct, int qstatslen, int qstatsbufct, int databuflen, int databufct, int alignmt) {}//: //ts(nullptr), cov(nullptr), xcorr(nullptr), xind(nullptr), mu(nullptr), invn(nullptr), df(nullptr), dx(nullptr), q(nullptr), 
                                                                                                                    //qcov(nullptr), qcorr(nullptr), covbufs(nullptr), qmatch(nullptr), covtsks(nullptr){}
   struct buf_indexed<dtype> ts;
   struct buf_indexed<dtype> cov;
   struct buf_indexed<dtype> xcorr;
   struct buf_indexed<int> xind;
   struct buf_indexed<dtype> mu;
   struct buf_indexed<dtype> invn;
   struct buf_indexed<dtype> df;
   struct buf_indexed<dtype> dx;

   struct buf_strided<dtype> q;
   struct buf_strided<dtype> qcov;
   struct buf_strided<dtype> qcorr;
   struct buf_strided<dtype> covbufs;
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


   inline int __attribute__((always_inline)) isinitialized(){return (q.dat != nullptr) && (qcov.dat != nullptr) && (qcorr.dat != nullptr) && (covbufs.dat != nullptr) && (qmatch.dat != nullptr) && (covtsks != nullptr) && (ts.dat != nullptr) 
                                                                 && (cov.dat != nullptr) && (xcorr.dat != nullptr) && (xind.dat != nullptr) && (mu.dat != nullptr) && (invn.dat != nullptr) && (df.dat != nullptr) && (dx.dat != nullptr);} 

   inline int __attribute__((always_inline)) required_passes(){return (q.bcount/querycount) + ((q.bcount%querycount) ? 1 : 0);} 
};


struct query_stat{
   query_stat(double qcov, int qind) : qcov(qcov), qind(qind){}
   double qcov;
   int qind;
};



