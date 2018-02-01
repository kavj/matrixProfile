#include<cstdio>
#include "../arch/v256.h"
#include "../utils/xprec_math.hpp"

/*
template<typename vt>
static inline void xcorr_reduce_horizontal(double* __restrict__ buf, double* __restrict__ prof, long* __restrict__ profi, int diag, int offset){
   int stride = 4;
   vt ca = aload(buf,offset);
   vt cb = aload(buf,offset+stride);
   vt cc = aload(buf,offset+2*stride);
   vt cd = aload(buf,offset+3*stride);
   vt ce = aload(buf,offset+4*stride);
   vt cf = aload(buf,offset+5*stride);
//   for(int i = 0; i < tilesz; i++){
      
//   }
}*/

static inline void xcorr_reduce_vertical(double* __restrict__ buf, double* __restrict__ prof, long* __restrict__ profi, int diag, int offset){
   

}


static inline void xcorr_2(vt &ca, vt &cb, vt &cc, vt &cd, vt &ce, vt &cf, vt &cg, vt &ch, double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double*  __restrict__ output, int diag, int k) __attribute__ ((always_inline));

template<typename vt>
static inline void xcorr_2(vt &ca, vt &cb, vt &cc, vt &cd, vt &ce, vt &cf, vt &cg, vt &ch, double *dx, double *xf, double *s, double *buf, int diag, int offset){
   int p = diag+offset;
   int stride = 4;
   vt fbc = bcast(dx,offset);
   
   ca = vfma(fbc,aload(xf,p),ca);
   cb = vfma(fbc,aload(xf,p+stride),cb);
   cc = vfma(fbc,aload(xf,p+2*stride),cc);
   cd = vfma(fbc,aload(xf,p+3*stride),cd);
   ce = vfma(fbc,aload(xf,p+4*stride),ce);
   cf = vfma(fbc,aload(xf,p+5*stride),cf);
   cg = vfma(fbc,aload(xf,p+6*stride),cg);
   ch = vfma(fbc,aload(xf,p+7*stride),ch);

   fbc = bcast(xf,offset); 
   
   ca = vfma(fbc,aload(dx,p),ca);
   cb = vfma(fbc,aload(dx,p+stride),cb);
   cc = vfma(fbc,aload(dx,p+2*stride),cc);
   cd = vfma(fbc,aload(dx,p+3*stride),cd);
   ce = vfma(fbc,aload(dx,p+4*stride),ce);
   cf = vfma(fbc,aload(dx,p+5*stride),cf);
   cg = vfma(fbc,aload(dx,p+6*stride),cg);
   ch = vfma(fbc,aload(dx,p+7*stride),ch);

   vt da, db, dc, dd, de, df, dg, dh;
   fbc = bcast(s,offset);
   da = ca * fbc;
   db = cb * fbc;
   dc = cc * fbc;
   dd = cd * fbc;
   de = ce * fbc;
   df = cf * fbc;
   dg = cg * fbc;
   dh = ch * fbc;

   da *= aload(s,p);
   db *= aload(s,p+stride);
   dc *= aload(s,p+2*stride);
   dd *= aload(s,p+3*stride);
   de *= aload(s,p+4*stride);
   df *= aload(s,p+5*stride);
   dg *= aload(s,p+6*stride);
   dh *= aload(s,p+7*stride);
  
   int k = 4*p;

   strided_store_x8<vt,double,4>(da,db,dc,dd,de,df,dg,dh,buf,4*p);
 
}





static inline void xcorr_reduce(double* __restrict__ buf, double* __restrict__ output, long*  __restrict__ outputi, int diag, int k) __attribute__ ((always_inline));

static inline void xcorr_reduce(double* __restrict__ buf, double* __restrict__ output, long*  __restrict__ outputi, int diag, int k){
 // for(int i = 0; i < 32 

}

template<typename vt>
static inline void xcorr_kern(double* __restrict__ cx, const double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double* __restrict__ buf, int diag, int offset) __attribute__ ((always_inline));

template<typename vt>
static inline void xcorr_kern(double* __restrict__ cx, double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double* __restrict__ buf, int diag, int offset){
   for(int subdiag = diag; subdiag < diag + 2048; subdiag += 32){
      vt ca, cb, cc, cd, ce, cf, cg, ch;
     // strided_load_x8<vt,double,4>(ca,cb,cc,cd,ce,cf,cg,ch,cx,subdiag);
      
      for(int suboff = offset; suboff < offset+2048; suboff+=4){
         xcorr_2 <__m256d>(ca,cb,cc,cd,ce,cf,cg,ch, dx, xf, s, buf, subdiag, suboff);
      }
     // strided_store_x8<vt,double,4>(ca,cb,cc,cd,ce,cf,cg,ch,cx,subdiag);
   }
}

//#define tilesz 2048
template<typename dtype, typename itype>
void  accumTest4_7_10(dtype* cx,dtype* dx, dtype* xf, dtype*s, dtype* buf, dtype* output, itype* outputi,int n, int m){
   unsigned long t= 0;
   const int tilesz = 2048;
   const int stride = 4;
   #pragma omp parallel for  
   for(int diag = 0; diag < n-16*m; diag += 2048){
      for(int offset = 0; offset < n-diag-16*m; offset+=2048){
         t ++;
         xcorr_kern(cx,dx,xf,s,buf,diag,offset);
      }
   }   
   printf("iterations: %lu \n",t);
}


