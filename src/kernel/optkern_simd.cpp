#include<cstdio>
//#include "v256.h"
#include "../arch/avx2arith_2.hpp"
//#include "../utils/xprec_math.hpp"

// temporary 
typedef  __m256d vt;
using namespace vmth;

//template<typename vt>
static inline void xcorr_2(vt &ca, vt &cb, vt &cc, vt &cd, vt &ce, vt &cf, vt &cg, vt &ch, double* __restrict__ dx, double* __restrict__ xf, int diag, int offset) __attribute__ ((always_inline));

//template<typename vt>

static inline void xcorr_2(vt &ca, vt &cb, vt &cc, vt &cd, vt &ce, vt &cf, vt &cg, vt &ch, double *dx, double *xf, int diag, int offset){
   int p = diag+offset;
   int stride = 4;

   vt fbc = bcast(dx,offset);
   
   ca = fma(fbc,aload(xf,p),ca);
   cb = fma(fbc,aload(xf,p+stride),cb);
   cc = fma(fbc,aload(xf,p+2*stride),cc);
   cd = fma(fbc,aload(xf,p+3*stride),cd);
   ce = fma(fbc,aload(xf,p+4*stride),ce);
   cf = fma(fbc,aload(xf,p+5*stride),cf);
   cg = fma(fbc,aload(xf,p+6*stride),cg);
   ch = fma(fbc,aload(xf,p+7*stride),ch);

   fbc = bcast(xf,offset);   

   ca = fma(fbc,aload(dx,p),ca);
   cb = fma(fbc,aload(dx,p+stride),cb);
   cc = fma(fbc,aload(dx,p+2*stride),cc);
   cd = fma(fbc,aload(dx,p+3*stride),cd);
   ce = fma(fbc,aload(dx,p+4*stride),ce);
   cf = fma(fbc,aload(dx,p+5*stride),cf);
   cg = fma(fbc,aload(dx,p+6*stride),cg);
   ch = fma(fbc,aload(dx,p+7*stride),ch);


}
typedef __m256i vi;
static inline void xcorr_reduce(vt ca, vt cb, vt cc, vt cd, vt ce, vt cf, vt cg, vt ch, double *s, double *output, long *outputi, int diag, int offset){
   const int stride = 4;
   const int p = diag+offset;
   const int q = 4*p; 
   const int r = 4*offset;
   
   ca *= aload(s,p);
   cb *= aload(s,p+stride);
   cc *= aload(s,p+2*stride);
   cd *= aload(s,p+3*stride);
   ce *= aload(s,p+4*stride);
   cf *= aload(s,p+5*stride);
   cg *= aload(s,p+6*stride);
   ch *= aload(s,p+7*stride);
  
   vt fbc = bcast(s,offset);
   ca = ca * fbc;
   cb = cb * fbc;
   cc = cc * fbc;
   cd = cd * fbc;
   ce = ce * fbc;
   cf = cf * fbc;
   cg = cg * fbc;
   ch = ch * fbc;

   vt da, db, dc, dd;
   vi ga, gb, gc, gd;

   const vi four = bcast(4);
   const vi st   = set(0,1,2,3);
   vi st2 = st + bcast(offset);
   vi st3 = st2 + four;
   
   da = cb < ca;
   db = cc < cd;
   dc = ce < cf;
   dd = cg < ch;

   ga = blend(st2,st3,da);
   st2 = st3;
   st3 = st3 + four;
   gb = blend(st2,st3,db);
   st2 = st3;
   st3 = st3 + four;
   gc = blend(st2,st3,dc);
   st2 = st3;
   st3 = st3 + four;
   gd = blend(st2,st3,dd);

   ca = max(ca,cb);
   cc = max(cc,cd);
   cd = max(ce,cf);
   cg = max(cg,ch);

   cb = cc < ca;
   ce = cg < cd;

   ga = blend(ga,gb,cb);
   gc = blend(gc,gd,ce);
   
    
   

/*   ca = max(da,db);
   cc = max(dc,dd);
   ce = max(de,df);
   cg = max(dg,dh);

   ca = max(ca,cc);
   ce = max(ce,cg);
   ca = max(ca,ce);
   astore(ca,output,r);

   da = max(da,aload(output,q));
   db = max(db,aload(output,q+stride));
   dc = max(dc,aload(output,q+2*stride));
   dd = max(dd,aload(output,q+3*stride));
   de = max(de,aload(output,q+4*stride));
   df = max(da,aload(output,q+5*stride));
   dg = max(db,aload(output,q+6*stride));
   dh = max(dc,aload(output,q+7*stride));
 
   strided_store_x8<vt,4>(da,db,dc,dd,de,df,dg,dh,output,q 
*/
}



static inline void xcorr_kern(double* __restrict__ cx, const double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double* __restrict__ output, long* __restrict__ outputi, int diag, int offset) __attribute__ ((always_inline));

static inline void xcorr_kern(double* __restrict__ cx, double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double* __restrict__ output, long* __restrict__ outputi, int diag, int offset){
   for(int subdiag = diag; subdiag < diag + 2048; subdiag += 32){
      vt ca, cb, cc, cd, ce, cf, cg, ch;
      strided_load_x8<vt,4>(ca,cb,cc,cd,ce,cf,cg,ch,cx,subdiag);
      for(int suboff = offset; suboff < offset+2048; suboff+=4){
         xcorr_2(ca,cb,cc,cd,ce,cf,cg,ch,dx,xf,subdiag,suboff);
         xcorr_reduce(ca,cb,cc,cd,ce,cf,cg,ch,s,output,outputi,subdiag,suboff);
      }
      strided_store_x8<vt,4>(ca,cb,cc,cd,ce,cf,cg,ch,cx,subdiag);
   }
}

//#define tilesz 2048
//template<typename dtype, typename itype>
void  accumTest4_7_10(double* cx,double* dx, double* xf, double*s, double* buf, double* output, long* outputi,int n, int m){
   unsigned long t= 0;
   const int tilesz = 2048;
   const int stride = 4;
   
   #pragma omp parallel for  
   for(int diag = 0; diag < n-16*m; diag += 2048){
      for(int offset = 0; offset < n-diag-16*m; offset+=2048){
         t ++;
         xcorr_kern(cx,dx,xf,s,output,outputi,diag,offset);
      }
   }   
   printf("iterations: %lu \n",t);
}


