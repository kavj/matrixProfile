#include<cstdio>
#include "../utils/xprec_math.hpp"
//#include "../utils/blockarith.hpp"
#include "../arch/avx2arith_2.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
#define vwidth 4
typedef  __m256d vtf;
typedef  __m256i vti;

static inline void xcorr_2(double* __restrict__ cx, double* __restrict__ dx, double* __restrict__ xf, double* __restrict__ s, double*  __restrict__ output, int diag, int k) __attribute__ ((always_inline));
 
static inline void xcorr_2(double *cx, double *dx, double *xf, double *s, double *output, int diag, int k){
   int p = diag+k;
   int stride = 4;
   vtf ca,cb,cc,cd,ce,cf,cg,ch,ci,cj,ck,cl;    
   ca = aload(cx,diag);
   cb = aload(cx,diag+stride);
   cc = aload(cx,diag+2*stride);
   cd = aload(cx,diag+3*stride);
   ce = aload(cx,diag+4*stride);
   cf = aload(cx,diag+5*stride);
   cg = aload(cx,diag+6*stride);
   ch = aload(cx,diag+7*stride);

   vtf fbc = bcast(dx,k);
   ca = fma(fbc,aload(xf,p),ca);
   cb = fma(fbc,aload(xf,p+stride),cb);
   cc = fma(fbc,aload(xf,p+2*stride),cc);
   cd = fma(fbc,aload(xf,p+3*stride),cd);
   ce = fma(fbc,aload(xf,p+4*stride),ce);
   cf = fma(fbc,aload(xf,p+5*stride),cf);
   cg = fma(fbc,aload(xf,p+6*stride),cg);
   ch = fma(fbc,aload(xf,p+7*stride),ch);

   fbc = bcast(xf,k);
          
   ca = fma(fbc,aload(dx,p),ca);
   cb = fma(fbc,aload(dx,p+stride),cb);
   cc = fma(fbc,aload(dx,p+2*stride),cc);
   cd = fma(fbc,aload(dx,p+3*stride),cd);
   ce = fma(fbc,aload(dx,p+4*stride),ce);
   cf = fma(fbc,aload(dx,p+5*stride),cf);
   cg = fma(fbc,aload(dx,p+6*stride),cg);
   ch = fma(fbc,aload(dx,p+7*stride),ch);
              
   astore(ca,cx,diag);
   astore(cb,cx,diag+stride);
   astore(cc,cx,diag+2*stride);
   astore(cd,cx,diag+3*stride);
   astore(ce,cx,diag+4*stride);
   astore(cf,cx,diag+5*stride);
   astore(cg,cx,diag+6*stride);
   astore(ch,cx,diag+7*stride);

   fbc = bcast(dx,k);

   vtf da,db,dc,dd,de,df,dg,dh;
   da = aload(cx,diag);
   db = aload(cx,diag+stride);
            dc = aload(cx,diag+2*stride);
            dd = aload(cx,diag+3*stride);
            de = aload(cx,diag+4*stride);
            df = aload(cx,diag+5*stride);
            dg = aload(cx,diag+6*stride);
            dh = aload(cx,diag+7*stride);
/*
            ca = mul(ca,aload(s,p));
            cb = mul(cb,aload(s,p+stride));
            cc = mul(cc,aload(s,p+2*stride));
            cd = mul(cd,aload(s,p+3*stride));
            ce = mul(ce,aload(s,p+4*stride));
            cf = mul(cf,aload(s,p+5*stride));
            cg = mul(cg,aload(s,p+6*stride));
            ch = mul(ch,aload(s,p+7*stride));
*/
            ca = mul(ca,aload(s,p));
            cb = mul(cb,aload(s,k+stride));
            cc = mul(cc,aload(s,k+2*stride));
            cd = mul(cd,aload(s,k+3*stride));
            ce = mul(ce,aload(s,k+4*stride));
            cf = mul(cf,aload(s,k+5*stride));
            cg = mul(cg,aload(s,k+6*stride));
            ch = mul(ch,aload(s,k+7*stride));

            da = fma(fbc,aload(dx,p),ca);
            db = fma(fbc,aload(dx,p+stride),cb);
            dc = fma(fbc,aload(dx,p+2*stride),cc);
            dd = fma(fbc,aload(dx,p+3*stride),cd);
            de = fma(fbc,aload(dx,p+4*stride),ce);
            df = fma(fbc,aload(dx,p+5*stride),cf);
            dg = fma(fbc,aload(dx,p+6*stride),cg);
            dh = fma(fbc,aload(dx,p+7*stride),ch);
            

            da = fma(fbc,aload(xf,p),ca);
            db = fma(fbc,aload(xf,p+stride),cb);
            dc = fma(fbc,aload(xf,p+2*stride),cc);
            dd = fma(fbc,aload(xf,p+3*stride),cd);
            de = fma(fbc,aload(xf,p+4*stride),ce);
            df = fma(fbc,aload(xf,p+5*stride),cf);
            dg = fma(fbc,aload(xf,p+6*stride),cg);
            dh = fma(fbc,aload(xf,p+7*stride),ch);
/*
            da = mul(da,aload(s,p));
            db = mul(db,aload(s,p+stride));
            dc = mul(dc,aload(s,p+2*stride));
            dd = mul(dd,aload(s,p+3*stride));
            de = mul(de,aload(s,p+4*stride));
            df = mul(df,aload(s,p+5*stride));
            dg = mul(dg,aload(s,p+6*stride));
            dh = mul(dh,aload(s,p+7*stride));
  */        
            da = mul(da,aload(s,p));
            db = mul(db,aload(s,k+stride));
            dc = mul(dc,aload(s,k+2*stride));
            dd = mul(dd,aload(s,k+3*stride));
            de = mul(de,aload(s,k+4*stride));
            df = mul(df,aload(s,k+5*stride));
            dg = mul(dg,aload(s,k+6*stride));
            dh = mul(dh,aload(s,k+7*stride));
          

            ca = max(ca,cb);
            cc = max(cc,cd);
            ce = max(ce,cf);
            cg = max(cg,ch);
            da = max(da,db);
            dc = max(dc,dd);
            de = max(de,df);
            dg = max(dg,dh);
            ca = max(ca,cc);
            ce = max(ce,cg);
            de = max(da,db);
            df = max(dc,dd);
            ca = max(ca,ce);
            dg = max(de,df);
            ca = max(ca,dg); 
            astore(ca,output,diag+8*k);
}


/*
static inline void xcorr_3(double *cx, double *dx, double *xf, double *s, double *output, int diag, int k) __attribute__ ((always_inline));


static inline void xcorr_3(double *cx, double *dx, double *xf, double *s, double *output, int diag, int k){
   int stride = 4;
   int p = offset+k;
   vtf ca = aload(cx,diag);
   vtf cb = aload(cx,diag+stride);
   vtf cc = aload(cx,diag+2*stride);
   vtf cd = aload(cx,diag+3*stride);
   vtf f1 = bcast(dx,k);
   for(int i = diag; i < diag+512; i += 4){

      ca = fma(f1,aload(xf,p));
      cb = fma(f1,aload(xf,p+stride));
      cc = fma(f1,aload(xf,p+2*stride));
      cd = fma(f1,aload(xf,p+3*stride));

      f1 = bcast(xf,k);
      ca = fma(f1,aload(xf,p));
      cb = fma(f1,aload(xf,p+stride));
      cc = fma(f1,aload(xf,p+2*stride));
      cd = fma(f1,aload(xf,p+3*stride));

      ca = mul(ca,aload(s,p));
      cb = mul(cb,aload(s,p+stride));
      cc = mul(cc,aload(s,p+2*stride));
      cd = mul(cc,aload(s,p+3*stride));

      f1 = bcast(s,k);   
      ca = mul(ca,f1);
      cb = mul(cb,f1);
      cc = mul(cc,f1);
      cd = mul(cd,f1);

      
   }
   
}
*/


void  accumTest4_7_10(double* cx,double* dx, double* xf, double*s, double* buf, double* output, long* outputi,int n, int m){
   unsigned long t= 0;
   //printf("n
   #pragma omp parallel for  
   for(int diag = 0; diag < n-4*m; diag += 32){
     // for(int offset = 0; offset < n-4*m; offset+=2048){
         int z = diag;
         for(int k = 0; k < n-4*m; k+=4){
            xcorr_2(cx, dx, xf, s, output, diag, k);
            t += 64;
         }
     //}
   }     
   printf("iterations: %lu\n",t);
}


