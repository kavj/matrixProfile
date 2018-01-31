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
   ci = aload(cx,diag+8*stride);
   cj = aload(cx,diag+9*stride);
   ck = aload(cx,diag+10*stride);
   cl = aload(cx,diag+11*stride);

   for(int i = p; i < p+1024; i+=4){ 
      vtf fbc = bcast(dx,i);
      ca = fma(fbc,aload(xf,p),ca);
      cb = fma(fbc,aload(xf,p+stride),cb);
      cc = fma(fbc,aload(xf,p+2*stride),cc);
      cd = fma(fbc,aload(xf,p+3*stride),cd);
      ce = fma(fbc,aload(xf,p+4*stride),ce);
      cf = fma(fbc,aload(xf,p+5*stride),cf);
      cg = fma(fbc,aload(xf,p+6*stride),cg);
      ch = fma(fbc,aload(xf,p+7*stride),ch);
      ci = fma(fbc,aload(xf,p+8*stride),ci);
      cj = fma(fbc,aload(xf,p+9*stride),cj);
      ck = fma(fbc,aload(xf,p+10*stride),ck);
      cl = fma(fbc,aload(xf,p+11*stride),cl);
      fbc = bcast(xf,i);
          
      ca = fma(fbc,aload(dx,p),ca);
      cb = fma(fbc,aload(dx,p+stride),cb);
      cc = fma(fbc,aload(dx,p+2*stride),cc);
      cd = fma(fbc,aload(dx,p+3*stride),cd);
      ce = fma(fbc,aload(dx,p+4*stride),ce);
      cf = fma(fbc,aload(dx,p+5*stride),cf);
      cg = fma(fbc,aload(dx,p+6*stride),cg);
      ch = fma(fbc,aload(dx,p+7*stride),ch);
      ci = fma(fbc,aload(dx,p+8*stride),ci);
      cj = fma(fbc,aload(dx,p+9*stride),cj);
      ck = fma(fbc,aload(dx,p+10*stride),ck);
      cl = fma(fbc,aload(dx,p+11*stride),cl);

      vtf da,db,dc,dd,de,df,dg,dh,di,dj,dk,dl;

      fbc = bcast(dx,i);

      da = mul(ca,aload(s,p));
      db = mul(cb,aload(s,p+stride));
      dc = mul(cc,aload(s,p+2*stride));
      dd = mul(cd,aload(s,p+3*stride));
      de = mul(ce,aload(s,p+4*stride));
      df = mul(cf,aload(s,p+5*stride));
      dg = mul(cg,aload(s,p+6*stride));
      dh = mul(ch,aload(s,p+7*stride));
      di = mul(ci,aload(s,p+8*stride));
      dj = mul(cj,aload(s,p+9*stride));
      dk = mul(ck,aload(s,p+10*stride));
      dl = mul(cl,aload(s,p+11*stride));

      da = mul(da,aload(s,p));
      db = mul(db,aload(s,k+stride));
      dc = mul(dc,aload(s,k+2*stride));
      dd = mul(dd,aload(s,k+3*stride));
      de = mul(de,aload(s,k+4*stride));
      df = mul(df,aload(s,k+5*stride));
      dg = mul(dg,aload(s,k+6*stride));
      dh = mul(dh,aload(s,k+7*stride));
      di = mul(di,aload(s,k+8*stride));
      dj = mul(dj,aload(s,k+9*stride));
      dk = mul(dk,aload(s,k+10*stride));
      dl = mul(dl,aload(s,k+11*stride));

      astore(da,output,4*p);
      astore(db,output,4*p+stride);
      astore(dc,output,4*p+2*stride);
      astore(dd,output,4*p+3*stride);
      astore(de,output,4*p+4*stride);
      astore(df,output,4*p+5*stride);
      astore(dg,output,4*p+6*stride);
      astore(dh,output,4*p+7*stride);
      astore(di,output,4*p+8*stride);
      astore(dj,output,4*p+9*stride);
      astore(dk,output,4*p+10*stride);
      astore(dl,output,4*p+11*stride);

   }
   astore(ca,cx,diag);
   astore(cb,cx,diag+stride);
   astore(cc,cx,diag+2*stride);
   astore(cd,cx,diag+3*stride);
   astore(ce,cx,diag+4*stride);
   astore(cf,cx,diag+5*stride);
   astore(cg,cx,diag+6*stride);
   astore(ch,cx,diag+7*stride);
   astore(ci,cx,diag+8*stride);
   astore(cj,cx,diag+9*stride);
   astore(ck,cx,diag+10*stride);
   astore(cl,cx,diag+11*stride);

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


static inline void reduce(double* __restrict__ buf, double*  __restrict__ output, double* __restrict__ outputi, int diag, int k) __attribute__ ((always_inline));


static inline void reduce(double* __restrict__ buf, double*  __restrict__ output, long* __restrict__  outputi, int diag, int k){ 
   int stride = 4;
   for(int i = 0; i < 12*1024; i+= 8*stride){
       vtf ca = aload(buf,i);
       vtf cb = aload(buf,i+stride);
       vtf cc = aload(buf,i+2*stride);
       vtf cd = aload(buf,i+3*stride);
       vtf ce = aload(buf,i+4*stride);
       vtf cf = aload(buf,i+5*stride);
       vtf cg = aload(buf,i+6*stride);
       vtf ch = aload(buf,i+7*stride);
       ca = max(ca,cb);
       cc = max(cc,cd);
       ce = max(ce,cf);
       cg = max(cg,ch);
       vti ones = bcast(1);
       int p = diag+k;
       vti base = set(p,p+1,p+2,p+3);
       cb = cmpgtr(cb,ca);
       cd = cmpgtr(cd,ca);
       cf = cmpgtr(ce,cf);
       ch = cmpgtr(cg,ch);

       vti da = set(p,p+1,p+2,p+3);
       vti db = da + ones;
       da = blend(da,db,cb);
       vti dc = db + ones;
       vti dd = dc + ones;
       dc = blend(dc,dd,cd);
       vti de = dd + ones;
       vti df = de + ones;
       de = blend(de,df,cf);
       vti dg = df + ones;
       vti dh = dg + ones;
       dh = blend(dg,dh,ch);
       
       ca = max(ca,cc);
       ce = max(ce,cg);
       cb = cmpgtr(ca,cc);
       cc = cmpgtr(ce,cg);

       ca = max(ca,ce);
       cd = cmpgtr(ca,ce);
       astore(ca,output,i);
       astore(dh,outputi,i);
   }

} 


void  accumTest4_7_10(double* cx,double* dx, double* xf, double*s, double* buf, double* output, long* outputi,int n, int m){
   unsigned long t= 0;
   //printf("n
   #pragma omp parallel for  
   for(int diag = 0; diag < n-4*m; diag += 48){
     // for(int offset = 0; offset < n-4*m; offset+=2048){
         int z = diag;
         for(int k = 0; k < n-4*m; k+=1024){
            xcorr_2(cx, dx, xf, s, buf, diag, k);
            t += 49152;
            //reduce(buf,output,outputi,diag,k);
         }
     //}
   }     
   printf("iterations: %lu\n",t);
}


