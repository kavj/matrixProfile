#include "../mp/avx2arith_2.hpp"
#include <cstdio>

using namespace vmth;
typedef __m256d vtf;
#define wid 4


//template <typename dtype>
//static inline void ssinit_1x1(dtype *a, dtype *mu, dtype *s, int base, int w){
vtf vcss_1x1(double* a, vtf mu, int base, int w){
   vtf s = setzero();
   for(int i = base; i < base+w; i++){
      vtf b = sub(uload(a,i),mu);
      s = fma(b,b,s);
   }
   return s;
}

//template <typename dtype>
//static inline void muinit_1x1(dtype *a, dtype *mu, int base, int w){
 vtf vsum_1x1(double* a, int base, int w){
   vtf s = setzero();
   for(int i = base+32; i < base+w; i+=4){
      s = add(s,uload(a,i));
   }
   return s;
}

//template <typename dtype>
//static inline void ssinit_8x1(dtype *a, dtype *mu, dtype *s, int base, int w){
static void vcss_8x1(double* a, double*mu, double* s, int base, int w){
   vtf ma = aload(mu,base);
   vtf mb = aload(mu,base+wid);
   vtf mc = aload(mu,base+2*wid);
   vtf md = aload(mu,base+3*wid);
   vtf me = aload(mu,base+4*wid);
   vtf mf = aload(mu,base+5*wid);
   vtf mg = aload(mu,base+6*wid);
   vtf mh = aload(mu,base+6*wid);
   

   vtf sa = vcss_1x1(a,ma,base,7*wid);
   vtf sb = vcss_1x1(a,mb,base+wid,6*wid);
   vtf sc = vcss_1x1(a,mc,base+2*wid,5*wid);
   vtf sd = vcss_1x1(a,md,base+3*wid,4*wid);
   vtf se = vcss_1x1(a,me,base+4*wid,3*wid);
   vtf sf = vcss_1x1(a,mf,base+5*wid,2*wid);
   vtf sg = vcss_1x1(a,mg,base+6*wid,wid);
   vtf sh = setzero();

   for(int i = base+7*wid; i < base+w-7*wid; i++){
      vtf ba = uload(a,i);

      vtf ca = sub(ba,ma);
      vtf cb = sub(ba,mb);
      vtf cc = sub(ba,mc);
      vtf cd = sub(ba,md);
      vtf ce = sub(ba,me);
      vtf cf = sub(ba,mf);
      vtf cg = sub(ba,mg);
      vtf ch = sub(ba,mh);

      sa = fma(ca,ca,sa);
      sb = fma(cb,cb,sb);
      sc = fma(cc,cc,sc);
      sd = fma(cd,cd,sd);
      se = fma(ce,ce,se);
      sf = fma(cf,cf,sf);
      sg = fma(cg,cg,sg);
      sh = fma(ch,ch,sh);
   }
   
   astore(sa,s,base);
   astore(sb,s,base+wid);
   astore(sc,s,base+2*wid);
   astore(sd,s,base+3*wid);
   astore(se,s,base+4*wid);
   astore(sf,s,base+5*wid);
   astore(sg,s,base+6*wid);
   astore(sh,s,base+7*wid);
}

//template <typename dtype>
//static inline void ssinit_8x1(dtype *a, dtype *mu, dtype *s, int base, int w){
static void vcsp_8x1(double* a, double*mu, double* s, int base, int offset, int w){
   vtf sa = uload(a,base);
   vtf sb = uload(a,base+wid);
   vtf sc = uload(a,base+2*wid);
   vtf sd = uload(a,base+3*wid);
   vtf se = uload(a,base+4*wid);
   vtf sf = uload(a,base+5*wid);
   vtf sg = uload(a,base+6*wid);
   vtf sh = uload(a,base+7*wid);

   for(int i = base; i < base+w; i++){
      vtf ba = bcast(a,i);

      vtf ca = sub(ba,uload(mu,i));
      vtf cb = sub(ba,uload(mu,i+wid));
      vtf cc = sub(ba,uload(mu,i+2*wid));
      vtf cd = sub(ba,uload(mu,i+3*wid));
      vtf ce = sub(ba,uload(mu,i+4*wid));
      vtf cf = sub(ba,uload(mu,i+5*wid));
      vtf cg = sub(ba,uload(mu,i+6*wid));
      vtf ch = sub(ba,uload(mu,i+7*wid));

      sa = fma(ca,ca,sa);
      sb = fma(cb,cb,sb);
      sc = fma(cc,cc,sc);
      sd = fma(cd,cd,sd);
      se = fma(ce,ce,se);
      sf = fma(cf,cf,sf);
      sg = fma(cg,cg,sg);
      sh = fma(ch,ch,sh);
   }
   
   astore(sa,s,base);
   astore(sb,s,base+wid);
   astore(sc,s,base+2*wid);
   astore(sd,s,base+3*wid);
   astore(se,s,base+4*wid);
   astore(sf,s,base+5*wid);
   astore(sg,s,base+6*wid);
   astore(sh,s,base+7*wid);
}


//template <typename dtype>
//static inline void muinit_8x1(dtype *a, dtype *mu, int base, int w){
static inline void winmean_8x1(double* a, double* mu, int base, int w){
   vtf sa = vsum_1x1(a,base,7*wid);
   vtf sb = vsum_1x1(a,base+wid,6*wid);
   vtf sc = vsum_1x1(a,base+2*wid,5*wid);
   vtf sd = vsum_1x1(a,base+3*wid,4*wid);
   vtf se = vsum_1x1(a,base+4*wid,3*wid);
   vtf sf = vsum_1x1(a,base+5*wid,2*wid);
   vtf sg = vsum_1x1(a,base+6*wid,wid);
   vtf sh = setzero();

   for(int i = base+7*wid; i < base+w-7*wid; i++){
      vtf ba = uload(a,i);

      sa = add(sa,ba);
      sb = add(sb,ba);
      sc = add(sc,ba);
      sd = add(sd,ba);
      se = add(se,ba);
      sf = add(sf,ba);
      sg = add(sg,ba);
      sh = add(sh,ba);
   }
   vtf v = bcast((double)w); 

   sb = add(sb,vsum_1x1(a,base,wid));
   sc = add(sc,vsum_1x1(a,base,2*wid));
   sd = add(sd,vsum_1x1(a,base,3*wid));
   se = add(se,vsum_1x1(a,base,4*wid));
   sf = add(sf,vsum_1x1(a,base,5*wid));
   sg = add(sg,vsum_1x1(a,base,6*wid)); 
   sh = add(sh,vsum_1x1(a,base,7*wid));

   astore(div(sa,v),mu,base);
   astore(div(sb,v),mu,base+wid);
   astore(div(sc,v),mu,base+2*wid);
   astore(div(sd,v),mu,base+3*wid);
   astore(div(se,v),mu,base+4*wid);
   astore(div(sf,v),mu,base+5*wid);
   astore(div(sg,v),mu,base+6*wid);
   astore(div(sh,v),mu,base+7*wid);

}

static inline vtf winmean_gen(double* a,int base, int w){
   vtf s = setzero();
   for(int i = base; i < base+w; i++){
      s = add(s,uload(a,i));
   }
   return div(s,bcast(static_cast<double>(w)));
}


//template <typename dtype>
//void initstats(dtype *a, dtype *mu, dtype *s, int n, int w){
void initstats(double *a, double *mu, double *s, int n, int w){
   printf("in init stats %d %d\n",n,w);
   #pragma omp parallel for 
   for(int i = 0; i < n-w-31; i+=32){
      winmean_8x1(a,mu,i,w);
      vcss_8x1(a,mu,s,i,w);
   }
}

void initxy(double* a, double* b, double* ma, double* mb, double* c, int n, int w){
   #pragma omp parallel for
   for(int i = 0; i < n-w-31; i+=32){
      
   }
} 


//void initstats(double* a, double* mu, 
