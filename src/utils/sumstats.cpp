#include "../mp/avx2arith_2.hpp"
#include <cstdio>

using namespace vmth;
typedef __m256d vtf;
#define vwidth 4



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
   vtf mb = aload(mu,base+vwidth);
   vtf mc = aload(mu,base+2*vwidth);
   vtf md = aload(mu,base+3*vwidth);
   vtf me = aload(mu,base+4*vwidth);
   vtf mf = aload(mu,base+5*vwidth);
   vtf mg = aload(mu,base+6*vwidth);
   vtf mh = aload(mu,base+6*vwidth);
   

   vtf sa = vcss_1x1(a,ma,base,28);
   vtf sb = vcss_1x1(a,mb,base+4,24);
   vtf sc = vcss_1x1(a,mc,base+8,20);
   vtf sd = vcss_1x1(a,md,base+12,16);
   vtf se = vcss_1x1(a,me,base+16,12);
   vtf sf = vcss_1x1(a,mf,base+20,8);
   vtf sg = vcss_1x1(a,mg,base+24,4);
   vtf sh = setzero();

   for(int i = base; i < base+w+28; i++){
      vtf ba = uload(a,i);

      vtf ca = sub(ba,uload(mu,i));
      vtf cb = sub(ba,uload(mu,i+vwidth));
      vtf cc = sub(ba,uload(mu,i+2*vwidth));
      vtf cd = sub(ba,uload(mu,i+3*vwidth));
      vtf ce = sub(ba,uload(mu,i+4*vwidth));
      vtf cf = sub(ba,uload(mu,i+5*vwidth));
      vtf cg = sub(ba,uload(mu,i+6*vwidth));
      vtf ch = sub(ba,uload(mu,i+7*vwidth));

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
   astore(sb,s,base+vwidth);
   astore(sc,s,base+2*vwidth);
   astore(sd,s,base+3*vwidth);
   astore(se,s,base+4*vwidth);
   astore(sf,s,base+5*vwidth);
   astore(sg,s,base+6*vwidth);
   astore(sh,s,base+7*vwidth);
}

//template <typename dtype>
//static inline void ssinit_8x1(dtype *a, dtype *mu, dtype *s, int base, int w){
static void vcsp_8x1(double* a, double*mu, double* s, int base, int offset, int w){
   vtf sa = uload(a,base);
   vtf sb = uload(a,base+4);
   vtf sc = uload(a,base+8);
   vtf sd = uload(a,base+12);
   vtf se = uload(a,base+16);
   vtf sf = uload(a,base+20);
   vtf sg = uload(a,base+24);
   vtf sh = uload(a,base+28);

   for(int i = base; i < base+w; i++){
      vtf ba = bcast(a,i);

      vtf ca = sub(ba,uload(mu,i));
      vtf cb = sub(ba,uload(mu,i+vwidth));
      vtf cc = sub(ba,uload(mu,i+2*vwidth));
      vtf cd = sub(ba,uload(mu,i+3*vwidth));
      vtf ce = sub(ba,uload(mu,i+4*vwidth));
      vtf cf = sub(ba,uload(mu,i+5*vwidth));
      vtf cg = sub(ba,uload(mu,i+6*vwidth));
      vtf ch = sub(ba,uload(mu,i+7*vwidth));

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
   astore(sb,s,base+vwidth);
   astore(sc,s,base+2*vwidth);
   astore(sd,s,base+3*vwidth);
   astore(se,s,base+4*vwidth);
   astore(sf,s,base+5*vwidth);
   astore(sg,s,base+6*vwidth);
   astore(sh,s,base+7*vwidth);
}




//template <typename dtype>
//static inline void muinit_8x1(dtype *a, dtype *mu, int base, int w){
static inline void winmean_8x1(double* a, double* mu, int base, int w){
   vtf sa = vsum_1x1(a,base,28);
   vtf sb = vsum_1x1(a,base+4,24);
   vtf sc = vsum_1x1(a,base+8,20);
   vtf sd = vsum_1x1(a,base+12,16);
   vtf se = vsum_1x1(a,base+16,12);
   vtf sf = vsum_1x1(a,base+20,8);
   vtf sg = vsum_1x1(a,base+24,4);
   vtf sh = setzero();

   for(int i = base+28; i < base+w; i++){
      vtf ba = bcast(a,i);

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

   sb = add(sb,vsum_1x1(a,base,4));
   sc = add(sc,vsum_1x1(a,base,8));
   sd = add(sd,vsum_1x1(a,base,12));
   se = add(se,vsum_1x1(a,base,16));
   sf = add(sf,vsum_1x1(a,base,20));
   sg = add(sg,vsum_1x1(a,base,24)); 
   sh = add(sh,vsum_1x1(a,base,28));

   astore(div(sa,v),mu,base);
   astore(div(sb,v),mu,base+4);
   astore(div(sc,v),mu,base+8);
   astore(div(sd,v),mu,base+12);
   astore(div(se,v),mu,base+16);
   astore(div(sf,v),mu,base+20);
   astore(div(sg,v),mu,base+24);
   astore(div(sh,v),mu,base+28);

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
