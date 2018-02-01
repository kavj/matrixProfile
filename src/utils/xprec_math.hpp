#include<cmath>
#include<stdint.h>
#include "../arch/avx2arith_2.hpp"
//#include <cmath>
// Collection of in place 2 operand arithmetic and sliding window sum and average. These are based on 
// ACCURATE SUM AND DOT PRODUCT Ogita et al

// The first argument is always updated with the resulting sum or product and the second is updated to reflect the error term.
//using namespace vmth;


template<typename dtype>
static inline void xadd(dtype &a, dtype &b){
   dtype c = a + b;
   dtype d = c - a;
   b = ((a - (c - d)) + (b - d));
   a = c;
}

static inline void xsplit(double &a,double &b){
   double f = a*((double)(2<<27));
   b = (f - (f - a));
   a = a - f;
}

// this should probably move to a separate header or the top or something
template<typename dtype>
static inline void xmul(dtype &a, dtype &b){
   dtype c = a*b;
   b = fma(a,b,-1*c);
   a = c;
}

#ifdef __FMA__
static inline void xmul(__m256d &a,__m256d &b){
   __m256d c = a*b;
   b = vmth::fms(a,b,c);
   a = c;
}
#endif


void xsum(float *a, float &s, float &e, int len){
   for (int i = 0; i < len; i++){
      float b = a[i];
      xadd(s,b);
      e += b;
   }
}

void xsum(double *a, double &s, double &e, int len){
   for (int i = 0; i < len; i++){
      double b = a[i];
      xadd(s,b);
      e += b;
   }
}


template<typename dtype, typename vtype, typename indextype>
static inline vtype xsInv_simd(dtype *a, dtype *mu, dtype *sI, int offset, indextype winlen){
   vtype z = vmth::uload(mu,offset);
   vtype p = vmth::uload(a,offset) - z;
   vtype s = p;
   xmul(p,s);
 //  dtype u = static_cast<dtype>(winlen);
 //  vtype t = bcast(u);
   for(int j = offset+1; j < offset+winlen; j++){
       vtype h = vmth::uload(a,j) - z;
       vtype r = h;
       xmul(h,r);
       xadd(p,h);
       s += (h+r);
   }
   p += s;
   p = p/vmth::bcast(static_cast<dtype>(winlen));
   p = sqrt(p); 
   return p;
}

template<typename dtype,typename indextype>
static inline dtype xsInv_scalar(dtype *a, dtype *mu, dtype *sI, int offset, indextype winlen){
   dtype z = mu[offset];
   dtype p = a[offset] - z;
   dtype s = p;
   xmul(p,s);
   for(int j = offset+1; j < offset+winlen; j++){
       dtype h = a[j] - z;
       dtype r = h;
       xmul(h,r);
       xadd(p,h);
       s += (h+r);
   }
   p += s;
   p = p/(static_cast<dtype>(winlen));
   p = vmth::sqrt(p); 
   return p;
}


// This file is way too messy. I should just split the vector and 
template<typename dtype>
void xsInv(dtype *a, dtype *mu, dtype *sI, int len, int winlen){
   int k = len - winlen + 1;
   long z = static_cast<long>(winlen);
   int m = k - k%8;
   #pragma omp parallel for
   for(int i = 0; i < m; i+=8){
      
      vmth::ustore(xsInv_simd<dtype,__m256d>(a,mu,sI,i,z),sI,i);
      vmth::ustore(xsInv_simd<dtype,__m256d>(a,mu,sI,i+4,z),sI,i+4);
      
   }
   for(int i = m; i < k; i++){
      sI[i] = xsInv_scalar(a,mu,sI,i,winlen);
   }
}


template<typename dtype>
void init_dfdx(dtype *a, dtype* mu, dtype *df, dtype *dx, int w, int n){
   df[0] = 0; 
   dx[0] = 0;
   for(int i = 0; i < n-w; i++){
      df[i+1] = (a[i+w] - mu[i+1]) + (a[i]-mu[i]);
      dx[i+1] = (1/2)*(a[i+w] - a[i]);
   }
}


template<typename dtype>
void xmean_windowed(dtype *a, dtype *output, int len, int winlen){
   int mlen = len - winlen + 1;
   dtype rsum = 0;
   dtype rerr = 0;
   xsum(a,rsum,rerr,winlen);
   output[0] = (rsum + rerr)/winlen;
   for (int i = winlen; i < len; i++){
      dtype addrem = -1*a[i-winlen];
      xadd(rsum,addrem);
      rerr += addrem; 
      addrem = a[i];
      xadd(rsum,addrem);
      rerr += addrem;
      output[i-winlen+1] = (rsum + addrem)/winlen;
   }
}


