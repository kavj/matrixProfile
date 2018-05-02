#include<cmath>
#include<stdint.h>
#pragma once

// Collection of in place 2 operand arithmetic and sliding window sum and average. These are based on 
// ACCURATE SUM AND DOT PRODUCT Ogita et al


// Todo: unwind the inline functions in place as I don't care for the reference passing here


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

template<typename dtype>
static inline void xmul(dtype &a, dtype &b){
   dtype c = a*b;
   b = fma(a,b,-1*c);
   a = c;
}

static inline void xsum(double *a, double &s, double &e, int len){
   for (int i = 0; i < len; i++){
      double b = a[i];
      xadd(s,b);
      e += b;
   }
}

template<typename dtype>
void fast_invcn(dtype* __restrict__ invn, const dtype* __restrict__ ts, const dtype* __restrict__ mu, int len, int sublen){
   dtype a = 0;
   for(int i = 0; i < sublen; i++){
      dtype t = ts[i] - mu[0];
      a += t*t;
   }
   invn[0] = 1.0/a;
   for(int i = 1; i < len-sublen+1; i++){
      dtype b = ts[i+sublen-1];
      dtype c = ts[i-1];
      a += ((b - mu[i+sublen-1]) + (c - mu[i-1])) * (b - c); 
      invn[i-sublen+1] = sqrt(1.0/a);
   }
}


template<typename dtype>
dtype invn(dtype *a, dtype *mu, dtype *sI, int offset, int winlen){
   dtype z = mu[offset];
   dtype p = a[offset] - z;
   dtype s = p;
   xmul(p,s);
   for(int j = offset + 1; j < offset+winlen; j++){
       dtype h = a[j] - z;
       dtype r = h;
       xmul(h,r);
       xadd(p,h);
       s += (h+r);
   }
   p += s;
   p = 1.0/sqrt(p); 
   return p;
}


template<typename dtype>
void xsInv(dtype *a, dtype *mu, dtype *sI, int len, int winlen){
   int k = len - winlen + 1;
   #pragma omp parallel for
   for(int i = 0; i < k; i++){
      sI[i] = invn(a,mu,sI,i,winlen);
   }
}

// move to cython layer 
template<typename dtype>
void init_dfdx(dtype *a, dtype* mu, dtype *df, dtype *dx, int w, int n){
   df[0] = 0; 
   dx[0] = 0;
   for(int i = 0; i < n-w; i++){
      df[i+1] = (a[i+w] - mu[i+1]) + (a[i]-mu[i]);
      dx[i+1] = (a[i+w] - a[i])/2.0;
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

