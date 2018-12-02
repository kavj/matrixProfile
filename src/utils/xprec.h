#include<math.h>
#pragma once

// Collection of in place 2 operand arithmetic and sliding window sum and average. These are based on 
// ACCURATE SUM AND DOT PRODUCT Ogita et al


static inline void xadd (double a, double b){
   double c = a + b;
   double d = c - a;
   b = ((a - (c - d)) + (b - d));
   a = c;
}


static inline void xmul (double a, double b){
   double c = a * b;
   b = fma(a, b, -1 * c);
   a = c;
}


void fast_invcn(double* restrict invn, const double* restrict ts, const double* restrict mu, int len, int sublen){
   double a = 0;
   for(int i = 0; i < sublen; i++){
      double t = ts[i] - mu[0];
      a += t * t;
   }
   invn[0] = sqrt(1.0/a);
   for(int i = 1; i < len - sublen + 1; i++){
      double b = ts[i + sublen - 1];
      double c = ts[i - 1];
      a += ((b - mu[i + sublen - 1]) + (c - mu[i - 1])) * (b - c); 
      invn[i - sublen + 1] = 1.0/sqrt(a);
   }
}


double invcn(const double* restrict ts, const double* restrict sI, double z, int winlen){
   double p = ts[0] - z;
   double s = p;
   xmul(p,s);
   for(int i = 1; i < winlen; i++){
       double h = ts[i] - z;
       double r = h;
       xmul(h, r);
       xadd(p, h);
       s += (h + r);
   }
   p += s;
   p = 1.0/sqrt(p); 
   return p;
}


void xInvn(const double* restrict ts, const double* restrict mu, double* restrict sI, int len, int winlen){
   int mlen = len - winlen + 1;
   #pragma omp parallel for
   for(int i = 0; i < mlen; i++){
      sI[i] = invcn(ts + i, sI + i, mu[i], winlen);
   }
}


void xmean_windowed(const double* restrict ts, double* restrict mu, int len, int winlen){
   double u, v, w, x, y;
   u = ts[0];
   v = 0;
   for(int i = 1; i < winlen; i++){
      w = u;
      y = ts[i];
      u += y;
      x = u - w;
      v += (w - (u - x)) + (y - x);
   }
   mu[0] = (u + v)/winlen;
   for(int i = winlen; i < len; i++){
      w = u; 
      y = -1 * ts[i - winlen];
      u += y;
      x = u - w;
      v += (w - (u - x)) + (y - x);

      w = u;
      y = ts[i];
      u += y;
      x = u - w;
      v += (w - (u - x)) + (y - x);
      mu[i - winlen + 1] = (u + v)/winlen;
   }
}


void sum_windowed(const double* restrict ts, double* restrict s, int len, int winlen){
   int mlen = len - winlen + 1;
   double u, v, w, x, y;
   u = ts[0];
   v = 0;
   for(int i = 1; i < winlen; i++){
      w = u;
      y = ts[i];
      u += y;
      x = u - w;
      v += (w - (u - x)) + (y - x);
   }
   s[0] = u+v;
   for(int i = winlen; i < len; i++){
      w = u; 
      y = -1 * ts[i - winlen];
      u += y;
      x = u - w;
      v += (w - (u - x)) + (y - x);

      w = u;
      y = ts[i];
      u += y;
      x = u - w;
      v += (w - (u - x)) +(y - x);
      s[i - winlen + 1] = u + v;
   }
}
