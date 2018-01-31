


#include "refpwed.h" 



void mpkern(double* __restrict__ cx, const double* __restrict__ dx, const double* __restrict__ xf, double* __restrict__ s, double* __restrict__ buf, int diag, int offset){
   cx = (double*)__builtin_assume_aligned(cx,64);
   dx = (double*)__builtin_assume_aligned(dx,64);
   xf = (double*)__builtin_assume_aligned(xf,64);
   s  = (double*)__builtin_assume_aligned(s,64);
   buf =(double*)__builtin_assume_aligned(buf,64);
   
   const int tilesz = 32;
   int m = diag+offset;
   for(int i = offset; i < offset+tilesz; i++){
      double s1 = s[i];
      double dxi = dx[i];
      double xfi = xf[i];      
      for(int j = diag; j < diag+tilesz; j++){
         int k = i+j;
         double c = cx[i];
         c += dxi*xf[k];
         c += dx[k]*xfi;
         cx[j] = c;
         double d = c*s1;
         d *= s[k];
         buf[k-offset] = d;
      }
   }
}

/*
void mpkern2(double* __restrict__ cx, const double* __restrict__ dx, const double* __restrict__ xf, double* __restrict__ s, double* __restrict__ buf, int diag, int offset){
   cx = (double*)__builtin_assume_aligned(cx,64);
   dx = (double*)__builtin_assume_aligned(dx,64);
   xf = (double*)__builtin_assume_aligned(xf,64);
   s  = (double*)__builtin_assume_aligned(s,64);
   buf =(double*)__builtin_assume_aligned(buf,64);
   
   const int tilesz = 512;
   for(int j = diag; j < diag+tilesz; j++){
      for(int i = offset; i < offset+tilesz; i+=10){
         int k = i+j;
         int m = k-offset;
         double ca = cx[i];
         double cb = cx[i+1];
         double cc = cx[i+2];
         double cd = cx[i+3];
         double ce = cx[i+4];
         double cf = cx[i+5];
         double cg = cx[i+6];
         double ch = cx[i+7];
         double ci = cx[i+8];
         double cj = cx[i+9];

         double bcst = dx[i]; 

         ca += bcst*xf[k];
         cb += bcst*xf[k+1];
         cc += bcst*xf[k+2];
         cd += bcst*xf[k+3];
         ce += bcst*xf[k+4];
         cf += bcst*xf[k+5];
         cg += bcst*xf[k+6];
         ch += bcst*xf[k+7];
         ci += bcst*xf[k+8];
         cj += bcst*xf[k+9]; 
 
         bcst = xf[i];

         ca += bcst*dx[k];
         cb += bcst*dx[k+1];
         cc += bcst*dx[k+2];
         cd += bcst*dx[k+3];
         ce += bcst*dx[k+4];
         cf += bcst*dx[k+5];
         cg += bcst*dx[k+6];
         ch += bcst*dx[k+7];
         ci += bcst*dx[k+8];
         cj += bcst*dx[k+9];

         cx[i] = ca;
         cx[i+1] = cb;
         cx[i+2] = cc;
         cx[i+3] = cd;
         cx[i+4] = ce;
         cx[i+5] = cf;
         cx[i+6] = cg;
         cx[i+7] = ch; 
         cx[i+8] = ci;
         cx[i+9] = cj;

         bcst = s[i]; 

         ca *= bcst;
         cb *= bcst;
         cc *= bcst;
         cd *= bcst;
         ce *= bcst;
         cf *= bcst;
         cg *= bcst;
         ch *= bcst;
         ci *= bcst;
         cj *= bcst;
         
         ca *= s[k];
         cb *= s[k+1];
         cc *= s[k+2];
         cd *= s[k+3];
         ce *= s[k+4];
         cf *= s[k+5];
         cg *= s[k+6];
         ch *= s[k+7];
         ci *= s[k+8];
         cj *= s[k+9];

         buf[m] = ca;
         buf[m+1] = cb;
         buf[m+2] = cc;
         buf[m+3] = cd;
         buf[m+4] = ce;
         buf[m+5] = cf;
         buf[m+6] = cg;
         buf[m+7] = ch;
         buf[m+8] = ci;
         buf[m+9] = cj;

      }
   }
}*/


