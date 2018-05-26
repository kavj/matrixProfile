#ifndef CHECK_ARRAY
#define CHECK_ARRAY

#include<cstdio>
#include<cmath>

inline bool chk(double* a, int n){
   for(int i = 0; i < n; i++){
      if(isnan(a[i])){
         printf("a[%d]:%lf\n",a[i]);
      }
   }
}

inline bool chk(double a, int i){
   if(isnan(a)){
      printf("a[%d]:%lf\n",i,a);
   }
}

inline bool chk(double a){
   return isnan(a);
}

#endif


