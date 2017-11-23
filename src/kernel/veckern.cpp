#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 16 
#define tileSz 256
#define alignedInner (innerWid/4-1)

using namespace vmth;
#define ksz 64
typedef  __m256d vtf;
typedef  __m256i vti;
/*
void accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    
    unsigned long t3 = 0;
    for(int i = 0; i < n; i +=tileSz){
        for(int j = i; j < n;j+=4){
            for(int k = 0; k < tileSz; k+=innerWid){
                vtf a1 = bcast(dF,i+k); 
                vtf a2 = bcast(dF,j+k);
                for(int l = 0; l < innerWid; l+=4){
                   storea(mult_add(a2,loada(dX,j+l),mult_add(a1,loada(dX,i+l),loada(Cxy,i+l))),Cxy,i+l);
                }
                t3+=40;
          }
       }
    }
    printf("%lu %d \n",t3,n);
}*/


/*
void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    double a[innerWid];
    double b[innerWid];
    double e[innerWid];
    double f[innerWid];
    unsigned long t3 = 0;
    for(int i = 0; i < n; i +=tileSz){
        for(int j = i; j < n;j+=4){
            for(int k = 0; k < tileSz; k+=innerWid){
                vtf a1 = bcast(dF,i+k); 
                vtf a2 = bcast(dF,j+k);
                vtf s2 = bcast(s,i+k);
                for(int l = 0; l < innerWid; l++){
                    storea(mult_add(a1,loada(dX,i+k+l),loada(Cxy,i+k+l)),a,l);
                    vtf c1 = mult_add(a2,loada(dX,j+k+l),mult_add(a1,loada(dX,i+k+l),loada(Cxy,i+k+l)));
                    storea(c1,output,i+k+l); 
                    storea(c1,a,l);
                        
                    storea(mult_add(a2,loada(dX,j+k+l),loada(a,l)),output,i+k+l);
                }
                for(int l = 0; l < innerWid; l++){
                    vtf c1 = mult(s2,loadu(s,j+k+l));
                }
           }
            
        }
    }
    printf("%lu %d \n",t3,n);
}
*/

static inline void  kernelTest(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int k){
    double a[tileSz*innerWid];
    double b[tileSz*innerWid];
    int t3 = 0;
    vtf c1 = loada(Cxy,0);
    for(int l = 0; l < innerWid; l+=4){
       vtf c2 = loada(dF,l);
       vtf c3 = loada(dX,l);
       vtf c4 = bcast(dF,l);
       vtf c5 = bcast(dX,l);
       vtf c6 = bcast(s,l);
       vtf c7 = loada(s,l);
       c1 = mult_add(c4,c5,mult_add(c2,c3,c1));
       vtf c8 = mult(c7,mult(c6,c1));
       storea(c8,a,l);
       t3 += 4;
    }
    storea(c1,Cxy,0);
    for(int l = 0; l < innerWid; l+=4){
          vtf c2 = loada(output,l);
          vtf c3 = loada(output,l);
          vti c4 = loada(outputI,l);
          vti c5 = loada(outputI,l);
          vtf c6 = loada(a,k+l);
          vti c7 = set(l,l+1,l+2,l+3);
          vti c8 = bcast(l);
    }
    for(int l = 0; l < innerWid; l+=4){
       storea(loada(a,l),output,k);
    }
    //printf("%lu %i \n",t3,k);
}


void callKern(double* Cxy, double* dX, double* dF, double* s, double* output, long* outputI, int n){
   int t3 = 0;
   for(int i = 0; i < n; i+=tileSz){
      for(int j = i; j < n; j+=tileSz){
         kernelTest(Cxy+j,dX+i,dF+j,s,output,outputI,i);
         t3 += 4;
      }
   }
   printf("%lu\n",t3);
}


void  accumTest6(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    double a[tileSz*innerWid*4];
    double b[tileSz*innerWid*4];
    int t3 = 0;
    for(int i = 0; i < n/2; i += tileSz){
        for(int j = i; j < n; j+=tileSz){
            for(int k = 0; k < tileSz; k+=4){
                vtf c1 = loada(Cxy,j+k*innerWid);
                vtf F1 = loada(dF,i+k*innerWid);
                vtf X1 = loada(dX,i+k*innerWid);
                vtf F2 = bcast(dF,j+k*innerWid);
                vtf X2 = bcast(dX,j+k*innerWid);
                vtf F3 = preshuffle(F2,loada(dF,i+k*innerWid+4));
                vtf X3 = preshuffle(X2,loada(dX,j+k*innerWid+4));
                for(int l = 0; l < innerWid; l+=4){
                    c1 = mult_add(F2,X2,mult_add(F1,X1,c1)); 
                    storea(c1,a,l+k*innerWid);
                    F1 = shiftOnly1(F1,F3);
                    X1 = shiftOnly1(X1,X3);
                    F2 = shiftOnly1(F1,F3);
                    X2 = shiftOnly1(X2,X3);
                    t3 += 4;
                }
                storea(c1,Cxy,j+k);
            }
            for(int k = 0; k < tileSz; k+=4){
               for(int l = 0; l < innerWid; l+=4){
                    vtf c2 = loada(output,i+k+l);
                    vtf c3 = loada(output,j+k+l);
                    vti c4 = loada(outputI,i+k+l);
                    vti c5 = loada(outputI,j+k+l);
                    vtf c6 = loada(a,k+l);
                    vti c7 = set(i+k+l,i+k+l+1,i+k+l+2,i+k+l+3);
                    vti c8 = bcast(j+k+l);

                }
                storea(loada(a,k),output,k);
            }
        }

    }
    printf("%lu %d \n",t3,n);
}

double* recursiveTest(double* Cxy, double* dX, double* dF, double* s, double* output, long* outputI, int n,int m){
    double* a;
    double* b;
    double* c;
    long t3 = 0;
    posix_memalign((void**)&a,64,4*n*sizeof(double));
    posix_memalign((void**)&b,64,4*n*sizeof(double));
    posix_memalign((void**)&c,64,4*n*sizeof(double));

    n -= n%tileSz; //chopping unaligned portion for now for tesst reasons
    for(int i = 0; i < n-m; i+= m){
      for(int j = i; j < n-m; j += m){
         for(int k = 0; k < 4*m; k+=4){
            vtf c1 = bcast(Cxy,i+k);
            vtf c2 = bcast(Cxy,j+k);
            for(int l = 0; l < innerWid; l+=4){
               vtf c2 = loada(dF,i+k*innerWid+l);
               vtf c3 = bcast(dX,j+k*innerWid+l);
               vtf c4 = loada(dF,i+k*innerWid+l);
               vtf c5 = loada(dX,j+k*innerWid+l);
               c1 = mult_add(c2,c3,c1);
               c2 = mult_add(c4,c5,c2);
               storea(c1,a,k+l);
               storea(c2,b,k+l);
               t3 += 4;
            }
         }
         for(int k = 0; k < 4*m; k+=4){
            storea(add(loada(a,k),loada(b,k)),output,i+k);
         }
      }
   } 
   printf("%lu\n",t3);
   return a;
}


void  accumTest5(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    double a[tileSz*innerWid];
    double b[tileSz*innerWid];
    int t3 = 0;
    for(int i = 0; i < n; i += tileSz){
        for(int j = i; j < n; j+=tileSz){
            for(int k = 0; k < tileSz; k+=4){
                vtf c1 = loada(Cxy,j+k);
                for(int l = 0; l < innerWid; l+=4){
                    vtf c2 = loadu(dF,i+k*innerWid+l);
                    vtf c3 = loadu(dX,i+k*innerWid+l);
                    vtf c4 = bcast(dF,j+k*innerWid+l);
                    vtf c5 = bcast(dX,j+k*innerWid+l);
                    vtf c6 = bcast(s,i+k*innerWid+l);
                    vtf c7 = loadu(s,j+k*innerWid+l);
                    c1 = mult_add(c4,c5,mult_add(c2,c3,c1));
                    vtf c8 = mult(c7,mult(c6,c1));
                    storea(c8,a,l+k*innerWid);
                    t3 += 4;
                }
                storea(c1,Cxy,j+k);
            }
            for(int k = 0; k < tileSz; k+=4){
               for(int l = 0; l < innerWid; l+=4){
                    vtf c2 = loada(output,i+k+l);
                    vtf c3 = loada(output,j+k+l);
                    vti c4 = loada(outputI,i+k+l);
                    vti c5 = loada(outputI,j+k+l);
                    vtf c6 = loada(a,k+l);
                    vti c7 = set(i+k+l,i+k+l+1,i+k+l+2,i+k+l+3);
                    vti c8 = bcast(j+k+l);
                    
                }
                storea(loada(a,k),output,k);
            }
        }

    }
    printf("%lu %d \n",t3,n);
}



/*
static inline void kernTest(double* Cxy, double* dX, double* dF, double* s, double* output, long* outputI, int cInd, int tileRow, int tileCol){
   double a[tileSz];
   double b[tileSz];
   for(int k = 0; k < tileSz; k+=innerWid){  // blocking over parallel iterations
      vtf c = loada(Cxy,cInd+k);
      vtf d1 = loada(output,tileRow,
      for(int l = 0; l < innerWid; l+=simdWid){
         
      }
   }
}*/


// bad
void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    double a[tileSz];
    double b[tileSz];
    unsigned long t3 = 0;

    // Variable i should determine the index of Cxy, For each of those, we need to actually consider all corresponding segments
    //int tCnt = n/tileSz;  // aligned tile count
    //
    // i selects the diagonal index as an index of 
    for(int j = 0; j < n; j +=tileSz){
        for(int i = 0; i < j;i+=tileSz){
            // k and l operate over the innermost loop only
            for(int k = 0; k < tileSz; k+=innerWid){
                vtf c1 = loada(Cxy,i);
                for(int l = 0; l < innerWid; l+=4){
                    c1 = mult_add(bcast(dX,j+k+l),loadu(dF,i+k+l),c1);
                    c1 = mult_add(bcast(dF,j+k+l),loadu(dX,i+k+l),c1);
                    storea(c1,a,l+k);
                }
                storea(c1,Cxy,j+k);
                vtf s2 = bcast(s,j+k);
                for(int l = 0; l < innerWid; l+=4){
                    vtf c1 = loada(a,l);
                    storea(mult(loada(a,l),mult(s2,loadu(s,i+k+l))),b,l);
                }
                vtf outi = bcast(output,i+k);
                vti outIi = loada(outputI,i+k);
                vti out = set(i+k,i+k+1,i+k+2,i+k+3);
                for(int l = 0; l < innerWid; l+=4){
                    vtf c3 = loada(b,l);
                    vtf c2 = cmpgtr(outi,c3);
                    outi = blend(outi,c3,c2);
                    outIi = blend(outIi,out,c2);
                }
                storea(outi,output,i+k);
                storea(outIi,outputI,i+k);
                outi = loada(output,j+k);
                outIi = loada(outputI,j+k);
                for(int l = 4; l < innerWid-4; l+=4){
                for(int r = 1; r < 4; r++){
                    int offs = r*innerWid + r;
                    vtf c3 = loadu(b,offs);
                    vti c4 = set(i+k,i+k+1,i+k+2,i+k+3);
                    vtf c2 = cmpgtr(outi,c3);
                    outi = blend(outi,c3,c2);
                    outIi = blend(outIi,c4,c2);
                }
                storea(outi,output,j+k);
                storea(outIi,outputI,j+k);
                }
                t3 += innerWid;
            }
        }    
    }
    printf("%lu %d \n",t3,n);
}



/*

void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
    double a[innerWid];
    double b[innerWid];
    double e[innerWid];
    double f[innerWid];
    unsigned long t3 = 0;
    for(int i = 0; i < n; i +=tileSz){
        for(int j = i; j < n;j+=4){
            for(int k = 0; k < tileSz; k+=innerWid){
                vtf c1 = loada(Cxy,j+k);
                for(int l = 0; l < innerWid; l+=4){
                    c1 = mult_sub(bcast(dX,j+k+l),loadu(dF,i+k+l),c1);
                    c1 = mult_add(bcast(dF,j+k+l),loadu(dX,i+k+l),c1);
                    storea(c1,a,l);
                }
                storea(c1,Cxy,j+k);
                vtf s2 = bcast(s,j+k);
                for(int l = 0; l < innerWid; l+=4){
                    vtf c1 = loada(a,l);
                    storea(mult(loada(a,l),mult(s2,loadu(s,i+k+l))),b,l);
                }
                vtf outi = bcast(output,i+k);
                vti outIi = loada(outputI,i+k);
                vti out = set(i+k,i+k+1,i+k+2,i+k+3);
                for(int l = 0; l < innerWid; l+=4){
                    vtf c3 = loada(b,l);
                    vtf c2 = cmpgtr(outi,c3);
                    outi = blend(outi,c3,c2);
                    outIi = blend(outIi,out,c2);
                }
                storea(outi,output,i+k);
                storea(outIi,outputI,i+k);
                outi = loada(output,j+k);
                outIi = loada(outputI,j+k);
                for(int l = 4; l < innerWid-4; l+=4){
                for(int r = 1; r < 4; r++){
                    int offs = r*innerWid + r;
                    vtf c3 = loadu(b,offs);
                    vti c4 = set(i+k,i+k+1,i+k+2,i+k+3);
                    vtf c2 = cmpgtr(outi,c3);
                    outi = blend(outi,c3,c2);
                    outIi = blend(outIi,c4,c2);
                }
                storea(outi,output,j+k);
                storea(outIi,outputI,j+k);
                }
                storea(outi,output,i+k);
                storea(outIi,outputI,i+k);
                t3 += innerWid;
            }

        }    
    }
    printf("%lu %d \n",t3,n);
}

*/
  /* for(int r = 4; r < innerWid-4; r+=4){
                    vtf outj = loada(output,j+k);
                    vti outIj = loada(outputI,j+k);
                    // need to pass over simd width here
                    // it goes 1 2 3 4
                    //           2 3 4
                    //             3 4
                    //               4
                    // but this only applies to the first in a series
                    // and the rest are complete. Rather than repeatedly blend, it might make sense to reduce the partials in their own pass
                    //
                    // the reason for this is that in a longer pass it becomes, 
                    //    1 2 3 4 5 6 7 8
                    //      2 3 4 5 6 7 8 9
                    //        3 4 5 6 7 8 9 10
                    //          4 5 6 7 8 9 10 11
                    //  it relies on unaligned loads, but only the first segment is actually incomplete 
                    //  This can be 8, 16, or 32 in width and everything in the interior can form a contiguous chunk.
                    
                }
                storea(outi,output,i+k);
                storea(outIi,outputI,i+k);
               */ 


/*
void  accumTest5(double* Cxy, double* dX, double* dF, double*s, double* output,long* outputI,int n,int m){
    double a[256];
    unsigned long t3 = 0;
    int iters = 1;
    for(int i = 0; i < (n-m); i +=32){
        for(int j = 0; j < i;j++){
            vtf c1 = loada(Cxy,j);
            vtf indI = set(i);
            vtf indJ = set(j,j+1,j+2,j+3);
            vtf inc = set(1);
            for(int k = 0; k < 32; k+=4){
                vtf a1 = bcast(dF,i+k);
                vtf a2 = bcast(dX,j+k);
                vtf b1 = loadu(dX,i+k);
                vtf b2 = loadu(dF,j+k);
       
                                
 
                __m256d op1 = _mm256_broadcast_sd(dF+i+k);
                __m256d op2 = _mm256_loadu_pd(dX+i+k);
                __m256d op3 = _mm256_broadcast_sd(dF+j+k+m);
                __m256d op4 = _mm256_loadu_pd(dX+j+k+m);
                __m256i ind1 = _mm256_loadu_si256((__m256i*)outputI);
                __m256i ind2 = _mm256_loadu_si256((__m256i*)outputI+m);
                c = _mm256_fmadd_pd(op1,op2,c);
                c = _mm256_fmadd_pd(op3,op4,c);
                _mm256_store_pd(a+k,c);
                __m256d g = _mm256_mul_pd(_mm256_loadu_pd(s+j+k+m),_mm256_mul_pd(_mm256_loadu_pd(a+k),_mm256_loadu_pd(s+j+k)));
                __m256d f = _mm256_loadu_pd(output+j);
                __m256i t = _mm256_loadu_si256((__m256i*)outputI+j);
                __m256i t2 = _mm256_loadu_si256((__m256i*)outputI+j+m);
                __m256d h = _mm256_loadu_pd(output+k+m);
                __m256d d1 = _mm256_cmp_pd(g,f,_CMP_GT_OQ);
                __m256d d2 = _mm256_cmp_pd(g,h,_CMP_GT_OQ);          
                g = _mm256_max_pd(g,f);
                h = _mm256_max_pd(g,h);
                ind1 = _mm256_blendv_epi8(ind3,ind1,(__m256i)d1);
                ind2 = _mm256_blendv_epi8(ind3,ind2,(__m256i)d2);
                _mm256_store_pd(output+k,g);
                _mm256_store_pd(output+k+m,h);
                t3+=4;
            }
        }
        _mm256_store_pd(Cxy+i,c);
    }
    printf("%lu \n",t3,n,n-m+1);
}

*/
