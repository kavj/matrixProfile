#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 64 
using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;

/*
void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n-4096; diag += 32){
      vtf z = setzero();
      vtf c1 = loada(Cxy,diag);
      for(int offset = 0; offset < diag; offset++){
         for(int suboff = 0; suboff < 1024; suboff+=4){
         }
         t3+= 32;
         storeu(c1,Cxy,offset);
      }
   }
   printf("%lu %d \n",t3,n);
}*/




void  accumTest4_7_9(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n, int m){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n-4096; diag += 32){
      vtf z = setzero();
      for(int offset = 0; offset < n-diag-32*128; offset+=4){
         vtf t1 = _mm256_set1_pd(0x7fffffff);
         for(int i = 0; i < 32; i += 4){
            vtf c1 = setzero();
            for(int suboff = 0; suboff < 128; suboff+=4){
               vtf b1 = bcast(dx,suboff + offset);
               vtf d1 = max(z,flabs(loada(df,suboff+offset+i)-loada(dx,diag+suboff+i))*loada(s,suboff+offset+i));
               c1 = mult_add(d1,d1,c1);
            }
            storea(c1,output,diag+offset+i);
         }
         t3+= 32;
      }
   }
   printf("%lu %d \n",t3,n);
}


/*
void  accumtest4_7_9(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n-4096; diag += 32){
      vtf z = setzero();
  for(int offset = 0; offset < n-diag-32*128; offset+=4){
          vtf c1 = setzero();
          vtf c2 = setzero();
          vtf c3 = setzero();
          vtf c4 = setzero();
          vtf c5 = setzero();
          vtf c6 = setzero();
          vtf c7 = setzero();
          vtf c8 = setzero();
         for(int suboff = 0; suboff < 128; suboff+=4){
            vtf t1 = _mm256_set1_pd(0x7fffffff);
            vtf b1 = bcast(dx,suboff + offset);
            vtf d1 = max(z,flabs(loada(df,suboff+offset)-loada(dx,diag+suboff))*loada(s,suboff+offset));
            vtf d2 = max(z,flabs(loada(df,suboff+offset+4)-loada(dx,diag+suboff+4))*loada(s,suboff+offset+4));
            vtf d3 = max(z,flabs(loada(df,suboff+offset+8)-loada(dx,diag+suboff+8))*loada(s,suboff+offset+8));
            vtf d4 = max(z,flabs(loada(df,suboff+offset+12)-loada(dx,diag+suboff+12))*loada(s,suboff+offset+12));
            vtf d5 = max(z,flabs(loada(df,suboff+offset+16)-loada(dx,diag+suboff+16))*loada(s,suboff+offset+16));
            vtf d6 = max(z,flabs(loada(df,suboff+offset+20)-loada(dx,diag+suboff+20))*loada(s,suboff+offset+20));
            vtf d7 = max(z,flabs(loada(df,suboff+offset+24)-loada(dx,diag+suboff+24))*loada(s,suboff+offset+24));
            vtf d8 = max(z,flabs(loada(df,suboff+offset+28)-loada(dx,diag+suboff+28))*loada(s,suboff+offset+28));

            c1 = mult_add(d1,d1,c1);
         //   d1 = max(z,b1 - (loada(df,suboff+offset+4)-loada(df,diag+suboff+4))*loada(s,suboff+offset+4));
            c2 = mult_add(d2,d2,c2);
         //   d1 = max(z,b1 - (loada(df,suboff+offset+8)-loada(df,diag+suboff+8))*loada(s,suboff+offset+8));
            c3 = mult_add(d3,d3,c3);
         //   d1 = max(z,b1 - (loada(df,suboff+offset+12)-loada(df,diag+suboff+12))*loada(s,suboff+offset+12));
            c4 = mult_add(d4,d4,c4);
         //   d1 = max(z,b1 - (loada(df,suboff+offset+16)-loada(df,diag+suboff+16))*loada(s,suboff+offset+16));
            c5 = mult_add(d5,d5,c5);
         //   d1 = max(z,b1 - (loada(df,suboff+offset+20)-loada(df,diag+suboff+20))*loada(s,suboff+offset+20));
            c6 = mult_add(d6,d6,c6);
         //   d1 = max(z,b1 - (loada(df,suboff+offset+24)-loada(df,diag+suboff+24))*loada(s,suboff+offset+24));
            c7 = mult_add(d7,d7,c7);
        //    d1 = max(z,b1 - (loada(df,suboff+offset+28)-loada(df,diag+suboff+28))*loada(s,suboff+offset+28));
	    c8 = mult_add(d8,d8,c8);
         }
         t3+= 32;
      
         storea(c1,output,offset);
         storea(c2,output,offset+4);
         storea(c3,output,offset+8);
         storea(c4,output,offset+12);
         storea(c5,output,offset+16);
         storea(c6,output,offset+20);
         storea(c7,output,offset+24);
         storea(c8,output,offset+28);
     }
   }
   printf("%lu %d \n",t3,n);
}*/


/*
void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n-4096; diag += 32){
      vtf z = setzero();
      vtf c1 = setzero();
      vtf c2 = setzero();
      vtf c3 = setzero();
      vtf c4 = setzero();
      vtf c5 = setzero();
      vtf c6 = setzero();
      vtf c7 = setzero();
      vtf c8 = setzero();
   for(int offset = 0; offset < n-diag-32*128; offset+=4){
         for(int suboff = 0; suboff < 128; suboff+=4){
            vtf t1 = _mm256_set1_pd(0x7FFFFFFF);
            vtf b1 = bcast(dX,suboff + offset);
            for(int k = 0; k < 8; k+=4){
               vtf d1 = max(z,flabs(loada(dF,suboff+offset+k)-loada(dX,diag+suboff+k))*loada(s,suboff+offset+k));
              
            }
            vtf t1 = _mm256_set1_pd(0x7FFFFFFF);
            vtf b1 = bcast(dX,suboff + offset);
            vtf d1 = max(z,flabs(loada(dF,suboff+offset)-loada(dX,diag+suboff))*loada(s,suboff+offset));
            vtf d2 = max(z,flabs(loada(dF,suboff+offset+4)-loada(dX,diag+suboff+4))*loada(s,suboff+offset+4));
            vtf d3 = max(z,flabs(loada(dF,suboff+offset+8)-loada(dX,diag+suboff+8))*loada(s,suboff+offset+8));
            vtf d4 = max(z,flabs(loada(dF,suboff+offset+12)-loada(dX,diag+suboff+12))*loada(s,suboff+offset+12));
            vtf d5 = max(z,flabs(loada(dF,suboff+offset+16)-loada(dX,diag+suboff+16))*loada(s,suboff+offset+16));
            vtf d6 = max(z,flabs(loada(dF,suboff+offset+20)-loada(dX,diag+suboff+20))*loada(s,suboff+offset+20));
            vtf d7 = max(z,flabs(loada(dF,suboff+offset+24)-loada(dX,diag+suboff+24))*loada(s,suboff+offset+24));
            vtf d8 = max(z,flabs(loada(dF,suboff+offset+28)-loada(dX,diag+suboff+28))*loada(s,suboff+offset+28));

            c1 = mult_add(d1,d1,c1);
            c2 = mult_add(d2,d2,c2);
            c3 = mult_add(d3,d3,c3);
            c4 = mult_add(d4,d4,c4);
            c5 = mult_add(d5,d5,c5);
            c6 = mult_add(d6,d6,c6);
            c7 = mult_add(d7,d7,c7);
	    c8 = mult_add(d8,d8,c8);
         }
         t3+= 32;
      
         storea(c1,Cxy,offset);
         storea(c2,Cxy,offset+4);
         storea(c3,Cxy,offset+8);
         storea(c4,Cxy,offset+12);
         storea(c5,Cxy,offset+16);
         storea(c6,Cxy,offset+20);
         storea(c7,Cxy,offset+24);
         storea(c8,Cxy,offset+28);
     }
   }
   printf("%lu %d \n",t3,n);
}*/




/*
void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n-4096; diag += 32){
      vtf z = setzero();
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      vtf c5 = loada(Cxy,diag+16);
      vtf c6 = loada(Cxy,diag+20);
      vtf c7 = loada(Cxy,diag+24);
      vtf c8 = loada(Cxy,diag+28);
      for(int offset = 0; offset < diag; offset+=4){
         for(int suboff = 0; suboff < 1024; suboff++){
            vtf b1 = bcast(dX,suboff + offset);
            vtf d1 = b1 - loadu(dX,suboff+offset);
            c1 = mult_add(d1,d1,c1);
            d1 = b1 - loadu(dX,suboff+offset+4);
            c2 = mult_add(d1,d1,c2);
            d1 = b1 - loadu(dX,suboff+offset+8);
            c3 = mult_add(d1,d1,c3);
            d1 = b1 - loadu(dX,suboff+offset+12);
            c4 = mult_add(d1,d1,c4);
            d1 = b1 - loadu(dX,suboff+offset+16);
            c5 = mult_add(d1,d1,c5);
            d1 = b1 - loadu(dX,suboff+offset+20);
            c6 = mult_add(d1,d1,c6);
            d1 = b1 - loadu(dX,suboff+offset+24); 
            c7 = mult_add(d1,d1,c7);
            d1 = b1 - loadu(dX,suboff+offset+28);
            c8 = mult_add(d1,d1,c8);
         }
         t3+= 32;
      
         storea(c1,Cxy,offset);
         storea(c2,Cxy,offset+4);
         storea(c3,Cxy,offset+8);
         storea(c4,Cxy,offset+12);
         storea(c5,Cxy,offset+16);
         storea(c6,Cxy,offset+20);
         storea(c7,Cxy,offset+24);
         storea(c8,Cxy,offset+28);
      }
   }
   printf("%lu %d \n",t3,n);
}
*/

void  accumTest4_7_13(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
  // double a[innerWid*innerWid];
for(int i = 0; i < 4; i++){
   for(int diag = 0; diag < n; diag += 16){
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      vtf dF1 = loada(dF,diag);
      vtf dF2 = loada(dF,diag+4);
      vtf dF3 = loada(dF,diag+8);
      vtf dF4 = loada(dF,diag+12); 
      vtf dX1 = loada(dX,diag);
      vtf dX2 = loada(dX,diag+4);
      vtf dX3 = loada(dX,diag+8);
      vtf dX4 = loada(dX,diag+12);

      for(int offset = 0; offset < diag; offset+=16){
         t3 += 64;
         vtf dFY = bcast(dF,offset);
      //   vtf sY = bcast(s,offset);
         c1 = mult_add(dFY,dX1,c1);
         c2 = mult_add(dFY,dX2,c2); 
         c3 = mult_add(dFY,dX3,c3);
         c4 = mult_add(dFY,dX4,c4);


         vtf dXY = bcast(dX,offset);
         c1 = mult_add(dXY,dF1,c1);
         c2 = mult_add(dXY,dF2,c2);
         c3 = mult_add(dXY,dF3,c3);
         c4 = mult_add(dXY,dF4,c4);

      /*   storea(c1,output,diag+offset);
         storea(c2,output,diag+offset+4);
         storea(c3,output,diag+offset+8);
         storea(c4,output,diag+offset+12);
     */    
              
         dX1 = loada(dX,diag+4);
         dFY = bcast(dF,offset+1);
         c1 = mult_add(dFY,dX2,c1);
         c2 = mult_add(dFY,dX3,c2);
         c3 = mult_add(dFY,dX4,c3);
         c4 = mult_add(dFY,dX1,c4);
 
         dXY = bcast(dX,offset+1);
         dF1 = loada(dF,diag+4);
         c1  = mult_add(dXY,dF2,c1);
         c2  = mult_add(dXY,dF3,c2);
         c3  = mult_add(dXY,dF4,c3);
         c4  = mult_add(dXY,dF1,c4);

         storea(c1,output,diag+offset+20);
         storea(c2,output,diag+offset+24);
         storea(c3,output,diag+offset+28);
         storea(c4,output,diag+offset+32);
           
        
         dFY = bcast(dF,offset+2);
         dX2 = loada(dX,diag+8);
         c1  = mult_add(dFY,dX3,c1);
         c2  = mult_add(dFY,dX4,c2);
         c3  = mult_add(dFY,dX1,c3);
         c4  = mult_add(dFY,dX2,c4);
         
         dXY = bcast(dX,offset+2);
         dF2 = loada(dF,diag+8);
         c1  = mult_add(dXY,dF3,c1);
         c2  = mult_add(dXY,dF4,c2);
         c3  = mult_add(dXY,dF1,c3);
         c4  = mult_add(dXY,dF2,c4);
   
            
         storea(c1,output,diag+offset+40);
         storea(c2,output,diag+offset+44);
         storea(c3,output,diag+offset+48);
         storea(c4,output,diag+offset+52);
           
      }
      storeu(c1,Cxy,diag);
      storeu(c2,Cxy,diag+4);
      storeu(c3,Cxy,diag+8);
      storeu(c4,Cxy,diag+12);
      
      }
   }

printf("%lu %d \n",t3,n);
}






void  accumTest4_7_3(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*16];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset+diag < n; offset += 4){
         vtf dX1 = bcast(dX,offset);
         vtf dX2 = bcast(dX,offset+1);
         vtf dX3 = bcast(dX,offset+2);
         vtf dX4 = bcast(dX,offset+3);
         vtf dF1 = bcast(dF,offset);
         vtf dF2 = bcast(dF,offset+1);
         vtf dF3 = bcast(dF,offset+2);
         vtf dF4 = bcast(dF,offset+3);
         vtf mp1 = bcast(output,offset);
         vtf mp2 = bcast(output,offset+1);
         vtf mp3 = bcast(output,offset+2);
         vtf mp4 = bcast(output,offset+3);
         vti mpI1 =bcast(outputI[offset]);
         vti mpI2 =bcast(outputI[offset+1]);
         vti mpI3 =bcast(outputI[offset+2]);
         vti mpI4 =bcast(outputI[offset+3]);
         for(int subdiag = diag+offset; subdiag <  diag+offset+innerWid; subdiag += 4){

            vti ind = _mm256_set_epi64x(subdiag,subdiag+1,subdiag+2,subdiag+3);
            vti inc = bcast(1);
            vtf dX12 = loada(dX,subdiag);
            vtf dX22 = preshuffle(dX12,loada(dX,subdiag+4));
            vtf dX32 = shift2(dX12,dX22);
                dX22 = shift1(dX12,dX22);
            vtf dX42 = loadu(dX,subdiag+3);

            vtf dF12 = loada(dF,subdiag);
            vtf dF22 = preshuffle(dF12,loada(dF,subdiag+4));
            vtf dF32 = shift2(dF12,dF22);
                dF22 = shift1(dF12,dF22);
            vtf dF42 = loadu(dF,subdiag+3);

            vtf s12  = loada(s,subdiag);
            vtf s22  = preshuffle(s12,loada(s,subdiag+4));
            vtf s32  = shift2(s12,s22);
                s22  = shift1(s12,s22);
            vtf s42  = loadu(s,subdiag+3);

            vtf c1 = mult_add(dF12,dX1,loada(Cxy,subdiag));
                c1 = mult_add(dX12,dF1,c1);
            vtf sx = mult(c1,s12);
            vtf cmp2 = cmpgtr(mp1,c1);
                mp1 = blend(mp1,sx,cmp2);
                mpI1= blend(mpI1,ind,cmp2);
                 
                c1 = mult_add(dF22,dX2,c1);
                c1 = mult_add(dF2,dX22,c1);
                sx = mult(c1,s22);
                cmp2 = cmpgtr(mp2,sx);
                mp2 = blend(mp2,sx,cmp2); 
                mpI2= blend(mpI2,ind+1,cmp2); 

                c1   = mult_add(dF32,dX3,c1);
                c1   = mult_add(dF3,dX32,c1);
                sx   = mult(c1,s32);
                cmp2 = cmpgtr(mp3,sx);
                mp3  = blend(mp3,sx,cmp2);
                mpI3 = blend(mpI3,ind+2,cmp2);

                c1 = mult_add(dF42,dX4,c1);
                c1 = mult_add(dF4,dX42,c1);
                sx = mult(c1,s42);
                cmp2 = cmpgtr(mp4,sx);
                mp4 = blend(mp4,sx,cmp2);
                mpI4 = blend(mpI4,ind+3,cmp2);
                t3 += 16;
                storea(c1,Cxy,subdiag);
         }
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}

static inline double hmax(__m256d a){
   __m256d b = _mm256_permute2f128_pd(a,a,1);
   __m256d c = max(a,b);
   c = max(c,_mm256_permute_pd(c,5));
   return c[0]; // adjust later
}

/*static inline double hmax(__m256d a){
   __m128d b = _mm_max_pd(_mm256_extractf128_pd(a,1);
  // a  = _mm256_extractf128_pd(a,1);
   // use zero upper
   
}*/


// also reasonable quite possibly the fastest  
void  accumTest4_7_11(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset+diag < n; offset += 4){
    int col = diag+offset;
    vtf dX1 = loada(dX,col);
    vtf dX2 = loadu(dX,col+1);
    vtf dX3 = loadu(dX,col+2);
    vtf dX4 = loadu(dX,col+3);
    vtf dF1 = loada(dF,col);
    vtf dF2 = loadu(dF,col+1);
    vtf dF3 = loadu(dF,col+2);
    vtf dF4 = loadu(dF,col+3);
    vtf s1  = loada(s,col);
    vtf s2  = loadu(s,col+1);
    vtf s3  = loadu(s,col+2);
    vtf s4  = loadu(s,col+3);
    vtf mp1 = loada(output,col);
    vtf mp2 = loadu(output,col+1);
    vtf mp3 = loadu(output,col+2);
    vtf mp4 = loadu(output,col+3);
    vti mpI1 = loada(outputI,col);
    vti mpI2 = loadu(outputI,col+1);
    vti mpI3 = loadu(outputI,col+2);
    vti mpI4 = loadu(outputI,col+3);
    for(int subdiag = 0; subdiag < innerWid; subdiag += 4){
         vtf dX12  = bcast(dX,offset+subdiag);
         vtf dX22  = bcast(dX,offset+subdiag+1);
         vtf dX32  = bcast(dX,offset+subdiag+2);
         vtf dX42  = bcast(dX,offset+subdiag+3);
         vtf dF12  = bcast(dF,offset+subdiag);
         vtf dF22  = bcast(dF,offset+subdiag+1);
         vtf dF32  = bcast(dF,offset+subdiag+2);
         vtf dF42  = bcast(dF,offset+subdiag+3);
         vtf s12   = bcast(s,offset+subdiag);
         vtf s22   = bcast(s,offset+subdiag+1);
         vtf s32   = bcast(s,offset+subdiag+2);
         vtf s42   = bcast(s,offset+subdiag+3);
     /*    vtf mp12  = bcast(output,offset+subdiag);
         vtf mp22  = bcast(output,offset+subdiag+1);
         vtf mp32  = bcast(output,offset+subdiag+2);
         vtf mp42  = bcast(output,offset+subdiag+3);
*/
       // definitely a bit faster with reversed order, which means that we have to do the prior iteration of maxes on the current pass.
       // This can of course be tuned  a bit further for the broadcast dimension.

              // vtf c1 = loada(Cxy,subdiag);
              // c1    = mult_add(dF12,dX1,c1);
                vtf c1   = mult_add(dF12,dX1,loada(Cxy,subdiag));
                c1   = mult_add(dX12,dF1,c1);
                vtf s1 = mult(c1,s12);
                vtf cmp1 = cmpgtr(mp1,c1);
                mp1 = blend(mp1,s1,cmp1);
                mp1 = max(mp1,s1);
                mpI1 = blend(mpI1,set(subdiag+diag,subdiag+diag+1,subdiag+diag+2,subdiag+diag+3),cmp1);
                output[subdiag+offset] = hmax(s1);
               
                c1   = mult_add(dF22,dX2,c1);
                c1   = mult_add(dF2,dX22,c1);
                s1  = mult(c1,s22);
                
                mp2 = max(mp2,s1);
                cmp1 = cmpgtr(mp2,s1);
                mp2  = blend(mp2,s1,cmp1);
                mpI2 = blend(mpI2,set(subdiag+diag+1,subdiag+diag+2,subdiag+diag+3,subdiag+diag+4),cmp1);
                output[subdiag+offset+1] = hmax(s1);

                c1   = mult_add(dF32,dX3,c1);
                c1   = mult_add(dF3,dX32,c1);
                s1   = mult(c1,s32);
                mp3  = max(mp3,mult(c1,s32));
                cmp1  = cmpgtr(mp3,s1);
                mp3   = blend(mp3,s1,cmp1);
                mpI3  = blend(mpI3,set(subdiag+diag+2,subdiag+diag+3,subdiag+diag+4,subdiag+diag+5),cmp1);
                output[subdiag+offset+2] = hmax(s1);


                c1   = mult_add(dF42,dX4,c1);
                c1   = mult_add(dF4,dX42,c1);
                s1   = mult(c1,s42);
                mp4  = max(mp4,mult(c1,s42));
                cmp1 = cmpgtr(mp4,s1);
                mp4  = blend(mp4,s1,cmp1);
                mpI4 = blend(mpI4,set(subdiag+diag+3,subdiag+diag+4,subdiag+diag+5,subdiag+diag+6),cmp1);
                output[subdiag+offset+3] = hmax(s1); 
               
                t3 += 16;                
                storea(c1,Cxy,subdiag);

         }
         
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
         storea(mpI1,outputI,diag+offset);
         storea(mpI2,outputI,diag+offset+4);
         storea(mpI3,outputI,diag+offset+8);
         storea(mpI4,outputI,diag+offset+12);

      }
   }
   printf("%lu %d \n",t3,n);
}

void  accumTest4_7_16(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*4];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset+diag < n; offset += 4){
    int col = diag+offset;
    vtf mp1 = loada(output,col);
    vtf mp2 = loada(output,col+4);
    vtf mp3 = loada(output,col+8);
    vtf mp4 = loada(output,col+12);
    vtf s1  = loada(s,col);
    vtf s2  = loada(s,col+4);
    vtf s3  = loada(s,col+8);
    vtf s4  = loada(s,col+12);
    vtf X1  = loada(dX,col);
    vtf X2  = loada(dX,col+4);
    vtf X3  = loada(dX,col+8);
    vtf X4  = loada(dX,col+12);
    vtf F1  = loada(dF,col);
    vtf F2  = loada(dF,col+4);
    vtf F3  = loada(dF,col+8);
    vtf F4  = loada(dF,col+12);

    for(int subdiag = 0; subdiag < innerWid; subdiag += 4){
            vtf c1   = mult_add(F1,bcast(dX,offset),loada(Cxy,col));
                c1   = mult_add(X1,bcast(dF,offset),  c1);
                mp1 = max(mp1,mult(bcast(s,offset),mult(c1,s1)));
                c1   = mult_add(F2,bcast(dX,offset+4),c1);
                c1   = mult_add(X2,bcast(dF,offset+4), c1);
                mp2   = max(mp2,mult(bcast(s,offset+4),mult(c1,s2)));
                c1   = mult_add(F3,bcast(dX,offset+8),c1);
                c1   = mult_add(X3,bcast(dF,offset+8), c1);
                mp3 = max(mp3,mult(bcast(s,offset+8),mult(c1,s3)));
                c1  = mult_add(F4,bcast(dX,offset+12),c1);
                c1  = mult_add(X4,bcast(dF,offset+12), c1);
                s1  = max(mp4,mult(bcast(s,offset+12),mult(c1,s4)));
                storea(c1,Cxy,subdiag);
         }
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}




// this one might be a winner
void  accumTest4_7_15(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*4];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset+diag < n; offset += 4){
    int col = diag+offset;
    vtf mp1 = loada(output,col);
    vtf mp2 = loada(output,col+4);
    vtf mp3 = loada(output,col+8);
    vtf mp4 = loada(output,col+12);
    
    for(int subdiag = 0; subdiag < innerWid; subdiag += 16){
            col = subdiag+diag;
            vtf c1   = mult_add(loada(dF,col),bcast(dX,offset)  ,loada(Cxy,col));
            vtf c2   = mult_add(loada(dF,col),bcast(dX,offset+4),loada(Cxy,col+4));
            vtf c3   = mult_add(loada(dF,col),bcast(dX,offset+8),loada(Cxy,col+8));
            vtf c4   = mult_add(loada(dF,col),bcast(dX,offset+12),loada(Cxy,col+12));

                c1   = mult_add(loada(dX,col),bcast(dF,offset),  c1);
                c2   = mult_add(loada(dX,col),bcast(dF,offset+4),c2);
                c3   = mult_add(loada(dX,col),bcast(dF,offset+8),c3);
                c4   = mult_add(loada(dX,col),bcast(dF,offset+12),c4);
     
                vtf s1 = mult(loada(s,col),mult(c1,bcast(s,offset)));
                vtf s2 = mult(loada(s,col),mult(c2,bcast(s,offset+4)));
                vtf s3 = mult(loada(s,col),mult(c3,bcast(s,offset+8)));
                vtf s4 = mult(loada(s,col),mult(c4,bcast(s,offset+12)));

                mp1 = max(mp1,max(max(s1,s2),max(s3,s4)));
              
                col += 4;
                c1   = mult_add(loada(dF,col),bcast(dX,offset+1),c1);
                c2   = mult_add(loada(dF,col),bcast(dX,offset+5),c2);
                c3   = mult_add(loada(dF,col),bcast(dX,offset+9),c3);
                c4   = mult_add(loada(dF,col),bcast(dX,offset+13),c4);
            
                c1   = mult_add(loada(dX,col),bcast(dF,offset+1), c1);
                c2   = mult_add(loada(dX,col),bcast(dF,offset+5), c2);
                c3   = mult_add(loada(dX,col),bcast(dF,offset+9),c3);
                c4   = mult_add(loada(dX,col),bcast(dF,offset+13),c4);
     
                s1 = mult(loada(s,col),mult(c1,bcast(s,offset+1)));
                s2 = mult(loada(s,col),mult(c2,bcast(s,offset+5)));
                s3 = mult(loada(s,col),mult(c3,bcast(s,offset+9)));
                s4 = mult(loada(s,col),mult(c4,bcast(s,offset+13)));

                mp2 = max(mp2,max(max(s1,s2),max(s3,s4)));

                col += 4;
                c1   = mult_add(loada(dF,col),bcast(dX,offset+2),c1);
                c2   = mult_add(loada(dF,col),bcast(dX,offset+6),c2);
                c3   = mult_add(loada(dF,col),bcast(dX,offset+10),c3);
                c4   = mult_add(loada(dF,col),bcast(dX,offset+14),c4);
            
                c1   = mult_add(loada(dX,col),bcast(dF,offset+2), c1);
                c2   = mult_add(loada(dX,col),bcast(dF,offset+6), c2);
                c3   = mult_add(loada(dX,col),bcast(dF,offset+10),c3);
                c4   = mult_add(loada(dX,col),bcast(dF,offset+14),c4);
     
                s1 = mult(loada(s,col),mult(c1,bcast(s,offset+2)));
                s2 = mult(loada(s,col),mult(c2,bcast(s,offset+6)));
                s3 = mult(loada(s,col),mult(c3,bcast(s,offset+10)));
                s4 = mult(loada(s,col),mult(c4,bcast(s,offset+14)));

                mp3 = max(mp3,max(max(s1,s2),max(s3,s4)));

                col += 4;
                c1   = mult_add(loada(dF,col),bcast(dX,offset+3),c1);
                c2   = mult_add(loada(dF,col),bcast(dX,offset+7),c2);
                c3   = mult_add(loada(dF,col),bcast(dX,offset+11),c3);
                c4   = mult_add(loada(dF,col),bcast(dX,offset+15),c4);
            
                c1   = mult_add(loada(dX,col),bcast(dF,offset+3), c1);
                c2   = mult_add(loada(dX,col),bcast(dF,offset+7),c2);
                c3   = mult_add(loada(dX,col),bcast(dF,offset+11),c3);
                c4   = mult_add(loada(dX,col),bcast(dF,offset+15),c4);
     
                s1 = mult(loada(s,col),mult(c1,bcast(s,offset+3)));
                s2 = mult(loada(s,col),mult(c2,bcast(s,offset+7)));
                s3 = mult(loada(s,col),mult(c3,bcast(s,offset+11)));
                s4 = mult(loada(s,col),mult(c4,bcast(s,offset+15)));

                mp4 = max(mp4,max(max(s1,s2),max(s3,s4)));


                t3 += 64;                

                storea(c1,Cxy,subdiag);
         }
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}



void  accumTest4_7_2(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*16];
   for(int diag = 0; diag < n-4*innerWid; diag += innerWid){
      for(int offset = 0; offset+diag+innerWid < n-4*innerWid; offset += 4){
         vtf dX1  = bcast(dX,offset);
         vtf dX2  = bcast(dX,offset+1);
         vtf dX3  = bcast(dX,offset+2);
         vtf dX4  = bcast(dX,offset+3);
         vtf dF1  = bcast(dF,offset);
         vtf dF2  = bcast(dF,offset+1);
         vtf dF3  = bcast(dF,offset+2);
         vtf dF4  = bcast(dF,offset+3);
         vtf s1   = bcast(s,offset);
         vtf s2   = bcast(s,offset+1);
         vtf s3   = bcast(s,offset+2);
         vtf s4   = bcast(s,offset+3);
         vtf mp1  = bcast(output,offset);
         vtf mp2  = bcast(output,offset+1);
         vtf mp3  = bcast(output,offset+2);
         vtf mp4  = bcast(output,offset+3);
         vti mpI1 = bcast(outputI[offset]);
         vti mpI2 = bcast(outputI[offset+1]);
         vti mpI3 = bcast(outputI[offset+2]);
         vti mpI4 = bcast(outputI[offset+3]);
         for(int subdiag = diag+offset; subdiag < diag+offset+innerWid; subdiag += 4){
            
            vti ind  = _mm256_set_epi64x(subdiag,subdiag+1,subdiag+2,subdiag+3);
            vti inc  = bcast(1);
            
            vtf dX12 = loada(dX,subdiag);
            vtf dX22 = preshuffle(dX12,loada(dX,subdiag+4));
            vtf dX32 = shift2(dX12,dX22);
                dX22 = shift1(dX12,dX22);
            
            vtf dF12 = loada(dF,subdiag);
            vtf dF22 = preshuffle(dF12,loada(dF,subdiag+4));
            vtf dF32 = shift2(dF12,dF22);
                dF22 = shift1(dF12,dF22);

            vtf s12  = loada(s,subdiag);
            vtf s22  = preshuffle(s12,loada(s,subdiag+4));
            vtf s32  = mult(s3,shift2(s12,s22));
                s22  = mult(s2,shift1(s12,s22));
                s12  = mult(s2,s12);
            vtf s42  = mult(s4,loadu(s,subdiag+3));

            vtf c1   = mult_add(dF12,dX1,loada(Cxy,subdiag));
                c1   = mult_add(dX12,dF1,c1);
                s12  = mult(c1,s12);
//                storea(s12,a,subdiag-diag-offset);
            vtf cmp2 = cmpgtr(mp1,c1);
                mp1  = blend(mp1,s12,cmp2);
                mpI1 = blend(mpI1,ind,cmp2);
                 
                c1   = mult_add(dF22,dX2,c1);
                c1   = mult_add(dF2,dX22,c1);
                s22  = mult(c1,s22);
      //          storea(s22,a,subdiag-diag-offset+4);
                cmp2 = cmpgtr(mp2,s22);
                mp2  = blend(mp2,s22,cmp2); 
                mpI2 = blend(mpI2,ind+1,cmp2); 

                c1   = mult_add(dF32,dX3,c1);
                c1   = mult_add(dF3,dX32,c1);
                s32  = mult(c1,s32);
    //            storea(s32,a,subdiag-diag-offset+8);
                cmp2 = cmpgtr(mp3,s32);
                mp3  = blend(mp3,s32,cmp2);
                mpI3 = blend(mpI3,ind+2,cmp2);

                c1   = mult_add(loadu(dF,subdiag+3),dX4,c1);
                c1   = mult_add(dF4,loadu(dX,subdiag+3),c1);
                s42  = mult(c1,s42);
  //              storea(s42,a,subdiag-diag-offset+12);
                cmp2 = cmpgtr(mp4,s42);
                mp4  = blend(mp4,s42,cmp2);
                mpI4 = blend(mpI4,ind+3,cmp2);
                t3 += 16;                

                storea(c1,Cxy,subdiag);
         }
/*
       for(int i = 0; i < innerWid; i+=4){
            vti e1 = set(1,2,3,4);
            vti e2 = set(5,6,7,8);
            vti e3 = set(6,7,8,9);
            vti e4 = set(10,11,12,13);
            vtf f1 = loadu(output,diag);
            vtf d1 = loada(a,i);
            vtf d2 = loada(a,4*i+4);
            vtf d3 = loada(a,4*i+8);
            vtf d4 = loada(a,4*i+12);
            vtf cmp1 = cmpgtr(d1,d2);
            vtf cmp2 = cmpgtr(d3,d4);
                d1   = blend(d1,d2,cmp1);
                e1   = blend(e1,e2,cmp1);
                d2   = blend(d3,d4,cmp2);
                e2   = blend(e3,e4,cmp2);
                cmp1 = cmpgtr(d1,d2);
                d1   = blend(d1,d2,cmp1);
                e1   = blend(e1,e2,cmp2);
                storea(d1,output,diag+offset);
                storea(d2,output,diag+offset+4);
                storea(e1+bcast(offset),outputI,diag+offset);
                storea(e2+bcast(offset),outputI,diag+offset+4);
         }
*/
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}

