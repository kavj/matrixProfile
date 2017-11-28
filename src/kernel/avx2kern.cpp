#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 128 
#define alignedInner (innerWid/4-1)

using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;
/*
void  accumTest4_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   vti nums = set(0,1,2,3);
//   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset++){
            vtf s1  = loadu(s,diag+offset); 
            vtf dF1 = loadu(dF,diag+offset);
            vtf dX1 = loadu(dX,diag+offset);
            vtf mp1 = loadu(output,diag+offset);
            vti mpI1 =loadu(outputI,diag+offset);
            vti mpI2 = add(bcast(offset),nums);
         for(int subdiag = 0; subdiag <  innerWid; subdiag += 4){
            vtf c1 = mult_add(bcast(dF,offset),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset),dF1,c1);
            vtf sx = mult(s1,bcast(s,offset));
                sx = mult(sx,c1);
            vtf cmp1 = cmpgtr(sx,mp1);
                mp1 = blend(mp1,sx,cmp1);
                mpI1 = blend(mpI1,mpI2,cmp1);
                t3 += 4;
         }         
         storeu(mp1,output,diag+offset);
         storeu(mpI1,outputI,diag+offset);
      }
   }
   printf("%lu %d \n",t3,n);
}*/



void  accumTest4_8(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   double b[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += 4){
            vtf s1  = loada(s,diag+offset);
            vtf dX1 = loada(dX,diag+offset);
            vtf dF1 = loada(dF,diag+offset);
 
            vtf dF2 = preshuffle(dF1,loada(dF,diag+offset+4));
            vtf dX2 = preshuffle(dX1,loada(dX,diag+offset+4));
            vtf s2  = preshuffle(s1,loada(s,diag+offset+4));

            vtf dF3 = shift2(dF1,dF2);
            vtf dX3 = shift2(dX1,dX2);
            vtf s3  = shift2(s1,s2);
      
                dF2 = shift1(dF1,dF2);
                dX2 = shift1(dX1,dX2);
                 s2 = shift1(s1,s2);
 
            vtf s4  = loadu(s,diag+offset+3); 
            vtf dF4 = loadu(dF,diag+offset+3);
            vtf dX4 = loadu(dX,diag+offset+3);
         //   vtf mp4 = loadu(output,diag+offset+3);

         for(int subdiag = 0; subdiag <  innerWid; subdiag += 4){

            vtf c1 = mult_add(bcast(dF,offset),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset),dF1,c1);
            vtf sx = mult(s1,bcast(s,offset));
                sx = mult(sx,c1);
                storea(sx,a,4*subdiag);
                

                c1 = mult_add(bcast(dF,offset+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+1),c1);
                sx = mult(s2,bcast(s,offset+1));
                sx = mult(sx,c1);
                storea(sx,a,4*subdiag+4);
 
                c1 = mult_add(bcast(dF,offset+2),shift2(dX1,dX2),c1);
                c1 = mult_add(dF3,bcast(dX,offset+2),c1);
                sx = mult(s3,bcast(s,offset+2));
                sx = mult(sx,c1);  
                storea(sx,a,4*subdiag+8);

                c1 = mult_add(bcast(dF,offset+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+3),dF4,c1);
                sx = mult(s4,bcast(s,offset+3));
                sx = mult(sx,c1);
                storea(sx,a,4*subdiag+12);

                t3 += 16;
         }         
         // reminder: decided it's easier to do many aligned loads, then merge
         // this means we need a scalar for the row.... somehow
         
         //vti inc = set(1,2,3,4);
         vtf mp1 = loada(output,diag+offset);
         vtf mp2 = loadu(output,diag+offset+1);
         vtf mp3 = loadu(output,diag+offset+2);
         vtf mp4 = loadu(output,diag+offset+3);
         //vti mpICmp = add(inc,mp);
            for(int subdiag = 0; subdiag < innerWid; subdiag += 4){
               vtf mpCmp1 = loada(a,4*subdiag);
               vtf mpCmp2 = loada(a,4*(subdiag+1));
               vtf mpCmp3 = loada(a,4*(subdiag+2));
               vtf mpCmp4 = loada(a,4*(subdiag+3));
               vtf cmp1 = cmpgtr(mp1,mpCmp1);
                   mp1  = blend(mp1,mpCmp1,cmp1);
               vtf cmp2 = cmpgtr(mp2,mpCmp2);
                   mp2  = blend(mp2,mpCmp2,cmp2);
               vtf cmp3 = cmpgtr(mp3,mpCmp3);
                   mp3  = blend(mp3,mpCmp3,cmp3);
               vtf cmp4 = cmpgtr(mp4,mpCmp4);    
                   mp4  = blend(mp4,mpCmp4,cmp4);                     
            }
            storeu(mp1,output,offset);
            storeu(mp2,output,offset+4);
            storeu(mp3,output,offset+8);
            storeu(mp4,output,offset+12);

         }
   }
   printf("%lu %d \n",t3,n);
}



void  accumTest4_7(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += 4){
            vtf dX1 = loada(dX,diag+offset);
            vtf dF1 = loada(dF,diag+offset);
            vtf mp1 = loada(output,diag+offset);
            vti mpI1= loada(outputI,diag+offset); 
 
            vtf dF2 = preshuffle(dF1,loada(dF,diag+offset+4));
            vtf dX2 = preshuffle(dX1,loada(dX,diag+offset+4));
            vtf mp2 = preshuffle(mp1,loada(output,diag+offset+4));
            vti mpI2= preshuffle(mpI2,loada(outputI,diag+offset+4));

            vtf dF3 = shift2(dF1,dF2);
            vtf dX3 = shift2(dX1,dX2);
            vtf mp3 = shift2(mp1,mp2);
            vti mpI3 =shift2(mpI1,mpI2);
      
                dF2 = shift1(dF1,dF2);
                dX2 = shift1(dX1,dX2);
                mp2 = shift1(mp1,mp2);
                mpI2= shift1(mpI1,mpI2);        
 
            vtf dF4 = loadu(dF,diag+offset+3);
            vtf dX4 = loadu(dX,diag+offset+3);
            vtf mp4 = loadu(output,diag+offset+3);
            vti mpI4= loadu(outputI,diag+offset+3);

         for(int subdiag = 0; subdiag <  innerWid; subdiag += 4){

            vtf c1 = mult_add(bcast(dF,offset),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);
                       

                c1 = mult_add(bcast(dF,offset+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+1),c1);
                sx = mult(c1,bcast(s,offset+1));
                cmp2 = cmpgtr(mp2,sx);
                mp2 = blend(mp2,sx,cmp2);               
 
                c1 = mult_add(bcast(dF,offset+2),shift2(dX1,dX2),c1);
                c1 = mult_add(dF3,bcast(dX,offset+2),c1);
                sx = mult(c1,bcast(s,offset+2));
                cmp2 = cmpgtr(mp3,sx);
                mp3 = blend(mp3,sx,cmp2);

                c1 = mult_add(bcast(dF,offset+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+3),dF4,c1);
                sx = mult(c1,bcast(s,offset+3));
                cmp2 = cmpgtr(mp4,sx);
                mp4 = blend(mp4,sx,cmp2);

                t3 += 16;
         }
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}


