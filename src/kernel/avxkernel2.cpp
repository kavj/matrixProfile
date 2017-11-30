#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 256 

using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;

/*

void  accumTest4_7_13(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double* a = (double*)outputI;

//   double a[innerWid*40];
for(int i = 0; i < 4; i++){
for(int superoff = 0; superoff < n; superoff+=innerWid){
  for(int superdiag = 0; superdiag < superoff; superdiag += innerWid){
   for(int diag = superdiag; diag < superdiag+innerWid; diag += 40){
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      vtf c5 = loada(Cxy,diag+16);
      vtf c6 = loada(Cxy,diag+20);
      vtf c7 = loada(Cxy,diag+24);
      vtf c8 = loada(Cxy,diag+28);
      vtf c9 = loada(Cxy,diag+32);
      vtf c10= loada(Cxy,diag+36);
      for(int offset = superoff; offset < superoff+innerWid; offset+=4){
         t3 += 40;
         vtf dXY = bcast(dX,offset);
         vtf dFY = bcast(dF,offset);
         vtf sY = bcast(s,offset);
 
         c1 = mult_add(dXY,loada(dF,offset+diag),c1);
         c2 = mult_add(dXY,loada(dF,offset+diag)+4,c2);
         c3 = mult_add(dXY,loada(dF,offset+diag+4),c3);
         c4 = mult_add(dXY,loada(dF,offset+diag+8),c4);
         c5 = mult_add(dXY,loada(dF,offset+diag+12),c5);
         c6 = mult_add(dXY,loada(dF,offset+diag+16),c6);
         c7 = mult_add(dXY,loada(dF,offset+diag+20),c7);
         c8 = mult_add(dXY,loada(dF,offset+diag+24),c8);
         c9 = mult_add(dXY,loada(dF,offset+diag+28),c9);
         c10= mult_add(dXY,loada(dF,offset+diag+32),c10);
         
         vtf s1 = mult(c1,loada(s,offset+diag));
         vtf s2 = mult(c2,loada(s,offset+diag+4));
         vtf s3 = mult(c3,loada(s,offset+diag+8));
         vtf s4 = mult(c4,loada(s,offset+diag+12));
         vtf s5 = mult(c5,loada(s,offset+diag+16));
         vtf s6 = mult(c6,loada(s,offset+diag+20));
         vtf s7 = mult(c7,loada(s,offset+diag+24));
         vtf s8 = mult(c8,loada(s,offset+diag+28));
         vtf s9 = mult(c9,loada(s,offset+diag+32));
         vtf s10= mult(c10,loada(s,offset+diag+36));

         vtf mp1 = max(s1,bcast(output,offset));
         mp1 = max(mp1,s2);
         vtf mp2 = max(s3,s4);
         vtf mp3 = max(s4,s5);
         vtf mp4 = max(s6,s7);
         vtf mp5 = max(s8,s9);
         mp1 = max(mp1,s10);
         mp2 = max(mp2,mp3);
         mp4 = max(mp4,mp5);
         mp1 = max(mp1,mp2);
         mp1 = max(mp1,mp4); 

         storea(max(s1,loada(a,40*(offset-superoff))),a,40*(offset-superoff));
         storea(max(s2,loada(a,40*(offset-superoff+1))),a,40*(offset-superoff+1));
         storea(max(s3,loada(a,40*(offset-superoff+2))),a,40*(offset-superoff+2));
         storea(max(s4,loada(a,40*(offset=superoff+3))),a,40*(offset-superoff+3));
      /*   storea(max(s5,loada(a,40*(offset-superoff+4))),a,40*(offset-superoff+4));
         storea(max(s6,loada(a,40*(offset-superoff+5))),a,40*(offset-superoff+5));
         storea(max(s7,loada(a,40*(offset-superoff+6))),a,40*(offset-superoff+6));
         storea(max(s8,loada(a,40*(offset=superoff+7))),a,40*(offset-superoff+7));
         storea(max(s9,loada(a,40*(offset-superoff+8))),a,40*(offset-superoff+8));
         storea(max(s10,loada(a,40*(offset-superoff+9))),a,40*(offset-superoff+9));

         storea(mp1,output,diag+offset);
      }
      storea(c1,Cxy,diag);
      storea(c2,Cxy,diag+4);
      storea(c3,Cxy,diag+8);
      storea(c4,Cxy,diag+12);
      storea(c5,Cxy,diag+16);
      storea(c6,Cxy,diag+20);
      storea(c7,Cxy,diag+24);
      storea(c8,Cxy,diag+28);
      storea(c9,Cxy,diag+32);
      storea(c10,Cxy,diag+36);
      
      }
 }
}
*/

void accumTest4_7_8(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double* a = (double*)outputI;
//   double a[innerWid*40];
for(int i = 0; i < 4; i++){
   for(int diag = 0; diag < n; diag += 40){
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      vtf c5 = loada(Cxy,diag+16);
      vtf c6 = loada(Cxy,diag+20);
      vtf c7 = loada(Cxy,diag+24);
      vtf c8 = loada(Cxy,diag+28);
      vtf c9 = loada(Cxy,diag+32);
      vtf c10= loada(Cxy,diag+36);
      for(int offset = 0; offset < diag; offset+=4){
         t3 += 40;
         vtf dXY = bcast(dX,offset);
         vtf dFY = bcast(dF,offset);
         vtf sY = bcast(s,offset);
 
         c1 = mult_add(dXY,loada(dF,offset+diag),c1);
         c2 = mult_add(dXY,loada(dF,offset+diag)+4,c2);
         c3 = mult_add(dXY,loada(dF,offset+diag+4),c3);
         c4 = mult_add(dXY,loada(dF,offset+diag+8),c4);
         c5 = mult_add(dXY,loada(dF,offset+diag+12),c5);
         c6 = mult_add(dXY,loada(dF,offset+diag+16),c6);
         c7 = mult_add(dXY,loada(dF,offset+diag+20),c7);
         c8 = mult_add(dXY,loada(dF,offset+diag+24),c8);
         c9 = mult_add(dXY,loada(dF,offset+diag+28),c9);
         c10= mult_add(dXY,loada(dF,offset+diag+32),c10);
         
         vtf s1 = mult(c1,loada(s,offset+diag));
         vtf s2 = mult(c2,loada(s,offset+diag+4));
         vtf s3 = mult(c3,loada(s,offset+diag+8));
         vtf s4 = mult(c4,loada(s,offset+diag+12));
         vtf s5 = mult(c5,loada(s,offset+diag+16));
         vtf s6 = mult(c6,loada(s,offset+diag+20));
         vtf s7 = mult(c7,loada(s,offset+diag+24));
         vtf s8 = mult(c8,loada(s,offset+diag+28));
         vtf s9 = mult(c9,loada(s,offset+diag+32));
         vtf s10= mult(c10,loada(s,offset+diag+36));

         vtf mp1 = max(s1,bcast(output,offset));
         mp1 = max(mp1,s2);
         vtf mp2 = max(s3,s4);
         vtf mp3 = max(s4,s5);
         vtf mp4 = max(s6,s7);
         vtf mp5 = max(s8,s9);
         mp1 = max(mp1,s10);
         mp2 = max(mp2,mp3);
         mp4 = max(mp4,mp5);
         mp1 = max(mp1,mp2);
         mp1 = max(mp1,mp4); 

         storea(max(s1,loada(a,40*(offset))),a,40*(offset));
         storea(max(s2,loada(a,40*(offset+1))),a,40*(offset+1));
         storea(max(s3,loada(a,40*(offset+2))),a,40*(offset+2));
         storea(max(s4,loada(a,40*(offset+3))),a,40*(offset+3));
         storea(max(s5,loada(a,40*(offset+4))),a,40*(offset+4));
         storea(max(s6,loada(a,40*(offset+5))),a,40*(offset+5));
         storea(max(s7,loada(a,40*(offset+6))),a,40*(offset+6));
         storea(max(s8,loada(a,40*(offset+7))),a,40*(offset+7));
         storea(max(s9,loada(a,40*(offset+8))),a,40*(offset+8));
         storea(max(s10,loada(a,40*(offset+9))),a,40*(offset+9));

         storea(mp1,output,diag+offset);
      }
      storea(c1,Cxy,diag);
      storea(c2,Cxy,diag+4);
      storea(c3,Cxy,diag+8);
      storea(c4,Cxy,diag+12);
      storea(c5,Cxy,diag+16);
      storea(c6,Cxy,diag+20);
      storea(c7,Cxy,diag+24);
      storea(c8,Cxy,diag+28);
      storea(c9,Cxy,diag+32);
      storea(c10,Cxy,diag+36);
      
      }
   }

printf("%lu %d \n",t3,n);
}

void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
for(int i = 0; i < 4; i++){
   for(int diag = 0; diag < n; diag += 40){
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      vtf c5 = loada(Cxy,diag+16);
      vtf c6 = loada(Cxy,diag+20);
      vtf c7 = loada(Cxy,diag+24);
      vtf c8 = loada(Cxy,diag+28);
      vtf c9 = loada(Cxy,diag+32);
      vtf c10= loada(Cxy,diag+36);
      for(int offset = 0; offset < diag; offset+=4){
         t3 += 32;
         vtf dXY = bcast(dX,offset);
         vtf dFY = bcast(dF,offset);
         vtf sY = bcast(s,offset);
 
         c1 = mult_add(dXY,loada(dF,offset+diag),c1);
         c2 = mult_add(dXY,loada(dF,offset+diag)+4,c2);
         c3 = mult_add(dXY,loada(dF,offset+diag+4),c3);
         c4 = mult_add(dXY,loada(dF,offset+diag+8),c4);
         c5 = mult_add(dXY,loada(dF,offset+diag+12),c5);
         c6 = mult_add(dXY,loada(dF,offset+diag+16),c6);
         c7 = mult_add(dXY,loada(dF,offset+diag+20),c7);
         c8 = mult_add(dXY,loada(dF,offset+diag+24),c8);
         c9 = mult_add(dXY,loada(dF,offset+diag+28),c9);
         c10= mult_add(dXY,loada(dF,offset+diag+32),c10);
         
         c1 = mult_add(dFY,loada(dX,offset+diag),c1);
         c2 = mult_add(dFY,loada(dX,offset+diag)+4,c2);
         c3 = mult_add(dFY,loada(dX,offset+diag+4),c3);
         c4 = mult_add(dFY,loada(dX,offset+diag+8),c4);
         c5 = mult_add(dFY,loada(dX,offset+diag+12),c5);
         c6 = mult_add(dFY,loada(dX,offset+diag+16),c6);
         c7 = mult_add(dFY,loada(dX,offset+diag+20),c7);
         c8 = mult_add(dFY,loada(dX,offset+diag+24),c8);
         c9 = mult_add(dFY,loada(dX,offset+diag+28),c9);
         c10= mult_add(dFY,loada(dX,offset+diag+32),c10);
         
          
         vtf s1 = mult(sY,mult(c1,loada(s,offset+diag)));
         vtf s2 = mult(sY,mult(c2,loada(s,offset+diag+4)));
         vtf s3 = mult(sY,mult(c3,loada(s,offset+diag+8)));
         vtf s4 = mult(sY,mult(c4,loada(s,offset+diag+12)));
         vtf s5 = mult(sY,mult(c5,loada(s,offset+diag+16)));
         vtf s6 = mult(sY,mult(c6,loada(s,offset+diag+20)));
         vtf s7 = mult(sY,mult(c7,loada(s,offset+diag+24)));
         vtf s8 = mult(sY,mult(c8,loada(s,offset+diag+28)));
         vtf s9 = mult(sY,mult(c9,loada(s,offset+diag+32)));
         vtf s10= mult(sY,mult(c10,loada(s,offset+diag+36)));

         s1 = max(s1,bcast(output,offset));
         s2 = max(s2,s3);
         s4 = max(s4,s5);
         s6 = max(s6,s7);
         s8 = max(s8,s9);
         s1 = max(s1,s10);
         s1 = max(s1,s2);
         s1 = max(s1,s4);
         s1 = max(s1,s6);
         s1 = max(s1,s8);
         
         storea(mult(sY,s1),output,offset);
         storea(max(loadu(output,offset+diag),s1),output,offset+diag);
         storea(max(loadu(output,offset+diag),s2),output,offset+diag+4);
         storea(max(loadu(output,offset+diag),s3),output,offset+diag+8);
         storea(max(loadu(output,offset+diag),s4),output,offset+diag+12);
         storea(max(loadu(output,offset+diag),s5),output,offset+diag+16);
         storea(max(loadu(output,offset+diag),s6),output,offset+diag+20);
         storea(max(loadu(output,offset+diag),s7),output,offset+diag+24);
         storea(max(loadu(output,offset+diag),s8),output,offset+diag+28);
         storea(max(loadu(output,offset+diag),s9),output,offset+diag+32);
         storea(max(loadu(output,offset+diag),s10),output,offset+diag+36);
         
      }
      storea(c1,Cxy,diag);
      storea(c2,Cxy,diag+4);
      storea(c3,Cxy,diag+8);
      storea(c4,Cxy,diag+12);
      storea(c5,Cxy,diag+16);
      storea(c6,Cxy,diag+20);
      storea(c7,Cxy,diag+24);
      storea(c8,Cxy,diag+28);
      storea(c9,Cxy,diag+32);
      storea(c10,Cxy,diag+36);
      
      }
   }

printf("%lu %d \n",t3,n);
}


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








void  accumTest4_7_2(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
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
         vtf s1  = bcast(s,offset);
         vtf s2  = bcast(s,offset+1);
         vtf s3  = bcast(s,offset+2);
         vtf s4  = bcast(s,offset+3);
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
            
            vtf dX12 = loada(dX,diag+offset);
            vtf dX22 = preshuffle(dX12,loada(dX,subdiag+4));
            vtf dX32 = shift2(dX12,dX22);
                dX22 = shift1(dX12,dX22);
            
            vtf dF12 = loada(dF,diag+offset);
            vtf dF22 = preshuffle(dF12,loada(dF,subdiag+4));
            vtf dF32 = shift2(dF12,dF22);
                dF22 = shift1(dF12,dF22);

            vtf s12  = loada(s,diag+offset);
            vtf s22  = preshuffle(s12,loada(s,subdiag+4));
            vtf s32  = mult(s3,shift2(s12,s22));
                s22  = mult(s2,shift1(s12,s22));
                s12  = mult(s2,s12);
            vtf s42   = mult(s4,loadu(s,subdiag+3));

            vtf c1 = mult_add(dF12,dX1,loada(Cxy,subdiag));
                c1 = mult_add(dX12,dF1,c1);
                s12 = mult(c1,s12);
            storea(s12,a,subdiag-diag-offset);
            vtf cmp2 = cmpgtr(mp1,c1);
                mp1 = blend(mp1,s12,cmp2);
                mpI1= blend(mpI1,ind,cmp2);
                 
                c1 = mult_add(dF22,dX2,c1);
                c1 = mult_add(dF2,dX22,c1);
                s22 = mult(c1,s22);
                storea(s22,a,subdiag-diag-offset+4);
                cmp2 = cmpgtr(mp2,s22);
                mp2 = blend(mp2,s22,cmp2); 
                mpI2= blend(mpI2,ind+1,cmp2); 

                c1   = mult_add(dF32,dX3,c1);
                c1   = mult_add(dF3,dX32,c1);
                s32   = mult(c1,s32);
                storea(s32,a,subdiag-diag-offset+8);
                cmp2 = cmpgtr(mp3,s32);
                mp3  = blend(mp3,s32,cmp2);
                mpI3 = blend(mpI3,ind+2,cmp2);

                c1 = mult_add(loadu(dF,subdiag+3),dX4,c1);
                c1 = mult_add(dF4,loadu(dX,subdiag+3),c1);
                s42 = mult(c1,s42);
                storea(s42,a,subdiag-diag-offset+12);
                cmp2 = cmpgtr(mp4,s42);
                mp4 = blend(mp4,s42,cmp2);
                mpI4 = blend(mpI4,ind+3,cmp2);
                t3 += 16;

                storea(c1,Cxy,subdiag);
         }
         for(int i = 0; i < innerWid; i+=4){
            vti e1 = set(1,2,3,4);
            vti e2 = set(5,6,7,8);
            vti e3 = set(6,7,8,9);
            vti e4 = set(10,11,12,13);
            vtf f1 = loada(output,diag);
            vtf f2 = loada(output,diag+4);
            vtf f3 = loada(output,diag+8);
            vtf f4 = loada(output,diag+12);
            vtf d1 = loada(a,i);
            vtf d2 = loada(a,i+4);
            vtf d3 = loada(a,i+8);
            vtf d4 = loada(a,i+12);
            vtf cmp1 = cmpgtr(f1,d1);
            vtf cmp2 = cmpgtr(f2,d2);
            vtf cmp3 = cmpgtr(f3,d3);
            vtf cmp4 = cmpgtr(f4,d4);
            d1 = blend(f1,d1,cmp1);
            e1 = blend(e1,loada(outputI,diag),cmp1);
            d2 = blend(f2,d2,cmp2);
            e2 = blend(e2,loada(outputI,diag+4),cmp2);
            d3 = blend(f3,d3,cmp3);
            e3 = blend(e3,loada(outputI,diag+8),cmp3);
            d4 = blend(f4,d4,cmp4);
            e4 = blend(e4,loada(outputI,diag+12),cmp4);
         
            storea(d1,output,diag+offset);
            storea(d2,output,diag+offset+4);
            storea(d3,output,diag+offset+8);
            storea(d4,output,diag+offset+12);
            storea(e1,outputI,diag+offset);
            storea(e2,outputI,diag+offset+4);
            storea(e3,outputI,diag+offset+8);
            storea(e4,outputI,diag+offset+12);
         }

         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}

