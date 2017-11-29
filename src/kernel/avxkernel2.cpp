#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 65536 

using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;



/*
void  accumTest4_11(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
  // long b[4] = {0,1,2,3};
   //const vti inc = loada(b);
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset++){
         vtf dX1 = loadu(dX,diag+offset);
         vtf dF1 = loadu(dF,diag+offset);
         vtf mp1 = loadu(output,diag+offset);
         vti mpI1= loadu(outputI,diag+offset); 

         for(int subdiag = 0; subdiag < innerWid; subdiag += 4){
            vtf c1 = mult_add(bcast(dF,offset+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+subdiag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);
                mpI1= blend(mpI1,bcast(offset+diag+subdiag),cmp2);
            
            vtf mpX = bcast(output,offset);
            vti mpXI= bcast(outputI[offset]);
             
            vtf y = _mm256_permute2f128_pd(sx,sx,1);
            vti cX = bcast(offset+diag+subdiag);
            vtf cX = cmpgtr(sx,y1);
            sx = blend(c1,y1,c3);
            
            vtf mX = _mm256_permute_pd(c1,5);
            c3 = cmpgtr(sx,mX);
            sx = blend(sx,mX,c3);
            c2 = blend(c2,mp2I,c3);  //this needs to eventually get appropriate integer permute or better yet, just generate permuted constants using broadcast + compile time stuff given that we should have a free add port
            output[offset+diag+subdiag] = _mm256_extract_epi64(c2,1);
            outputI[offset+diag+subdiag] = _mm256_extract_epi64((__m256i)sx,1);
            storea(c1,Cxy,subdiag+diag);
            t3 += 4;
         }
          
         storeu(mp1,output,diag+offset);
      }
   }
   printf("%lu %d \n",t3,n);
}
*/


/*
void accumTest4_7_4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI, int n){
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag+= 16){
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      for(int offset = 0; offset < diag; offset += 4){
            vti ind = _mm256_set_epi64x(subdiag,subdiag+1,subdiag+2,subdiag+3);
            vti inc = bcast(1);
            vtf dX1 = loada(dX,diag+offset);
            vtf dF1 = loada(dF,diag+offset);
            q23vtf mp1 = loada(output,diag+offset);
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
             
            
      }
   }
}
*/
void  accumTest4_7_8(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
  // double a[innerWid*innerWid];
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

         storea(mp1,output,offset);
      }
      storeu(c1,Cxy,diag);
      storeu(c2,Cxy,diag+4);
      storeu(c3,Cxy,diag+8);
      storeu(c4,Cxy,diag+12);
      storeu(c5,Cxy,diag+16);
      storeu(c6,Cxy,diag+20);
      storeu(c7,Cxy,diag+24);
      storeu(c8,Cxy,diag+28);
      storeu(c9,Cxy,diag+32);
      storeu(c10,Cxy,diag+36);
      
      }
   }
}

}
printf("%lu %d \n",t3,n);
}

void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
  // double a[innerWid*innerWid];
for(int i = 0; i < 4; i++){
for(int superoff = 0; superoff < n; superoff+=innerWid){
  for(int superdiag = 0; superdiag < superoff; superdiag += innerWid){
   for(int diag = superdiag; diag < superdiag+innerWid; diag += 32){
      vtf c1 = loada(Cxy,diag);
      vtf c2 = loada(Cxy,diag+4);
      vtf c3 = loada(Cxy,diag+8);
      vtf c4 = loada(Cxy,diag+12);
      vtf c5 = loada(Cxy,diag+16);
      vtf c6 = loada(Cxy,diag+20);
      vtf c7 = loada(Cxy,diag+24);
      vtf c8 = loada(Cxy,diag+28);
  //    vtf c9 = loada(Cxy,diag+32);
  //    vtf c10= loada(Cxy,diag+36);
      for(int offset = superoff; offset < superoff+innerWid; offset+=4){
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
      //   c9 = mult_add(dXY,loada(dF,offset+diag+28),c9);
      //   c10= mult_add(dXY,loada(dF,offset+diag+32),c10);
         
         c1 = mult_add(dFY,loada(dX,offset+diag),c1);
         c2 = mult_add(dFY,loada(dX,offset+diag)+4,c2);
         c3 = mult_add(dFY,loada(dX,offset+diag+4),c3);
         c4 = mult_add(dFY,loada(dX,offset+diag+8),c4);
         c5 = mult_add(dFY,loada(dX,offset+diag+12),c5);
         c6 = mult_add(dFY,loada(dX,offset+diag+16),c6);
         c7 = mult_add(dFY,loada(dX,offset+diag+20),c7);
         c8 = mult_add(dFY,loada(dX,offset+diag+24),c8);
      //   c9 = mult_add(dFY,loada(dX,offset+diag+28),c9);
      //   c10= mult_add(dFY,loada(dX,offset+diag+32),c10);
         

          
         vtf s1 = mult(c1,loada(s,offset+diag));
         vtf s2 = mult(c2,loada(s,offset+diag+4));
         vtf s3 = mult(c3,loada(s,offset+diag+8));
         vtf s4 = mult(c4,loada(s,offset+diag+12));
         vtf s5 = mult(c5,loada(s,offset+diag+16));
         vtf s6 = mult(c6,loada(s,offset+diag+20));
         vtf s7 = mult(c7,loada(s,offset+diag+24));
         vtf s8 = mult(c8,loada(s,offset+diag+28));
      //   vtf s9 = mult(c9,loada(s,offset+diag+32));
      //   vtf s10= mult(c10,loada(s,offset+diag+36));

         vtf mp1 = max(s1,bcast(output,offset));
         mp1 = max(mp1,s2);
         vtf mp2 = max(s3,s4);
         vtf mp3 = max(s4,s5);
         vtf mp4 = max(s6,s7);
       /*  vtf mp5 = max(s8,s9);
         mp1 = max(mp1,s10);
         mp2 = max(mp2,mp3);
         mp4 = max(mp4,mp5);
         mp1 = max(mp1,mp2);
         mp1 = max(mp1,mp4); 
        
         mp1 = max(mp1,mp8);
         mp2 = max(mp2,mp3);
         mp1 = max(mp1,mp4);
         mp1 = max(mp1,mp2);
       */   storeu(mult(sY,mp1),output,offset);
          storeu(max(loadu(output,offset+diag),s1),output,offset+diag);
      /*   storeu(max(loadu(output,offset+diag),s2),output,offset+diag+4);
         storeu(max(loadu(output,offset+diag),s3),output,offset+diag+8);
         storeu(max(loadu(output,offset+diag),s4),output,offset+diag+12);
         storeu(max(loadu(output,offset+diag),s5),output,offset+diag+16);
         storeu(max(loadu(output,offset+diag),s6),output,offset+diag+20);
         storeu(max(loadu(output,offset+diag),s7),output,offset+diag+24);
         storeu(max(loadu(output,offset+diag),s8),output,offset+diag+28);
        */// storeu(max(loadu(output,offset+diag),s9),output,offset+diag+32);
        // storeu(max(loadu(output,offset+diag),s10),output,offset+diag+36);
         
      }
      storeu(c1,Cxy,diag);
      storeu(c2,Cxy,diag+4);
      storeu(c3,Cxy,diag+8);
      storeu(c4,Cxy,diag+12);
      storeu(c5,Cxy,diag+16);
      storeu(c6,Cxy,diag+20);
      storeu(c7,Cxy,diag+24);
      storeu(c8,Cxy,diag+28);
     // storeu(c9,Cxy,diag+32);
     // storeu(c10,Cxy,diag+36);
      
      }
   }
}

}
printf("%lu %d \n",t3,n);
}
/*

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
            vtf mp11 = loada(output,subdiag+diag+offset);
            vtf mp12 = preshuffle(mp11,loada(output,subdiag+diag+offset));
            vtf mp13 = shift2(mp11,mp12);
            vtf mp12 = shift1(mp11,mp12);
            vtf mp14 = loadu(output,subdiag+diag+offset+3);
           
            vti ind = set(subdiag,subdiag+1,subdiag+2,subdiag+3); 
            vtf c1 = mult_add(bcast(dF,offset+diag+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+diag+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+diag+subdiag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);
                mpI1= blend(mpI1,ind,cmp2);
                cmp2 = cmpgtr(mp11,sx);
                mp11 = blend(mp11,sx,cmp2);
               
        
                c1 = mult_add(bcast(dF,offset+subdiag+diag+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+subdiag+diag+1),c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+1));

                cmp2 = cmpgtr(mp2,sx);
                mp2 = blend(mp2,sx,cmp2); 
                mpI2= blend(mpI2,ind+1,cmp2); 
                storea(sx,a,subdiag+4);

                c1   = mult_add(bcast(dF,offset+subdiag+diag+2),dX3,c1);
                c1   = mult_add(dF3,bcast(dX,offset+subdiag+diag+2),c1);
                sx   = mult(c1,bcast(s,offset+subdiag+diag+2));
                cmp2 = cmpgtr(mp3,sx);
                mp3  = blend(mp3,sx,cmp2);
                mpI3 = blend(mpI3,ind+2,cmp2);
                storea(sx,a,subdiag+8);

                c1 = mult_add(bcast(dF,offset+subdiag+diag+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+subdiag+diag+3),dF4,c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+3));
                cmp2 = cmpgtr(mp4,sx);
                mp4 = blend(mp4,sx,cmp2);
                mpI4 = blend(mpI4,ind+3,cmp2);
                storea(c1,Cxy,subdiag+diag);
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
*/



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
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8); 
         storea(mp4,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}


/*
void  accumTest4_10(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += 4){
            vtf dX1 = loada(dX,diag+offset);
            vtf dF1 = loada(dF,diag+offset);
            vtf dF2 = preshuffle(dF1,loada(dF,diag+offset+4));
            vtf dX2 = preshuffle(dX1,loada(dX,diag+offset+4));
            vtf dF3 = shift2(dF1,dF2);
            vtf dX3 = shift2(dX1,dX2);
                dF2 = shift1(dF1,dF2);
                dX2 = shift1(dX1,dX2);
            vtf dF4 = loadu(dF,diag+offset+3);
            vtf dX4 = loadu(dX,diag+offset+3);

         for(int subdiag = 0; subdiag <  innerWid; subdiag += 4){
            vtf c1 = mult_add(bcast(dF,offset+diag+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+diag+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+diag+subdiag));
           // vtf cmp2 = cmpgtr(mp1,sx);
           //     mp1 = blend(mp1,sx,cmp2);
                storea(sx,a,16*subdiag); 
                 
                c1 = mult_add(bcast(dF,offset+subdiag+diag+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+subdiag+diag+1),c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+1));
          //      cmp2 = cmpgtr(mp2,sx);
          //      mp2 = blend(mp2,sx,cmp2); 
                storea(sx,a,16*subdiag+4);              
 
                c1   = mult_add(bcast(dF,offset+subdiag+diag+2),dX3,c1);
                c1   = mult_add(dF3,bcast(dX,offset+subdiag+diag+2),c1);
                sx   = mult(c1,bcast(s,offset+subdiag+diag+2));
        //        cmp2 = cmpgtr(mp3,sx);
        //        mp3  = blend(mp3,sx,cmp2);
                storea(sx,a,16*subdiag+8);

                c1 = mult_add(bcast(dF,offset+subdiag+diag+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+subdiag+diag+3),dF4,c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+3));
          //      cmp2 = cmpgtr(mp4,sx);
          //      mp4 = blend(mp4,sx,cmp2);
                storea(sx,a,16*subdiag+12);

                t3 += 16;
         }
         vtf mp1 = loada(output,diag+offset);
         vti mpI1= loada(outputI,diag+offset); 
 
         vtf mp2 = preshuffle(mp1,loada(output,diag+offset+4));
         vti mpI2= preshuffle(mpI2,loada(outputI,diag+offset+4));
         vtf mp3 = shift2(mp1,mp2);
         vti mpI3 =shift2(mpI1,mpI2);
         mp2 = shift1(mp1,mp2);
         mpI2= shift1(mpI1,mpI2);        
         vtf mp4 = loadu(output,diag+offset+3);
         vti mpI4= loadu(outputI,diag+offset+3);
         for(int subdiag = 0; subdiag < innerWid; subdiag +=4){
            vtf mpx1 = loada(a,4*subdiag);
            vtf mpx2 = loada(a,4*(subdiag+1));
            vtf mpx3 = loada(a,4*(subdiag+2));
            vtf mpx4 = loada(a,4*(subdiag+3));
            vtf cmp1 = cmpgtr(mp1,mpx1);
            vtf cmp2 = cmpgtr(mp2,mpx2);
            vtf cmp3 = cmpgtr(mp3,mpx3);
            vtf cmp4 = cmpgtr(mp4,mpx4);
            mp1 = blend(mp1,mpx1,cmp1);
            mp2 = blend(mp2,mpx2,cmp2);
            mp3 = blend(mp3,mpx3,cmp3);
            mp4 = blend(mp4,mpx4,cmp4);
           
          //  mp1 = max(mp1,loada(a,4*(subdiag)));
          //  mp2 = max(mp2,loada(a,4*(subdiag+1)));
          //  mp3 = max(mp3,loada(a,4*(subdiag+2)));
          //  mp4 = max(mp4,loada(a,4*(subdiag+3)));
           

         }
         storea(mp1,output,diag+offset);
         storea(mp2,output,diag+offset+4);
         storea(mp3,output,diag+offset+8);
         storea(mp4,output,diag+offset+12);
*/
/*        for(int subdiag = 0; subdiag < innerWid; subdiag += 1){
            vtf c1 = loada(a,4*subdiag);
            vtf y1 = _mm256_permute2f128_pd(c1,c1,1);
            c1 = max(c1,y1);
            vtf m1 = _mm256_permute_pd(c1,5);
            c1 = max(c1,m1);
            storea(c1,output,diag+4*subdiag);
         }
*/          
     //    storea(mp1,output,diag+offset);
      //   storea(mp2,output,diag+offset+4);
      //   storea(mp3,output,diag+offset+8); 
      //   storea(mp4,output,diag+offset+12);
/*      }
   }
   printf("%lu %d \n",t3,n);
}



void  accumTest4_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
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
            vtf c1 = mult_add(bcast(dF,offset+diag+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+diag+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+diag+subdiag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);

                c1 = mult_add(bcast(dF,offset+subdiag+diag+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+subdiag+diag+1),c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+1));
                cmp2 = cmpgtr(mp2,sx);
                mp2 = blend(mp2,sx,cmp2);               
 
                c1   = mult_add(bcast(dF,offset+subdiag+diag+2),dX3,c1);
                c1   = mult_add(dF3,bcast(dX,offset+subdiag+diag+2),c1);
                sx   = mult(c1,bcast(s,offset+subdiag+diag+2));
                cmp2 = cmpgtr(mp3,sx);
                mp3  = blend(mp3,sx,cmp2);

                c1 = mult_add(bcast(dF,offset+subdiag+diag+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+subdiag+diag+3),dF4,c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+3));
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
}*/
/*
void  accumTest4_8(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += 4){
            vtf dX1 = loada(dX,diag+offset);
            vtf dF1 = loada(dF,diag+offset);
            vtf mp1 = loada(output,diag+offset);
            vti mpI1= loada(outputI,diag+offset); 
         for(int subdiag = 0; subdiag <  innerWid; subdiag += 4){
            vtf c1 = mult_add(bcast(dF,offset+diag+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+diag+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+diag+subdiag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);
         }
         storea(mp1,output,diag+offset);
         vtf dF3 = preshuffle(dF1,loada(dF,diag+offset+4));
         vtf dX3 = preshuffle(dX1,loada(dX,diag+offset+4));
         vtf mp3 = preshuffle(mp1,loada(output,diag+offset+4));
         vti mpI3= preshuffle(mpI2,loada(outputI,diag+offset+4));

         vtf dF2 = shift2(dF1,dF3);
         vtf dX2 = shift2(dX1,dX3);
         vtf mp2 = shift2(mp1,mp3);
         vti mpI2 =shift2(mpI1,mpI3);
      
         dF1 = shift1(dF1,dF3);
         dX1 = shift1(dX1,dX3);
         mp1 = shift1(mp1,mp3);
         mpI1= shift1(mpI1,mpI3);
         
         for(int subdiag = 1; subdiag < innerWid+1; subdiag += 4){
            vtf c1 = mult_add(bcast(dF,offset+subdiag+diag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(dF1,bcast(dX,offset+subdiag+diag),c1);
            vtf sx = mult(c1,bcast(s,offset+subdiag+diag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);               
 
                c1   = mult_add(bcast(dF,offset+subdiag+diag+1),dX2,c1);
                c1   = mult_add(dF2,bcast(dX,offset+subdiag+diag+1),c1);
                sx   = mult(c1,bcast(s,offset+subdiag+diag+1));
                cmp2 = cmpgtr(mp2,sx);
                mp2  = blend(mp2,sx,cmp2);
         }
         storea(mp1,output,diag+offset+4);
         storea(mp2,output,diag+offset+8);
         vtf dF1 = loadu(dF,diag+offset+3);
         vtf dX1 = loadu(dX,diag+offset+3);
         vtf mp1 = loadu(output,diag+offset+3);
         vti mpI1= loadu(outputI,diag+offset+3);
 
         for(int subdiag = 3; subdiag < innerWid+3; subdiag += 4){
                c1 = mult_add(bcast(dF,offset+subdiag+diag),dX1,c1);
                c1 = mult_add(bcast(dX,offset+subdiag+diag),dF1,c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag));
                cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp4,sx,cmp2);

                t3 += 16;
         }
         storea(mp1,output,diag+offset+12);
      }
   }
   printf("%lu %d \n",t3,n);
}*/
