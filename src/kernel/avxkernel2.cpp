#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 128 

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

         for(int subdiag = diag+offset; subdiag <  diag+offset+innerWid; subdiag += 4){
     
            vtf c1 = mult_add(bcast(dF,offset+diag+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+diag+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+diag+subdiag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);
                mpI1= blend(mpI1,ind,cmp2);
                 
                c1 = mult_add(bcast(dF,offset+subdiag+diag+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+subdiag+diag+1),c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+1));
                cmp2 = cmpgtr(mp2,sx);
                mp2 = blend(mp2,sx,cmp2); 
                mpI2= blend(mpI2,ind+1,cmp2); 

                c1   = mult_add(bcast(dF,offset+subdiag+diag+2),dX3,c1);
                c1   = mult_add(dF3,bcast(dX,offset+subdiag+diag+2),c1);
                sx   = mult(c1,bcast(s,offset+subdiag+diag+2));
                cmp2 = cmpgtr(mp3,sx);
                mp3  = blend(mp3,sx,cmp2);
                mpI3 = blend(mpI3,ind+2,cmp2);

                c1 = mult_add(bcast(dF,offset+subdiag+diag+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+subdiag+diag+3),dF4,c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+3));
                cmp2 = cmpgtr(mp4,sx);
                mp4 = blend(mp4,sx,cmp2);
                mpI4 = blend(mpI4,ind+3,cmp2);

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






void  accumTest4_7_2(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += 4){
         vtf dX1 = bcast(offset);
         vtf dX2 = bcast(offset+1);
         vtf dX3 = bcast(offset+2);
         vtf dX4 = bcast(offset+3);
         vtf mp1 = bcast(output,offset);
         vtf mp2 = bcast(output,offset+1);
         vtf mp3 = bcast(output,offset+2);
         vtf mp4 = bcast(output,offset+3);
         vtf mpI1 =bcast(outputI,offset);
         vtf mpI2 =bcast(outputI,offset+1);
         vtf mpI3 =bcast(outputI,offset+2);
         vtf mpI4 =bcast(outputI,offset+3); 
/*            vtf dX1 = loada(dX,diag+offset);
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
*/
         for(int subdiag = diag+offset; subdiag <  diag+offset+innerWid; subdiag += 4){
            vti ind = _mm256_set_epi64x(offset+subdiag,offset+subdiag+1,offset+subdiag+2,offset+subdiag+3);
            vti inc = bcast(1);           
            vtf dX12 = loada(dX,subdiag);
            vtf dX22 = loadu(dX,subdiag+1);
            vtf dX23 = loadu(dX,subdiag+2);
            vtf dX24 = loadu(dX,subidag+3);
            vtf dF12 = loada(dF,subdiag);
            vtf dF22 = loadu(dF,subdiag+1);
            vtf dF32 = loadu(dF,subdiag+2);
            vtf dF42 = loadu(dF,subdiag+3);
     
            vtf c1 = mult_add(bcast(dF,offset+diag+subdiag),dX1,loada(Cxy,subdiag+diag));
                c1 = mult_add(bcast(dX,offset+diag+subdiag),dF1,c1);
            vtf sx = mult(c1,bcast(s,offset+diag+subdiag));
            vtf cmp2 = cmpgtr(mp1,sx);
                mp1 = blend(mp1,sx,cmp2);
                mpI1= blend(mpI1,ind,cmp2);
                 
                c1 = mult_add(bcast(dF,offset+subdiag+diag+1),dX2,c1);
                c1 = mult_add(dF2,bcast(dX,offset+subdiag+diag+1),c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+1));
                cmp2 = cmpgtr(mp2,sx);
                mp2 = blend(mp2,sx,cmp2); 
                mpI2= blend(mpI2,ind+1,cmp2); 

                c1   = mult_add(bcast(dF,offset+subdiag+diag+2),dX3,c1);
                c1   = mult_add(dF3,bcast(dX,offset+subdiag+diag+2),c1);
                sx   = mult(c1,bcast(s,offset+subdiag+diag+2));
                cmp2 = cmpgtr(mp3,sx);
                mp3  = blend(mp3,sx,cmp2);
                mpI3 = blend(mpI3,ind+2,cmp2);

                c1 = mult_add(bcast(dF,offset+subdiag+diag+3),dX4,c1);
                c1 = mult_add(bcast(dX,offset+subdiag+diag+3),dF4,c1);
                sx = mult(c1,bcast(s,offset+subdiag+diag+3));
                cmp2 = cmpgtr(mp4,sx);
                mp4 = blend(mp4,sx,cmp2);
                mpI4 = blend(mpI4,ind+3,cmp2);

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
