#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 32 
#define tileSz 128 
#define alignedInner (innerWid/4-1)

using namespace vmth;
#define ksz 64
typedef  __m256d vtf;
typedef  __m256i vti;



void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   double a[tileSz*simdWid];
   double b[tileSz*simdWid];
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += tileSz){
      for(int offset = 0; offset < diag; offset += tileSz){
         for(int suboff = offset; suboff < offset + tileSz; suboff++){
            vtf s1 = bcast(s,suboff);
            vtf dFi = bcast(dF,suboff);
            vtf dXi = bcast(dX,suboff);
            vtf mpi = bcast(output,suboff);
            vtf mpiI = bcast(output,suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = loada(Cxy,subdiag);
               c1 = mult_add(loadu(dF,suboff+subdiag),dXi,mult_add(dFi,loadu(dX,suboff+subdiag),c1));
               vtf c2 = mult(c1,s1);
               c2 = mult(c2,loadu(s,suboff+subdiag));
               mpi = max(c2,mpi);
               storea(c1,Cxy,subdiag-diag);
               t3 += 4;
            }
            storeu(mpi,output,suboff);
         }
      }   
   }
   printf("%lu %d \n",t3,n);
}




void  accumTest4_2(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   double a[tileSz*simdWid];
   double b[tileSz*simdWid];
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
            vtf c1 = loada(Cxy,subdiag);
            for(int suboff = offset; suboff < offset + innerWid; suboff++){
               c1 = mult_add(bcast(dF,suboff),loadu(dX,subdiag),c1);
               c1 = mult_add(bcast(dX,suboff),loadu(dF,subdiag),c1);
               //storea(c1,a,simdWid*(suboff-offset));
               vtf s1 = mult(bcast(s,subdiag),loadu(s,suboff));
               vtf s2 = mult(c1,s1);
               storea(s2,b,simdWid*(suboff-offset));
               t3 += 4;
            }
            storea(c1,Cxy,subdiag);
            c1 = loada(output,subdiag);
            
            for(int i = 0; i < innerWid; i+=simdWid){
               c1 = max(c1,loada(b,i));
            }
            storea(c1,output,subdiag);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}


void  accumTest4_3(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
            vtf c1 = loada(Cxy,subdiag);
            vtf c2 = loada(output,subdiag);
            for(int suboff = 0; suboff < innerWid; suboff++){
               c1 = mult_add(bcast(dF,suboff+offset),loadu(dX,subdiag),c1);
               c1 = mult_add(bcast(dX,suboff+offset),loadu(dF,subdiag),c1);
               vtf s1 = mult(bcast(s,subdiag),loadu(s,suboff+offset));
               vtf s2 = mult(c1,s1);
               c2 = max(s2,c2);
               storea(c2,output,4*suboff+subdiag);
               t3 += 4;
            }
            storea(c1,Cxy,subdiag);

         }
      }
   }
   printf("%lu %d \n",t3,n);
}


void  accumTest4_7(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff++){
            vtf cmx = bcast(output,suboff);
            vtf s1 = bcast(s,suboff);
            vtf dx1 = bcast(dX,suboff);
            vtf df1 = bcast(dF,suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = mult_add(df1,loadu(dX,suboff),loada(Cxy,subdiag));
              // vtf c2 = mult_add(df
               c1 = mult_add(dx1,loadu(dF,suboff),c1);
               vtf s2 = mult(s1,loadu(s,suboff+offset));
               s2 = mult(c1,s2);
            //   c2 = max(s2,c2);
               t3 += 4;
               storea(s2,a,innerWid*suboff+(subdiag-diag));
               storea(c1,Cxy,subdiag);
            }
           // storea(c2,output,offset+suboff);
         }
         for(int i = 0; i < innerWid; i+= simdWid){
            vtf c1 = loada(output,diag);
            c1 = max(c1,loada(a,i));
            c1 = max(c1,loadu(a,i+innerWid+1));
            c1 = max(c1,loadu(a,i+2*(innerWid+1)));
            c1 = max(c1,loadu(a,i+3*(innerWid+1)));
            storea(c1,output,diag);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}


void  accumTest4_4(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff++){
            vtf c2 = bcast(output,offset+suboff);
            vtf s1 = bcast(s,offset+suboff);
            vtf dx1 = bcast(dX,offset+suboff);
            vtf df1 = bcast(dF,offset+suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = mult_add(df1,loadu(dX,offset+suboff),loada(Cxy,subdiag));
               c1 = mult_add(dx1,loadu(dF,offset+suboff),c1);
               vtf s2 = mult(s1,loadu(s,offset+suboff));
               s2 = mult(c1,s2);
               c2 = max(s2,c2);
               t3 += 4;
               storea(s2,a,innerWid*suboff+(subdiag-diag));
               storea(c1,Cxy,subdiag);
            }
            storea(c2,output,offset+4*suboff);
         }
         for(int i = 0; i < innerWid; i+= simdWid){
            vtf c1 = loada(output,diag);
            c1 = max(c1,loada(a,i));
            c1 = max(c1,loadu(a,i+innerWid+1));
            c1 = max(c1,loadu(a,i+2*(innerWid+1)));
            c1 = max(c1,loadu(a,i+3*(innerWid+1)));
            storea(c1,output,diag);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}


/* shuffles are really really annoying
 * I'm testing an approach where I pre-align the rows of the implicit pairwise distance matrix,
 * Since rows are broadcast, they aren't a big deal. Now the unaligned portion consumes only one set of registers.
 * */

/*
void  accumTest4_8(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += 4){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff+=4){
            vtf dx1 = bcast(dX,suboff+offset);
            vtf df1 = bcast(dF,suboff+offset);
            vtf s   = bcast(s,suboff+offset);
            for(int subd = 0; subd < innerWid; subd += 4){
               vtf dff1  = loada(dF,suboff+offset+diag);
               vtf dxx1  = loada(dX,suboff+offset+diag);
               vtf dxx2_ = preshuffle(dxx1,loada(dX,suboff+offset+diag+4));
               vtf dff1  = loada(dF,suboff+offset+diag);
               vtf dff2_ = preshuffle(dff1,loada(dF,suboff+offset+diag+4));
               vtf ss1   = loada(s,suboff+offset+diag);
               vtf ss2_   = preshuffle(s,loada(s,suboff+offset+diag+4));               

	       vtf c1 = loada(Cxy,diag+subd);
               
}

void  accumTest4_6(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff+=4){
            vtf c2  = bcast(output,suboff+offset);
            vtf s1  = bcast(s,suboff+offset);
            vtf s2  = bcast(s,suboff+offset+1);
            vtf s3  = bcast(s,suboff+offset+2);
            vtf s4  = bcast(s,suboff+offset+3);
            vtf dx1 = bcast(dX,suboff+offset);
            vtf dx2 = bcast(dX,suboff+offset+1);
            vtf dx3 = bcast(dX,suboff+offset+2);
            vtf dx4 = bcast(dX,suboff+offset+3);
            vtf df1 = bcast(dF,suboff+offset);
            vtf df2 = bcast(dF,suboff+offset+1);
            vtf df3 = bcast(dF,suboff+offset+2);
            vtf df4 = bcast(dF,suboff+offset+3);
            vtf cmx1 = bcast(output,suboff+offset);
            vtf cmx2 = bcast(output,suboff+offset+1);
            vtf cmx3 = bcast(output,suboff+offset+2);
            vtf cmx4 = bcast(output,suboff+offset+3);
          for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf dxx1 = loada(dX,subdiag);
               vtf dxx2 = preshuffle(dxx1,loada(dX,subdiag+simdWid));
               vtf dff1 = loada(dF,subdiag);
               vtf dff2 = preshuffle(dff1,loada(dF,subdiag+simdWid));
               vtf ss1  = loada(s,subdiag);
               vtf ss2  = preshuffle(ss1,loada(s,subdiag+simdWid));

               vtf c1 = mult_add(dx1,dff1,loada(Cxy,subdiag));
               c1     = mult_add(df1,dxx1,c1);
               vtf s5 = mult(c1,mult(s1,ss1));
               storea(s5,a,subdiag-diag);
               cmx1 = max(cmx1,s5);

               c1 = mult_add(dx2,shift1(dff1,dff2),c1); 
               c1 = mult_add(df2,shift1(dxx1,dxx2),c1);
               s5 = mult(c1,mult(s2,shift1(ss1,ss2)));
               storea(s5,a,subdiag-diag+4);
               cmx2 = max(cmx2,s5);
    
               c1 = mult_add(dx3,shift2(dff1,dff2),c1);
               c1 = mult_add(df3,shift2(dxx1,dxx2),c1);
               s5 = mult(c1,mult(s3,shift2(ss1,ss2)));
               storea(s5,a,subdiag-diag+8);
               cmx3 = max(cmx3,s5);
      
               storea(c1,Cxy,subdiag);
               t3 += 16;
            }
            storea(c2,output,suboff+offset);
            storea(cmx1,output,suboff+offset);
            storea(cmx2,output,suboff+offset+4);
            storea(cmx3,output,suboff+offset+8);
            storea(cmx4,output,suboff+offset+12);
         }
         for(int i = 0; i < innerWid; i+= simdWid){
            vtf c1 = loada(output,diag);
            c1 = max(c1,loada(a,i));
            c1 = max(c1,loadu(a,i+(innerWid-1)));
            c1 = max(c1,loadu(a,i+2*(innerWid-1)));
            c1 = max(c1,loadu(a,i+3*(innerWid-1)));
            storea(c1,output,diag+i);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}



*/


void  accumTest4_6(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   double a[innerWid*innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += 4){

            // this might increase the number of extract operations, which is unfortunate. It's arguably best to do this with an unpack or shuffle
            // to get upper and lower separately, but it's actually a lot of instructions either way. This does cut down on unaligned loads considerably

            // it might be worth batching the maximizing phase due to the potential unaligned writebacks


            vtf s1  = loada(s,diag+offset);
            vtf dX1 = loada(dX,diag+offset);
            vtf dF1 = loada(dF,diag+offset);
            vtf mp1 = loada(output,diag+offset);
            vti mpI1= loada(outputI,diag+offset); 
 
            // we preload the second ones because this removes 2 unaligned acccesses
            vtf dF2 = preshuffle(dF1,loada(dF,diag+offset+4));
            vtf dX2 = preshuffle(dX1,loada(dX,diag+offset+4));
            vtf s2  = preshuffle(s1,loada(s,diag+offset+4));
            vtf mp2 = preshuffle(mp1,loada(output,diag+offset+4));
            //vti mpI2= preshuffle(mpI2,loada(output,diag+offset+4));
         
         for(int subdiag = 0; subdiag <  innerWid; subdiag += 4){
            vtf c1 = mult_add(bcast(dF,offset),dX1,loada(Cxy,subdiag));
                c1 = mult_add(bcast(dX,offset),dF1,c1);
            vtf sx = mult(s1,bcast(s,offset));
                sx = mult(c1,sx);
                storea(sx,a,subdiag);
                c1 = mult_add(bcast(dF,offset+1),shift1(dX1,dX2),c1);
                c1 = mult_add(shift1(dF1,dF2),bcast(dX,offset+1),c1);
                sx = mult(shift1(s1,s2),bcast(s,offset+1));
                sx = mult(c1,sx);
                storea(sx,a,subdiag+4);
                c1 = mult_add(bcast(dF,offset+2),shift2(dX1,dX2),c1);
                c1 = mult_add(shift2(dF1,dF2),bcast(dX,offset+2),c1);
                sx = mult(shift1(s1,s2),bcast(s,offset+2));
                sx = mult(c1,sx);
                storea(sx,a,subdiag+8);
                c1 = mult_add(bcast(dF,offset+3),loadu(dX,diag+offset+3),c1);
                c1 = mult_add(bcast(dX,offset+3),loadu(dX,diag+offset+3),c1);
                sx = mult(loadu(s,diag+offset+3),bcast(s,offset+3));
                sx = mult(c1,sx);
                storea(sx,a,subdiag+12);
         
         }
         
         for(int subdiag = 0; subdiag < innerWid; subdiag += 16){
                mp1 = max(mp1,loada(a,subdiag));
                mp1 = max(mp1,loadu(a,subdiag+3));
                mp1 = max(mp1,loadu(a,subdiag+2));
                mp1 = max(mp1,loadu(a,subdiag+1));
         }
         t3 += innerWid*4;
      }
   }
   printf("%lu %d \n",t3,n);
}



void  accumTest4_5(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;

   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset < diag; offset += innerWid){
         for(int suboff = 0; suboff < innerWid; suboff++){
            vtf c2 = bcast(output,suboff);
            vtf s1 = bcast(s,suboff);
            vtf dx1 = bcast(dX,suboff);
            vtf df1 = bcast(dF,suboff);
            for(int subdiag = diag; subdiag < diag + innerWid; subdiag += simdWid){
               vtf c1 = mult_add(df1,loadu(dX,suboff),loada(Cxy,subdiag));
               c1 = mult_add(dx1,loadu(dF,suboff),c1);
               vtf s2 = mult(s1,loadu(s,suboff+offset));
               s2 = mult(c1,s2);
               c2 = max(s2,c2);
               t3 += 4;
               storea(c1,Cxy,subdiag);
            }
            storea(c2,output,4*suboff);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}



