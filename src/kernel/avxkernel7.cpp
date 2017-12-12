#include<cstdio>
#include "../mp/avx2arith.hpp"
#define simdWid 4
#define innerWid 128 
using namespace vmth;
typedef  __m256d vtf;
typedef  __m256i vti;


static inline double hmax(__m256d a){
   __m256d b = _mm256_permute2f128_pd(a,a,1);
   __m256d c = max(a,b);
   c = max(c,_mm256_permute_pd(c,5));
   return c[0]; // adjust later
}


// blend casts the mask since most of our comparisons are done in floating point
/*static inline void blend4x1(vti &mpi1, vti &mpi2, vti &mpi3, vti &mpi4,
                            vti &s1,  vti &s2,  vti &s3,  vti &s4,
                            vtf &cmp1,vtf &cmp2,vtf &cmp3,vtf &cmp4){
   mpi1 = blend(mpi1,s1,cmp1);
   mpi2 = blend(mpi2,s2,cmp2);
   mpi3 = blend(mpi3,s3,cmp3);
   mpi4 = blend(mpi4,s4,cmp4);
}*/



static inline void symprodsum4x1(vtf &c1, vtf &c2, vtf &c3, vtf &c4,
                              double* f, double* x, int i, int j){

   c1 = mult_add(bcast(f,i),loada(x,j),mult_add(bcast(x,i),loada(f,j),c1));
   c2 = mult_add(bcast(f,i+4),loada(x,j+4),mult_add(bcast(x,i+4),loada(f,j+4),c2));
   c3 = mult_add(bcast(f,i+8),loada(x,j+8),mult_add(bcast(x,i+8),loada(f,j+8),c3));
   c4 = mult_add(bcast(f,i+12),loada(x,j+12),mult_add(bcast(x,i+12),loada(f,j+12),c4));
}

static inline void symprod4x1(vtf &c1, vtf &c2, vtf &c3, vtf &c4, 
                              vtf &s1, vtf &s2, vtf &s3, vtf &s4,
                              double* s, int i, int j){
   s1 = mult(mult(c1,bcast(s,i)),loada(s,j));
   s2 = mult(mult(c2,bcast(s,i+4)),loada(s,j+4));
   s3 = mult(mult(c3,bcast(s,i+8)),loada(s,j+8));
   s4 = mult(mult(c4,bcast(s,i+12)),loada(s,j+12));
}




// also reasonable quite possibly the fastest  
void  accumTest4_7_11(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   vtf a[innerWid];
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset+diag < n; offset += innerWid){
         for(int suboff = offset; suboff < offset + innerWid; suboff += 4){
            int col  = diag+offset;
            vtf mp1  = loada(output,col);
            vtf mp3  = preshuffle(loada(output,col+4),mp1);
            vtf mp2  = shift1(mp1,mp3);
                mp3  = shift2(mp1,mp3);
            vtf mp4  = loadu(output,col+3);
            vti mpI1 = loada(outputI,col);
            vti mpI3 = preshuffle(mpI1,loada(outputI,col+4));
            vti mpI2 = shift1(mpI1,mpI3);
                mpI3 = shift2(mpI1,mpI3);
            vti mpI4 = loadu(outputI,col+3);
            vtf c2 = loada(Cxy,diag+4);
            vtf c3 = loada(Cxy,diag+8);
            vtf c4 = loada(Cxy,diag+12);
            vtf c1;
            
            for(int subdiag = diag; subdiag < diag+innerWid; subdiag += 4){
                c1 = mult_add(bcast(dF,suboff),loada(dX,subdiag+suboff),c1);
                c2 = mult_add(bcast(dF,suboff+4),loada(dX,subdiag+suboff+4),c2);
                c3 = mult_add(bcast(dF,suboff+8),loada(dX,subdiag+suboff+8),c3);
                c4 = mult_add(bcast(dF,suboff+12),loada(dX,subdiag+suboff+12),c4);
           
                c1 = mult_add(bcast(dX,suboff),loada(dF,subdiag+suboff),loada(Cxy,subdiag));
                c2 = mult_add(bcast(dX,suboff+4),loada(dF,subdiag+suboff+4),c2);
                c3 = mult_add(bcast(dX,suboff+8),loada(dF,subdiag+suboff+8),c3);
                c4 = mult_add(bcast(dX,suboff+12),loada(dF,subdiag+suboff+12),c4);
             
                vtf s1 = mult(bcast(s,suboff),loada(s,subdiag+suboff));
                vtf s2 = mult(bcast(s,suboff+4),loada(s,subdiag+suboff));
                vtf s3 = mult(bcast(s,suboff+8),loada(s,subdiag+suboff));
                vtf s4 = mult(bcast(s,suboff+12),loada(s,subdiag+suboff)); 

                s1 = mult(s1,c1); 
                s2 = mult(s2,c2);
                s3 = mult(s3,c3);
                s4 = mult(s4,c4);               
                a[subdiag-diag] = max(a[subdiag-diag],s1);
                a[subdiag-diag+1]=max(a[subdiag-diag+1],s2);
                a[subdiag-diag+2]=max(a[subdiag-diag+2],s3);
                a[subdiag-diag+3]=max(a[subdiag-diag+3],s4);
               
/*                storea(max(loada(a,subdiag-diag),s1),a,subdiag-diag);
                storea(max(loada(a,subdiag-diag+4),s2),a,subdiag-diag);
                storea(max(loada(a,subdiag-diag+8),s3),a,subdiag-diag);

  */    
           /*     storea(max(loada(output,offset),s1),output,offset);
                storea(max(loada(output,offset+4),s2),output,offset);
                storea(max(loada(output,offset+8),s3),output,offset);
                storea(max(loada(output,offset+12),s4),output,offset); 
             *//* output[offset]    = hmax(s1);
                output[offset+4]  = hmax(s2);
                output[offset+8]  = hmax(s3);
                output[offset+12] = hmax(s4);
                */
                mp1  = max(mp1,s1);
                mp2  = max(mp2,s2);
                mp3  = max(mp3,s3);
                mp4  = max(mp4,s4);
     
                s1   = cmpgtr(s1,mp1);
                s2   = cmpgtr(s2,mp2);
                s3   = cmpgtr(s3,mp3);
                s4   = cmpgtr(s4,mp4);
        
                mpI1 = blend(mpI1, bcast(subdiag),mp1);
                mpI2 = blend(mpI2, bcast(subdiag)+4,mp2);
                mpI3 = blend(mpI3, bcast(subdiag+8),mp3);
                mpI4 = blend(mpI4, bcast(subdiag+12),mp4);
                
                t3 += 16;                
                storea(c4,Cxy,subdiag+12);
                // not supported on all compilers, can make adjustments later
                c4 = c3;
                c3 = c2;
                c2 = c1;
            }
            
            storea(c1,Cxy,diag+innerWid-12);
            storea(c2,Cxy,diag+innerWid-8);
            storea(c3,Cxy,diag+innerWid-4);
            storea(mp1,output,diag+offset);
            storea(mp2,output,diag+offset+4);
            storea(mp3,output,diag+offset+8); 
            storea(mp4,output,diag+offset+12);
            storea(mpI1,outputI,diag+offset);
            storea(mpI2,outputI,diag+offset+4);
            storea(mpI3,outputI,diag+offset+8);
            storea(mpI4,outputI,diag+offset+12);
            for(int i = 0; i < innerWid; i+= 4){
               output[offset+i] = max(hmax(a[i]),output[offset+i]);
            }
         }
      }
   }
   printf("%lu %d \n",t3,n);
}

void  accumTest4_7_11_alt(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n){
   unsigned long t3 = 0;
   #pragma omp parallel for schedule(static)
   for(int diag = 0; diag < n; diag += innerWid){
      for(int offset = 0; offset+diag < n; offset += innerWid){
         for(int suboff = offset; suboff < offset + innerWid; suboff += 4){
            int col  = diag+offset;
            vtf mp1  = loada(output,col);
            vtf mp3  = preshuffle(loada(output,col+4),mp1);
            vtf mp2  = shift1(mp1,mp3);
                mp3  = shift2(mp1,mp3);
            vtf mp4  = loadu(output,col+3);
            vti mpi1 = loada(outputI,col);
            vti mpi3 = preshuffle(mpi1,loada(outputI,col+4));
            vti mpi2 = shift1(mpi1,mpi3);
                mpi3 = shift2(mpi1,mpi3);
            vti mpi4 = loadu(outputI,col+3);


            for(int subdiag = diag; subdiag < diag+innerWid; subdiag += 4){
            vtf c1 = mult_add(bcast(dX,suboff),loada(dF,subdiag+suboff),loada(Cxy,subdiag));
                c1 = mult_add(bcast(dF,suboff),loada(dX,subdiag+suboff),c1);
            vtf s1 = mult(c1,bcast(s,suboff));
                mp1  = max(mp1,s1);
                s1   = cmpgtr(s1,mp1);
                mpi1 = blend(mpi1, bcast(subdiag),s1);
                
                c1 = mult_add(bcast(dX,suboff+4),loada(dF,subdiag+suboff+4),c1);
                c1 = mult_add(bcast(dF,suboff+4),loada(dX,subdiag+suboff+4),c1);
                s1 = mult(c1,bcast(s,suboff+4));
                output[subdiag+offset+4] = hmax(s1);
                mp2 = max(mp2,s1);
                s1 = cmpgtr(s1,mp2);
                mpi2 = blend(mpi2,bcast(subdiag+4),s1);
             
                c1 = mult_add(bcast(dX,suboff+8),loada(dF,subdiag+suboff+8),c1);
                c1 = mult_add(bcast(dF,suboff+8),loada(dX,subdiag+suboff+8),c1);
                s1 = mult(c1,bcast(s,suboff+8));
                output[subdiag+offset+8] = hmax(s1);
                mp3 = max(mp3,s1);
                s1 = cmpgtr(s1,mp3);
                mpi3 = blend(mpi3,bcast(subdiag+4),s1);

                c1 = mult_add(bcast(dX,suboff+12),loada(dF,subdiag+suboff+12),c1);
                c1 = mult_add(bcast(dF,suboff+12),loada(dX,subdiag+suboff+12),c1);
                s1 = mult(c1,bcast(s,suboff+12));
                output[subdiag+offset+12] = hmax(s1); 
                mp4 = max(mp4,s1);
                s1 = cmpgtr(s1,mp4);       
                mpi4 = blend(mpi4,bcast(subdiag+12),s1);

          //      t3 += 16;                
                storea(c1,Cxy,subdiag);
            }
            storea(mp1,output,diag+offset);
            storea(mp2,output,diag+offset+4);
            storea(mp3,output,diag+offset+8); 
            storea(mp4,output,diag+offset+12);
            storea(mpi1,outputI,diag+offset);
            storea(mpi2,outputI,diag+offset+4);
            storea(mpi3,outputI,diag+offset+8);
            storea(mpi4,outputI,diag+offset+12);
         //   printf(" diag: %d suboff: %d diag+offset:%d \n",diag,suboff,diag+offset);
         }
      }
   }
   printf("%lu %d \n",t3,n);
}



