#include <stdio.h>


//temporary file to help check inputs and outputs. 


int importnum(const char* filename, double* d, int linecount){
    FILE* fp = fopen(filename,"r");
    if(fp == NULL){
        perror("fopen");
        return -1;
    }
    for(int i = 0; i < n; i++){
        fscanf(fp,"%lf\n",&T[i]);
    }    
    fclose(fp);
    return 0;
}

int importnum(const char* filename, long* I, int linecount){
   FILE* fp = fopen(filename,"r");
   if(fp == NULL){
      perror("fopen");
      return -1;
   }
   for(int i = 0; i < n; i++){
       fscanf(fp,"%d\n",&I[i];
   }
}


int exportnum(const char* filename, double* d, int n){
    FILE* fp = fopen(filename,"w");
    if(f == NULL){
        perror("fopen");  
        return -1;
    }
    for(int i = 0; i < n; i++){
        fprintf(fp,"%.15lf\n",d[i]);
    }
    fclose(fp);
}

int exportnum(const char* filename, int* t, int n){
    FILE* fp = fopen(filename,"w");
    if(fp == NULL){
        perror("fopen");  
        return -1;
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%d\n",t[i]);
    }
    fclose(fp);
}


