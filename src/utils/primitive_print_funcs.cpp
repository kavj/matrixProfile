#ifndef DUMB_PRINT
#define DUMB_PRINT

#include<cstdio>
#include<cstdlib>
#include<unistd.h>

void writeDoubles(const char* name, const double* t, const int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        fprintf(stderr,"error opening %s:  ",name);
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%.15lf\n",t[i]);
    }
    fclose(f);
}

void appendDoubles(const char* name, const double* t, const int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%.15lf\n",t[i]);
    }
    int k = fclose(f);
    if(k != 0){
       perror("fclose");
       exit(1);
    }
}

void writeInts(const char* name,const  int* t, const int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%d\n",t[i]);
    }
    fclose(f);
}

void writeLongs(const char* name,const long long* t, const int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%d\n",t[i]);
    }
    fclose(f);
}

#endif


