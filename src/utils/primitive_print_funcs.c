#ifndef DUMB_PRINT
#define DUMB_PRINT
#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS 1
#endif
#include<stdio.h>
#include<stdlib.h>
/*
#ifdef _WIN32
void writeDoubles(const char* name, const double* t, const int n){
    FILE* f;
    fopen_s(&f, name, "w");
    if(f == NULL){
        fprintf_s(stderr, "error opening %s:  ", name);
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf_s(f, "%.15lf\n", t[i]);
    }
    fclose(f);
}

void appendDoubles(const char* name, const double* t, const int n){
    FILE* f;
    fopen_s(&f, name, "w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf_s(f, "%.15lf\n", t[i]);
    }
    int k = fclose(f);
    if(k != 0){
       perror("fclose");
       exit(1);
    }
}

void writeInts(const char* name, const  int* t, const int n){
    FILE* f;
    fopen_s(&f, name, "w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf_s(f, "%d\n", t[i]);
    }
    fclose(f);
}

void writeLongs(const char* name,const long long* t, const int n){
    FILE* f;
    fopen_s(&f, name, "w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf_s(f, "%lld\n", t[i]);
    }
    fclose(f);
}

#else
*/
void writeDoubles(const char* name, const double* t, const int n){
    FILE* f = fopen(name, "w");
    if(f == NULL){
        fprintf(stderr, "error opening %s:  ", name);
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f, "%.15lf\n", t[i]);
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
        fprintf(f, "%.15lf\n", t[i]);
    }
    int k = fclose(f);
    if(k != 0){
       perror("fclose");
       exit(1);
    }
}

void writeInts(const char* name, const  int* t, const int n){
    FILE* f = fopen(name, "w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f, "%d\n", t[i]);
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
        fprintf(f, "%lli\n", t[i]);
    }
    fclose(f);
}

//#endif
#endif

