

#CC = gcc
CC = ~/gcc-versions/7.3.0/bin/g++


CFLAGS =  -O3 -march=native -std=c++11 -mavx2 -mfma    

LFLAGS =  -lm  

ASMFLAGS = -S  -masm=intel -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti -fverbose-asm

SRC = src/simple_driver.cpp src/solvers/ref_tiled_pearson.cpp  src/solvers/pearson.cpp  src/utils/primitive_print_funcs.cpp  

opt:
	$(CC) $(SRC)   $(CFLAGS) $(LFLAGS) 

asm:
	$(CC) $(SRC)   $(CFLAGS) $(ASMFLAGS) 

debug:
	$(CC) $(SRC)  -g -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti

ompbuild: 
	$(CC) $(SRC)  $(CFLAGS) -O3 -fopenmp $(LFLAGS)

