

CC = ~/gcc-versions/8.3/bin/g++ 

CPPFLAGS = -I$(SRC)

CFLAGS =  -O3 -march=native -funroll-loops -fprefetch-loop-arrays

CXXFLAGS = $(CFLAGS) -std=c++11 

DIRS = -Isrc -Isrc/solvers/ -Isrc/utils/ 

LFLAGS =  -lm  

ASMFLAGS = -S   -masm=intel -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti -fverbose-asm

SRC = src/utils/moments.cpp src/main.cpp  src/solvers/pearson.cpp src/utils/alloc.cpp src/utils/primitive_print_funcs.cpp  

opt:
	$(CC) $(DIRS) $(SRC) -o mpSolve $(CFLAGS) $(LFLAGS) 

asm:
	$(CC) $(SRC)   $(CFLAGS) $(ASMFLAGS) 

debug:
	$(CC) $(SRC)  -g -fopenmp -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti

vdebug:
	$(CC) $(SRC)  -g -O2 -mavx2 -mfma -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti
 

ompbuild: 
	$(CC) $(SRC)  $(CFLAGS) -O3 -fopenmp $(LFLAGS)

