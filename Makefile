

CC = g++

CFLAGS =  -O3 -march=native -funroll-loops -std=c++11  

LFLAGS =  -lm  

ASMFLAGS = -S   -masm=intel -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti -fverbose-asm

SRC = src/utils/alloc.cpp src/utils/cov.cpp src/simple_driver.cpp  src/solvers/pearson.cpp src/utils/primitive_print_funcs.cpp  


opt:
	$(CC) -o scamp $(SRC)   $(CFLAGS) -fprefetch-loop-arrays $(LFLAGS) 

asm:
	$(CC) $(SRC)   $(CFLAGS) $(ASMFLAGS) 

debug:
	$(CC) $(SRC)  -g -fopenmp -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti

vdebug:
	$(CC) $(SRC)  -g -O2 -mavx2 -mfma -fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti
 

ompbuild: 
	$(CC) -o scamp $(SRC)  $(CFLAGS) -O3 -fopenmp $(LFLAGS)

