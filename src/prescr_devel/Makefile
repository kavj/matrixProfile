
CC = ~/gcc-versions/7.3.0/bin/g++
#CC = g++

CFLAGS =  -O3 -march=haswell  -std=c++11 -mavx2 -mfma -funroll-loops -fprefetch-loop-arrays  #-DMKL_ILP64 -I ~/anaconda3/include   

#LFLAGS = -lm
LFLAGS =  -lm #-L ~/anaconda3/lib/ -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl 

ASMFLAGS = -S  -masm=intel #-fno-asynchronous-unwind-tables -fno-exceptions -fno-rtti -fverbose-asm



SRC = tiled_Pearson.cpp ../utils/cov.cpp exact_Pearson.h prescr_driver.cpp  




#	$(CC) $(SRC) $(CFLAGS) $(LFLAGS) 

opt:
	$(CC) $(SRC)   $(CFLAGS) $(LFLAGS) 

asm:
	$(CC) $(SRC)   $(CFLAGS) $(ASMFLAGS) 

ompbuild: 
	$(CC) $(SRC)  $(CFLAGS) -O3 -fopenmp $(LFLAGS)

mskylake:
	$(CC) $(SRC) $(ASMFLAGS) -O3 -march=skylake-avx512 -mfma -funroll-loops -fprefetch-loop-arrays



all:    opt
