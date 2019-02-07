ARCH=61

DEBUG_CPU= -O3
DEBUG_GPU= -O3

CXX=g++
CPU= $(DEBUG_CPU) -funroll-loops -ffinite-math-only -fno-signaling-nans -fno-math-errno -fno-signed-zeros -frename-registers -march=native -fPIC -L ../../shared -I ../../include -lproteinManager

BASIC_LINE= nvcc -lineinfo $(DEBUG_GPU) -ccbin="$(CXX)" -Xcompiler="$(CPU)" -src-in-ptx -gencode arch=compute_$(ARCH),code=sm_$(ARCH) -x cu -std=c++11 --expt-relaxed-constexpr

all:
	$(BASIC_LINE)  -L ../../shared -I ../../include -lproteinManager forceFieldFittingTest.cuh

PF:
	$(BASIC_LINE)  -L ../../shared -I ../../include -lproteinManager forceFieldFitting_PF_Test.cuh
