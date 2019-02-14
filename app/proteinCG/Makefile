ARCH=61

DEBUG_CPU= -O3
DEBUG_GPU= -O3

CXX=g++
CPU= $(DEBUG_CPU) -funroll-loops -ffinite-math-only -fno-signaling-nans -fno-math-errno -fno-signed-zeros -frename-registers -march=native -fPIC   -lproteinManager

BASIC_LINE= nvcc -lineinfo $(DEBUG_GPU) -ccbin="$(CXX)" -Xcompiler="$(CPU)" -src-in-ptx -gencode arch=compute_$(ARCH),code=sm_$(ARCH) -x cu -std=c++11 --expt-relaxed-constexpr

proteinCG:
	$(BASIC_LINE) -lproteinManager proteinCG_test.cpp
