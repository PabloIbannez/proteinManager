/*Pablo Ibáñez Freire 2018, pablo.ibannez@uam.es*/

/*
 * How each integrator is initialized depends on the algorithm used.
 * In the initialization is where the region of integration is indicated.
 * Once initialized to calculate the integral on the specified region, the following method is used:
 * 
 * real computeIntegral<function>(function f)
 * 
 *      -"function" is a class that must have the operator () overloaded as follows:
 *          __host__ __device__ real operator()(real3 point), which returns the value of the function at the specified point.
 * 
 *      -"f" is an instance of that class.
 */

#ifndef INTEGRATOR_3D_CUH
#define INTEGRATOR_3D_CUH

#include <iostream>
#include <string>
#include <sstream>

#include <cuda.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <cub/cub.cuh>

//#define DEBUG

#include <proteinManager.hpp>

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

#include "integrator3D_grid.cuh"
#include "integrator3D_MC.cuh"

#endif
