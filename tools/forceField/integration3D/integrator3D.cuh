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

#include "vectorData.cuh"
#include "grid.cuh"
#include "saruprng.cuh"

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
   if (code != cudaSuccess) 
   {
      fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
      if (abort) exit(code);
   }
}

using namespace vectorData;

namespace integrator{

namespace integrator3D_grid_ns{
        
    template<class function>
    __global__ void functionGrid3DEvaluation2vector(real4* grid_ptr,
                                                    real*  values_ptr,
                                                    int N,
                                                    function f){
                                                        
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if(i >= N) return;
        
        real3 point = make_real3(grid_ptr[i]);
        values_ptr[i] = f(point);
        
        //printf("%f %f %f, %f\n",point.x,point.y,point.z,values_ptr[i]);
        
        return;
    }
    
    template<class function>
    __global__ void functionGrid3D_PF_Evaluation2vector(real4* grid_ptr,
                                                    real*  values_ptr,
                                                    int N,
                                                    function f){
                                                        
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if(i >= N) return;
        
        real3 point = make_real3(grid_ptr[i]);
        values_ptr[i] = f(point,grid_ptr[i].w);
        
        //printf("%f %f %f %f, %f\n",point.x,point.y,point.z,grid_ptr[i].w,values_ptr[i]);
        
        return;
    }
    
    class integrator3D_grid{
        
        private:
            
            grid::grid3D g3D;
            thrust::device_vector<real> functionValuesOverGrid;
            
            bool initialized = false;
            bool precompFunct = false;
            
            //CUDA variables
            
            cudaStream_t integratorStream;
            
            //Reduction variables
            
            real* sumResult;
            
            int      N_temp = 0;
            void*    d_temp_storage = NULL;
            size_t   temp_storage_bytes = 0;
            
        public:
            
            integrator3D_grid(){
                cudaStreamCreate(&integratorStream);
				cudaMallocManaged(&sumResult,sizeof(real));
			}
			
			~integrator3D_grid(){
				
				if(d_temp_storage != NULL){
					cudaFree(&d_temp_storage);
				}
				
				cudaFree(&sumResult);
				cudaStreamDestroy(integratorStream);
			}
            
            void free(){
                
                if(d_temp_storage != NULL){
					cudaFree(&d_temp_storage);
                    d_temp_storage = NULL;
                    N_temp = 0;
				}
                
                functionValuesOverGrid.clear();
                g3D.free();
                
                initialized  = false;
                precompFunct = false;
            }
            
            void init(real3 boxMin,real3 boxMax,int3 gridSize){
                this->init_fixedBox(boxMin,boxMax,gridSize);
            }
            
            void init(real3 boxMin,real3 boxMax,int gridSize){
                this->init_fixedBox(boxMin,boxMax,{gridSize,gridSize,gridSize});
            }
            
            void init_fixedBox(real3 boxMin,real3 boxMax,int3 gridSize){
                g3D.setUpGrid_fixedBox(boxMin,boxMax,gridSize);
                functionValuesOverGrid.resize(g3D.getSize());
                initialized = true;
            }
            
            void init_fixedCellSize(real3 boxMin,real3 boxMax, real3 cellSize){
                g3D.setUpGrid_fixedCellSize(boxMin,boxMax,cellSize);
                functionValuesOverGrid.resize(g3D.getSize());
                initialized = true;
            }
            
            void init_fixedCellSize(real3 boxMin,real3 boxMax, real cellSize){
                this->init_fixedCellSize(boxMin,boxMax,{cellSize,cellSize,cellSize});
            }
            
            void init_fixedBox(real3 boxMin,real3 boxMax,int gridSize){
                this->init_fixedBox(boxMin,boxMax,{gridSize,gridSize,gridSize});
            }
            
            void init_precompFunct(std::string inputFilePath, real coordFactor = real(1.0) ,real precompFunctFactor = real(1.0)){
                g3D.input_CUBE(inputFilePath, coordFactor,precompFunctFactor);
                functionValuesOverGrid.resize(g3D.getSize());
                initialized = true;
                precompFunct = true;
            }
            
            template<class function>
            real computeIntegral(function& f){
                
                std::stringstream ss;
                
                if(initialized == false){
                    ss.clear();
                    ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                    throw std::runtime_error(ss.str());
                }
                    
                real4* grid_ptr = g3D.getGPU_raw();
				real*  functionValuesOverGrid_ptr = thrust::raw_pointer_cast(functionValuesOverGrid.data());
                
                int N = g3D.getSize();
				
				///////////////////////////////////////////////////////
				
				int Nthreads=128;
				int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
                
                functionGrid3DEvaluation2vector<function><<<Nblocks, Nthreads,0,integratorStream>>>(grid_ptr,
                                                                                                    functionValuesOverGrid_ptr,
                                                                                                    N,
                                                                                                    f);
				cudaStreamSynchronize(integratorStream);
				
				///////////////////////////////////////////////////////
				
				//reduction
				if(N > N_temp){
					N_temp = N;
					if(d_temp_storage != NULL){
						cudaFree(&d_temp_storage);
					}
					d_temp_storage = NULL;
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverGrid_ptr, sumResult, N,integratorStream);
					cudaMalloc(&d_temp_storage, temp_storage_bytes);
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverGrid_ptr, sumResult, N,integratorStream);
				} else {
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverGrid_ptr, sumResult, N,integratorStream);
				}
				cudaStreamSynchronize(integratorStream);
                
                real3 cellSize = g3D.getCellSize();
				return sumResult[0]*cellSize.x*cellSize.y*cellSize.z;
            }
            
            template<class function>
            real computeIntegral_PF(function& f){
                
                std::stringstream ss;
                
                if(initialized == false){
                    ss.clear();
                    ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                    throw std::runtime_error(ss.str());
                }
                
                if(precompFunct == false){
                    ss.clear();
                    ss << "ERROR: " << __FUNCTION__ << " can not be used if a precomputed function has not been given";
                    throw std::runtime_error(ss.str());
                    
                }
                    
                real4* grid_ptr = g3D.getGPU_raw();
				real*  functionValuesOverGrid_ptr = thrust::raw_pointer_cast(functionValuesOverGrid.data());
                
                int N = g3D.getSize();
				
				///////////////////////////////////////////////////////
				
				int Nthreads=128;
				int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
                
                functionGrid3D_PF_Evaluation2vector<function><<<Nblocks, Nthreads,0,integratorStream>>>(grid_ptr,
                                                                                                        functionValuesOverGrid_ptr,
                                                                                                        N,
                                                                                                        f);
                    
				//gpuErrchk( cudaPeekAtLastError() );
                //gpuErrchk( cudaDeviceSynchronize() );
                
                cudaStreamSynchronize(integratorStream);
				
				///////////////////////////////////////////////////////
				
				//reduction
				if(N > N_temp){
					N_temp = N;
					if(d_temp_storage != NULL){
						cudaFree(&d_temp_storage);
					}
					d_temp_storage = NULL;
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverGrid_ptr, sumResult, N,integratorStream);
					cudaMalloc(&d_temp_storage, temp_storage_bytes);
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverGrid_ptr, sumResult, N,integratorStream);
				} else {
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverGrid_ptr, sumResult, N,integratorStream);
				}
				cudaStreamSynchronize(integratorStream);
                
                real3 cellSize = g3D.getCellSize();
				return sumResult[0]*cellSize.x*cellSize.y*cellSize.z;
            }
    };
    
}

namespace integrator3D_MC_ns{
    
    template<class function>
    __global__ void functionMC3DEvaluation2vector(real*  values_ptr,
                                                  real3 boxMin,
                                                  real3 boxMax,
                                                  int N,
                                                  int iterationNumber,
                                                  unsigned int prngSeed,
                                                  function f){
                                                        
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if(i >= N) return;
        
        Saru saruPRNG(i,iterationNumber,prngSeed);
        
        real3 point = {saruPRNG.f(boxMin.x,boxMax.x),saruPRNG.f(boxMin.y,boxMax.y),saruPRNG.f(boxMin.z,boxMax.z)};
        
        values_ptr[i] = f(point);
        
        //printf("%f %f %f, %f\n",point.x,point.y,point.z,values_ptr[i]);
        
        return;
    }
    
    class integrator3D_MC{
        
        private:
            
            thrust::device_vector<real> functionValuesOverUniformDistribution;
            
            //Box
            
            real3 boxMin_;
            real3 boxMax_;
            
            ////////////////////////////////////////////
            
            bool initialized = false;
            
            //CUDA variables
            
            cudaStream_t integratorStream;
            
            //Reduction variables
            
            real* sumResult;
            
            int      N_temp = 0;
            void*    d_temp_storage = NULL;
            size_t   temp_storage_bytes = 0;
            
            //PRNG
            
            int N_ = 0; //Number of random points where the function is evalued
            int iterationNumber_ = 0;
            unsigned int prngSeed_ = 0xFEEDC0DE;
            
        public:
            
            integrator3D_MC(){
                cudaStreamCreate(&integratorStream);
				cudaMallocManaged(&sumResult,sizeof(real));
			}
			
			~integrator3D_MC(){
				
				if(d_temp_storage != NULL){
					cudaFree(&d_temp_storage);
				}
				
				cudaFree(&sumResult);
				cudaStreamDestroy(integratorStream);
			}
            
            void free(){
                
                if(d_temp_storage != NULL){
					cudaFree(&d_temp_storage);
                    d_temp_storage = NULL;
                    N_temp = 0;
				}
                
                functionValuesOverUniformDistribution.clear();
                initialized = false;
            }
            
            void init(real3 boxMin, real3 boxMax, int N){
                boxMin_ = boxMin;
                boxMax_ = boxMax;
                
                N_=N;
                functionValuesOverUniformDistribution.resize(N_);
                
                initialized = true;
            }
            
            void setUpPRNG(unsigned int prngSeed){
                prngSeed_ = prngSeed;
            }
            
            template <class function>
            real computeIntegral(function f){
                
                std::stringstream ss;
                
                if(initialized == false){
                    ss.clear();
                    ss << "ERROR: " << __FUNCTION__ << " can not be used until box has been initialized";
                    throw std::runtime_error(ss.str());
                }
                
				real*  functionValuesOverUniformDistribution_ptr = thrust::raw_pointer_cast(functionValuesOverUniformDistribution.data());
                
                iterationNumber_ ++;
				
				///////////////////////////////////////////////////////
				
				int Nthreads=128;
				int Nblocks=N_/Nthreads + ((N_%Nthreads)?1:0);
				
				functionMC3DEvaluation2vector<function><<<Nblocks, Nthreads,0,integratorStream>>>(functionValuesOverUniformDistribution_ptr,
                                                                                                  boxMin_,
                                                                                                  boxMax_,
                                                                                                  N_,
                                                                                                  iterationNumber_,
                                                                                                  prngSeed_,
                                                                                                  f);
				cudaStreamSynchronize(integratorStream);
				
				///////////////////////////////////////////////////////
				
				//reduction
				if(N_ > N_temp){
					N_temp = N_;
					if(d_temp_storage != NULL){
						cudaFree(&d_temp_storage);
					}
					d_temp_storage = NULL;
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverUniformDistribution_ptr, sumResult, N_, integratorStream);
					cudaMalloc(&d_temp_storage, temp_storage_bytes);
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverUniformDistribution_ptr, sumResult, N_, integratorStream);
				} else {
					cub::DeviceReduce::Sum(d_temp_storage, temp_storage_bytes, functionValuesOverUniformDistribution_ptr, sumResult, N_, integratorStream);
				}
				cudaStreamSynchronize(integratorStream);
                
                real V = (boxMax_-boxMin_).x*(boxMax_-boxMin_).y*(boxMax_-boxMin_).z;
                
				return sumResult[0]*(V/N_);
                
            }
            
            template <class function>
            real2 computeIntegralAverage(function f,int itrNum){
                
                std::vector<real> ensemble;
                
                for(int i =0; i < itrNum; i++){
                    ensemble.push_back(this->computeIntegral<function>(f));
                }
                
                real mean = 0;
                for(int i =0; i < itrNum; i++){mean += ensemble[i];}
                mean /= itrNum;
            
                real s = 0;
                for(int i =0; i < itrNum; i++){s += std::pow(ensemble[i]-mean,2);}
                s = std::sqrt(s/(itrNum-1));
                
                return {mean,s/std::sqrt(real(itrNum))};
                
            }
        
        
    };
    
    
    
}
}

#endif
