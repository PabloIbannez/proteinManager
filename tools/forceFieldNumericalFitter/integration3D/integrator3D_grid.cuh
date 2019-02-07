#ifndef INTEGRATOR_3D_GRID_CUH
#define INTEGRATOR_3D_GRID_CUH

#include "grid.cuh"

namespace proteinManager{
namespace integrator{
namespace integrator3D_grid_ns{
        
    template<class function>
    __global__ void functionGrid3DEvaluation2vector(real4* grid_ptr,
                                                    real*  values_ptr,
                                                    int N,
                                                    function f){
                                                        
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if(i >= N) return;
        
        real3 point = {grid_ptr[i].x,grid_ptr[i].y,grid_ptr[i].z};
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
        
        real3 point = {grid_ptr[i].x,grid_ptr[i].y,grid_ptr[i].z};
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
            
            void init_precompFunct(std::string inputFilePath){
                g3D.inputIndex_Value(inputFilePath);
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
    
}}}

#endif
