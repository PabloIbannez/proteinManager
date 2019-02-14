#ifndef INTEGRATOR_3D_MC_CUH
#define INTEGRATOR_3D_MC_CUH

#include "saruprng.cuh"

namespace proteinManager{
namespace integrator{
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
					cudaFree(d_temp_storage);
				}
				
				cudaFree(sumResult);
				cudaStreamDestroy(integratorStream);
			}
            
            void free(){
                
                if(d_temp_storage != NULL){
					cudaFree(d_temp_storage);
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
    
    
    
}}}

#endif
