#ifndef FIELD_COMPUTING_CUH
#define FIELD_COMPUTING_CUH

#include <iostream> 
#include <string>
#include <fstream>
#include <sstream>

#include <proteinManager/proteinManager.hpp>
#include <integration3D/src/pIntegrator.cuh>

namespace proteinManager {
namespace fieldComputing{
    
    template <class function>
    __global__ void functionGrid3DEvaluation(vectorData::real4* grid_ptr,
                                             int N,
                                             function f){
                                            
        int i = blockIdx.x*blockDim.x + threadIdx.x;
        if(i >= N) return;
        
        grid_ptr[i].w += f({grid_ptr[i].x,grid_ptr[i].y,grid_ptr[i].z});
        
        //printf("%f\n",grid_ptr[i].w);
                                            
    }
                                        
    
    class fieldComputing{
        
        private:
        
            integrator::grid::grid3D grid;
            
            bool initialized_ = false;
            
            cudaStream_t fieldComputingStream;
            
        public:
        
            fieldComputing(){
                cudaStreamCreate(&fieldComputingStream);
            }
        
            fieldComputing(real3 boxMin, real3 boxMax, real cellSize){
                this->init(boxMin,boxMax,cellSize);
                cudaStreamCreate(&fieldComputingStream);
            }
            
            ~fieldComputing(){
                cudaStreamDestroy(fieldComputingStream);
            }
            
            void init(real3 boxMin, real3 boxMax, real cellSize){
                grid.setUpGrid_fixedCellSize({boxMin.x,boxMin.y,boxMin.z},
                                             {boxMax.x,boxMax.y,boxMax.z},
                                             {cellSize,cellSize,cellSize});
                initialized_ = true;
            }
            
            template<class potential>
            void computeField(proteinManager::STRUCTURE& structIn, potential& pot){
                
                std::stringstream ss;
                if(initialized_ == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                }
                
                int N  = grid.getSize();
                vectorData::real4* grid_ptr = grid.getGPU_raw();
                
                /////////////////////////////////////////////////////////////////////////////////
                
                int Nthreads=128;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
                
                for(proteinManager::MODEL& md : structIn.model())   {
                for(proteinManager::CHAIN& ch : md.chain())         {
                for(proteinManager::RESIDUE& res : ch.residue())    {
                for(proteinManager::ATOM& atm : res.atom())         {
                    
                    auto potAP = pot.getPotentialAtPoint({atm.getAtomCoord().x/10,atm.getAtomCoord().y/10,atm.getAtomCoord().z/10},
                                                          potential::getInteractionParm(atm));
                                                                                  
                    functionGrid3DEvaluation<<<Nblocks, Nthreads,0,fieldComputingStream>>>(grid_ptr,
                                                                                           N,
                                                                                           potAP);
                    cudaStreamSynchronize(fieldComputingStream);
                    
                }}}}
                
            }
            
            template<class potential>
            void computeRefField(proteinManager::STRUCTURE& structIn, potential& pot){
                
                std::stringstream ss;
                if(initialized_ == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                }
                
                int N  = grid.getSize();
                vectorData::real4* grid_ptr = grid.getGPU_raw();
                
                /////////////////////////////////////////////////////////////////////////////////
                
                int Nthreads=128;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
                
                for(proteinManager::MODEL& md : structIn.model())   {
                for(proteinManager::CHAIN& ch : md.chain())         {
                for(proteinManager::RESIDUE& res : ch.residue())    {
                for(proteinManager::ATOM& atm : res.atom())         {
                    
                    auto potAP = pot.getPotentialAtPointRef({atm.getAtomCoord().x/10,atm.getAtomCoord().y/10,atm.getAtomCoord().z/10},
                                                             potential::getInteractionParm(atm));
                                                                                  
                    functionGrid3DEvaluation<<<Nblocks, Nthreads,0,fieldComputingStream>>>(grid_ptr,
                                                                                           N,
                                                                                           potAP);
                    cudaStreamSynchronize(fieldComputingStream);
                    
                }}}}
                
            }
            
            template<class potential>
            void computeObjField(proteinManager::STRUCTURE& structIn, potential& pot){
                
                std::stringstream ss;
                if(initialized_ == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                }
                
                int N  = grid.getSize();
                vectorData::real4* grid_ptr = grid.getGPU_raw();
                
                /////////////////////////////////////////////////////////////////////////////////
                
                int Nthreads=128;
                int Nblocks=N/Nthreads + ((N%Nthreads)?1:0);
                
                for(proteinManager::MODEL& md : structIn.model())   {
                for(proteinManager::CHAIN& ch : md.chain())         {
                for(proteinManager::RESIDUE& res : ch.residue())    {
                for(proteinManager::ATOM& atm : res.atom())         {
                    
                    auto potAP = pot.getPotentialAtPointObj({atm.getAtomCoord().x/10,atm.getAtomCoord().y/10,atm.getAtomCoord().z/10},
                                                             potential::getInteractionParm(atm));
                                                                                  
                    functionGrid3DEvaluation<<<Nblocks, Nthreads,0,fieldComputingStream>>>(grid_ptr,
                                                                                           N,
                                                                                           potAP);
                    cudaStreamSynchronize(fieldComputingStream);
                    
                }}}}
                
            }
            
            void setFieldValue(real value){
                grid.setValue(value);
            }
            
            void outputIndex_Field(std::ostream& out){
                grid.outputIndex_Value(out);
            }
            
            void output_CUBE(std::ostream& out,std::string comment1 = " ", std::string comment2 = " ", real lFactor = real(1.0), real fFactor = real(1.0)){
                grid.output_CUBE(out,comment1,comment2,lFactor,fFactor);
            }
            
        
    };
}
}

#endif
