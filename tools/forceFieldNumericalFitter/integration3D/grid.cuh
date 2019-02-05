/*Pablo Ibáñez Freire 2018, pablo.ibannez@uam.es*/

/*
 * Data structures where spatial points are stored.
 * 
 * A grid has to be able to handle both CPU and GPU data.
 * A grid class has to ensure that the data from both devices are synchronized
 * and provide an interface to access them.
 */

#ifndef GRID_CUH
#define GRID_CUH

#include <iostream>
#include <fstream>
#include <sstream>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

#include <proteinManager.hpp>

namespace proteinManager{
namespace integrator{
namespace grid{

    class grid3D{
            
            private:
            
                //Data
                bool CPUneedUpdate = false;
                bool GPUneedUpdate = false;
                
                thrust::host_vector<real4> gridDataCPU;
                thrust::device_vector<real4> gridDataGPU;
                
                enum device {CPU,GPU};
                enum operation {read,write,read_write};
                
                //Box
                
                real3 boxMin_;
                real3 boxMax_;
                real3 boxSize_;
                    
                //Grid
                
                real3 cellSize_;
                int3 gridSize_;
                
                bool gridInitialized_ = false;
                
                //Memory management
                void updateState(device dev,operation opt){
                    
                    std::stringstream ss;
                    
                    switch(opt){
                        
                        ///////////////////////////////////////////////////
                        case read:
                            switch(dev){
                            
                                case CPU:
                                    if(CPUneedUpdate) gridDataCPU = gridDataGPU;
                                    CPUneedUpdate = false;
                                    GPUneedUpdate = false;
                                    break;
                                case GPU:
                                    if(GPUneedUpdate) gridDataGPU = gridDataCPU;
                                    CPUneedUpdate = false;
                                    GPUneedUpdate = false;
                                    break;
                                default:
                                    ss.clear();
                                    ss << "A correct device has to be specified. There are two options: CPU or GPU";
                                    throw std::runtime_error(ss.str());
                            } break;
                        ////////////////////////////////////////////////////
                        case write:
                            switch(dev){
                                case CPU:
                                    CPUneedUpdate = false;
                                    GPUneedUpdate = true;
                                    break;
                                case GPU:
                                    CPUneedUpdate = true;
                                    GPUneedUpdate = false;
                                    break;
                                default:
                                    ss.clear();
                                    ss << "A correct device has to be specified. There are two options: CPU or GPU";
                                    throw std::runtime_error(ss.str());
                            } break;
                        ////////////////////////////////////////////////////
                        case read_write:
                            switch(dev){
                                case CPU:
                                    if(CPUneedUpdate) gridDataCPU = gridDataGPU;
                                    CPUneedUpdate = false;
                                    GPUneedUpdate = true;
                                    break;
                                case GPU:
                                    if(GPUneedUpdate) gridDataGPU = gridDataCPU;
                                    CPUneedUpdate = true;
                                    GPUneedUpdate = false;
                                    break;
                                default:
                                    ss.clear();
                                    ss << "A correct device has to be specified. There are two options: CPU or GPU";
                                    throw std::runtime_error(ss.str());
                            } break;
                        ////////////////////////////////////////////////////
                        default:
                            ss.clear();
                            ss << "A correct operation has to be specified. There are three options: read, write or read_write";
                            throw std::runtime_error(ss.str());
                    }
                }
            
            public:
            
                void free(){
                    
                    gridDataCPU.clear();
                    gridDataGPU.clear();
                    
                    gridInitialized_ = false;
                    
                }
                
                //Data access
                
                real4* getCPU_raw(){
                    updateState(CPU,read_write);
                    return gridDataCPU.data();
                }
                
                real4* getGPU_raw(){
                    updateState(GPU,read_write);
                    return thrust::raw_pointer_cast(gridDataGPU.data());
                }
                
                //Get class variables
                
                real3 getBoxMin(){
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return boxMin_;
                    }
                }
                    
                real3 getBoxMax(){
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return boxMax_;
                    }
                }
                
                real3 getCellSize(){
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return cellSize_;
                    }
                }
                    
                int3 getGridSize(){
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return gridSize_;
                    }
                }
                
                size_t getSize(){
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return gridDataCPU.size();
                    }
                }
                
                bool isGridInitialized(){return gridInitialized_;}
                                
                //Value accessing
                
                void setVoxelValue(int x,int y,int z, real4 value){
                    updateState(CPU,write);
                    
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        gridDataCPU[x*gridSize_.y*gridSize_.z + y*gridSize_.z + z] = value;
                    }
                }
                
                void setValue(int x,int y,int z, real value){
                    updateState(CPU,write);
                    
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        gridDataCPU[x*gridSize_.y*gridSize_.z + y*gridSize_.z + z].w = value;
                    }
                }
                
                void setValue(real value){
                    for(int i=0; i < gridSize_.x; i++){
                    for(int j=0; j < gridSize_.y; j++){
                    for(int k=0; k < gridSize_.z; k++){
                            this->setValue(i,j,k,value);
                    }}}
                }
                
                real4 getVoxelValue(int x,int y,int z){
                    updateState(CPU,read);
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return gridDataCPU[x*gridSize_.y*gridSize_.z + y*gridSize_.z + z];
                    }
                }
                
                real getValue(int x,int y,int z){
                    updateState(CPU,read);
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << " can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        return gridDataCPU[x*gridSize_.y*gridSize_.z + y*gridSize_.z + z].w;
                    }
                }
                
                //Box and grid setting up
                void setUpGrid_fixedCellSize(real3 boxMin,real3 boxMax, real3 cellSize){
                    updateState(CPU,write);
                    
                    boxMin_ = boxMin;
                    boxMax_ = boxMax;
                    
                    cellSize_ = cellSize;
                    
                    boxSize_ = boxMax_-boxMin_;
                    
                    boxMax_.x += cellSize_.x*(std::ceil(boxSize_.x/cellSize_.x) - boxSize_.x/cellSize_.x);
                    boxMax_.y += cellSize_.y*(std::ceil(boxSize_.y/cellSize_.y) - boxSize_.y/cellSize_.y);
                    boxMax_.z += cellSize_.z*(std::ceil(boxSize_.z/cellSize_.z) - boxSize_.z/cellSize_.z);
                    
                    boxSize_ = boxMax_-boxMin_;
                    
                    gridSize_.x= boxSize_.x/cellSize_.x;
                    gridSize_.y= boxSize_.y/cellSize_.y;
                    gridSize_.z= boxSize_.z/cellSize_.z;
                    
                    gridDataCPU.resize(gridSize_.x*gridSize_.y*gridSize_.z);
                    
                    gridInitialized_ = true;
                    
                    for(int i=0; i < gridSize_.x; i++){
                    for(int j=0; j < gridSize_.y; j++){
                    for(int k=0; k < gridSize_.z; k++){
                            this->setVoxelValue(i,j,k,{boxMin_.x+i*cellSize_.x+cellSize_.x/real(2.0),
                                                       boxMin_.y+j*cellSize_.y+cellSize_.y/real(2.0),
                                                       boxMin_.z+k*cellSize_.z+cellSize_.z/real(2.0),0});
                    }}}
                    
                    #ifdef DEBUG
                    std::cerr << "Box min: "   << boxMin_ << std::endl;
                    std::cerr << "Box max: "   << boxMax_ << std::endl;
                    std::cerr << "Box size: "  << boxSize_ << std::endl;
                    std::cerr << "Cell size: " << cellSize_ << std::endl;
                    std::cerr << "Grid size: " << gridSize_.x << " " << gridSize_.y << " " << gridSize_.z << std::endl;
                    std::cerr << "Total grid elements: " << gridDataCPU.size() << std::endl;
                    std::cerr << "Grid size (Mb): "      << gridSize_.x*gridSize_.y*gridSize_.z*sizeof(real)*3.0/1024/1024 << std::endl;
                    #endif
                    
                    
                }
                
                //Box and grid setting up
                void setUpGrid_fixedBox(real3 boxMin,real3 boxMax, int3 gridSize){
                    updateState(CPU,write);
                    
                    boxMin_ = boxMin;
                    boxMax_ = boxMax;
                    boxSize_ = boxMax_-boxMin_;
                    
                    gridSize_ = gridSize;
            
                    cellSize_.x= boxSize_.x/gridSize_.x;
                    cellSize_.y= boxSize_.y/gridSize_.y;
                    cellSize_.z= boxSize_.z/gridSize_.z;
                    
                    gridDataCPU.resize(gridSize_.x*gridSize_.y*gridSize_.z);
                    
                    gridInitialized_ = true;
                    
                    for(int i=0; i < gridSize_.x; i++){
                    for(int j=0; j < gridSize_.y; j++){
                    for(int k=0; k < gridSize_.z; k++){
                            this->setVoxelValue(i,j,k,{boxMin_.x+i*cellSize_.x,
                                                       boxMin_.y+j*cellSize_.y,
                                                       boxMin_.z+k*cellSize_.z,0});
                    }}}
                    
                    #ifdef DEBUG
                    std::cerr << "Box min: "   << boxMin_ << std::endl;
                    std::cerr << "Box max: "   << boxMax_ << std::endl;
                    std::cerr << "Box size: "  << boxSize_ << std::endl;
                    std::cerr << "Cell size: " << cellSize_ << std::endl;
                    std::cerr << "Grid size: " << gridSize_.x << " " << gridSize_.y << " " << gridSize_.z << std::endl;
                    std::cerr << "Total grid elements: " << gridDataCPU.size() << std::endl;
                    std::cerr << "Grid size (Mb): "      << gridSize_.x*gridSize_.y*gridSize_.z*sizeof(real)*3.0/1024/1024 << std::endl;
                    #endif
                    
                    
                    
                }
                
                //Output
                
                void outputIndex_Value(std::ostream& out){
                    updateState(CPU,read);
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << "can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                    
                        for(int i=0;i<gridSize_.x;i++){
                        for(int j=0;j<gridSize_.y;j++){
                        for(int k=0;k<gridSize_.z;k++){
                            out << i << " " << j << " " << k << " " << gridDataCPU[i*gridSize_.y*gridSize_.z + j*gridSize_.z + k].w << std::endl;
                        }}}
                    }
                }
                
                void output_CUBE(std::ostream& out,std::string comment1 = " ", std::string comment2 = " ", real lFactor = real(1.0) ,real fFactor = real(1.0)){
                    updateState(CPU,read);
                    std::stringstream ss;
                    if(this->isGridInitialized() == false){
                        ss.clear();
                        ss << "ERROR: " << __FUNCTION__ << "can not be used until grid has been initialized";
                        throw std::runtime_error(ss.str());
                    } else {
                        
                        real coeff=0.5291772108; // Bohr radius Armstrongs
                        
                        out << comment1 << std::endl;
                        out << comment2 << std::endl;
                        
                        out << "1 " << boxMin_*lFactor/coeff << std::endl;
                        
                        out << gridSize_.x << " " << cellSize_.x*lFactor/coeff << " " << 0 << " " << 0 << std::endl;
                        out << gridSize_.y << " " << 0 << " " << cellSize_.y*lFactor/coeff << " " << 0 << std::endl;
                        out << gridSize_.z << " " << 0 << " " << 0 << " " << cellSize_.z*lFactor/coeff << std::endl;
                        
                        out << "1 " << 0 << " " << 0 << " " << 0 << " " << 0 << std::endl;
                        
                        for(int i=0;i<gridSize_.x;i++) {
                            for(int j=0;j<gridSize_.y;j++) {
                                for(int k=0;k<gridSize_.z;k++) {
                                    out << fFactor*this->getValue(i,j,k) << " ";
                                    
                                    if (k % 6 == 5) {out << std::endl;}
                                }
                                
                                out << std::endl;
                            }
                        }
                        
                        
                    }
                }
                
                void input_CUBE(std::string inputFilePath, real lFactor = real(1.0) ,real fFactor = real(1.0)){
                    
                    std::stringstream ss;
                    
                    //Check if file exists
                    std::ifstream inputFile(inputFilePath);
                    if(!inputFile){
                        ss.clear();
                        ss << "File not found: " << inputFilePath;
                        throw std::runtime_error(ss.str());
                    }
                    
                    real3 boxMin;
                    real3 boxMax;
                    int3 gridSize;
                    real3 cellSize;
                    
                    real coeff=0.5291772108; // Bohr radius Armstrongs
                    
                    //Processing file
                    std::string line;
                    
                    int intBuffer;
                    double doubleBuffer1;
                    double doubleBuffer2;
                    double doubleBuffer3;
                    double doubleBuffer4;
                    
                    //First two lines are ignored
                    std::getline(inputFile,line);
                    std::getline(inputFile,line);
                    
                    //////////////////////////////////////////////////////////////////
                    std::getline(inputFile,line);
                    ss.clear();
                    ss.str(line);
                    ss >> intBuffer >> doubleBuffer1 >> doubleBuffer2 >> doubleBuffer3;
                    
                    if(intBuffer != 1){
                        ss.clear();
                        ss << "Format error in " << line << " . First element should be 1";
                        throw std::runtime_error(ss.str());
                    }
                    
                    boxMin = {real(doubleBuffer1*lFactor*coeff),
                              real(doubleBuffer2*lFactor*coeff),
                              real(doubleBuffer3*lFactor*coeff)};
                    
                    //////////////////////////////////////////////////////////////////
                    std::getline(inputFile,line);
                    ss.clear();
                    ss.str(line);
                    ss >> intBuffer >> doubleBuffer1 >> doubleBuffer2 >> doubleBuffer3;
                    
                    if( doubleBuffer2 != 0 or doubleBuffer3 != 0){
                        ss.clear();
                        ss << "Format error in " << line << " . Expected a rectangular box";
                        throw std::runtime_error(ss.str());
                    } else {
                        gridSize.x = intBuffer;
                        cellSize.x = doubleBuffer1*lFactor*coeff;
                    }
                    
                    std::getline(inputFile,line);
                    ss.clear();
                    ss.str(line);
                    ss >> intBuffer >> doubleBuffer1 >> doubleBuffer2 >> doubleBuffer3;
                    
                    if( doubleBuffer1 != 0 or doubleBuffer3 != 0){
                        ss.clear();
                        ss << "Format error in " << line << " . Expected a rectangular box";
                        throw std::runtime_error(ss.str());
                    } else {
                        gridSize.y = intBuffer;
                        cellSize.y = doubleBuffer2*lFactor*coeff;
                    }
                    
                    std::getline(inputFile,line);
                    ss.clear();
                    ss.str(line);
                    ss >> intBuffer >> doubleBuffer1 >> doubleBuffer2 >> doubleBuffer3;
                    
                    if( doubleBuffer1 != 0 or doubleBuffer2 != 0){
                        ss.clear();
                        ss << "Format error in " << line << " . Expected a rectangular box";
                        throw std::runtime_error(ss.str());
                    } else {
                        gridSize.z = intBuffer;
                        cellSize.z = doubleBuffer3*lFactor*coeff;
                    }
                    
                    //////////////////////////////////////////////////////////////////
                    std::getline(inputFile,line);
                    ss.clear();
                    ss.str(line);
                    ss >> intBuffer >> doubleBuffer1 >> doubleBuffer2 >> doubleBuffer3 >> doubleBuffer4;
                    
                    if( intBuffer != 1 or doubleBuffer1 != 0 or doubleBuffer2 != 0 or doubleBuffer3 != 0){
                        ss.clear();
                        ss << "Format error in " << line << " .Expected data 1 0 0 0";
                        throw std::runtime_error(ss.str());
                    }
                    
                    //////////////////////////////////////////////////////////////////
                    
                    boxMax = {boxMin.x+gridSize.x*cellSize.x,
                              boxMin.y+gridSize.y*cellSize.y,
                              boxMin.z+gridSize.z*cellSize.z};
                    
                    this->setUpGrid_fixedBox(boxMin,boxMax,gridSize);
                    
                    //////////////////////////////////////////////////////////////////
                    
                    std::getline(inputFile,line);
                    ss.clear();
                    ss.str(line);
                    
                    for(int i=0;i<gridSize.x;i++) {
                        for(int j=0;j<gridSize.y;j++) {
                            for(int k=0;k<gridSize.z;k++) {
                                ss >> doubleBuffer1;
                                
                                this->setValue(i,j,k,doubleBuffer1*fFactor);
                                
                                if (k % 6 == 5) {
                                    std::getline(inputFile,line);
                                    ss.clear();
                                    ss.str(line);
                                }
                            }
                            std::getline(inputFile,line);
                            ss.clear();
                            ss.str(line);
                        }
                    }
                    
                }
        };
}
}
}

#endif
