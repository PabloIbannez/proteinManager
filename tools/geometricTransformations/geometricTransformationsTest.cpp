#include "geometricTransformations.hpp"

#include "../centers/centroid.hpp"
#include <random>

#define real double

#define angleMax M_PI/16.0

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");
    
    //proteinManager::geometricTransformations::uniformScaling(pdb,0.1);
    
    //proteinManager::real3 center = proteinManager::computeCentroid(pdb);
    //proteinManager::geometricTransformations::rotation(pdb,center,{1,0,0},M_PI/4.0);
    //std::cerr << center << std::endl;
    //std::cout << pdb    << std::endl;
    
    std::random_device rd;
    std::mt19937 gen(rd());

    std::uniform_real_distribution<real> angleDist(-angleMax,angleMax);
    
    proteinManager::real3 center = proteinManager::computeCentroid(pdb);
    
    for(int i = 0; i<1000; i++){
        //proteinManager::geometricTransformations::randomRotation(pdb.model(0),center,gen);
        proteinManager::geometricTransformations::randomRotation(pdb.model(0),center,angleDist(gen),gen);

        auto atomVector = pdb.model(0).atom();

        std::cout << "#" << std::endl;
        for(auto& atm : atomVector){
            std::cout << atm.getAtomCoord() << " 3.8 0" << std::endl;
        }
    }
    
    
}
