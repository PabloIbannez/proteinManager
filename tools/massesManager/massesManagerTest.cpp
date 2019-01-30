#include <proteinManager.hpp>

#define DEBUG

#include "massesManager.hpp"
#include "../centroid/centerOfMass.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1akiClear.pdb");
    
    proteinManager::massesManager mM;
    
    mM.loadMassesData("massesData/atomMasses.dat");
    mM.applyMassesData(pdb);
    
    std::cout << proteinManager::computeCenterOfMass(pdb) << std::endl;
}
