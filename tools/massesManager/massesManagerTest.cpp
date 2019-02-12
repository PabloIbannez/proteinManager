#include <proteinManager/proteinManager.hpp>

#define DEBUG

#include "massesManager.hpp"
#include "../centers/centerOfMass.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");
    
    proteinManager::massesManager mM;
    
    mM.loadMassesData("massesData/atomMasses.dat");
    mM.applyMassesData(pdb.model()[0]);
    
    std::cout << proteinManager::computeCenterOfMass(pdb.model()[0]) << std::endl;
}
