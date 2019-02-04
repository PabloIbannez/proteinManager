#include "geometricTransformations.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");
    
    proteinManager::geometricTransformations::uniformScaling(pdb,0.1);
    
    std::cout << pdb << std::endl;
}
