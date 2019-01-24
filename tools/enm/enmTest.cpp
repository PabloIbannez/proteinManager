#include "enm.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.setSerialOverWrite(true);
    pdb.setResSeqOverWrite(true);
    
    pdb.loadPDB("1aki.pdb");
    
    
    proteinManager::enm<proteinManager::enm_models::caOrellana> enmTest;
    
    enmTest.computeENM(pdb);
    
    //std::cout << pdb << std::endl;
    
    return EXIT_SUCCESS;
}
