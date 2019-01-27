#include "enm.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.setSerialOverWrite(true);
    pdb.setResSeqOverWrite(true);
    
    pdb.loadPDB("streptavidin_1mk5_CG.pdb");
    
    
    proteinManager::enm<proteinManager::enm_models::caOrellana> enmTest;
    
    enmTest.computeENM(pdb);
    
    //std::cout << pdb << std::endl;
    
    std::cout << enmTest << std::endl;
    
    //std::cout << pdb << std::endl;
    
    return EXIT_SUCCESS;
}
