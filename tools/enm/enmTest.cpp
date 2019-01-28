#include "enm.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.setSerialOverWrite(true);
    pdb.setResSeqOverWrite(true);
    
    pdb.loadPDB("./examples/streptavidin_4jo6_CG.pdb");
    
    
    proteinManager::enm<proteinManager::enm_models::REACH> enmTest;
    
    enmTest.computeENM(pdb);
    
    pdb.setOutputFormat(proteinManager::DATA_FORMAT::XYZ);
    
    std::cout << pdb << std::endl;
    std::cout << enmTest << std::endl;
    
    
    return EXIT_SUCCESS;
}
