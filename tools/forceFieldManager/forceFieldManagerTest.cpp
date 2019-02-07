#define DEBUG

#include "forceFieldManager.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");
    
    proteinManager::ffManager::forceFieldManager fM;
    
    fM.loadForceFieldData("./forceFieldModels/gromos.ff");
    
    fM.applyForceFieldData(pdb);
    
    return EXIT_SUCCESS;
}
