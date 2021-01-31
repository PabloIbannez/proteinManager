#include <proteinManager/proteinManager.hpp>
#include "bonds.hpp"

int main(int argc,char *argv[]){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");

    auto bondPairVector  = proteinManager::bonds::residuePairBonds(pdb);
    auto angleVector = proteinManager::bonds::residueAngleBonds(pdb);
    auto dihedralVector = proteinManager::bonds::residueDihedralBonds(pdb);

    return EXIT_SUCCESS;
}
