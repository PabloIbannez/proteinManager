#include "analyticChargeFitting.hpp"

int main(){
    
    proteinManager::STRUCTURE pdbRef;
    proteinManager::STRUCTURE pdbObj;
    
    pdbRef.loadPDB("./examples/1aki.pqr");
    pdbObj.loadPDB("./examples/1akiCG_MARRINK.pqr");
    
    proteinManager::chargeFitting::analyticChargeFitting aChgF(pdbRef,pdbObj);
    
    aChgF.fitChargesAndPositions(1000);
    
    std::ofstream out("1akiChgPosFit.pqr");
    out << pdbObj << std::endl;
}
