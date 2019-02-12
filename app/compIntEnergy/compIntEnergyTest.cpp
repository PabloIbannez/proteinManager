#include <proteinManager/proteinManager.hpp>

#include "compIntEnergy.hpp"

using namespace proteinManager;

int main(){
    
    STRUCTURE p1;
    STRUCTURE p2;
    
    p1.loadPDB("");
    p2.loadPDB("");
    
    compInt::sasaPot pot(INFINITY,78,0.7);
    
    real energySep = 0;
    energySep  += compInt::compIntEnergy(p1.model(0),pot);
    energySep  += compInt::compIntEnergy(p2.model(0),pot);
    
    real energyComp = 0;
    energyComp += compInt::compIntEnergy(p1.model(0),p2.model(0),pot);
    
    
    
    return EXIT_SUCCESS;
}
