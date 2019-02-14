#include "./proteinCG.hpp"

#include "../compIntEnergy/compIntEnergy.hpp"
#include "../compIntEnergy/potential.hpp"

#include "../../tools/enm/enm.hpp"

using namespace proteinManager;

int main(int argc, char *argv[]){
        
    STRUCTURE proteinIn;
    
    proteinIn.loadPDB("./proteins/barnase_barstar/barnase_barstar.pdrs");
    proteinIn.renumber();
    
    geometricTransformations::uniformScaling(proteinIn,0.1);
    
    coarseGrainedTop cgTop;
    
    cgTop.loadStructure(proteinIn);
    
    cgTop.loadChargesSurf(1,"./proteins/barnase/barnase.charge");
    cgTop.loadChargesSurf(2,"./proteins/barstar/barstar.charge");
    
    cgTop.fit_C6(1,1,0.606962);
    cgTop.fit_C12(1,1,0.779078);
    
    cgTop.fit_C6(2,1,0.606962);
    cgTop.fit_C12(2,1,0.779078);
    
    std::ofstream out("barnase_barstar.top");
    
    cgTop.output(out);
    
    ////////////////////////////////////////////////////////////////////
    
    //proteinManager::enm<proteinManager::enm_models::REACH> ENM;    
    //
    //geometricTransformations::uniformScaling(proteinIn,10);
    //geometricTransformations::uniformScaling(proteinOut,10);
    //
    //ENM.computeENM(proteinOut.model(1));
    //ENM.computeENM(proteinOut.model(2));
    
    return EXIT_SUCCESS;
}
