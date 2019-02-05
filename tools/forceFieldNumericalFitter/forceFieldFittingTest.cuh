#include "../forceFieldManager/forceFieldManager.hpp"
#include "../geometricTransformations/geometricTransformations.hpp"

#include "potential.cuh"
#include "forceFieldNumericalFitter.cuh"

#include "../fieldComputing/fieldComputing.cuh"

int main(){
    
    
    proteinManager::STRUCTURE pdbRef;
    proteinManager::STRUCTURE pdbObj;
    
    pdbRef.loadPDB("./examples/1aki.pqr");
    pdbObj.loadPDB("./examples/1akiCG.pqr");
    
    /////////////////////////////////////////////////////////////////
    
    proteinManager::ffManager::forceFieldManager ffM;
    
    ffM.loadForceFieldData("../forceFieldManager/forceFieldModels/gromos.ff");
    ffM.applyForceFieldData(pdbRef);
    
    /////////////////////////////////////////////////////////////////
    
    proteinManager::geometricTransformations::uniformScaling(pdbRef,0.1);
    proteinManager::geometricTransformations::uniformScaling(pdbObj,0.1);
    
    std::cout << "scaled" << std::endl;
    
    /////////////////////////////////////////////////////////////////
    
    using potentialRef = proteinManager::ffManager::attractive6;
    using potentialObj = proteinManager::ffManager::attractive6;
    
    potentialRef potRef(1,0.606962);
    potentialObj potObj(1,0.606962);
    
    /*
    proteinManager::ffManager::forceFieldFitter<potentialRef,potentialObj>::MonteCarloIntegration mc;
    mc.pointsPerIntegral = 100000;
    mc.samplesPerIntegral = 10;
    
    proteinManager::ffManager::forceFieldFitter<potentialRef,potentialObj> test(pdbRef,pdbObj,mc,potRef,potObj);
    */
    
    
    proteinManager::ffManager::forceFieldFitter<potentialRef,potentialObj>::GridIntegration gInt;
    gInt.cellSize = 0.05;
    
    proteinManager::ffManager::forceFieldFitter<potentialRef,potentialObj> test(pdbRef,pdbObj,gInt,potRef,potObj);
    
    
    /*
    proteinManager::ffManager::forceFieldFitter<potential>::Grid_PF_Integration gridPF;
    gridPF.inputFilePath = "./examples/phimap1aki.cube";
    gridPF.lFactor = 0.1;
    gridPF.fFactor = 0.593;
    
    proteinManager::ffManager::forceFieldFitter<potential> ffF(pdbRef,pdbObj,gridPF,pot);
    */
    
    test.computeNewParametersTotalChargeConstraint();
	
    
	/////////////////////////////////////////////////////////////////
    
    potRef = potentialRef(1,0);
    potObj = potentialObj(1,0);
    
    proteinManager::fieldComputing::fieldComputing fC;
    
    std::ofstream ref("refField.dat");
    std::ofstream obj("objField.dat");
    
    fC.init({-0.0356,-0.6309,-2.6852},{5.7644,5.6191,2.7898},0.05);
    
    fC.computeField<potentialRef>(pdbRef,potRef);
    fC.outputIndex_Field(ref);                                        
    
    fC.setFieldValue(0);
    
    fC.computeField<potentialObj>(pdbObj,potObj);
    fC.outputIndex_Field(obj);

    return EXIT_SUCCESS;
}
