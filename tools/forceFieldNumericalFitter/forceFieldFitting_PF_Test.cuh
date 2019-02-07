#define DEBUG

#include "../geometricTransformations/geometricTransformations.hpp"

#include "potential.cuh"
#include "forceFieldNumericalFitter.cuh"

#include "../fieldComputing/fieldComputing.cuh"

int main(){
    
    proteinManager::STRUCTURE pdbObj;
    
    pdbObj.loadPDB("./examples/1akiCG.pqr");
    
    /////////////////////////////////////////////////////////////////
    
    using potentialObj = proteinManager::ff_fitting::debye;
    
    potentialObj potObj(1,5.22651,78.000,10,0.05);
    
    proteinManager::ff_fitting::forceFieldFitter_PrecomputedField<potentialObj>::GridIntegrationPF gridPF;
    gridPF.precompFieldFilePath = "1akiField.xyz";
    
    proteinManager::ff_fitting::forceFieldFitter_PrecomputedField<potentialObj> test(pdbObj,gridPF,potObj);

    test.computeNewParameters();
    
    ////////////////////////////////////////////////////////////////////
    
    potObj = potentialObj(1,5.22651,78.000,INFINITY,0);
    
    proteinManager::fieldComputing::fieldComputing fC;
    
    std::ofstream obj("objField_PF.dat");
    
    fC.init({-1.77600,-7.78300,-27.88000},{58.98400,57.39500,28.70200},{0.4746875,0.5092031,0.442046});
    
    fC.computeField<potentialObj>(pdbObj,potObj);
    fC.outputIndex_Field(obj);
    
    ////////////////////////////////////////////////////////////////////
    
    std::ofstream ref("refField_PF.dat");
    
    proteinManager::integrator::grid::grid3D gr;
    
    gr.inputIndex_Value("1akiField.xyz");
    gr.outputIndex_Value(ref);

    return EXIT_SUCCESS;
}
