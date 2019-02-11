#include <proteinManager/proteinManager.hpp>

#include "../../tools/massesManager/massesManager.hpp"
#include "../../tools/forceFieldManager/forceFieldManager.hpp"

#include "../../tools/coarseGrained/coarseGrainedManager.hpp"
#include "../../tools/coarseGrained/coarseGrainedMappingSchemes.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"

#include "../../tools/forceFieldNumericalFitter/forceFieldNumericalFitter.cuh"
#include "../../tools/forceFieldNumericalFitter/potential.cuh"

#include "../../tools/enm/enm.hpp"

using namespace proteinManager;

void loadChargesSurf(MODEL& mdlIn,std::string ChgSurfFilePath){
    
    struct chainRes{
        std::string chainId;
        int resSeq;
    };
    
    std::map<int,chainRes> oneChainMap;
    chainRes chRbuffer;
    
    int i = 1;
    for(CHAIN&   ch  : mdlIn.chain())  {
    for(RESIDUE& res : ch.residue())   {
        chRbuffer.chainId = ch.getChainId();
        chRbuffer.resSeq  = res.getResSeq();
        oneChainMap[i]=chRbuffer;
        i++;
    }}
    
    ////////////////////////////////////////////////////////////////////
    std::ifstream inputFile(ChgSurfFilePath);
    
    std::string line;
    std::stringstream ss;
    
    int    resSeqBuffer;
    double chgBuffer;
    
    for(CHAIN&   ch  : mdlIn.chain())  {
    for(RESIDUE& res : ch.residue())   {
        
        if(res.atom().size() > 1){
            ss.clear();
            ss << "ERROR. Expected one atom per risude only";
            throw std::runtime_error(ss.str());
        }
        res.atom()[0].setAtomCharge(0);
        res.atom()[0].setAtomSurf(false);
    }}
    
    while(std::getline(inputFile,line)){
        ss.clear();
        ss.str(line);
        
        ss >> resSeqBuffer >> chgBuffer;
        mdlIn.chain(oneChainMap[resSeqBuffer].chainId).residue(oneChainMap[resSeqBuffer].resSeq).atom()[0].setAtomCharge(chgBuffer);
        mdlIn.chain(oneChainMap[resSeqBuffer].chainId).residue(oneChainMap[resSeqBuffer].resSeq).atom()[0].setAtomSurf(true);
    }
    
}

int main(int argc, char *argv[]){
    
    STRUCTURE proteinIn;
    STRUCTURE proteinOut;
    proteinIn.loadPDB("./proteins/barnase_barstar/barnase_barstar.pdrs");
    proteinIn.renumber();
    
    /////////////////////Load masses and ff/////////////////////////////
    
    //masses
    massesManager massesM;
    
    massesM.loadMassesData("../../tools/massesManager/massesData/atomMasses.dat");
    massesM.applyMassesData(proteinIn);
    
    //ff
    
    ffManager::forceFieldManager ffM;
    
    ffM.loadForceFieldData("../../tools/forceFieldManager/forceFieldModels/gromos.ff");
    ffM.applyForceFieldData(proteinIn);
    
    //Apply CG
    
    coarseGrainedManager::coarseGrainedGenerator cg;
    
    cg.loadCGmodel("../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/aminoAcid2bead_RES2BEAD.map", \
                   "../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/bead2atom_RES2BEAD.map");
                   
    cg.applyCoarseGrainedMap<coarseGrainedManager::coarseGrainedMappingSchemes::sasa>(proteinIn,proteinOut);
    
    //std::cout << proteinIn << std::endl;
    //std::cout << proteinOut << std::endl;
    
    ////////////////////////////////////////////////////////////////////

    loadChargesSurf(proteinOut.model(1),"./proteins/barnase/barnase.charge");
    loadChargesSurf(proteinOut.model(2),"./proteins/barstar/barstar.charge");
    
    //////////////////////////Fit C6 C12////////////////////////////////
    
    geometricTransformations::uniformScaling(proteinIn,0.1);
    geometricTransformations::uniformScaling(proteinOut,0.1);
    
    //C6
    
    using potentialRef_C6 = ff_fitting::attractive6;
    using potentialObj_C6 = ff_fitting::attractive6;
    
    potentialRef_C6 potRef_C6(1,0.606962);
    potentialObj_C6 potObj_C6(1,0.606962);
    
    ff_fitting::forceFieldFitter<potentialRef_C6,potentialObj_C6>::MonteCarloIntegration mc_C6;
    mc_C6.pointsPerIntegral = 100000;
    mc_C6.samplesPerIntegral = 10;
    
    ff_fitting::forceFieldFitter<potentialRef_C6,potentialObj_C6> model1_C6(proteinIn.model(1),proteinOut.model(1),mc_C6,potRef_C6,potObj_C6);
    model1_C6.computeNewParametersTotalChargeConstraint();
    
    ff_fitting::forceFieldFitter<potentialRef_C6,potentialObj_C6> model2_C6(proteinIn.model(2),proteinOut.model(2),mc_C6,potRef_C6,potObj_C6);
    model2_C6.computeNewParametersTotalChargeConstraint();
    
    //C12
    
    using potentialRef_C12 = ff_fitting::repulsive12;
    using potentialObj_C12 = ff_fitting::repulsive12;
    
    potentialRef_C12 potRef_C12(1,0.779078);
    potentialObj_C12 potObj_C12(1,0.779078);
    
    ff_fitting::forceFieldFitter<potentialRef_C12,potentialObj_C12>::MonteCarloIntegration mc_C12;
    mc_C12.pointsPerIntegral = 100000;
    mc_C12.samplesPerIntegral = 10;
    
    ff_fitting::forceFieldFitter<potentialRef_C12,potentialObj_C12> model1_C12(proteinIn.model(1),proteinOut.model(1),mc_C12,potRef_C12,potObj_C12);
    model1_C12.computeNewParametersTotalChargeConstraint();
    
    ff_fitting::forceFieldFitter<potentialRef_C12,potentialObj_C12> model2_C12(proteinIn.model(2),proteinOut.model(2),mc_C12,potRef_C12,potObj_C12);
    model2_C12.computeNewParametersTotalChargeConstraint();
    
    ////////////////////////////////////////////////////////////////////
    
    std::ofstream out("barnase_barstar.top");
    
    int atomCount = 1;
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()  ){
    for(RESIDUE& res : ch.residue() ){
    for(ATOM&    atm : res.atom()   ){
        out << atomCount << " " << res.getResName() << " " << mdl.getModelId() << " " << atm.getAtomSASA() << " " << ((atm.getAtomSurf()==true)?1:0) << " " << mdl.getModelId() << " " 
            << atm.getAtomMass() << " " << atm.getAtomCoord() << " " << atm.getAtomC12() << " " << atm.getAtomC6() << " "
            << atm.getAtomCharge() << " " << atm.getAtomSolvE() << std::endl;
        atomCount ++;
    }}}}
    
    
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
