#define DEBUG

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
    
    //cgTop.fit_C6(1,1,0.9);
    //cgTop.fit_C12(1,1,0.779078);
    //
    //cgTop.fit_C6(2,1,0.9);
    //cgTop.fit_C12(2,1,0.779078);
    
    std::shared_ptr<STRUCTURE> strOut = cgTop.getStructureOut();
    
    for(MODEL& mdl   : strOut->model()) {
    for(CHAIN& ch    : mdl.chain())     {
    for(RESIDUE& res : ch.residue())    {
    for(ATOM&    atm : res.atom())      {
        atm.scaleAtomSASA(0.01);
    }}}}
    
    std::ofstream out("barnase_barstar.top");
    
    cgTop.output(out);
    
    ////////////////////////////////////////////////////////////////////
    
    //compInt::sasaPot pot(INFINITY,78,0.7);
    //std::cout << "E: " << compInt::compIntEnergy(strOut->model(1),strOut->model(2),pot) << std::endl;
    
    ////////////////////////////////////////////////////////////////////
    
    std::ofstream out_enm("barnase_barstar.bond");
    
    //proteinManager::enm<proteinManager::enm_models::REACH_nm> ENM;    
    //proteinManager::enm<proteinManager::enm_models::caOrellana_nm> ENM;    
    proteinManager::enm<proteinManager::enm_models::basic_nm> ENM;    
    
    ENM.computeENM(strOut->model(1));
    ENM.computeENM(strOut->model(2));
    
    out_enm << ENM << std::endl;
    
    return EXIT_SUCCESS;
}
