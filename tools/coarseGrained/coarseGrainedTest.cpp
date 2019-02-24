#define DEBUG
#include "coarseGrainedManager.hpp"
#include "../massesManager/massesManager.hpp"
#include "../forceFieldManager/forceFieldManager.hpp"
#include "../geometricTransformations/geometricTransformations.hpp"

int main(){
    
    proteinManager::coarseGrainedManager::coarseGrainedGenerator cg;
    
    cg.loadCGmodel("./coarseGrainedModels/RES2BEAD/aminoAcid2bead_RES2BEAD.map","./coarseGrainedModels/RES2BEAD/bead2atom_RES2BEAD.map");
    //cg.loadCGmodel("./coarseGrainedModels/MARRINK/aminoAcid2bead_MARRINK.map","./coarseGrainedModels/MARRINK/bead2atom_MARRINK.map");
    
    ////////////////////////////////////////////
    
    std::string pdbInputPath = "./examples/barnase_barstar.pdrs";
    std::string pdbOutputPath = "./examples/barnase_barstar_CG.pdb";
    
    proteinManager::STRUCTURE pdbInput;
    proteinManager::STRUCTURE pdbOutput;
    
    ////////////////////////////////////////////
    
    proteinManager::geometricTransformations::uniformScaling(pdbInput,0.1);
    
    ////////////////////////////////////////////
    
    pdbInput.loadPDB(pdbInputPath);
    
    proteinManager::ffManager::forceFieldManager ffM;
    
    ffM.loadForceFieldData("../forceFieldManager/forceFieldModels/gromos.ff");
    ffM.applyForceFieldData(pdbInput);
    
    proteinManager::massesManager massesM;
    
    massesM.loadMassesData("../massesManager/massesData/atomMasses.dat");
    massesM.applyMassesData(pdbInput);
    
    cg.applyCoarseGrainedMap<proteinManager::coarseGrainedManager::coarseGrainedMappingSchemes::sasaFitting>(pdbInput,pdbOutput);
    
    ////////////////////////////////////////////
    
    std::ofstream outPutFile;
    outPutFile.open(pdbOutputPath);
    pdbOutput.setOutputFormat(proteinManager::DATA_FORMAT::PDB);
    
    outPutFile << pdbOutput << std::endl;
    
    return EXIT_SUCCESS;
}

