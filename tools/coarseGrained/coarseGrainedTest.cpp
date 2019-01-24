#define DEBUG
#include "coarseGrainedManager.hpp"

int main(){
    
    proteinManager::coarseGrainedManager::coarseGrainedGenerator cg;
    
    cg.loadElementsData("./elements/element_Types.prm");
    cg.loadCGmodel("./coarseGrainedModels/MARRINK/aminoAcid2bead_MARRINK.map","./coarseGrainedModels/MARRINK/bead2atom_MARRINK.map");
    
    ////////////////////////////////////////////
    
    std::string pdbInputPath = "./examples/3apgP.pqr";
    std::string pdbOutputPath = "./examples/3apgCG_MARRINK.pqr";
    
    proteinManager::STRUCTURE pdbInput;
    proteinManager::STRUCTURE pdbOutput;
    
    ////////////////////////////////////////////
    
    pdbInput.loadPDB(pdbInputPath);
    
    cg.applyCoarseGrainedMap(pdbInput,pdbOutput);
    
    ////////////////////////////////////////////
    
    std::ofstream outPutFile;
    outPutFile.open(pdbOutputPath);
    pdbOutput.setOuputFormat(proteinManager::DATA_FORMAT::PQR);
    
    outPutFile << pdbOutput << std::endl;
    
    return EXIT_SUCCESS;
}

