#define DEBUG
#include "coarseGrainedManager.hpp"

int main(){
    
    proteinManager::coarseGrainedManager::coarseGrainedGenerator cg;
    
    cg.loadElementsData("./elements/element_Types.prm");
    cg.loadCGmodel("./coarseGrainedModels/RES2BEAD/aminoAcid2bead_RES2BEAD.map","./coarseGrainedModels/RES2BEAD/bead2atom_RES2BEAD.map");
    
    ////////////////////////////////////////////
    
    std::string pdbInputPath = "./examples/streptavidin_1mk5.pqr";
    std::string pdbOutputPath = "./examples/streptavidin_1mk5_CG.pdb";
    
    proteinManager::STRUCTURE pdbInput;
    proteinManager::STRUCTURE pdbOutput;
    
    ////////////////////////////////////////////
    
    pdbInput.loadPDB(pdbInputPath);
    
    cg.applyCoarseGrainedMap(pdbInput,pdbOutput);
    
    ////////////////////////////////////////////
    
    std::ofstream outPutFile;
    outPutFile.open(pdbOutputPath);
    pdbOutput.setOuputFormat(proteinManager::DATA_FORMAT::PDB);
    
    outPutFile << pdbOutput << std::endl;
    
    return EXIT_SUCCESS;
}

