#include "SASA.hpp"
#include "../../tools/coarseGrained/coarseGrainedCA.hpp"

int main(){

    proteinManager::STRUCTURE pdbIn;

    pdbIn.loadPDB("1egl.pdb");
    
    proteinManager::STRUCTURE pdbOut;
    proteinManager::coarseGrainedCA::coarseGrainedCA(pdbIn,pdbOut);

    proteinManager::SASA sasa;

    //sasa.addSASArndCoil(pdbOut);
    sasa.addPCASA(pdbOut);

    std::ofstream outFile("1egl.test");
    
    for(proteinManager::MODEL&   mdl : pdbOut.model()){
    for(proteinManager::CHAIN&   ch  : mdl.chain()   ){
    for(proteinManager::RESIDUE& res : ch.residue()  ){
    for(proteinManager::ATOM&    atm : res.atom()    ){

        outFile << atm.getAtomSerial() << " "
                << atm.getAtomName()   << " "
                << atm.getAtomSASA()   << std::endl;

    }}}}
}
