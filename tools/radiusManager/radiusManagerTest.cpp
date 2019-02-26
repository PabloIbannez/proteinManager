#include <proteinManager/proteinManager.hpp>

#define DEBUG

#include "radiusManager.hpp"

int main(){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");
    
    proteinManager::radiusManager rM;
    
    rM.loadRadiusData("radiusData/atomRadius.dat");
    rM.applyRadiusData(pdb.model()[0]);
    
    for(proteinManager::MODEL& md : pdb.model()) {
    for(proteinManager::CHAIN& ch : md.chain()) {
    for(proteinManager::RESIDUE& res : ch.residue()) {
    for(proteinManager::ATOM& atm : res.atom()) {
        
        std::cout << atm.getAtomName() << " " << atm.getAtomRadius() << std::endl;
        
    }}}}
}
