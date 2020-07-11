#include <iostream>

#include <proteinManager/proteinManager.hpp>

int main(int argc, char *argv[]) {

    proteinManager::STRUCTURE pdb1;
    
    pdb1.setResSeqOverWrite(true);
    pdb1.setSerialOverWrite(true);
    pdb1.loadPDB(argv[1]);
    
    std::cout << pdb1 << std::endl;
    
    std::cout << pdb1.model(0).chain("A").residue(1).atom("CA").getAtomCoord() << std::endl;
    std::cout << pdb1.model(0).chain("A").residue(1).atom("CA").getChainId() << std::endl;
    
    auto& a = pdb1.model(0).chain("A").residue(1).atom("CA");

    std::cout << a.getChainId() << std::endl;
    pdb1.model(0).chain("A").residue(1).atom("CA").getParentChain().setChainId("Z");
    std::cout << a.getChainId() << std::endl;

    return EXIT_SUCCESS;
}
