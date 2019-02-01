#include <iostream>

#include <proteinManager.hpp>

int main(int argc, char *argv[]) {

    proteinManager::STRUCTURE pdb1;
    
    pdb1.setResSeqOverWrite(true);
    pdb1.setSerialOverWrite(true);
    pdb1.loadPDB(argv[1]);

    std::cout << pdb1 << std::endl;
    
    return EXIT_SUCCESS;
}
