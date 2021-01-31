#include <proteinManager/proteinManager.hpp>


int main(int argc, char* argv[]){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB(argv[1]);

}
