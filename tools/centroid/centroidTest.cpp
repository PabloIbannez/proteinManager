#include <proteinManager.hpp>
#include "centroid.hpp"
#include "centerOfMass.hpp"

int main(int argc,char *argv[]){
    
    proteinManager::STRUCTURE pdb;
    
    pdb.loadPDB("1aki.pdb");
    
    std::cout << computeCentroid(pdb) << std::endl;
    //std::cout << computeCenterOfMass(pdb) << std::endl;
    
    return EXIT_SUCCESS;
}
