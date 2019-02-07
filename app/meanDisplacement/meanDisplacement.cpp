#include <proteinManager.hpp>

#include "../massesManager/massesManager.hpp"
#include "../centroid/centerOfMass.hpp"

int main(int argc, char* argv[]){
    
    int equiFrames = 300;
    
    proteinManager::STRUCTURE pdb;
    proteinManager::massesManager mM;
    
    mM.loadMassesData("../massesManager/massesData/aminoacidsMasses.dat");
    
    pdb.loadPDB(argv[1]);
    mM.applyMassesData(pdb);
    
    proteinManager::real3 meanPosition = {0,0,0};
    
    int i = 0;
    int framesCount = 0;
    for(proteinManager::MODEL& md : pdb.model()){
	if(i>equiFrames){
	    meanPosition += proteinManager::computeCenterOfMass(md);
	    framesCount ++;
	}
	i++;
    }
    
    meanPosition /= framesCount;
    
    proteinManager::real meanDev = 0;
    i = 0;
    framesCount = 0;
    for(proteinManager::MODEL& md : pdb.model()){
	if(i>equiFrames){
	    proteinManager::real3 dr = proteinManager::computeCenterOfMass(md)-meanPosition;
	    meanDev += sqrt(proteinManager::dot(dr,dr));
	    framesCount ++;
	}
	i++;
    }
    
    std::cout << meanDev/framesCount << std::endl;
}
