#include <proteinManager/proteinManager.hpp>

#include "../../tools/coarseGrained/coarseGrainedCA.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"
#include "../../tools/centers/centroid.hpp"

#include "../../tools/geometric/geometric.hpp"
#include "../../tools/neighbourList/neighbourList.hpp"
#include "../../tools/bonds/bonds.hpp"

#include "./nativeContactSOP.hpp" 

using namespace proteinManager;

struct bond : public bonds::pair{};

int main(int argc, char *argv[]){

    real cutOff=8.0; //angstrom
    
    real eIntra=1.26;
    real eInter=1.1;

    bool merge_models = true;
    
    STRUCTURE proteinIn;
    
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2];

    proteinIn.loadPDB(inputFileName);
    proteinIn.renumber();

    ////////////////////////////////////////////////
    
    real3 centroid = computeCentroid(proteinIn);
    //geometricTransformations::translation(proteinIn,real(-1.0)*centroid);
    //geometricTransformations::rotation(proteinIn,centroid,{1,0,0},34.0*(M_PI/180.0));
    //geometricTransformations::rotation(proteinIn,centroid,{0,1,0},13.0*(M_PI/180.0));
    
    ////////////////////////////////////////////////
    
    std::cerr << "Generating CG model" << std::endl;

    STRUCTURE proteinOut;
    
    coarseGrainedCA::coarseGrainedCA(proteinIn,proteinOut);

    ////////////////////////////////////////////////
    
    std::cerr << "Generating bonds" << std::endl;

    std::vector<bond> bondVector; 

    bonds::residuePairBonds(proteinOut,bondVector);
    
    ////////////////////////////////////////////////
    
    auto& atomOutVector  = proteinOut.atom();

    nativeContactSOP::Parameters par;

    par.eIntra = eIntra;
    par.eInter = eInter;

    par.cutOff = cutOff;

    par.merge_models = merge_models;

    nativeContactSOP  ncSOP(atomOutVector,par);

    auto nc = ncSOP.getNativeContactList();
    
    ////////////////////////////////////////////////

    std::cerr << "Generating exclusion list"  << std::endl;

    std::map<int,std::vector<int>> exclusionsList;

    for(auto& b : bondVector){
        exclusionsList[b.i].push_back(b.j);
        exclusionsList[b.j].push_back(b.i);
    }
    
    for(auto& c : nc){
        exclusionsList[c.iSerial].push_back(c.jSerial);
        exclusionsList[c.jSerial].push_back(c.iSerial);
    }

    ////////////////////////////////////////////////
    
    std::cerr << "Writing coord file" << std::endl;

    std::ofstream coord(outputFileName+".coord");
    
    for(uint i = 0  ;i<atomOutVector.size();i++){
        
        coord << std::left << std::setw(6) 
              << atomOutVector[i].getAtomSerial() << " " 
              << atomOutVector[i].getAtomCoord()  << std::endl;
    }
    
    ////////////////////////////////////////////////
    
    std::ofstream topology(outputFileName+".top");
    
    ////////////////////////////////////////////////

    topology << "[STRUCTURE]" << std::endl;

    //std::map<std::string,int> chain2int;
    //int chainCount = 0;
    int chainCount = 1;
    std::string prevChain = atomOutVector[0].getChainId();
    for(uint i = 0  ;i<atomOutVector.size();i++){
        
        if(prevChain != atomOutVector[i].getChainId()){
            prevChain = atomOutVector[i].getChainId();
            chainCount ++;
        }

        /*
        if(chain2int.count(atomOutVector[i].getChainId())==0){
           chain2int[atomOutVector[i].getChainId()]=chainCount;
           chainCount++;
        }*/

        topology << std::left << std::setw(6) 
                 << atomOutVector[i].getAtomSerial()         << " " 
                 << std::left << std::setw(6) 
                 << atomOutVector[i].getAtomName()           << " " 
                 << std::left << std::setw(6) 
                 << atomOutVector[i].getResSeq()             << " " 
                 << std::left << std::setw(6) 
                 //<< chain2int[atomOutVector[i].getChainId()] << " " 
                 << chainCount << " " ;
                 if(merge_models){
                    topology << std::left << std::setw(6) 
                             << 0         << std::endl;
                 } else {
                    topology << std::left << std::setw(6) 
                             << atomOutVector[i].getModelId()            
                             << std::endl; 
                 }
    }
    
    topology << "[SOP_BOND]" << std::endl;

    for(auto& b : bondVector){
        
        topology << std::left << std::setw(6) 
                 << b.i  << " " 
                 << std::left << std::setw(6) 
                 << b.j  << " "
                 << std::fixed << std::setprecision(6)
                 << b.r0 <<  std::endl;
    }
    
    topology << "[SOP_NATIVE_CONTACT]" << std::endl;
    for(auto& c : nc){

        topology << std::left  << std::setw(6)
                 << c.iSerial  << " " 
                 << std::left  << std::setw(6) 
                 << c.jSerial  << " "
                 << std::fixed << std::setprecision(6)
                 << c.r0       << " "
                 << std::fixed << std::setprecision(6)
                 << c.e        << std::endl;
        
    }                              
    
    uint nExclusions=0;
    uint maxExclusion=0;
    for(auto& E : exclusionsList){
        std::sort(E.second.begin(),E.second.end());
        E.second.erase(std::unique(E.second.begin(),E.second.end()),E.second.end());
        nExclusions+=E.second.size();
        if(E.second.size()>maxExclusion){
            maxExclusion=E.second.size();
        }
    }
    
    topology << "[SOP_EXCLUSIONS]" << std::endl;
    for(auto& E : exclusionsList){
        topology << E.first << " ";
        
        for(auto& e : E.second){
            topology << std::left << std::setw(6) << e << " ";
        }

        topology << std::endl;
    }
    
}
