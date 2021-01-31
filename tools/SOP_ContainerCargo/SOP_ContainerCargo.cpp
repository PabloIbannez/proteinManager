#include <proteinManager/proteinManager.hpp>

#include "../../tools/coarseGrained/coarseGrainedCA.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"
#include "../../tools/centers/centroid.hpp"

#include "../../tools/geometric/geometric.hpp"
#include "../../tools/neighbourList/neighbourList.hpp"
#include "../../tools/bonds/bonds.hpp"

#include "../SOP/nativeContactSOP.hpp" 

using namespace proteinManager;

struct bond : public bonds::pair{};

int main(int argc, char *argv[]){
    
    real cutOff=8.0; //angstrom
    
    real eIntraContainer=1.28;
    real eInterContainer=1.05;
    
    real eIntraCargo=1.28;
    real eInterCargo=1.05;

    bool merge_models = true;
    
    STRUCTURE containerIn;
    STRUCTURE cargoIn;
    
    std::string inputContainerFileName  = argv[1];
    std::string inputCargoFileName      = argv[2];
    std::string outputFileName          = argv[3];

    containerIn.loadPDB(inputContainerFileName);
    containerIn.renumber();
    
    cargoIn.loadPDB(inputCargoFileName);
    cargoIn.renumber();

    ////////////////////////////////////////////////
    
    //geometricTransformations::uniformScaling(proteinIn,0.1);

    real3 centroid = computeCentroid(containerIn);

    geometricTransformations::translation(containerIn,real(-1.0)*centroid);
    geometricTransformations::translation(cargoIn,real(-1.0)*centroid);
    
    //geometricTransformations::rotation(proteinIn,centroid,{1,0,0},34.0*(M_PI/180.0));
    //geometricTransformations::rotation(proteinIn,centroid,{0,1,0},13.0*(M_PI/180.0));
    
    ////////////////////////////////////////////////
    
    std::cerr << "Generating CG model" << std::endl;
    STRUCTURE containerOut;
    coarseGrainedCA::coarseGrainedCA_SASA(containerIn,containerOut);
    
    std::cerr << "Generating CG model cargo" << std::endl;
    STRUCTURE cargoOut;
    coarseGrainedCA::coarseGrainedCA_SASA(cargoIn,cargoOut);
    //
    
    int serialOffset;

    for(MODEL&   mdl : containerOut.model()){
    for(CHAIN&   ch  : mdl.chain()       ){
    for(RESIDUE& res : ch.residue()      ){
    for(ATOM&   atm : res.atom()         ){
        serialOffset = atm.getAtomSerial();
    }}}}
    
    for(MODEL&   mdl : cargoOut.model() ){
    for(CHAIN&   ch  : mdl.chain()      ){
    for(RESIDUE& res : ch.residue()     ){
    for(ATOM&   atm : res.atom()        ){
        atm.setAtomSerial(serialOffset+atm.getAtomSerial()+1);
    }}}}

    ////////////////////////////////////////////////
    
    std::cerr << "Generating bonds" << std::endl;

    std::vector<bond> bondVector; 

    bonds::residuePairBonds(containerOut,bondVector);
    bonds::residuePairBonds(cargoOut,bondVector);
    
    ////////////////////////////////////////////////
    
    //Native contacts
    
    //Container
    auto& atomContainerOutVector  = containerOut.atom();

    nativeContactSOP::Parameters parContainer;

    parContainer.eIntra = eIntraContainer;
    parContainer.eInter = eInterContainer;

    parContainer.cutOff = cutOff;

    parContainer.merge_models = merge_models;

    nativeContactSOP  ncContainerSOP(atomContainerOutVector,parContainer);

    auto nc = ncContainerSOP.getNativeContactList();
    
    //Cargo
    auto& atomCargoOutVector  = cargoOut.atom();

    nativeContactSOP::Parameters parCargo;

    parCargo.eIntra = eIntraCargo;
    parCargo.eInter = eInterCargo;

    parCargo.cutOff = cutOff;

    parCargo.merge_models = merge_models;

    nativeContactSOP  ncCargoSOP(atomCargoOutVector,parCargo);

    auto ncCargo = ncCargoSOP.getNativeContactList();

    nc.insert(nc.end(),ncCargo.begin(),ncCargo.end());
    
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
    
    std::cerr << "Writing sp file" << std::endl;
    
    std::map<std::string,int> type2int;
    int typeCount = 0;

    std::ofstream sp(outputFileName+".sp");
    
    for(uint i = 0  ;i<atomContainerOutVector.size();i++){
        
        if(type2int.count(atomContainerOutVector[i].getAtomName())==0){
           type2int[atomContainerOutVector[i].getAtomName()]=typeCount;
           typeCount++;
        }
        
        sp << std::left << std::setw(6) 
           << atomContainerOutVector[i].getAtomCoord() << " " 
           << real(3.8) << " "
           << type2int[atomContainerOutVector[i].getAtomName()] << std::endl; 
    }
    
    for(uint i = 0  ;i<atomCargoOutVector.size();i++){
        
        if(type2int.count(atomCargoOutVector[i].getAtomName())==0){
           type2int[atomCargoOutVector[i].getAtomName()]=typeCount;
           typeCount++;
        }
        
        sp << std::left << std::setw(6) 
           << atomCargoOutVector[i].getAtomCoord() << " " 
           << real(3.8) << " "
           << type2int[atomCargoOutVector[i].getAtomName()] << std::endl; 
    }
    
    ////////////////////////////////////////////////
    
    std::cerr << "Writing coord file" << std::endl;

    std::ofstream coord(outputFileName+".coord");
    
    for(uint i = 0  ;i<atomContainerOutVector.size();i++){
        
        coord << std::left << std::setw(6) 
              << atomContainerOutVector[i].getAtomSerial()<< " " 
              << atomContainerOutVector[i].getAtomCoord() << std::endl;
    }
    
    for(uint i = 0  ;i<atomCargoOutVector.size();i++){
        
        coord << std::left << std::setw(6) 
              << atomCargoOutVector[i].getAtomSerial()    << " " 
              << atomCargoOutVector[i].getAtomCoord()     << std::endl;
    }
    
    ////////////////////////////////////////////////
    
    std::ofstream topology(outputFileName+".top");
    
    ////////////////////////////////////////////////
    
    topology << "[STRUCTURE]" << std::endl;
    int chainCount = 0;
    std::string prevChain = atomContainerOutVector[0].getChainId();
    for(uint i = 0  ;i<atomContainerOutVector.size();i++){
        
        if(prevChain != atomContainerOutVector[i].getChainId()){
            prevChain=atomContainerOutVector[i].getChainId();
            chainCount++;
        }
        
        topology  << std::left << std::setw(6) 
                  << atomContainerOutVector[i].getAtomSerial()  << " " 
                  << std::left << std::setw(6) 
                  << atomContainerOutVector[i].getAtomName()    << " " 
                  << std::left << std::setw(6) 
                  << atomContainerOutVector[i].getResSeq()      << " " 
                  << std::left << std::setw(6) 
                  << chainCount                                 << " " 
                  << std::left << std::setw(6) 
                  << 0                                          << std::endl;
    }
    
    chainCount = 0;
    prevChain = atomCargoOutVector[0].getChainId();

    int modelCount = 1;
    int prevModel  = atomCargoOutVector[0].getModelId();
    
    for(uint i = 0  ;i<atomCargoOutVector.size();i++){
        
        if(prevChain != atomCargoOutVector[i].getChainId()){prevChain=atomCargoOutVector[i].getChainId();
                                                            chainCount++;}
        
        
        if(prevModel != atomCargoOutVector[i].getModelId()){prevModel=atomCargoOutVector[i].getModelId();
                                                            modelCount++;
                                                            chainCount=0;}
        
        topology  << std::left << std::setw(6) 
                  << atomCargoOutVector[i].getAtomSerial()         << " " 
                  << std::left << std::setw(6) 
                  << atomCargoOutVector[i].getAtomName()           << " " 
                  << std::left << std::setw(6) 
                  << atomCargoOutVector[i].getResSeq()             << " " 
                  << std::left << std::setw(6) 
                  << chainCount                                    << " " 
                  << std::left << std::setw(6) 
                  << modelCount                                    << std::endl;
    }
    
    topology << "[SASA]" << std::endl;
    for(uint i = 0 ; i < atomContainerOutVector.size(); i++){
        topology  << std::left << std::setw(6)
                  << atomContainerOutVector[i].getAtomSerial() << " " 
                  << std::fixed << std::setprecision(6)
                  << atomContainerOutVector[i].getAtomSASA()   << std::endl;
    }
    
    for(uint i = 0 ; i < atomCargoOutVector.size(); i++){
        topology  << std::left << std::setw(6)
                  << atomCargoOutVector[i].getAtomSerial() << " " 
                  << std::fixed << std::setprecision(6)
                  << atomCargoOutVector[i].getAtomSASA()   << std::endl;
    }
    
    ////////////////////////////////////////////////
    
    topology << "[SOP_BOND]" << std::endl;
    for(auto& b : bondVector){
        topology << std::left << std::setw(6) 
                 << b.i  << " " 
                 << std::left << std::setw(6) 
                 << b.j  << " "
                 << std::fixed << std::setprecision(6)
                 << b.r0 <<  std::endl;
    }
    ////////////////////////////////////////////////
    
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
    
    ////////////////////////////////////////////////
    
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
