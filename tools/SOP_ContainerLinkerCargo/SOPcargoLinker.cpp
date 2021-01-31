#include <chrono>

#include <proteinManager/proteinManager.hpp>

#include "../../tools/coarseGrained/coarseGrainedCA.hpp"

#include "../../tools/coarseGrained/coarseGrainedManager.hpp"
#include "../../tools/coarseGrained/coarseGrainedMappingSchemes.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"
#include "../../tools/centers/centroid.hpp"

#include "../geometric/geometric.hpp"

#include "../../tools/bonds/bonds.hpp"

#include "../SASA/SASA.hpp"

#include "../SOP/nativeContactSOP.hpp" 

#define cargoInitRes 0
#define cargoEndRes 64

using namespace proteinManager;
struct bond     : public bonds::pair{};
struct angle    : public bonds::angle{};
struct dihedral : public bonds::dihedral{};

using namespace geometric;

#include "cargoLinkerMinimization.hpp"

struct nativeContact{
    int iMol,jMol;
    std::string chi,chj;
    int resi,resj;
    int iSerial,jSerial;
    real r0;
    real e;
};

std::vector<std::string> linkerRes = {"GLU","LEU","PRO",
                                      "TRP","LEU","VAL",
                                      "PRO","ARG","GLY",
                                      "SER","CYS"};
    
bool checkIfExcluded(int serial_i,int serial_j,std::map<int,std::vector<int>>& exclusionsList){

    bool j_in_i = false;
    bool i_in_j = false;

    auto& neig_i = exclusionsList[serial_i];
    if(std::find(neig_i.begin(), neig_i.end(), serial_j) != neig_i.end()){
        j_in_i = true;
    }
    
    auto& neig_j = exclusionsList[serial_j];
    if(std::find(neig_j.begin(), neig_j.end(), serial_i) != neig_j.end()){
        i_in_j = true;
    }

    if(j_in_i and !i_in_j){
        std::cerr << "ERROR in exclusions list, not symmetric" << std::endl;
    }
    
    if(i_in_j and !j_in_i){
        std::cerr << "ERROR in esclusions list, not symmetric" << std::endl;
    }

    return (j_in_i and i_in_j);
}
                                      
int main(int argc, char *argv[]){     
    
    real cutOff=8.0; //angstrom
    
    real eIntraContainer=1.28;
    real eInterContainer=1.05;
    
    real eIntraCargo=1.;
    real eInterCargo=1.;

    bool merge_models = true;
                                      
    STRUCTURE containerIn;              
    STRUCTURE cargoIn;                
    
    std::string inputProteinFileName = argv[1];
    std::string inputCargoFileName   = argv[2];
    
    int nCargo=std::atoi(argv[3]);
    
    std::string outputFileName       = argv[4];

    containerIn.loadPDB(inputProteinFileName);
    containerIn.renumber();

    cargoIn.loadPDB(inputCargoFileName);
    cargoIn.renumber();

    ////////////////////////////////////////////////
    
    //geometricTransformations::uniformScaling(containerIn,0.1);

    real3 centroid = computeCentroid(containerIn);
    
    geometricTransformations::translation(containerIn,real(-1.0)*centroid);
    
    ////////////////////////////////////////////////
    
    std::cerr << "Generating CG model" << std::endl;

    STRUCTURE containerOut;
    coarseGrainedCA::coarseGrainedCA(containerIn,containerOut);
    
    real3 centroidProteinCG = computeCentroid(containerOut);
    geometricTransformations::translation(containerOut,real(-1.0)*centroidProteinCG);
          
    centroidProteinCG = computeCentroid(containerOut);
    
    ///////////////////
    
    std::cerr << "Generating CG model cargo" << std::endl;
    STRUCTURE cargoOut;
    coarseGrainedCA::coarseGrainedCA(cargoIn,cargoOut);
    //coarseGrainedCA::coarseGrainedCA_SASA(cargoIn,cargoOut);
    
    real3 origCargoCG = cargoOut.atom()[0].getAtomCoord();
    geometricTransformations::translation(cargoOut,real(-1.0)*origCargoCG);

    ///////////////////
    
    real containerInnerRadius = INFINITY;
    
    for(MODEL&   mdl : containerOut.model()){
    for(CHAIN&   ch  : mdl.chain()         ){
    for(RESIDUE& res : ch.residue()        ){
    for(ATOM&    atm : res.atom()          ){
        
        real atmCenterDst =  sqrt(dot(atm.getAtomCoord(),atm.getAtomCoord()));
        containerInnerRadius = std::min(containerInnerRadius,atmCenterDst);

    }}}}

    std::cerr<< "Container inner radius: " << containerInnerRadius << std::endl;

    ////////////////////
    
    //Add linker and cargo
    
    STRUCTURE linkersCargos;
   
    int linkerCargoOffSet = 0;

    for(MODEL&   mdl : containerOut.model()){
    for(CHAIN&   ch  : mdl.chain()         ){
    for(RESIDUE& res : ch.residue()        ){
    for(ATOM&   atm : res.atom()           ){
        if(atm.getAtomSerial() > linkerCargoOffSet){
            linkerCargoOffSet = atm.getAtomSerial();
        }
    }}}}

    //offset
    linkerCargoOffSet++;
    
    std::cerr<< "Linker cargo offset: " << linkerCargoOffSet << std::endl;

    /////////////////////////////

    std::random_device rd;
    std::mt19937 rndg(rd());
    //std::mt19937 rndg(1000);
    
    //model chain list for random choose
    struct mdlCh{
        int modelId;
        std::string chainId;
    };
    
    std::vector<mdlCh> allMdlCh;
    
    for(MODEL&   mdl : containerOut.model()){
    for(CHAIN&   ch  : mdl.chain()       ){
        allMdlCh.push_back({mdl.getModelId(),ch.getChainId()});
    }}

    std::shuffle(allMdlCh.begin(), allMdlCh.end(), rndg);
    
    /////////////////////////////

    int linkerCargoCount = 1;
    int linkerCargoAtomCount = linkerCargoOffSet;
    
    std::vector<bond>     linkersBondVector;
    std::vector<angle>    linkersAngleVector;
    std::vector<dihedral> linkersDihedralVector;
    
    for(int cargN=0;cargN<nCargo;cargN++){

        std::cerr << "\r" << "Adding cargo : " << cargN+1 << "/" << nCargo << " Max: " << allMdlCh.size() << std::flush;

        CHAIN& ch = containerOut.model(allMdlCh[cargN].modelId).chain(allMdlCh[cargN].chainId);

        real3 linkerStartPos = ch.residue().front().atom()[0].getAtomCoord();

        real3 linkerVector = centroidProteinCG - linkerStartPos;
              linkerVector = linkerVector/sqrt(dot(linkerVector,linkerVector));

        linkerStartPos = linkerStartPos + real(3.8)*linkerVector;

        //Add linker
        linkersCargos.addModel(linkerCargoCount);
        linkersCargos.model(linkerCargoCount).addChain("L");
       
        //Capsid linker bond
        {
            bond bondBuffer;
            bondBuffer.i  = ch.residue().front().atom()[0].getAtomSerial();
            bondBuffer.j  = linkerCargoAtomCount;
            bondBuffer.r0 = real(3.8);
        
            linkersBondVector.push_back(bondBuffer);
        }
        
        {
            angle angleBuffer;
            angleBuffer.i  = ch.residue().front().atom()[0].getAtomSerial();
            angleBuffer.j  = linkerCargoAtomCount;
            angleBuffer.k  = linkerCargoAtomCount+1;
        
            linkersAngleVector.push_back(angleBuffer);
        }
        
        {
            dihedral dihedralBuffer;
            dihedralBuffer.i  = ch.residue().front().atom()[0].getAtomSerial();
            dihedralBuffer.j  = linkerCargoAtomCount;
            dihedralBuffer.k  = linkerCargoAtomCount+1;
            dihedralBuffer.l  = linkerCargoAtomCount+2;
        
            linkersDihedralVector.push_back(dihedralBuffer);
        }
        
        int linkerResCount = 1;
        for(auto& linkerComp : linkerRes){
            linkersCargos.model(linkerCargoCount).chain("L").addResidue(linkerComp,linkerResCount," ");
            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).addAtom(linkerCargoAtomCount,linkerComp);

            //Old
            //real3 p = linkerStartPos+real(3.8)*(linkerResCount - 1)*linkerVector;
            
            real3 p;
            if(linkerResCount%2==1){
                    p = linkerStartPos+(real(3.8)*(linkerResCount - 1)*linkerVector)/real(2.0); 
            } else {
                    p = linkerStartPos+(real(3.8)*(linkerResCount - 1)*linkerVector)/real(2.0); 
            
                    std::uniform_real_distribution<real> distr(0.0,1.0);

                    real3 rndAxis;

                    real theta = real(2.0)*M_PI*distr(rndg);
                    real phi   = acos(real(1.0)-real(2.0)*distr(rndg));

                    rndAxis.x = sin(phi)*cos(theta);
                    rndAxis.y = sin(phi)*sin(theta);
                    rndAxis.z = cos(phi);

                    real3 u = rndAxis-(dot(rndAxis,linkerVector)*linkerVector);
                          u = u/sqrt(dot(u,u));

                    p = p+u*real(sqrt(3.0)/2.0)*real(3.8);

            }

            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).atom(linkerComp).setAtomCoord(p);
            
            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).atom(linkerComp).setAtomAltLoc(" ");
            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).atom(linkerComp).setAtomOccupancy(1);
            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).atom(linkerComp).setAtomTempFactor(0);
            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).atom(linkerComp).setAtomElement("");
            linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount).atom(linkerComp).setAtomCharge(0);
            
            linkerCargoAtomCount++;
            linkerResCount++;
        }
        
        //Add cargo
        
        //Linker cargo bond
        {
            bond bondBuffer;
            bondBuffer.i  = linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount-1).atom().back().getAtomSerial();
            bondBuffer.j  = linkerCargoAtomCount;
            bondBuffer.r0 = real(3.8);
        
            linkersBondVector.push_back(bondBuffer);
        
        }
        
        {
            angle angleBuffer;
            angleBuffer.i  = linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount-2).atom().back().getAtomSerial();
            angleBuffer.j  = linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount-1).atom().back().getAtomSerial();
            angleBuffer.k  = linkerCargoAtomCount;
        
            linkersAngleVector.push_back(angleBuffer);
        }
        
        {
            dihedral dihedralBuffer;
            dihedralBuffer.i  = linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount-3).atom().back().getAtomSerial();
            dihedralBuffer.j  = linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount-2).atom().back().getAtomSerial();
            dihedralBuffer.k  = linkersCargos.model(linkerCargoCount).chain("L").residue(linkerResCount-1).atom().back().getAtomSerial();
            dihedralBuffer.l  = linkerCargoAtomCount;
        
            linkersDihedralVector.push_back(dihedralBuffer);
        }

        linkersCargos.model(linkerCargoCount).addChain("C");
        
        auto cargoAtoms = cargoOut.atom();

        real3 cargoInit = cargoAtoms[cargoInitRes].getAtomCoord();
        real3 cargoEnd  = cargoAtoms[cargoEndRes].getAtomCoord();

        real3 cargoVector = cargoEnd-cargoInit;
              cargoVector = cargoVector/sqrt(dot(cargoVector,cargoVector));
            
        real  angle = acos(dot(linkerVector,cargoVector));
        real3 axis  = cross(linkerVector,cargoVector);
              axis  = axis/sqrt(dot(axis,axis));
        
        real3 orig = linkerStartPos+(real(3.8)*(linkerResCount - 1)*linkerVector)/real(2.0);
        
        for(ATOM& atm : cargoAtoms){
        
            linkersCargos.model(linkerCargoCount).chain("C").addResidue(atm.getResName(),atm.getResSeq(),atm.getResInsCode());
            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).addAtom(linkerCargoAtomCount,atm.getResName());

            real3 p = atm.getAtomCoord() + orig;

            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).atom(atm.getResName()).setAtomCoord(p);
            
            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).atom(atm.getResName()).setAtomAltLoc(" ");
            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).atom(atm.getResName()).setAtomOccupancy(1);
            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).atom(atm.getResName()).setAtomTempFactor(0);
            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).atom(atm.getResName()).setAtomElement("");
            linkersCargos.model(linkerCargoCount).chain("C").residue(atm.getResSeq()).atom(atm.getResName()).setAtomCharge(0);

            linkerCargoAtomCount++;
        }
        
        geometricTransformations::rotation(linkersCargos.model(linkerCargoCount).chain("C"),orig,axis,-angle);
        
        linkerCargoCount++;
    }

    std::cerr << std::endl;

    //removeCargoClashes(rndg,linkersCargos,centroidProteinCG,containerInnerRadius,real(3.8));
    minimizeCargoLost(rndg,linkersCargos,centroidProteinCG,containerInnerRadius,real(7.0));

    std::ofstream lcFile("linkersCargos.sp");

    std::map<std::string,int> type2int; int typeCount = 0;

    for(MODEL&   mdl : linkersCargos.model()){
    for(CHAIN&   ch  : mdl.chain()          ){
    for(RESIDUE& res : ch.residue()         ){
    for(ATOM&    atm : res.atom()           ){
        
        if(type2int.count(atm.getAtomName())==0){
            type2int[atm.getAtomName()]=typeCount;
            typeCount++;
        }

        lcFile << atm.getAtomCoord()  << " " << 
                  real(3.8)           << " " <<  
                  type2int[atm.getAtomName()] << std::endl;
    }}}}

    //Separate linkers and cargos
    
    STRUCTURE linkers;
    STRUCTURE cargos;

    for(MODEL&   mdl : linkersCargos.model()){
        linkers.addModel(mdl.getModelId());
        cargos.addModel(mdl.getModelId());
    for(CHAIN&   ch  : mdl.chain()       ){

        if        (ch.getChainId() == "L"){
            linkers.model(mdl.getModelId()).addChain(ch.getChainId());
            
            for(RESIDUE& res : ch.residue()      ){
                linkers.model(mdl.getModelId()).chain(ch.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
            for(ATOM&   atm : res.atom()         ){

                linkers.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).addAtom(atm.getAtomSerial(),atm.getAtomName());
                linkers.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).atom(atm.getAtomName()).setAtomCoord(atm.getAtomCoord());

            }}

        } else if (ch.getChainId() == "C"){
            cargos.model(mdl.getModelId()).addChain(ch.getChainId());
            
            for(RESIDUE& res : ch.residue()      ){
                cargos.model(mdl.getModelId()).chain(ch.getChainId()).addResidue(res.getResName(),res.getResSeq(),res.getResInsCode());
            for(ATOM&   atm : res.atom()         ){

                cargos.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).addAtom(atm.getAtomSerial(),atm.getAtomName());
                cargos.model(mdl.getModelId()).chain(ch.getChainId()).residue(res.getResSeq()).atom(atm.getAtomName()).setAtomCoord(atm.getAtomCoord());

            }}

        } else {

        }
    
    }}

    /*
    linkers.setOutputFormat(DATA_FORMAT::SP);
    std::ofstream linkSP("link.sp");
    linkSP << linkers << std::endl;
    
    cargos.setOutputFormat(DATA_FORMAT::SP);
    std::ofstream cargoSP("cargo.sp");
    cargoSP << cargos << std::endl;
    */

    ////////////////////////////////////////////////
    
    std::cerr << "Generating bonds" << std::endl;

    std::vector<bond> bondVector;

    bonds::residuePairBonds(containerOut,bondVector);
    bonds::residuePairBonds(cargos,bondVector);

    bonds::residuePairBonds(linkers,linkersBondVector);
    bonds::residueAngleBonds(linkers,linkersAngleVector);
    bonds::residueDihedralBonds(linkers,linkersDihedralVector);

    ////////////////////////////////////////////////
    
    std::cerr << "Adding SASA" << std::endl;
    
    proteinManager::SASA sasa;

    sasa.addPCASA(containerOut,true);
    sasa.addPCASA(cargos,false);
    //sasa.addSASArndCoil(containerOut);
    //sasa.addSASArndCoil(cargos);
    sasa.addSASArndCoil(linkers);
    
    ////////////////////////////////////////////////
    
    //Native contacts
    
    //Container
    auto& atomContainerOutVector  = containerOut.atom();

    nativeContactSOP::Parameters parContainer;

    parContainer.eIntra = eIntraContainer;
    parContainer.eInter = eInterContainer;

    parContainer.cutOff = cutOff;

    parContainer.merge_models = true;

    nativeContactSOP  ncContainerSOP(atomContainerOutVector,parContainer);

    auto nc = ncContainerSOP.getNativeContactList();
    
    //Cargo
    auto& atomCargosVector  = cargos.atom();

    nativeContactSOP::Parameters parCargo;

    parCargo.eIntra = eIntraCargo;
    parCargo.eInter = eInterCargo;

    parCargo.cutOff = cutOff;

    parCargo.merge_models = false;

    nativeContactSOP  ncCargoSOP(atomCargosVector,parCargo);

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
    
    std::cerr << "Generating exclusion list for linkers"  << std::endl;
    
    for(auto& b : linkersBondVector){
        exclusionsList[b.i].push_back(b.j);
        exclusionsList[b.j].push_back(b.i);
    }
    
    for(auto& a : linkersAngleVector){
        exclusionsList[a.i].push_back(a.j);
        exclusionsList[a.i].push_back(a.k);
        
        exclusionsList[a.j].push_back(a.i);
        exclusionsList[a.j].push_back(a.k);
        
        exclusionsList[a.k].push_back(a.i);
        exclusionsList[a.k].push_back(a.j);
    }
    
    for(auto& d : linkersDihedralVector){
        exclusionsList[d.i].push_back(d.j);
        exclusionsList[d.i].push_back(d.k);
        exclusionsList[d.i].push_back(d.l);
        
        exclusionsList[d.j].push_back(d.i);
        exclusionsList[d.j].push_back(d.k);
        exclusionsList[d.j].push_back(d.l);
        
        exclusionsList[d.k].push_back(d.i);
        exclusionsList[d.k].push_back(d.j);
        exclusionsList[d.k].push_back(d.l);
        
        exclusionsList[d.l].push_back(d.i);
        exclusionsList[d.l].push_back(d.j);
        exclusionsList[d.l].push_back(d.k);
    }

    auto& atomLinkersVector  = linkers.atom();
    
    std::vector<bond> linkersKimHummerVector;
    
    for(uint i=0 ; i<atomLinkersVector.size() ; i++){
        for(uint j=i+1 ; j<atomLinkersVector.size() ; j++){
            int mdl_i = atomLinkersVector[i].getModelId();
            int mdl_j = atomLinkersVector[j].getModelId();

            if(mdl_i == mdl_j){
                int serial_i = atomLinkersVector[i].getAtomSerial();
                int serial_j = atomLinkersVector[j].getAtomSerial();

                if(!checkIfExcluded(serial_i,serial_j,exclusionsList)){
                    bond bondBuffer;
                    bondBuffer.i  = serial_i;
                    bondBuffer.j  = serial_j;
        
                    linkersKimHummerVector.push_back(bondBuffer);
                }
            }
    }}

    
    ////////////////////////////////////////////////
    
    auto& atomLinkersCargosVector = linkersCargos.atom();

    std::cerr << "Writing coord file" << std::endl;

    std::ofstream coord(outputFileName+".coord");
    
    for(uint i = 0  ;i<atomContainerOutVector.size();i++){
        
        coord << std::left << std::setw(6) 
              << atomContainerOutVector[i].getAtomSerial() << " " 
              << atomContainerOutVector[i].getAtomCoord()  << std::endl;
    }
    
    for(uint i = 0  ;i<atomLinkersCargosVector.size();i++){
        
        coord << std::left << std::setw(6) 
              << atomLinkersCargosVector[i].getAtomSerial() << " " 
              << atomLinkersCargosVector[i].getAtomCoord()  << std::endl;
    }
    
    ////////////////////////////////////////////////
    
    std::ofstream topology(outputFileName+".top");
    
    ////////////////////////////////////////////////

    topology << "[STRUCTURE]" << std::endl;

    //std::map<std::string,int> chain2int;
    //int chainCount = 0;
    int chainCount = 1;
    std::string prevChain = atomContainerOutVector[0].getChainId();
    for(uint i = 0  ;i<atomContainerOutVector.size();i++){
        
        if(prevChain != atomContainerOutVector[i].getChainId()){
            prevChain = atomContainerOutVector[i].getChainId();
            chainCount ++;
        }

        topology << std::left << std::setw(6) 
                 << atomContainerOutVector[i].getAtomSerial()         << " " 
                 << std::left << std::setw(6) 
                 << atomContainerOutVector[i].getAtomName()           << " " 
                 << std::left << std::setw(6) 
                 << atomContainerOutVector[i].getResSeq()             << " " 
                 << std::left << std::setw(6) 
                 << chainCount                                        << " "                   
                 << std::left << std::setw(6) 
                 << 0         << std::endl;
    }

    int modelCount = 1;
        chainCount = 1;
    int prevModel = atomCargosVector[0].getModelId();
        prevChain = atomCargosVector[0].getChainId();
    for(uint i = 0  ;i<atomCargosVector.size();i++){
        
        if(prevChain != atomCargosVector[i].getChainId()){
            prevChain = atomCargosVector[i].getChainId();
            chainCount ++;
        }
        
        if(prevModel != atomCargosVector[i].getModelId()){
            prevModel = atomCargosVector[i].getModelId();
            modelCount ++;
            chainCount = 1;
        }

        topology << std::left  << std::setw(6) 
                 << atomCargosVector[i].getAtomSerial() << " " 
                 << std::left  << std::setw(6) 
                 << atomCargosVector[i].getAtomName()   << " " 
                 << std::left  << std::setw(6) 
                 << atomCargosVector[i].getResSeq()     << " " 
                 << std::left  << std::setw(6) 
                 << chainCount                          << " "                   
                 << std::left  << std::setw(6) 
                 << modelCount << std::endl;
    }
    
    chainCount = 1;
    prevModel = atomLinkersVector[0].getModelId();
    prevChain = atomLinkersVector[0].getChainId();
    for(uint i = 0  ;i<atomLinkersVector.size();i++){
        
        if(prevChain != atomLinkersVector[i].getChainId()){
            prevChain = atomLinkersVector[i].getChainId();
            chainCount ++;
        }
        
        if(prevModel != atomLinkersVector[i].getModelId()){
            prevModel = atomLinkersVector[i].getModelId();
            modelCount ++;
            chainCount = 1;
        }

        topology << std::left  << std::setw(6) 
                 << atomLinkersVector[i].getAtomSerial() << " " 
                 << std::left  << std::setw(6) 
                 << atomLinkersVector[i].getAtomName()   << " " 
                 << std::left  << std::setw(6) 
                 << atomLinkersVector[i].getResSeq()     << " " 
                 << std::left  << std::setw(6) 
                 << chainCount                           << " "                   
                 << std::left  << std::setw(6) 
                 << modelCount << std::endl;
    }
    
    topology << "[SASA]" << std::endl;

    //std::map<std::string,int> chain2int;
    //int chainCount = 0;
    for(uint i = 0  ;i<atomContainerOutVector.size();i++){
        
        topology << std::left << std::setw(6) 
                 << atomContainerOutVector[i].getAtomSerial() << " " 
                 << std::left << std::setw(6) 
                 << atomContainerOutVector[i].getAtomSASA()   << std::endl;
    }
    
    for(uint i = 0  ;i<atomCargosVector.size();i++){
        
        topology << std::left << std::setw(6) 
                 << atomCargosVector[i].getAtomSerial()       << " " 
                 << std::left << std::setw(6) 
                 << atomCargosVector[i].getAtomSASA()         << std::endl;
    }
    
    for(uint i = 0  ;i<atomLinkersVector.size();i++){
        
        topology << std::left << std::setw(6) 
                 << atomLinkersVector[i].getAtomSerial()      << " " 
                 << std::left << std::setw(6) 
                 << atomLinkersVector[i].getAtomSASA()        << std::endl;
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
    
    topology << "[SOP_BOND]" << std::endl;

    for(auto& b : bondVector){
        
        topology << std::left << std::setw(6) 
                 << b.i  << " " 
                 << std::left << std::setw(6) 
                 << b.j  << " "
                 << std::fixed << std::setprecision(6)
                 << b.r0 <<  std::endl;
    }
    
    topology << "[LINKER_BOND]" << std::endl;

    for(auto& b : linkersBondVector){
        
        topology << std::left << std::setw(6) 
                 << b.i  << " " 
                 << std::left << std::setw(6) 
                 << b.j  << " "
                 << std::fixed << std::setprecision(6)
                 << b.r0 <<  std::endl;
    }
    
    topology << "[LINKER_ANGLE]" << std::endl;

    for(auto& a : linkersAngleVector){
        
        topology << std::left << std::setw(6) 
                 << a.i  << " " 
                 << std::left << std::setw(6) 
                 << a.j  << " "
                 << std::left << std::setw(6) 
                 << a.k <<  std::endl;
    }
    
    topology << "[LINKER_DIEHEDRAL]" << std::endl;

    for(auto& d : linkersDihedralVector){
        
        topology << std::left << std::setw(6) 
                 << d.i  << " " 
                 << std::left << std::setw(6) 
                 << d.j  << " "
                 << std::left << std::setw(6) 
                 << d.k  << " "
                 << std::left << std::setw(6) 
                 << d.l <<  std::endl;
    }
    
    topology << "[LINKER_KIM_HUMMER]" << std::endl;

    for(auto& b : linkersKimHummerVector){
        
        topology << std::left << std::setw(6) 
                 << b.i  << " " 
                 << std::left << std::setw(6) 
                 << b.j  <<  std::endl;
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
