#include <proteinManager/proteinManager.hpp>

#include "../../tools/coarseGrained/coarseGrainedManager.hpp"
#include "../../tools/coarseGrained/coarseGrainedMappingSchemes.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"

//#define DEBUG
#define EXCL 3

using namespace proteinManager;

real beadCharge(std::string const & beadName){

    if(beadName.substr(0,3) == "LYS" or beadName.substr(0,3) == "ARG"){
        return 1.0;
    }
    
    if(beadName.substr(0,3) == "ASP" or beadName.substr(0,3) == "GLU"){
        return -1.0;
    }
    
    if(beadName.substr(0,3) == "HIS" or beadName.substr(0,3) == "HID" or beadName.substr(0,3) == "HIE" or beadName.substr(0,3) == "HIP" ){
        return 0.5;
    }

    return 0.0;
}

void toHIS(RESIDUE& res){
    if(res.getResName().substr(0,3) == "HID" or res.getResName().substr(0,3) == "HIE" or res.getResName().substr(0,3) == "HIP" ){
        res.setResName("HIS");
    }
}

struct KH{
                                
    void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
        
        ////////////////////////////////////////////////
        
        real3 pos = resIn.atom("CA").getAtomCoord();
        
        ////////////////////////////////////////////////
        
        resOut.atom(beadName).setAtomCoord(pos);
        resOut.atom(beadName).setAtomCharge(beadCharge(beadName));

        //Common properties
        resOut.atom(beadName).setAtomAltLoc(" ");
        resOut.atom(beadName).setAtomOccupancy(1);
        resOut.atom(beadName).setAtomTempFactor(0);
        resOut.atom(beadName).setAtomElement("");
    }
                                
};

std::map<std::tuple<std::string,std::string>,std::vector<real>> diehMap(){
    
    std::string line;
    std::stringstream ss;
    
    std::map<std::tuple<std::string,std::string>,std::vector<real>> KBdiehParamMap; 

    std::ifstream KBdiehParamFile("./KaranicolasBrooksDiheParm.dat");
    if(!KBdiehParamFile){
        ss.clear();
        ss << "File not found: " << "./KaranicolasBrooksDiheParm.dat";
        throw std::runtime_error(ss.str());
    }
    
    std::string atm1,atm2;
    int   n;
    real  V;
    real  theta0;

    while(std::getline(KBdiehParamFile,line)){
        
        //Empty lines or lines which starts with # are ignored
        if(line[0]=='#' or line.empty()) continue;
        
        //Process line
        ss.clear();
        ss.str(line);
        ss >> atm1 >> atm2 >> V >> n >> theta0;

        std::tuple<std::string,std::string> key=std::make_tuple(atm1,atm2);

        KBdiehParamMap[key].resize(8);
        KBdiehParamMap[key][n-1]=V;
        KBdiehParamMap[key][n-1+4]=theta0;
    }
    
#ifdef DEBUG
    for(auto& k : KBdiehParamMap){

        std::cout << std::get<0>(k.first) << " " << std::get<1>(k.first) << " " ;

        for(auto& v : k.second){
            if(v==0){
                ss.clear();
                ss << "diehMap: no value added for " << std::get<0>(k.first) << " " << std::get<1>(k.first);
                throw std::runtime_error(ss.str());
            }
            std::cout << v << " ";
        }

        std::cout << std::endl;
    }
#endif
    
    return KBdiehParamMap; 
}

bool isNonHydrogenSideChainAtom(ATOM& atm){
    std::vector<std::string> backboneAtoms = {"CA", "N", "H", "O", "C"};
    std::string atmName = atm.getAtomName();

    if(!(std::find(backboneAtoms.begin(), backboneAtoms.end(), atmName) != backboneAtoms.end())){
        if(atmName[0] != 'H'){
            return true;            
        }
    }

    return false;
}

bool isNonHydrogen(ATOM& atm){
    std::string atmName = atm.getAtomName();

    if(atmName[0] != 'H'){
        return true;            
    }

    return false;
}


real atomsDst(ATOM& atmi,ATOM& atmj){
    real3 dr = atmj.getAtomCoord()-atmi.getAtomCoord();
    return sqrt(dot(dr,dr));
}

RESIDUE& getAllAtomResidue(STRUCTURE& structAllAtom,ATOM& bead){
    return structAllAtom.model(bead.getModelId()).chain(bead.getChainId()).residue(bead.getResSeq());
}

bool areResiduesInContact(RESIDUE& resi,RESIDUE& resj){

    for(ATOM& atmi : resi.atom()){
    for(ATOM& atmj : resj.atom()){
        if(isNonHydrogenSideChainAtom(atmi) and isNonHydrogenSideChainAtom(atmj)){
            if(atomsDst(atmi,atmj)<real(4.5)){
                return true;
            }
        }
    }}

    return false;
}

bool isFirstChainResidue(RESIDUE& res){
    return res.getParentChain().residue()[0].getResSeq()==res.getResSeq();
}

bool isLastChainResidue(RESIDUE& res){
    return res.getParentChain().residue().back().getResSeq()==res.getResSeq();
}

bool areResiduesHydrogenBondend(RESIDUE& resi,RESIDUE& resj){ //Not symmetric !!!
    
    if(resj.getResName()=="PRO"){return false;}
    if(isFirstChainResidue(resj) or isLastChainResidue(resj)){return false;}
    
    real q1 = 0.42;
    real q2 = 0.20;
    real f =   332;
    real F=q1*q2*f;

    real r_ON = atomsDst(resi.atom("O"),resj.atom("N"));
    real r_CH = atomsDst(resi.atom("C"),resj.atom("H"));
    real r_OH = atomsDst(resi.atom("O"),resj.atom("H"));
    real r_CN = atomsDst(resi.atom("C"),resj.atom("N"));

    real e = F*(1.0/r_ON+1.0/r_CH-1.0/r_OH-1.0/r_CN);

    if(e<real(-0.5)){
        return true;
    } else {
        return false;
    }
}

int main(int argc, char *argv[]){
    
    std::map<std::tuple<std::string,std::string>,std::vector<real>> KBdiehMap = diehMap();

    ////////////////////////////////////////////////

    STRUCTURE proteinIn;
    STRUCTURE proteinOut;
    
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2];

    proteinIn.loadPDB(inputFileName);
    proteinIn.renumber();

    ////////////////////////////////////////////////
    
    //geometricTransformations::uniformScaling(proteinIn,0.1);

    ////////////////////////////////////////////////
    
    coarseGrainedManager::coarseGrainedGenerator cg;
    
    cg.loadCGmodel("../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/aminoAcid2bead_RES2BEAD.map","../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/bead2atom_RES2BEAD.map");
    cg.applyCoarseGrainedMap<KH>(proteinIn,proteinOut);
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()       ){
    for(RESIDUE& res : ch.residue()      ){
        toHIS(res);
    }}}
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()       ){
    for(RESIDUE& res : ch.residue()      ){
    for(ATOM&   atm : res.atom()         ){
        atm.setAtomName(atm.getResName());
    }}}}

    //std::cout << proteinOut << std::endl;

    ////////////////////////////////////////////////
    
    std::ofstream bonds(outputFileName+".bond");
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()  ){
    for(uint i=0 ; i<ch.residue().size()-1 ; i++ ){
        
        real3 r1 = ch.residue()[i  ].atom()[0].getAtomCoord();
        real3 r2 = ch.residue()[i+1].atom()[0].getAtomCoord();
        
        real3 r21 = r2-r1;

        real r = sqrt(dot(r21,r21));;

        bonds << ch.residue()[i  ].atom()[0].getAtomSerial() << " " 
              << ch.residue()[i+1].atom()[0].getAtomSerial() << " "
              << r << std::endl;
    }}}
    
    ////////////////////////////////////////////////
    
    std::ofstream angles(outputFileName+".angle");
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()  ){
    for(uint i=0 ; i<ch.residue().size()-2 ; i++){
        
        real3 r1 = ch.residue()[i  ].atom()[0].getAtomCoord();
        real3 r2 = ch.residue()[i+1].atom()[0].getAtomCoord();
        real3 r3 = ch.residue()[i+2].atom()[0].getAtomCoord();

        real3 dr21     = r1-r2;
        real  r21      = sqrt(dot(dr21,dr21));

        dr21/=r21;
        
        real3 dr32     = r3-r2;
        real  r32    = sqrt(dot(dr32,dr32));

        dr32/=r32;

        real theta0 = acos(dot(dr21,dr32))*(180.0/M_PI);
        
        angles << ch.residue()[i  ].atom()[0].getAtomSerial() << " " 
               << ch.residue()[i+1].atom()[0].getAtomSerial() << " "
               << ch.residue()[i+2].atom()[0].getAtomSerial() << " "
               << theta0 << std::endl;
    }}}
    
    ////////////////////////////////////////////////
    
    std::ofstream diehs(outputFileName+".dieh");
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()  ){
    for(uint i=0 ; i<ch.residue().size()-3 ; i++){
        
        std::string atm1 = ch.residue()[i+1].atom()[0].getAtomName();
        std::string atm2 = ch.residue()[i+2].atom()[0].getAtomName();

        atm1=atm1.substr(0,3);
        atm2=atm2.substr(0,3);
        
        if(KBdiehMap.count(std::make_tuple(atm1,atm2))!=1){
            std::stringstream ss;
            ss << "Error generating .dieh file, no diehdral value for " << atm1 << " " << atm2;
            throw std::runtime_error(ss.str());
        } 
        
        std::vector<real> diehValues = KBdiehMap[std::make_tuple(atm1,atm2)];

        diehs << ch.residue()[i  ].atom()[0].getAtomSerial() << " " 
              << ch.residue()[i+1].atom()[0].getAtomSerial() << " "
              << ch.residue()[i+2].atom()[0].getAtomSerial() << " "
              << ch.residue()[i+3].atom()[0].getAtomSerial() << " ";

        for(real& v : diehValues){
            diehs << v << " ";
        }

        diehs << std::endl;

    }}}

    ////////////////////////////////////////////////
    
    //Native contacts
    
    std::map<std::tuple<int,int>,std::string> nc;

    auto& atomOutVector  = proteinOut.atom();

    for(uint i = 0  ;i<atomOutVector.size();i++){
    for(uint j = i+1;j<atomOutVector.size();j++){

        int iSerial = atomOutVector[i].getAtomSerial();
        int jSerial = atomOutVector[j].getAtomSerial();

        RESIDUE& resi = getAllAtomResidue(proteinIn,atomOutVector[i]); 
        RESIDUE& resj = getAllAtomResidue(proteinIn,atomOutVector[j]);

        std::string chi = resi.getChainId();
        std::string chj = resj.getChainId();
        
        int resi_seq = resi.getResSeq();
        int resj_seq = resj.getResSeq();

        int c=0;
        int h=0;

        int nexcl;

        if(chi == chj){
            nexcl = EXCL;
        } else {
            nexcl = 0;
        }
            
        if(abs(resi_seq-resj_seq)>=nexcl){
            if(areResiduesInContact(resi,resj)){
                nc[std::make_tuple(iSerial,jSerial)].append("S");
                c++;
            }

            if(areResiduesHydrogenBondend(resi,resj)){ //ij
                nc[std::make_tuple(iSerial,jSerial)].append("H");
                h++;
            }

            if(areResiduesHydrogenBondend(resj,resi)){ //ji
                nc[std::make_tuple(iSerial,jSerial)].append("H");
                h++;
            }

            if(c+h>=2){

                std::string bType;

                if(h  ==2){bType = "O";}
                if(c+h==2){bType = "O";}
                if(c+h==3){bType = "OO";}

                if(atomOutVector[i].getParentChain().isRes(resi_seq-1)){ //i-1,j
                    if(abs(resi_seq-1-resj_seq)>=nexcl){
                        nc[std::make_tuple(atomOutVector[i].getParentChain().residue(resi_seq-1).atom()[0].getAtomSerial(),jSerial)].append(bType);
                    }
                }
                if(atomOutVector[j].getParentChain().isRes(resj_seq-1)){ //i,j-1
                    if(abs(resi_seq-(resj_seq-1))>=nexcl){
                        nc[std::make_tuple(iSerial,atomOutVector[j].getParentChain().residue(resj_seq-1).atom()[0].getAtomSerial())].append(bType);
                    }
                }

                if(atomOutVector[i].getParentChain().isRes(resi_seq+1)){ //i+1,j
                    if(abs(resi_seq+1-resj_seq)>=nexcl){
                        nc[std::make_tuple(atomOutVector[i].getParentChain().residue(resi_seq+1).atom()[0].getAtomSerial(),jSerial)].append(bType);
                    }
                }
                if(atomOutVector[j].getParentChain().isRes(resj_seq+1)){ //i,j+1
                    if(abs(resi_seq-(resj_seq+1))>=nexcl){
                        nc[std::make_tuple(iSerial,atomOutVector[j].getParentChain().residue(resj_seq+1).atom()[0].getAtomSerial())].append(bType);
                    }
                }
            }
        }
    }}
    
    std::ofstream Qlist(outputFileName+".Qlist");

    for(auto& v : nc){
        int i = std::get<0>(v.first);
        int j = std::get<1>(v.first);
        std::string details = v.second;
        
        Qlist << atomOutVector[i].getAtomName() << " " << atomOutVector[i].getResSeq() << " " 
              << atomOutVector[j].getAtomName() << " " << atomOutVector[j].getResSeq() << " "
              << details    << std::endl;
    }

    ////////////////////////////////////////////////
    
    std::vector<real> repRadius(atomOutVector.size(),std::numeric_limits<real>::max());

    for(uint i = 0  ;i<atomOutVector.size();i++){
    for(uint j = i+1;j<atomOutVector.size();j++){
        
        std::string chi = atomOutVector[i].getChainId();
        std::string chj = atomOutVector[j].getChainId();
        
        int resi_seq = atomOutVector[i].getResSeq();
        int resj_seq = atomOutVector[j].getResSeq();

        int nexcl;

        if(chi == chj){
            nexcl = EXCL+1;
        } else {
            nexcl = 0;
        }
        
        if(abs(resi_seq-resj_seq)>=nexcl){
            if(nc.count(std::make_tuple(atomOutVector[i].getAtomSerial(),atomOutVector[j].getAtomSerial()))==0){
                real r = atomsDst(atomOutVector[i],atomOutVector[j]);
                if(repRadius[i] > r){
                    //std::cout << "i " << atomOutVector[i].getAtomSerial() << " " << atomOutVector[j].getAtomSerial() << " " << chi << " " << chj << " " << atomOutVector[i].getAtomName() << " " << atomOutVector[i].getResSeq() << " " << r << " " << r*std::pow(2.0,1.0/6.0)/2.0 <<std::endl;
                    repRadius[i]=r;
                    //std::cin.get();
                }
                if(repRadius[j] > r){
                    //std::cout << "j " << atomOutVector[i].getAtomSerial() << " " << atomOutVector[j].getAtomSerial() << " " << chi << " " << chj << " " << atomOutVector[j].getAtomName() << " " << atomOutVector[j].getResSeq() << " " << r << " " << r*std::pow(2.0,1.0/6.0)/2.0 <<std::endl;
                    repRadius[j]=r;
                    //std::cin.get();
                }
            }
        }

    }}
    
    std::ofstream wca(outputFileName+".rep");
    
    for(uint i = 0 ; i < repRadius.size(); i++){
        wca << atomOutVector[i].getAtomSerial() << " " << repRadius[i]*std::pow(2.0,1.0/6.0)/2.0 << std::endl;
    }
}
