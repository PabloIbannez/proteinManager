#include <proteinManager/proteinManager.hpp>

#include "../../tools/coarseGrained/coarseGrainedManager.hpp"
#include "../../tools/coarseGrained/coarseGrainedMappingSchemes.hpp"

#include "../../tools/geometricTransformations/geometricTransformations.hpp"
#include "../../tools/centers/centroid.hpp"

//#define DEBUG

#define MERGE_MODELS 

#define EXCL 3
#define Tf 350.0
#define cutOff 25 //angstrom

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

struct CA{
                                
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

struct CA_SASA{
                                
    void mappingScheme(RESIDUE& resIn, RESIDUE& resOut, std::string const & beadName,std::vector<std::string>& beadComponents){
        
        ////////////////////////////////////////////////
        
        real3 pos = resIn.atom("CA").getAtomCoord();
        
        ////////////////////////////////////////////////
        
        resOut.atom(beadName).setAtomCoord(pos);
        resOut.atom(beadName).setAtomCharge(beadCharge(beadName));
        
        real SASA = 0;
        for(ATOM& atm : resIn.atom()){
            SASA+=atm.getAtomSASA();
        }

        resOut.atom(beadName).setAtomSASA(SASA);

        //Common properties
        resOut.atom(beadName).setAtomAltLoc(" ");
        resOut.atom(beadName).setAtomOccupancy(1);
        resOut.atom(beadName).setAtomTempFactor(0);
        resOut.atom(beadName).setAtomElement("");
    }
                                
};

std::map<std::tuple<std::string,std::string>,std::vector<real>> diheMap(){
    
    std::string line;
    std::stringstream ss;
    
    std::map<std::tuple<std::string,std::string>,std::vector<real>> KBdiheParamMap; 

    std::ifstream KBdiheParamFile("../../data/KaranicolasBrooksDiheParm.dat");
    if(!KBdiheParamFile){
        ss.clear();
        ss << "File not found: " << "../../data/KaranicolasBrooksDiheParm.dat";
        throw std::runtime_error(ss.str());
    }
    
    std::string atm1,atm2;
    int   n;
    real  V;
    real  theta0;

    while(std::getline(KBdiheParamFile,line)){
        
        //Empty lines or lines which starts with # are ignored
        if(line[0]=='#' or line.empty()) continue;
        
        //Process line
        ss.clear();
        ss.str(line);
        ss >> atm1 >> atm2 >> V >> n >> theta0;

        std::tuple<std::string,std::string> key=std::make_tuple(atm1,atm2);

        KBdiheParamMap[key].resize(8);
        KBdiheParamMap[key][n-1]=V;
        KBdiheParamMap[key][n-1+4]=theta0;
    }
    
#ifdef DEBUG
    for(auto& k : KBdiheParamMap){

        std::cout << std::get<0>(k.first) << " " << std::get<1>(k.first) << " " ;

        for(auto& v : k.second){
            if(v==0){
                ss.clear();
                ss << "diheMap: no value added for " << std::get<0>(k.first) << " " << std::get<1>(k.first);
                throw std::runtime_error(ss.str());
            }
            std::cout << v << " ";
        }

        std::cout << std::endl;
    }
#endif
    
    return KBdiheParamMap; 
}

std::map<std::tuple<std::string,std::string>,real> interactionEnergyMap(){
    
    std::string line;
    std::stringstream ss;
    
    std::map<std::tuple<std::string,std::string>,real> interactionMap; 

    std::ifstream interactionParamFile("../../data/MiyazawaJernigaEnergyParm.dat");
    if(!interactionParamFile){
        ss.clear();
        ss << "File not found: " << "../../data/MiyazawaJernigaEnergyParm.dat";
        throw std::runtime_error(ss.str());
    }
    
    std::string atm1,atm2;
    real  e;
    
    real E=0;
    int  eValues=0;

    while(std::getline(interactionParamFile,line)){
        
        //Empty lines or lines which starts with # are ignored
        if(line[0]=='#' or line.empty()) continue;
        
        //Process line
        ss.clear();
        ss.str(line);
        ss >> atm1 >> atm2 >> e;

        std::tuple<std::string,std::string> key=std::make_tuple(atm1,atm2);
        interactionMap[key]=e;
        E+=e;
        eValues++;

        if(atm1!=atm2){
            std::tuple<std::string,std::string> key=std::make_tuple(atm2,atm1);
            interactionMap[key]=e;
        }
    }
    
#ifdef DEBUG
    for(auto& e : interactionMap){
        std::cout << std::get<0>(e.first) << " " << std::get<1>(e.first) << " " << e.second << std::endl;
    }
#endif

    //Renormalization
    E=E/eValues;
    for(auto& e : interactionMap){
        e.second=real(0.5)*e.second/E;
    }
    
    return interactionMap; 
}

real atomsDst(ATOM& atmi,ATOM& atmj){
    real3 dr = atmj.getAtomCoord()-atmi.getAtomCoord();
    return sqrt(dot(dr,dr));
}

real3 cross(real3& ri,real3& rj){
    return {ri.y*rj.z - ri.z*rj.y, (-ri.x*rj.z + ri.z*rj.x), ri.x*rj.y - ri.y*rj.x};
}

real sign(real x){
    if(x>real(0)){return real( 1);}
    if(x<real(0)){return real(-1);}
    return real(0);
}

real dihedral(ATOM& atmi,ATOM& atmj,ATOM& atmk,ATOM& atml){
    real3 ri=atmi.getAtomCoord();
    real3 rj=atmj.getAtomCoord();
    real3 rk=atmk.getAtomCoord();
    real3 rl=atml.getAtomCoord();

    real3 rij=rj-ri;
    real3 rjk=rk-rj;
    real3 rkl=rl-rk;
    
    real3 va=cross(rij,rjk);
    real3 vb=cross(rjk,rkl);
    
    real vam = sqrt(dot(va,va));
    real vbm = sqrt(dot(vb,vb));
    
    va/=vam;
    vb/=vbm;

    real s = dot(cross(rij,rjk),rkl);
         s = sign(s);

    return s*acos(dot(va,vb))*180/M_PI;
}

RESIDUE& getAllAtomResidue(STRUCTURE& structAllAtom,ATOM& bead){
    return structAllAtom.model(bead.getModelId()).chain(bead.getChainId()).residue(bead.getResSeq());
}

bool isNonHydrogen(ATOM& atm){
    std::string atmName = atm.getAtomName();

    if(atmName[0] != 'H'){
        return true;            
    }

    return false;
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


bool isFirstChainResidue(RESIDUE& res){
    return res.getParentChain().residue()[0].getResSeq()==res.getResSeq();
}

bool isLastChainResidue(RESIDUE& res){
    return res.getParentChain().residue().back().getResSeq()==res.getResSeq();
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

template<class atomVector>
std::vector<std::vector<int>> generateNeighbourList(atomVector& atmV,real rCut){

    std::vector<std::vector<int>> neigList(atmV.size());
    
    //std::cerr << "Generating neighbour list (parallel) ..." << std::endl;
    #pragma omp parallel for schedule(dynamic)
    for(uint i = 0  ;i<atmV.size();i++){
    if(omp_get_thread_num()==0){
        std::cerr << "\r" << "Generating neighbour list (parallel) " << i+1 << "/" << atmV.size()  << " (number of threads: " << omp_get_num_threads() << ")"<<std::flush;
    }
    for(uint j = i+1;j<atmV.size();j++){
        if(atomsDst(atmV[i],atmV[j])<cutOff){
            neigList[i].push_back(j);
        }
    }}
    std::cerr << "\r\e[K" << "Generating neighbour list " << atmV.size() << "/" << atmV.size() << std::endl;

    return neigList;
}

struct bond{
    int i,j;
    real r0;
};

struct angle{
    int i,j,k;
    real theta0;
};

struct dihe{
    int i,j,k,m;
    real K1,K2,K3,K4;
    real a1,a2,a3,a4;
};

struct nativeContact{
    ATOM* atmi;
    ATOM* atmj;
    std::string type;
    real e=0;
    real eScaled=0;
};

int main(int argc, char *argv[]){
    
    std::map<std::tuple<std::string,std::string>,std::vector<real>> KBdiheMap = diheMap();
    std::map<std::tuple<std::string,std::string>,real>              intEneMap = interactionEnergyMap();

    ////////////////////////////////////////////////

    STRUCTURE proteinIn;
    STRUCTURE proteinOut;
    
    std::string inputFileName  = argv[1];
    std::string outputFileName = argv[2];

    proteinIn.loadPDB(inputFileName);
    proteinIn.renumber();

    ////////////////////////////////////////////////
    
    //geometricTransformations::uniformScaling(proteinIn,0.1);

    real3 centroid = computeCentroid(proteinIn);
    geometricTransformations::translation(proteinIn,real(-1.0)*centroid);
    
    ////////////////////////////////////////////////
    
    std::cerr << "Generating CG model" << std::endl;

    coarseGrainedManager::coarseGrainedGenerator cg;
    
    cg.loadCGmodel("../../../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/aminoAcid2bead_RES2BEAD.map","../../../../tools/coarseGrained/coarseGrainedModels/RES2BEAD/bead2atom_RES2BEAD.map");
    //cg.applyCoarseGrainedMap<CA>(proteinIn,proteinOut);
    cg.applyCoarseGrainedMap<CA_SASA>(proteinIn,proteinOut);
    
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
    
    std::cerr << "Generating bonds" << std::endl;

    std::vector<bond> bondVector;
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()  ){
    for(uint i=0 ; i<ch.residue().size()-1 ; i++ ){
        
        real3 r1 = ch.residue()[i  ].atom()[0].getAtomCoord();
        real3 r2 = ch.residue()[i+1].atom()[0].getAtomCoord();
        
        real3 r21 = r2-r1;

        real r = sqrt(dot(r21,r21));;

        bond bondBuffer;
        bondBuffer.i = ch.residue()[i  ].atom()[0].getAtomSerial();
        bondBuffer.j = ch.residue()[i+1].atom()[0].getAtomSerial();
        bondBuffer.r0=r;

        bondVector.push_back(bondBuffer);

    }}}
    
    ////////////////////////////////////////////////
    
    std::cerr << "Generating angular bonds" << std::endl;
    
    std::vector<angle> angleVector;
    
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

        real theta0 = acos(dot(dr21,dr32));
        
        angle angleBuffer;

        angleBuffer.i      = ch.residue()[i  ].atom()[0].getAtomSerial();
        angleBuffer.j      = ch.residue()[i+1].atom()[0].getAtomSerial();
        angleBuffer.k      = ch.residue()[i+2].atom()[0].getAtomSerial();
        angleBuffer.theta0 = theta0;

        angleVector.push_back(angleBuffer);
    }}}
    
    ////////////////////////////////////////////////
    
    std::cerr << "Generating dihedral" << std::endl;
    
    std::vector<dihe> diheVector;
    
    for(MODEL&   mdl : proteinOut.model()){
    for(CHAIN&   ch  : mdl.chain()  ){
    for(uint i=0 ; i<ch.residue().size()-3 ; i++){
        
        std::string atm1 = ch.residue()[i+1].atom()[0].getAtomName();
        std::string atm2 = ch.residue()[i+2].atom()[0].getAtomName();

        atm1=atm1.substr(0,3);
        atm2=atm2.substr(0,3);
        
        if(KBdiheMap.count(std::make_tuple(atm1,atm2))!=1){
            std::stringstream ss;
            ss << "Error generating .dihe file, no dihedral value for " << atm1 << " " << atm2;
            throw std::runtime_error(ss.str());
        } 
        
        std::vector<real> diheValues = KBdiheMap[std::make_tuple(atm1,atm2)];
        
        dihe diheBuffer;

        diheBuffer.i = ch.residue()[i  ].atom()[0].getAtomSerial();
        diheBuffer.j = ch.residue()[i+1].atom()[0].getAtomSerial();
        diheBuffer.k = ch.residue()[i+2].atom()[0].getAtomSerial();
        diheBuffer.m = ch.residue()[i+3].atom()[0].getAtomSerial();
        
        diheBuffer.K1 = diheValues[0];
        diheBuffer.K2 = diheValues[1];
        diheBuffer.K3 = diheValues[2];
        diheBuffer.K4 = diheValues[3];
        
        diheBuffer.a1 = diheValues[4]*real(M_PI/180);
        diheBuffer.a2 = diheValues[5]*real(M_PI/180);
        diheBuffer.a3 = diheValues[6]*real(M_PI/180);
        diheBuffer.a4 = diheValues[7]*real(M_PI/180);

        diheVector.push_back(diheBuffer);
        
    }}}

    ////////////////////////////////////////////////
    
    //Native contacts
    
    std::map<std::tuple<int,int>,nativeContact> nc;

    auto& atomOutVector  = proteinOut.atom();

    auto nL = generateNeighbourList(atomOutVector,cutOff);

    for(uint i = 0  ;i<atomOutVector.size();i++){
    std::cerr << "\r" << "Generating native contacts " << i+1 << "/" << atomOutVector.size()  << std::flush;
    //for(uint j = i+1;j<atomOutVector.size();j++){
    for(uint j : nL[i]){

        #ifdef MERGE_MODELS
            
        #else
            if(atomOutVector[i].getModelId()!=atomOutVector[j].getModelId()){continue;}
        #endif

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
                nc[std::make_tuple(iSerial,jSerial)].atmi=&atomOutVector[i];
                nc[std::make_tuple(iSerial,jSerial)].atmj=&atomOutVector[j];
                nc[std::make_tuple(iSerial,jSerial)].type.append("S");
                nc[std::make_tuple(iSerial,jSerial)].e+=intEneMap[std::make_tuple(atomOutVector[i].getAtomName(),atomOutVector[j].getAtomName())];
                c++;
            }

            if(areResiduesHydrogenBondend(resi,resj)){ //ij
                nc[std::make_tuple(iSerial,jSerial)].atmi=&atomOutVector[i];
                nc[std::make_tuple(iSerial,jSerial)].atmj=&atomOutVector[j];
                nc[std::make_tuple(iSerial,jSerial)].type.append("H");
                h++;
            }

            if(areResiduesHydrogenBondend(resj,resi)){ //ji
                nc[std::make_tuple(iSerial,jSerial)].atmi=&atomOutVector[i];
                nc[std::make_tuple(iSerial,jSerial)].atmj=&atomOutVector[j];
                nc[std::make_tuple(iSerial,jSerial)].type.append("H");
                h++;
            }

            if(c==0 and h>0){
                nc[std::make_tuple(iSerial,jSerial)].e+=real(1.0);
            }

            if(c+h>=2){

                std::string bType;
                int         nadd;

                if(h  ==2){bType = "O";nadd=1;}
                if(c+h==2){bType = "O";nadd=1;}
                if(c+h==3){bType = "OO";nadd=2;}

                int orientationalHydrogenBondCount = 0;
                if(atomOutVector[i].getParentChain().isRes(resi_seq-1)){ //i-1,j
                    if(abs(resi_seq-1-resj_seq)>=nexcl){
                        int i_nei = atomOutVector[i].getParentChain().residue(resi_seq-1).atom()[0].getAtomSerial();
                        int j_nei = jSerial;
                        nc[std::make_tuple(i_nei,j_nei)].atmi=&atomOutVector[i].getParentChain().residue(resi_seq-1).atom()[0];
                        nc[std::make_tuple(i_nei,j_nei)].atmj=&atomOutVector[j];
                        orientationalHydrogenBondCount++;
                    }
                }
                if(atomOutVector[j].getParentChain().isRes(resj_seq-1)){ //i,j-1
                    if(abs(resi_seq-(resj_seq-1))>=nexcl){
                        int i_nei = iSerial;
                        int j_nei = atomOutVector[j].getParentChain().residue(resj_seq-1).atom()[0].getAtomSerial();
                        nc[std::make_tuple(i_nei,j_nei)].atmi=&atomOutVector[i];
                        nc[std::make_tuple(i_nei,j_nei)].atmj=&atomOutVector[j].getParentChain().residue(resj_seq-1).atom()[0];
                        orientationalHydrogenBondCount++;
                    }
                }

                if(atomOutVector[i].getParentChain().isRes(resi_seq+1)){ //i+1,j
                    if(abs(resi_seq+1-resj_seq)>=nexcl){
                        int i_nei = atomOutVector[i].getParentChain().residue(resi_seq+1).atom()[0].getAtomSerial();
                        int j_nei = jSerial;
                        nc[std::make_tuple(i_nei,j_nei)].atmi=&atomOutVector[i].getParentChain().residue(resi_seq+1).atom()[0];
                        nc[std::make_tuple(i_nei,j_nei)].atmj=&atomOutVector[j];
                        orientationalHydrogenBondCount++;
                    }
                }
                if(atomOutVector[j].getParentChain().isRes(resj_seq+1)){ //i,j+1
                    if(abs(resi_seq-(resj_seq+1))>=nexcl){
                        int i_nei = iSerial;
                        int j_nei = atomOutVector[j].getParentChain().residue(resj_seq+1).atom()[0].getAtomSerial();
                        nc[std::make_tuple(i_nei,j_nei)].atmi=&atomOutVector[i];
                        nc[std::make_tuple(i_nei,j_nei)].atmj=&atomOutVector[j].getParentChain().residue(resj_seq+1).atom()[0];
                        orientationalHydrogenBondCount++;
                    }
                }

                real oe=real(nadd)/orientationalHydrogenBondCount;

                if(atomOutVector[i].getParentChain().isRes(resi_seq-1)){ //i-1,j
                    if(abs(resi_seq-1-resj_seq)>=nexcl){
                        int i_nei = atomOutVector[i].getParentChain().residue(resi_seq-1).atom()[0].getAtomSerial();
                        int j_nei = jSerial;
                        nc[std::make_tuple(i_nei,j_nei)].type.append(bType);
                        nc[std::make_tuple(i_nei,j_nei)].e+=oe;
                    }
                }
                if(atomOutVector[j].getParentChain().isRes(resj_seq-1)){ //i,j-1
                    if(abs(resi_seq-(resj_seq-1))>=nexcl){
                        int i_nei = iSerial;
                        int j_nei = atomOutVector[j].getParentChain().residue(resj_seq-1).atom()[0].getAtomSerial();
                        nc[std::make_tuple(i_nei,j_nei)].type.append(bType);
                        nc[std::make_tuple(i_nei,j_nei)].e+=oe;
                    }
                }

                if(atomOutVector[i].getParentChain().isRes(resi_seq+1)){ //i+1,j
                    if(abs(resi_seq+1-resj_seq)>=nexcl){
                        int i_nei = atomOutVector[i].getParentChain().residue(resi_seq+1).atom()[0].getAtomSerial();
                        int j_nei = jSerial;
                        nc[std::make_tuple(i_nei,j_nei)].type.append(bType);
                        nc[std::make_tuple(i_nei,j_nei)].e+=oe;
                    }
                }
                if(atomOutVector[j].getParentChain().isRes(resj_seq+1)){ //i,j+1
                    if(abs(resi_seq-(resj_seq+1))>=nexcl){
                        int i_nei = iSerial;
                        int j_nei = atomOutVector[j].getParentChain().residue(resj_seq+1).atom()[0].getAtomSerial();
                        nc[std::make_tuple(i_nei,j_nei)].type.append(bType);
                        nc[std::make_tuple(i_nei,j_nei)].e+=oe;
                    }
                }




            }
        }
    }}


    for(auto& v : nc){
        if(std::get<0>(v.first) > std::get<1>(v.first)){
            std::stringstream ss;
            ss << "Error in native contact list, this special case has to be checked";
            throw std::runtime_error(ss.str());
        }
    }

    std::cerr << std::endl;

    ////////////////////////////////////////////////
    
    std::vector<real> innerRadius(atomOutVector.size(),std::numeric_limits<real>::max());

    for(uint i = 0  ;i<atomOutVector.size();i++){
    std::cerr << "\r" << "Generating inner radius " << i+1 << "/" << atomOutVector.size()  << std::flush;
    //for(uint j = i+1;j<atomOutVector.size();j++){
    for(uint j : nL[i]){
        
        #ifdef MERGE_MODELS
            
        #else
            if(atomOutVector[i].getModelId()!=atomOutVector[j].getModelId()){continue;}
        #endif
        
        std::string chi = atomOutVector[i].getChainId();
        std::string chj = atomOutVector[j].getChainId();
        
        int resi_seq = atomOutVector[i].getResSeq();
        int resj_seq = atomOutVector[j].getResSeq();

        int nexcl;

        if(chi == chj){
            nexcl = EXCL+1; //!!!!!!!!!!!!!!
        } else {
            nexcl = 0;
        }
        
        if(abs(resi_seq-resj_seq)>=nexcl){
            if(nc.count(std::make_tuple(atomOutVector[i].getAtomSerial(),atomOutVector[j].getAtomSerial()))==0){
                real r = atomsDst(atomOutVector[i],atomOutVector[j]);
                if(innerRadius[i] > r){
                    innerRadius[i]=r;
                }
                if(innerRadius[j] > r){
                    innerRadius[j]=r;
                }
            }
        }
    }}

    for(real& r : innerRadius){
        r/=real(2.0);
        r*=std::pow(2.0,1.0/6.0);
    }
    
    std::cerr << std::endl;

    ////////////////////////////////////////////////

    std::cerr << "Generating exclusion list"  << std::endl;

    std::map<int,std::vector<int>> exclusionsList;

    for(auto& b : bondVector){
        exclusionsList[b.i].push_back(b.j);
        exclusionsList[b.j].push_back(b.i);
    }
    
    for(auto& a : angleVector){
        exclusionsList[a.i].push_back(a.j);
        exclusionsList[a.i].push_back(a.k);
        
        exclusionsList[a.j].push_back(a.i);
        exclusionsList[a.j].push_back(a.k);
        
        exclusionsList[a.k].push_back(a.i);
        exclusionsList[a.k].push_back(a.j);
    }
    
    for(auto& d : diheVector){
        exclusionsList[d.i].push_back(d.j);
        exclusionsList[d.i].push_back(d.k);
        exclusionsList[d.i].push_back(d.m);
        
        exclusionsList[d.j].push_back(d.i);
        exclusionsList[d.j].push_back(d.k);
        exclusionsList[d.j].push_back(d.m);
        
        exclusionsList[d.k].push_back(d.i);
        exclusionsList[d.k].push_back(d.j);
        exclusionsList[d.k].push_back(d.m);
        
        exclusionsList[d.m].push_back(d.i);
        exclusionsList[d.m].push_back(d.j);
        exclusionsList[d.m].push_back(d.k);
    }
    
    for(auto& v : nc){
        exclusionsList[v.second.atmi->getAtomSerial()].push_back(v.second.atmj->getAtomSerial());
        exclusionsList[v.second.atmj->getAtomSerial()].push_back(v.second.atmi->getAtomSerial());
    }

    ////////////////////////////////////////////////
    
    real epsRes = Tf*0.0054;

    #ifdef MERGE_MODELS
        
        real totalNCenergy=0;
        for(auto& v : nc){
            totalNCenergy+=v.second.e;
        }
    
        int N = atomOutVector.size();
        real epsNC = epsRes*N/totalNCenergy;
        
        for(auto& v : nc){
            v.second.eScaled=v.second.e*epsNC;
        }

    #else
        
        std::map<int,real> totalNCenergy;
        for(MODEL& md : proteinOut.model()){
            totalNCenergy[md.getModelId()]=0;
        }

        for(auto& v : nc){
            ATOM* atmi = v.second.atmi;
            ATOM* atmj = v.second.atmj;
            if(atmi->getModelId()!=atmj->getModelId()){
                std::stringstream ss;
                ss << "Error in native contact list (1), it should not have happen";
                throw std::runtime_error(ss.str());
            }
            totalNCenergy[atmi->getModelId()]+=v.second.e;
        }

        for(auto& tNCe : totalNCenergy){
            int mId = tNCe.first;
            int N = proteinOut.model(mId).atom().size();
            tNCe.second = epsRes*N/tNCe.second;
        }
        
        for(auto& v : nc){
            ATOM* atmi = v.second.atmi;
            ATOM* atmj = v.second.atmj;
            if(atmi->getModelId()!=atmj->getModelId()){
                std::stringstream ss;
                ss << "Error in native contact list (2), it should not have happen";
                throw std::runtime_error(ss.str());
            }
            v.second.eScaled=v.second.e*totalNCenergy[atmi->getModelId()];
        }
    
#endif

    ////////////////////////////////////////////////
    
    std::ofstream bonds(outputFileName+".bond");
    
    bonds << bondVector.size() << std::endl;
    for(auto& b : bondVector){
        bonds << std::left << std::setw(6) 
              << b.i  << " " 
              << std::left << std::setw(6) 
              << b.j  << " "
              << std::fixed << std::setprecision(6)
              << b.r0 << " " << real(200)*epsRes << std::endl;
    }
    
    std::ofstream angles(outputFileName+".angle");
    
    angles << angleVector.size() << std::endl;
    for(auto& a : angleVector){
        angles << std::left << std::setw(6)
               << a.i      << " " 
               << std::left << std::setw(6) 
               << a.j      << " "
               << std::left << std::setw(6) 
               << a.k      << " "
               << std::fixed << std::setprecision(6)
               << a.theta0 << " " << real(40)*epsRes << std::endl;
    }
    
    std::ofstream dihes(outputFileName+".dihe");
    
    dihes << diheVector.size() << std::endl;
    for(auto& d : diheVector){
        
        dihes << std::left << std::setw(6)
              << d.i  << " " 
              << std::left << std::setw(6) 
              << d.j  << " "
              << std::left << std::setw(6) 
              << d.k  << " "
              << std::left << std::setw(6) 
              << d.m  << " "
              << std::fixed << std::setprecision(6)
              << d.K1*real(0.4)*epsRes << " "
              << std::fixed << std::setprecision(6)
              << d.K2*real(0.4)*epsRes << " "
              << std::fixed << std::setprecision(6)
              << d.K3*real(0.4)*epsRes << " "
              << std::fixed << std::setprecision(6)
              << d.K4*real(0.4)*epsRes << " "
              << std::fixed << std::setprecision(6)
              << d.a1 << " "
              << std::fixed << std::setprecision(6)
              << d.a2 << " "
              << std::fixed << std::setprecision(6)
              << d.a3 << " "
              << std::fixed << std::setprecision(6)
              << d.a4 << std::endl;
    }
    
    std::ofstream Qdetails(outputFileName+".Qdetails");
    std::ofstream Qlist(outputFileName+".Qlist");
    
    Qlist<< nc.size() << std::endl;

    for(auto& v : nc){
        ATOM* atmi = v.second.atmi;
        ATOM* atmj = v.second.atmj;
        std::string details = v.second.type;
        
        /*
        Qdetails << atmi->getAtomName() << " " << atmi->getAtomSerial() << " " 
                 << atmi->getResSeq()   << " " << atmi->getChainId()    << " " 
                 << atmi->getModelId()  << " " <<
                    atmj->getAtomName() << " " << atmj->getAtomSerial() << " " 
                 << atmj->getResSeq()   << " " << atmj->getChainId()    << " " 
                 << atmj->getModelId()  << " " << details << std::endl;
        */
        
        //Qdetail is in A and kcal/mol !!!!!
        Qdetails << std::left << std::setw(3)
                 << atmi->getAtomName() << " "
                 << std::left << std::setw(3)
                 << atmi->getResSeq()   << " " 
                 << std::left << std::setw(3)
                 << atmj->getAtomName() << " "  
                 << std::left << std::setw(3)
                 << atmj->getResSeq()   << " "
                 << std::fixed << std::setprecision(6)
                 << -v.second.e         << " "
                 << v.second.type       << std::endl;

        Qlist << std::left << std::setw(6)
              << atmi->getAtomSerial() << " " 
              << std::left << std::setw(6) 
              << atmj->getAtomSerial() << " "
              << std::fixed << std::setprecision(6)
              << atomsDst(*atmi,*atmj) << " "
              << std::fixed << std::setprecision(6)
              << v.second.eScaled      << std::endl;
    }                              
                                   
    std::ofstream iRad(outputFileName+".iRad");
    
    for(uint i = 0 ; i < innerRadius.size(); i++){
        iRad << std::left << std::setw(6)
             << atomOutVector[i].getAtomSerial() << " " 
             << std::fixed << std::setprecision(6)
             << innerRadius[i] << std::endl;
    }

    std::ofstream excl(outputFileName+".excl");
    
    int nExclusions=0;
    int maxExclusion=0;
    for(auto& E : exclusionsList){
        std::sort(E.second.begin(),E.second.end());
        E.second.erase(std::unique(E.second.begin(),E.second.end()),E.second.end());
        nExclusions+=E.second.size();
        if(E.second.size()>maxExclusion){
            maxExclusion=E.second.size();
        }
    }
    
    excl << nExclusions << " " << maxExclusion << std::endl;
    for(auto& E : exclusionsList){
        excl << E.first << " ";
        
        for(auto& e : E.second){
            excl << std::left << std::setw(6) << e << " ";
        }

        excl << std::endl;
    }

    ////////////////////////////////////////////////
    
    std::cerr << "Generating SP file" << std::endl;

    std::ofstream sp(outputFileName+".sp");
    
    std::map<std::string,int> chain2int;
    int chainCount = 0;
    for(uint i = 0  ;i<atomOutVector.size();i++){
        
        if(chain2int.count(atomOutVector[i].getChainId())==0){
           chain2int[atomOutVector[i].getChainId()]=chainCount;
           chainCount++;
        }
        
        sp << std::left << std::setw(3) 
           << atomOutVector[i].getAtomName()           << " " 
           << std::left << std::setw(6) 
           << atomOutVector[i].getAtomSerial()         << " " 
           << std::left << std::setw(6) 
           << atomOutVector[i].getResSeq()             << " " 
           << std::left << std::setw(6) 
           << chain2int[atomOutVector[i].getChainId()] << " " 
           #ifdef MERGE_MODELS
           << std::left << std::setw(6) 
           << atomOutVector[i].getModelId()            << " " 
           #else
           << std::left << std::setw(6) 
           << 0                                        << " "
           #endif
           << atomOutVector[i].getAtomCoord()          << "        "
           << innerRadius[i]                           << " "
           //<< std::endl;
           << atomOutVector[i].getAtomSASA()           << std::endl;
    }

    std::cout << "epsRes: " << epsRes;     
    
    #ifdef MERGE_MODELS
        
        std::cout << " epsNC: " << epsNC << std::endl;

    #else
        
        std::cout << " epsNC: " ;

        for(auto& tNCe : totalNCenergy){
            std::cout << tNCe.second << " ";
        }
        
        std::cout << std::endl;

    #endif
}
