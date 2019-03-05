#ifndef ENM_HPP
#define ENM_HPP

#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    struct bond{
            
            std::shared_ptr<ATOM> ptr1;
            std::shared_ptr<ATOM> ptr2;
            
            double r0;
            double k;
            
    };           
    
    template<class enmModel>
    class enm;
    
    template <class T>
    std::ostream& operator<<(std::ostream& os, const enm<T>& enmOut);
        
    template<class enmModel>
    class enm{
        
            friend std::ostream& operator<< <>( std::ostream& os, const enm& enmOut );
        
        private:
            
            std::shared_ptr<enmModel> ENM;
            
            std::vector<bond> network;
        
        public:
        
            void computeENMbetweenModels(STRUCTURE& strIn){
                
                ENM = std::make_shared<enmModel>();
                
                if(!ENM->init(strIn)){
                    std::stringstream ss;
                    ss << "The current structure is not compatible with the selected elastic network model";
                    throw std::runtime_error(ss.str());
                }
                
                //Model vector
                
                std::vector<int> modelCount;
                
                for(MODEL&   mdl : strIn.model()){ 
                    std::cout << "Adding model: " << mdl.getModelId() << std::endl;
                    modelCount.push_back(mdl.getModelId());
                }
                
                for(int i = 0;i<modelCount.size();i++){
                    
                    std::cout << "Current model i:" << modelCount[i] << std::endl;
                    
                    //Buffer atom vector
                    std::vector<std::shared_ptr<proteinManager::ATOM>> atomVectorMdl_i;
                    std::vector<std::shared_ptr<proteinManager::ATOM>> atomVectorMdl_j;
                    
                    for(CHAIN&   ch  : strIn.model(modelCount[i]).chain()  ){
                    for(RESIDUE& res : ch.residue() ){
                    for(ATOM&    atm : res.atom()   ){
                        atomVectorMdl_i.push_back(std::make_shared<proteinManager::ATOM>(atm));
                    }}}
                    
                    
                    for(int j = i+1;j<modelCount.size();j++){
                        
                        std::cout << " ---- Current model j:" << modelCount[j] << std::endl;
                        
                        for(CHAIN&   ch  : strIn.model(modelCount[j]).chain()  ){
                        for(RESIDUE& res : ch.residue() ){
                        for(ATOM&    atm : res.atom()   ){
                            atomVectorMdl_j.push_back(std::make_shared<proteinManager::ATOM>(atm));
                        }}}
                        
                        for(int k=0; k < atomVectorMdl_i.size(); k++){
                        for(int t=0; t < atomVectorMdl_j.size(); t++){
                            
                            bond bd = ENM->computeBond(atomVectorMdl_i[k],atomVectorMdl_j[t]);
                            
                            if(bd.k > 0){
                                network.push_back(bd);
                            }
                        }}
                        
                    }
                }
                
            }
            
            void computeENM(STRUCTURE& strIn){
                
                ENM = std::make_shared<enmModel>();
                
                if(!ENM->init(strIn)){
                    std::stringstream ss;
                    ss << "The current structure is not compatible with the selected elastic network model";
                    throw std::runtime_error(ss.str());
                }
                
                //Buffer atom vector
                std::vector<std::shared_ptr<proteinManager::ATOM>> atomVector;
                
                
                for(MODEL&   mdl : strIn.model()){
                for(CHAIN&   ch  : mdl.chain()  ){
                for(RESIDUE& res : ch.residue() ){
                for(ATOM&    atm : res.atom()   ){
                    atomVector.push_back(std::make_shared<proteinManager::ATOM>(atm));
                }}}}
                
                for(int i=0; i < atomVector.size(); i++){
                for(int j=i+1; j < atomVector.size(); j++){
                    
                    bond bd = ENM->computeBond(atomVector[i],atomVector[j]);
                    
                    if(bd.k > 0){
                        network.push_back(bd);
                    }
                }}
                
            }
            
            void computeENM(MODEL& mdlIn){
                
                ENM = std::make_shared<enmModel>();
                
                if(!ENM->init(mdlIn)){
                    std::stringstream ss;
                    ss << "The current structure is not compatible with the selected elastic network model";
                    throw std::runtime_error(ss.str());
                }
                
                //Buffer atom vector
                std::vector<std::shared_ptr<proteinManager::ATOM>> atomVector;
                
                for(CHAIN&   ch  : mdlIn.chain()  ){
                for(RESIDUE& res : ch.residue()  ){
                for(ATOM&    atm : res.atom()    ){
                    atomVector.push_back(std::make_shared<proteinManager::ATOM>(atm));
                }}}
                
                for(int i=0; i < atomVector.size(); i++){
                for(int j=i+1; j < atomVector.size(); j++){
                    
                    bond bd = ENM->computeBond(atomVector[i],atomVector[j]);
                    
                    if(bd.k > 0){
                        network.push_back(bd);
                    }
                }}
                
            }
            
            void computeENM(CHAIN& chIn){
                
                ENM = std::make_shared<enmModel>();
                
                if(!ENM->init(chIn)){
                    std::stringstream ss;
                    ss << "The current structure is not compatible with the selected elastic network model";
                    throw std::runtime_error(ss.str());
                }
                
                //Buffer atom vector
                std::vector<std::shared_ptr<proteinManager::ATOM>> atomVector;
                
                for(RESIDUE& res : chIn.residue()  ){
                for(ATOM&    atm : res.atom()    ){
                    atomVector.push_back(std::make_shared<proteinManager::ATOM>(atm));
                }}
                
                for(int i=0; i < atomVector.size(); i++){
                for(int j=i+1; j < atomVector.size(); j++){
                    
                    bond bd = ENM->computeBond(atomVector[i],atomVector[j]);
                    
                    if(bd.k > 0){
                        network.push_back(bd);
                    }
                }}
                
            }
            
            void computeENM(RESIDUE& resIn){
                
                ENM = std::make_shared<enmModel>();
                
                if(!ENM->init(resIn)){
                    std::stringstream ss;
                    ss << "The current structure is not compatible with the selected elastic network model";
                    throw std::runtime_error(ss.str());
                }
                
                //Buffer atom vector
                std::vector<std::shared_ptr<proteinManager::ATOM>> atomVector;
                
                for(ATOM&    atm : resIn.atom()    ){
                    atomVector.push_back(std::make_shared<proteinManager::ATOM>(atm));
                }
                
                for(int i=0; i < atomVector.size(); i++){
                for(int j=i+1; j < atomVector.size(); j++){
                    
                    bond bd = ENM->computeBond(atomVector[i],atomVector[j]);
                    
                    if(bd.k > 0){
                        network.push_back(bd);
                    }
                }}
                
            }
    };
    
    template <class T>
    std::ostream& operator<<(std::ostream& os, const enm<T>& enmOut){
        
        if(enmOut.network.size() <= 0){
            std::stringstream ss;
            ss << "No network has been generated";
            throw std::runtime_error(ss.str());
        }
        
        os << enmOut.ENM->info;
        os << enmOut.network.size() << std::endl;
        
        for(const bond& bd : enmOut.network){
            
            os << std::setw(10) << bd.ptr1->getAtomSerial() <<
                  std::setw(10) << bd.ptr2->getAtomSerial() <<
                  std::setprecision(6)                     <<
                  std::setw(12) << bd.k                    <<
                  std::setw(12) << bd.r0                   << std::endl;
            
            /*
            ///////////////////////////////////////////
            
            os << std::left << std::fixed            <<
            std::setw(6) << "ATOM"                   <<
            std::right                               <<
            std::setw(5) << bd.ptr1->getAtomSerial() <<
            " "                                   ;
        
            if(bd.ptr1->getAtomName().size() < 4) {
                os << std::left << std::fixed <<" "   <<
                std::setw(3) << bd.ptr1->getAtomName() ;
            } else {
                os << std::left << std::fixed         <<
                std::setw(4) << bd.ptr1->getAtomName() ;
            }
        
            os << std::left << std::fixed             <<
            std::setw(1) << bd.ptr1->getAtomAltLoc()  <<
            std::setw(3) << bd.ptr1->getResName()     <<
            " "                                       <<
            std::setw(1) << bd.ptr1->getChainId()     <<
            std::right                                <<
            std::setw(4) << bd.ptr1->getResSeq()      <<
            std::setw(1) << bd.ptr1->getResInsCode()  <<
            "   "                                     <<
            std::setprecision(3)                      <<
            std::setw(8) << bd.ptr1->getAtomCoord().x <<
            std::setw(8) << bd.ptr1->getAtomCoord().y <<
            std::setw(8) << bd.ptr1->getAtomCoord().z <<
            "   "                                     ;
            
            ///////////////////////////////////////////
            
            os << std::left << std::fixed            <<
            std::setw(6) << "ATOM"                   <<
            std::right                               <<
            std::setw(5) << bd.ptr2->getAtomSerial() <<
            " "                                   ;
            
            if(bd.ptr2->getAtomName().size() < 4) {
                os << std::left << std::fixed <<" "   <<
                std::setw(3) << bd.ptr2->getAtomName() ;
            } else {
                os << std::left << std::fixed         <<
                std::setw(4) << bd.ptr2->getAtomName() ;
            }
            
            os << std::left << std::fixed             <<
            std::setw(1) << bd.ptr2->getAtomAltLoc()  <<
            std::setw(3) << bd.ptr2->getResName()     <<
            " "                                       <<
            std::setw(1) << bd.ptr2->getChainId()     <<
            std::right                                <<
            std::setw(4) << bd.ptr2->getResSeq()      <<
            std::setw(1) << bd.ptr2->getResInsCode()  <<
            "   "                                     <<
            std::setprecision(3)                      <<
            std::setw(8) << bd.ptr2->getAtomCoord().x <<
            std::setw(8) << bd.ptr2->getAtomCoord().y <<
            std::setw(8) << bd.ptr2->getAtomCoord().z <<
            "   "                                     ;
            
            ////////////////////////////////////////////
            
            os << std::setprecision(6)         <<
                  std::setw(12) << bd.r0       <<
                  std::setw(12) << bd.k        << std::endl;
            */
            
        }
        
        return os;
    }
}

#include "enm_models.hpp"

#endif
