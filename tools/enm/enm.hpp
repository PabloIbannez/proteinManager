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
            
            std::vector<bond> network;
        
        public:
            
            void computeENM(STRUCTURE& strIn){
                
                enmModel ENM;
                
                if(!ENM.init(strIn)){
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
                    
                    bond bd = ENM.computeBond(atomVector[i],atomVector[j]);
                    
                    //std::cout << bond.ptr1.getAtomName() << std::endl;
                    //std::cin.get();
                    
                    if(bd.k > 0){
                        network.push_back(bd);
                    }
                }}
                
            }
    };
    
    template <class T>
    std::ostream& operator<<(std::ostream& os, const enm<T>& enmOut){
        
        os << enmOut.network.size() << std::endl;
        
        for(const bond& bd : enmOut.network){
            
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
            "   "                                     ;
            
            ////////////////////////////////////////////
            
            os << std::setprecision(6)         <<
                  std::setw(12) << bd.r0       <<
                  std::setw(12) << bd.k        << std::endl;
            
        }
        
        return os;
    }
}

#include "enm_models.hpp"

#endif
