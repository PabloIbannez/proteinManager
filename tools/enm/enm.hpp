#ifndef ENM_HPP
#define ENM_HPP

#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    struct bondInfo{
            
            ATOM& bondedParticle1;
            ATOM& bondedParticle2;
            double k;
    };
        
    template<class enmModel>
    class enm{
        
        private:
            
            std::vector<bondInfo> network;
        
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
                for(int j=0; j < atomVector.size(); j++){
                    
                    if (i != j){
                        double currentK = ENM.springConstant(*atomVector[i],*atomVector[j]);
                        
                        if(currentK > 0){
                            network.push_back({*atomVector[i],*atomVector[j],currentK});
                        }
                    }
                }}
                
            }
            
            
        
    };
    
    namespace enm_models{
        
        struct caOrellana{
            
            double rcut;
            
            int M = 3;
            
            double Cseq = 60;
            double Ccart = 6;
            
            int n_seq = 2;
            int n_cart = 6;
            
            const char* info = "units are ...";
            
            bool init(STRUCTURE& str){ 
                
                rcut = 17;
                return true;
            }
            bool init(MODEL&     str){ return true;}
            bool init(CHAIN&     str){ return true;}
            bool init(RESIDUE&   str){ return true;}
            
            double springConstant(ATOM& atm1, ATOM& atm2){
                
                if(atm1.getAtomName() != "CA" or atm2.getAtomName() != "CA"){ return 0;}
                
                int S12 = abs(atm1.getParentResidue()->getResSeq()-atm2.getParentResidue()->getResSeq());
                
                if(S12 <= M){
                    
                    return double(Cseq)/(pow(S12,n_seq));
                    
                } else {
                    
                    real3 r12 = atm1.getAtomCoord()-atm2.getAtomCoord();
                    double r = sqrt(dot(r12,r12));
                    
                    if(r <= rcut){
                        //check pow(Ccart/r,n_cart)
                        std::cout << double(Ccart)/(pow(r,n_cart)) << std::endl;
                        return double(Ccart)/(pow(r,n_cart));
                    } else {
                        return 0;
                    }
                }
            }
        };
        
    }

}

#endif
