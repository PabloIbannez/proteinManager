#include <proteinManager/proteinManager.hpp>

namespace proteinManager{
    
    namespace misc{

        namespace misc_ns{
            
            void toHIS(RESIDUE& res){
                if(res.getResName().substr(0,3) == "HID" or 
                   res.getResName().substr(0,3) == "HIE" or 
                   res.getResName().substr(0,3) == "HIP" ){
                    res.setResName("HIS");
                }
            }
            
            void toCYS(RESIDUE& res){
                if(res.getResName().substr(0,3) == "CYX"){
                    res.setResName("CYS");
                }
            }
        }

        void fixHIS(STRUCTURE& structure){
            for(MODEL&   mdl : structure.model()){
            for(CHAIN&   ch  : mdl.chain()       ){
            for(RESIDUE& res : ch.residue()      ){
                misc_ns::toHIS(res);
            }}}
        }
        
        void fixCYS(STRUCTURE& structure){
            for(MODEL&   mdl : structure.model()){
            for(CHAIN&   ch  : mdl.chain()       ){
            for(RESIDUE& res : ch.residue()      ){
                misc_ns::toCYS(res);
            }}}
        }

        real beadCharge(std::string const & beadName){
        
            if(beadName.substr(0,3) == "LYS" or 
               beadName.substr(0,3) == "ARG"){
                return 1.0;
            }
            
            if(beadName.substr(0,3) == "ASP" or 
               beadName.substr(0,3) == "GLU"){
                return -1.0;
            }
            
            if(beadName.substr(0,3) == "HIS" or 
               beadName.substr(0,3) == "HID" or 
               beadName.substr(0,3) == "HIE" or 
               beadName.substr(0,3) == "HIP" ){
                return 0.5;
            }
        
            return 0.0;
        }
        
        void setAtomsChargeEqualToResCharge(STRUCTURE& structure){
            for(MODEL&   mdl : structure.model()){
            for(CHAIN&   ch  : mdl.chain()       ){
            for(RESIDUE& res : ch.residue()      ){
            for(ATOM&    atm : res.atom()         ){
                atm.setAtomCharge(beadCharge(res.getResName()));
            }}}}
        }
        
        void setAtomsNameEqualToResName(STRUCTURE& structure){
            for(MODEL&   mdl : structure.model()){
            for(CHAIN&   ch  : mdl.chain()       ){
            for(RESIDUE& res : ch.residue()      ){
            for(ATOM&    atm : res.atom()         ){
                atm.setAtomName(atm.getResName());
            }}}}
        }
    }
}
