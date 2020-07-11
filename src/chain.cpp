/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/

#include "chain.hpp"

namespace proteinManager {

    CHAIN::CHAIN(std::string chainID,MODEL* parentModel)
    {
            this->setChainId(chainID);
            parentModel_ = parentModel;

    }
    
    CHAIN::CHAIN(std::string chainID) {
        CHAIN(chainID,nullptr);
    }
    
    RESIDUE& CHAIN::residue(int resSeq) {
        for(RESIDUE& res : residueVector) {
            if(res.getResSeq() == resSeq) {
                return res;
            }
        }
    
        std::stringstream ss;
        ss << "No residue added with the resSeq \"" << resSeq << "\"";
        throw std::runtime_error(ss.str());
    }
    
    boost::ptr_vector<RESIDUE>& CHAIN::residue() {
        return residueVector;
    }
    
    boost::ptr_vector<ATOM>& CHAIN::atom(){
        
        atomVector.clear();

        for(RESIDUE& res : residueVector){
            atomVector.insert(atomVector.end(),res.atom().begin(),res.atom().end());
        }

        return atomVector;
    }
    
    bool CHAIN::isRes(int resSeq){
        
        for(RESIDUE& res : residueVector){
            if(res.getResSeq() == resSeq){
                return true;
            }
        }
        
        return false;
    }
    
    STRUCTURE& CHAIN::getParentStructure() const {
        return this->getParentModel().getParentStructure();    
    }
    
    MODEL& CHAIN::getParentModel() const {
        if(parentModel_ == nullptr) {
            std::stringstream ss;
            ss << "ERROR (CHAIN). No parent model has been provided" ;
            throw std::runtime_error(ss.str());
        } else {
            return *parentModel_;
        }
    }
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_MDL_PROPERTY_IMPL_T(type,Name,name)  GET_MDL_PROPERTY_IMPL_R(type,Name,name)
    #define GET_MDL_PROPERTY_IMPL_R(type,Name,name)  type CHAIN::getModel##Name() const{ return this->getParentModel().getModel##Name();}
    #define GET_MDL_PROPERTY_IMPL(r, data, tuple) GET_MDL_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
    MDL_PROPERTY_LOOP(GET_MDL_PROPERTY_IMPL)
    
    #undef GET_MDL_PROPERTY_IMPL_T
    #undef GET_MDL_PROPERTY_IMPL_R
    #undef GET_MDL_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_CHAIN_PROPERTY_IMPL_T(type,Name,name)  GET_CHAIN_PROPERTY_IMPL_R(type,Name,name)
    #define GET_CHAIN_PROPERTY_IMPL_R(type,Name,name)  type CHAIN::getChain##Name() const{ \
                                                                                         try { return boost::any_cast<type>(chainProperties.at(#name));} \
                                                                                         catch (const std::out_of_range& e){ \
                                                                                             std::stringstream ss; \
                                                                                             ss << "ERROR ( " << e.what() << " ). The property \"" << #name << "\" has not been added previously." ; \
                                                                                             throw std::runtime_error(ss.str());}\
                                                                                       }
    #define GET_CHAIN_PROPERTY_IMPL(r, data, tuple) GET_CHAIN_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
    CHAIN_PROPERTY_LOOP(GET_CHAIN_PROPERTY_IMPL)
    
    #undef GET_CHAIN_PROPERTY_IMPL_T
    #undef GET_CHAIN_PROPERTY_IMPL_R
    #undef GET_CHAIN_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define SET_CHAIN_PROPERTY_IMPL_T(type,Name,name)  SET_CHAIN_PROPERTY_IMPL_R(type,Name,name)
    #define SET_CHAIN_PROPERTY_IMPL_R(type,Name,name)  void CHAIN::setChain##Name( type name ) { chainProperties[#name] = name;}
    #define SET_CHAIN_PROPERTY_IMPL(r, data, tuple) SET_CHAIN_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
    CHAIN_PROPERTY_LOOP(SET_CHAIN_PROPERTY_IMPL)
    
    #undef SET_CHAIN_PROPERTY_IMPL_T
    #undef SET_CHAIN_PROPERTY_IMPL_R
    #undef SET_CHAIN_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    void CHAIN::addResidue(std::string resName,int resSeq,std::string iCode) {
    
        for(RESIDUE& res : residueVector) {
            if(res.getResName()    == resName and
               res.getResSeq()     == resSeq  and
               res.getResInsCode() == iCode 
               ) {
                std::stringstream ss;
                ss << "Error (CHAIN). A residue with the chain identifier \"" 
                   << resName << "\" "  
                   << resSeq  << "\" " 
                   << iCode   << "\" " 
                   << "has been added before.";
                throw std::runtime_error(ss.str());
            }
        }
    
        this->residueVector.push_back(new RESIDUE(resName,resSeq,iCode,this));
    }
    
    void CHAIN::addAtom(std::string resName, int resSeq, std::string iCode,
                        int serial, std::string name) {
    
    
        for(RESIDUE& res: residueVector) {
            if(resName == res.getResName() and
               resSeq  == res.getResSeq()  and
               iCode   == res.getResInsCode()) {
                res.addAtom(serial,name);
                return;
            }
        }
    
        this->addResidue(resName,resSeq,iCode);
        residueVector.back().addAtom(serial,name);
    
    }
    
    std::ostream& operator<<(std::ostream& os, const CHAIN& chain) {
    
        DATA_FORMAT outputFormat = chain.getParentModel().getParentStructure().outputFormat;
    
        if(chain.residueVector.size() ==0 ) {
            return os;
        }
    
        size_t i;
    
            
        switch(outputFormat) {
    
        case PDB:
        case PDBQ:
        case PQR:
        case PDRS:
            for(i = 0; i < chain.residueVector.size(); i++) {
                os << chain.residueVector[i] << std::endl;
            }
            
            {
                int serial = chain.residueVector.back().atomVector.back().getAtomSerial()+1;
                
                os << std::left << std::fixed <<
                   std::setw(6) << "TER"   <<
                   std::right              <<
                   std::setw(5) << ((serial>99999)?99999:serial)           <<
                   "      "                <<
                   std::setw(3) << chain.residueVector[i-1].getResName()   <<
                   " "                     <<
                   std::setw(1) << chain.getChainId()  <<
                   std::setw(4) << chain.residueVector[i-1].getResSeq()    <<
                   std::setw(1) << chain.residueVector[i-1].getResInsCode()<<
                   "                                                     " ;
            }

            break;
        case SP:
        case SPQ:
        case XYZ:
    
            for(i = 0; i < chain.residueVector.size()-1 ; i++) {
                os << chain.residueVector[i] << std::endl;
            }
    
            os << chain.residueVector[chain.residueVector.size()-1];
    
            break;
        }
    
        return os;
    }
}

