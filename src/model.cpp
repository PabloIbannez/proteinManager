#include "model.hpp"


namespace proteinManager {

    MODEL::MODEL(int modelID,STRUCTURE* parentStructure){
            
            this->setModelId(modelID);
            parentStructure_ = parentStructure;
        }
    
    MODEL::MODEL(int modelID) {
        MODEL(modelID,nullptr);
    }
    
    CHAIN& MODEL::chain(std::string chainID) {
        for(CHAIN& ch : chainVector) {
            if(ch.getChainId() == chainID) {
                return ch;
            }
        }
    
        std::stringstream ss;
        ss << "No chain added with the chainID \"" << chainID << "\"" << std::endl;
        throw std::runtime_error(ss.str());
    }
    
    boost::ptr_vector<CHAIN>& MODEL::chain() {
        return chainVector;
    }
    
    bool MODEL::isChain(std::string chainID){
        
        for(CHAIN& ch : chainVector){
            if(ch.getChainId() == chainID){
                return true;
            }
        }
        
        return false;
    }
    
    STRUCTURE* MODEL::getParentStructure() const {
        if(parentStructure_ == nullptr) {
            std::stringstream ss;
            ss << "ERROR (MODEL). No parent structure has been provided";
            throw std::runtime_error(ss.str());
        } else {
            return parentStructure_;
        }
    }
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_MDL_PROPERTY_IMPL_T(type,Name,name)  GET_MDL_PROPERTY_IMPL_R(type,Name,name)
    #define GET_MDL_PROPERTY_IMPL_R(type,Name,name)  type MODEL::getModel##Name() const{ \
                                                                                         try { return boost::any_cast<type>(modelProperties.at(#name));} \
                                                                                         catch (const std::out_of_range& e){ \
                                                                                             std::stringstream ss; \
                                                                                             ss << "ERROR ( " << e.what() << " ). The property \"" << #name << "\" has not been added previously." ; \
                                                                                             throw std::runtime_error(ss.str());}\
                                                                                       }
    #define GET_MDL_PROPERTY_IMPL(r, data, tuple) GET_MDL_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        MDL_PROPERTY_LOOP(GET_MDL_PROPERTY_IMPL)
    
    #undef GET_MDL_PROPERTY_IMPL_T
    #undef GET_MDL_PROPERTY_IMPL_R
    #undef GET_MDL_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define SET_MDL_PROPERTY_IMPL_T(type,Name,name)  SET_MDL_PROPERTY_IMPL_R(type,Name,name)
    #define SET_MDL_PROPERTY_IMPL_R(type,Name,name)  void MODEL::setModel##Name( type name ) { modelProperties[#name] = name;}
    #define SET_MDL_PROPERTY_IMPL(r, data, tuple) SET_MDL_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
    MDL_PROPERTY_LOOP(SET_MDL_PROPERTY_IMPL)
    
    #undef SET_MDL_PROPERTY_IMPL_T
    #undef SET_MDL_PROPERTY_IMPL_R
    #undef SET_MDL_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    void MODEL::addChain(std::string chainID) {
    
        for(CHAIN& ch : chainVector) {
            if(chainID == ch.getChainId()) {
                std::stringstream ss;
                ss << "Error (MODEL). A chain with the chain identifier \"" << chainID << "\" has been added before.";
                throw std::runtime_error(ss.str());
            }
        }
    
        chainVector.push_back(new CHAIN(chainID,this));
    
    }
    
    void MODEL::addResidue(std::string chainID,
                           std::string resName,int resSeq,std::string iCode) {
    
        for(CHAIN& ch : chainVector) {
            if(chainID == ch.getChainId()) {
                ch.addResidue(resName,resSeq,iCode);
                return;
            }
        }
    
        this->addChain(chainID);
        chainVector.back().addResidue(resName,resSeq,iCode);
    
    }
    
    void MODEL::addAtom(std::string chainID,
                        std::string resName, int resSeq, std::string iCode,
                        int serial, std::string name) {
    
        for(CHAIN& ch : chainVector) {
            if(chainID == ch.getChainId()) {
                ch.addAtom(resName,resSeq,iCode,
                           serial,name);
                return;
            }
        }
    
        this->addChain(chainID);
        chainVector.back().addAtom(resName,resSeq,iCode,
                                   serial,name);
    }
    
    std::ostream& operator<<(std::ostream& os, const MODEL& md) {
    
        DATA_FORMAT outputFormat = md.getParentStructure()->outputFormat;
    
        switch(outputFormat) {
    
        case PDB:
        case PDBQ:
        case PQR:
        case SP:
        case SPQ:
        case XYZ:
        case PDRS:
    
            for(size_t i = 0; i != md.chainVector.size(); i++) {
                if( i != md.chainVector.size()-1) {
                    os << md.chainVector[i] << std::endl;
                } else {
                    os << md.chainVector[i] ;
                }
            }
    
            break;
        }
    
        return os;
    }

}

