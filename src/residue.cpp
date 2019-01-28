#include "residue.hpp"

namespace proteinManager {

    RESIDUE::RESIDUE(std::string resName,int resSeq,std::string iCode,CHAIN* parentChain){
        
        this->setResName(resName);
        this->setResSeq(resSeq);
        this->setResInsCode(iCode);
        
        parentChain_ = parentChain;
    }
    
    RESIDUE::RESIDUE(std::string resName, int resSeq, std::string iCode) {
        RESIDUE(resName,resSeq,iCode,nullptr);
    }
    
    ATOM& RESIDUE::atom(std::string atomName) {
        if(atomMap.count(atomName) < 1) {
            std::stringstream ss;
            ss << "ERROR. No atom called \"" << atomName << "\" has been added to the residue \"" << this->getResName()
                                                         << "\" (" << this->getResSeq() << ")." << std::endl;
            throw std::runtime_error(ss.str());
        } else {
            return atomVector[atomMap[atomName]];
        }
    }   
    
    boost::ptr_vector<ATOM>& RESIDUE::atom() {
        return atomVector;
    }
    
    bool RESIDUE::isAtom(std::string atomName) {
    
        if(atomMap.count(atomName) < 1) {
            return false;
        } else {
            return true;
        }
    }
    
    CHAIN* RESIDUE::getParentChain() const {
        if(parentChain_ == nullptr) {
            std::stringstream ss;
            ss << "ERROR (RESIDUE). No parent chain has been provided" ;
            throw std::runtime_error(ss.str());
        } else {
            return parentChain_;
        }
    }
    
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_MDL_PROPERTY_IMPL_T(type,Name,name)  GET_MDL_PROPERTY_IMPL_R(type,Name,name)
    #define GET_MDL_PROPERTY_IMPL_R(type,Name,name)  type RESIDUE::getModel##Name() const{ return this->getParentChain()->getModel##Name();}
    #define GET_MDL_PROPERTY_IMPL(r, data, tuple) GET_MDL_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
    MDL_PROPERTY_LOOP(GET_MDL_PROPERTY_IMPL)
    
    #undef GET_MDL_PROPERTY_IMPL_T
    #undef GET_MDL_PROPERTY_IMPL_R
    #undef GET_MDL_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_CHAIN_PROPERTY_IMPL_T(type,Name,name)  GET_CHAIN_PROPERTY_IMPL_R(type,Name,name)
    #define GET_CHAIN_PROPERTY_IMPL_R(type,Name,name)  type RESIDUE::getChain##Name() const{ return this->getParentChain()->getChain##Name();}
    #define GET_CHAIN_PROPERTY_IMPL(r, data, tuple) GET_CHAIN_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
    CHAIN_PROPERTY_LOOP(GET_CHAIN_PROPERTY_IMPL)
    
    #undef GET_CHAIN_PROPERTY_IMPL_T
    #undef GET_CHAIN_PROPERTY_IMPL_R
    #undef GET_CHAIN_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_RES_PROPERTY_IMPL_T(type,Name,name)  GET_RES_PROPERTY_IMPL_R(type,Name,name)
    #define GET_RES_PROPERTY_IMPL_R(type,Name,name)  type RESIDUE::getRes##Name() const{ \
                                                                                         try { return std::any_cast<type>(residueProperties.at(#name));} \
                                                                                         catch (const std::out_of_range& e){ \
                                                                                             std::stringstream ss; \
                                                                                             ss << "ERROR ( " << e.what() << " ). The property \"" << #name << "\" has not been added previously." ; \
                                                                                             throw std::runtime_error(ss.str());}\
                                                                                       }
    #define GET_RES_PROPERTY_IMPL(r, data, tuple) GET_RES_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        RES_PROPERTY_LOOP(GET_RES_PROPERTY_IMPL)
    
    #undef GET_RES_PROPERTY_IMPL_T
    #undef GET_RES_PROPERTY_IMPL_R
    #undef GET_RES_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define SET_RES_PROPERTY_IMPL_T(type,Name,name)  SET_RES_PROPERTY_IMPL_R(type,Name,name)
    #define SET_RES_PROPERTY_IMPL_R(type,Name,name)  void RESIDUE::setRes##Name( type name ) { residueProperties[#name] = name;}
    #define SET_RES_PROPERTY_IMPL(r, data, tuple) SET_RES_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        RES_PROPERTY_LOOP(SET_RES_PROPERTY_IMPL)
    
    #undef SET_RES_PROPERTY_IMPL_T
    #undef SET_RES_PROPERTY_IMPL_R
    #undef SET_RES_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    
    void RESIDUE::addAtom(int serial,std::string name) {
    
        atomVector.push_back(new ATOM(serial,name,this));
        atomMap[name]=atomVector.size()-1;
    }
    
    std::ostream& operator<<(std::ostream& os, const RESIDUE& res) {
    
        DATA_FORMAT outputFormat = res.getParentChain()->getParentModel()->getParentStructure()->outputFormat;
    
        size_t i;
    
        switch(outputFormat) {
    
        case PDB:
        case PDBQ:
        case PQR:
        case PDRS:
        case SP:
        case SPQ:
        case XYZ:
    
            for(i = 0; i < res.atomVector.size(); i++) {
                if( i != res.atomVector.size()-1) {
                    os << res.atomVector[i] << std::endl;
                } else {
                    os << res.atomVector[i];
                }
            }
    
            break;
        }
    
        return os;
    }
}
