#include "atom.hpp"

namespace proteinManager {
    
    /*
    #define SFINAE_DEFINE_HAS_MEMBER(name, X)                     \
    template <typename T>  \
    class has_##name  \
    {  \
            using one = char;  \
            using two = long;  \
        \
            template <class C> static constexpr inline one test( decltype(&C::X) ) ;    \
            template <class C> static constexpr inline two test(...);  \
        \
        public:  \
            static constexpr bool value = std::is_same<decltype(test<T>(0)),one>::value; \
    };  \
    
    SFINAE_DEFINE_HAS_MEMBER(mult, operator*)
    */
    
    template<class T, bool active>
    struct hasMultDelegator;
    
    template<class T>
    struct hasMultDelegator<T, true>{
        static T call(T prop,real factor, std::string propName){    
            return prop*factor;
        }
    };
    
    template<class T>
        struct hasMultDelegator<T, false>{
            static T call(T prop,real factor, std::string propName){    
            std::stringstream ss;
            ss << "The proterty " << propName << " con not been scaled.";
            throw std::runtime_error(ss.str());
        }
    };
    
    template<class T>
    T applyScale(T prop,real factor,std::string propName){
        //return hasMultDelegator<T, has_mult<T>::value>::call(prop,factor,propName);
        return hasMultDelegator<T, std::is_same<T,real>::value or std::is_same<T,real3>::value>::call(prop,factor,propName);
    }

    ATOM::ATOM(int serial,std::string name,RESIDUE* parentResidue){

        this->setAtomSerial(serial);
        this->setAtomName(name);
        
        parentResidue_ = parentResidue;
    }
    
    
    ATOM::ATOM(int serial,std::string name) {
        ATOM(serial,name,nullptr);
    }
    
    RESIDUE* ATOM::getParentResidue() const {
        if(parentResidue_ == nullptr) {
            std::stringstream ss;
            ss << "ERROR (RESIDUE). No parent residue has been provided" << std::endl;
            throw std::runtime_error(ss.str());
        } else {
            return parentResidue_;
        }
    }
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_MDL_PROPERTY_IMPL_T(type,Name,name)  GET_MDL_PROPERTY_IMPL_R(type,Name,name)
    #define GET_MDL_PROPERTY_IMPL_R(type,Name,name)  type ATOM::getModel##Name() const{ return this->getParentResidue()->getModel##Name();}
    #define GET_MDL_PROPERTY_IMPL(r, data, tuple) GET_MDL_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        MDL_PROPERTY_LOOP(GET_MDL_PROPERTY_IMPL)
    
    #undef GET_MDL_PROPERTY_IMPL_T
    #undef GET_MDL_PROPERTY_IMPL_R
    #undef GET_MDL_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_CHAIN_PROPERTY_IMPL_T(type,Name,name)  GET_CHAIN_PROPERTY_IMPL_R(type,Name,name)
    #define GET_CHAIN_PROPERTY_IMPL_R(type,Name,name)  type ATOM::getChain##Name() const{ return this->getParentResidue()->getChain##Name();}
    #define GET_CHAIN_PROPERTY_IMPL(r, data, tuple) GET_CHAIN_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        CHAIN_PROPERTY_LOOP(GET_CHAIN_PROPERTY_IMPL)
    
    #undef GET_CHAIN_PROPERTY_IMPL_T
    #undef GET_CHAIN_PROPERTY_IMPL_R
    #undef GET_CHAIN_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_RES_PROPERTY_IMPL_T(type,Name,name)  GET_RES_PROPERTY_IMPL_R(type,Name,name)
    #define GET_RES_PROPERTY_IMPL_R(type,Name,name)  type ATOM::getRes##Name() const{ \
                                                                                     return this->getParentResidue()->getRes##Name(); \
                                                                                    }
    #define GET_RES_PROPERTY_IMPL(r, data, tuple) GET_RES_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        RES_PROPERTY_LOOP(GET_RES_PROPERTY_IMPL)
    
    #undef GET_RES_PROPERTY_IMPL_T
    #undef GET_RES_PROPERTY_IMPL_R
    #undef GET_RES_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define GET_ATOM_PROPERTY_IMPL_T(type,Name,name)  GET_ATOM_PROPERTY_IMPL_R(type,Name,name)
    #define GET_ATOM_PROPERTY_IMPL_R(type,Name,name)  type ATOM::getAtom##Name() const{ \
                                                                                        try { return boost::any_cast<type>(atomProperties.at(#name));} \
                                                                                        catch (const std::out_of_range& e){ \
                                                                                            std::stringstream ss; \
                                                                                            ss << "ERROR ( " << e.what() << " ). The property \"" << #name << "\" has not been added previously." ; \
                                                                                            throw std::runtime_error(ss.str());}\
                                                                                      }
    #define GET_ATOM_PROPERTY_IMPL(r, data, tuple) GET_ATOM_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        ATOM_PROPERTY_LOOP(GET_ATOM_PROPERTY_IMPL)
    
    #undef GET_ATOM_PROPERTY_IMPL_T
    #undef GET_ATOM_PROPERTY_IMPL_R
    #undef GET_ATOM_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define SET_ATOM_PROPERTY_IMPL_T(type,Name,name)  SET_ATOM_PROPERTY_IMPL_R(type,Name,name)
    #define SET_ATOM_PROPERTY_IMPL_R(type,Name,name)  void ATOM::setAtom##Name( type name ) { atomProperties[#name] = name;}
    #define SET_ATOM_PROPERTY_IMPL(r, data, tuple) SET_ATOM_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        ATOM_PROPERTY_LOOP(SET_ATOM_PROPERTY_IMPL)
    
    #undef SET_ATOM_PROPERTY_IMPL_T
    #undef SET_ATOM_PROPERTY_IMPL_R
    #undef SET_ATOM_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    #define SCALE_ATOM_PROPERTY_IMPL_T(type,Name,name)  SCALE_ATOM_PROPERTY_IMPL_R(type,Name,name)
    #define SCALE_ATOM_PROPERTY_IMPL_R(type,Name,name)  void ATOM::scaleAtom##Name( real scaleFactor ) {atomProperties[#name] = applyScale<type>(boost::any_cast<type>(atomProperties[#name]),scaleFactor,#Name);}
    #define SCALE_ATOM_PROPERTY_IMPL(r, data, tuple) SCALE_ATOM_PROPERTY_IMPL_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple),PROPNAME(tuple))
    
        ATOM_PROPERTY_LOOP(SCALE_ATOM_PROPERTY_IMPL)
    
    #undef SCALE_ATOM_PROPERTY_IMPL_T
    #undef SCALE_ATOM_PROPERTY_IMPL_R
    #undef SCALE_ATOM_PROPERTY_IMPL
    
    ////////////////////////////////////////////////////////////////////
    
    std::ostream& operator<<(std::ostream& os, const ATOM& a) {
        
        DATA_FORMAT outputFormat = a.getParentResidue()->getParentChain()->getParentModel()->getParentStructure()->outputFormat;
    
        switch(outputFormat) {
        
        case PDB:
        
            os << std::left << std::fixed            <<
            std::setw(6) << "ATOM"                   <<
            std::right                               <<
            std::setw(5) << a.getAtomSerial()        <<
            " "                                   ;
        
            if(a.getAtomName().size() < 4) {
                os << std::left << std::fixed <<" "  <<
                std::setw(3) << a.getAtomName()   ;
            } else {
                os << std::left << std::fixed        <<
                std::setw(4) << a.getAtomName()   ;
            }
        
            os << std::left << std::fixed            <<
            std::setw(1) << a.getAtomAltLoc()        <<
            std::setw(3) << a.getResName()           <<
            " "                                      <<
            std::setw(1) << a.getChainId()           <<
            std::right                               <<
            std::setw(4) << a.getResSeq()            <<
            std::setw(1) << a.getResInsCode()        <<
            "   "                                    <<
            std::setprecision(3)                     <<
            std::setw(8) << a.getAtomCoord().x       <<
            std::setw(8) << a.getAtomCoord().y       <<
            std::setw(8) << a.getAtomCoord().z       <<
            std::setprecision(2)                     <<
            std::setw(6) << a.getAtomOccupancy()     <<
            std::setw(6) << a.getAtomTempFactor()    <<
            "          "                             <<
            std::setw(2) << a.getAtomElement()    ;
        
            if(int(a.getAtomCharge())==0) {
                os << std::left << std::fixed        <<
                std::setw(2) << "";
            } else {
                os << std::left << std::fixed        <<
                std::setw(2) << int(a.getAtomCharge());
            }
            break;
        
        case PDBQ:
        
            os << std::left << std::fixed            <<
            std::setw(6) << "ATOM"                   <<
            std::right                               <<
            std::setw(5) << a.getAtomSerial()        <<
            " "                                   ;
        
            if(a.getAtomName().size() < 4) {
                os << std::left << std::fixed <<" "  <<
                std::setw(3) << a.getAtomName()   ;
            } else {
                os << std::left << std::fixed        <<
                std::setw(4) << a.getAtomName()   ;
            }
        
            os << std::left << std::fixed            <<
            std::setw(1) << a.getAtomAltLoc()        <<
            std::setw(3) << a.getResName()           <<
            " "                                      <<
            std::setw(1) << a.getChainId()           <<
            std::right                               <<
            std::setw(4) << a.getResSeq()            <<
            std::setw(1) << a.getResInsCode()        <<
            "   "                                    <<
            std::setprecision(3)                     <<
            std::setw(8) << a.getAtomCoord().x       <<
            std::setw(8) << a.getAtomCoord().y       <<
            std::setw(8) << a.getAtomCoord().z       <<
            std::setprecision(2)                     <<
            std::setw(6) << a.getAtomOccupancy()     <<
            std::setw(6) << a.getAtomTempFactor()    <<
            "    "                                   <<
            std::setprecision(3)                     <<
            std::setw(6) << a.getAtomCharge()     ;
            break;
        
        case PQR:
        
            os << std::left << std::fixed            <<
            std::setw(6) << "ATOM"                   <<
            std::right                               <<
            std::setw(5) << a.getAtomSerial()        <<
            " "                                   ;
        
            if(a.getAtomName().size() < 4) {
                os << std::left << std::fixed <<" "  <<
                std::setw(3) << a.getAtomName()   ;
            } else {
                os << std::left << std::fixed <<
                std::setw(4) << a.getAtomName()   ;
            }
        
            os << std::left << std::fixed            <<
            std::setw(1) << a.getAtomAltLoc()        <<
            std::setw(3) << a.getResName()           <<
            " "                                      <<
            std::setw(1) << a.getChainId()           <<
            std::right                               <<
            std::setw(4) << a.getResSeq()            <<
            std::setw(1) << a.getResInsCode()        <<
            "   "                                    <<
            std::setprecision(3)                     <<
            std::setw(8) << a.getAtomCoord().x       <<
            std::setw(8) << a.getAtomCoord().y       <<
            std::setw(8) << a.getAtomCoord().z       <<
            std::setprecision(4)                     <<
            std::setw(8) << a.getAtomCharge()        <<
            std::setw(8) << a.getAtomRadius()     ;
            break;
        case PDRS:
            
            os << std::left << std::fixed            <<
            std::setw(6) << "ATOM"                   <<
            std::right                               <<
            std::setw(5) << a.getAtomSerial()        <<
            " "                                   ;
        
            if(a.getAtomName().size() < 4) {
        	os << std::left << std::fixed <<" "      <<
        	std::setw(3) << a.getAtomName()   ;
            } else {
        	os << std::left << std::fixed            <<
        	std::setw(4) << a.getAtomName()   ;
            }
        
            os << std::left << std::fixed            <<
            std::setw(1) << a.getAtomAltLoc()        <<
            std::setw(3) << a.getResName()           <<
            " "                                      <<
            std::setw(1) << a.getChainId()           <<
            std::right                               <<
            std::setw(4) << a.getResSeq()            <<
            std::setw(1) << a.getResInsCode()        <<
            "   "                                    <<
            std::setprecision(3)                     <<
            std::setw(8) << a.getAtomCoord().x       <<
            std::setw(8) << a.getAtomCoord().y       <<
            std::setw(8) << a.getAtomCoord().z       <<
            std::setprecision(2)                     <<
            std::setw(6) << a.getAtomRadius()        <<
            std::setw(6) << a.getAtomSASA()          <<
            "          "                             ;
            break;
        case SP:
            os << std::left << std::fixed << std::setprecision(6) <<
            a.getAtomCoord();
            break;
        case SPQ:
            os << std::left << std::fixed << std::setprecision(6) <<
            a.getAtomCoord() << " " <<
            a.getAtomCharge() ;
            break;
        case XYZ:
            os << a.getAtomSerial() << " " << 
            std::left << std::fixed << std::setprecision(6) <<
            a.getAtomCoord();
            break;
    }
    
    
	return os;
    }
}



