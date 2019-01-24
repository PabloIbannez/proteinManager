/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/

#ifndef ATOM_HPP
#define ATOM_HPP

#include "proteinManager.hpp"

namespace proteinManager {

    class STRUCTURE;
    class MODEL;
    class CHAIN;
    class RESIDUE;
    
    class ATOM {
    
        friend STRUCTURE;
        friend MODEL;
        friend CHAIN;
        friend RESIDUE;
    
        friend std::ostream& operator<<(std::ostream& os, const STRUCTURE& structure);
        friend std::ostream& operator<<(std::ostream& os, const MODEL& md);
        friend std::ostream& operator<<(std::ostream& os, const CHAIN& res);
        friend std::ostream& operator<<(std::ostream& os, const RESIDUE& res);
        friend std::ostream& operator<<(std::ostream& os, const ATOM& a);
    
      private:
      
        std::map<std::string, std::any> atomProperties;
    
        RESIDUE* parentResidue_;
    
      public:
    
        ATOM(int serial,std::string name,RESIDUE* parentResidue);
        ATOM(int serial,std::string name);
        
        /*
        ~ATOM(){
            std::cerr << "Calling atom destructor" << std::endl;
        }*/
        
        RESIDUE* getParentResidue() const;
        
        ////////////////////////////////////////////////////////////////////
        
        // A set of functions are declared to access the properties
        // of the model to which the atom belongs.
        
        #define GET_MDL_PROPERTY_T(type,Name)  GET_MDL_PROPERTY_R(type,Name)
        #define GET_MDL_PROPERTY_R(type,Name)  type getModel ## Name() const;
        #define GET_MDL_PROPERTY(r, data, tuple) GET_MDL_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
        
            MDL_PROPERTY_LOOP(GET_MDL_PROPERTY)
        
        #undef GET_MDL_PROPERTY_T
        #undef GET_MDL_PROPERTY_R
        #undef GET_MDL_PROPERTY
        
        ////////////////////////////////////////////////////////////////////
        
        // A set of functions are declared to access the properties
        // of the chain to which the atom belongs.
    
        #define GET_CHAIN_PROPERTY_T(type,Name)  GET_CHAIN_PROPERTY_R(type,Name)
        #define GET_CHAIN_PROPERTY_R(type,Name)  type getChain ## Name() const;
        #define GET_CHAIN_PROPERTY(r, data, tuple) GET_CHAIN_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
        
            CHAIN_PROPERTY_LOOP(GET_CHAIN_PROPERTY)
        
        #undef GET_CHAIN_PROPERTY_T
        #undef GET_CHAIN_PROPERTY_R
        #undef GET_CHAIN_PROPERTY
        
        ////////////////////////////////////////////////////////////////////
        
        // A set of functions are declared to access the properties
        // of the residue to which the atom belongs.
    
        #define GET_RES_PROPERTY_T(type,Name)  GET_RES_PROPERTY_R(type,Name)
        #define GET_RES_PROPERTY_R(type,Name)  type getRes ## Name() const;
        #define GET_RES_PROPERTY(r, data, tuple) GET_RES_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
        
        RES_PROPERTY_LOOP(GET_RES_PROPERTY)
        
        #undef GET_RES_PROPERTY_T
        #undef GET_RES_PROPERTY_R
        #undef GET_RES_PROPERTY
        
        ////////////////////////////////////////////////////////////////////
        
        // A set of functions are declared to access the properties
        // of the current atom.
    
        #define GET_ATOM_PROPERTY_T(type,Name)  GET_ATOM_PROPERTY_R(type,Name)
        #define GET_ATOM_PROPERTY_R(type,Name)  type getAtom ## Name() const;
        #define GET_ATOM_PROPERTY(r, data, tuple) GET_ATOM_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
        
            ATOM_PROPERTY_LOOP(GET_ATOM_PROPERTY)
        
        #undef GET_ATOM_PROPERTY_T
        #undef GET_ATOM_PROPERTY_R
        #undef GET_ATOM_PROPERTY
        
        ////////////////////////////////////////////////////////////////////
        
        // A set of functions are declared to set the properties
        // of the current atom.
        
        #define SET_ATOM_PROPERTY_T(Name,type,name)  SET_ATOM_PROPERTY_R(Name,type,name)
        #define SET_ATOM_PROPERTY_R(Name,type,name)  void setAtom ## Name( type name ) ;
        #define SET_ATOM_PROPERTY(r, data, tuple) SET_ATOM_PROPERTY_T(PROPNAME_CAPS(tuple),PROPTYPE(tuple),PROPNAME(tuple))
        
            ATOM_PROPERTY_LOOP(SET_ATOM_PROPERTY)
        
        #undef SET_ATOM_PROPERTY_T
        #undef SET_ATOM_PROPERTY_R
        #undef SET_ATOM_PROPERTY
    
    };
}

#endif
