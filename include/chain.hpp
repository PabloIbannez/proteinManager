/*Pablo Ibáñez Freire, pablo.ibannez@uam.es*/

#ifndef CHAIN_HPP
#define CHAIN_HPP

#include <vector>

#include "proteinManager.hpp"

namespace proteinManager {

    class STRUCTURE;
    class MODEL;
    class RESIDUE;
    class ATOM;
    
    class CHAIN {
    
        friend STRUCTURE;
        friend MODEL;
        friend RESIDUE;
        friend ATOM;
    
        friend std::ostream& operator<<(std::ostream& os, const STRUCTURE& structure);
        friend std::ostream& operator<<(std::ostream& os, const MODEL& md);
        friend std::ostream& operator<<(std::ostream& os, const CHAIN& res);
        friend std::ostream& operator<<(std::ostream& os, const RESIDUE& res);
        friend std::ostream& operator<<(std::ostream& os, const ATOM& a);
    
      private:
      
        std::map<std::string, std::any> chainProperties;
    
        MODEL* parentModel_;
    
        boost::ptr_vector<RESIDUE> residueVector;
    
      public:
    
        CHAIN(std::string chainID,MODEL* parentModel);
        CHAIN(std::string chainID);
    
        /*
         ~CHAIN(){
        std::cerr << "Calling chain destructor" << std::endl;
        }*/
    
        RESIDUE& residue(int resSeq);
        boost::ptr_vector<RESIDUE>& residue();
        
        bool isRes(int resSeq);
    
        MODEL* getParentModel() const ;
        
        ////////////////////////////////////////////////////////////////
    
        #define GET_MDL_PROPERTY_T(type,Name)  GET_MDL_PROPERTY_R(type,Name)
        #define GET_MDL_PROPERTY_R(type,Name)  type getModel ## Name() const;
        #define GET_MDL_PROPERTY(r, data, tuple) GET_MDL_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
        
            MDL_PROPERTY_LOOP(GET_MDL_PROPERTY)
        
        #undef GET_MDL_PROPERTY_T
        #undef GET_MDL_PROPERTY_R
        #undef GET_MDL_PROPERTY
        
        ////////////////////////////////////////////////////////////////
        
        #define GET_CHAIN_PROPERTY_T(type,Name)  GET_CHAIN_PROPERTY_R(type,Name)
        #define GET_CHAIN_PROPERTY_R(type,Name)  type getChain ## Name() const;
        #define GET_CHAIN_PROPERTY(r, data, tuple) GET_CHAIN_PROPERTY_T(PROPTYPE(tuple),PROPNAME_CAPS(tuple))
        
            CHAIN_PROPERTY_LOOP(GET_CHAIN_PROPERTY)
        
        #undef GET_CHAIN_PROPERTY_T
        #undef GET_CHAIN_PROPERTY_R
        #undef GET_CHAIN_PROPERTY
    
        ////////////////////////////////////////////////////////////////
    
        #define SET_CHAIN_PROPERTY_T(Name,type,name)  SET_CHAIN_PROPERTY_R(Name,type,name)
        #define SET_CHAIN_PROPERTY_R(Name,type,name)  void setChain ## Name( type name ) ;
        #define SET_CHAIN_PROPERTY(r, data, tuple) SET_CHAIN_PROPERTY_T(PROPNAME_CAPS(tuple),PROPTYPE(tuple),PROPNAME(tuple))
        
            CHAIN_PROPERTY_LOOP(SET_CHAIN_PROPERTY)
        
        #undef SET_CHAIN_PROPERTY_T
        #undef SET_CHAIN_PROPERTY_R
        #undef SET_CHAIN_PROPERTY
        
        ////////////////////////////////////////////////////////////////
    
        void addResidue(std::string resName,int resSeq,std::string iCode);
    
        void addAtom(std::string resName, int resSeq, std::string iCode,
                     int serial, std::string name);
    };
}

#endif
